from ..utils import load_api_config, greek2english
from ..ETM_API.references import Reference,DocMine
from ..ETM_API.references import AUTHORS, INSTITUTIONS,GRANT_APPLICATION,JOURNAL,SENTENCE,RELEVANCE
from scibite_toolkit.scibite_search import SBSRequestBuilder as s
import re
from collections import defaultdict
from time import sleep
from datetime import datetime

DATASET2IDTYPE = {'Medline':'PMID','PMC':'PMC','US Clinical Trials':'NCT ID'}

def remove_ids(text:str):
    pattern = r"\([^)]*\)"
    new_text = re.sub(pattern, "", text)
    new_text = new_text.replace('[','').replace(']','')
    return new_text


def parse_authors(authors:list):
    names = []
    affiliations = []
    for author in authors:
        # Get author names from name_normalized
        names.append(author.get('name_normalized', ''))
        # Get affiliations
        orgs = [affil for affil in author.get('affiliations', [])]
        affiliations += orgs # Append affiliations if they exist
        
    return list(filter(None, names)), list(filter(None, affiliations))


class SBSjson(DocMine):
    @staticmethod
    def dump_fname(search_name): return search_name + '.json'

    def __init__(self,doc:dict):
        id_type, native_id = self.__parse_ids(doc)
        super().__init__(id_type, native_id)

        try:
            title = remove_ids(doc['title'])
            self._set_title(title)
        except KeyError: pass

        try:
            last_modified = doc["last_modified_date"]
            year = last_modified[:4]
            self.set_date(year)
        except KeyError: pass

        authors = doc.get('authors',[])
        if authors:
            author_names,institutions = parse_authors(authors)
            if author_names: self[AUTHORS] = author_names
            if institutions: self[INSTITUTIONS] = institutions

        journal = doc.get('journal_title',GRANT_APPLICATION)
        self[JOURNAL] = journal

        issn = doc.get('issn','')
        if issn:
            self.append_property('ISSN', issn)

        # PUBLISHER?
        article_type = doc.get('article_type')
        if article_type:
            self['categories'] = [article_type]
        text_ref = self._make_textref()
        snippets = doc['_snippets'].values()
        self.add_sentence_props(text_ref, SENTENCE, snippets)
        self[RELEVANCE] = [doc['_score']]


    @staticmethod
    def __parse_ids(doc:dict):
        dataset, native_id = str(doc['dataset']), str(doc['native_id'])
        id_type = DATASET2IDTYPE.get(dataset,'')
        if not id_type:
            doi = str(doc.get('doi',''))
            if doi:
               id_type = 'PII' if doi.startswith('P0') else 'DOI' 
            else:
               id_type = dataset
        return id_type,native_id


class SBSapi():
  def __init__(self,*args,**kwargs):
      '''
      kwargs:
          api_type - defaults to "search"
      '''
      super().__init__()
      self.APIconfig = dict(args[0]) if args else load_api_config()
      self.limit = kwargs.get('limit',5)
      self.searchterm2id = dict()
      self.SBSsearch = self.SBSearch()


  def __try2getdocs(self,query='',markup=True,limit=20,offset=0,additional_fields=[]):
      sleep_time = 5
      for attempt in range(1,11):
          try:
              json_response = self.SBSsearch.get_docs(query,markup,limit,offset,additional_fields)
              if json_response:
                  return json_response
              else:
                  print('SBS server response is empty.  Will make {attempt} attempt in {sleep_time}')
                  sleep(sleep_time)
                  continue
          except Exception as ex:
              print(f'{attempt} attempt to retreive data for SBS {query} failed with {ex} exception')
              sleep(sleep_time)
              continue
      raise ex


  def SBSearch(self):  # This is actually regenerating a new builder request object.
      s.set_url(self.APIconfig['SBSurl'])
      s.set_auth_url(self.APIconfig['SBStoken_url'])
      s.set_oauth2(self.APIconfig['SBSuser'], self.APIconfig['SBSecret'])
      self.timestamp = datetime.now()
      return s  # This is actually regenerating a new builder request object.
  

  def get_id(self,search_term:str, in_vocabs:list=['GENETREE','EMTREE']):
      def __getid(prefix:str,vocab:str):
          try:
              return self.searchterm2id[search_term]
          except KeyError:
            response = self.SBSsearch.get_entities(suggest_prefix=prefix,include_vocabularies=[vocab],limit=1)
            data = response['data']   
            _id = data[0]['id'] if data else ''
            self.searchterm2id[search_term] = _id
            return _id
                          
      for vocab in in_vocabs:
          _id = __getid(search_term,vocab)
          if _id: return _id
      return ''


  def __hit_count(self,query:str):
      docs = self.__try2getdocs(query=query,limit=0)
      return docs['pagination']['totalItems']
  

  @staticmethod
  def deduplicate(refs:list[Reference]):
    def normalize(title:str):
      normal_title = re.sub("[^a-zA-Z]", "", title)
      return greek2english(normal_title).lower()
        
    #original_len = len(refs)
    normtitle2ref = dict()
    for ref in refs:
        title = ref.title()
        if title:
          normtitle = normalize(title)
          try:
              normtitle2ref[normtitle]._merge(ref)
          except KeyError:
              normtitle2ref[normtitle] = ref

    deduplicated_refs = list(normtitle2ref.values())
    #print(f'{original_len-len(deduplicated_refs)} references deduplicated')
    return deduplicated_refs


  def search(self,entity1:str,entity2:str,add2query:list)->tuple[int,dict[str,list],list[Reference]]:
    '''
    output tuple:
        [0] hit_count - TOTAL number of reference found by ETM basic search 
        [1] ref_ids = {id_type:[identifiers]}, where len(ref_ids) == ETMsearch.params['limit']\n
        id_type is from [PMID, DOI, 'PII', 'PUI', 'EMBASE','NCT ID']\n
        [2] references = [ref] list of Reference objects sorted by ETM relevance. len(references) == ETMsearch.params['limit'] 
        Relevance score is stored in ref['Relevance'] for every reference
    '''
    if (datetime.now() - self.timestamp).total_seconds() > 3000:
        self.SBSsearch = self.SBSearch()

    hit_count = 0
    references = list()
    ref_ids = defaultdict(list)

    fields = ['doi',"article_type","journal_title","issn","last_modified_date",'authors']
    step = 100
    entity1_id = self.get_id(entity1)
    if entity1_id:
        entity2_id = self.get_id(entity2)
        if entity2_id:
            query = f'{entity1_id} AND {entity2_id}'
            docs = self.__try2getdocs(query=query,limit=step,additional_fields=fields)
            hit_count =  docs['pagination']['totalItems']
            
            for offset in range(step,hit_count,step):
              for doc in docs['data']:
                  ref = SBSjson(doc)
                  references.append(ref)
                  id_type,identifier = ref.get_doc_id()
                  ref_ids[id_type].append(identifier)
              docs = self.__try2getdocs(query=query,limit=step,markup=False,additional_fields=fields,offset=offset)

            references = self.deduplicate(references)
    return len(references),dict(ref_ids),references

'''
if __name__ == "__main__":
    my_s = SBSapi()
    response_getdocs = s.__try2getdocs(query='abstract ~ INDICATION$D008175 AND DRUG$*',limit=100,markup=True)
    print(response_getdocs['data'])
'''
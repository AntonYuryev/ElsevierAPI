from ..utils import load_api_config, greek2english
from ..ETM_API.references import Reference,DocMine, Author
from ..ETM_API.references import AUTHORS, GRANT_APPLICATION,JOURNAL,SENTENCE,RELEVANCE,AUTHORS_STR
from scibite_toolkit.scibite_search import SBSRequestBuilder as s
import re
from collections import defaultdict
from time import sleep
from datetime import datetime

DATASET2IDTYPE = {'Medline':'PMID',
                  'PMC':'PMC',
                  'US Clinical Trials':'NCT ID',
                  'Embase':'EMBASE'}

def remove_ids(text:str):
    pattern = r"\([^)]*\)"
    new_text = re.sub(pattern, "", text)
    new_text = new_text.replace('[','').replace(']','')
    return new_text


def parse_authors(sbs_authors:list)->list[Author]:
  authors = list()
  names = list()
  for sbs_author in sbs_authors:
    # Get author names from name_normalized
    name = str(sbs_author.get('name_normalized', ''))
    if name:
      names.append(name)
      author = Author.fromStr(name)
      for affil in sbs_author.get('affiliations', []):
        institution_name = str(affil.get('source_text',''))
        if institution_name:
          address = str(affil.get('address',''))
          author.affiliations[institution_name] = address.replace('/',',')
        authors.append(author)
  return authors,names


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
        authors,names = parse_authors(authors)
        self[AUTHORS] = authors
        self[AUTHORS_STR] = ';'.join(names)


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
      self[RELEVANCE] = [float(doc['_score'])]


    @staticmethod
    def parse_doi(doi:str):
      if doi[:2] == 'P0':
        return 'PII',doi   
      s_pos = doi.find('/S')
      if s_pos > 0:
        return 'PII', doi[s_pos+1:]
      else:
        return 'DOI', doi


    @staticmethod
    def __parse_ids(doc:dict):
      dataset, native_id = str(doc['dataset']), str(doc['native_id'])
      try:
        return DATASET2IDTYPE[dataset],native_id
      except KeyError:
        try:
          return SBSjson.parse_doi(doc['doi'])
        except KeyError:
          return dataset, native_id
      

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
    self.SBSsearch = self.token_refresh()
    self.normtitle2ref = dict() # contains dictionary of normalized titles to refs for deduplication


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
        self.SBSsearch = self.token_refresh()
        continue
    raise ex


  def token_refresh(self):  # This is actually regenerating a new builder request object.
    sbs = s()
    sbs.set_url(self.APIconfig['SBSurl'])
    sbs.set_auth_url(self.APIconfig['SBStoken_url'])
    sbs.set_oauth2(self.APIconfig['SBSclientID'], self.APIconfig['SBSecret'])
    self.timestamp = datetime.now()
    return sbs
  

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
  

  def add_refs(self,refs:list[Reference])->dict[str,Reference]:
    '''
    output:
      {normalized_title:ref}
    '''
    def normalize(title:str):
      normal_title = re.sub("[^a-zA-Z]", "", title)
      return greek2english(normal_title).lower()
    
    normalized_titles = list()
    for ref in refs:
      title = ref.title()
      if title:
        normtitle = normalize(title)
        normalized_titles.append(normtitle)
        try:
          self.normtitle2ref[normtitle]._merge(ref)
        except KeyError:
          self.normtitle2ref[normtitle] = ref

    return {k:v for k,v in self.normtitle2ref.items() if k in normalized_titles}


  @staticmethod
  def process_chunk(docs:dict,refids:defaultdict[str,list[str]]):
    chunk_refs = list()
    for doc in docs['data']:
      ref = SBSjson(doc)
      chunk_refs.append(ref)
      id_type,identifier = ref.get_doc_id()
      if id_type:
        refids[id_type].append(identifier)
    return chunk_refs


  def search(self,entity1:str,entity2:str,add2query:list=[])->tuple[list[Reference],dict[str,list[str]]]:
    '''
    output:
      references found by search "entity1 AND entity2" deduplicated by titles and are added to self.normtitle2ref
      ref_ids = {id_type:[identifiers]}, where len(ref_ids) == ETMsearch.params['limit']   
    '''
    if (datetime.now() - self.timestamp).total_seconds() > 3000:
        self.SBSsearch = self.token_refresh()
        print(f'SBS token was refreshed at {datetime.now()}')

    hit_count = 0
    ref_ids = defaultdict(list)
    references = dict() # {mormalized_titles:ref}

    fields = ['doi',"article_type","journal_title","issn","last_modified_date",'authors']
    step = 100
    entity1_id = self.get_id(entity1)
    hit_count = 0
    if entity1_id:
      entity2_id = self.get_id(entity2)
      if entity2_id:
        query = f'{entity1_id} AND {entity2_id}'
        docs = self.__try2getdocs(query=query,limit=step,additional_fields=fields,markup=False)
        hit_count =  docs['pagination']['totalItems']
        for offset in range(step,hit_count,step):
          references.update(self.add_refs(self.process_chunk(docs,ref_ids)))
          docs = self.__try2getdocs(query=query,limit=step,markup=False,additional_fields=fields,offset=offset)
        references.update(self.add_refs(self.process_chunk(docs,ref_ids))) # to load last chunk
    # have to use len(references) instead of hit_count because of document duplication in SBS
    return list(references.values()),dict(ref_ids)
  

  def multiple_search(self,for_entity:str, and_concepts:list|set, add2query:list=[], ref_limit=10)->tuple[int,dict[str,list[str]],list[Reference]]:
    '''
    supports progressive hierachical concept search: leading concepts on the list have priority over trailing concepts
    '''
    references = list()
    refidtype2ids = defaultdict(set)
    for concept in and_concepts:
      if for_entity != concept: 
        refs,ref_ids = self.search(for_entity,concept,add2query)
        [refidtype2ids[idtype].update(ids) for idtype,ids in ref_ids.items()]
        [references.append(ref) for ref in refs if ref not in references] # preserve reference ranking
        if len(references) >= ref_limit:
          break

    refidtype2ids = {k:list(v) for k,v in refidtype2ids.items()}
    return len(references), dict(refidtype2ids), list(references)

'''
if __name__ == "__main__":
    s = SBSapi()
    #response_getdocs = s.__try2getdocs(query='abstract ~ INDICATION$D008175 AND DRUG$*',limit=100,markup=True)
    s.search('EGFR','cancer')
    #print(response_getdocs['data'])
'''

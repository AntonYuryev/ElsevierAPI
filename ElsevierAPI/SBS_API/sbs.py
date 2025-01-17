from ..utils import load_api_config, greek2english
from ..ETM_API.references import Reference,DocMine, Author
from ..ETM_API.references import AUTHORS, GRANT_APPLICATION,JOURNAL,SENTENCE,RELEVANCE,AUTHORS_STR
from scibite_toolkit.scibite_search import SBSRequestBuilder as s
from concurrent.futures import ThreadPoolExecutor, as_completed
import re,threading
from time import sleep
from datetime import datetime

DATASET2IDTYPE = {'Medline':'PMID',
                  'PMC':'PMC',
                  'SD CE':'PII',
                  'US Clinical Trials':'NCT ID',
                  'EU Clinical Trials':'EudraCT Number',
                  'Embase':'EMBASE'}

BM25SCORE = 'BM25score'
MAX_SBS_SESSIONS = 3
SBS_ID = 'sbs_id'


class SBSRef(DocMine):
  def __init__(self,id_type:str, native_id:str):
    super().__init__(id_type, native_id)


  def sbsid(self):
    return str(self.Identifiers[SBS_ID])


  def key(self):
    try: 
      return self['__key__']
    except KeyError:
      title = self.title()
      if title:
        key = (SBSRef.__normalize(title),self.pubyear())
      else:
        key = ('No title', 1812)

      self['__key__'] = key
      return key


  def __hash__(self):
    return hash(self.key())
  

  def __eq__(self, other:"SBSRef"):
    return self.__hash__() == other.__hash__()


  def _merge(self,other:"SBSRef"):
    other.pop('__key__','')
    super()._merge(other)


  @staticmethod
  def __docids(doc:dict):
    dataset, native_id = str(doc['dataset']), str(doc['native_id'])
    try:
      return DATASET2IDTYPE[dataset],native_id
    except KeyError:
      try:
        return SBSRef.parse_doi(doc['doi'])
      except KeyError:
        return dataset, native_id


  @staticmethod
  def __normalize(title:str):
    normal_title = re.sub("[^a-zA-Z]", "", title)
    return greek2english(normal_title).lower()


  @staticmethod
  def remove_markup(text:str): # use markup=False option to remove markup
    pattern = r"\([^)]*\)"
    new_text = re.sub(pattern, "", text)
    new_text = new_text.replace('[','').replace(']','')
    return new_text


  @staticmethod
  def parse_authors(sbs_authors:list)->tuple[list[Author],list[str]]:
    authors = []
    names = []
    for sbs_author in sbs_authors:
      # Get author names from name_normalized
      name = str(sbs_author.get('name_normalized', '')).strip()
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


  @classmethod
  def from_doc(cls,doc:dict,relevance_rank=1):
    doc_id = doc['id']
    sbsj = SBSRef(SBS_ID, doc_id)
    id_type, native_id = cls.__docids(doc)
    sbsj.Identifiers[id_type] = native_id

    try:
      sbsj._set_title(doc['title'])
    except KeyError: pass

    try:
      last_modified = doc["last_modified_date"]
      year = last_modified[:4]
      sbsj.set_date(year)
    except KeyError: pass

    authors = doc.get('authors',[])
    if authors:
      author_objs,author_names = SBSRef.parse_authors(authors)
      sbsj[AUTHORS] = author_objs
      sbsj[AUTHORS_STR] = [';'.join(author_names)]


    journal = doc.get('journal_title',GRANT_APPLICATION)
    sbsj[JOURNAL] = journal

    issn = doc.get('issn','')
    if issn:
      sbsj.append_property('ISSN', issn)

    # PUBLISHER?
    article_type = doc.get('article_type')
    if article_type:
      sbsj['categories'] = [article_type] # marker for loaded bibliography 
    
    try:
      sbsj[BM25SCORE] = [float(doc['_score'])/relevance_rank]
    except KeyError:
      pass

    return sbsj


  @classmethod
  def from_sent(cls,sent:dict,relevance_rank=1):
    dataset = str(sent['dataset'])
    native_id = str(sent['documentId'])[37:]
    id_type = DATASET2IDTYPE.get(dataset,dataset)
    field = sent['field']
    sentence = sent['content']

    sbsj = SBSRef(id_type, native_id)
    try:
      title = sent['title']
      sbsj._set_title(title)
    except KeyError: pass

    try:
      pybyear = sent['publishDate'][:4]
      sbsj.set_date(pybyear)
    except KeyError: pass

    sbsj[BM25SCORE] = [float(sent['_score'])/relevance_rank]
    
    sent_id = str(sent['id'])
    last_hyphen_index = sent_id.rfind('-')
    doc_id = sent_id[:last_hyphen_index]
    sbsj.Identifiers[SBS_ID] = doc_id

    textref = f'info:{id_type}/{native_id}#{field}:{sent_id[last_hyphen_index+1:]}'
    sbsj.add_sentence_prop(textref,SENTENCE,sentence)
    return sbsj


  def is_valid(self):
    try:
      authors = self[AUTHORS]
      return len(authors) > 0
    except KeyError:
      return any(w in self.Identifiers for w in ['NCT ID', 'EudraCT Number'])


  def has_bibliography(self):
    return not {'categories',AUTHORS,JOURNAL}.isdisjoint(self)
  

  @staticmethod
  def parse_doi(doi:str):
    if doi[:2] == 'P0':
      return 'PII',doi   
    s_pos = doi.find('/S')
    if s_pos > 0:
      return 'PII', doi[s_pos+1:]
    else:
      return 'DOI', doi


class SBSapi():
  def __init__(self,*args,**kwargs):
    '''
    kwargs:
      api_type - defaults to "search"
    '''
    super().__init__()
    self.APIconfig = dict(args[0]) if args else load_api_config()
    self.term2id = dict() #contains dictionary of search terms to ontology IDs
    self.SBSsearch = self.token_refresh()
    self.refCache = dict() # {ref_key:SBSRef}


  def clone(self):
    new_session = SBSapi(self.APIconfig)
    new_session.searchterm2id = self.term2id
    new_session.SBSsearch = self.token_refresh()
    new_session.refCache = self.refCache
    return new_session


  def _try2getdocs(self,query='',**kwargs):
    my_kwargs = {'markup':False,
                 'limit':self.limit,
                 'offset':0,
                 'maxSnippets':0,
                 'fields':[]}
    
    my_kwargs.update(kwargs)
    
    sleep_time = 5
    for attempt in range(1,11):
      try:
        json_response = self.SBSsearch.get_docs2(query,**my_kwargs)
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
  

  def __sent2doc(self,sent_ref:SBSRef):
    assert (not sent_ref.has_bibliography())
    doc = self.SBSsearch.get_document(sent_ref.sbsid(),markup=False)
    ref = SBSRef.from_doc(doc['data'])
    sent_ref._merge(ref)
    assert(self.refCache[sent_ref.key()].has_bibliography())
    return sent_ref if ref.is_valid() else None


  def sents2docs(self,sent_refs:list[SBSRef]):
    '''
    output:
      list of clean references with bibliography
    '''
    futures = []
    with ThreadPoolExecutor(max_workers=MAX_SBS_SESSIONS, thread_name_prefix='__sent2doc') as e:
      [futures.append(e.submit(self.__sent2doc,r)) for r in sent_refs]

    refs = []
    for f in as_completed(futures):
      merged_ref = f.result()
      if merged_ref is not None:
        refs.append(merged_ref)

    return refs
  

  def _try2getsents(self,query:str,**kwargs)->tuple[int,list[SBSRef]]:
    '''
    output:
      hit_count - to continue retreival with offset
      liste of retreived SBSRef objects
    '''
    my_kwargs = {'markup':False,
                 'limit':100, # cannot be more than 100
                 'offset':0,
                 'hit_count':0,
                 'relevance_rank':1}
    
    my_kwargs.update(kwargs)
    hit_count = my_kwargs.pop('hit_count',0)
    sleep_time = 5
    for attempt in range(1,11):
      try:
        json_response = self.SBSsearch.get_sentences2(query,**my_kwargs)
        if 'data' not in json_response: # bug in SBS
          limit = my_kwargs['limit']
          if hit_count: # case when hit count is known but is incorrect
            leftover = hit_count-my_kwargs['offset']
            if leftover < limit:
              limit = leftover
            while 'error' in json_response: #usually the real hit count is close to incorrect hit count
              limit -= 1
              my_kwargs['limit'] = limit
              json_response = self.SBSsearch.get_sentences2(query,**my_kwargs)
                 
          while 'error' in json_response: # case when hit count is not known
            limit -= 10
            my_kwargs['limit'] = limit
            json_response = self.SBSsearch.get_sentences2(query,**my_kwargs)
            
          while 'error' not in json_response:
            limit += 5
            my_kwargs['limit'] = limit
            json_response = self.SBSsearch.get_sentences2(query,**my_kwargs)

          limit -= 5
          while 'error' in json_response:
            limit -= 1
            my_kwargs['limit'] = limit
            json_response = self.SBSsearch.get_sentences2(query,**my_kwargs)
            
        data = json_response['data']
        if data: # data may be empty
          hit_count = json_response['pagination']['totalItems']
          sent_refs = [SBSRef.from_sent(json_response['data'][0])]
          for i in range(1,len(json_response['data'])):
            new_ref = SBSRef.from_sent(json_response['data'][i])
            if sent_refs[-1] == new_ref:
              sent_refs[-1]._merge(new_ref)
            else:
              sent_refs.append(new_ref)
          return hit_count, sent_refs
        else:
          return 0,[]
      except Exception as ex:
        print(f'{attempt} attempt to retreive data for SBS query {query} failed with {ex} exception',flush=True)
        sleep(sleep_time)
        self.SBSsearch = self.token_refresh()
        continue
    raise ex
  

  def get_id(self,search_term:str, in_vocabs:list=['GENETREE','EMTREE'])->str:
    def __getid(prefix:str,vocab:str):
      try:
          return self.term2id[search_term]
      except KeyError:
        response = self.SBSsearch.get_entities(suggest_prefix=prefix,include_vocabularies=[vocab],limit=1)
        data = response['data']   
        _id = data[0]['id'] if data else ''
        self.term2id[search_term] = _id
        return str(_id)
                        
    with threading.Lock():
      for vocab in in_vocabs:
        _id = __getid(search_term,vocab)
        if _id: return str(_id)
    return ''
  

  def __search_term(self,search_term:str, in_vocabs:list=['GENETREE','EMTREE']):
    term_id = self.get_id(search_term,in_vocabs)
    return term_id if term_id else f'"{search_term}"'
  

  def __map2terms(self,terms:list[str],in_vocabs:list=['GENETREE','EMTREE']):
    return list(filter(None,(map(lambda t: self.__search_term(t, in_vocabs),terms))))


  def __hit_count(self,query:str):
      docs = self._try2getdocs(query=query,limit=0)
      return docs['pagination']['totalItems']
  
  
  def add_refs(self,refs:list[SBSRef],from_sent=False)->list[SBSRef]:
    with threading.Lock(): # to prevent lock during multithreading
      added_refs = set()
      new_refs = []
      for ref in refs:
        ref_key = ref.key()
        if ref_key in self.refCache:
          exist_ref = self.refCache[ref_key]
          exist_ref._merge(ref)
          added_refs.add(exist_ref)
        else:
          self.refCache[ref_key] = ref
          added_refs.add(ref)
          new_refs.append(ref)

      if from_sent:
        self.sents2docs(new_refs)
      
      return list(added_refs)
  

  def bm25_Relevance(self):
    for ref in self.refCache.values():
      ref[RELEVANCE] = [ref.relevance(BM25SCORE)*(1.0+ref.number_of_sentences()/2.0)]


  def conjunct2lists(self,entities1:list,entities2:list,add2query=[]):
    entity1_ids = self.__map2terms(entities1)
    if entity1_ids:
      entity2_ids = self.__map2terms(entities2)
      if entity2_ids:
        str1 = ' OR '.join(entity1_ids)
        str2 = ' OR '.join(entity2_ids)
        query = f'({str1}) AND ({str2})'
        if add2query:
          query += ' AND '.join(self.__map2terms(add2query))
        return query
    return ''
  

  def search_sents(self,query:str,relevance_rank=1)->tuple[int,list[SBSRef]]:
    if (datetime.now() - self.timestamp).total_seconds() > 3000:
        self.SBSsearch = self.token_refresh()
        print(f'SBS token was refreshed at {datetime.now()}')

    if query:
      kwargs = {'relevance_rank':relevance_rank,'limit':100}
      hit_count,sent_refs = self._try2getsents(query,**kwargs)
      return hit_count,self.add_refs(sent_refs,from_sent=True)
    else:
      return 0,[]


  def __cooc4list(self,entities:list[str], and_concepts:list[str]|set[str],
                add2query:list[str]=[])->dict[str,tuple[int,list[SBSRef]]]:
    results = dict()
    for term in entities:
      term_concepts = [x for x in and_concepts if x != term]
      if term not in term_concepts:
        query = self.conjunct2lists([term],term_concepts,add2query)
        hit_count,refs = self.search_sents(query)
        results[term] = (hit_count,refs)
    return results
  

  def cooc4list(self,entities:list[str], link2concepts:list[str]|set[str],
    add2query:list[str]=[],multithread=True)->dict[str,tuple[int,list[Reference]]]:
    '''
    output:
      {entity:(hit_count,[SBSRef])}
    '''
    entity2rfks = dict()
    if multithread:
      chunk_size = int(len(entities)/MAX_SBS_SESSIONS)
      with ThreadPoolExecutor(max_workers=MAX_SBS_SESSIONS, thread_name_prefix='SBS') as e:
        futures = list()
        new_session = self.clone()
        for chunk in range(0,len(entities),chunk_size):
          fut_entities = entities[chunk:chunk+chunk_size]
          future = e.submit(new_session.__cooc4list,fut_entities,link2concepts,add2query)
          futures.append(future)
        
        for f in as_completed(futures):
          results = f.result()
          entity2rfks.update(results)
    else:
      entity2rfks = self.__cooc4list(entities,link2concepts,add2query)

    self.bm25_Relevance()

    e2refs = dict()
    for entity in entities:
      hit_count,refs = entity2rfks.get(entity,(0,[]))
      clean_refs = [r for r in refs if r.is_valid()]
      clean_refs.sort(key=lambda x:x.relevance(),reverse=True)
      e2refs[entity] = (hit_count,clean_refs)

    return e2refs


  @staticmethod
  def process_chunk(docs:dict,relevance_rank=1,from_sent=False)->list[SBSRef]:
    chunk_refs = list()
    make_ref = SBSRef.from_sent if from_sent else SBSRef.from_doc
    try:
      data = docs['data'] # sometime server error is returned
      for doc in data:
        ref = make_ref(doc,relevance_rank)
        chunk_refs.append(ref)
      return chunk_refs
    except KeyError:
      print(f'Error in {docs}')
      return []


  def search_docs(self,query:str,relevance_rank=1)->dict[tuple[str,int],Reference]:
    '''
    output:
      references found by search "entity1 AND entity2" deduplicated by titles and are added to self.refCache
      ref_ids = {id_type:[identifiers]}, where len(ref_ids) == ETMsearch.params['limit']   
    '''
    if (datetime.now() - self.timestamp).total_seconds() > 3000:
        self.SBSsearch = self.token_refresh()
        print(f'SBS token was refreshed at {datetime.now()}')

    refkeys = set() # {(normalized_title,pubyear)}
    fields = ['doi',"article_type","journal_title","issn","last_modified_date",'authors']
    if query:
      step = 100
      docs = self._try2getdocs(query,limit=step,fields=fields)
      hit_count =  docs['pagination']['totalItems']
      for offset in range(step,hit_count,step):
        step_refs = self.process_chunk(docs,relevance_rank)
        step_refkeys = self.add_refs(step_refs)
        refkeys.update(step_refkeys)
        docs = self._try2getdocs(query,limit=step,fields=fields,offset=offset)
      refkeys.update(self.add_refs(self.process_chunk(docs))) # to load last chunk
    return refkeys

'''
if __name__ == "__main__":
    s = SBSapi()
    #response_getdocs = s._try2getdocs(query='abstract ~ INDICATION$D008175 AND DRUG$*',limit=100,markup=True)
    s.search('EGFR','cancer')
    #print(response_getdocs['data'])
'''

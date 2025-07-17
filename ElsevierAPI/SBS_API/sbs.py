from ..utils import load_api_config, greek2english,urlencode,json,multithread,execution_time
from ..ETM_API.references import Reference,DocMine, Author
from ..ETM_API.references import AUTHORS,_AUTHORS_,GRANT_APPLICATION,JOURNAL,SENTENCE,RELEVANCE
from scibite_toolkit.scibite_search import SBSRequestBuilder as s
from concurrent.futures import ThreadPoolExecutor, as_completed
import re,threading
from time import sleep, time
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
      sbsj[_AUTHORS_] = author_objs
      sbsj[AUTHORS] = [';'.join(author_names)]


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
    '''
    checks:
      if self has _AUTHORS_ or valid identifiers: 'NCT ID','EudraCT Number','PMC','PMID','PII','DOI'
    '''
    try:
      authors = self[_AUTHORS_]
      return len(authors) > 0
    except KeyError:
      return any(w in self.Identifiers for w in ['NCT ID', 'EudraCT Number','PMC','PMID','PII','DOI'])


  def has_bibliography(self):
    return not {'categories',_AUTHORS_,JOURNAL}.isdisjoint(self)
  

  @staticmethod
  def parse_doi(doi:str):
    if doi[:2] == 'P0':
      return 'PII',doi   
    s_pos = doi.find('/S')
    if s_pos > 0:
      return 'PII', doi[s_pos+1:]
    else:
      return 'DOI', doi


####################### SBSapi ############# SBSapi ###################### SBSapi ##############
class SBSapi():
  def __init__(self,api_config:dict=dict()):
    self.APIconfig = api_config.copy() if api_config else load_api_config()
    self.term2id = dict() #contains dictionary of search terms to ontology IDs
    self.refCache = dict() # {ref_key:SBSRef}
    self.timestamp = datetime.now()
    self.SBSsearch = self.__get_token()
    

  def clone(self):
    new_session = SBSapi(self.APIconfig)
    new_session.term2id = self.term2id
    new_session.refCache = self.refCache
    return new_session


  def __get_token(self):  # This is actually regenerating a new builder request object.
    sbs = s()
    sbs.set_url(self.APIconfig['SBSurl'])
    sbs.set_auth_url(self.APIconfig['SBStoken_url'])
    sbs.set_oauth2(self.APIconfig['SBSclientID'], self.APIconfig['SBSecret'])
    self.timestamp = datetime.now()
    return sbs
  

  def token_refresh(self):
    if (datetime.now() - self.timestamp).total_seconds() > 3000:
      self.SBSsearch = self.__get_token()
      print(f'SBS token was refreshed at {datetime.now()}')
  

  def get_id(self,search_term:str, in_vocabs:list=['GENETREE','EMTREE'])->str:
    with threading.Lock():
      for vocab in in_vocabs:
        response = self.SBSsearch.get_entities(suggest_prefix=search_term,include_vocabularies=[vocab],limit=1)
        if 'data' not in response: continue
        data = response['data']
        if not data: continue
        _id = str(data[0]['id'])
        self.term2id[search_term] = _id
        return _id
    return ''
  

  def __search_term(self,search_term:str, in_vocabs:list=['GENETREE','EMTREE'])->str:
    # quoted search_term is a signal to search it as is without ontology synonym expansion
    # + term expansion sign is stripped from the end of search_term
    if search_term[0] == '"' and search_term[-1] == '"':
      if len(search_term) > 2:
        return search_term
      else:
        print(f'{search_term} is empty')
        return ''
    else:
      if search_term.endswith('+'):
        noplus_term = search_term.rstrip(' +')
        expand = ' +'
      else:
        noplus_term = search_term
        expand = ''

      if noplus_term in self.term2id:
        return self.term2id[noplus_term]+expand
      else:
        term_id = self.get_id(noplus_term,in_vocabs)
        if not term_id:
          term_id = f'"{noplus_term}"'
          self.term2id[noplus_term] = term_id
          return term_id
        return term_id + expand
  

  def __map2terms(self,terms:list[str],in_vocabs:list=['GENETREE','EMTREE']):
    '''
    handles search terms with quotes and with + signs at the end
    '''
    return list(filter(None,(map(lambda t: self.__search_term(t, in_vocabs),terms))))
  
  
  def add_refs(self,fetched_refs:list[SBSRef],from_sent=False)->list[SBSRef]:
    '''
    output:
      equivalent to fetched_refs [SBSRef] from self.refCache with updated bibliography and all sentences
    '''
    with threading.Lock(): # to prevent lock during multithreading
      cache_equivalent_refs = set()
      new_refs = []
      for ref in fetched_refs:
        ref_key = ref.key()
        if ref_key in self.refCache:
          exist_ref = self.refCache[ref_key]
          exist_ref._merge(ref)
          cache_equivalent_refs.add(exist_ref)
        else:
          self.refCache[ref_key] = ref
          cache_equivalent_refs.add(ref)
          new_refs.append(ref)

      if from_sent:
        self.sents2docs(new_refs)
      
      return list(cache_equivalent_refs)
  

  def bm25_Relevance(self):
    '''
    adjusts ref[RELEVANCE] in self.refCache by ref.number_of_sentences()
    '''
    for ref in self.refCache.values():
      ref[RELEVANCE] = [ref.relevance(BM25SCORE)*(1.0+ref.number_of_sentences()/2.0)]


  def search_url(self,entities:list[str], link2concepts:list[str]|set[str])->str:
    #https://psgscibitesearch.lifesciencepsg.com/documents?document_schema=journal&ssql=[sptbn1](GENETREE$GENETREE_6711) [AND](AND) [dopamine](GENETREE$GENETREE_1003345)
    url = self.APIconfig['SBSurl'] + '/documents?document_schema=journal&ssql='
    entity_ids = self.__map2terms(entities)
    concept_ids = self.__map2terms(link2concepts)
    entities_str = ' [OR](OR) '.join([f'[{entity}]({entity_id})' for entity,entity_id in zip(entities,entity_ids)])
    concepts_str = ' [OR](OR) '.join([f'[{concept}]({concept_id})' for concept,concept_id in zip(link2concepts,concept_ids)])
    return urlencode(url + entities_str + ' [AND](AND) ' + concepts_str)
  

  def conjunct2lists(self,entities1:list[str],entities2:list[str],add2query:list[str]=[])->tuple[str,str]:
    '''
    output:
      query - conjunct search query for SBS
      encoded_search_url - url for SBS search
    '''
    entity1_ids = self.__map2terms(entities1)
    entity2_ids = self.__map2terms(entities2)
    str1 = ' OR '.join(entity1_ids)
    str2 = ' OR '.join(entity2_ids)

    entities1_search_str = ' [OR](OR) '.join([f'[{entity}]({entity_id})' for entity,entity_id in zip(entities1,entity1_ids)])
    entities2_search_str = ' [OR](OR) '.join([f'[{concept}]({concept_id})' for concept,concept_id in zip(entities2,entity2_ids)])
    
    search_url = f'[(](() {entities1_search_str} [)]()) [AND](AND) [(](() {entities2_search_str} [)]())'
    query = f'({str1}) AND ({str2})'

    if add2query:
      add2query_ids = self.__map2terms(add2query)
      query += ' AND (' + 'OR'.join(add2query_ids)+')'
      add2query_search_str = ' [OR](OR) '.join([f'[{entity}]({entity_id})' for entity,entity_id in zip(add2query,add2query_ids)])
      search_url += f' [AND](AND) [(](() {add2query_search_str} [)]())'

    base_url = self.APIconfig['SBSurl'] + '/documents?document_schema=journal&ssql='
    encoded_search_url = base_url+urlencode(search_url)
    return query,encoded_search_url


  def join2query(self,entity:str,concepts:list[str],add2query:list[str]=[])->tuple[str,str]:
    '''
      output:
        query - conjunct search query for SBS
        encoded_search_url - url for SBS search
      if entity is in double quotes it will be searched as is without ontology synonym expansion
    '''
    my_concepts = [c for c in concepts if c != entity]
    return self.conjunct2lists([entity],my_concepts,add2query)


################# SENTENCE SEARCH ############### SENTENCE SEARCH ##########################
  def __sent2doc(self,sent_ref:SBSRef):
    '''
    loads bibliography for sent_ref
    '''
    assert (not sent_ref.has_bibliography())
    senref_sbsid = sent_ref.sbsid()
    doc = self.SBSsearch.get_document(senref_sbsid,markup=False)
    if 'data' in doc:
      ref = SBSRef.from_doc(doc['data'])
      sent_ref._merge(ref)
      assert(self.refCache[sent_ref.key()].has_bibliography())
      return sent_ref if ref.is_valid() else None
    else:
      print(f'Invalid document: {str(self.APIconfig['SBSurl']).rstrip('/')+'/documents/'+senref_sbsid}\n')
      return None


  def sents2docs(self,sent_refs:list[SBSRef]):
    '''
    input:
      SBS sentence documents
    output:
      list of clean references with bibliography
    '''
    refs = multithread(sent_refs,self.__sent2doc,max_workers=MAX_SBS_SESSIONS)
    refs = list(filter(None,refs))
    return refs
  

  def __get_sentences2__(self,query,**my_kwargs):
    try:
      return self.SBSsearch.get_sentences2(query,**my_kwargs)
    except json.JSONDecodeError: # get_sentences2 raises error only if token is expired
      self.token_refresh()
      return self.__get_sentences2__(query,**my_kwargs)


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
              if limit <= 0: break
              my_kwargs['limit'] = limit
              json_response = self.SBSsearch.get_sentences2(query,**my_kwargs)
              
          while 'error' in json_response: # case when hit count is not known
            limit -= 10
            if limit <= 0: break
            my_kwargs['limit'] = limit
            json_response = self.SBSsearch.get_sentences2(query,**my_kwargs)
            
          while 'error' not in json_response:
            limit += 5
            if limit > 100: break
            my_kwargs['limit'] = limit
            json_response = self.SBSsearch.get_sentences2(query,**my_kwargs)

          limit -= 5
          while 'error' in json_response:
            limit -= 1
            if limit <= 0: break
            my_kwargs['limit'] = limit
            json_response = self.SBSsearch.get_sentences2(query,**my_kwargs)
            
        if 'data' in json_response:
          data = json_response['data']
          if data: # data may be empty
            sentence_count = json_response['pagination']['totalItems']
            sent_refs = [SBSRef.from_sent(json_response['data'][0])]
            for i in range(1,len(data)):
              new_ref = SBSRef.from_sent(json_response['data'][i])
              if sent_refs[-1] == new_ref: # sentence is from the same document
                sent_refs[-1]._merge(new_ref)
              else:
                sent_refs.append(new_ref)
            return sentence_count, sent_refs
          else:
            return 0,[]
        else:
          return 0,[]
      except Exception as ex:
        print(f'{attempt} attempt to retreive data for SBS query {query} failed with {ex} exception',flush=True)
        sleep(5)
        self.token_refresh()
        continue
    raise ex
  
  
  def search_sents(self,query:str,relevance_rank=1)->tuple[int,list[SBSRef]]:
    '''
    output:
        total_sentence_count, [SBSRef] for first 100 sentences
    '''
    self.token_refresh()

    if query:
      kwargs = {'relevance_rank':relevance_rank,'limit':100}
      sentence_count,sent_refs = self._try2getsents(query,**kwargs)
      return sentence_count,self.add_refs(sent_refs,from_sent=True)
    else:
      return 0,[]


  def __sentcooc4list(self,entities:list[str], and_concepts:list[str]|set[str],
                add2query:list[str]=[])->dict[str,tuple[str,int,list[SBSRef]]]:
    '''
    input:
      quoted entities are searched as is without ontology synonym expansion
    output:
      references for each search "entity and_concepts" as {entity:(sentence_count,[SBSRef])}
    '''
    results = dict()
    for entity in entities:
      query,search_url = self.join2query(entity,and_concepts,add2query)
      if entity == entities[0]:
        print(f'Will find sentence co-occurence for {len(entities)} entities')
        print(f'Sample query for sentence co-occurence: {query}\n')
      sentence_count,refs = self.search_sents(query)
      results[entity] = (search_url,sentence_count,refs)
    return results
  

  def sentcooc4list(self,entities:list[str], link2concepts:list[str]|set[str],
    add2query:list[str]=[],multithread=True)->dict[str,tuple[str,int,list[Reference]]]:
    '''
    Entry function for RefStat.add_refs_sbs\n
    output:
      {entity:(sentence_count,[SBSRef])}, where RELEVANCE is adjusted by number of sentences in every SBSRef
    '''
    start = time()
    entity2rfks = dict()
    if multithread:
      chunk_size = int(len(entities)/MAX_SBS_SESSIONS)
      with ThreadPoolExecutor(max_workers=MAX_SBS_SESSIONS, thread_name_prefix='SBS') as e:
        futures = list()
        new_session = self.clone()
        for chunk in range(0,len(entities),chunk_size):
          fut_entities = entities[chunk:chunk+chunk_size]
          future = e.submit(new_session.__sentcooc4list,fut_entities,link2concepts,add2query)
          futures.append(future)
        
        for f in as_completed(futures):
          results = f.result()
          entity2rfks.update(results)
    else:
      entity2rfks = self.__sentcooc4list(entities,link2concepts,add2query)

    self.bm25_Relevance()

    e2refs = dict()
    rows_counter = 0
    for entity in entities:
      search_url,sentence_count,refs = entity2rfks.get(entity,('',0,[]))
      rows_counter += bool(sentence_count)
      clean_refs = [r for r in refs if r.is_valid()]
      clean_refs.sort(key=lambda x:x.relevance(),reverse=True)
      e2refs[entity] = (search_url,sentence_count,clean_refs)

    print(f'Found sentences in SBS for {rows_counter} rows in {execution_time(start)}')
    return e2refs

################# DOCUMENT SEARCH ####################### DOCUMENT SEARCH ###############
  def _try2getdocs(self,query='',**kwargs)->tuple[int,list[SBSRef]]:
    '''
    output:
      hit_count, references
    '''
    my_kwargs = {'markup':False,
                 'limit':self.limit,
                 'offset':0,
                 'maxSnippets':0,
                 'fields':[]}
    my_kwargs.update(kwargs)
    
    json_response = self.get_docs2(query,**my_kwargs)
    if json_response and 'data' in json_response:
      data = json_response['data']
      if data: # data may be empty
        hit_count = int(data['pagination']['totalItems'])
        refs = [SBSRef.from_doc(data[i]) for i in range(0,len(data))]
        # refs are between offset and offset+limit
        return hit_count, refs
    return 0,[]


  def search_docs(self,query:str,relevance_rank=1)->list[SBSRef]:
    '''
    output:
      references found by search "entity1 AND entity2" deduplicated by titles and are added to self.refCache
      ref_ids = {id_type:[identifiers]}, where len(ref_ids) == ETMsearch.params['limit']   
    '''
    self.token_refresh()
    fields = ['doi',"article_type","journal_title","issn","last_modified_date",'authors']
    if query:
      step = 100
      hit_count,_1strefs = self._try2getdocs(query,fields=fields,limit=step)
      all_refs = set(self.add_refs(_1strefs))
      for offset in range(step,hit_count,step):
        hit_count,step_refs = self._try2getdocs(query,fields=fields,offset=offset,limit=step)
        all_refs.update(self.add_refs(step_refs))
    return all_refs


  def get_docs2(self,query,**kwargs):
    try:
      return self.SBSsearch.get_docs2(query,**kwargs)
    except json.JSONDecodeError: # get_docs2 raises error only if token expires
      self.token_refresh()
      return self.get_docs2(query,**kwargs)


  def __abscooc4list(self,entities:list[str],concepts:list[str])->dict[str,tuple[int,float]]:
    '''
    input:
      quoted entities are searched as is without ontology synonym expansion
    output:
      references found by search "entity1 AND entity2" deduplicated by titles and are added to self.refCache
      ref_ids = {id_type:[identifiers]}, where len(ref_ids) == ETMsearch.params['limit']   
    '''
    self.token_refresh()
    entity2abscount = dict()
    number_of_entities = len(entities)
    #abs_query = 'dataset = "medline" AND field = "abstract" AND content ~ ({q})'
    abs_query = 'dataset = "medline" AND abstract ~ ({q})'
    for entity in entities:
      jq,_ = self.join2query(entity,concepts)
      if jq:
        limit = 100
        query = abs_query.format(q=jq)
        if entity == entities[0]:
          print(f'\nWill find abstract coocurence for {number_of_entities} entities')
          print(f'Sample query for abstract co-ocurence: {query}\n')
        kwargs = {'markup':False,
                 'limit':limit,
                 'offset':0,
                 'maxSnippets':0,
                 'fields':[]}
        abstract_count = 0
        abstract_relevance = 0.0
        json_response = self.get_docs2(query,**kwargs)
        if json_response and 'data' in json_response:
          data = json_response['data']
          if data:
            hit_count = int(json_response['pagination']['totalItems'])
            abstract_count = hit_count
            abstract_relevance = sum([a['_score'] for a in data])
          
            for offset in range(limit,hit_count,limit):
              kwargs['offset'] = offset
              json_response = self.get_docs2(query,**kwargs)
              if json_response and 'data' in json_response :
                abstract_relevance += sum([a['_score'] for a in json_response['data']])
        else:
          print(f'No abstract co-occurence for query {query}')
          continue

      entity2abscount[entity]=(abstract_count,abstract_relevance)
    return entity2abscount
  

  def abscooc4list(self,entities:list[str],link2concepts:list[str]|set[str],
                   multithread=True)->dict[str,tuple[int,float]]:
    '''
    Entry function for RefStat.add_abs_cooc\n
    '''
    start = time()
    entity2stat = dict()
    if multithread:
      chunk_size = int(len(entities)/MAX_SBS_SESSIONS)
      with ThreadPoolExecutor(max_workers=MAX_SBS_SESSIONS, thread_name_prefix='SBS') as e:
        futures = list()
        new_session = self.clone()
        for chunk in range(0,len(entities),chunk_size):
          fut_entities = entities[chunk:chunk+chunk_size]
          future = e.submit(new_session.__abscooc4list,fut_entities,link2concepts)
          futures.append(future)
        
        for f in as_completed(futures):
          results = f.result()
          entity2stat.update(results)
    else:
      entity2stat = self.__abscooc4list(entities,link2concepts)

    print(f'Abstract co-occurence for {len(entities)} entities finished in {execution_time(start)}')
    return entity2stat


if __name__ == "__main__":
    s = SBSapi()
    #response_getdocs = s._try2getdocs(query='abstract ~ INDICATION$D008175 AND DRUG$*',limit=100,markup=True)
    SBSsearch = s.token_refresh()
    SBSsearch.get_sentences(query='GENETREE$GENETREE_3000051 AND GENETREE$GENETREE_12813384', markup=False, limit=100, offset=0)
    #s.search('EGFR','cancer')
    #print(response_getdocs['data'])


from ..pandas.panda_tricks import df,pd
from ..utils import execution_time,sortdict,load_api_config
from ..NCBI.pubmed import medlineTA2issn
from ..NCBI.PMC import docid2pmid
from ..SBS_API.sbs import SBSapi
from .etm import ETMsearch
from ..ScopusAPI.scopus import Scopus,AuthorSearch
from .references import Reference,DocMine,pubmed_hyperlink,make_hyperlink,pmc_hyperlink,pii_hyperlink,doi_hyperlink
from ..ScopusAPI.scopus import SCOPUS_AUTHORIDS,SCOPUS_CITESCORE,SCOPUS_SJR,SCOPUS_SNIP
from .references import AUTHORS,INSTITUTIONS,JOURNAL,PUBYEAR,RELEVANCE,ETM_CITATION_INDEX,IN_OPENACCESS,PUBLISHER,GRANT_APPLICATION
import math,time,json,os,zipfile
from datetime import datetime,date
from concurrent.futures import ThreadPoolExecutor
from urllib.error import HTTPError
from collections import defaultdict


IDENTIFIER_COLUMN = 'Document identifier: PMID or DOI'
MAX_TM_SESSIONS = 10 #ETM performance deteriorates if number of concurrent sessions is > 5

class RefStats:
  """
  self.ref_counter = {str(id_type+':'+identifier):(ref,count)}
  """
  def __init__(self,APIconfig:dict,**kwargs):
      """
        self.ref_counter = {(id_type,identifier):(ref,count)}
      """
      my_kwargs = {'load_medlineTA': True,'limit':5}
      my_kwargs.update(kwargs)
      self.ref_counter = dict() # {(id_type,identifier):(ref,count)}
      self.refcols = set() # stores refcount columns names for formatting by SemanticSearch.clean_df()
      self.doi_columns = set()

      self.APIconfig = APIconfig if APIconfig else load_api_config()
      self.AuthorSearch = AuthorSearch(self.APIconfig)
      self.Scopus = Scopus(self.APIconfig)
      self.ref_limit = my_kwargs.get('limit',5)
      self.NormalizedAffs = dict()

      if my_kwargs['load_medlineTA']:
        self.abbrev2journal, self.journal2issn = medlineTA2issn() # {journal_title:"[issn_print,issn_online]}
      

  def clone(self):
    kwargs = {'load_medlineTA':False,'limit':self.ref_limit}
    refstat = RefStats(self.APIconfig,**kwargs)
    refstat.abbrev2journal = self.abbrev2journal.copy()
    refstat.journal2issn = self.journal2issn.copy()
    refstat.NormalizedAffs = self.NormalizedAffs.copy()
    return refstat
  

  def _limit(self):
      return self.ref_limit


  def references(self)->set[DocMine]:
      return set([x[0] for x in self.ref_counter.values()])


  def _add2counter(self, ref:DocMine):
    '''
    updates:
      self.ref_counter - {str(id_type+':'+identifier):(ref,count)}
      ref[RELEVANCE] with ref.relevance()

    output:
      id_type,identifier of added ref
    '''
    counter_key = ref.get_doc_id()
    if counter_key in self.ref_counter:
      count_exist = self.ref_counter[counter_key][1]
      self.ref_counter[counter_key] = (ref, count_exist+1)
      self.ref_counter[counter_key][0][RELEVANCE][0] += ref.relevance()
    else:
      self.ref_counter[counter_key] = (ref,1)


  def counter2df(self, use_relevance=True):
    """
    Creates DataFrame from self.ref_counter = {identifier:(ref,ref_count)}
    used to count references from ETM
    """
    if not self.ref_counter: 
      print('No TM references found')
      return df()
    
    first_col = RELEVANCE if use_relevance else ETM_CITATION_INDEX
    header = [first_col,PUBYEAR,'Identifier type',IDENTIFIER_COLUMN,'Citation']

    # normalization constant for relevance score between 0 and 100
    log_scores = []
    for ref, refcount in self.ref_counter.values():
      boosted_relevance = ref.relevance() * int(refcount)
      ref[RELEVANCE] = [boosted_relevance]
      log_scores.append(math.log(boosted_relevance,10))
    
    #normalization_constant = statistics.median(relevance_vec)
    min_log = min(log_scores)
    max_log = max(log_scores)
    maxmin = max_log - min_log

    table_rows = set()
    for (id_type,identifier), (ref,refcount) in self.ref_counter.items():
      assert(isinstance(ref,Reference))
      if use_relevance:
        score = 100 * (math.log(ref[RELEVANCE][0],10)-min_log)/maxmin
      else:
        score = refcount
        
      biblio_str,_,_ = ref._biblio_tuple()
      if id_type == 'PMID':
          identifier = pubmed_hyperlink([identifier],identifier)
      elif id_type == 'DOI':
          identifier = make_hyperlink(identifier,'http://dx.doi.org/')
      elif id_type == 'PMC':
          identifier = pmc_hyperlink([identifier],identifier)
      elif id_type == 'PII':
          identifier =  pii_hyperlink(identifier,identifier)
      table_rows.add(tuple([score,ref.pubyear(),id_type,identifier,biblio_str]))

    return_df = df.from_rows(list(table_rows),header)
    if use_relevance:
      return_df[RELEVANCE] = return_df[RELEVANCE].astype(float).round(2)
    else:
      return_df[ETM_CITATION_INDEX] = return_df[ETM_CITATION_INDEX].astype(int)

    return_df[PUBYEAR] = pd.to_numeric(return_df[PUBYEAR], errors='coerce')
    return_df = return_df.sortrows(by=first_col)
    return_df.add_column_format('Citation','width',150)
    return_df.add_column_format('Citation','wrap_text',True)
    return_df.make_header_horizontal()
    return_df.set_hyperlink_color([IDENTIFIER_COLUMN])
    return return_df


  @staticmethod
  def external_counter2pd(ref_counter:set[Reference], stat_prop=RELEVANCE):
      '''
      use for external ref_counter = {Reference} where Reference objects are annotated as Reference[stat_prop]
      used to count references in ResnetGraph()
      '''
      table_rows = set()
      for ref in ref_counter:
          biblio_str, id_type, identifier = ref._biblio_tuple()
          if id_type == 'PMID':
              identifier = pubmed_hyperlink([identifier],identifier)
          elif id_type == 'DOI':
              identifier = make_hyperlink(identifier,'http://dx.doi.org/')
          elif id_type == 'NCT ID':
              identifier = make_hyperlink(identifier,'https://clinicaltrials.gov/ct2/show/')
          elif id_type == 'PMC':
              identifier = pmc_hyperlink([identifier],identifier)
          table_rows.add(tuple([ref[stat_prop][0],id_type,identifier,biblio_str]))

      header = [stat_prop,'Identifier type',IDENTIFIER_COLUMN,'Citation']
      return_df = df.from_rows(list(table_rows),header)
      if stat_prop == RELEVANCE:
          return_df[RELEVANCE] = return_df[RELEVANCE].astype(float).round(2)
      else:
          return_df[stat_prop] = return_df[stat_prop].astype(int)
      return_df = return_df.sortrows(by=stat_prop)
      return_df.add_column_format('Citation','width',150)
      return_df.add_column_format('Citation','wrap_text',True)
      return_df.make_header_horizontal()
      return_df.set_hyperlink_color([IDENTIFIER_COLUMN])
      return return_df
  

  def __normalize_journal(self,ref:Reference):
      '''
      updates:
          ref[JOURNAL] by normalized journal name from self.abbrev2journal
          ref['ISSN'] from self.journal2issn.get(norm_journal,[])
      '''
      maybe_abbr = str(ref.get_prop(JOURNAL))
      norm_journal = ref.journal()
      try:
          norm_journal = str(self.abbrev2journal[maybe_abbr])
      except KeyError:
          try:
              norm_journal = str(self.abbrev2journal[norm_journal])
          except KeyError:
              maybe_abbr_no_dots = maybe_abbr.replace('. ',' ')
              try:
                  norm_journal = str(self.abbrev2journal[maybe_abbr_no_dots])
              except KeyError: 
                  pass

      ref[JOURNAL] = [norm_journal]

      my_issns = self.journal2issn.get(norm_journal,[])
      ref.update_with_list('ISSN',my_issns)
      return norm_journal


  def get_publisher(self,ref:Reference):
      '''
      annotates _4ref with fields:
          PUBLISHER, SCOPUS_CITESCORE,SCOPUS_SJR,SCOPUS_SNIP
      '''
      self.__normalize_journal(ref)
      self.Scopus.scopus_stats4(ref)
      return
  

  def citation_index(self,ref:Reference)->Reference:
    '''
    returns citation index for ref from Scopus
    '''
    self.Scopus.scopus_stats4(ref)
    return ref


  @staticmethod
  def count_refs(ref_counter:set, references:list):
      '''
      updates "ref_counter" with "references"
      updates 'Citation index' for reference in "ref_counter" by 1 if reference is in "references"
      '''
      ref_counter.update(references)
      counter_refs = ref_counter.intersection(set(references))
      for ref in counter_refs:
          try:
              count = ref['Citation index'][0]
              ref['Citation index'] = [count + 1]
          except KeyError:
              ref['Citation index'] = [1]
      return
  

  def _2hyperlink(self,references:list[Reference],total_hits:int,default_url='') -> str:
    references.sort(key=lambda x:float(x[RELEVANCE][0]),reverse=True)
    my_limit = self._limit()

    pmids = list()
    dois = list()
    pmcs = list()
    piis = list()
    linkable_id_counter = 0
    for ref in references:
      id_type,identifier = ref.get_doc_id()
      if id_type == 'PMID':
        pmids.append(identifier)
        linkable_id_counter += 1
      elif id_type == 'DOI':
        dois.append(identifier)
      elif id_type == 'PII':
        piis.append(identifier)
      elif id_type == 'PMC':
        pmcs.append(identifier)
        linkable_id_counter += 1
      else:
          continue
      
      if linkable_id_counter == my_limit:
          break

    # remap PMC ID to PMIDs using NCBI web service:
    if pmcs:
      pmids += docid2pmid(pmcs).values() # docid2pmid works only with PMCs or DOIs from PMC!!

    if pmids:
      return pubmed_hyperlink(pmids,total_hits)
    elif piis:
      return pii_hyperlink(piis[0],total_hits)
    elif dois:
      return doi_hyperlink(dois[0],total_hits)
    elif total_hits and default_url and len(default_url) < 255:
      return f'=HYPERLINK("{default_url}","{total_hits}")'
    else:
      return ''


  @staticmethod
  def doi_column(between_column:str, and_concepts:str|list):
      if isinstance(and_concepts,str):
          return 'DOIs' + ' between '+between_column+' and '+and_concepts
      else:
          return 'DOIs' + ' between '+between_column+' and '+','.join(and_concepts)



  @staticmethod
  def refcounter2tsv(fname:str, ref_counter:set[Reference], use_relevance=False,include_idtype=False):
      to_sort = list(ref_counter)
      if use_relevance:
          for ref in to_sort:
              try:
                  relevance_index = float(ref[ETM_CITATION_INDEX][0])*float(ref[RELEVANCE][0])
              except KeyError: relevance_index = 0.0
              ref['Relevance index'] = [float(relevance_index)]
          to_sort.sort(key=lambda x: x['Relevance index'][0], reverse=True)
          with open(fname, 'w', encoding='utf-8') as f:
              f.write('Relevance\tCitation\tPMID or DOI\n')
              for ref in to_sort:
                  biblio_tup = ref._biblio_tuple()
                  if not biblio_tup[2]: 
                      continue
                  if include_idtype:
                      biblio_str = biblio_tup[0]+'\t'+biblio_tup[1]+':'+biblio_tup[2]
                  else:
                      biblio_str = biblio_tup[0]+'\t'+biblio_tup[2]

                  f.write(str(ref['Relevance index'][0])+'\t'+biblio_str+'\n')
      else:        
          to_sort.sort(key=lambda x: x['Citation index'][0], reverse=True)
          with open(fname, 'w', encoding='utf-8') as f:
              f.write('Citation index\tCitation\tPMID or DOI\n')
              for ref in to_sort:
                  biblio_tup = ref._biblio_tuple()
                  if not biblio_tup[2]: 
                      continue #to remove reference with no ID type
                  biblio_tup = ref._biblio_tuple()
                  if include_idtype:
                      biblio_str = biblio_tup[0]+'\t'+biblio_tup[1]+':'+biblio_tup[2]
                  else:
                      biblio_str = biblio_tup[0]+'\t'+biblio_tup[2]
                  f.write(str(ref['Citation index'][0])+'\t'+biblio_str+'\n')


class SBSstats(RefStats):
  def __init__(self, APIconfig, **kwargs):
      super().__init__(APIconfig, **kwargs)
      self.SBSsearch = SBSapi(self.APIconfig)


  def reflinks(self,to_df:df,between_names_in_col:str,and_concepts:list[str],
               add2query:list[str]=[],multithread=False,expand=True)->dict[str,str]:
    '''
    output:
      {name in to_df[between_names_in_col]:sentence_coocurence hyperlinked to top 10 references in pubmed}
    '''
    start_time = time.time()
    names = to_df[between_names_in_col].to_list()
    if expand:
      names = [x+'+' for x in names if ':' not in x]
      self.SBSsearch.multithread = multithread 
      name2refsX = self.SBSsearch.sentcooc4list(names,and_concepts,add2query)
      name2refs = {k[:-1]:v for k,v in name2refsX.items()}
    else:
      names = [x for x in names if ':' not in x] # ':' is invalid character in URL parameter
      name2refs = self.SBSsearch.sentcooc4list(names,and_concepts,add2query,multithread)
    
    names2hyperlinks = dict()
    for name, (search_url,count, refs) in name2refs.items():
      names2hyperlinks[name] = self._2hyperlink(refs,count,search_url)
      [self._add2counter(ref) for ref in refs]
    print(f'Added SBS references to {len(name2refs)} out of {len(to_df)} rows in worksheet "{to_df._name_}" in {execution_time(start_time)}')
    return names2hyperlinks
  

  def add_reflinks(self,names2hyperlinks:dict[str,str], _2col:str,in_df:df, map2col:str)->df:
    '''
    input:
      names2hyperlinks - {name in to_df[map_col]:sentence_coocurence hyperlinked to top 10 references in pubmed}
      _2col - new column name in in_df with hyperlinks
      map_col - column name in in_df used for mapping
    '''
    if not names2hyperlinks: return in_df

    in_df[_2col] = in_df[map2col].map(names2hyperlinks)
    in_df.add_column_format(_2col,'align','center')
    in_df.set_hyperlink_color([_2col])
    self.refcols.add(_2col)
    return in_df
  

  def abs_cooc(self,to_df:df,between_names_in_col:str,and_concepts:list[str],
                   multithread=False)->tuple[dict[str,int],dict[str,float]]:
    '''
    output:
    {name:abs_coocurence},{name:relevance},
      where name is cell value from "between_names_in_col" column
    '''
    names = to_df[between_names_in_col].to_list()
    names = [n for n in names if ':' not in n]
    self.SBSsearch.multithread = multithread
    name2abs_cooc = self.SBSsearch.abscooc4list(names,and_concepts)
    name2relevance = dict()
    name2abscooc = dict()
    for name, (cooc, relevance) in name2abs_cooc.items():
      name2abscooc[name] = cooc
      name2relevance[name] = relevance
    return name2abscooc, name2relevance


  def search_docs(self,query:str):
    refs = self.SBSsearch.search_docs(query)
    pmc_ids = [ref.Identifiers['PMC'] for ref in refs if 'PMC' in ref.Identifiers]
    pmc2pmid = docid2pmid(pmc_ids)
    for r in refs:
      if 'PMC' in r.Identifiers:
        pmc_id = r.Identifiers['PMC'].upper()
        if pmc_id in pmc2pmid:
          r.Identifiers['PMID'] = str(pmc2pmid[pmc_id])
    return refs

##################  DEPRICATED ###################### DEPRICATED #################### DEPRICATED ##############
class ETMStats(RefStats):
  def __init__(self, APIconfig, **kwargs):
    super().__init__(APIconfig, **kwargs)
    self.ETMsearch = ETMsearch(self.APIconfig,**kwargs)

  def __etm42columns(self,in_df:df,between_col:str,and_col:str,my_query,add2query=[]):
      my_df = in_df.dfcopy()
      refcols = self.refcount_column(between_col,and_col)
      doi_ref_column_name = self.doi_column(between_col,and_col)
      for i in in_df.index:
          col1 = my_df.loc[i][between_col]
          col2 = my_df.loc[i][and_col]
          refcountcol,hyperlink2doi = self.__get_refs(col1,[col2],my_query,add2query)
          my_df.at[i,refcols] = refcountcol
          my_df.at[i,doi_ref_column_name] = hyperlink2doi

      return my_df,self.ref_counter


  def refs42columns(self,in_df:df,between_col:str,and_col:str,my_query,add2query=[],max_row=0):
      start_time = time.time()
      if max_row:
          df2annotate = df.from_pd(in_df.iloc[:max_row])
          unannoated_rows = df(in_df.iloc[max_row:])
      else:
          df2annotate = in_df.dfcopy()
          unannoated_rows = df()

      partition_size = 100
      thread_name = f'Retrieve {partition_size} ETM references'
      dfs2concat = list()
      with ThreadPoolExecutor(max_workers=MAX_TM_SESSIONS, thread_name_prefix=thread_name) as e:    
          futures = list()
          for i in range(0,len(df2annotate),partition_size):
              df_part = df(df2annotate.iloc[i:i+partition_size])
              new_session = self.clone() # need to clone here to avoid self.ref_counter mutation
              futures.append(e.submit(new_session.__etm42columns,df_part,between_col,and_col,my_query,add2query))

          for f in futures: # cannot use as_completed here to ensure the same order for concatenation
              session_df, session_refcounter = f.result()
              dfs2concat.append(session_df)
              [self._add2counter(t[0]) for t in session_refcounter.values()] # combine references from all futures 

          dfs2concat.append(unannoated_rows)

      annotated_df = df.concat_df(dfs2concat,in_df._name_)
      
      refcols = self.refcount_column(between_col,and_col)
      doi_column = self.doi_column(between_col,and_col)
      self.refcols.add(refcols)
      self.doi_columns.add(doi_column)

      annotated_df.add_column_format(refcols,'align','center')
      annotated_df.set_hyperlink_color([refcols,doi_column])
      print('Annotated %d rows from %s with ETM references in %s' % 
              (len(in_df),in_df._name_,execution_time(start_time)))
      return annotated_df
  

  def add_refs_etm(self,to_df:df,between_names_in_col:str,and_concepts:list,
               use_query:str,add2query=[],max_row=0):
      """
      input:
          use_query - function to generate query from "between_names_in_col" and each concept in "and_concepts"
          "use_query" must have 3 arguments: my_query(entity1:str, entity2:str, add2query:list)\n
          where add2query - list of additinal keywords used for all pairs "between_names_in_col" and "and_concepts"

      output:
          copy of "to_df" with added columns. Added reference count columns are listed in self.refcols\n
          self.refcols - []\n
          self.doi_columns - []\n
      """
      if to_df.empty: return to_df
      
      if max_row:
          df2annotate = df.from_pd(to_df.iloc[:max_row])
          unannoated_rows = df.from_pd(to_df.iloc[max_row:])
      else:
          df2annotate = to_df.dfcopy()
          unannoated_rows = df()
      row_count = len(to_df)
      
      if row_count > 9:
          thread_name = f'{use_query} refs 4 {len(df2annotate)} rows in {MAX_TM_SESSIONS} thrds'
          partition_size = 100
          with ThreadPoolExecutor(max_workers=MAX_TM_SESSIONS, thread_name_prefix=thread_name) as e:
              futures = list()
              for i in range(0,len(df2annotate),partition_size):
                  df_part = df.from_pd(df2annotate.iloc[i:i+partition_size])
                  new_session = self.clone() # need to clone here to avoid self.ref_counter mutation
                  futures.append(e.submit(new_session.__add_refs1,df_part,between_names_in_col,and_concepts,use_query,add2query))

              dfs2concat = list()
              for f in futures: # cannot use as_completed here to ensure the same order for concatenation
                  session_df, session_refcounter = f.result()
                  dfs2concat.append(session_df)
                  [self._add2counter(t[0]) for t in session_refcounter.values()] # combine references from all futures
              
              dfs2concat.append(unannoated_rows)
              annotated_df = df.concat_df(dfs2concat,to_df._name_)
      else:
          annotated_df = self. __add_refs1(to_df,between_names_in_col,and_concepts,use_query,add2query)[0]

      return annotated_df
  
  def __add_refs1(self,to_df:df,between_names_in_col:str,and_concepts:list,use_query,add2query=[]):
      '''
          applies search results and creates self.refcount_column and self.doi_columns
      '''
      my_df = to_df.dfcopy()
      refcountcol = self.refcount_column(between_names_in_col,and_concepts)
      doi_ref_column_name = self.doi_column(between_names_in_col,and_concepts)
      my_df[[refcountcol,doi_ref_column_name]] = my_df[between_names_in_col].apply(lambda row: self.__get_refs(row,and_concepts,use_query,add2query)).apply(pd.Series)
      return my_df,self.ref_counter # need to return self.ref_counter to concatenate ref_counters after cloning

  def dispatch2TMsoft(self,search_name:str,for_entity:str,concept:str,add2query:list)->tuple[int,dict[str,list],list[Reference]]:
      '''
      input:
          search_name in ['ETMbasicSearch','ETMadvancedSearch','ETMrel','SBSsearch']
      output tuple:
          [0] hit_count - TOTAL number of reference found by ETM basic search 
          [1] ref_ids = {id_type:[identifiers]}, where len(ref_ids) == ETMsearch.params['limit']\n
          id_type is from [PMID, DOI, 'PII', 'PUI', 'EMBASE','NCT ID']\n
          [2] references = [ref] list of Reference objects sorted by ETM relevance. len(references) == ETMsearch.params['limit'] 
          Relevance score is stored in ref['Relevance'] for every reference
      ''' 
      if search_name == 'ETMbasicSearch':
          search4concepts = [for_entity,concept]+add2query
          return self.ETMsearch.basic_search(search4concepts)
      elif search_name == 'ETMadvancedSearch':
          return self.ETMsearch.advanced_query(for_entity,concept,add2query)
      elif search_name == 'ETMrel':
          return self.ETMsearch.advanced_query_rel(for_entity,concept,add2query)
      else:
          return None
  

  def __multiple_search(self,for_entity:str,and_concepts:list,tm_soft:str,
                        add2query=[],getScopusInfo=False)->tuple[list[Reference],int]:
    '''
      performs search for_entity and each concept in "and_concepts"
      input:
        tm_soft - type of TM serach. Options: ETMbasicSearch, ETMadvancedSearch, ETMrel, SBSsearch
    '''
    assert(isinstance(and_concepts,list|set))
    references = set()
    total_hits = 0
    for concept in and_concepts:
      if for_entity != concept:
        hit_count,_,refs = self.dispatch2TMsoft(tm_soft,for_entity,concept,add2query)
        total_hits += hit_count
        references.update(refs)

    [self._add2counter(ref) for ref in references]
    if getScopusInfo:
      [self.Scopus.scopus_stats4(ref) for ref in references]

    return references, total_hits
  

  def __get_refs(self, entity_name:str, concepts2link:list,use_query,add2query=[]):
      references = set()
      references, total_hits = self.__multiple_search(entity_name,concepts2link,use_query,add2query)
      return self._2hyperlink(references,total_hits),''


########## ETMcache class ########## ETMcache class ########## ETMcache class ##########
# ETMcache is used to download and cache ETM results
ETM_CACHE_DIR = os.path.join(os.getcwd(),'ENTELLECT_API/ElsevierAPI/ETM_API/__etmcache__')

class ETMcache (RefStats):
    """
    self.statistics = {prop_name:{prop_value:count}}
    """
    pass
    etm_results_dir = ETM_CACHE_DIR
    etm_stat_dir = ETM_CACHE_DIR
    
    def __init__(self,query,search_name,APIconfig=dict(),etm_dump_dir='',etm_stat_dir='',add_params={},use_scopus=False):
        super().__init__(APIconfig)
        self.ETMsearch.request_type = '/search/advanced?'
        self.ETMsearch.params.update(add_params)
        self.ETMsearch.params.pop('limit')
        self.ETMsearch._set_query(query)
        try:
            self.articles, self.hit_count = self.ETMsearch._get_articles()
            self.can_connect2server = True
        except HTTPError as er:
            print(er)
            print(f'Cannot connect to {APIconfig["ETMURL"]}.  Will use only cache files')
            self.can_connect2server = False
            pass
        self.search_name = search_name
        self.etm_results_dir = etm_dump_dir if etm_dump_dir else ETM_CACHE_DIR
        self.etm_stat_dir = etm_stat_dir if etm_stat_dir else ETM_CACHE_DIR
        self.statistics = dict() # {prop_name:{prop_value:count}}
        self.scopusAPI = AuthorSearch(APIconfig) if use_scopus else None


    def set_query(self,query:str,search_name:str):
        super().ETMsearch._set_query(query)
        self.search_name = search_name


    def set_stat_props(self,stat_props:list):
        for p in stat_props:
            self.statistics[p] = dict()
            if p == AUTHORS: 
                self.statistics[SCOPUS_AUTHORIDS] = dict() 
    

    def __download(self,dumpfile_path:str,cached_articles:list=[],dump_date_modified=''):
        '''
        Input
        -----
        dump_date_create format: 2003-03-19\n
        if dump_date_modified is empty will replace entire cache file if today hit count > len(cached_articles)
        '''
        def __update(articles:list,hit_count:int):
            with ThreadPoolExecutor(max_workers=2, thread_name_prefix='backup ETM cache') as e:
                pathname,_ = os.path.splitext(dumpfile_path)
                try:
                    with zipfile.ZipFile(pathname+'.zip', 'w',zipfile.ZIP_DEFLATED) as z:
                        z.write(dumpfile_path,os.path.basename(dumpfile_path)) 
                except FileNotFoundError:
                    print(f'Cannot backup to {pathname}.zip')
                    pass

            start = time.time()
            for page_start in range(self.ETMsearch.page_size,hit_count,self.ETMsearch.page_size):
                more_articles,_ = self.ETMsearch._get_articles(page_start=page_start)
                articles += more_articles
                print(f"Downloaded {len(articles)} out of {hit_count} hits in {execution_time(start)}")
            
            e.shutdown()

        if dump_date_modified: # update cache with articles published after dump_date_create
            if datetime.strptime(dump_date_modified, "%Y-%m-%d").date() < date.today(): 
                self.ETMsearch.params.update({'so_d':dump_date_modified+'<'}) # dump_date_create format: 2003-03-19
                update_articles, update_hit_count = self.ETMsearch._get_articles()
                print(f'Updating cache modified on {dump_date_modified} with {update_hit_count} recent articles')
                __update(update_articles,update_hit_count)
                self.ETMsearch.params.pop('so_d','')
                self.articles = cached_articles+update_articles
                return True
            else:
                self.articles = cached_articles
                return False
        elif self.hit_count > len(cached_articles): #replace or create new cache.  
            print('Today ETM search finds %d results. %d more than in cache' %(self.hit_count, self.hit_count-len(cached_articles)))
            __update(self.articles,self.hit_count)
            return True
        else:
            print('No new articles found. Will use articles from cache')
            return False


    def __dumps(self, articles:list, into_file:str):
        my_file = into_file if into_file else self.search_name+'.json'
        print(f'Writing ETM data into {my_file} file')
        dump = open(my_file, "w", encoding='utf-8')
        dump.write(json.dumps(articles,indent=1))
 

    def load_from_json(self,use_cache=True):
        '''
        Loads articles from ETM cache if use_cache = True\n
        Downloads additional articles from ETM if necessary\n
        articles in self.artciles are annotated with RELEVANCE equal to ETM score
        '''
        dump_dir = os.path.join(self.etm_results_dir,'')    
        dumpfile_path = dump_dir+self.search_name+'.json'
        if use_cache:  
            try:
                # attepts to find json file with ETM results saved after ETM API call below
                cached_articles = json.load(open(dumpfile_path,'r'))
                fname = os.path.basename(dumpfile_path)
                print(f'Found file "{fname}" in "{dump_dir}" directory with {len(cached_articles)} articles.\nWill use it to load cached results')
                if self.can_connect2server:
                    d =  datetime.strptime(time.ctime(os.path.getmtime(dumpfile_path)), "%c")
                    dump_modified_date =  d.strftime('%Y-%m-%d')
                    if self.__download(dumpfile_path,cached_articles,dump_modified_date):
                        self.__dumps(self.articles,dumpfile_path)
                else:
                    self.articles = cached_articles
            except (FileNotFoundError,json.JSONDecodeError):
                #if json dump file is not found or it is corrupted new ETM search is initiated
                self.__download(dumpfile_path)
                self.__dumps(self.articles,dumpfile_path)
        else: # if not use cache
                self.__download(dumpfile_path)
                self.__dumps(self.articles,dumpfile_path)
        

    def scopus_annotation(self):
        def set_oa_status(ref:Reference):
            ref_doi = ref.doi()
            if ref_doi:
                try:
                    ref[IN_OPENACCESS] = [bool(doi2oa[ref_doi])]
                except KeyError:
                    is_in_oa = self.Scopus.is_in_open_access(ref)
                    ref[IN_OPENACCESS] = [is_in_oa]
                    doi2oa[ref_doi] = is_in_oa
            else:
                print(f'#{i} article "{etm_ref.title()}" out of {len(self.articles)} has no DOI!!!')

        dump_dir = os.path.join(self.etm_results_dir,'')
        do12oa_dump = dump_dir+self.search_name+'_do12oa.json'
        try:
            doi2oa = json.load(open(do12oa_dump, 'r',encoding='utf-8'))
        except FileNotFoundError:
            doi2oa = dict()

        print(f'Reannotating articles from "{self.search_name}" query with Scopus data')
        for i,article in enumerate(self.articles):
            etm_ref = ETMsearch.article2ref(article)
            if etm_ref: # etm_ref can be empty if it is conference proceedings
            # scopus_authors = self.AuthorSearch.get_authors(etm_ref)
                if etm_ref.journal() != GRANT_APPLICATION:
                    with ThreadPoolExecutor(max_workers=3, thread_name_prefix='Scopus annotation') as e:
                        #e.submit(self.AuthorSearch.normalize_institution,etm_ref)
                        e.submit(set_oa_status,etm_ref)
                        e.submit(self.get_publisher,etm_ref)
                    e.shutdown()

                relevance_score = float(article['score'])
                etm_ref[RELEVANCE] = [relevance_score]
                self._add2counter(etm_ref)

        self.AuthorSearch.close()
        self.Scopus.close()
        with open(do12oa_dump, 'w',encoding='utf-8') as f:
            json.dump(doi2oa,f,indent=1)


    def get_statistics(self, stat_prop_list):
        self.term2refs = dict() # {term:{ref}}
        self.keywords2ref = dict() #{keyword:{ref}}
        self.scopusid2name = dict()
        references = self.references()
        for ref in references:
            [ref.count_property(self.statistics[p], p) for p in stat_prop_list]

            if hasattr(ref,'term_ids'):
                for term in ref.term_ids: #term_ids = {term_id+'\t'+term_name}
                    try:
                        self.term2refs[term].add(ref)
                    except KeyError:
                        self.term2refs[term] = {ref}
            
            if hasattr(ref,'keywords'):
                for k in ref.keywords:
                    try:
                        self.term2refs[k].add(ref)
                    except KeyError:
                        self.term2refs[k] = {ref}


    def _org2address_stats(self):
        return dict(sorted(self.statistics[INSTITUTIONS].items(), key=lambda item: item[1],reverse=True))


    def count_oa_stats(self):
        journal_oa_counts = defaultdict(int)
        for ref in self.references():
            if ref.get_prop(IN_OPENACCESS,0,''):
                journal_oa_counts[ref.journal()] += 1
        return dict(journal_oa_counts)


    def count_affiliations(self,affils:set):
        journal_aff_counts = dict()
        normalized_affils = set(self.AuthorSearch.normalize_affiliations(list(affils)))
        for ref in self.references():
            ref_affils = set(ref.get(INSTITUTIONS,[]))
            if not ref_affils.isdisjoint(normalized_affils):
                journal_aff_counts[ref.journal()] += 1
        return journal_aff_counts


    def to_excel(self, stat_props:list,_4affiliations:set={}):
      if self.scopusid2name:
        replaceid4name = dict()
        scopus_stats = dict(self.statistics[SCOPUS_AUTHORIDS])
        for k,v in scopus_stats.items():
          au_name = self.scopusid2name[k]
          replaceid4name[au_name] = v
        self.statistics[SCOPUS_AUTHORIDS] = replaceid4name
      
      excel_file = os.path.join(self.etm_stat_dir,self.search_name+'_ETMstats.xlsx')
      writer = pd.ExcelWriter(excel_file, engine='xlsxwriter')
      print(f'Writing ETM statistics into {excel_file} file')
      for p in stat_props:
          try:
              dic = dict(self.statistics[p])
              by_key = True if p == PUBYEAR else False
              #prop_stat = self._org2address_stats() if p == INSTITUTIONS else self._sort_dict(dic,sort_by_value)
              prop_stat = sortdict(dic,by_key)
              stat_df = df.from_dict2(prop_stat,p,'#References')
              if p == JOURNAL:
                  oa_stats = self.count_oa_stats()
                  stat_df = stat_df.merge_dict(oa_stats,'Articles in open access',JOURNAL)

                  if _4affiliations:
                      aff_stats = self.count_affiliations(_4affiliations)
                      stat_df = stat_df.merge_dict(aff_stats,'Affiliation count',JOURNAL)
                  
                  citescore_dict = {v[0]:v[2] for k,v in self.AuthorSearch.JournalInfo.items()}
                  stat_df = stat_df.merge_dict(citescore_dict,SCOPUS_CITESCORE,JOURNAL)
                  sjr_dict = {v[0]:v[3] for k,v in self.AuthorSearch.JournalInfo.items()}
                  stat_df = stat_df.merge_dict(sjr_dict,SCOPUS_SJR,JOURNAL)
                  csnip_dict = {v[0]:v[4] for k,v in self.AuthorSearch.JournalInfo.items()}
                  stat_df = stat_df.merge_dict(csnip_dict,SCOPUS_SNIP,JOURNAL)
                  publ_dict = {v[0]:v[1] for k,v in self.AuthorSearch.JournalInfo.items()}
                  stat_df = stat_df.merge_dict(publ_dict,PUBLISHER,JOURNAL)
              stat_df.df2excel(writer,p)
          except KeyError: continue

      ref_df = self.counter2df()
      ref_df.df2excel(writer,'References')
      writer.close()

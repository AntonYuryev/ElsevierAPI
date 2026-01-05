import psycopg2
from ..utils import load_api_config, plot_distribution,print_error_info
from ..pandas.panda_tricks import df,pd
from concurrent.futures import ThreadPoolExecutor,as_completed
from ..ETM_API.references import AUTHORS,JOURNAL,MEDLINETA,SENTENCE,PUBYEAR,TITLE,Reference
from collections import defaultdict

SNIPPET_ID = 'unique_id'
RELATION_ID = 'id'
REFID2ATTR = {'doi':'DOI','pmid':'PMID','embase':'EMBASE','pii':'PII', 'pui':'PUI','nct_id':'NCT ID'}

SENTENCE_PROPS = {'msrc':SENTENCE,'organism':'Organism','source':'Source','textmods':'TextMods','organ':'Organ','tissue':'Tissue',
  'biomarkertype':'BiomarkerType', 'celllinename':'CellLineName','celltype':'CellType', 'px':'pX','quantitativetype':'QuantitativeType',
  'changetype':'ChangeType','collaborator':'Collaborator','company':'Company','condition':'Condition',
  'experimental_system':'Experimental System','intervention':'Intervention','percent':'Percent',
  'phase':'Phase','start':'Start', 'studytype':'StudyType','trialstatus':'TrialStatus','url':'URL'}

COLUMN2ATTR = {'authors':AUTHORS,'title':TITLE,'pubyear':PUBYEAR,'journal':JOURNAL,'medlineta':MEDLINETA,
               'issn':'ISSN','essn':'ESSN','id':'Postgres ID'}

SCOPUS_DATA = {'citation_type':'Article type', 'citation_count':'Citations', 'fwci':'FWCI', 
'fwci_perc':f'FWCI %ile', 'citation_count_ns':'Non-self citations', 'fwci_ns':'Non-self FWCI', 
'fwci_perc_ns':f'Non-self FWCI %ile', 'citescore2024':'CiteScore2024', 
'min_asjc_citescore_percentile_raw':f'CiteScore %ile', 'patent_citation_count':'Patent citations', 
'corporate':'Corporate', 'num_refs':'References', 'independent_ref_count':'Independent References', 
'document_score':'Document score', 'relation_score':'Relation score'}

class PostgreSQL:
  def __init__(self, APIconfig:dict={}, resnetV='resnet18'):
    self.resnet_version = resnetV
    self.rel2refDict = dict() # {int(relid):[Reference]}
    self.executor = ThreadPoolExecutor(thread_name_prefix='postgres')
    self.futures = [] # futures of reference retrieval
    if not APIconfig:
      APIconfig = load_api_config()

    try:
        self.db = psycopg2.connect(
            database = APIconfig['postgreSQLdb'],
            host = APIconfig['postgreSQLhost'],
            user = APIconfig['postgreSQLuser'],
            password = APIconfig['postgreSQLpswd'],
            port = APIconfig['postgreSQLport']
          )
        print('Connected to Postgres')
    except Exception as e:
        self.db = None
        print("Error connecting to PostgreSQL:", e)

  def close(self):
    """
    Closes the connection to PostgreSQL.
    """
    if self.db:
      self.db.close()
      print('Postgres connection closed')


  def full_table_name(self,table:str):
    return f'{self.resnet_version}.{table}'


  def get_relprops(self,relations_id:list[str]):
    relid_str = ','.join(map(str, relations_id))
    sql = f"SELECT * FROM {self.resnet_version}.control WHERE {self.resnet_version}.control.id IN ({relid_str})"

    with self.db.cursor() as cur:
        cur.execute(sql)
        data = cur.fetchall()
        colnames = [desc[0] for desc in cur.description]

    return df(data, columns=colnames)
  

  def get_stats(self,table:str, columns:list[str],filter:dict[str:tuple[str,int]]=dict()):
    '''
    input:
      filter = {columns_name:(sign,value)}
    '''
    count_str = ','.join([f'count({c})' for c in columns])
    count_str += ',count(*)'
    sql = f'SELECT {count_str} from {self.resnet_version}.{table}'
    if filter:    
      clauses = [f'{colname} {value[0]} {value[1]}' for colname,value in filter.items()]
      sql += ' WHERE '+ ' AND '.join(clauses)
    with self.db.cursor() as cur:
      cur.execute(sql)
      data = cur.fetchall()
      counts = list(data[0])
      row_count = counts[-1]
      result_str = f'"{table}" table has {row_count} rows\n'
      for i,count in enumerate(counts[:-1]):
        result_str += f'{columns[i]}: {count}\n'
      print('Counts:\n'+result_str)
      return
    

  def averages(self,table:str, columns:list[str],min_values:dict[str:tuple[str,int]]=dict()):
    '''
    input:
      filter = {columns_name:(sign,value)}
    '''
    count_str = ','.join([f'AVG({c})' for c in columns])
    sql = f'SELECT {count_str} from {self.resnet_version}.{table}'
    if min_values:    
      clauses = [f'{colname} {value[0]} {value[1]}' for colname,value in min_values.items()]
      sql += ' WHERE '+ ' AND '.join(clauses)
    with self.db.cursor() as cur:
      cur.execute(sql)
      data = cur.fetchall()
      averages = list(data[0])
      result_str = ''
      for i,average in enumerate(averages):
        result_str += f'Average of {columns[i]}: {average:.2f}\n'
      print(result_str)
      return averages
    
  
  def execute_sql(self, sql:str):
    with self.db.cursor() as cur:
      cur.execute(sql)
      data = cur.fetchall()
      return list(data)


  def plot_distribution(self,table:str,columns:list[str],outdir=''):
    for col in columns:
      sql = f'SELECT {col} FROM {self.resnet_version}.{table}'
      with self.db.cursor() as cur:
        cur.execute(sql)
        distribution = cur.fetchall()
        distribution = [t[0] for t in distribution]
        print(f'Will plot distribution for {table}.{col} with {len(distribution)} values')
        plot_distribution([{f'{table}.{col}':distribution}],outdir=outdir)


  def get_refs(self,relations_id:list[str]):
    new_relids = set(relations_id).difference(self.rel2refDict)
    newrelid_str = ','.join(map(str, new_relids))
    sql = f"SELECT * FROM {self.resnet_version}.reference WHERE {self.resnet_version}.reference.id IN ({newrelid_str})"
    with self.db.cursor() as cur:
      cur.execute(sql)
      rows = cur.fetchall()
      colnames = [desc[0] for desc in cur.description]
    rows = [list(r) for r in rows]
    ref_pd = pd.DataFrame(rows,columns=colnames)
    return ref_pd
    

  def submit_refs(self, relations_ids:list[str]):
    """
      submits reference fetching job to ThreadPoolExecutor future that is added to self.futures
    """
    if relations_ids:
      self.futures.append(self.executor.submit(self.get_refs,relations_ids))


  def scopus_data(self, refids:list[int]):
    '''
      scopus_data.reference_id is joined with reference.unique_id
    '''
    refid_str = ','.join(map(str, refids))
    sql = f"SELECT * FROM {self.resnet_version}.scopus_data WHERE {self.resnet_version}.scopus_data.reference_id IN ({refid_str})"
    with self.db.cursor() as cur:
      cur.execute(sql)
      rows = cur.fetchall()
      colnames = [desc[0] for desc in cur.description]
    rows = [list(r) for r in rows]
    scopus_pd = pd.DataFrame(rows,columns=colnames).set_index('reference_id')
    return scopus_pd
  

  def __rows2refs(self,ref_pd:pd.DataFrame,scopus_pd:pd.DataFrame)->dict[str,list[Reference]]:
    '''
    output:
      {relation_id:[Reference]}
      scopus_data.reference_id is joined with reference.unique_id
    '''
    relid2refs = defaultdict(list)
    for refpd_idx in ref_pd.index:
      relid = ref_pd.at[refpd_idx,RELATION_ID]
      textref = ref_pd.at[refpd_idx,'textref']
      ref = dict()
      ref_idtypes = list(REFID2ATTR.keys())
      for idtype_idx, idtype in enumerate(ref_idtypes):
        refid = ref_pd.at[refpd_idx,idtype]
        if not pd.isna(refid):
          ref = Reference(REFID2ATTR[idtype],refid)
          for idt in ref_idtypes[idtype_idx+1:]:
            id = ref_pd.at[refpd_idx,idt]
            if not pd.isna(id):
              ref.Identifiers[REFID2ATTR[idt]] = id
          break

      if isinstance(ref,Reference):
        for col, attr in SENTENCE_PROPS.items():
          attr_val = ref_pd.at[refpd_idx,col]
          if not pd.isna(attr_val):
            ref.add_sentence_prop(textref,attr,attr_val)

        for col, attr in COLUMN2ATTR.items():
          attr_val = ref_pd.at[refpd_idx,col]
          if not pd.isna(attr_val):
            if attr in [PUBYEAR]:
              attr_val = int(attr_val)
            ref[attr] = [attr_val]

        snippet_id = int(ref_pd.at[refpd_idx,SNIPPET_ID])
        if snippet_id in scopus_pd.index:
          ref_scopus_data = scopus_pd.loc[snippet_id].to_dict()
          ref.update({SCOPUS_DATA[k]:[v] for k,v in ref_scopus_data.items() if k in SCOPUS_DATA})
        else:
          refid = ref.doi_or_id()
          if not refid.startswith('NCT'):
            print(f'Reference {refid} has no Scopus data')
      
        ref.toAuthors()
        relid2refs[int(relid)].append(ref)

    return dict(relid2refs)


  def load_refs(self):
    """
    output:
      {embio_relation_id:[Reference]}
    """
    if self.futures:
      processed_futures = []
      print(f'Got {len(self.futures)} Postgres futures to process')
      for get_refs_future in as_completed(self.futures):
        try:
          ref_pd = get_refs_future.result()
          scopus_pd = self.scopus_data(set(ref_pd[SNIPPET_ID].to_list()))
          self.rel2refDict.update(self.__rows2refs(ref_pd,scopus_pd))
          processed_futures.append(get_refs_future)
        except Exception as e:
          print_error_info(e,f'Error loading references from Postgres with SQL {e.cursor.query.decode()}')
          self.db.rollback()
      self.futures = [f for f in  self.futures if f not in processed_futures]
      #print(f'Processed {len(processed_futures)} out of {self.futures} postgres futures')
      print(f'Cached references for {len(self.rel2refDict)} relations from Postgres')
      #print(f'{len(self.futures)} futures remains in cue')
    return self.rel2refDict
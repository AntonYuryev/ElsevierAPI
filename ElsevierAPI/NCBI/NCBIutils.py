
#C:Windows> py -m pip install entrezpy --user
import json, datetime,os
import xml.etree.ElementTree as ET
from time import sleep
from urllib.parse import urlencode
from collections import defaultdict
from titlecase import titlecase
from ..utils import attempt_request4,remove_duplicates,sortdict
from ..ETM_API.references import Reference,JOURNAL,TITLE,AUTHORS,PUBYEAR


RETMAX = 10000

class NCBIeutils:
  baseURL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
  
  def __init__(self,db:str,cache_path:str,retmode='xml'):
      #database names are in # from https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly
      self.params = {'db':db,'retmode':retmode}
      self.cache_path = cache_path
      self.query = ''
  
  def _esearch_url(self,params:dict):
      return self.baseURL+'esearch.fcgi?'+urlencode(params)
  
  def _efetch_url(self,params:dict):
      return self.baseURL+'efetch.fcgi?'+urlencode(params)
  
  def mydb(self):
    return self.params['db']
  
  def myquery(self):
    return self.query
  

  def get_count(self, query:str=''):
    q = query if query else self.myquery()
    self.params = {'db':self.mydb(),'term':q}
    count_param = dict(self.params)
    count_param.update({'rettype':'count','retmode':'json'})
    my_url = self._esearch_url(count_param)
    http_response = attempt_request4(url=my_url)
    data = dict(json.loads(http_response.data.decode('utf-8')))
    return int(data["esearchresult"]['count'])
  

  def _retmax_uids(self,params:dict={})->list[int]:
    '''
    Return
    ------
    [PMIDs] with size < RETMAX sorted in ascending order
    '''
    my_params = {'db':self.mydb(),'term':self.myquery(),'retstart':0,'retmax':RETMAX,'sort':'pub_date','retmode':'json'}
    my_params.update(params)
    my_url = self._esearch_url(my_params)
    http_response = attempt_request4(my_url)
    data = dict(json.loads(http_response.data.decode('utf-8')))
    ids = list(map(int,data["esearchresult"]['idlist'])) # NCBI returns PMIDs list in descending year order
    reversed_list = []
    [reversed_list.append(ids[i]) for i in range(len(ids) - 1, -1, -1)]
    return reversed_list
  

  def get_uids(self,query_name:str)->list[int]:
      '''
      Return
      ------
      List of PMIDs sorted by PDAT
      '''
      json_id_dump = query_name+'PMIDs.json'
      json_id_path = os.path.join(self.cache_path,json_id_dump)
      try:
        all_ids = json.load(open(json_id_path,'r'))
        print(f'Loaded {len(all_ids)} IDs from {json_id_dump}')
        return all_ids
      except FileNotFoundError:
        count = self.get_count()
        all_ids = list()
        if count > RETMAX: # spliting downlaods by year
          current_year = datetime.date.today().year
          for year in range(1965, current_year,1):
            year_params = dict(self.params)
            year_params['term'] += f' AND ({year}/01/01[PDat]:{year}/12/31[PDat])'
            all_ids += self._retmax_uids(year_params)
            sleep(0.5) # 0.1 - too many requests
          all_ids = remove_duplicates(all_ids)
          json.dump(all_ids, open(json_id_path,'w'), indent = 2)
          return all_ids
        else:
          return self._retmax_uids(self.params)


  def ids2records(self,ids:list[int]):
    params = {'db':self.params['db'],'rettype':'XML'}
    stepSize = 200
    for i in range(0, len(ids), stepSize):
      str_ids = ','.join(str(s) for s in ids[i:i+stepSize])
      params.update({'id':str_ids})
      my_url = self._efetch_url(params)
      http_response = attempt_request4(my_url)
      yield http_response.data.decode()


  def fetch(self,query_name:str):
    allids = self.get_uids(query_name)
    return self.ids2records(allids)


  def path2cache(self,query_name:str,extension='xml'):
    return self.cache_path+query_name+"."+extension
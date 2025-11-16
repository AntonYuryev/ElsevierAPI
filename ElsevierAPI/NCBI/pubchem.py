
#C:Windows> py -m pip install entrezpy --user
import urllib.request, json, os
import xml.etree.ElementTree as ET
from .NCBIutils import NCBIeutils

NCBI_CACHE = os.path.join(os.getcwd(),'ENTELLECT_API/ElsevierAPI/NCBI/__ncbipubchemcache__/')

class Pubchem(NCBIeutils):
  def __init__(self):
    self.params = {'db':'pccompound'} # pcassay,pccompound
    # from https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly
    self.cache_path = NCBI_CACHE


  def get_uids(self,query:tuple[str,str]):
      '''
      Return
      ------
      List of PMIDs sorted by PDAT
      '''
      query_name = query[1]
      json_id_dump = query_name+'.json'
      json_id_path = self.cache_path+json_id_dump
      try:
          all_ids = json.load(open(json_id_path,'r'))
          print(f'Loaded {len(all_ids)} IDs from {json_id_dump}')
          return all_ids
      except FileNotFoundError:
          count = self.get_count(query)
          all_ids = list()
          return self._retmax_uids(self.params)
                    

  def fetch(self,query:tuple[str,str]):
    allids = self.get_uids(query)
    params = {'db':self.params['db'],'rettype':'XML'}
    stepSize = 200
    for i in range(0, len(allids), stepSize):
      ids = ','.join(str(s) for s in allids[i:i+stepSize])
      params.update({'id':ids})
      my_url = self._efetch_url(params)
      req = urllib.request.Request(url=my_url)
      xml_str = urllib.request.urlopen(req).read()
      yield xml_str


  def download_pubmed(self,query:tuple[str,str]):
      """
      Output
      ------
      query_name.json: has list of pmids found by query\n
      query_name.xml: pubmed abstracts
      """
      query_name = query[1]
      fpath = self.path2cache(query_name)
      result_counter = 0
      with open(fpath, "w", encoding='utf-8') as file_result:
        file_result.write('<PubmedArticleSet>\n')
        for xml_str in self.fetch(query):
          abstractTree = ET.fromstring(xml_str)
          articles = abstractTree.findall('PubmedArticle')
          result_counter += len(articles)
          for article in articles:
            file_result.write(ET.tostring(article, encoding="unicode"))
#          sleep(1.0)
        file_result.write('</PubmedArticleSet>\n')
        print(f'Downloaded {result_counter} pubmed abstracts')
        print(f'Downloaded abstracts are in "{fpath}"')

        fstatpath = os.path.join(self.cache_path,query_name+'_stats.tsv')
        with open(fstatpath, "w", encoding='utf-8') as f:
            [f.write(f'{k}\t{v}\n') for k,v in self.journalCounter.items()]
        print(f'Statistics is in "{fstatpath}"')


if __name__ == "__main__":
  pubchem = Pubchem('Selegiline')
  pubchem.download_pubmed('Selegiline')


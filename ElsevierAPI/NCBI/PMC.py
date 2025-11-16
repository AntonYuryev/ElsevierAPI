
#C:Windows> py -m pip install entrezpy --user
import json,os
import xml.etree.ElementTree as ET
from ..utils import attempt_request4,list2chunks_generator
from .NCBIutils import NCBIeutils
from ..pandas.panda_tricks import df


RETMAX = 10000
PMC_CACHE = os.path.join(os.getcwd(),'ENTELLECT_API/ElsevierAPI/NCBI/__ncbipmccache__/')


class PMC(NCBIeutils): 
  def __init__(self,query:str,retmode='xml'):
    super().__init__('pmc',PMC_CACHE,retmode)
    self.query = query
    self.params.update({'filter':'availability.pmc_public'})

  def download(self,query_name:str):
    fpath = self.path2cache(query_name)
    counter = 0
    with open(fpath, "w", encoding='utf-8') as result:
      result.write('<pmc-articleset>\n')
      for xml_str in self.fetch(query_name):
        articles = ET.fromstring(xml_str).findall('article')
        [result.write(ET.tostring(a, encoding="unicode")) for a in articles]
        counter += len(articles)
      result.write('</pmc-articleset>\n')
    print(f'Downloaded {counter} PMC articles')


  @staticmethod
  def pub_ids(article:ET.Element)->list[str]:
    pubids = dict()
    for id_elem in article.findall('front/article-meta/article-id'):
      pubid_type = id_elem.attrib.get('pub-id-type','')
      if pubid_type:
        pubids[pubid_type] = id_elem.text
    return pubids


  def download_by_ids(self,docids:list[str],query_name:str):
    id_fname = self.path2cache(query_name,'tab')
    xml_fname = self.path2cache(query_name,'xml')
    downloaded_ids = set()
    pubids = list()
    with open(xml_fname,"w",encoding='utf-8') as xmlfile, open(id_fname,'w') as idfile:
      xmlfile.write('<pmc-articleset>\n')
      for i, id_chunk in list2chunks_generator(docids,chunk_size=75):
        id_query = ' OR '.join([f'"{docid}"' for docid in id_chunk])
        id_query  += ' [lid]'
        self.query = id_query
        for xml_str in self.fetch(query_name):
          articles = ET.fromstring(xml_str).findall('article')
          for a in articles:
            a_ids = self.pub_ids(a)
            downloaded_ids.update([v for idtype,v in a_ids.items() if idtype in ['pmcid','doi','pii']])
            pubids.append(a_ids)
            xmlfile.write(ET.tostring(a, encoding="unicode"))
      xmlfile.write('</pmc-articleset>\n')
    id_df = df(pubids)
    id_df.to_csv(id_fname,sep='\t',index=False)
    missed_ids = set(docids) - downloaded_ids
    if missed_ids:
      missedids_fname = self.path2cache(query_name,'missed_ids.txt')
      with open(missedids_fname,'w',encoding='utf-8') as mf:
        [mf.write(f'{mid}\n') for mid in missed_ids]

    print(f'Downloaded {len(pubids)} PMC articles into {xml_fname}\nCould not find {len(missed_ids)} IDs')
    print(f'IDs for downloaded articles are in {id_fname}\nMissed IDs are in {missedids_fname}')


  def download_by_idfile(self,idfile:str,query_name:str):
    with open(idfile,'r',encoding='utf-8') as f:
      docids = [line.strip() for line in f if line.strip()]
    self.download_by_ids(docids,query_name)


def docid2pmid(docids:list[str])->dict[str,str]:
  '''
  input:
    docids - homogeneous list of PMC or DOI 
    only DOIs for PMC articles will produce valid response
  '''
  docid2pmid = dict()
  for chunk in range(0,len(docids),100):
    idstr = ','.join(docids[chunk:chunk+100])
    requesturl = 'https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids='+idstr+'&format=json'
    http_response = attempt_request4(requesturl)
    if http_response:
      data = json.loads(http_response.data.decode('utf-8'))
    if data:
      for record in data["records"]:
        if 'pmcid' in record: 
          if 'pmid' in record:
            pmid = record['pmid']
            pmcid = record['pmcid']
            doi = record.get('doi','')
            doc_id = doi if doi in docids else pmcid
            docid2pmid[doc_id] = pmid
          else:
            print(f'No PMID in {record}')
  return docid2pmid


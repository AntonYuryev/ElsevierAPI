
#C:Windows> py -m pip install entrezpy --user
import urllib.request, urllib.parse, json, datetime,os
import xml.etree.ElementTree as ET
from time import sleep
from collections import defaultdict
from titlecase import titlecase

PUBMED_URL = 'https://pubmed.ncbi.nlm.nih.gov/?'
PMC_URL = 'https://www.ncbi.nlm.nih.gov/pmc/?'
DOI_URL = 'http://dx.doi.org/'
RETMAX = 10000
NCBI_CACHE = os.path.join(os.getcwd(),'ENTELLECT_API/ElsevierAPI/NCBI/__ncbipubmedcache__/')
#query = ['Aquilegia OR Arabidopsis OR Hordeum OR Brassica OR Theobroma OR Phaseolus OR Gossypium OR "Vitis vinifera" OR Mesembryanthemum OR Lactuca OR "Zea mays" OR Medicago OR "Nicotiana benthamiana" OR Allium OR Capsicum OR "Petunia hybrida" OR Pinus OR Strobus OR Solanum OR Oryza OR Secale OR "Sorghum bicolor" OR Picea OR Saccharum OR Helianthus OR Festuca OR Lycopersicon OR Triphysaria OR Triticum OR Populus OR anthocyanin OR "abscisic acid" OR brassinolide OR brassinosteroid OR brassinosteroids OR cytokinin OR cytokinins OR gibberellin OR gibberellins OR "gibberellic acid" OR "jasmonic acid" OR jasmonate OR jasmonates OR cultivar OR cultivars OR endosperm AND ((2021/08/07[PDat]:3019/12/31[PDat])','ArabidopsisUpdate']
#query = ['(rice OR oryza OR sativa) AND (2021/08/07[PDat]:3019/12/31[PDat])','RiceUpdate']
#query= ['(corn OR maize OR zea OR mays) AND (2021/08/07[PDat]:3019/12/31[PDat])','CornUpdate']
#query= '(COVID19 OR SARS OR coronavirus) AND ("2020/09/17"[Date - Create] : "3000"[Date - Create])'
#query= "cytolethal distending toxin"
#query = '"inflammaging" OR "Inflamm-aging" OR "inflamm-ageing" OR "inflammageing" OR "age-increased inflammation" OR "age-associated inflammation" OR "age increased inflammation" OR "age associated inflammation"'
#query= '("lewy body" OR "lewy bodies" OR protofibrils OR protofibril)'
#query='(infection OR infectious OR pathogen OR parasite) AND (1989 [dp] OR 1994 [dp] OR 1999 [dp] OR 2004 [dp] OR 2009 [dp] OR 2014 [dp] OR 2019 [dp]) AND aug NOT (virus OR viral)'
#query = 'MAIT OR "Mucosa-associated invariant T" OR "Mucosa-associated invariant T cells" OR "Mucosa-associated invariant T cell"'
#query = 'bacillus OR subtilis AND ("2021/10/21"[Date - Create] : "3000"[Date - Create])'
#qname = 'Bsubtilis_since10212021'
#query = 'Clostridium AND ("2021/10/21"[Date - Create] : "3000"[Date - Create])'
#qname = 'Clostridium_since10212021'
#query = 'Cupriavidus OR necator AND ("2021/10/21"[Date - Create] : "3000"[Date - Create])'
#qname = 'Cnecator_since10212021'
#query = 'escherichia OR coli AND ("2021/10/21"[Date - Create] : "3000"[Date - Create])'
#qname = 'Ecoli_since10212021'
#query = 'H3K27M OR (K27M AND histone)'
#query = '(mutant OR mutated OR mutation OR mutations) AND (p53 OR tp53)'
#query = '("2019/11/29"[Date - Entry] : "3000"[Date - Entry])'
#query = '"Perineuronal nets" OR "Perineuronal net" OR PNNs'
#query = 'paxalisib OR GDC0084 OR "GDC-0084" OR P5DKZ70636 OR "EX-A1019"'
#query = '(sulfoquinovosyl OR monogalactosylmonoacylglycerol OR digalactosyldiacylglycerol OR monogalactosyldiacylglycerol OR digalactosylmonoacylglycerol OR monogalactosylmonoacylglycerol OR digalactosyldiacylglycerols OR monogalactosyldiacylglycerols OR digalactosylmonoacylglycerols) AND (human OR sapiens OR mammal OR rat OR mouse)'
#query = 'serine OR threonine OR tyrosine OR histidine'
#query = ['"retracted article" OR  "retracted publication"','retracted articles']
#query = ["Altria[Affil] OR Cronos Group[Affil] OR Lexaria Bioscience[Affil] OR Micreos[Affil] OR Philip Morris International[Affil] OR Softhale[Affil] OR Biovotion[Affil] OR Fertin Pharma[Affil] OR Swedish Match[Affil] OR OtiTopic[Affil] OR Biognosys[Affil] OR Syqe Medical[Affil] OR Vectura Group[Affil] OR Biofourmis[Affil] OR Japan Tobacco[Affil] OR Tularik[Affil] OR Akros Pharma[Affil] OR Torii Pharmaceutical[Affil] OR Imperial Brands[Affil] OR Auxly Cannabis[Affil] OR Compania de Distribucion[Affil] OR Oxford Cannabinoid Technologies[Affil] OR British American Tobacco[Affil] OR Philter Labs[Affil] OR Open Book Extracts[Affil] OR Velo[Affil] OR The Kanvas Co[Affil] OR Cannopeia[Affil] OR Purissima[Affil] OR Hesperos[Affil] OR Organigram[Affil] OR Trait Biosciences[Affil] OR Ajna Biosciences[Affil] OR Niconovum[Affil] OR Sanity Group[Affil] OR PlatoScience[Affil]' OR Charlotte's Web[Affil] OR KBio[Affil]",'TobaccoCompaniesPubs']
query = ["(Altria[Affil] OR Cronos Group[Affil] OR Lexaria Bioscience[Affil] OR Micreos[Affil] OR Philip Morris International[Affil] OR Softhale[Affil] OR Biovotion[Affil] OR Fertin Pharma[Affil] OR Swedish Match[Affil] OR OtiTopic[Affil] OR Biognosys[Affil] OR Syqe Medical[Affil] OR Vectura Group[Affil] OR Biofourmis[Affil] OR Japan Tobacco[Affil] OR Tularik[Affil] OR Akros Pharma[Affil] OR Torii Pharmaceutical[Affil] OR Imperial Brands[Affil] OR Auxly Cannabis[Affil] \
OR Compania de Distribucion[Affil] OR Oxford Cannabinoid Technologies[Affil] OR British American Tobacco[Affil] OR Philter Labs[Affil] OR Open Book Extracts[Affil] OR Velo[Affil] OR The Kanvas Co[Affil] OR Cannopeia[Affil] OR Purissima[Affil] OR Hesperos[Affil] OR Organigram[Affil] \
OR Trait Biosciences[Affil] OR Ajna Biosciences[Affil] OR Niconovum[Affil] OR Sanity Group[Affil] OR PlatoScience[Affil]' OR Charlotte's Web[Affil] OR KBio[Affil]) AND (2023[Date - Create] : 3000[Date - Create])",'TobaccoCompaniesPubsSince2023']


def removeThe(t:str):
    return t[4:] if t.startswith('The ') else t

class NCBIeutils:
  baseURL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
  
  def __init__(self,query:str):
      self.params = {'db':'pubmed','term':query}
      self.journalCounter = defaultdict(int)
      self.cache_path = NCBI_CACHE
      self.query = query
  
  def _esearch_url(self,params:dict):
      return self.baseURL+'esearch.fcgi?'+urllib.parse.urlencode(params)
  
  def _efetch_url(self,params:dict):
      return self.baseURL+'efetch.fcgi?'+urllib.parse.urlencode(params)
  
  def mydb(self):
      return self.params['db']
  
  def myquery(self):
      return self.query[0]
  

  def get_count(self):
    count_param = self.params = {'db':self.mydb(),'term':self.myquery()}
    count_param.update({'rettype':'count'})
    my_url = self._esearch_url(count_param)
    req = urllib.request.Request(url=my_url)
    response = ET.fromstring(urllib.request.urlopen(req).read())
    return int(str(response.find('./Count').text))
  

  def _retmax_uids(self,params:dict={}):
    '''
    Return
    ------
    [PMIDs] with size < RETMAX sorted in ascending order
    '''
    my_params = {'db':self.mydb(),'term':self.myquery(),'retstart':0,'retmax':RETMAX,'sort':'pub_date'}
    my_url = self._esearch_url(my_params)
    req = urllib.request.Request(url=my_url)
    xml_str = urllib.request.urlopen(req).read()
    response = ET.fromstring(xml_str)
    ids = response.findall('./IdList/Id')
    ids2return = [int(x.text) for x in ids]
    # sort returns PMIDs list in descending order
    reversed_list = []
    [reversed_list.append(ids2return[i]) for i in range(len(ids2return) - 1, -1, -1)]
    return reversed_list
  

  def get_uids(self,query_name:str):
      '''
      Return
      ------
      List of PMIDs sorted by PDAT
      '''
      json_id_dump = query_name+'.json'
      json_id_path = self.cache_path+json_id_dump

      def remove_duplicates(uids:list):
        '''
        keeps uids order in uids
        '''
        seen = set()
        return [id for id in uids if id not in seen and not seen.add(id)]
      
      try:
          all_ids = json.load(open(json_id_path,'r'))
          print(f'Loaded {len(all_ids)} IDs from {json_id_dump}')
          return all_ids
      except FileNotFoundError:
          count = self.get_count()
          all_ids = list()
          if count > RETMAX:
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
                    

  def fetch(self,query_name:str):
    ids = self.get_uids(query_name)
    params = {'db':self.params['db'],'rettype':'XML'}
    stepSize = 200
    for i in range(0, len(ids), stepSize):
      ids = ','.join(str(s) for s in ids[i:i+stepSize])
      params.update({'id':ids})
      my_url = self._efetch_url(params)
      req = urllib.request.Request(url=my_url)
      xml_str = urllib.request.urlopen(req).read()
      yield xml_str


  def path2cache(self,query_name:str):
      return self.cache_path+query_name+".xml"

  def download_pubmed(self,query_name:str):
      """
      Output
      ------
      query_name.json: has list of pmids found by query\n
      query_name.xml: pubmed abstracts
      """
      fpath = self.path2cache(query_name)
      result_counter = 0
      with open(fpath, "a", encoding='utf-8') as file_result:
          file_result.write('<PubmedArticleSet>\n')
          for xml_str in self.fetch(query_name):
              abstractTree = ET.fromstring(xml_str)
              articles = abstractTree.findall('PubmedArticle')
              result_counter += len(articles)
              for article in articles:
                  journal = str(article.find('MedlineCitation/Article/Journal/Title').text)
                  jnames = normalize_journal(journal)
                  for j in jnames:
                      self.journalCounter[j] += 1
                  file_result.write(ET.tostring(article, encoding="unicode"))
          file_result.write('</PubmedArticleSet>\n')
          print(f'Downloaded {result_counter} pubmed abstracts')
          print(f'Downloaded abstracts are in "{fpath}"')

          fstatpath = os.path.join(self.cache_path,query_name+'_stats.tsv')
          with open(fstatpath, "w", encoding='utf-8') as f:
              [f.write(f'{k}\t{v}\n') for k,v in self.journalCounter.items()]
          print(f'Statistics is in "{fstatpath}"')


def pubmed_hyperlink(pmids:list,display_str='',as_excel_formula=True):
    if as_excel_formula:
        ids4hyperlink = pmids[:20] 
        # hyperlink in Excel does not work with long URLs, 
        params = {'term':','.join(ids4hyperlink)}
        if display_str:
            to_display = str(display_str)
        elif len(pmids) == 1:
            to_display = str(pmids[0])
        else:
            to_display = str(len(pmids))

        data = urllib.parse.urlencode(params, quote_via=urllib.parse.quote)
        return '=HYPERLINK("'+PUBMED_URL+data+'",\"{}\")'.format(to_display)
    else:
        ids4hyperlink = pmids[:100]
        params = {'term':','.join(ids4hyperlink)}
        data = urllib.parse.urlencode(params, quote_via=urllib.parse.quote)
        return PUBMED_URL+data


def pmc_hyperlink(pmcs:list,display_str='',as_excel_formula=True):
    if as_excel_formula:
        ids4hyperlink = pmcs[:20] 
        # hyperlink in Excel does not work with long URLs
        terms_string = '+OR+'.join(ids4hyperlink)
        query_string = f'term={terms_string}'

        #params = {'term':ids4hyperlink}
        if display_str:
            to_display = str(display_str)
        elif len(pmcs) == 1:
            to_display = str(pmcs[0])
        else:
            to_display = str(len(pmcs))
        
        return '=HYPERLINK("'+PMC_URL+query_string+'",\"{}\")'.format(to_display)
    else:
        ids4hyperlink = pmcs[:100]
        terms_string = '+OR+'.join(ids4hyperlink)
        query_string = f'term={terms_string}'
        return PMC_URL+query_string


def normalize_journal(journal_title:str):
    jname = removeThe(titlecase(journal_title)).replace(('. '),' ')
    jname = jname.split(' : ')[0]
    jnames = jname.split(' = ')
    all_names = list()
    for j in jnames:
        all_names += j.split(': ')
    return all_names


def medlineTA2issn()->tuple[dict[str,str],dict[str,str]]:
    abbrev2journal = dict()
    journal2issns = defaultdict(list)
    path2cache = NCBI_CACHE+'J_Medline.txt'
    try: 
        m = open(path2cache,'r',encoding='utf-8')
        line = m.readline().strip()
        journal = str()
        j_abbrev = str()
        IsoAbbr = str()
        issn_print = str()
        issn_online = str()
        while line:
            if line.startswith('JrId'):
                if journal:
                    jnames = normalize_journal(journal)
                    norm_journal = jnames[0]
                    abbrev2journal[j_abbrev] = norm_journal
                    if IsoAbbr != j_abbrev:
                        abbrev2journal[IsoAbbr] = norm_journal
                    if issn_print:
                        journal2issns[norm_journal].append(issn_print)
                    if issn_online:
                        journal2issns[norm_journal].append(issn_online)

                    if len(jnames) > 1:
                        for j in jnames:
                            abbrev2journal[j] = norm_journal
                        
                journal = ''
                j_abbrev = ''
                issn_print = ''
                issn_online = ''
                IsoAbbr = ''
            elif line.startswith('JournalTitle'):
                journal = line[len('JournalTitle')+2:]
            elif line.startswith('MedAbbr'):
                j_abbrev = line[len('MedAbbr')+2:]
            elif line.startswith('ISSN (Print)'):
                issn_print = line[len('ISSN (Print)')+2:]
            elif line.startswith('ISSN (Online)'):
                issn_online = line[len('ISSN (Online)')+2:]
            elif line.startswith('IsoAbbr'):
                IsoAbbr = line[len('IsoAbbr')+2:]
            line = m.readline().strip()
    except FileNotFoundError:
        print(f'Cannot find {path2cache}')
        pass
    return abbrev2journal, dict(journal2issns)


def docid2pmid(docids:list[str]):
  '''
  input:
    docids - homogeneous list of PMC or DOI 
    only DOIs for PMC articles will produce valid response
  '''
  docid2pmid = dict()
  idstr = ','.join(docids)
  requesturl = 'https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids='+idstr+'&format=json'
  req = urllib.request.Request(url=requesturl)
  response_dict = json.loads(urllib.request.urlopen(req).read())
  for record in response_dict["records"]:
    pmid = record['pmid']
    pmcid = record['pmcid']
    doi = record.get('doi','')
    if doi in docids:
      docid2pmid[doi] = pmid
    else:
      docid2pmid[pmcid] = pmid
  return docid2pmid


if __name__ == "__main__":
  ncbi = NCBIeutils(query[0])
  ncbi.download_pubmed(query[1])


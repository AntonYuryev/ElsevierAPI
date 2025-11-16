
#C:Windows> py -m pip install entrezpy --user
import os
import xml.etree.ElementTree as ET
from collections import defaultdict
from titlecase import titlecase
from ..utils import sortdict
from ..ETM_API.references import Reference,JOURNAL,TITLE,AUTHORS,PUBYEAR
from .NCBIutils import NCBIeutils


PUBMED_CACHE = os.path.join(os.getcwd(),'ENTELLECT_API/ElsevierAPI/NCBI/__ncbipubmedcache__/')

def removeThe(t:str):
    return t[4:] if t.startswith('The ') else t

class Pubmed(NCBIeutils):
  def __init__(self,query:str,retmode='xml'):
    super().__init__('pubmed',PUBMED_CACHE,retmode)
    self.journalCounter = defaultdict(int)
    self.query = query

  
  @staticmethod
  def pubmed2ref(article:ET.Element)->Reference:
    def parse_authors(AuthorList:ET.Element):
      authors = []
      for author_elem in AuthorList.findall('Author'):
          last_name = author_elem.find('LastName').text
          fore_name = author_elem.find('ForeName').text
          initials = author_elem.find('Initials').text
          # Assuming Initials is always available and accurate for the format:
          formatted_name = f"{last_name} {initials}{fore_name}"
          authors.append(formatted_name)
      return ';'.join(authors)
    
    def parse_year(journal:ET.Element):
      year = journal.find('JournalIssue/PubDate/Year')
      if year is not None:
        return year.text
      else:
        medline_date = journal.find('JournalIssue/PubDate/MedlineDate')
        if medline_date is not None:
          return medline_date.text.split(' ')[0]
        else:
          return ''
       
      
    medline_citation = article.find('MedlineCitation')
    pmid_elem = medline_citation.find('PMID')
    ref = Reference('PMID', pmid_elem.text)
    article_elem = medline_citation.find('Article')
    journal_elem = article_elem.find('Journal')
    ref['ISSN'] = [journal_elem.find('ISSN').text]
    ref[JOURNAL] = [journal_elem.find('ISOAbbreviation').text]
    ref[PUBYEAR] = [parse_year(journal_elem)]
    ref[TITLE] = [article_elem.find('ArticleTitle').text]
    ref[AUTHORS] = [parse_authors(article_elem.find('AuthorList'))]
    return ref
    

  def download_refs(self,pmids:list[int])->list[Reference]:
    refs = []
    for xml_str in self.ids2records(pmids):
      articles = ET.fromstring(xml_str).findall('PubmedArticle')
      [refs.append(NCBIeutils.pubmed2ref(article))for article in articles]
    return refs


  def download(self,query_name:str)->list[Reference]:
      fpath = self.path2cache(query_name)
      refs = []
      with open(fpath, "w", encoding='utf-8') as result:
        result.write('<PubmedArticleSet>\n')
        for xml_str in self.fetch(query_name):
          abstractTree = ET.fromstring(xml_str)
          articles = abstractTree.findall('PubmedArticle')
          for article in articles:
            refs.append(self.pubmed2ref(article))
            journal = str(article.find('MedlineCitation/Article/Journal/Title').text)
            jnames = normalize_journal(journal)
            for j in jnames:
                self.journalCounter[j] += 1
            result.write(ET.tostring(article, encoding="unicode"))
#          sleep(1.0)
        result.write('</PubmedArticleSet>\n')
        print(f'Downloaded {len(refs)} pubmed abstracts')
        print(f'Downloaded abstracts are in "{fpath}"')

      fstatpath = os.path.join(self.cache_path,query_name+'_stats.tsv')
      self.journalCounter = sortdict(self.journalCounter,by_key=False,reverse=True)
      with open(fstatpath, "w", encoding='utf-8') as f:
          [f.write(f'{k}\t{v}\n') for k,v in self.journalCounter.items()]
      print(f'Download statistics is in "{fstatpath}"')
      return refs


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
    path2cache = PUBMED_CACHE+'J_Medline.txt'
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


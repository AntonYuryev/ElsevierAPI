from builtins import len
from .medscan import MedScan
import xlsxwriter,re,time,json,unicodedata
from datetime import timedelta
from urllib.parse import urlencode,quote
from titlecase import titlecase
from ..utils import list2str,sortdict
from collections import defaultdict
from typing import Generator


AUTHORS = 'Authors'
_AUTHORS_ = 'AuthorsObject'
INSTITUTIONS = 'Institutions'
JOURNAL = 'Journal'
MEDLINETA = 'MedlineTA'
PUBYEAR = 'PubYear'
PUBMONTH = 'PubMonth'
PUBDAY = 'PubDay'
SENTENCE = 'msrc'
TITLE = 'Title'
ABSTRACT = 'Abstract'
CLAIMS = 'Claims'
PATENT_APP_NUM = 'Patent Application Number'
PATENT_GRANT_NUM = 'Patent Grant Number'
RELEVANCE = 'Relevance'
MEASUREMENT = 'Measurement'
ETM_CITATION_INDEX = 'ETM Citation index'
SBS_CITATION_INDEX = 'SBS Citation index'
PS_CITATION_INDEX = 'Graph Citation index'
SCOPUS_CI = 'Scopus Citation index'
LOINCID = 'LOINC ID'
THRESHOLD = 'Threshold'
hGRAPHID = 'hGraph ID'
EDMID = 'EDM ID'
IN_OPENACCESS = 'is_openaccess'
PUBLISHER = 'Publisher'
GRANT_APPLICATION = 'Grant Application'
EMAIL = re.compile(r"\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}\b", flags=re.IGNORECASE)
CLINVAR_ID = 'Clinvar RCV ID'
CLINVAR_ACC = 'Clinvar RCV Accession'

INT_PROPS = {PS_CITATION_INDEX,PUBYEAR,ETM_CITATION_INDEX,SBS_CITATION_INDEX}
FLOAT_PROPS = {RELEVANCE}
ARTICLE_ID_TYPES = ['PMID', 'DOI', 'PII', 'PUI', 'EMBASE']
CLINVAR_ID_TYPES = ['CLINVAR','Clinvar',CLINVAR_ACC,CLINVAR_ID]
PS_REFIID_TYPES = ARTICLE_ID_TYPES + ['NCT ID','NCTID'] + CLINVAR_ID_TYPES
#keep PS_ID_TYPES as list for efficient identifier sort.  ID types are ordered by frequency in Resnet
ETM_ID_TYPES = ['ELSEVIER','PMC','REPORTER','GRANTNUMREPORTER']
PATENT_ID_TYPES = [PATENT_APP_NUM, PATENT_GRANT_NUM]
CLINTRIAL_PROPS = {'TrialStatus','Phase','StudyType','Start','Intervention','Condition','Company','Collaborator'}

JOURNAL_PROPS = {JOURNAL,'ISSN','ESSN',MEDLINETA}
PS_BIBLIO_PROPS_ALL = {PUBYEAR,AUTHORS,TITLE,'PubMonth','PubDay','PubTypes','Start'}|JOURNAL_PROPS # 'Start' = 'PubYear'
PS_BIBLIO_PROPS = {TITLE, PUBYEAR, AUTHORS, JOURNAL, MEDLINETA,'Start','ISSN'} # contains only props necessary for reference::__biblio_tuple

BIBLIO_PROPS = PS_BIBLIO_PROPS | {INSTITUTIONS,RELEVANCE}
REF_ID_TYPES = PS_REFIID_TYPES+ETM_ID_TYPES+PATENT_ID_TYPES

PS_SENTENCE_PROPS = [SENTENCE,'Organism','CellType','CellLineName','Organ','Tissue','Source','Percent',THRESHOLD,
                     'pX','Phase','Start','TrialStatus','URL','Experimental System',CLINVAR_ID,'TextRef',
                     CLINVAR_ACC,'Clinvar ID','TextMods','BiomarkerType','ChangeType','QuantitativeType']
SENTENCE_PROPS = PS_SENTENCE_PROPS + ['Evidence','msrc','mref','Similarity']
# SENTENCE_PROPS needs to be a list for ordered printing
#also TextRef - used as key in Reference.snippets

PS_REFERENCE_PROPS = list(CLINTRIAL_PROPS)+PS_REFIID_TYPES+list(PS_BIBLIO_PROPS_ALL)+PS_SENTENCE_PROPS

REFERENCE_PROPS = list(PS_BIBLIO_PROPS_ALL)+list(CLINTRIAL_PROPS)+REF_ID_TYPES+SENTENCE_PROPS

NOT_ALLOWED_IN_SENTENCE='[\t\r\n\v\f]' # regex to clean up special characters in sentences, titles, abstracts

IDENTIFIER_PREFIXES = [
        ('PMID', 'info:pmid/{}'),
        ('DOI', 'info:doi/{}'),
        ('PII', 'info:pii/{}'),
        ('PUI', 'info:pui/{}'),
        ('EMBASE', 'info:embase/{}'),
        ('ELSEVIER', 'info:pii/{}'),
        ('PMC', 'info:pmc/{}'),
        (CLINVAR_ID, 'info:clinvar/{}'),
        ('Clinvar', 'info:clinvar/{}'),
        (CLINVAR_ACC, 'info:clinvar/{}'),
        
      ]

NONSTANDARD_IDENTIFIER_PREFIXES = [
        ('NCT ID', 'info:nctid/{}'),
        ('TextRef', '{}'), # No prefix needed for TextRef itself
        ('REPORTER', 'info:nihreporter/{}'),
        ('GRANTNUMREPORTER', 'info:nihgrant/{}'),
        (PATENT_APP_NUM, 'info:patentapp/{}'),
        (PATENT_GRANT_NUM, 'info:patentgrant/{}'),
        (TITLE, 'info:title/{}')
      ]

PUBMED_URL = 'https://pubmed.ncbi.nlm.nih.gov/?'
PMC_URL = 'https://www.ncbi.nlm.nih.gov/pmc/?'
DOI_URL = 'http://dx.doi.org/'

def pubmed_hyperlink(pmids:list,display_str='',as_excel_formula=True):
    if as_excel_formula:
      ids4hyperlink = list(map(str,pmids[:20]))
      # hyperlink in Excel does not work with long URLs, 
      params = {'term':','.join(ids4hyperlink)}
      if display_str:
        to_display = str(display_str)
      elif len(pmids) == 1:
        to_display = str(pmids[0])
      else:
        to_display = str(len(pmids))

      data = urlencode(params, quote_via=quote)
      return '=HYPERLINK("'+PUBMED_URL+data+'",\"{}\")'.format(to_display)
    else:
      ids4hyperlink = pmids[:100]
      params = {'term':','.join(ids4hyperlink)}
      data = urlencode(params, quote_via=quote)
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


def make_hyperlink(identifier:str,url:str,display_str=''):
        # hyperlink in Excel does not work with long URLs, 
        display_str = display_str if display_str else identifier
        return '=HYPERLINK("'+url+identifier+'",\"{}\")'.format(display_str)


def pii_hyperlink(identifier:str,display_str=''):
  # hyperlink in Excel does not work with long URLs, 
  display_str = display_str if display_str else identifier
  identif = identifier.replace('-','').replace('_','')
  url = 'https://www.sciencedirect.com/science/article/pii/'
  return '=HYPERLINK("'+url+identif+'",\"{}\")'.format(display_str)


def doi_hyperlink(identifier:str,display_str=''):
  # hyperlink in Excel does not work with long URLs, 
  display_str = display_str if display_str else identifier
  url = 'http://doi.org/'
  return '=HYPERLINK("'+url+identifier+'",\"{}\")'.format(display_str)


class Reference(dict):
  '''
  Reference{BIBLIO_PROPS[i]:[values]};\n
  Reference.Identifiers{REF_ID_TYPES[i]:identifier};\n
  Reference.Sentences{TextRef:{SENTENCE_PROPS[i]:{Value}}}\n
  '''
  pass

  def __init__(self, idType:str, ID:str):
    super().__init__(dict()) # self{BIBLIO_PROPS[i]:[values]};
    self.Identifiers = {idType:ID} #from REF_ID_TYPES
    self.snippets = defaultdict(lambda: defaultdict(set)) # {TextRef:{PropID:{Values}}} PropID is from SENTENCE_PROPS contains sentences marked up by NLP
    # PropID:{Values} = defaultdict(set)

  def copy_ref(self):
      '''
      Return
      ------
      Reference object copy if self has valid identifiers\n
      otherwise returns empty Reference() object
      '''
      for i in REF_ID_TYPES:
        try:
          id_value = self.Identifiers[i]
          my_copy = Reference(i,id_value)
          my_copy.update(self)
          my_copy.Identifiers = self.Identifiers.copy()
          my_copy.snippets = self.snippets.copy()
          return my_copy
        except KeyError: continue
      return Reference()
  

  @staticmethod
  def _textref2id(textref:str)->tuple[str,str]:
    '''
    output:
      (idtype,id) if textref valid else ('TextRef', textref)
    '''
    #TexRef example: 'info:doi/10.1016/j.gendis.2015.05.001#body:49'
    prefix = textref[:textref.find(':')]
    if prefix == 'info':
      tr = textref[5:]
      slash_pos = tr.find('/')
      if slash_pos > 0:
        id_type = tr[:slash_pos]
        identifier_start = slash_pos+1
        identifier_end = tr.rfind('#',identifier_start)
        if identifier_end < 0: identifier_end = len(tr)
        identifier = tr[identifier_start:identifier_end]
        return id_type.upper(),identifier
      else:
        return 'TextRef', textref
    else:
      return 'TextRef', textref


  @classmethod
  def from_textref(cls, textref:str):
      id_type, identifier = cls._textref2id(textref)
      return cls(id_type,identifier)
  
  
  @classmethod
  def from_iddict(cls, idtype2id:dict):
      '''
      Input
      -----
      {id_type:id}
      '''
      id_type, identifier = next(iter(idtype2id.items()))
      new_ref = cls(id_type, identifier)
      new_ref.Identifiers.update(idtype2id)
      return new_ref
  

  def __key(self):
      for id_type in REF_ID_TYPES:
          try: return self.Identifiers[id_type]
          except KeyError: continue
      
      try: return self.Identifiers['TextRef']
      #if reference has non-canonical TextRef that cannot be parsed by _textref2id
      # self.Identifiers has ['TextRef'] value
      except KeyError: return NotImplemented


  def __hash__(self):#__hash__ needs __eq__ to work properly
    return hash(self.__key())
  

  def __eq__(self, other:"Reference"):#__hash__ needs __eq__ to work properly
    for id_type in REF_ID_TYPES:
      try:
        return self.Identifiers[id_type] == other.Identifiers[id_type]
      except KeyError:
        continue
 
  
  def same_as(self, other:"Reference"):
    for id_type,identifier in self.Identifiers.items():
      if other.Identifiers[id_type] == identifier:
        return True
    return False
  

  def append_property(self, PropId, PropValue):
      try:
          self[PropId].append(PropValue)
      except KeyError:
          self[PropId] = [PropValue]


  def get_equivalent(self, ref_list:set):
      for ref in ref_list:
          if ref == self: return ref
      return dict()


  def my_sentence_props(self)->set[str]:
      prop_names = set()
      for prop2value in self.snippets.values():
          prop_names.update(prop2value.keys())
      return prop_names


  def update_with_value(self, PropId, PropValue:int|str):
      clean_prop = PropValue.strip(' .\n') if isinstance(PropValue,str) else PropValue
      try:
          my_props = set(self[PropId])
          my_props.add(clean_prop)
          self[PropId] = list(my_props)
      except KeyError:
          self[PropId] = [clean_prop]


  @staticmethod
  def __clean_vals(prop:str,values:list):
    if prop == _AUTHORS_:
       return set(values)
    elif prop in INT_PROPS:
      return set(map(int,values))
    elif prop in FLOAT_PROPS:
      return set(map(float,values))
    else:
      stripped_vals = {str(v).strip(' .') for v in values if v is not None}
      return set(filter(None,stripped_vals))


  def update_with_list(self, prop:str, with_values:list):
    clean_vals = self.__clean_vals(prop,with_values)
    try:
      self[prop] = list(set(self[prop])|clean_vals)
    except KeyError:
      self[prop] = list(clean_vals)


  def add_sentence_prop(self, text_ref:str, propID:str, prop_value:str):
    if propID == SENTENCE:
      prop_value = prop_value.strip(' .')

    self.snippets[text_ref][propID].add(prop_value)



  def number_of_sentences(self):
    count = 0
    for snippet in self.snippets.values():
      count += len(snippet.get(SENTENCE,{})) #
    return count


  def add_snippet(self,textref:str,snippet:dict[str,set]):
    '''
    input:
      snippet = {prop_name:{values}}
    '''
    [self.add_sentence_props(textref,k,list(v)) for k,v in snippet.items()]
  
    
  def get_sentence(self,textref:str):
    if textref in self.snippets:
      snippets = self.snippets[textref]
      return next(iter(snippets[SENTENCE])) if SENTENCE in snippets else ''
    return ''


  def __is_new(self,sentence:str):
    clean_sent = sentence.strip(' .')
    sent_no_white = clean_sent.replace(" ", "")
    for prop2vals in self.snippets.values():
      try:
        my_sentences = prop2vals[SENTENCE]
        for sentence in my_sentences:
          if sent_no_white == sentence.replace(" ", ""): ## sometime whitespaces are screwedup
            return ''
      except KeyError: continue
    return clean_sent


  def add_sentence_props(self, TextRef:str, propID:str, prop_values:list):
    if propID == SENTENCE:
      prop_values = list(filter(None,[self.__is_new(x) for x in prop_values if x]))
               
    if prop_values:
      self.snippets[TextRef][propID].update(prop_values)


  def has_property(self, prop_name:str):
    for prop2values in self.snippets.values():
      if prop_name in prop2values.keys():
        return True
    return prop_name in self.keys()
  

  def get_values(self,prop_name:str):
    try:
      return list(self[prop_name])
    except KeyError: 
      prop_vals = set()
      for prop2values in self.snippets.values():
        if prop_name in prop2values:
          prop_vals.update(prop2values[prop_name])
      return list(prop_vals)

  
  def has_values_in(self, in_prop2values:dict,case_sensitive=False):
    '''
    Input
    -----
    prop2values = {prop_name:[values]}
    '''
    # now seaching among ref properties
    for prop, values in self.items():
      try:
        match_values = set(in_prop2values[prop])
        if case_sensitive:
          search_set = values
        else:
          match_values =set(map(lambda x: x.lower(),match_values))
          search_set = set(map(lambda x: x.lower(),values))

        if not match_values.isdisjoint(search_set): 
          return True
        else: continue
      except KeyError: continue

    # if nothing found seach among snippet properties 
    for my_prop2values in self.snippets.values():
      for prop, self_snippet_values in my_prop2values.items():
        try:
          match_values = set(in_prop2values[prop])
          if case_sensitive:
            search_set = values
          else:
            match_values =set(map(lambda x: str(x).lower(),match_values))
            search_set = set(map(lambda x: str(x).lower(),self_snippet_values))
          if not match_values.isdisjoint(search_set): 
              return True
          else: continue
        except KeyError: continue

    return False
  

  def remove_props(self, prop_names:list):
      my_copy = self.copy_ref()
      if my_copy:
        for prop2values in my_copy.snippets.values():
          [prop2values.pop(p,'') for p in prop_names]
        
        [my_copy.pop(p,'') for p in prop_names]
        [my_copy.Identifiers.pop(p,'') for p in prop_names]
        return my_copy
      return Reference()
  

  def rename_prop(self, old_prop_name:str, new_prop_name:str):
      '''
      Rename property in self and in self.snippets
      '''
      if old_prop_name in self:
        self[new_prop_name] = self.pop(old_prop_name)
        return True
      
      was_renamed = False
      for snippet_props in self.snippets.values():
        if old_prop_name in snippet_props:
          snippet_props[new_prop_name] = snippet_props.pop(old_prop_name)
          was_renamed = True

      return was_renamed
  

  def get_doc_id(self):
      '''
      output:
        tuple(id_type, identifier) for the first id type from REF_ID_TYPES\n
        tuple('','') if reference does not have id_type in REF_ID_TYPES
      '''
      for id_type in REF_ID_TYPES:
        if id_type in self.Identifiers:
          return  id_type, self.Identifiers[id_type]
      return str(),str()
  

  def doi_or_id(self):
    if 'DOI' in self.Identifiers:
      return self.Identifiers['DOI']
    else:
      id_type, identifier = self.get_doc_id()
      return id_type  +':'+identifier if id_type  else ''

  
  @staticmethod
  def identifiers_str(id_type:str,identifier:str):
      return id_type + ':' + identifier
  

  def _identifiers_str(self,id_type=''):
      '''
      output:
        if id_type provided and exists in self.Identifiers: id_type:id_value\n
        else: first id_type:id_value from REF_ID_TYPES
      '''
      if id_type in self.Identifiers:
        return id_type+':'+self.Identifiers[id_type]
      else:
        id_type, identifier = self.get_doc_id() 
        return id_type  +':'+identifier if id_type  else ''
  

  def identifier(self,identifier_type):
      try:
          return self.Identifiers[identifier_type]
      except KeyError:
          return ''
  

  def biblioprop2str(self,propid:str):
    if propid in self:
      map_func = Author.tostr if propid == _AUTHORS_ else str
      return ';'.join(list(map(map_func,self[propid])))
    else:
      return ''


  def to_list(self,id_types=list(),print_snippets=False,biblio_props=list(),other_props=list(),with_hyperlinks=False):
      '''
      Return
      ------
      order of properties in return list: other_props, id_types, biblio_props, snippets\n
      reference identifiers for PMID and DOI are hyperlinked if with_hyperlinks is True
      snippets are printed as json dump in one column
      '''
      row = list()
      id_types = id_types if isinstance(id_types,list) else ['PMID']
      for p in other_props:
          try:
              prop_values_str = ';'.join(list(map(str,self.get_props(p))))
              row.append(prop_values_str)
          except KeyError:
              row.append('')

      if id_types:
        for t in id_types:
          try:
            identifier = self.Identifiers[t]
            if with_hyperlinks:
                if t == 'PMID':
                    identifier = pubmed_hyperlink([identifier])
                elif t == 'DOI':
                    identifier = doi_hyperlink(identifier)
                elif t == 'PMC':
                    identifier = pmc_hyperlink([identifier])
                elif t == PATENT_APP_NUM:
                    identifier = make_hyperlink(identifier,'https://patents.google.com/patent/')
            row.append(identifier)
          except KeyError:
            row.append('')
      else:
          row.append(self._identifiers_str())
              
      for prop_id in biblio_props:
        try:
          prop_values_str = self.biblioprop2str(prop_id)
          if prop_id in ['Title', 'Abstract']:
            prop_values_str = re.sub(NOT_ALLOWED_IN_SENTENCE,' ',prop_values_str)
        except KeyError:
            prop_values_str = ''
        row.append(prop_values_str)

      if print_snippets:
        list_snippets = {k:{p:list(l)} for k,v in self.snippets.items() for p,l in v.items()}
        sentence_props = json.dumps(list_snippets)
        sentence_props = re.sub(NOT_ALLOWED_IN_SENTENCE,' ',sentence_props)
        row.append(sentence_props)

      return row
  

  def toAuthors(self):
    """
      creates self[_AUTHORS_] from self[AUTHORS]
    """
    if _AUTHORS_ not in self:
      author_strs = list(filter(None,self.author_list()))
      self[_AUTHORS_] = list(map(Author.fromStr,author_strs))
    return
    

  def to_str(self,id_types=list(),col_sep='\t',print_snippets=False,biblio_props=[],other_props=[],with_hyperlinks=False):
    '''
    Return
    ------
    order of properties in return list: other_props, id_types, biblio_props, snippets\n
    reference identifiers for PMID and DOI are hyperlinked if with_hyperlinks is True
    snippets are printed as json dump in one column
    '''
    row = self.to_list(id_types,print_snippets,biblio_props,other_props,with_hyperlinks)
    return col_sep.join(row)
  

  def pubyear(self):
      if PUBYEAR in self:
        return int(self[PUBYEAR][0])
      else:
        if 'Start' in self: # Clinical trials case
          year = str(self['Start'][0]) 
          if len(year) >= 4 and year[-4:].isdigit(): # format: 20-Apr-2020; August 2004 
              return int(year[-4:])
          elif year[-2:].isdigit():
              return int('20'+year[-2:]) # format: 20-Apr-20
          elif year[:2].isdigit(): 
              return int('20'+year[:2]) # format 20-May
          else:
              print(f'Unknown Clinical trial "Start" format: {year}')
              return 1812
        else: return 1812 # No PubYear case


  def title(self):
    return self.get(TITLE,[''])[0]

  
  def pmid(self):
    return self.Identifiers.get('PMID','')
  

  def pubmed_link(self):
    if 'PMID' in self.Identifiers:
      pmid = self.Identifiers['PMID']
      return pubmed_hyperlink([pmid],pmid)
    else:
      return ''


  def doi(self):
    return self.Identifiers.get('DOI','')

  
  def doi_link(self):
    if 'DOI' in self.Identifiers:
      doi = self.Identifiers['DOI']
      return doi_hyperlink(doi,doi)
    else:
      return ''


  def number_of_snippets(self):
      return len(self.snippets)


  def textrefs(self):
      return list(self.snippets.keys())


  def journal(self): 
    return self.get(JOURNAL,['No journal name'])[0]
                
      
  def relevance(self,score_name:str=RELEVANCE):
    try:
      return max(map(float,self[score_name]))
    except KeyError:
      return 0.0
  

  def is_clinical_trial(self):
      try:
          return self.Identifiers['NCT ID']
      except KeyError:
          return ''
      

  def author_list(self)->list[str]:
    '''
    output:
      sorted list of authors from self[AUTHORS]
    '''
    if AUTHORS in self:
      individual_authors = list({au for au_list in self[AUTHORS] for au in str(au_list).split(';') if au})
      individual_authors.sort()
      self[AUTHORS] = individual_authors
      return individual_authors
    else:
      return []


  def _biblio_tuple(self):
      """
      Returns
      --------
      tuple str(title+' ('+pubyear+'). journal.'+authors), identifier_type, identifier 
      """
      try:
          title = self[TITLE][0]
          if title[-1] != '.': title += '.'
      except KeyError:
          title = 'No title available'

      year = self.pubyear()
      year = 'year unknown' if year == 1812 else str(year)

      try:
        authors = self[_AUTHORS_]
        authors_list = [x.tostr() for x in authors if isinstance(x,Author)]
        authors_str = ','.join(authors_list[:3]) if authors_list else 'unknown authors.'
        if len(authors_list) > 3:
          authors_str = authors_list[0]+' et al.'
        else:
          if authors_str[-1] != '.': authors_str += '.'
      except KeyError:
        authors_str = 'unknown authors.'

      journal = self.journal()

      biblio_str = title+' ('+year+'). '+journal+'. '+authors_str
      if self.get(IN_OPENACCESS,False):
          biblio_str += ' [Open access]'
      identifier_type, identifier  = self.get_doc_id()
      return biblio_str, identifier_type, identifier 

  
  def get_biblio_str(self, sep='\t'):
      biblio, id_type,identifier = self._biblio_tuple()
      return biblio+sep+id_type+':'+identifier
  

  def _merge(self, other:"Reference"):
    for prop, values in other.items():
      self.update_with_list(prop,values)

    self.Identifiers.update(other.Identifiers)

    for textref, other_snippet_p2v in other.snippets.items():
      clean_snippet = {k:{str(p).strip(' .') for p in v} for k,v in other_snippet_p2v.items()}
      [self.add_sentence_props(textref,p,v) for p,v in clean_snippet.items()]
    del other
    return


  def is_from_abstract(self):
      for textref in self.snippets.keys():
          try:
              return bool(str(textref).rindex('#abs',-8,-3))
          except ValueError:
              continue
      return False


  def _make_standard_textref(self):
    for identifier, prefix_format in IDENTIFIER_PREFIXES:
      value = self.Identifiers.get(identifier,'')
      if value:
        return prefix_format.format(value)
    return ''


  def _make_non_standard_textref(self):
    for identifier, prefix_format in NONSTANDARD_IDENTIFIER_PREFIXES:
      value = self.Identifiers.get(identifier,'')
      if value:
        return prefix_format.format(value)
    return NotImplemented


  def _make_textref(self):
      textref = self._make_standard_textref()
      return textref if textref else self._make_non_standard_textref()


  def set_weight(self, weight:float,weight_name='weight'):
      '''
          create property "weight_name" with value = weight
          old value is replaced if it is smaller than input "weight_name"
      '''
      try:
          w = float(self[weight_name][0])
          if w < weight:
              self[weight_name] = [weight]
      except KeyError:
          self[weight_name] = [weight]


  def  add_weight(self, weight:float,weight_name='weight'):
      '''
      adds:
          weight to property "weight_name"
      '''
      assert(weight <= 1.00)
      try:
          w = float(self[weight_name][0])
          self[weight_name] = [w+weight]
      except KeyError:
          self[weight_name] = [weight]


  def get_weight(self,weight_name='weight'):
      try:
          return self[weight_name][0]
      except KeyError:
          return 0.0 # allows to exclude certain references from counting if they are not annotated by weight_name


  def _sort_key(self, by_property, is_numerical=True):
      if by_property == PUBYEAR: return self.pubyear()
      else:
          try:
              return float(self[by_property][0]) if is_numerical else str(self[by_property][0])
          except KeyError:
              return 0.0 if is_numerical else '0'


  def sentences(self)->Generator[tuple[str,str], None, None]:
    for textref, snippet in self.snippets.items():
      sentences = snippet.get(SENTENCE,{''})
      for sentence in sentences:
        yield textref,sentence


  def _snippets(self):
    for textref, snippet in self.snippets.items():
      sentences = snippet.get(SENTENCE,set())
      if len(sentences) > 1:
        for i,sentence in enumerate(sentences):
          new_snippet = snippet.copy()
          new_snippet[SENTENCE] = {sentence}
          yield f'{textref}:{i}',new_snippet
      else:
        yield textref,snippet


  def get_snippet_prop(self,prop_name:str):
      """
      output:
        {textref:[prop_vals]}
      """
      prop_values = dict()
      for textref, sentence_props in self.snippets.items():
        if prop_name in sentence_props:
          prop_values[textref] = sentence_props[prop_name]
      return prop_values
  

  def get_props(self,prop_name:str)->list[str]:
    '''
    input:
        prop_name can be either in self.snippet or self.Identifiers or self
    '''
    value = self.get(prop_name, None)
    if value is not None:
      return list(value)

    # First check failed, try .Identifiers
    value = self.Identifiers.get(prop_name,None)
    if value is not None:
      return [value]

    # Both checks failed, run the fallback logic
    snippet_prop_dic = self.get_snippet_prop(prop_name)
    props = set()
    [props.update(prop_vals) for prop_vals in snippet_prop_dic.values()]
    props = list(props)
    props.sort()
    return props


  def get_prop(self,prop_name:str,value_index=0,if_missing_return='')->str|int|bool:
    '''
    input:
        prop_name can be either in self.snippet or self.Identifiers or self
    '''
    my_props = self.get_props(prop_name)
    if my_props:
      try:
        return my_props[value_index]
      except IndexError:
        return if_missing_return
    return if_missing_return
  

  def todict(self,as_str2str=False,relname='')->tuple[str,dict[str,list]]:
    dic = dict(self)
    dic.update({k:[v] for k,v in self.Identifiers.items()})

    add2snippet = relname+'\n' if relname else ''
    my_snippets = sortdict(dict(self.snippets))
    for i, textref2props in enumerate(my_snippets.items()):
        textref = textref2props[0]
        props = textref2props[1]
        props_str = str()
        for prop, values in props.items():
          values_str = ' '.join(values).strip(' .\n,')
          props_str += prop +': '+ values_str + '. '
        dic[f'[{i+1}]Snippet'] = [add2snippet+textref+': '+props_str.strip()]

    return self.doi_or_id(), list2str(dic) if as_str2str else self.doi_or_id(), dic


class ReferenceEncoder(json.JSONEncoder):
  def default(self, obj):
      if isinstance(obj, Reference):
          # Convert CustomObject to a dictionary
          dump_dict = dict(obj)
          dump_dict.update({'Identifiers':obj.Identifiers})
          if obj.snippets:
              dump_dict.update({'snippets':obj.snippets})
          dump_dict['hook'] = 'reference'
          return {dump_dict} # creat hook for ReferenceDecoder
      else:
          super().default(obj)


class ReferenceDecoder(json.JSONDecoder):
  def __init__(self, *args, **kwargs):
      json.JSONDecoder.__init__(self, object_hook=self.object_hook, *args, **kwargs)

  def object_hook(self, obj):
      hook = obj.pop('hook','')
      if hook == 'reference':
        Identifiers = obj.pop('Identifiers')
        snippets = obj.pop('snippets',dict())
        addresses = obj.pop('addresses',dict())

        id_type,id = next(iter(obj.items()))
        ref = Reference(id_type,id)
        ref.update(obj)
        ref.Identifiers.update(Identifiers)
        ref.snippets.update(snippets)
        return ref
      return obj

#########################################DocMine#############################################DocMine##################

INSTITUTION_KEYWORDS = {'institute', 'institut', 'instituto', 'istituto','clinic', 'klinik','genetics','hospital', 'university', 'universitat', 'universiti', 'université','università','universidade','universidad','universitï¿½','universitaria','centre', 'center', 'centro', 'inc', 'colleges',
                      'ltd','gmbh', 'llc', 'school','politecnic','politecnico', 'college', 'department', 'departamento','division', 'council', 'academy','faculty', 'co', 'corp','ministry','campus','group','corporation', 'consulting', 'pharma','union',
                      'laboratory', 'laboratoire', 'laboratories','society','project', 'association', 'commission','program','drug', 'labs', 'biolabs', 'biolab', 'lab', 'unit', 'pharmaceutical', 'policlinico', 'cátedra', 'sciences', 'r&d', 'research', 'trust', 'fund', 'plant', 'biomed','salud', 'health','foundation','federation','service','cluster','fondazione','bio','pharmaceuticals', 'farmaceutico'}

INSTITUTION_KEYWORDS = {unicodedata.normalize('NFKD', word).casefold() for word in INSTITUTION_KEYWORDS} 

def removeThe(t:str):
  return t[4:] if t.startswith('The ') else t


class Author:
  def __init__(self,LastName:str,_1stName='',middle_name=''):
    self.LastName = LastName
    self.MiddleName = middle_name
    self._1stName = _1stName
    self.affiliations = dict() # {organization:address}
    self.email = ''

  def institutions(self):
    return list(self.affiliations.keys())
  
  def __key(self):
    return tuple([self._1stName,self.MiddleName,self.LastName] + self.institutions())

  def __hash__(self):#__hash__ needs __eq__ to work properly
    return hash(self.__key())
  

  def __eq__(self, other:"Author"):#__hash__ needs __eq__ to work properly
    my_key = self.__key()
    other_key = tuple([other._1stName,other.MiddleName,other.LastName] + other.institutions())
    return my_key == other_key


  def name(self):
    return self._1stName[0].upper() +self.LastName.capitalize()
  

  def tostr(self):
    return self._1stName[0]+' '+self.LastName if self._1stName else self.LastName
     
     
  @classmethod
  def fromStr(cls,author:str):
    '''
    assumes "author" has format "FirstName LastName"
    '''
    author_str = author.strip(' .,')
    whitepos = author_str.find(' ')
    if whitepos > 0:
      return Author(author[whitepos+1:],author[:whitepos])
    else:
      last_comma_pos = author_str.rfind(',')
      if last_comma_pos > 0:
        return Author(author_str[:last_comma_pos], author_str[last_comma_pos+1:])
    
    return Author(author_str)
    

class DocMine (Reference):
    '''
    use this class to normalize annotations from references in disparate sources, e.g. ETM,USPO,EPO,Pubmed\n
    DocMine stores text in sections: Abstract, Results,Discussion, Claims, Descriptions
    '''
    def __init__(self, doc_id_type:str, doc_id:str):
        super().__init__(doc_id_type, doc_id)
        self.sections = dict() #{{abstract:[text]}, {claims:[claims]}. Texts to be markerd up by NLP.
        #self.authors = list() # [Author]
        self.addresses = list() # {orgname:adress}


    @staticmethod
    def __textref_suffix(section_name:str):
        name2suffix = {TITLE:'title', ABSTRACT:'abs',CLAIMS:'claims'}
        try:
            return name2suffix[section_name]
        except KeyError:
            return 'cont'


    @staticmethod
    def normalize_journal(journal_title:str):
        return removeThe(titlecase(journal_title)).replace(('. '),' ')


    def add2section(self,section_name:str, paragraph:str): 
        if not paragraph: 
            #print('%s in %s is empty' % (section_name, self.get_title()))
            return
        if not isinstance(paragraph,str):
            #print('%s has no section %s' % (self.get_title(),section_name))
            return
        try:
            self.sections[section_name].append(paragraph)
        except KeyError:
            self.sections[section_name] = [paragraph]


    def set_date (self, year, month='', day=''): 
        self[PUBYEAR] = [year]
        if month: self[PUBMONTH] = [month]
        if day: self[PUBDAY] = [day]


    def _set_title(self, title:str):
        self[TITLE] = [title]
        self.add2section(TITLE,title)

    def medscan_annotate(self,medscan:MedScan):
        base_text_ref = self._make_textref()
        for secname, paragraphs in self.sections.items():
            textref_suf = self.__textref_suffix(secname)
            sentence_idx = 1
            for paragraph in paragraphs:
                paragraph_annotation = medscan.find_concepts(paragraph) # paragraph_annotation = {snippet:{id_range:{id:obj_name}}}
                for sentence_markup, range2dict in paragraph_annotation.items():
                    if range2dict:
                        text_ref = base_text_ref+'#'+textref_suf+':'+str(sentence_idx)
                        self.add_sentence_prop(text_ref,SENTENCE,sentence_markup)
                        for msid_range, concept_dict in range2dict.items():
                            prop_name = medscan.get_concept_type(msid_range)
                            self.add_sentence_props(text_ref,prop_name,list(concept_dict.values()))

                    sentence_idx +=1


    def get_annotations(self, prop_name:str):
        if prop_name == JOURNAL:
            return [self.journal()]
        
        try:
            return list(self[prop_name])
        except KeyError:
            return []


    def journal(self): 
        '''
        Return
        ------
        self[JOURNAL][0] normalized by titlecase and removeThe
        '''
        try:
            journal_name = self.normalize_journal(self[JOURNAL][0])
        except KeyError:
            try:
                journal_name = self.normalize_journal(self[MEDLINETA][0])
            except KeyError:
                return str('No journal name')
                  
        return journal_name
    

    def journal_publisher(self):
        '''
        Return
        ------
        self.journal() (publisher)
        '''
        j_name = self.journal()
        try:
            publisher = str(self[PUBLISHER][0])
            return j_name +'('+publisher+')'
        except KeyError:
            return j_name
        

    def organizations(self):
        return {a.organization for a in self[_AUTHORS_] for i in a.organization}
    
    
    def count_property(self,counter:dict, prop_name:str):
        prop_values = self.get_annotations(prop_name)
        for v in prop_values:
            try:
                current_count = counter[v]
                counter[v] = current_count+1
            except KeyError:
                counter[v] = 1


    @staticmethod
    def dict2worksheet(workbook:xlsxwriter.Workbook, worksheet_name, dict2print:dict, header:list, key_col_width=50):
        worksheet = workbook.add_worksheet(worksheet_name)
        col = 0
        for h in header:
            worksheet.write(0, col, h)
            col +=1
            
        row = 1
        for k,v in dict2print.items():
            worksheet.write(row, 0, k)
            worksheet.write(row, 1, v)
            row += 1

        worksheet.set_column('A:A', key_col_width)

    @staticmethod
    def has_institution_keyword(name:str):
        name_no_punctuation = name.replace('.', ' ').replace(',', ' ')
        name_no_punctuation = name_no_punctuation.strip(' ()[]')
        name_words = name_no_punctuation.split(' ')
        for w in name_words:
            if unicodedata.normalize('NFKD',w).casefold() in INSTITUTION_KEYWORDS: return True
        return False

    @staticmethod
    def execution_time(execution_start):
        return "{}".format(str(timedelta(seconds=time.time() - execution_start)))

            
def reflist2dict(refs:list[Reference])->dict[str,Reference]:
  '''
  output:
    {'id_type:identifier':Reference}
  '''
  refdict = dict()
  for ref in refs:
    identifiers = [ref.identifiers_str(k,v) for k,v in ref.Identifiers.items()]
    was_merged = False
    for i in identifiers:
      try:
        refdict[i]._merge(ref)
        was_merged = True
        break
      except KeyError:
        continue
    if not was_merged:
      refdict.update({i:ref for i in identifiers})
    
  return refdict



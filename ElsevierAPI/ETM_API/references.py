from .medscan import MedScan
import xlsxwriter
import re
from datetime import timedelta
import time
import json


AUTHORS = 'Authors'
INSTITUTIONS = 'Institutions'
JOURNAL = 'Journal'
MEDLINETA = 'MedlineTA'
PUBYEAR = 'PubYear'
PUBMONTH = 'PubMonth'
PUBDAY = 'PubDay'
SENTENCE = 'Sentence'
TITLE = 'Title'
ABSTRACT = 'Abstract'
CLAIMS = 'Claims'
PATENT_APP_NUM = 'Patent Application Number'
PATENT_GRANT_NUM = 'Patent Grant Number'
EFFECT = 'Effect'
RELEVANCE = 'Relevance'
EMAIL = re.compile(r"\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}\b", flags=re.IGNORECASE)

PS_ID_TYPES = ['PMID', 'DOI', 'PII', 'PUI', 'EMBASE','NCT ID'] 
#keep PS_ID_TYPES as list for efficient identifier sort.  ID types are ordered by frequency in Resnet
ETM_ID_TYPES = ['ELSEVIER','PMC','REPORTER','GRANTNUMREPORTER']
PATENT_ID_TYPES = [PATENT_APP_NUM, PATENT_GRANT_NUM]
CLINTRIAL_PROPS = {'TrialStatus','Phase','StudyType','Start','Intervention','Condition','Company','Collaborator'}
JOURNAL_PROPS = {JOURNAL,'ISSN','ESSN',MEDLINETA}
PS_BIBLIO_PROPS = {PUBYEAR,AUTHORS,TITLE,'PubMonth','PubDay','PubTypes','Start'}|JOURNAL_PROPS # 'Start' = 'PubYear'
BIBLIO_PROPS = PS_BIBLIO_PROPS | {INSTITUTIONS,RELEVANCE}
REF_ID_TYPES = PS_ID_TYPES+ETM_ID_TYPES+PATENT_ID_TYPES

SENTENCE_PROPS = [SENTENCE,'Organism','CellType','CellLineName','Organ','Tissue','Source','Percent','Evidence']
# SENTENCE_PROPS needs to be a list for ordered printing
#also TextRef - used as key in Reference.Sentences
RELATION_PROPS = {EFFECT,'Mechanism','ChangeType','BiomarkerType','QuantitativeType'}


NOT_ALLOWED_IN_SENTENCE='[\t\r\n\v\f]' # regex to clean up special characters in sentences, titles, abstracts

class Reference(dict):  
# self{BIBLIO_PROPS[i]:[values]}; Identifiers{REF_ID_TYPES[i]:identifier}; Sentences{TextRef:{SENTENCE_PROPS[i]:Value}}
    pass

    def __init__(self, idType:str, ID:str):
        super().__init__(dict()) # self{BIBLIO_PROPS[i]:[values]};
        self.Identifiers = {idType:ID} #from REF_ID_TYPES
        self.snippets = dict() # {TextRef:{PropID:[Values]}} PropID is from SENTENCE_PROPS contains sentences marked up by NLP
        self.addresses = dict() # {orgname:adress}

    @classmethod
    def copy(cls, other):
        if isinstance(other, Reference):
            for i in REF_ID_TYPES:
                try:
                    id_value = other.Identifiers[i]
                    cls(i,id_value)
                    cls.Identifiers.update(other.Identifiers)
                    cls.update(other)
                    cls.snippets.update(other.snippets)
                    cls.addresses.update(other.addresses)
                    return cls
                except KeyError: continue

        return dict()

        
    @staticmethod
    def __parse_textref(textref:str):
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
        id_type, identifier = cls.__parse_textref(textref)
        return cls(id_type,identifier)

    def __key(self):
        for id_type in REF_ID_TYPES:
            try: return self.Identifiers[id_type]
            except KeyError: continue
        
        try: return self.Identifiers['TextRef']
        #if reference has non-canonical TextRef that cannot be parsed by __parse_textref
        # self.Identifiers has ['TextRef'] value
        except KeyError: return NotImplemented

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other):
        if isinstance(other, Reference):
            for id_type in REF_ID_TYPES:
                try:
                    return self.Identifiers[id_type] == other.Identifiers[id_type]
                except KeyError:
                    continue
    
    def append_property(self, PropId, PropValue):
        try:
            self[PropId].append(PropValue)
        except KeyError:
            self[PropId] = [PropValue]

    def get_equivalent(self, ref_list:set):
        for ref in ref_list:
            if ref == self: return ref

        return dict()

    def update_with_value(self, PropId, PropValue:str):
        clean_prop = PropValue.strip(' .')
        self_clean_props = set(map(lambda x: str(x).strip(' .'), PropValue))
        try:
            self[PropId] = list(set(self_clean_props|{clean_prop}))
        except KeyError:
            self[PropId] = [clean_prop]

    def update_with_list(self, prop_id, prop_values: list):
        try:
            self[prop_id] = list(set(self[prop_id] + prop_values))
        except KeyError:
            self[prop_id] = list(set(prop_values))

    def add_sentence_prop(self, text_ref:str, propID:str, prop_value:str):
        try:
            exist_sentence = self.snippets[text_ref]
            try:
                exist_sentence[propID].append(prop_value)
            except KeyError:
                exist_sentence[propID] = [prop_value]
            self.snippets[text_ref] = exist_sentence
        except KeyError:
            self.snippets[text_ref] = {propID:[prop_value]}

    def add_sentence_props(self, text_ref:str, propID:str, prop_values:list):
        try:
            exist_sentence = self.snippets[text_ref]
            try:
                vals = set(exist_sentence[propID]+ prop_values)
                exist_sentence[propID] = list(vals)
            except KeyError:
                exist_sentence[propID] = prop_values
            self.snippets[text_ref] = exist_sentence
        except KeyError:
            self.snippets[text_ref] = {propID:prop_values}


    def has_properties(self, prop_names2values:dict,case_sensitive=False):
        # prop_names2values = {prop_name:[values]}
        for prop2values in self.snippets.values():
            for prop, values in prop2values.items():
                try:
                    match_values = set(prop_names2values[prop])
                    # assumes match_values are in lowercase for case insensitive search
                    if case_sensitive:
                        search_set = values
                    else:
                        search_set = set(map(lambda x: x.lower(),values))
                    if not match_values.isdisjoint(search_set): 
                        return True
                    else: continue
                except KeyError: continue

        for prop, values in self.items():
            try:
                match_values = set(prop_names2values[prop])
                if case_sensitive:
                    search_set = values
                else:
                    search_set = set(map(lambda x: x.lower(),values))
                if not match_values.isdisjoint(search_set): 
                    return True
                else: continue
            except KeyError: continue

        return False
           

    def __identifiers_str(self,col_sep='\t'):
        for i in REF_ID_TYPES:
            try:
                return i+':'+self.Identifiers[i]
            except KeyError: 
                continue
        return 'no identifier available'
    
    def to_str(self,id_types:list=None,col_sep='\t',print_snippets=False,biblio_props=[]):
        id_types = id_types if isinstance(id_types,list) else ['PMID']
        row = str()
        if isinstance(id_types,list):
            for t in id_types:
                try:
                    row = row + t+':'+self.Identifiers[t]+col_sep
                except KeyError:
                    row = row + col_sep
        else:
            row = self.__identifiers_str() + col_sep


        row = row[:len(row)-1] #remove last separator

        for prop_id in biblio_props:
            try:
                prop_values_str = ';'.join(self[prop_id])
                if prop_id in ['Title', 'Abstract']:
                    prop_values_str = re.sub(NOT_ALLOWED_IN_SENTENCE, ' ', prop_values_str)
            except KeyError:
                prop_values_str = ''
            
            row = row + col_sep + prop_values_str
            #row = row + col_sep + prop_id + ':' + prop_values_str

        if print_snippets:
            sentence_props =  json.dumps(self.snippets)
            sentence_props = re.sub(NOT_ALLOWED_IN_SENTENCE, ' ', sentence_props)
            row += col_sep + sentence_props

        return row


    def get_biblio_str(self):
        try:
            title = self[TITLE][0]
            if title[-1] != '.': title += '.'
        except KeyError:
            title = 'No title available'
        try:
            year = str(self[PUBYEAR][0])
        except KeyError:
            year = 'year unknown'
        try:
            authors_list = self[AUTHORS][:3]
            authors_str = ','.join(authors_list) if authors_list else 'unknown authors.'
            if len(self[AUTHORS]) > 3:
                authors_str = authors_list[0]+' et al.'
            else:
                if authors_str[-1] != '.': authors_str += '.'
        except KeyError:
            authors_str = 'unknown authors.'
        return title+' ('+year+'). '+authors_str, self.__identifiers_str()


    def _merge(self, other):
        if isinstance(other, Reference):
            for prop, values in other.items():
                if isinstance(values[0], str):
                    clean_vals = [v.strip(' .') for v in values]
                else:
                    clean_vals = values
                self.update_with_list(prop,clean_vals)

            self.Identifiers.update(other.Identifiers)

            for textref, prop2val in other.snippets.items():
                for prop, values in prop2val.items():
                    self.add_sentence_props(textref,prop,values)
            self.addresses.update(other.addresses)

    def is_from_abstract(self):
        for textref in self.snippets.keys():
            try:
                return bool(str(textref).rindex('#abs', len(textref) - 8, len(textref) - 3))
            except ValueError:
                continue
        return False

    def _make_standard_textref(self):
        try:
            return 'info:pmid/'+ self.Identifiers['PMID']
        except KeyError:
            try:
                return 'info:doi/'+ self.Identifiers['DOI']
            except KeyError:
                try:
                    return 'info:pii/'+ self.Identifiers['PII']
                except KeyError:
                    try:
                        return 'info:pui/'+ self.Identifiers['PUI']
                    except KeyError:
                        try:
                            return 'info:embase/'+ self.Identifiers['EMBASE']
                        except KeyError:
                            try:
                                return 'info:pii/'+ self.Identifiers['ELSEVIER']
                            except KeyError:
                                try:
                                    return 'info:pmc/'+ self.Identifiers['PMC']
                                except KeyError: return ''


    def _make_non_standard_textref(self):
        try:
            return 'info:nctid/'+ self.Identifiers['NCT ID']
        except KeyError:
            try:
                return self.Identifiers['TextRef']
            except KeyError:
                try:
                    return 'info:nihreporter/'+ self.Identifiers['REPORTER']
                except KeyError:
                    try:
                        return 'info:nihgrant/'+ self.Identifiers['GRANTNUMREPORTER']
                    except KeyError:
                        try:
                            return 'info:patentapp/'+ self.Identifiers[PATENT_APP_NUM]
                        except KeyError:
                            try:
                                return 'info:patentgrant/'+ self.Identifiers[PATENT_GRANT_NUM]
                            except KeyError:
                                return NotImplemented


    def _make_textref(self):
        textref = self._make_standard_textref()
        return textref if textref else self._make_non_standard_textref()


    def set_weight(self, weight:float):
        try:
            w = self['weight']
            if w < weight:
                self['weight'] = weight
        except KeyError:
            self['weight'] = weight



#########################################DocMine#############################################DocMine##################

INSTITUTION_KEYWORDS = {'institute', 'institut', 'clinic', 'hospital', 'university', 'universitat', 'universiti', 'centre', 'center', 'inc', 'colleges',
                        'ltd','gmbh', 'school','politecnic','politecnico', 'college', 'department', 'division', 'council', 'academy','faculty', 'co', 'corp',
                        'laboratory', 'labs', 'biolabs', 'biolab', 'lab'}



class DocMine (Reference):
    def __init__(self, doc_id_type, doc_id):
        super().__init__(doc_id_type, doc_id)
        self.sections = dict() #{{abstract:[text]}, {claims:[claims]}. Texts to be markerd up by NLP.

    @staticmethod
    def __textref_suffix(section_name:str):
        name2suffix = {TITLE:'title', ABSTRACT:'abs',CLAIMS:'claims'}
        try:
            return name2suffix[section_name]
        except KeyError:
            return 'cont'

    def get_title(self): return self['Title'][0]

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
                        #self.snippets[text_ref] = {SENTENCE: [sentence_markup]}
                        for msid_range, concept_dict in range2dict.items():
                            prop_name = medscan.get_concept_type(msid_range)
                            self.add_sentence_props(text_ref,prop_name,list(concept_dict.values()))
                            #self.snippets[text_ref]= {prop_name: list(concept_dict.values())}

                    sentence_idx +=1

    def get_annotations(self, prop_name:str):
        try:
            return list(self[prop_name])
        except KeyError:
            if prop_name == INSTITUTIONS:
                return list(self.addresses.keys())
            else: 
                return []

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
        name_words = name.split(' ')
        for w in name_words:
            if w.lower() in INSTITUTION_KEYWORDS: return True
        return False

    @staticmethod
    def execution_time(execution_start):
        return "{}".format(str(timedelta(seconds=time.time() - execution_start)))



            
            


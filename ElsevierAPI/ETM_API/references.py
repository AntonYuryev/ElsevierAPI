#from ElsevierAPI.ResnetAPI.NetworkxObjects import PSObject
from textblob import TextBlob
from ElsevierAPI.ETM_API.medscan import MedScan
import xlsxwriter
import re

AUTHORS = 'Authors'
INSTITUTIONS = 'Institutions'
JOURNAL = 'Journal'
PUBYEAR = 'PubYear'
SENTENCE = 'Sentence'
EMAIL = re.compile(r"\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}\b", flags=re.IGNORECASE)

PS_ID_TYPES = {'PMID', 'DOI', 'PII', 'PUI', 'EMBASE','NCT ID'}
ETM_ID_TYPES = {'ELSEVIER','PMC','REPORTER','GRANTNUMREPORTER'}
OTHER_ID_TYPES = {'Patent ID'}
CLINTRIAL_PROPS = {'TrialStatus','Phase','StudyType','Start','Intervention','Condition','Company','Collaborator'}
PS_BIBLIO_PROPS = {PUBYEAR,AUTHORS,JOURNAL,'MedlineTA','Title','PubMonth','PubDay'}
BIBLIO_PROPS = set(PS_BIBLIO_PROPS)
BIBLIO_PROPS.add(INSTITUTIONS)
REF_ID_TYPES = PS_ID_TYPES | ETM_ID_TYPES |OTHER_ID_TYPES

SENTENCE_PROPS = {SENTENCE,'Organism','CellType','CellLineName','Organ','Tissue','Source','Percent'}
#also TextRef - used as key in Reference.Sentences
RELATION_PROPS = {'Effect','Mechanism','ChangeType','BiomarkerType','QuantitativeType'}

NOT_ALLOWED_IN_SENTENCE='[\t\r\n\v\f]' # regex to clean up special characters in sentences

class Reference(dict):  
# self{BIBLIO_PROPS[i]:value}; Identifiers{REF_ID_TYPES[i]:identifier}; Sentences{TextRef:{SENTENCE_PROPS[i]:Value}}
    pass

    def __init__(self, idType:str, ID:str):
        super().__init__(dict()) # self{BIBLIO_PROPS[i]:value};
        self.Identifiers = {idType:ID} #from REF_ID_TYPES
        self.snippets = dict() # {TextRef:{PropID:[Values]}} from SENTENCE_PROPS. contains senetences marked up by NLP
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

    def update_with_value(self, PropId, PropValue: str):
        try:
            self[PropId] = list(set(self[PropId]) | {PropValue})
        except KeyError:
            self[PropId] = [PropValue]

    def update_with_list(self, PropId, PropValues: list):
        try:
            self[PropId] = list(set(self[PropId] + PropValues))
        except KeyError:
            self[PropId] = list(set(PropValues))

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

    def to_str(self, id_types: list=None, col_sep='\t'):
        id_types = id_types if isinstance(id_types,list) else ['PMID']
        row = str()
        for t in id_types:
            try:
                row = row + t+':'+self.Identifiers[t]+col_sep
            except KeyError:
                row = row + col_sep

        row = row[:len(row)-1] #remove last separator
        for prop_id, prop_values in self.items():
            prop_val = ';'.join(prop_values)
            if prop_id in ['Title', 'Abstract']:
                prop_val = re.sub(NOT_ALLOWED_IN_SENTENCE, ' ', prop_val)
            row = row + col_sep + prop_id + ':' + prop_val

        for text_ref, prop in self.snippets.items():
            for prop_id, prop_values in dict(prop).items():
                prop_value = '.'.join(prop_values)
                if prop_id == 'Sentence':
                    prop_value = re.sub(NOT_ALLOWED_IN_SENTENCE, ' ', prop_value)
                row = row + col_sep + prop_id + ' ('+text_ref+')' + ':' + prop_value
        return row


    def _merge(self, other):
        if isinstance(other, Reference):
            self.update(other)
            self.Identifiers.update(other.Identifiers)
            self.snippets.update(other.snippets)

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
                            return 'info:patent/'+ self.Identifiers['PATENT ID']
                        except KeyError:
                            return NotImplemented


    def _make_textref(self):
        textref = self._make_standard_textref()
        return textref if textref else self._make_non_standard_textref()



INSTITUTION_KEYWORDS = {'institute', 'institut', 'clinic', 'hospital', 'university', 'universitat', 'universiti', 'centre', 'center', 'inc', 'colleges',
                        'ltd','gmbh', 'school','politecnic','politecnico', 'college', 'department', 'division', 'council', 'academy','faculty', 'co', 'corp',
                        'laboratory', 'labs', 'biolabs', 'biolab', 'lab'}



class DocMine (Reference):
    def __init__(self, doc_id_type, doc_id):
        super().__init__(doc_id_type, doc_id)
        self.sections = dict() #{{abstract:[sentences]}, {claims:[sentences]}. Sentences to be markerd up by NLP

    def get_title(self): return self['Title']
    def add2section(self,section_name:str, sentence_list:str):
        for s in sentence_list:
            sentences = list(map(str,TextBlob(s).sentences))
            try:
                section_sentences = self.sections[section_name] + sentences
            except KeyError: 
                section_sentences = sentences

            self.sections[section_name] = section_sentences


    def set_date (self, year, month='', day=''): 
        self['PubYear'] = [year]
        if month: self['PubMonth'] = [month]
        if day: self['PubDay'] = [day]


    def medscan_annotate (self,id_ranges:list, org_name:str, license:str):
        medscan = MedScan(license,org_name)
        for secname, sentences in self.sections.items():
            sec_annotation = dict()
            for s in sentences:
                sec_annotation.update(medscan.find_concepts(s,id_ranges))
            self.section_annotations[secname] = sec_annotation

    def count_str(self,counter:dict, prop_name:str):
        try:
            prop_values = list(self[prop_name])
        except KeyError:
            if prop_name == INSTITUTIONS:
                prop_values = [k for k,v in self.addresses.items()]
            else: return

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

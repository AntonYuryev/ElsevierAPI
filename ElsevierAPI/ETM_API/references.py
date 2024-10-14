from builtins import len
from .medscan import MedScan
import xlsxwriter,re,time,json,unicodedata
from datetime import timedelta
from ..NCBI.pubmed import pubmed_hyperlink
from titlecase import titlecase

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
EMAIL = re.compile(r"\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}\b", flags=re.IGNORECASE)

INT_PROPS = {RELEVANCE,PS_CITATION_INDEX,PUBYEAR,ETM_CITATION_INDEX,SBS_CITATION_INDEX}

ARTICLE_ID_TYPES = ['PMID', 'DOI', 'PII', 'PUI', 'EMBASE']
PS_REFIID_TYPES = ARTICLE_ID_TYPES + ['NCT ID','NCTID']
#keep PS_ID_TYPES as list for efficient identifier sort.  ID types are ordered by frequency in Resnet
ETM_ID_TYPES = ['ELSEVIER','PMC','REPORTER','GRANTNUMREPORTER']
PATENT_ID_TYPES = [PATENT_APP_NUM, PATENT_GRANT_NUM]
CLINTRIAL_PROPS = {'TrialStatus','Phase','StudyType','Start','Intervention','Condition','Company','Collaborator'}

JOURNAL_PROPS = {JOURNAL,'ISSN','ESSN',MEDLINETA}
PS_BIBLIO_PROPS_ALL = {PUBYEAR,AUTHORS,TITLE,'PubMonth','PubDay','PubTypes','Start'}|JOURNAL_PROPS # 'Start' = 'PubYear'
PS_BIBLIO_PROPS = {TITLE, PUBYEAR, AUTHORS, JOURNAL, MEDLINETA,'Start'} # contains only props necessary for reference::__biblio_tuple

BIBLIO_PROPS = PS_BIBLIO_PROPS | {INSTITUTIONS,RELEVANCE}
REF_ID_TYPES = PS_REFIID_TYPES+ETM_ID_TYPES+PATENT_ID_TYPES

PS_SENTENCE_PROPS = [SENTENCE,'Organism','CellType','CellLineName','Organ','Tissue','Source','Percent',THRESHOLD,'pX','Phase','Start','TrialStatus','URL','Experimental System']
SENTENCE_PROPS = PS_SENTENCE_PROPS + ['Evidence','msrc','mref','Similarity']
# SENTENCE_PROPS needs to be a list for ordered printing
#also TextRef - used as key in Reference.Sentences
RELATION_PROPS = {EFFECT,'Mechanism','ChangeType','BiomarkerType','QuantitativeType'}

PS_REFERENCE_PROPS = list(CLINTRIAL_PROPS)+PS_REFIID_TYPES+list(PS_BIBLIO_PROPS_ALL)+PS_SENTENCE_PROPS+['TextRef']
ALL_PSREL_PROPS = list(RELATION_PROPS)+PS_REFERENCE_PROPS

REFERENCE_PROPS = list(BIBLIO_PROPS)+list(CLINTRIAL_PROPS)+REF_ID_TYPES+SENTENCE_PROPS

NOT_ALLOWED_IN_SENTENCE='[\t\r\n\v\f]' # regex to clean up special characters in sentences, titles, abstracts

def make_hyperlink(identifier:str,url:str,display_str=''):
        # hyperlink in Excel does not work with long URLs, 
        display_str = display_str if display_str else identifier
        return '=HYPERLINK("'+url+identifier+'",\"{}\")'.format(display_str)


class Reference(dict):
    '''
    Reference{BIBLIO_PROPS[i]:[values]};\n
    Reference.Identifiers{REF_ID_TYPES[i]:identifier};\n
    Reference.Sentences{TextRef:{SENTENCE_PROPS[i]:Value}}\n
    '''
    pass

    def __init__(self, idType:str, ID:str):
        super().__init__(dict({})) # self{BIBLIO_PROPS[i]:[values]};
        self.Identifiers = {idType:ID} #from REF_ID_TYPES
        self.snippets = dict() # {TextRef:{PropID:[Values]}} PropID is from SENTENCE_PROPS contains sentences marked up by NLP           

    def copy_ref(self):
        '''
        Return
        ------
        Reference object copy of self if self has valid identifiers\n
        otherwise returns dict()
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
        return dict({})
    

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
        #if reference has non-canonical TextRef that cannot be parsed by __parse_textref
        # self.Identifiers has ['TextRef'] value
        except KeyError: return NotImplemented


    def __hash__(self):
        #__hash__ needs __eq__ to work properly
        return hash(self.__key())
    

    def __eq__(self, other):
        #__hash__ needs __eq__ to work properly
        if isinstance(other, Reference):
            for id_type in REF_ID_TYPES:
                try:
                    return self.Identifiers[id_type] == other.Identifiers[id_type]
                except KeyError:
                    continue
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


    def update_with_list(self, prop_id, with_values:list):
        clean_vals = set(map(lambda x: str(x).strip(' .'),with_values))
        try:
            self[prop_id] = list(set(self[prop_id])|clean_vals)
        except KeyError:
            self[prop_id] = list(clean_vals)


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


    def has_property(self, prop_name:str):
        # prop_names2values = {prop_name:[values]}
        for prop2values in self.snippets.values():
            if prop_name in prop2values.keys():
                return True
        return prop_name in self.keys()
    

    def get_values(self,prop_name:str):
        for prop2values in self.snippets.values():
            try:
                return list(prop2values[prop_name])
            except KeyError: continue
        try:
            return list(self[prop_name])
        except KeyError: 
            return []

    
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
            for prop, values in my_prop2values.items():
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

        return False
    

    def remove_props(self, prop_names:list):
        my_copy = self.copy_ref()
        if isinstance(my_copy,Reference): #copy_ref
            for prop2values in my_copy.snippets.values():
                [prop2values.pop(p,'') for p in prop_names]
            
            [my_copy.pop(p,'') for p in prop_names]
            [my_copy.Identifiers.pop(p,'') for p in prop_names]
            return my_copy
        return dict({})

    
    def get_doc_id(self):
        '''
        Return
        ------
        tuple(id_type, identifier) for the first id type from REF_ID_TYPES\n
        tuple('','') if reference does not have id_type in REF_ID_TYPES
        '''
        for id_type in REF_ID_TYPES:
            try:
                return id_type, self.Identifiers[id_type]
            except KeyError: continue
        return str(),str()
    

    def _identifiers_str(self):
        '''
        Return
        ------
        id_type:id_value
        '''
        id_type, identifier = self.get_doc_id()
        return id_type  +':'+identifier
    

    def identifier(self,identifier_type):
        try:
            return self.Identifiers[identifier_type]
        except KeyError:
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
                            identifier = make_hyperlink(identifier,'http://dx.doi.org/')
                        elif t == PATENT_APP_NUM:
                            identifier = make_hyperlink(identifier,'https://patents.google.com/patent/')
                    row.append(identifier)
                except KeyError:
                    row.append('')
        else:
            row.append(self._identifiers_str())
                
        for prop_id in biblio_props:
            try:
                prop_values_str = ';'.join(list(map(str,self[prop_id])))
                if prop_id in ['Title', 'Abstract']:
                    prop_values_str = re.sub(NOT_ALLOWED_IN_SENTENCE,' ',prop_values_str)
            except KeyError:
                prop_values_str = ''
            row.append(prop_values_str)

        if print_snippets:
            sentence_props = json.dumps(self.snippets)
            sentence_props = re.sub(NOT_ALLOWED_IN_SENTENCE,' ',sentence_props)
            row.append(sentence_props)

        return row


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
        try:
            return int(self[PUBYEAR][0])
        except KeyError:
            try:
                year = str(self['Start'][0]) # Clinical trials case
                if year[-4:].isdigit(): # format: 20-Apr-2020; August 2004 
                    return int(year[-4:])
                elif year[-2:].isdigit():
                    return int('20'+year[-2:]) # format: 20-Apr-20
                elif year[:2].isdigit(): 
                    return int('20'+year[:2]) # format 20-May
                else:
                    print(f'Unknown Clinical trial "Start" format: {year}')
                    return 1812
            except KeyError: return 1812 # No PubYear case


    def title(self): 
        try:
            return self['Title'][0]
        except KeyError:
            return ''
        
    
    def pmid(self):
        try:
            return self.Identifiers['PMID']
        except KeyError:
            return ''


    def doi(self):
        try:
            return self.Identifiers['DOI']
        except KeyError:
            return ''


    def number_of_snippets(self):
        return len(self.snippets)


    def textrefs(self):
        return list(self.snippets.keys())


    def journal(self): 
        '''
        Return
        ------
        self[JOURNAL][0] 
        '''
        try:
            return self[JOURNAL][0]
        except KeyError:
            return str('No journal name')
                  
        
    def relevance(self):
        try:
            return float(self[RELEVANCE][0])
        except KeyError:
            return float(0.0)
    

    def is_clinical_trial(self):
        try:
            return self.Identifiers['NCT ID']
        except KeyError:
            return ''
        

    def author_list(self)->list[str]:
        '''
        Return
        ------
        sorted list of authors from self[AUTHORS]
        '''
        try:
            authors2split = self[AUTHORS]
            individual_authors = set()
            for au_list in authors2split:
                authors = au_list.split(';')
                individual_authors.update(authors)
            
            individual_authors = list(individual_authors)
            individual_authors.sort()
            self[AUTHORS] = individual_authors
            return individual_authors
        except KeyError:
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
            authors_list = self[AUTHORS][:3]
            authors_str = ','.join(authors_list) if authors_list else 'unknown authors.'
            if len(self[AUTHORS]) > 3:
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


    def _merge(self, other):
        if isinstance(other, Reference):
            for prop, values in other.items():
                if prop in INT_PROPS:
                    clean_vals = set(map(int,values))
                else:
                    clean_vals = {str(v).strip(' .') for v in values if v is not None}

                self.update_with_list(prop,list(clean_vals))

            self.Identifiers.update(other.Identifiers)

            for textref, prop2val in other.snippets.items():
                try:
                    exist_textref_props = self.snippets[textref]
                    if exist_textref_props == prop2val:
                        continue
                    else:
                        new_text_ref = textref+':'+str(len(self.snippets))
                        [self.add_sentence_props(new_text_ref,p,v) for p,v in prop2val.items()]
                except KeyError:
                    [self.add_sentence_props(textref,p,v) for p,v in prop2val.items()]
        return


    def is_from_abstract(self):
        for textref in self.snippets.keys():
            try:
                return bool(str(textref).rindex('#abs',-8,-3))
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
                                try:
                                    return 'info:title/'+self.Identifiers[TITLE] #.replace(' ','%20')
                                except KeyError:
                                    return 'NotImplemented'


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
            return 0.0 # allows to exclude cerain references from counting


    def _sort_key(self, by_property, is_numerical=True):
        if by_property == PUBYEAR: return self.pubyear()
        else:
            try:
                return float(self[by_property][0]) if is_numerical else str(self[by_property][0])
            except KeyError:
                return 0.0 if is_numerical else '0'


    def get_snippet_prop(self,prop_name:str):
        """
        Returns
        -------
        {textref:[prop_vals]}
        """
        prop_values = dict()
        for textref, sentence_props in self.snippets.items():
            try:
                prop_vals = sentence_props[prop_name]
                prop_values[textref] = prop_vals
            except KeyError:
                continue

        return prop_values
    

    def get_props(self,prop_name:str)->list[str]:
        '''
        input:
            prop_name can be either in self.snippet or self.Identifiers or self
        '''
        try:
            return list(self[prop_name])
        except KeyError:
            try:
                return [self.Identifiers[prop_name]]
            except KeyError:
                snippet_prop_dic = self.get_snippet_prop(prop_name)
                props = set()
                for text_ref, prop_vals in snippet_prop_dic.items():
                    props.update(prop_vals)
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
            except ValueError:
                return if_missing_return
        return if_missing_return


class SerializableRef:
    '''
    Enables Reference serialization into json dump because
    classes derived from dict and list cannot have custom JSONEncoder
    '''

    def __init__(self,ref:Reference):
        self.reference = ref


class ReferenceEncoder(json.JSONEncoder):
    def default(self, sref):
        if isinstance(sref, SerializableRef):
            # Convert CustomObject to a dictionary
            ref = sref.reference
            dump_dict = dict(ref)
            dump_dict.update({'Identifiers':ref.Identifiers})
            if ref.snippets:
                dump_dict.update({'snippets':ref.snippets})
            return {'reference':dump_dict}
        else:
            return json.JSONEncoder.default(self, sref)


class ReferenceDecoder(json.JSONDecoder):
    def __init__(self, *args, **kwargs):
        json.JSONDecoder.__init__(self, object_hook=self.object_hook, *args, **kwargs)

    def object_hook(self, dct):
        my_dict = dct['reference']
        Identifiers = my_dict.pop('Identifiers')
        snippets = my_dict.pop('snippets',dict())
        addresses = my_dict.pop('addresses',dict())

        id_type,id = next(iter(my_dict.items()))
        ref = Reference(id_type,id)
        ref.update(my_dict)
        ref.Identifiers.update(Identifiers)
        ref.snippets.update(snippets)
        return ref

#########################################DocMine#############################################DocMine##################

INSTITUTION_KEYWORDS = {'institute', 'institut', 'instituto', 'istituto','clinic', 'klinik','genetics','hospital', 'university', 'universitat', 'universiti', 'université','università','universidade','universidad','universitï¿½','universitaria','centre', 'center', 'centro', 'inc', 'colleges',
                        'ltd','gmbh', 'llc', 'school','politecnic','politecnico', 'college', 'department', 'departamento','division', 'council', 'academy','faculty', 'co', 'corp','ministry','campus','group','corporation', 'consulting', 'pharma','union',
                        'laboratory', 'laboratoire', 'laboratories','society','project', 'association', 'commission','program','drug', 'labs', 'biolabs', 'biolab', 'lab', 'unit', 'pharmaceutical', 'policlinico', 'cátedra', 'sciences', 'r&d', 'research', 'trust', 'fund', 'plant', 'biomed','salud', 'health','foundation','federation','service','cluster','fondazione','bio','pharmaceuticals', 'farmaceutico'}

INSTITUTION_KEYWORDS = {unicodedata.normalize('NFKD', word).casefold() for word in INSTITUTION_KEYWORDS} 

def removeThe(t:str):
    return t[4:] if t.startswith('The ') else t


class Author:
    LastName=str()
    _1stName=str()
    organization = list()
    address = str()
    email = str()

    def __init__(self,LastName,_1stName='',organization=[],address=str(),email=''):
        self.LastName = LastName
        self._1stName = _1stName
        self.organization = organization
        self.address = address
        self.email = email

    def institution(self):
        return self.organization[-1] if self.organization else []


class DocMine (Reference):
    '''
    use this class to normalize annotations from references in disparate sources, e.g. ETM,USPO,EPO,Pubmed\n
    DocMine stores text in sections: Abstract, Results,Discussion, Claims, Descriptions
    '''
    def __init__(self, doc_id_type, doc_id):
        super().__init__(doc_id_type, doc_id)
        self.sections = dict() #{{abstract:[text]}, {claims:[claims]}. Texts to be markerd up by NLP.
        self.authors = list() # [Author]
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
                        #self.snippets[text_ref] = {SENTENCE: [sentence_markup]}
                        for msid_range, concept_dict in range2dict.items():
                            prop_name = medscan.get_concept_type(msid_range)
                            self.add_sentence_props(text_ref,prop_name,list(concept_dict.values()))
                            #self.snippets[text_ref]= {prop_name: list(concept_dict.values())}

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
        return {a.organization for a in self.authors for i in a.organization}
    
    
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



            
            


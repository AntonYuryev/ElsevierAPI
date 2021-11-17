import re
import json
import itertools

class PSObject(dict):  # {PropId:[values], PropName:[values]}
    def __init__(self, dic: dict):
        super().__init__(dic)

    @staticmethod
    def from_zeep_objects(ZeepObjectRef):
        return_dict = dict()
        zeep_iter = iter(ZeepObjectRef)
        while True:
            try:
                item = next(zeep_iter)
            except StopIteration:
                break  # Iterator exhausted: stop the loop
            else:
                return_dict[item] = [ZeepObjectRef[item]]
        return return_dict

    @classmethod
    def from_zeep(cls, ZeepObjectRef):
        return cls(cls.from_zeep_objects(ZeepObjectRef))

    def __hash__(self):
        return self['Id'][0]

    def set_property(self, PropId, PropValue: str):
        self[PropId] = [PropValue]

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


    def data2str(self, columnPropNames: list, col_sep='\t', cell_sep=';', endOfline='\n'):
        # assumes all properties in columnPropNames were fetched from Database otherwise will crash
        table_row = str()
        for propName in columnPropNames:
            try:
                values = self[propName]
                prop_val = cell_sep.join(values)
            except KeyError:
                prop_val = ''

            table_row = table_row + prop_val + col_sep

        return table_row[0:len(table_row) - 1] + endOfline


    def has_property(self,prop_name,prop_values:list,case_sensitive=False):
        try:
            if case_sensitive:
                search_in = self[prop_name]
                search_set = set(prop_values)   
            else:
                search_in = list(self[prop_name])
                search_in  = list(map(lambda x: x.lower(),search_in))
                search_set = set(map(lambda x: x.lower(),prop_values))
                
            return not search_set.isdisjoint(search_in)
        except KeyError: return False


REF_ID_TYPES = {'PMID', 'DOI', 'PII', 'PUI', 'EMBASE', 'NCT ID'}
CLINTRIAL_PROPS = {'TrialStatus','Phase','StudyType','Start','Intervention','Condition','Company','Collaborator'}
BIBLIO_PROPS = {'PubYear','Authors','Journal','MedlineTA','Title'}

SENTENCE_PROPS = {'Sentence','Organism','CellType','CellLineName','Organ','Tissue','Source','Percent'}
#also TextRef - used as key in Reference.Sentences
RELATION_PROPS = {'Effect','Mechanism','ChangeType','BiomarkerType','QuantitativeType'}

NOT_ALLOWED_IN_SENTENCE='[\t\r\n\v\f]' # regex to clean up special characters in sentences


class Reference(PSObject):  
# self{BIBLIO_PROPS[i]:value}; Identifiers{REF_ID_TYPES[i]:identifier}; Sentences{TextRef:{SENTENCE_PROPS[i]:Value}}
    pass

    def __init__(self, idType:str, ID:str):
        super().__init__(dict()) # self{BIBLIO_PROPS[i]:value};
        self.Identifiers = {idType:ID} #from REF_ID_TYPES
        self.Sentences = dict() # {TextRef:{PropID:Value}} from SENTENCE_PROPS


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
            if prop_id == 'Title':
                prop_val = re.sub(NOT_ALLOWED_IN_SENTENCE, ' ', prop_val)
            row = row + col_sep + prop_id + ':' + prop_val

        for text_ref, prop in self.Sentences.items():
            for prop_id, prop_value in dict(prop).items():
                if prop_id == 'Sentence':
                    prop_value = re.sub(NOT_ALLOWED_IN_SENTENCE, ' ', prop_value)
                row = row + col_sep + prop_id + ' ('+text_ref+')' + ':' + prop_value
        return row


    def _merge(self, other):
        if isinstance(other, Reference):
            self.update(other)
            self.Identifiers.update(other.Identifiers)
            self.Sentences.update(other.Sentences)

    def is_from_abstract(self):
        for textref in self.Sentences.keys():
            try:
                return bool(str(textref).rindex('#abs', len(textref) - 8, len(textref) - 3))
            except ValueError:
                continue
        return False

    def _make_textref(self):
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
                                return 'info:nctid/'+ self.Identifiers['NCT ID']
                            except KeyError:
                                try:
                                    return self.Identifiers['TextRef']
                                except KeyError:
                                    return NotImplemented


class PSRelation(PSObject):
    pass
    PropSetToProps = dict()
    Nodes = dict()
    References = dict()

    def __init__(self, dic:dict):
        super().__init__(dic)
        self.PropSetToProps = dict()  # {PropSetID:{PropID:[values]}}
        self.Nodes = dict()  # {"Regulators':[(entityID, Dir, effect)], "Targets':[(entityID, Dir, effect)]}
        self.References = dict()  # {refIdentifier (PMID, DOI, EMBASE, PUI, LUI, Title):Reference}

    def copy(self, other):
        if isinstance(other, PSRelation):
            self.update(other)
            self.PropSetToProps.update(other.PropSetToProps)
            self.Nodes.update(other.Nodes)
            self.References.update(other.References)

    def __hash__(self):
        return self['Id'][0]


    def __props2str(self, prop_id, cell_sep: str =';'):
        try:
            return cell_sep.join(map(str, self[prop_id]))
        except KeyError:
            prop_set_values = []
            for prop in self.PropSetToProps.values():
                try:
                    prop_set_values.append(cell_sep.join(map(str, prop[prop_id])))
                except KeyError:
                    continue
            
            to_return = cell_sep.join(prop_set_values)
            if prop_id in ['Sentence','Title']:
                to_return = re.sub(NOT_ALLOWED_IN_SENTENCE, ' ', to_return)
            return to_return

    def props2list(self, propID, cell_sep=';'):
        try:
            return self[propID]
        except KeyError:
            to_return = list()
            for prop in self.PropSetToProps.values():
                prop_set_val = str()
                for prop_id, values in prop.items():
                    if prop_id == propID:
                        prop_set_val = cell_sep.join(values)
                        break
                to_return.append(prop_set_val)
            return to_return


    def load_references(self):
        if self.References: return
        for propSet in self.PropSetToProps.values():
            my_reference_tuples = list()
            for ref_id_type in REF_ID_TYPES:
                try:
                    ref_id = propSet[ref_id_type][0]
                    my_reference_tuples.append((ref_id_type, ref_id))
                except KeyError:
                    continue

            if my_reference_tuples: # propSet is valid reference - resolving duplicates 
                # case when reference have valid identifiers
                article_identifiers = {x[1] for x in my_reference_tuples}
                existing_ref = [r for i, r in self.References.items() if i in article_identifiers]
                if not existing_ref:  # case when reference is new
                    id_type, id_value = my_reference_tuples[0][0], my_reference_tuples[0][1]
                    ref = Reference(id_type, id_value)
                    self.References[id_value] = ref
                    for Id in range(1, len(my_reference_tuples)):
                        id_type, id_value = my_reference_tuples[Id][0], my_reference_tuples[Id][1]
                        ref.Identifiers[id_type] = id_value
                        self.References[id_value] = ref
                else: 
                    ref = existing_ref[0]
                    if len(existing_ref) > 1:  # identifiers from one propset point to different references
                    # will merge all references from propset with the first one
                        for i in range(1, len(existing_ref)):
                            ref._merge(existing_ref[i])
                            for id_value in dict(existing_ref[i].Identifiers).values():
                                self.References[id_value] = ref
                        del existing_ref[1:]

            else: #propSet is not valid reference - trying to create one using Title or TexRef
                try:
                    # trying id reference by title as a last resort since it does not have valid identifiers
                    propset_title = propSet['Title'][0]
                    try:
                        ref = self.References[propset_title]
                    except KeyError:
                        ref = Reference('Title', propset_title)
                        self.References[propset_title] = ref
                except KeyError:
                    try:
                        txtref = propSet['TextRef'][0]
                        if txtref != 'Admin imported':
                            ref = Reference.from_textref(txtref)
                            for ref_id_type, ref_id in ref.Identifiers.items():
                                self.References[ref_id] = ref
                        else: continue #'Admin imported' references have no identifiers and ignored 
                    except KeyError: continue

            textref = None 
            sentence_props = dict()
            for propId, propValues in propSet.items():  # adding all other valid properties to Ref
                if propId in BIBLIO_PROPS:
                    ref.update_with_list(propId, propValues)
                elif propId in SENTENCE_PROPS:
                    sentence_props[propId] = propValues
                elif propId in CLINTRIAL_PROPS:
                    ref.update_with_list(propId, propValues)
                elif propId == 'TextRef': textref = propValues[0]

            if not isinstance(textref,str): 
                textref = ref._make_textref()
            elif textref[4:10] == 'hash::': # references from older Resnet versions can start with 'urn:hash::'
                textref = ref._make_textref()
            elif textref == 'Admin imported': 
                continue # ignore and move to the next PropSet
            
            if sentence_props: 
                ref.Sentences[textref] = sentence_props
            else: 
                ref.Sentences[textref] = {'Sentence':[]} #load empty dict for data consistency

        self.PropSetToProps.clear() # PropSetToProps are not needed anymore. Everything is loaded into References.


    def get_reference_count(self, count_abstracts=False):
        if count_abstracts:
            ref_from_abstract = set([x for x in self.References.values() if x.is_from_abstract()])
            return len(ref_from_abstract)
        else:        
            return len(self.References) if self.References else self['RelationNumberOfReferences'][0]



    def to_table_dict(self, columnPropNames:list, cell_sep:str=';', RefNumPrintLimit=0, add_entities=False):
        # assumes all properties in columnPropNames were fetched from Database otherwise will crash
        # initializing table
        col_count = len(columnPropNames) +2 if add_entities else len(columnPropNames)
        RelationNumberOfReferences = int(self['RelationNumberOfReferences'][0])

        rowCount = 1
        if RelationNumberOfReferences >= RefNumPrintLimit:
            rowCount = max(1, len(self.PropSetToProps))

        table = dict()
        for r in range(0, rowCount):
            table[r] = [''] * col_count

        if add_entities:
            regulatorIDs = str()
            targetIDs = str()
            for k, v in self.Nodes.items():
                if k == 'Regulators':
                    regulatorIDs = ','.join(map(str, [x[0] for x in v]))
                else:
                    targetIDs = ','.join(map(str, [x[0] for x in v]))

            for r in range(0, rowCount):
                table[r][col_count-2] = regulatorIDs
                table[r][col_count-1] = targetIDs

        for col in range(len(columnPropNames)):
            propId = columnPropNames[col]
            if propId in self.keys(): # filling columns with relation properties
                for row in range(0, rowCount):
                    propValue = self.__props2str(propId)
                    table[row][col] = propValue
            elif RelationNumberOfReferences >= RefNumPrintLimit: #filling columns with reference properties
                row = 0
                for propList in self.PropSetToProps.values():
                    if propId in propList.keys():
                        propValues = propList[propId]
                        cellValue = cell_sep.join(propValues)
                        table[row][col] = cellValue
                    row += 1

        return table

    def to_table_str(self, columnPropNames:list, col_sep:str='\t', 
                     cell_sep=';', RefNumPrintLimit=0, add_entities=False):

        table_dict = self.to_table_dict(columnPropNames,cell_sep,RefNumPrintLimit,add_entities)

        tableStr = str()
        for row in table_dict.values():
            tableStr = tableStr + col_sep.join(row) + '\n'
        return tableStr

    def to1row(self, columnPropNames:list, col_sep:str='\t', cell_sep:str=';', RefNumPrintLimit=0, add_entities=False):
        # assumes all properties in columnPropNames were fetched from Database otherwise will crash
        # initializing table
        col_count = len(columnPropNames) +2 if add_entities else len(columnPropNames)
        RelationNumberOfReferences = int(self['RelationNumberOfReferences'][0])

        table = ['']*col_count
        if add_entities:
            regulatorIDs = str()
            targetIDs = str()
            for k, v in self.Nodes.items():
                if k == 'Regulators':
                    regulatorIDs = ','.join(map(str, [x[0] for x in v]))
                else:
                    targetIDs = ','.join(map(str, [x[0] for x in v]))

                table[col_count-2] = regulatorIDs
                table[col_count-1] = targetIDs

        for col in range(len(columnPropNames)):
            propId = columnPropNames[col]
            if propId in self.keys(): # filling columns with relation properties
                propValues = self.__props2str(propId)
                table[col] = propValues.replace(cell_sep, ' ')
            elif RelationNumberOfReferences >= RefNumPrintLimit: #filling columns with reference properties
                for propList in self.PropSetToProps.values():
                    try:
                        propValues = propList[propId]
                        no_sep_values = [str(p).replace(cell_sep, ' ') for p in propValues]
                        add2cell = ' '.join(no_sep_values)
                        table[col] += add2cell + cell_sep
                    except KeyError: table[col] += cell_sep

        for c in range(0,len(table)):
            col = table[c]
            table[c] = col.strip(cell_sep)

        return col_sep.join(table)+'\n'

    def triple2str(self,columnPropNames:list, col_sep:str='\t', cell_sep:str=';', RefNumPrintLimit=0, add_entities=False, as1row=False):
        if as1row:
            return self.to1row(columnPropNames,col_sep,cell_sep,RefNumPrintLimit,add_entities)
        else: 
            return  self.to_table_str(columnPropNames,col_sep,cell_sep,RefNumPrintLimit,add_entities)

    def get_regulators_targets(self):
        if len(self.Nodes) > 1:
            return [(r[0],t[0]) for r in self.Nodes['Regulators'] for t in self.Nodes['Targets']]
        else:
            try:
                objIdList = [x[0] for x in self.Nodes['Regulators']]
                return itertools.combinations(objIdList, 2)
            except KeyError:
                objIdList = [x[0] for x in self.Nodes['Targets']]
                return itertools.combinations(objIdList, 2)

    def get_entities_ids(self):
        nodeIds = {[x[0] for x in self.Nodes['Regulators']] + [x[0] for x in self.Nodes['Targets']]}
        return list(nodeIds)

    def to_json(self):
        str1 = '{"Relation Properties": ' + json.dumps(self) + '}'
        strP = '{"Relation References": ' + json.dumps(self.PropSetToProps) + '}'
        strR = json.dumps(self.Nodes)
        return str1 + '\n' + strP + '\n' + strR + '\n'

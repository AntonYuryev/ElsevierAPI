import re
import json
import itertools
from ElsevierAPI.ETM_API.references import JOURNAL,PS_ID_TYPES,NOT_ALLOWED_IN_SENTENCE,PS_BIBLIO_PROPS,SENTENCE_PROPS,CLINTRIAL_PROPS,MEDLINETA


PROTEIN_TYPES = ['Protein','FunctionalClass','Complex']

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
            if case_sensitive or prop_name in {'Id','id','ID'}:
                search_in = self[prop_name]
                search_set = set(prop_values)   
            else:
                search_in = list(self[prop_name])
                search_in  = list(map(lambda x: x.lower(),search_in))
                search_set = set(map(lambda x: x.lower(),prop_values))
                
            return not search_set.isdisjoint(search_in)
        except KeyError: return False

    def has_value_in(self,prop2values:dict,case_sensitive=False):
        # prop2values = {propName:[values]}
        for prop_name, prop_values in prop2values.items():
            if self.has_property(prop_name,prop_values,case_sensitive): return True

        return False
    
from ElsevierAPI.ETM_API.references import Reference
class PSRelation(PSObject):
    pass
    PropSetToProps = dict() # {PropSetID:{PropID:[values]}}
    Nodes = dict() # {"Regulators':[(entityID, Dir, effect)], "Targets':[(entityID, Dir, effect)]}
    References = dict() # {refIdentifier (PMID, DOI, EMBASE, PUI, LUI, Title):Reference}

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
        #if self['Id'][0] == 216172782119007911:
        #    print('')
        for propSet in self.PropSetToProps.values():
            my_reference_tuples = list()
            for ref_id_type in PS_ID_TYPES:
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
                #if propValues[0] == 'info:pmid/16842844':
                #    print('')
                if propId in PS_BIBLIO_PROPS:
                    if propId == MEDLINETA:
                        ref.update_with_list(JOURNAL, propValues)
                    else:
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
                ref.snippets[textref] = sentence_props
            else: 
                ref.snippets[textref] = {'Sentence':[]} #load empty dict for data consistency

        self.PropSetToProps.clear() # PropSetToProps are not needed anymore. Everything is loaded into References.


    def add_references(self, references:list):#list of Reference objects
        self.load_references()
        for r in references:
            was_merged = False
            for id in PS_ID_TYPES:
                try:
                    self.References[r.Identifiers[id]]._merge(r)
                    was_merged = True
                    break
                except KeyError:
                    continue

            if not was_merged:
                for i in r.Identifiers.values():
                    self.References[i] = r


    def sort_references(self, by_property:str):
        all_refs = list(set([x for x in self.References.values() if by_property in x.keys()]))
        all_refs.sort(key=lambda x: x[by_property][0])
        return all_refs


    def filter_references(self, prop_names2values:dict):
        # prop_names2values = {prop_name:[values]}
        all_refs = set(self.References.values())
        ref2keep = {ref for ref in all_refs if ref.has_properties(prop_names2values)}
        self.References = {k:r for k,r in self.References.items() if r in ref2keep}
        return
        #self['RelationNumberOfReferences'] will NOT be changed to keep the original number of references for relation


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
            # for directional relations
            return [(r[0],t[0]) for r in self.Nodes['Regulators'] for t in self.Nodes['Targets']]
        else:
            try:
                # for non-directions relations
                objIdList = [x[0] for x in self.Nodes['Regulators']]
                return itertools.combinations(objIdList, 2)
            except KeyError:
                # for wiered cases
                objIdList = [x[0] for x in self.Nodes['Targets']]
                return itertools.combinations(objIdList, 2)

    def get_regulator_ids(self):
        nodeIds = [x[0] for x in self.Nodes['Regulators']]
        return list(nodeIds)

    def get_target_ids(self):
        try:
            nodeIds = [x[0] for x in self.Nodes['Targets']]
            return list(nodeIds)
        except KeyError: return []

    def get_entities_ids(self):
        return self.get_regulator_ids()+self.get_target_ids()

    def to_json(self):
        str1 = '{"Relation Properties": ' + json.dumps(self) + '}'
        strP = '{"Relation References": ' + json.dumps(self.PropSetToProps) + '}'
        strR = json.dumps(self.Nodes)
        return str1 + '\n' + strP + '\n' + strR + '\n'

    def _weight2ref (self, weight_by_prop_name:str, proval2weight:dict):
        try:
            values2weight = set(self[weight_by_prop_name])
            max_weight = 0.0
            for v in values2weight:
                try:
                    weight = proval2weight[v]
                    if weight > max_weight: max_weight = weight
                except KeyError: continue
            for ref in self.References.values():
                ref.set_weight(max_weight)
        except KeyError:
                for ref in self.References:
                    ref.set_weight(0.0)


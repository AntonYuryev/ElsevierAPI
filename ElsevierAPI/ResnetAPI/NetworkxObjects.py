import re
import json
import itertools
import pickle
import hashlib

from ..ETM_API.references import Reference, len
from ..ETM_API.references import JOURNAL,PS_REFIID_TYPES,NOT_ALLOWED_IN_SENTENCE,BIBLIO_PROPS,SENTENCE_PROPS,CLINTRIAL_PROPS
from ..ETM_API.references import MEDLINETA,EFFECT,PUBYEAR


PROTEIN_TYPES = ['Protein','FunctionalClass','Complex']
REGULATORS = 'Regulators'
TARGETS = 'Targets'
REFCOUNT = 'RelationNumberOfReferences'
CHILDS = 'Child Ids'

#enums for objectypes to avoid misspeling
GENETICVARIANT = 'GeneticVariant'
FUNC_ASSOC = 'FunctionalAssociation'


class PSObject(dict):  # {PropId:[values], PropName:[values]}
    pass
    def __init__(self, dic=dict()):
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

    def name(self):
        try:
            return self['Name'][0]
        except KeyError:
            raise KeyError

    def urn(self):
        try:
            return self['URN'][0]
        except KeyError:
            raise KeyError

    def id(self):
        try:
            return self['Id'][0]
        except KeyError:
            raise KeyError

    def objtype(self):
        try:
            return self['ObjTypeName'][0]
        except KeyError:
            raise KeyError


    def descr(self):
        try:
            return self['Description'][0]
        except KeyError:
            return ''

    
    def get_prop(self,prop_name:str,value_index=-1,missing_value=''):
        try:
            values = self[prop_name]
            if value_index >= 0:
                return values[value_index]
            else:
                return values
        except KeyError:
            return missing_value


    @staticmethod
    def unpack(list_of_lists:list,make_unique=True):
        flat_list = [item for sublist in list_of_lists for item in sublist]
        return list(set(flat_list)) if make_unique else flat_list


    def update_with_value(self, prop_id:str, new_value):
        try:
            values = list(self[prop_id])
            if new_value not in values:
                values.append(new_value)
                self[prop_id] = values
        except KeyError:
            self[prop_id] = [new_value]


    def update_with_list(self, prop_id:str, new_values:list):
        try:
            values = list(self[prop_id])
            [values.append(x) for x in new_values if x not in values]
            self[prop_id] = values
        except KeyError:
            self[prop_id] = new_values


    def merge_obj(self,other:'PSObject'):
        [self.update_with_list(prop_name,values) for prop_name,values in other.items()]


    def _prop2str(self, prop_id:str,cell_sep:str =';'):
        try:
            return cell_sep.join(self[prop_id])
        except KeyError:
            return ''

    def props2dict(self,prop_ids:list):
        return {k:self._prop2str(k) for k in self.keys() if k in prop_ids}

    def data2str(self, columnPropNames: list, col_sep='\t', cell_sep=';', endOfline='\n'):
        table_row = str()
        for propName in columnPropNames:
            try:
                values = self[propName]
                prop_val = cell_sep.join(values)
            except KeyError:
                prop_val = ''

            table_row = table_row + prop_val + col_sep
        return table_row[0:len(table_row) - 1] + endOfline


    def has_property(self, prop_name:str):
        return prop_name in self.keys()

    def has_propeties(self,prop_names:set):
        return not prop_names.isdisjoint(set(self.keys()))

    def is_annotated(self,with_prop:str,having_values:list,case_sensitive=False):
        try:
            if case_sensitive or with_prop in {'Id','id','ID'}:
                search_in = self[with_prop]
                search_set = set(with_prop)   
            else:
                search_in = list(self[with_prop])
                search_in  = set(map(lambda x: x.lower(),search_in))
                search_set = set(map(lambda x: x.lower(),having_values))
                
            return not search_set.isdisjoint(search_in)
        except KeyError: return False


    def has_value_in(self,prop2values:dict,case_sensitive=False):
        """
        Input
        -----
        prop2values = {propName:[values]}
        """
        for prop_name, prop_values in prop2values.items():
            if self.is_annotated(prop_name,prop_values,case_sensitive): return True
        return False


    def prop_values2str(self, prop_name:str, sep=','):
        try:
            prop_values = self[prop_name]
            return sep.join(map(str,prop_values))
        except KeyError:
            return ''

    def _merge(self,with_other:"PSObject"):
            for k, vlist in with_other.items():
                self.update_with_list(k,vlist)

    def dump(self,to_file:str):
        with open(to_file+'.pickle', "wb") as outfile:
            # "wb" argument opens the file in binary mode
            pickle.dump(self, outfile)


    @classmethod
    def load(cls,from_file:str):
        dump_name = from_file+'.pickle'
        try:
            f = open(dump_name, "rb") 
            o = cls(pickle.load(f))
            f.close()
            return o
        except FileNotFoundError:
            raise FileNotFoundError


    def merge_with_rel(self,rel:"PSRelation"):
        self._merge(rel)
        self.update_with_list('References',rel._get_refs())


    def transform(self, remap:dict, remove_unmapped_props=False):
        transformed_copy = PSObject(self)
        for old_prop, new_prop in remap.items():
            try:
                transformed_copy[new_prop] = transformed_copy.pop(old_prop)
            except KeyError:
                try:
                    references = transformed_copy['References']
                    ref_values = set()
                    for ref in references:
                        ref_props = ref.get_props(old_prop)
                        ref_values.update(ref_props)
                    if ref_values:
                        transformed_copy[new_prop] = list(ref_values)
                except KeyError:
                    continue
        
        if remove_unmapped_props:
            transformed_copy = {k:v for k,v in transformed_copy.items() if k in remap.values()}

        return transformed_copy


class PSRelation(PSObject):
    pass
    PropSetToProps = dict() # {PropSetID:{PropID:[values]}}
    Nodes = dict() # {"Regulators':[(entityID, 0, effect)], "Targets':[(entityID, 1, effect)]}
    # where 0,1 direction of linktype
    References = dict() # {refIdentifier (PMID, DOI, EMBASE, PUI, LUI, Title):Reference}

    def __init__(self, dic=dict()):
        super().__init__(dic)
        self.PropSetToProps = dict()  # {PropSetID:{PropID:[values]}}
        self.Nodes = dict()  # {"Regulators':[(entityID, Dir, effect)], "Targets':[(entityID, Dir, effect)]}
        self.References = dict()  # {refIdentifier (PMID, DOI, EMBASE, PUI, LUI, Title):Reference}


    def add_references(self, references:list):#list of Reference objects
        self.load_references()
        for r in references:
            was_merged = False
            for id in PS_REFIID_TYPES:
                try:
                    self.References[r.Identifiers[id]]._merge(r)
                    was_merged = True
                    break
                except KeyError:
                    continue

            if not was_merged:
                for i in r.Identifiers.values():
                    self.References[i] = r


    def effects(self):
        try:
            return self[EFFECT]
        except KeyError:
            return []


    def effect(self):
        try:
            return self[EFFECT][0]
        except KeyError:
            raise KeyError

    def mechanisms(self):
        try:
            return self['Mechanism']
        except KeyError:
            return []


    def is_from_rnef(self):
        try:
            self.id()
            return False
        except KeyError:
            return True


    def _get_refs(self, only_with_id_type='',sort_by=PUBYEAR,reverse=True,ref_limit=0):
        self.load_references()
        all_refs = list(set(self.References.values()))
        if only_with_id_type:
            filter_refs = [ref for ref in all_refs if only_with_id_type in ref.Identifiers.keys()]
            if sort_by:
                filter_refs.sort(key=lambda r: r._sort_key(sort_by), reverse=reverse)
            return filter_refs[:ref_limit] if ref_limit else filter_refs
        else:
            if sort_by:
                all_refs.sort(key=lambda r: r._sort_key(sort_by), reverse=reverse)
            return all_refs[:ref_limit] if ref_limit else all_refs


    def merge_rel(self, other:'PSRelation'):
        self.merge_obj(other)
        self.add_references(other._get_refs())
     #   self.PropSetToProps.update(other.PropSetToProps)
     #   self.Nodes.update(other.Nodes)
        

    def copy(self):
        copy = PSRelation(self)
        copy.PropSetToProps= dict(self.PropSetToProps)
        copy.Nodes= dict(self.Nodes)
        copy.References = dict(self.References)
        return copy
        

    @ classmethod
    def make_rel(cls,regulator:PSObject,target:PSObject,props:dict,refs=[],is_directional=True):
        # props = {prop_name:[prop_values]}
        new_rel = cls(props)
        try:
            effect_val = props[EFFECT]
        except KeyError:
            effect_val = str()

        new_rel.Nodes[REGULATORS] = [(regulator['Id'][0], '0', effect_val)]
        if is_directional:
            new_rel.Nodes[TARGETS] = [(target['Id'][0], '1', effect_val)]
        else:
            new_rel.Nodes[REGULATORS].append((target['Id'][0], '0', effect_val))

        new_rel.add_references(refs)
        if REFCOUNT not in new_rel.keys():
            new_rel[REFCOUNT] = [len(refs)]
        return new_rel


    def __hash__(self):
        try:
            return self['Id'][0]
        except KeyError:
            # for ResnetGraph imported from RNEF:
            my_hash = hashlib.md5(str(self.urn()).encode())
            return int(my_hash.hexdigest(),32)


    def is_directional(self):
        return len(self.Nodes) == 2


    def _prop2str(self, prop_id, cell_sep:str =';'):
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


    def props2dict(self, prop_ids:list,cell_sep:str =';'):
        return {k:self._prop2str(k,cell_sep) for k in self.keys() if k in prop_ids}
    
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
            for ref_id_type in PS_REFIID_TYPES:
                try:
                    ref_id = propSet[ref_id_type][0]
                    my_reference_tuples.append((ref_id_type, ref_id))
                except KeyError:
                    continue

            if my_reference_tuples: # propSet is valid reference - resolving duplicates 
                # case when reference have valid identifiers
                article_identifiers = {x[1] for x in my_reference_tuples}
                existing_ref = [r for i,r in self.References.items() if i in article_identifiers]
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
                        if txtref not in ['Admin imported','Customer imported']:
                            ref = Reference.from_textref(txtref)
                            for ref_id_type, ref_id in ref.Identifiers.items():
                                self.References[ref_id] = ref
                        else: continue #'Admin imported' references have no identifiers and ignored 
                    except KeyError: continue

            textref = None 
            sentence_props = dict()
            for propId, propValues in propSet.items():  # adding all other valid properties to Ref
              #  if propValues[0] == 'info:pmid/21180064#abs:4':
              #      print('')
                if propId in BIBLIO_PROPS:
                    ref.update_with_list(propId, propValues)
                elif propId in SENTENCE_PROPS:
                    sentence_props[propId] = propValues
                elif propId in CLINTRIAL_PROPS:
                    ref.update_with_list(propId, propValues)
                elif propId == 'TextRef': textref = propValues[0]

            try:
                ref[JOURNAL] = ref.pop(MEDLINETA)
            except KeyError: pass

            if not isinstance(textref,str): 
                textref = ref._make_textref()
            elif textref[4:10] == 'hash::': # references from older Resnet versions can start with 'urn:hash::'
                textref = ref._make_textref()
            elif textref in ['Admin imported','Customer imported','']: 
                print('Reference has no TextRef property and will be ignored!!!')
                continue # ignore and move to the next PropSet
            
            if sentence_props: 
                ref.snippets[textref] = sentence_props
            else: 
                ref.snippets[textref] = {'Sentence':[]} #load empty dict for data consistency


    def filter_references(self, keep_prop2values:dict):
        '''
        !!!! does not change self[REFCOUNT] to keep the original number of references !!!!\n
        can be used to filter references by journal name

        Input
        -----
        prop_names2values = {prop_name:[values]}
        '''
        # self.References = {identifier:ref} contains reference duplications. Reference list needs compression for speed
        all_refs = set(self.References.values())
        ref2keep = {ref for ref in all_refs if ref.has_properties(keep_prop2values)}
        self.References = {k:r for k,r in self.References.items() if r in ref2keep}
        return


    def remove_reference(self, remove_prop2values:dict):
        # self.References = {identifier:ref} contains reference duplications. Reference list needs compression for speed
        all_refs = set(self.References.values()) 
        ref2keep = {ref for ref in all_refs if not ref.has_properties(remove_prop2values)}
        self.References = {k:r for k,r in self.References.items() if r in ref2keep}


    def get_reference_count(self, count_abstracts=False):
        self.load_references()
        if self.References:
            ref_count = len(self._get_refs())
            self[REFCOUNT] = [ref_count]
            if count_abstracts:
                ref_from_abstract = set([x for x in self.References.values() if x.is_from_abstract()])
                return len(ref_from_abstract)
            return ref_count
        else:
            try:
                refcount = self[REFCOUNT]
                # case when REFCOUNT was loaded from RNEF dump that did not contain references
                # e.g. for loading network from __pscache__
                refcount2merge = list(map(int,refcount))
                max_refcount = max(refcount2merge)
                self[REFCOUNT] = [max_refcount]
                return max_refcount
            except KeyError:
                self[REFCOUNT] = [0]
                return 0


    def rel_prop_str(self, sep=':'):
        # returns string of relation properties with no references
        to_return = str()
        for prop_id, prop_values in self.items():
            to_return += prop_id + sep + ','.join(prop_values)+';'


    def to_table_dict(self, columnPropNames:list, cell_sep:str=';', RefNumPrintLimit=0, add_entities=False):
        # assumes all properties in columnPropNames were fetched from Database otherwise will crash
        # initializing table
        col_count = len(columnPropNames) +2 if add_entities else len(columnPropNames)
        RelationNumberOfReferences = int(self[REFCOUNT][0])

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
                if k == REGULATORS:
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
                    propValue = self._prop2str(propId)
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
        RelationNumberOfReferences = self.get_reference_count()

        table = ['']*col_count
        if add_entities:
            regulatorIDs = str()
            targetIDs = str()
            for k, v in self.Nodes.items():
                if k == REGULATORS:
                    regulatorIDs = ','.join(map(str, [x[0] for x in v]))
                else:
                    targetIDs = ','.join(map(str, [x[0] for x in v]))

                table[col_count-2] = regulatorIDs
                table[col_count-1] = targetIDs

        for col in range(len(columnPropNames)):
            propId = columnPropNames[col]
            if propId in self.keys(): # filling columns with relation properties
                propValues = self._prop2str(propId)
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
            return self.to_table_str(columnPropNames,col_sep,cell_sep,RefNumPrintLimit,add_entities)

    def get_regulators_targets(self):
        """
        Returns
        -------
        regulator,target pairs for directional self
        all possible pairwise combinations for non-directional self
        """
        if len(self.Nodes) > 1:
            # for directional relations
            return [(r[0],t[0]) for r in self.Nodes[REGULATORS] for t in self.Nodes[TARGETS]]
        else:
            try:
                # for non-directions relations
                objIdList = [x[0] for x in self.Nodes[REGULATORS]]
                return itertools.combinations(objIdList, 2)
            except KeyError:
                # for wiered cases
                objIdList = [x[0] for x in self.Nodes[TARGETS]]
                return itertools.combinations(objIdList, 2)

    def regulator_ids(self):
        nodeIds = [x[0] for x in self.Nodes[REGULATORS]]
        return list(nodeIds)

    def target_ids(self):
        try:
            nodeIds = [x[0] for x in self.Nodes[TARGETS]]
            return list(nodeIds)
        except KeyError: return []

    def entities_ids(self):
        return self.regulator_ids()+self.target_ids()

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

    def effect_sign(self): 
        try:
            eff = self['Effect'][0]
            if eff == 'positive':
                return 1
            elif eff == 'negative':
                return -1
            else: 
                return 0
        except KeyError: 
            return 0

    
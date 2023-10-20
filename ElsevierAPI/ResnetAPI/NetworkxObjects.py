import re
import json
import itertools
import pickle
import hashlib
import datetime
import math

from ..ETM_API.references import Reference, len
from ..ETM_API.references import JOURNAL,PS_REFIID_TYPES,NOT_ALLOWED_IN_SENTENCE,BIBLIO_PROPS,SENTENCE_PROPS,CLINTRIAL_PROPS
from ..ETM_API.references import MEDLINETA,EFFECT,PUBYEAR,SENTENCE,TITLE


PROTEIN_TYPES = ['Protein','FunctionalClass','Complex']
REGULATORS = 'Regulators'
TARGETS = 'Targets'
REFCOUNT = 'RelationNumberOfReferences'
CHILDS = 'Childs'
#enums for objectypes to avoid misspeling
GENETICVARIANT = 'GeneticVariant'
FUNC_ASSOC = 'FunctionalAssociation'
DBID = 'Id'
ENTITY_IDENTIFIERS = ["CAS ID","EC Number","Ensembl ID","LocusLink ID","GO ID","GenBank ID",
                      "HMDB ID","MedScan ID","Microarray ID","PharmaPendium ID","PubChem CID","miRBase ID"]

CHEMICAL_PROPS = ["CAS ID","HMDB ID","IUPAC Name","InChIKey","LMSD ID","Alias",
                        "Molecular Formula", "Molecular Weight","NCIm ID","PubChem CID","PubChem SID",
                        "Reaxys ID","Rotatable Bond Count","XLogP","XLogP-AA","Description","PharmaPendium ID"]

DIRECT_RELTYPES = {'DirectRegulation','Binding','ChemicalReaction','ProtModification','PromoterBinding'}
        

STATE = "state"
ACTIVATED = 1
REPRESSED = -1
UNKNOWN_STATE = 0
MINREF4DIRECTREL = 5
MINAFFINITY4DIRECT = 6.0
DIRECT = 1
INDIRECT = 0


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
    
    def urn(self):
        try:
            return str(self['URN'][0])
        except KeyError:
            raise KeyError
        
    @staticmethod
    def urn2uid(urn:str):
        #return hash(urn)
        my_hash = hashlib.md5(str(urn).encode())
        return int(my_hash.hexdigest(),32)
        
    def __hash__(self):
        #__hash__ needs __eq__ to work properly
        return self.urn2uid(self.urn())

    def __eq__(self, other):
        #__hash__ needs __eq__ to work properly
        return self['URN'][0] == other['URN'][0]

    def uid(self):
        return self.__hash__()

    def dbid(self):
        try:
            return self[DBID][0]
        except KeyError:
            return 0

    def set_property(self, PropId, PropValue: str):
        self[PropId] = [PropValue]

    
    def childs(self):
        try:
            return list(self[CHILDS])
        except KeyError:
            return list()

    def child_dbids(self):
        try:
            children = self[CHILDS]
            return [c.dbid() for c in children]
        except KeyError:
            return []
        
    def child_uids(self):
        try:
            children = self[CHILDS]
            return [c.uid() for c in children]
        except KeyError:
            return []


    def append_property(self, PropId, PropValue):
        try:
            self[PropId].append(PropValue)
        except KeyError:
            self[PropId] = [PropValue]


    def name(self):
        try:
            return str(self['Name'][0])
        except KeyError:
            raise KeyError


    def objtype(self):
        try:
            return str(self['ObjTypeName'][0])
        except KeyError:
            raise KeyError


    def descr(self):
        try:
            return str(self['Description'][0])
        except KeyError:
            return ''
        

    def notes(self):
        try:
            return str(self['Notes'][0])
        except KeyError:
            return ''
        
    
    def set_state(self, state:int):
        assert (state in [ACTIVATED,REPRESSED,UNKNOWN_STATE])
        try:
            self[STATE] += state
        except KeyError:
            self[STATE] = state


    def state(self):
        try:
            return self[STATE]
        except KeyError:
            return UNKNOWN_STATE


    def is_from_rnef(self):
        try:
            database_id = self[DBID[0]]
            return False
        except KeyError:
            return True
    

    def get_prop(self,prop_name:str,value_index=-1,if_missing_return=''):
        '''
        Return
        ------
        if value_index >= 0 returns self[prop_name][value_index]\n
        else returns self[prop_name] as list
        '''
        try:
            values = self[prop_name]
            if value_index >= 0:
                return str(values[value_index])
            else:
                return values
        except KeyError:
            return if_missing_return


    def get_props(self,prop_names:list):
        '''
        Return
        ------
        returns set of all values of prop_names
        '''
        my_props = set()
        for prop_name in prop_names:
            try:
                values = self[prop_name]
                my_props.update(values)
            except KeyError:
                continue

        return my_props


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


    def update_with_list(self, prop_id:str, new_values):
        try:
            values = self[prop_id]
            if isinstance(values,list): 
            # list preserves order 
                [values.append(x) for x in new_values if x not in values]
            elif isinstance(values,set):
                [values.update(x) for x in new_values if x not in values]
            else:
                pass
            self[prop_id] = values
        except KeyError:
            self[prop_id] = new_values


    def merge_obj(self,other:'PSObject',replace_unique=False):
        '''
        properties from "other" take precedent including URN
        '''
        [self.update_with_list(prop_name,values) for prop_name,values in other.items()]
        if replace_unique:
            self['URN'] = other['URN']
            self['Name'] = other['Name']
            self['ObjTypeName'] = other['ObjTypeName']
        return self
        

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


    def has_propeties(self,prop_names:set):
        return not prop_names.isdisjoint(set(self.keys()))


    def is_annotated(self,with_prop:str,having_values:list=[],case_sensitive=False):
        '''
        Input
        -----
        if "having_values" is empty will return True if self has any value in "with_prop"
        '''
        try:
            search_in = self[with_prop]
            if having_values:
                if case_sensitive or with_prop in {DBID,'id','ID'}:
                    search_set = set(having_values)   
                else:
                    search_in  = set(map(lambda x: x.lower(),search_in))
                    search_set = set(map(lambda x: x.lower(),having_values))      
                return not search_set.isdisjoint(search_in)
            else:
                return True
        except KeyError: return False


    def has_value_in(self,prop2values:dict,case_sensitive=False):
        """
        Input
        -----
        prop2values = {propName:[values]}
        """
        for prop_name, prop_values in prop2values.items():
            if self.is_annotated(prop_name,prop_values,case_sensitive): 
                return True
        return False


    def prop_values2str(self, prop_name:str, sep=','):
        try:
            prop_values = self[prop_name]
            return sep.join(list(map(str,prop_values)))
        except KeyError:
            return ''

 #   def _merge(self,with_other:"PSObject"):
 #           for k, vlist in with_other.items():
 #               self.update_with_list(k,vlist)

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


    def transform(self, remap:dict, remove_props_not_in_remap=False):
        '''
        renames
        -------
        properties according to "remap"
        '''
        transformed_copy = PSObject(self)
        for old_prop, new_prop in remap.items():
            try:
                transformed_copy[new_prop] = transformed_copy.pop(old_prop)
            except KeyError:
                try:
                    references = transformed_copy['RefDict']
                except KeyError:
                    try:
                        references = transformed_copy['references']
                    except KeyError:
                        continue

                ref_values = set()
                for ref in references:
                    ref_props = ref.get_props(old_prop)
                    ref_values.update(ref_props)
                if ref_values:
                    transformed_copy[new_prop] = list(ref_values)
                
        if remove_props_not_in_remap:
            transformed_copy = {k:v for k,v in transformed_copy.items() if k in remap.values()}

        return transformed_copy


    def remove_props(self,prop_names:list):
        my_copy = self.copy()
        [my_copy.pop(p,'') for p in prop_names]
        return my_copy


class PSRelation(PSObject):
    pass
    PropSetToProps = dict() # {PropSetID:{PropID:[values]}}
    Nodes = dict() # {"Regulators':[(entityID, 0, effect)], "Targets':[(entityID, 1, effect)]}
    # where 0,1 direction of linktype
    RefDict = dict() # {refIdentifier (PMID, DOI, EMBASE, PUI, LUI, Title):Reference}
    references = list() # list of unique references sorted by PUBYEAR


    def __init__(self, dic=dict()):
        super().__init__(dic)
        self.PropSetToProps = dict()  # {PropSetID:{PropID:[values]}}
        self.Nodes = dict()  # {"Regulators':[(entityID, Dir, effect)], "Targets':[(entityID, Dir, effect)]}
        self.RefDict = dict()  # {refIdentifier (PMID, DOI, EMBASE, PUI, LUI, Title):Reference}


    def __hash__(self):
        #__hash__ needs __eq__ to work properly
        urn = self.urn()
        my_hash = hashlib.md5(str(urn).encode())
        return int(my_hash.hexdigest(),32)


    def __eq__(self, other):
        #__hash__ needs __eq__ to work properly
        return self['URN'][0] == other['URN'][0]
    

    def __effects(self):
        try:
            return self[EFFECT]
        except KeyError:
            return ['unknown']


    def effect(self):
        my_effects = self.__effects()
        return my_effects[0]


    def mechanisms(self):
        try:
            return self['Mechanism']
        except KeyError:
            return []
        
    def mechanism(self):
        try:
            return self['Mechanism'][0]
        except KeyError:
            return ''
        

    def has_properties(self,prop_names:set):
        has_props = super().has_propeties(prop_names)
        if not has_props:
            my_refs = self.refs()
            for prop in prop_names:
                for ref in my_refs:
                    if ref.has_property(prop):
                        return True
                    
        return has_props
        

    def make_urn(self,reg_urns:list,tar_urns:list):
        reg_urns.sort()
        tar_urns.sort()
        rel_urn = 'urn:agi-'+self.objtype()+':'
        if tar_urns:
            rel_urn += 'out:'+'out:'.join(tar_urns)
            rel_urn += ':in-out:'+'in-out:'.join(reg_urns)
            try:
                rel_urn += ':'+ self[EFFECT][0]
            except KeyError: pass
            try:
                rel_urn += ':'+self['Mechanism'][0]
            except KeyError: pass
        else:
            rel_urn += 'in-out:'+'in-out:'.join([u for u in reg_urns])

        self.set_property('URN', rel_urn)
        return str(rel_urn)
    

    def __refDict2refs(self):
        my_refs = list(set(self.RefDict.values()))
        my_refs.sort(key=lambda r: r._sort_key(PUBYEAR), reverse=True)
        self.references = my_refs
        return  self.references
    

    def _add_refs(self, references:list):
        '''
        Input
        -----
        references = [Reference]

        Return
        ------
        new list of references [Reference].\nMerges refs from input to existing refs in self.RefDict if possible
        '''
        self.load_refdict()
        for r in references:
            was_merged = False
            for id in PS_REFIID_TYPES:
                try:
                    self.RefDict[r.Identifiers[id]]._merge(r)
                    was_merged = True
                    break
                except KeyError:
                    continue

            if not was_merged:
                for i in r.Identifiers.values():
                    self.RefDict[i] = r

        return self.__refDict2refs()


    def replace_refs(self, newrefs:list):
        newrel = PSRelation(self)
        newrel.Nodes = dict(self.Nodes)
        newrel._add_refs(newrefs)
        return newrel
        

    def __load_refs(self):
        '''
        Input
        -----
        value of "with_id_type" must be from references.PS_REFIID_TYPES
        "sort_by" value must be one of keys in Reference. Defaults to PUBYEAR

        Return
        -------
        [Reference] sorted by "sort_by" in the order="reverse"
        if max_refcount > 0 returns only first max_refcount references
        '''
        self.load_refdict()
        return self.__refDict2refs()
    

    def refs(self, ref_limit=0):
        '''
        Return
        -------
        self.references sorted by PUBYEAR
        '''
        if not self.references:
            self.references = list(self.__load_refs())

        return self.references[:ref_limit] if ref_limit else self.references
       

    def pubage(self):
        today = datetime.date.today()
        this_year = today.year
        my_refs = self.refs()
        last_pubyear = my_refs[0].pubyear()
        return  this_year - last_pubyear


    def qw(self):
        '''
        Return
        -------
        publication quality Qscore
        '''
        #beta = 0.23
        return math.exp(-0.23*(self.pubage()))


    def cr(self):
        '''
        Return
        -------
        Relation citation ratio
        '''
        my_refs = self.refs()
        return float(len(my_refs))/float(self.pubage()+1)
    

    def cw(self):
        '''
        Return
        ------
        citation weight
        '''
        b = math.sqrt(3)/3
        cr = self.cr()
        cw = b*cr/math.sqrt(b*b*cr*cr +1)
        return cw
    

    def merge_rel(self, other:'PSRelation'):
        self.merge_obj(other)
        self._add_refs(other.__load_refs())
        return self
     #  self.PropSetToProps.update(other.PropSetToProps)

    
    def rel2psobj(self):
        '''
        Return
        ------
        PSObject with properties added from self\n
        rel.references are added to self['references'] attribute
        '''
        new_psobj = PSObject(self)
        new_psobj['references'] = self.refs()
        return new_psobj

        
    def make_copy(self):
        my_copy = PSRelation(self)
        my_copy.PropSetToProps= dict(self.PropSetToProps)
        my_copy.Nodes= dict(self.Nodes)
        my_copy.RefDict = dict(self.RefDict)
        my_copy.references = list(self.references)
        return my_copy
    

    def remove_props(self,prop_names:list):
        my_copy = PSRelation(super().remove_props(prop_names))
        my_copy.Nodes= dict(self.Nodes)
        my_copy.PropSetToProps= dict(self.PropSetToProps)
        for _,props in my_copy.PropSetToProps.items():
            [props.pop(p,'') for p in prop_names]

        new_RefDict = dict(self.RefDict)
        for refid,ref in new_RefDict.items():
            new_ref = ref.remove_props(prop_names)
            if isinstance(new_ref,Reference):
                my_copy.RefDict[refid] = new_ref
        
        my_copy.refs()
        return my_copy
    

    @classmethod
    def make_rel(cls,regulator:PSObject,target:PSObject,props:dict,refs=[],is_directional=True):
        # props = {prop_name:[prop_values]}
        new_rel = cls(props)
        try:
            effect_val = props[EFFECT]
        except KeyError:
            effect_val = str()

        new_rel.Nodes[REGULATORS] = [(regulator.uid(), '0', effect_val)]
        if is_directional:
            new_rel.Nodes[TARGETS] = [(target.uid(), '1', effect_val)]
        else:
            new_rel.Nodes[REGULATORS].append((target.uid(), '0', effect_val))

        new_rel._add_refs(refs)
        if REFCOUNT not in new_rel.keys():
            new_rel[REFCOUNT] = [len(refs)]
        return new_rel
            

    def is_directional(self):
        return len(self.Nodes) == 2
    

    def regulator_uids(self):
        nodeIds = [x[0] for x in self.Nodes[REGULATORS]]
        return list(nodeIds)


    def target_uids(self):
        try:
            nodeIds = [x[0] for x in self.Nodes[TARGETS]]
            return list(nodeIds)
        except KeyError: return []


    def entities_uids(self):
        return self.regulator_uids()+self.target_uids()


    def _prop2str(self, prop_id, cell_sep:str =';'):
        try:
            return cell_sep.join(list(map(str, self[prop_id])))
        except KeyError:
            prop_set_values = []
            for prop in self.PropSetToProps.values():
                try:
                    prop_set_values.append(cell_sep.join(list(map(str, prop[prop_id]))))
                except KeyError:
                    continue
            
            to_return = cell_sep.join(prop_set_values)
            if prop_id in [SENTENCE,TITLE]:
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


    def load_refdict(self,refresh=False):
        if self.RefDict:
            if not refresh:
                return
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
                existing_ref = [r for i,r in self.RefDict.items() if i in article_identifiers]
                if not existing_ref:  # case when reference is new
                    id_type, id_value = my_reference_tuples[0][0], my_reference_tuples[0][1]
                    ref = Reference(id_type, id_value)
                    self.RefDict[id_value] = ref
                    for Id in range(1, len(my_reference_tuples)):
                        id_type, id_value = my_reference_tuples[Id][0], my_reference_tuples[Id][1]
                        ref.Identifiers[id_type] = id_value
                        self.RefDict[id_value] = ref
                else: 
                    ref = existing_ref[0]
                    if len(existing_ref) > 1:  # identifiers from one propset point to different references
                    # will merge all references from propset with the first one
                        for i in range(1, len(existing_ref)):
                            ref._merge(existing_ref[i])
                            for id_value in dict(existing_ref[i].Identifiers).values():
                                self.RefDict[id_value] = ref
                        del existing_ref[1:]

            else: #propSet is not valid reference - trying to create one using Title or TexRef
                try:
                    # trying id reference by title as a last resort since it does not have valid identifiers
                    propset_title = propSet['Title'][0]
                    try:
                        ref = self.RefDict[propset_title]
                    except KeyError:
                        ref = Reference('Title', propset_title)
                        self.RefDict[propset_title] = ref
                except KeyError:
                    try:
                        txtref = propSet['TextRef'][0]
                        if txtref not in ['Admin imported','Customer imported']:
                            ref = Reference.from_textref(txtref)
                            for ref_id_type, ref_id in ref.Identifiers.items():
                                self.RefDict[ref_id] = ref
                        else: continue #'Admin imported' references have no identifiers and ignored 
                    except KeyError: continue

            textref = None 
            sentence_props = dict()
            for propId, propValues in propSet.items():  # adding all other valid properties to Ref
                if propId in BIBLIO_PROPS:
                    ref.update_with_list(propId, propValues)
                elif propId in SENTENCE_PROPS:
                    sentence_props[propId] = propValues
                elif propId == 'msrc':
                    sentence_props[SENTENCE] = propValues
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
        # self.RefDict = {identifier:ref} contains reference duplications. Reference list needs compression for speed
        all_refs = set(self.RefDict.values())
        ref2keep = {ref for ref in all_refs if ref.has_properties(keep_prop2values)}
        self.RefDict = {k:r for k,r in self.RefDict.items() if r in ref2keep}
        return


    def remove_reference(self, remove_prop2values:dict):
        # self.RefDict = {identifier:ref} contains reference duplications. Reference list needs compression for speed
        all_refs = set(self.RefDict.values()) 
        ref2keep = {ref for ref in all_refs if not ref.has_properties(remove_prop2values)}
        self.RefDict = {k:r for k,r in self.RefDict.items() if r in ref2keep}


    def get_reference_count(self, count_abstracts=False):
        self.load_refdict()
        if self.RefDict:
            ref_count = len(self.refs())
            self[REFCOUNT] = [ref_count]
            if count_abstracts:
                ref_from_abstract = set([x for x in self.RefDict.values() if x.is_from_abstract()])
                return len(ref_from_abstract)
            return ref_count
        else:
            try:
                refcount = self[REFCOUNT]
                # case when REFCOUNT was loaded from RNEF dump that did not contain references
                # e.g. for loading network from __pscache__
                if len(refcount) > 1:
                    refcount2merge = list(map(int,refcount))
                    max_refcount = max(refcount2merge)
                    self[REFCOUNT] = [max_refcount] 
                return int(self[REFCOUNT][0])
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
                    regulatorIDs = ','.join([str(x[0]) for x in v])
                else:
                    targetIDs = ','.join([str(x[0]) for x in v])

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
                    regulatorIDs = ','.join([str(x[0]) for x in v])
                else:
                    targetIDs = ','.join([str(x[0]) for x in v])

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
        [(regulator_uid,target_uid)] pairs for directional self
        all possible pairwise combinations for non-directional self
        """
        if len(self.Nodes) > 1:
            # for directional relations
            return [(r[0],t[0]) for r in self.Nodes[REGULATORS] for t in self.Nodes[TARGETS]]
        else:
            try:
                # for non-directions relations
                objIdList = [x[0] for x in self.Nodes[REGULATORS]]
                pairs = list(itertools.combinations(objIdList, 2))
                plus_reverse = list(pairs)
                [plus_reverse.append((p[1],p[0])) for p in pairs] 
                # non-directional relations are added in both direction into MultiDiGraph
                return plus_reverse
            except KeyError:
                # for wiered cases
                objIdList = [x[0] for x in self.Nodes[TARGETS]]
                pairs = list(itertools.combinations(objIdList, 2))
                plus_reverse = list(pairs)
                [plus_reverse.append((p[1],p[0])) for p in pairs]
                return plus_reverse


    def to_json(self):
        str1 = '{"Relation Properties": ' + json.dumps(self) + '}'
        strP = '{"Relation RefDict": ' + json.dumps(self.PropSetToProps) + '}'
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
            for ref in self.RefDict.values():
                ref.set_weight(max_weight)
        except KeyError:
                for ref in self.RefDict:
                    ref.set_weight(0.0)


    def effect_sign(self): 
        try:
            eff = self['Effect'][0]
            if eff == 'positive':
                return ACTIVATED
            elif eff == 'negative':
                return REPRESSED
            else: 
                return UNKNOWN_STATE
        except KeyError: 
            return UNKNOWN_STATE
        

    def _affinity(self):
        try:
            return float(self['Affinity'][0])
        except KeyError:
            try:
               return float(self['pX'][0]) 
            except KeyError:
                return -1.0

    '''
    @staticmethod
    def is_direct(rel:'PSRelation'):
        my_affinity = rel._affinity()
        if my_affinity >= MINAFFINITY4DIRECT:
            return True
        elif rel.objtype() in DIRECT_RELTYPES and rel.get_reference_count() >= MINREF4DIRECTREL:
            return True
        else:
            return False


    @staticmethod
    def is_indirect(rel:'PSRelation'):
        indirect_reltypes={'Regulation','MolTransport','Expression','MolSynthesis'}
        if rel.objtype() in indirect_reltypes:
            return True
        elif rel.objtype() in DIRECT_RELTYPES and rel.get_reference_count() < MINREF4DIRECTREL:
            return True
        elif rel._affinity() < MINAFFINITY4DIRECT and rel._affinity() > 0.0:
            return True
        else:
            return False
        '''

    def isdirect(self):
        my_affinity = self._affinity()
        if my_affinity >= MINAFFINITY4DIRECT:
            return DIRECT
        elif my_affinity > 0.0:
            return INDIRECT
        
        objtype = self.objtype()
        if objtype in DIRECT_RELTYPES:
            return DIRECT if self.get_reference_count() >= MINREF4DIRECTREL else INDIRECT
        elif objtype in {'Regulation','MolTransport','Expression','MolSynthesis'}:
            return INDIRECT
        
        return -100


    def _refprop2rel(self,ref_prop:str,relprop:str,min_max=0):
        '''
        Input
        -----
        if min_max < 0 assigns single value to relprop that is the minimum of all ref_prop values
        if min_max > 0 assigns single value to relprop that is the maximum of all ref_prop values
        min_max works only if ref_prop values are numerical
        '''
        propvals = set()
        for ref in self.refs():
            try:
                prop_val = ref[ref_prop]
                propvals.update(prop_val)
            except KeyError:
                try:
                    prop_val = ref.Identifiers[ref_prop]
                    propvals.add(prop_val)
                except KeyError:
                    for prop2values in ref.snippets.values():
                        try:
                            prop_val = prop2values[ref_prop]
                            propvals.update(prop_val)
                        except KeyError:
                            continue
        if propvals:
            if min_max:
                propvals = list(map(float,propvals))
                propvals.sort()
           #     if len(propvals) > 1:
           #         print('')
                propvals = [propvals[0]] if min_max < 0 else [propvals[-1]]

            self.update_with_list(relprop, list(propvals))

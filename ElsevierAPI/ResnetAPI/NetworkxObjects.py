import re,json,itertools,pickle,hashlib,math,copy
from datetime import datetime
from collections import defaultdict

from ..utils import normalize
from ..ETM_API.references import Reference, len, reflist2dict,pubmed_hyperlink,make_hyperlink,pmc_hyperlink
from ..ETM_API.references import JOURNAL,PS_REFIID_TYPES,NOT_ALLOWED_IN_SENTENCE,BIBLIO_PROPS,SENTENCE_PROPS,CLINTRIAL_PROPS
from ..ETM_API.references import MEDLINETA,PUBYEAR,SENTENCE,TITLE,AUTHORS,_AUTHORS_,PS_REFERENCE_PROPS

OBJECT_TYPE = 'ObjTypeName'
PROTEIN_TYPES = ['Protein','FunctionalClass','Complex']
REGULATORS = 'Regulators'
TARGETS = 'Targets'
REFCOUNT = 'RelationNumberOfReferences'
CHILDS = 'Childs'
EFFECT = 'Effect'
MECHANISM = 'Mechanism'
RELATION_PROPS = {EFFECT,MECHANISM,'ChangeType','BiomarkerType','QuantitativeType'}
ALL_PSREL_PROPS = list(RELATION_PROPS)+PS_REFERENCE_PROPS
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
NONDIRECTIONAL_RELTYPES = {'Binding','CellExpression',FUNC_ASSOC,'Metabolization','Paralog'}

STATE = "state"
ACTIVATED = 1
REPRESSED = -1
UNKNOWN_STATE = 0
MINREF4DIRECTREL = 5
MINAFFINITY4DIRECT = 6.0
DIRECT = 1
INDIRECT = 0


class PSObject(defaultdict):  # {PropId:[values], PropName:[values]}
  pass
  def __init__(self, dic:dict=dict()):
      super().__init__(list)
      if isinstance(dic,dict): # closseness passes <class:type> for some reason
          self.update(dic)

  @classmethod
  def from_zeep(cls, ZeepObjectRef):
      return_dict = dict()
      zeep_iter = iter(ZeepObjectRef)
      while True:
          try:
              item = next(zeep_iter)
          except StopIteration:
              break  # Iterator exhausted: stop the loop
          else:
              return_dict[item] = [ZeepObjectRef[item]]
      return cls(return_dict)
  
  
  def get_prop(self,prop_name:str,value_index=0,if_missing_return:str|int|float='')->str|int|float:
      '''
      output:
        self[prop_name][value_index]
      '''
      if prop_name in self:
        return self[prop_name][value_index] if value_index < len(self[prop_name]) else if_missing_return
      else:
         return if_missing_return


  def urn(self):
      '''
      returns empty string if URN property does not exist
      '''
      return self.get_prop('URN')
      
  
  def active_urn(self):
      return self.urn()+'a'
  

  def repressed_urn(self):
      return self.urn()+'i'
  

  def make_active(self):
      activated_self = PSObject(self.copy())
      activated_self['URN'] = [activated_self.active_urn()]
      return activated_self
  

  def make_repressed(self):
      repressed_self = PSObject(self.copy())
      repressed_self['URN'] = [repressed_self.repressed_urn()]
      return repressed_self


  @staticmethod
  def urn2uid(urn:str):
      #return hash(urn)
      my_hash = hashlib.md5(str(urn).encode())
      return int(my_hash.hexdigest(),32)


  def __hash__(self):
      #__hash__ needs __eq__ to work properly
      return self.urn2uid(self.urn())


  def __eq__(self, other:"PSRelation"):
      #__hash__ needs __eq__ to work properly
      return self.urn() == other.urn()


  def uid(self)->int:
      return self.__hash__()
  

  def propvalues(self,prop_name:str):
      '''
      Return
      ------
      returns set of all values of prop_names
      '''
      return list(self[prop_name])  if prop_name in self.keys() else []
  

  def get_props(self,prop_names:list):
      '''
      Return
      ------
      returns set of all values of prop_names
      '''
      my_props = set()
      for prop_name in prop_names:
          if prop_name in self.keys():
              my_props.update(self[prop_name])

      return my_props
  

  def dbid(self)->int:
      '''
      output:
          0 if missing
      '''
      return self.get_prop(DBID,if_missing_return=0)


  def set_property(self, PropId, PropValue: str):
      self[PropId] = [PropValue]


  def rename_prop(self, old_prop_name:str, new_prop_name:str):
      '''
      renames property in self
      '''
      if old_prop_name in self.keys():
        self[new_prop_name] = self.pop(old_prop_name)
        return True
      else:
        return False
  

  def childs(self)->list['PSObject']:
      if CHILDS in self.keys():
          return self[CHILDS]
      return []


  def child_dbids(self):
      children = self.childs()
      return [c.dbid() for c in children if c.dbid()]


  def child_uids(self):
      children = self.childs()
      return [c.uid() for c in children]


  def update_with_value(self, prop_id:str, new_value):
      if new_value not in self[prop_id]:
          self[prop_id].append(new_value)


  def update_with_list(self, prop_id:str, new_values):
      [self[prop_id].append(x) for x in new_values if x not in self[prop_id]]


  def name(self):
      return self.get_prop('Name')


  def objtype(self):
      return self.get_prop(OBJECT_TYPE)


  def organism(self):
      return str(self.get_prop('Organism'))


  def descr(self):
      return str(self.get_prop('Description'))
      

  def notes(self):
      return str(self.get_prop('Notes'))

  
  def set_state(self, state:int):
      assert (state in [ACTIVATED,REPRESSED,UNKNOWN_STATE])
      if STATE in self:
          self[STATE][0] += state
      else:
          self[STATE] = [state]


  def state(self):
      return str(self.get_prop(STATE,if_missing_return=UNKNOWN_STATE))
  

  def is_from_rnef(self):
      return not self.dbid()
  

  def copy(self):
      return PSObject(self)


  def merge_obj(self,other:'PSObject',replace_identity=False):
      '''
      properties from "other" take precedent,
      if replace_identity URN, Name, ObjTypeName are also replaced
      '''
      my_copy = PSObject(self)
      [my_copy.update_with_list(prop_name,values) for prop_name,values in other.items()]
      if replace_identity:
        my_copy['URN'] = other['URN']
        my_copy['Name'] = other['Name']
        my_copy['ObjTypeName'] = other['ObjTypeName']
      del other
      return my_copy
      

  def _prop2str(self, prop_id:str,cell_sep:str =';'):
      prop_values = self.propvalues(prop_id)
      return cell_sep.join(map(str,prop_values))


  def props2dict(self,prop_ids:list):
      return {k:self._prop2str(k) for k in self.keys() if k in prop_ids}


  def data2str(self, columnPropNames: list, col_sep='\t', cell_sep=';', endOfline='\n'):
      table_row = str()
      for propName in columnPropNames:
          values = self.propvalues(propName)
          prop_val = cell_sep.join(values)
          table_row = table_row + prop_val + col_sep
      return table_row[0:len(table_row) - 1] + endOfline


  def has_properties(self,prop_names:set):
      return not prop_names.isdisjoint(set(self.keys()))


  def is_annotated(self,with_prop:str,having_values:list=[],case_sensitive=False):
      '''
      Input
      -----
      if "having_values" is empty will return True if self has any value in "with_prop"
      '''
      search_in = self.get_props([with_prop])
      if search_in:
          if having_values:
              if case_sensitive or with_prop in {DBID,REFCOUNT}:
                  search_set = set(having_values)   
              else:
                  search_in  = set(map(lambda x: x.lower(),search_in))
                  search_set = set(map(lambda x: x.lower(),having_values))      
              return not search_set.isdisjoint(search_in)
          else:
              return True
      else: 
          return False


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
      prop_values = self.propvalues(prop_name)
      return sep.join(list(map(str,prop_values)))


  def dump(self,to_file:str):
      with open(to_file+'.pickle', "wb") as outfile:
          # "wb" argument opens the file in binary mode
          pickle.dump(dict(self), outfile)


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
          transformed_copy[new_prop] = transformed_copy.pop(old_prop)
          references = transformed_copy['references']

          ref_values = set()
          [ref_values.update(r.get_props(old_prop)) for r in references]
          if ref_values:
              transformed_copy[new_prop] = list(ref_values)
              
      if remove_props_not_in_remap:
          transformed_copy = {k:v for k,v in transformed_copy.items() if k in remap.values()}

      return transformed_copy


  def remove_props(self,prop_names:list):
      my_copy = self.copy()
      [my_copy.pop(p,'') for p in prop_names]
      return my_copy


  def remap(self,mapdic:dict[str,dict[str,dict[str,"PSObject"]]],map_by:list,normalize_propvalues=True):
    '''
    input:
      use normalized mapdic as input for best mapping
      map_by - defines list of properties for mapping and their order 
    output:
      if URN is remapped self['URN'] = [newURN,oldURN] and self is flagged by 'was_mapped' property
    '''
    for prop in map_by:
      try:
        provals = self.get_props([prop])
        for val in provals:
          v = val.lower()
          try:
            db_obj = mapdic[self.objtype()][prop][v]
            new_urn = db_obj.urn()
            old_urn = self.urn()
            if new_urn != old_urn:
              self['URN'] = [new_urn,old_urn]
            self['was_mapped'] = [True]
            return True
          except KeyError:
            if normalize_propvalues:
              v = normalize(val)
              try:
                db_obj = mapdic[self.objtype()][prop][v]
                new_urn = db_obj.urn()
                old_urn = self.urn()
                if new_urn != old_urn:
                  self['URN'] = [new_urn,old_urn]
                self['was_mapped'] = [True]
                return True
              except KeyError:
                continue
            else:
              continue
      except KeyError:
        continue
    return False


  def is_mapped(self):
    return 'was_mapped' in self


class PSObjectEncoder(json.JSONEncoder):
  pass
  def default(self, obj):
    if isinstance(obj, (PSObject,defaultdict)):
      obj['hook'] = 'PSObject'
      obj_dict = {k:v for k,v in obj if not str(k).endswith('Id')}
      return dict(obj_dict)  # Convert PSObject to a regular dictionary
    return super().default(obj)  # Let the base class handle other types


class PSObjectDecoder(json.JSONDecoder):
  def __init__(self, *args, **kwargs):
    json.JSONDecoder.__init__(self, object_hook=self.object_hook, *args, **kwargs)

  def object_hook(self, obj):
    hook = obj.pop('hook','')
    if hook == 'PSObject':
      return PSObject(obj)
    return obj


class PSRelation(PSObject):
  '''
  PropSetToProps = {PropSetID:{PropID:[values]}}
  Nodes = {"Regulators':[PSObject], "Targets':[PSObject]}, 
  references - list of unique references sorted by PUBYEAR in descending order
  '''
  pass

  def __init__(self, dic=dict()):
      '''
      self.Nodes - {"Regulators':[PSObject], "Targets':[PSObject]}
      '''
      super().__init__(dic)
      self.PropSetToProps = defaultdict(dict)  # {PropSetID:{PropID:[values]}}
      self.Nodes = defaultdict(list)  # {"Regulators':[(entityID, Dir, effect)], "Targets':[(entityID, Dir, effect)]}
      self.references = list() # has to be list for sorting


  def __hash__(self):
      #__hash__ needs __eq__ to work properly
      urn = self.urn()
      my_hash = hashlib.md5(str(urn).encode())
      return int(my_hash.hexdigest(),32)


  def __eq__(self, other):
      #__hash__ needs __eq__ to work properly
      return self['URN'][0] == other['URN'][0]


  def effect(self,unknown='unknown')->str:
      '''
      output:
          'positive','negative',unknown
      '''
      return str(self.get_prop(EFFECT,if_missing_return=unknown))
  

  def flip_effect(self):
    '''
    flip effect from positive to negative and vice versa\n
    output:
        True if effect was flipped;
        False if effect was not changed
    '''
    if self.effect() == 'positive':
      self[EFFECT][0] = 'negative'
      return True
    elif self.effect() == 'negative':
      self[EFFECT][0] = 'positive'
      return True
    return False


  def mechanisms(self):
      return self.propvalues(MECHANISM)
    

  def mechanism(self):
      return str(self.get_prop(MECHANISM))
  

  def name(self):
    if 'Name' in self.keys():
      return self['Name'][0]
    else:
      regulators, targets = self.Nodes[REGULATORS],self.Nodes[TARGETS]
      targ_names = ','.join([t['Name'][0] for t in targets])
      if targ_names:
          reg_names = ','.join([r['Name'][0] for r in regulators])
          arrow = str('--->')
          effect = self.effect()
          if effect == 'positive': 
              arrow = '--+>'
          elif effect == 'negative': 
              arrow = '---|'
      else:
          reg_names = regulators[0]['Name'][0]
          targ_names =  ','.join([r['Name'][0] for r in regulators[1:]])
          arrow = '----'
      name = reg_names+'----'+self.objtype()+arrow+targ_names
      if MECHANISM in self.keys():
          name += ':'+self[MECHANISM][0]

      self['Name'] = [name]
      return name
    

  def urn(self,refresh=False):
    '''
    Creates URN for rel if it does not exist
    '''
    if refresh:
      regulators = self.regulators()
      targets = self.targets()
      return self.__make_urn(regulators, targets)
    else:
      urn = super().urn()
      if urn:
          return urn
      else:
          regulators = self.regulators()
          targets = self.targets()
          return self.__make_urn(regulators, targets)
          

  def has_properties(self,prop_names:set):
      has_props = super().has_properties(prop_names)
      if not has_props:
          my_refs = self.refs()
          for prop in prop_names:
              for ref in my_refs:
                  if ref.has_property(prop):
                      return True
                  
      return has_props
  

  def has_value_in(self,prop2values:dict[str,list],case_sensitive=False):
      has_values = super().has_value_in(prop2values,case_sensitive)
      if not has_values:
          for ref in self.refs():
              if ref.has_values_in(prop2values):
                  return True
      return has_values


  def __make_urn(self,regulators:list[PSObject],targets:list[PSObject]):
      '''
      sets:
          URN property to self
      output:
          calculated URN
      '''
      reg_urns = sorted([r.urn() for r in regulators])
      tar_urns = sorted([r.urn() for r in targets])
      rel_urn = 'urn:agi-'+self.objtype()+':'
      if tar_urns:
          rel_urn += 'out:'+'out:'.join(tar_urns) # out:urn:out:
          rel_urn += ':in-out:'+'in-out:'.join(reg_urns)
          effect = self.effect()
          if effect in ['positive','negative']:
              rel_urn += ':'+ effect
          if MECHANISM in self:
              rel_urn += ':'+self[MECHANISM][0]
      else:
          rel_urn += 'in-out:'+':in-out:'.join([u for u in reg_urns])

      self.set_property('URN', rel_urn)
      return str(rel_urn)
  

  def rename_prop(self, old_prop_name:str, new_prop_name:str):
    '''
    Renames property in self
    '''
    if not super().rename_prop(old_prop_name, new_prop_name):
        for ref in self.refs():
          ref.rename_prop(old_prop_name, new_prop_name)

    return
   

  def _add_refs(self, references:list[Reference])->list[Reference]:
    '''
    output:
      new self.refernces\nmerges refs from input to existing refs in self.reference
    '''
    self.refs()
    if len(self.references) > len(references):
      larger_list = self.references
      smaller_list = references
    else:
      larger_list = references
      smaller_list = self.references

    ref_dict = reflist2dict(larger_list)
    for ref in smaller_list:
      was_merged = False
      for id_type in PS_REFIID_TYPES:
        try:
          identifier = ref.Identifiers[id_type]
          id_str = ref.identifiers_str(id_type,identifier)
          try:
            if ref != ref_dict[id_str]:
              ref_dict[id_str]._merge(ref)
            was_merged = True
          except KeyError:
              continue
        except KeyError:
            continue
        
      if not was_merged:
        ref_dict.update({ref.identifiers_str(k,v):ref for k,v in ref.Identifiers.items()})
        
    self.references = list(set(ref_dict.values()))
    return self.references


  def replace_refs(self, newrefs:list):
      newrel = PSRelation(self) # copy only properties and Nodes. 
      # PropSetToProps, RedDict and references are skipped
      newrel.Nodes = self.Nodes.copy()
      newrel._add_refs(newrefs)
      return newrel
      

  def _1st_ref(self):
      '''
      Return
      ------
      the earliest reference or empty Reference if publication year is unavalialble 
      '''
      for ref in reversed(self.refs()):
          pubyear = ref.pubyear()
          if pubyear > 1812:
              return ref
      empty_ref = Reference('','')
      empty_ref.Identifiers.clear()
      return empty_ref


  def pubage(self):
      today = datetime.today()
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
    '''
    does not copy other.Nodes
    '''
    my_copy = self.copy()
    my_copy.merge_obj(other)
    my_copy._add_refs(other.refs())
    del other
    return my_copy

  
  def number_of_snippets(self):
    return sum([ref.number_of_snippets() for ref in self.refs()])
      

  def textrefs(self):
    my_refs = self.refs()
    return sum([r.textrefs() for r in my_refs],[])


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


  def copy(self):
    my_copy = PSRelation(self)
    my_copy.PropSetToProps = copy.deepcopy(self.PropSetToProps)
    my_copy.Nodes = copy.deepcopy(self.Nodes)
    my_copy.references = self.references.copy()
    return my_copy
  

  def remove_props(self,prop_names:list):
    my_copy = PSRelation(super().remove_props(prop_names))
    my_copy.Nodes= dict(self.Nodes)
    my_copy.PropSetToProps= dict(self.PropSetToProps)
    for _,props in my_copy.PropSetToProps.items():
        [props.pop(p,'') for p in prop_names]

    my_copy.references.clear()
    for ref in self.references:
        new_ref = ref.remove_props(prop_names)
        if new_ref:
          my_copy.references.append(new_ref)
    
    return my_copy
  

  @classmethod
  def make_rel(cls,regulator:PSObject,target:PSObject,props:dict[str,list],refs:list[Reference],is_directional=True):
    new_rel = cls(props)
    if is_directional:
      new_rel.Nodes[REGULATORS] = [regulator]
      new_rel.Nodes[TARGETS] = [target]
    else:
      new_rel.Nodes[REGULATORS] = [regulator,target]

    new_rel._add_refs(refs)
    if REFCOUNT not in new_rel.keys():
      new_rel[REFCOUNT] = [len(refs)]
    
    new_rel.__make_urn([regulator],[target])
    return new_rel
          

  def is_directional(self):
      return len(self.Nodes) == 2


  def regulator_uids(self):
      nodeIds = [x.uid() for x in self.Nodes[REGULATORS]]
      return list(nodeIds)
  
  
  def target_uids(self)->list[int]:    
      return [x.uid() for x in self.Nodes[TARGETS]] if TARGETS in self.Nodes else []


  def entities_uids(self):
      return self.regulator_uids()+self.target_uids()
  

  def get_props(self,propname:str)->list:
    if propname in self:
      return self[propname]
    else:
      prop_set_values = set()
      my_refs = self.refs()
      [prop_set_values.update(ref.get_values(propname)) for ref in my_refs]
      return list(prop_set_values)


  def _prop2str(self, propname:str, cell_sep:str =';'):
    propvals = self.get_props(propname)
    to_return = cell_sep.join(map(str,propvals))
    if propname in [SENTENCE,TITLE]:
      to_return = re.sub(NOT_ALLOWED_IN_SENTENCE, ' ', to_return)
    return to_return
  

  def props2dict(self, prop_ids:list,cell_sep:str =';'):
      return {k:self._prop2str(k,cell_sep) for k in self.keys() if k in prop_ids}
  

  def props2list(self, propID, cell_sep=';'):
      try:
          return self[propID]
      except (IndexError,KeyError):
          to_return = list()
          for prop in self.PropSetToProps.values():
              prop_set_val = str()
              for prop_id, values in prop.items():
                  if prop_id == propID:
                      prop_set_val = cell_sep.join(values)
                      break
              to_return.append(prop_set_val)
          return to_return


  def refs(self,refresh=False,ref_limit=0)->list[Reference]:
    '''
    output:
      self.references sorted by PUBYEAR
    '''
    if not self.references or refresh: # making self.references from self.PropSetToProps:
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
          article_identifiers = {t[1] for t in my_reference_tuples}
          existing_ref = [r for r in self.references if not article_identifiers.isdisjoint(r.Identifiers.keys())]
          if not existing_ref:  # case when reference is new
            first_tuple =  my_reference_tuples[0]
            ref = Reference(first_tuple[0], first_tuple[1])
            ref.Identifiers.update({t[0]:t[1] for t in my_reference_tuples[1:]})
            self.references.append(ref)
          else: # duplicate ref exists !
            ref = existing_ref[0]
            for i in range(1, len(existing_ref)):
              existref = existing_ref[i]
              ref._merge(existref)
            del existing_ref[1:] # deleting duplicate refs to save memory
        else: #propSet is not valid reference - trying to create one using Title or TexRef
            try:
              # trying id reference by title as a last resort since it does not have valid identifiers
              propset_title = propSet['Title'][0]
              refs_with_same_title = [r for r in self.references if r.title() == propset_title]
              if refs_with_same_title:
                  ref = ref[0]
              else:
                  ref = Reference('Title', propset_title)
                  self.references.append(ref)
            except KeyError: # no title :(
              try:
                  txtref = propSet['TextRef'][0]
                  if txtref not in ['Admin imported','Customer imported']:
                      ref = Reference.from_textref(txtref)
                      self.references.append(ref)
                  else: continue #'Admin imported' references have no identifiers and ignored 
              except KeyError: continue

        # adding all other valid properties to ref
        try:
          textref = propSet['TextRef'][0]
          if textref[4:10] == 'hash::': # references from older Resnet versions can start with 'urn:hash::'
            textref = ref._make_textref()
          elif textref in ['Admin imported','Customer imported','']: 
            print('Reference has no TextRef property and will be ignored!!!')
            continue # ignore and move to the next PropSet
        except KeyError:
          textref = ref._make_textref()
          
        for propId, propValues in propSet.items():  
          if propId in BIBLIO_PROPS or propId in CLINTRIAL_PROPS:
            ref.update_with_list(propId, propValues)
          elif propId in SENTENCE_PROPS:
            ref.add_sentence_props(textref,propId, propValues)
          elif propId == 'msrc':
            ref.add_sentence_props(SENTENCE, propValues)

        if MEDLINETA in ref:
          ref[JOURNAL] = ref.pop(MEDLINETA)

        if AUTHORS not in ref and _AUTHORS_ in ref:
          ref[AUTHORS] = [x.tostr() for x in ref[_AUTHORS_]]
    
    [x.toAuthors() for x in self.references]
    self.references.sort(key=lambda r: r._sort_key(by_property=PUBYEAR), reverse=True)
    return self.references[:ref_limit] if ref_limit else self.references

  
  def filter_references(self, keep_prop2values:dict,in_place=True):
      '''
      !!!! does not change self[REFCOUNT] to keep the original number of references !!!!\n
      can be used to filter references by journal name

      Input
      -----
      prop_names2values = {prop_name:[values]}
      '''
      all_refs = self.refs()
      ref2keep = [ref for ref in all_refs if ref.has_values_in(keep_prop2values)]
      if in_place:
        self.references = ref2keep
      return self.references


  def remove_references(self, with_prop2values:dict):
      '''
      Input
      -----
      prop2values = {prop_name:[values]}
      '''
      all_refs = self.refs()
      ref2keep = {ref for ref in all_refs if not ref.has_values_in(with_prop2values)}
      if ref2keep:
          if len(ref2keep) < len(all_refs):
              mycopy = self.copy()
              mycopy.references = self.references.copy()
              return mycopy
          else:
              return self
      else:
          return PSRelation()
      

  def is_from_abstract(self):
      for ref in self.refs():
          if ref.is_from_abstract():
              return True
      return False


  def count_refs(self, count_abstracts=False)->int:
      if self.references:
        self[REFCOUNT] = [len(self.references)]
        if count_abstracts:
          ref_from_abstract = set([x for x in self.references if x.is_from_abstract()])
          return len(ref_from_abstract)
      else:
          refcount = self[REFCOUNT]
          # case when REFCOUNT was loaded from RNEF dump as string without references
          # e.g. for loading network from __pscache__
          if len(refcount) > 1: # case if refcount has 2 or more values after relation merge
              refcount2merge = list(map(int,refcount))
              self[REFCOUNT] = [max(refcount2merge)]
          else:
              self[REFCOUNT] = [int(refcount[0])] if refcount else [0]
      
      return int(self[REFCOUNT][0])


  def rel_prop_str(self, sep=':'):
    # returns string of relation properties with no references
    to_return = str()
    for prop_id, prop_values in self.items():
      to_return += prop_id + sep + ','.join(prop_values)+';'


  def is_annotated(self,with_prop:str,having_values:list=[],case_sensitive=False):
    if not super().is_annotated(with_prop,having_values,case_sensitive):
      for ref in self.refs():
         if ref.has_values_in({with_prop:having_values}):
          return True
      return False
    else:
       return True


  def to_table_dict(self, columnPropNames:list, cell_sep:str=';', RefNumPrintLimit=0, add_entities=False)->dict[int,str]:
      '''
      Return
      ------
      {rownum:[column_values]}
      '''
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
                  regulatorIDs = ','.join([str(x.dbid()) for x in v])
              else:
                  targetIDs = ','.join([str(x.dbid()) for x in v])

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
      RelationNumberOfReferences = self.count_refs()

      table = ['']*col_count
      if add_entities:
          regulatorIDs = str()
          targetIDs = str()
          for k, v in self.Nodes.items():
              if k == REGULATORS:
                  regulatorIDs = ','.join([str(x.dbid()) for x in v])
              else:
                  targetIDs = ','.join([str(x.dbid()) for x in v])

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


  def regulators(self)->list[PSObject]:
    return self.Nodes[REGULATORS] if REGULATORS in self.Nodes.keys() else []
      
  def targets(self)->list[PSObject]:
    return self.Nodes[TARGETS] if TARGETS in self.Nodes.keys() else []
      

  def get_regulators_targets(self)->list[tuple[int,int]]:
      """
      Returns
      -------
      [(regulator_uid,target_uid)] pairs for directional self
      all possible pairwise combinations for non-directional self
      """
      if self.is_directional():
          return [(r.uid(),t.uid()) for r in self.Nodes[REGULATORS] for t in self.Nodes[TARGETS]]
      else:
          # for non-directions relations
          assert(len(self.Nodes) == 1)
          if REGULATORS in self.Nodes:
              uIdList = [x.uid() for x in self.Nodes[REGULATORS]]
          else:
              uIdList = [x.uid() for x in self.Nodes[TARGETS]]
  
          pairs = list(itertools.combinations(uIdList, 2))
          if pairs:
              plus_reverse = list(pairs)
              [plus_reverse.append((p[1],p[0])) for p in pairs] 
              # non-directional relations are added in both direction into MultiDiGraph
              return plus_reverse
          else:
              if len(uIdList) > 1: # case of self-loop
                  return [(uIdList[0],uIdList[1])]
              else:
                  return list()


  def to_json(self):
      str1 = '{"Relation Properties": ' + json.dumps(self) + '}'
      strP = '{"Relation references": ' +json.dumps(self.references)
      strR = json.dumps(self.Nodes)
      return str1 + '\n' + strP + '\n' + strR + '\n'
  

  def _set_weight2ref (self, weight:float,weight_name='weight'):
      for ref in self.refs():
          ref.set_weight(weight,weight_name)


  def set_weight2ref (self, weight_by_prop_name:str, proval2weight:dict[str,float],weight_name='weight'):
      try:
          values2weight = set(self[weight_by_prop_name])
          max_weight = 0.0
          for v in values2weight:
              try:
                  weight = proval2weight[v]
                  if weight > max_weight: max_weight = weight
              except KeyError: continue
          for ref in self.references:
              ref.set_weight(max_weight,weight_name)
      except KeyError:
              for ref in self.references:
                  ref.set_weight(0.0,weight_name)


  def effect_sign(self): 
      try:
          eff = self['Effect'][0]
          if eff == 'positive':
              return ACTIVATED
          elif eff == 'negative':
              return REPRESSED
          else: 
              return UNKNOWN_STATE
      except (IndexError,KeyError): 
          return UNKNOWN_STATE
      

  def pX(self):
      pXs = [float(pX) for ref in self.refs() for pX in ref.get_values('pX')]
      return max(pXs) if pXs else -1.0
  

  def set_affinity(self):
      pXs = [float(pX) for ref in self.refs() for pX in ref.get_values('pX')]
      if pXs:
          self['Affinity'].append(max(pXs))
      return
  

  def _affinity(self):
      if 'Affinity' in self:
          return float(self['Affinity'][0])
      elif 'pX' in self:
              return float(self['pX'][0]) 
      else:
          return -1.0


  def isprimarytarget(self):
      '''
      ouput:
        True, if pX > 6.0 or reltype in DIRECT_RELTYPES if refcount >=5
      '''
      my_affinity = self._affinity()
      if my_affinity >= MINAFFINITY4DIRECT:
          return DIRECT
      elif my_affinity > 0.0:
          return INDIRECT
      
      objtype = self.objtype()
      if objtype in DIRECT_RELTYPES:
          return DIRECT if self.count_refs() >= MINREF4DIRECTREL else INDIRECT
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


  def todict(self,include_refs=False)->tuple[dict[str,list],list[tuple[str,dict[str,list]]]]:
    dic = dict(self)
    if include_refs:
      relname = self.name()
      my_refs = self.refs()
      return dic,[x.todict(relname=relname) for x in my_refs] 
    else:
      return dic,[]
    

  def hyperlinked_refcount(self,reflimit:int=10):
    last_refs = self.refs(ref_limit=reflimit)
    rel_pmids = list()
    rel_dois = list()
    rel_pmcids = list()
    for ref in last_refs:
      id_type, identifier = ref.get_doc_id()
      if id_type == 'PMID':
        if identifier not in rel_pmids:
          rel_pmids.append(identifier)
      elif id_type == 'DOI':
        if identifier not in rel_dois:
          rel_dois.append(identifier)
      elif id_type == 'PMC':
        if identifier not in rel_pmcids:
          rel_pmcids.append(identifier)
  
    refcount = self.count_refs()
    if rel_pmids:
      return pubmed_hyperlink(list(rel_pmids),str(refcount)),refcount
    elif rel_dois:
      return make_hyperlink(rel_dois[0],'http://dx.doi.org/',str(refcount)),refcount
    elif rel_pmcids:
      return pmc_hyperlink(rel_pmcids,str(refcount)),refcount
    else:
      return '',refcount

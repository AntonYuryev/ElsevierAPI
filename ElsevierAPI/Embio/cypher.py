from ..ResnetAPI.NetworkxObjects import OBJECT_TYPE,REFCOUNT,PSObject,CONNECTIVITY,EFFECT
from ..ResnetAPI.ResnetGraph import PHYSICAL_INTERACTIONS, PROTEIN_TYPES

class Cypher:
  @staticmethod
  def __urns(nodes:list[PSObject]):
    return [n.urn() for n in nodes]


  @staticmethod
  def match_psobjs(objs:list[PSObject],letter='a'):
    list_name = f'{letter}_urnList'
    objtypes = {o.objtype() for o in objs}
    objtypes_str = '|'.join(objtypes)
    cypher = f"WITH ${list_name} AS urns\nUNWIND urns AS urn\nMATCH ({letter}:{objtypes_str}{{URN:urn}})\n"
    parameter = {f'{list_name}':[obj.urn() for obj in objs]}
    return cypher, parameter
  

  def match_nodes_byprop(propValues:list[str|int],propName:str,objtype:str,letter='a'):
    list_name = f'{letter}_propList'
    cypher = f"WITH ${list_name} AS props\nUNWIND props AS prop\nMATCH ({letter}:{objtype}{{{propName}:prop}})\n"
    parameter = {f'{list_name}':propValues}
    return cypher, parameter
  


  @staticmethod
  def match_drugs(letter='d'):
    '''
    Needs MATCH and WHERE
    '''
    return f""" (
            ({letter})-[:is_a*]->(:SmallMol {{Name:'plant medicinal product'}}) 
            OR ({letter})-[:is_a*]->(:SemanticConcept {{Name:'drugs'}})
            )
            AND NOT ({letter})-[:is_a]->(:SmallMol {{Name:'PAINS compounds'}})
            """


  @staticmethod
  def select_drugs(only_from:list[PSObject]=[]):
    if only_from:
      cypher ="""WITH $drug_urnList AS urns
                UNWIND urns AS urn
                MATCH (d:SmallMol {URN:urn}) WHERE """
      cypher += Cypher.match_drugs()
      cypher += "RETURN d"
      return cypher, {'drug_urnList':Cypher.__urns(only_from)}
    else:
      cypher = 'MATCH (d:SmallMol) WHERE '
      cypher += Cypher.match_drugs()
      cypher += "RETURN d"
      return cypher, dict()
  

  @staticmethod
  def select_drug_targets(_4targets:list[PSObject],only_from:list[PSObject]=[],
                relProps:dict[str,list[str|int|float]]=[]) -> tuple[str,dict[str,list]]:
    '''
    input:
      only_from - list of drugs to limit the graph
    '''
    targettypes = {o.objtype() for o in _4targets}
    targettypes_str = ','.join([f"'{t}'" for t in targettypes])

    cypher, params = Cypher.select_drugs(only_from)
    cypher = cypher[:-8] # remove "RETURN d"
    cypher += f'MATCH (t) WHERE any(label IN labels(t) WHERE label IN [{targettypes_str}])\n'
    cypher += f'AND t.URN IN $targeturnList\n'
    cypher +=  'MATCH (d)-[r]->(t)'
    if relProps:
      cypher = Cypher.add_relProps(cypher,relProps)
    cypher += '\nRETURN d,r,t'

    params.update({'targeturnList':[o.urn() for o in _4targets]})
    return cypher, params


  @staticmethod
  def create_group(group_name:str, memberUrns:list[str]):
    """
      Cypher query requires {memberUrns:[urns]} as parameter for session.run(cypher,parameter)
    """
    return f"""
      MERGE (p:Group {{Name: '{group_name}'}})
      WITH p, $memberUrns AS memberUrnsList
      UNWIND memberUrnsList AS currentUrn
      OPTIONAL MATCH (c:SmallMol {{URN: currentUrn}})
      WITH p, c
      WHERE c IS NOT NULL
      MERGE (c)-[:part_of]->(p)
      WITH DISTINCT p
      OPTIONAL MATCH (linked_member:SmallMol)-[:part_of]->(p)
      RETURN p, collect(linked_member.URN) AS LinkedUrns
  """, {'memberUrns':memberUrns}


  @staticmethod
  def create_ontology_group(group_name:str, concept_type:str,memberUrns:list[str]):
    """
      Cypher query requires {memberUrns:[urns]} as parameter for session.run(cypher,parameter)
    """
    return f"""
      MERGE (p:{concept_type} {{Name: '{group_name}'}})
      WITH p, $memberUrns AS memberUrnsList
      UNWIND memberUrnsList AS currentUrn
      OPTIONAL MATCH (c:SmallMol {{URN: currentUrn}})
      WITH p, c
      WHERE c IS NOT NULL
      MERGE (c)-[:is_a]->(p)
      WITH DISTINCT p
      OPTIONAL MATCH (linked_member:SmallMol)-[:is_a]->(p)
      RETURN p, collect(linked_member.URN) AS LinkedUrns
  """, {'memberUrns':memberUrns}


  @staticmethod
  def get_nodes(objtype:str,propName:str,propVals:list[str],with_connectivity=False):
    """
    input:
      objtype can be empty

    output:
      cypher query, parameters from propVals for session.run(cypher,parameter)
      if with_connectivity = True, cypher query returns list of (node,connectivity) tuples
      otherwise returns list of nodes
    """
    cypher = f"""
      WITH $values AS valueList
      UNWIND valueList AS value"""
    if objtype:
      cypher += f' MATCH (c:{objtype} {{{propName}:value}})'
    else:
      cypher += f' MATCH (c {{{propName}:value}})'

    cypher +=  ' RETURN c AS node'
    if with_connectivity:
      cypher += f', COUNT{{(c)-[]-()}} AS {CONNECTIVITY}'
    return cypher, {'values':propVals}

  
  @staticmethod
  def match_childs(parent:PSObject,letter='c',max_childs=10):
    return f"""MATCH ({letter})
      WHERE ({letter})-[:is_a*..{max_childs}]->(:{parent.objtype()}{{{'URN'}:'{parent.urn()}'}})\n"""
  

  @staticmethod
  def get_childs(parent:PSObject,max_childs:int=None):
    '''
      if max_childs is specified, only parents with less than max_childs childs are returned
      set max_childs to None to get all childs
    '''
    if max_childs:
      cypher = f"""
              MATCH (p:{parent.objtype()} {{URN:$urn}})
              OPTIONAL MATCH (child)-[:is_a*]->(p)
              WITH p, collect(DISTINCT child) AS all_children  
              WITH size(all_children) AS count, all_children
              RETURN count, 
                CASE WHEN count > 0 AND count <= {max_childs} THEN all_children ELSE null 
                END AS childs
            """
    else:
      cypher = f"""
                MATCH (p:{parent.objtype()} {{URN:$urn}})
                OPTIONAL MATCH (child)-[:is_a*]->(p)
                WITH p, collect(DISTINCT child) AS all_children 
                RETURN size(all_children) AS count, all_children AS childs
            """
    return cypher, {'urn':parent.urn()} # urn is returned as parametert to overcome URNs with ' sign
  

  @staticmethod
  def node_connectivity(nodes:list[PSObject]):
    objtypes = {p.objtype() for p in nodes}
    objtypes_str = '|'.join(objtypes)
    cypher = f"""WITH $urnList AS urns\nUNWIND urns AS urn\nMATCH (n:{objtypes_str}{{URN:urn}})
            """
    cypher += f'\nRETURN urn as urn, COUNT{{(n)-[]-()}} AS {CONNECTIVITY}'
    return cypher, {'urnList':[n.urn() for n in nodes]}


  @staticmethod ## TO DO 
  def _get_childs_(parents:list[PSObject],max_childs=10,letter='c'):
    objtypes = {p.objtype() for p in parents}
    objtypes_str = '|'.join(objtypes)
    cypher = f'WITH $urnList AS urns\nUNWIND urns AS urn\n'
    cypher += f'MATCH ({letter}) WHERE ({letter})-[:is_a* ..{max_childs}]->(:{objtypes_str}{{URN:urn}})\n'
    cypher += f'RETURN {letter}'
  

  @staticmethod
  def __list2str(values:list[str|int|float])->str:
    if isinstance(values[0], str):
      return ', '.join([f"'{v}'" for v in values])
    else:
      return ', '.join([str(v) for v in values])


  @staticmethod
  def add_relProps(cypher: str, relProps: dict[str, list[str | int | float]]) -> str:
    """
      Appends WHERE clauses to a Cypher query based on properties.
    """
    if not relProps:
      return cypher

    conditions = []
    for prop, values in relProps.items():
      if not values:
        continue
      if prop == OBJECT_TYPE:
          conditions.append(f"type(r) IN [{Cypher.__list2str(values)}]")
      elif prop == REFCOUNT:
          val = values[0] if isinstance(values, list) else values
          conditions.append(f"r.{prop} {val}")
      elif prop == EFFECT:
        unknown_effect_requested = 'unknown' in values
        standard_effects = [v for v in values if v not in ('unknown')]

        effect_logic_parts = []
        if standard_effects:
          effect_logic_parts.append(f"r.Effect IN [{Cypher.__list2str(standard_effects)}]")
        if unknown_effect_requested:
          effect_logic_parts.append("coalesce(r.Effect, '_') IN ['unknown', '_']")
        
        if effect_logic_parts:
          conditions.append(f"({' OR '.join(effect_logic_parts)})")
      else:
        conditions.append(f"r.{prop} IN [{Cypher.__list2str(values)}]")

    if conditions:
        return f"{cypher}\nWHERE {' AND '.join(conditions)}"
    
    return cypher

  '''
  @staticmethod
  def add_relPropsOLD(cypher:str, relProps:dict[str,list[str|int|float]]):
    """
    input:
      cypher MUST NOT have WHERE and RETURN clauses yet
      by_relProps = {propName:[propValue1,propValue2,...]},
      use OBJECT_TYPE string to specify filtering by relation type
    """
    #how2connect can send relProps with empty values
    my_relprops = {k:v for k,v in relProps.items() if v}
    if my_relprops:
      cypher += '\nWHERE '
      for prop, values in my_relprops.items():
        if prop == OBJECT_TYPE:
          cypher += f'type(r) IN [{Cypher.__list2str(values)}]\nAND '
        elif prop == REFCOUNT:
          cypher += f'r.{prop} {values}\nAND '
        elif prop == EFFECT:
          if 'uknown' in values:
            vals = values.remove('uknown')
            if vals:
              cypher += f'r.Effect IN [{Cypher.__list2str(vals)}]\nAND ' # adding 'positive','negative'
              cypher += "coalesce(r.Effect, '') IN ['unknown', '']\nAND"
        else:
          cypher += f'r.{prop} IN [{Cypher.__list2str(values)}]\nAND '
      return cypher[:-4] # remove last AND
    else:
      return cypher
    '''


  @staticmethod
  def expand_upstream(seeds_with_values:list[str],in_prop='Name', 
                    _2neighbor_types:list[str]=[],by_relProps:dict[str,list[str|int|float]]={}):
    '''
      by_relProps = {reltype:[propValue1,propValue2,...]},
      use OBJECT_TYPE string to specify filtering by relation type
    '''
    cypher = f'WITH $propList AS props\nUNWIND props AS prop\nMATCH (c{{{in_prop}:prop}})'
    if _2neighbor_types:
      neigbors = '|'.join(_2neighbor_types)
      cypher += f'MATCH (n:{neigbors})-[r]->(c)\n'
    else:
      cypher += 'MATCH (n)-[r]->(c)'

    parameters = {'propList':seeds_with_values}
    cypher = Cypher.add_relProps(cypher, by_relProps)
    cypher += '\nRETURN n,r,c'
    return cypher,parameters
  

  @staticmethod
  def expand_downstream(seeds_with_values:list[str],in_prop='Name', 
                    _2neighbor_types:list[str]=[],by_relProps:dict[str,list[str|int|float]]={}):
    '''
      by_relProps = {reltype:[propValue1,propValue2,...]},
      use OBJECT_TYPE string to specify filtering by relation type
    '''
    cypher = f'WITH $propList AS props\nUNWIND props AS prop\nMATCH (c{{{in_prop}:prop}})'
    if _2neighbor_types:
      neigbors = '|'.join(_2neighbor_types)
      cypher += f'MATCH (c)-[r]->(n:{neigbors})\n'
    else:
      cypher += 'MATCH (c)-[r]->(n)'

    parameters = {'propList':seeds_with_values}
    cypher = Cypher.add_relProps(cypher, by_relProps)
    cypher += '\nRETURN c,r,n'
    return cypher,parameters
  

  @staticmethod
  def expand(seeds_with_values:list[str],in_prop='Name', 
                    _2neighbor_types:list[str]=[],by_relProps:dict[str,list[str|int|float]]={},dir=''):
    '''
      by_relProps = {reltype:[propValue1,propValue2,...]},
      use OBJECT_TYPE string to specify filtering by relation type
      dir: '', 'upstream', 'downstream'
    '''
    if dir == 'upstream':
      return Cypher.expand_upstream(seeds_with_values,in_prop,_2neighbor_types,by_relProps)
    elif dir == 'downstream':
      return Cypher.expand_downstream(seeds_with_values,in_prop,_2neighbor_types,by_relProps)
    else:
      cypher = f'WITH $propList AS props\nUNWIND props AS prop\nMATCH (c{{{in_prop}:prop}})\n'
      if _2neighbor_types:
        neigbors = '|'.join(_2neighbor_types)
        cypher += f'MATCH (n:{neigbors})-[r]-(c)\n'
      else:
        cypher += 'MATCH (n)-[r]-(c)'
        
      parameters = {'propList':seeds_with_values}
      cypher = Cypher.add_relProps(cypher, by_relProps)
      cypher += '\nRETURN c,r,n'
      return cypher,parameters


  @staticmethod
  def connect(regulator_objtypes:list[str], regulator_props:list[str],regulator_propName:str,
    target_objtypes:list[str], target_props:list[str]=[],target_propName='',
    by_relProps:dict[str,list[str|int|float]]={}, dir=False):
  
    cypher = 'MATCH (a) WHERE any(label IN labels(a) WHERE label IN $reglabelList)\n'
    parameters = {'reglabelList':regulator_objtypes,'tarlabelList':target_objtypes}
    if regulator_props:
      cypher += f'AND a.{regulator_propName} IN $regpropList\n'
      parameters['regpropList'] = regulator_props

    cypher += 'MATCH (b) WHERE any(label IN labels(b) WHERE label IN $tarlabelList)\n'
    if target_props:
      cypher += f'AND b.{target_propName} IN $tarpropList\n'
      parameters['tarpropList'] = target_props

    if dir:
      cypher += 'MATCH (a)-[r]->(b)'
    else:
      cypher += 'MATCH (a)-[r]-(b)'
    cypher = Cypher.add_relProps(cypher, by_relProps)

    cypher += '\nRETURN startNode(r) AS Regulator, r AS Relation, endNode(r) AS Target'
    return cypher, parameters
  

  @staticmethod
  def ppi(interactors:list[PSObject], minref:int=2):
    '''
      interactors: list of PSObject
      output:
        cypher query, parameters for session.run(cypher,parameter)
    '''
    physical_interactions = '|'.join([v for v in PHYSICAL_INTERACTIONS])
    proteins = '|'.join([v for v in PROTEIN_TYPES])
    cypher = f'''
      WITH $urnList AS urns
      UNWIND urns AS urn
      MATCH (a:{proteins} {{URN:urn}})-[r:{physical_interactions}]-(b:{proteins})
      WHERE r.{REFCOUNT} >= {minref}
      RETURN a, r, b
    '''
    parameter = {'urnList':[obj.urn() for obj in interactors]}
    return cypher, parameter

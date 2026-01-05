from ..ResnetAPI.ResnetGraph import ResnetGraph, PSObject, PSRelation
from ..utils import execution_time, load_api_config, multithread, ThreadPoolExecutor,as_completed
from ..ResnetAPI.NetworkxObjects import PS_REFIID_TYPES,OBJECT_TYPE,NONDIRECTIONAL_RELTYPES,CHILDS,CONNECTIVITY,DBID
import logging,neo4j, time
from neo4j import GraphDatabase
from neo4j import ManagedTransaction as tx
from neo4j.exceptions import ServiceUnavailable
from .cypher import Cypher
from .postgres import PostgreSQL

RELATIONID = 'RelationID'
NODECOLUMN2ATTR = {'id':DBID,'urn':'URN'}

nondirectional_reltype = list(map(str.upper,NONDIRECTIONAL_RELTYPES))

REL_PROP_Neo4j = ['Name', 'Effect', 'Mechanism', 'Source', 'TextRef']
ENT_PROP_Neo4j = ['URN', 'Name', 'Description']
REL_PROPs = list(PS_REFIID_TYPES) + REL_PROP_Neo4j


class nx2neo4j(GraphDatabase):
  def __init__(self, APIconfig:dict={}):
    '''
    required kwargs: uriNeo4j, userNeo4j, password, database
    '''
    if not APIconfig:
      APIconfig = load_api_config()

    self.NODE_REGISTRY = set()
    self.uri = APIconfig['neo4juri']
    self.database = APIconfig['neo4jdb']
    self.user =  APIconfig['neo4juser']
    self.password = APIconfig['neo4jpswd']
    self.__driver__ = super().driver(self.uri, auth=(self.user, self.password))
    self.postgres = PostgreSQL()
  

  def session(self):
    return self.__driver__.session(database=self.database)


  def close(self):
      # Don't forget to close the driver connection when you are finished with it
      self.driver.close()


  @staticmethod
  def __record2psobj(node_record):
    psobj = PSObject({NODECOLUMN2ATTR.get(k,k):[v] for k,v in node_record._properties.items() if v not in ['_','']})
    psobj[OBJECT_TYPE] =  list(node_record.labels)
    return psobj


  def _record2psrel(self, triple:neo4j.Record)->PSRelation:
      '''
        triple: regulator-relation-target
      '''
      assert(len(triple) == 3), 'Only triples Regulator-relationship-target are considered'
      regulator = PSObject({k:[v] for k,v in triple[0]._properties.items() if v not in ['_','']})
      regulator[OBJECT_TYPE] =  list(triple[0].labels)
      regulator['URN'] =  regulator.pop('urn',regulator['URN'])
      target = PSObject({k:[v] for k,v in triple[2]._properties.items() if v not in ['_','']})
      target[OBJECT_TYPE] =  list(triple[2].labels)
      target['URN'] =  target.pop('urn',target['URN'])
      rel = triple[1]
      reldict = {k:[v] for k,v in rel._properties.items() if v not in ['_','']}
      reldict[OBJECT_TYPE] = [rel.type]
      is_directional = reldict[OBJECT_TYPE][0] not in nondirectional_reltype
      rel_obj = PSRelation.make_rel(regulator,target,reldict,[],is_directional)
      return rel_obj


  def fetch_graph(self,cypher:str,parameters=dict(),request_name='')->ResnetGraph:
    '''
    input:
      cypher must MATCH (regulator)-[relation]->(target) 
    '''
    fetch_refs = parameters.pop('with_references',True)
    with self.session() as session:
      psrels = []
      try:
        neo4j_result = list(session.run(cypher,parameters))
        psrels = [self._record2psrel(record) for record in neo4j_result]
        if fetch_refs:
          relation_ids = [n[RELATIONID][0] for n in psrels]
          self.postgres.submit_refs(relation_ids)
        if psrels:
          to_return = ResnetGraph.from_rels(psrels)
          if request_name:
            print(f'Cypher query "{request_name}" found data:')
          print(f"Loaded network with {len(to_return)} nodes and {to_return.number_of_edges()} edges")
          return to_return
        else:
          if request_name:
            print(f'Cypher query "{request_name}" did not fetch any data')
          return ResnetGraph()
      except Exception as e:
        print(f"Error during network retrival: {e}")
        raise

  
  def _connect_(self,regulator_objtypes:list[str], regulator_props:list[str],regulator_propName:str,
                target_objtypes:list[str], target_props:list[str],target_propName:str,
                by_relProps:dict[str,list[str|int|float]]={}, dir=False):
      '''
      input:
        if not dir connects regulators, targets in BOTH directions, otherwise connects regulator->target
        by_relProps = {reltype:[propValue1,propValue2,...]},
      '''
      cypher1, params1 = Cypher.connect(regulator_objtypes, regulator_props,regulator_propName,
                                        target_objtypes, target_props,target_propName,by_relProps,dir)
      return self.fetch_graph(cypher1, params1)



  def connect_objs(self,regulators:set[PSObject],targets:set[PSObject],
                   by_relProps:dict[str,list[str|int|float]]={}, dir=False)->ResnetGraph:
    regtypes = {n.objtype() for n in regulators}
    targtypes = {n.objtype() for n in targets}
    regurns = [n.urn() for n in regulators]
    targurns = [n.urn() for n in targets]
    return self._connect_(list(regtypes),regurns,'URN',
                          list(targtypes),targurns,'URN',
                          by_relProps,dir)
  
  
  def get_ppi(self,interactors:set[PSObject], minref=2,with_references=True)->ResnetGraph:
    cypher, params = Cypher.ppi(list(interactors), minref=minref)
    params.update({'with_references':with_references})
    return self.fetch_graph(cypher, params)
  

  def _neighborhood_(self,seedProps:list[str],propType='Name',_2neighbor_types:list[str]=[],
                    by_relProps:dict[str,list[str|int|float]]={},dir=''):
    '''
    input:
      by_relProps = {reltype:[propValue1,propValue2,...]},
      use OBJECT_TYPE string to specify filtering by relation type
      dir: '', 'upstream', 'downstream'
    '''
    cypher,param = Cypher.expand(seedProps,propType,_2neighbor_types,by_relProps,dir)        
    return self.fetch_graph(cypher,param)


  def create_group(self,group_name:str,link_type:str,members:list[PSObject]):
    '''
    if link_type = "part_of" , creates Group with members
    if link_type = "is_a" creates ontology concept with label equal to label of the first member with members linked by "is_a"
    '''
    memberUrns = ResnetGraph.urns(members)
    if link_type == "part_of":
      cypher,params = Cypher.create_group(group_name,memberUrns)
    else:
      label = members[0].objtype()
      cypher,params = Cypher.create_ontology_group(group_name,label,memberUrns)
      
    with self.session() as session:
      result = session.run(cypher,params)
      record = result.single()
      if record:
        p_name = record["p"]["Name"]
        linked_count = len(record["LinkedUrns"])
        print(f"✅ Success! Merged '{p_name}' and linked {linked_count} new compounds.")
        print(f"Linked URNs: {record['LinkedUrns']}")
      else:
        print("⚠️ Query ran, but no group was created.")


  def add_connectivity(self,to_nodes:list[PSObject]):
     cypher,params = Cypher.node_connectivity(to_nodes)
     with self.session() as session:
      result = list(session.run(cypher,params))
      urn2connectivity = dict()
      for record in result:
        urn = record['urn']
        connectivity = record[CONNECTIVITY]
        urn2connectivity[urn] = connectivity

      for n in to_nodes:
        n[CONNECTIVITY] = [urn2connectivity[n.urn()]]
      return to_nodes


  def load_children(self,parent:PSObject,max_childs:int=None,with_connectivity= False)->list[PSObject]:
    '''
    output:
      if max_childs is None or zero loads all children,
      otherwise loads children only for parents with number of children less than max_childs
      parents with number of children exceeding max_childs are annotated with 
      parent[CHILDS] = [PSObject()]*count
    '''
    if CHILDS in parent:
      return parent
    children = []
    with self.session() as session:
      cypher,params = Cypher.get_childs(parent, max_childs)
      record = session.run(cypher,params).single()
      if record:
        count = record['count']
        children_records = record['childs']
        if count > 0:
          if not max_childs or count <= max_childs:
            children = [self.__record2psobj(record) for record in children_records]
            if with_connectivity:
              children = self.add_connectivity(children)
          else:
            children = [PSObject()]*count
    parent[CHILDS] = children
    return parent
  

  def _load_children_(self,parents:list[PSObject],max_childs=None)->list[PSObject]:
    '''
    output:
      list of parent annotated with CHILDS attributed
      if max_childs is None or zero loads all children,
      otherwise loads children only for parents with number of children less than max_childs
      parents with number of children exceeding max_childs are annotated with 
      parent[CHILDS] = [PSObject()]*count
    '''
    def process_single(parent:PSObject):
      return self.load_children(parent,max_childs)
    
    results = []
    with ThreadPoolExecutor(max_workers=20) as executor:
      futures = executor.map(process_single,parents)
      for res in futures:
        results.append(res)
    return results
    

  def get_nodes(self,objtype:str,propName:str,propVals:list[str],
                with_childs=False,with_connectivity=False)->list[PSObject]:
    """
    input:
      objtype (label) can be empty, but the query will be slower
    """
    cypher,params = Cypher.get_nodes(objtype,propName,propVals,with_connectivity)
    with self.session() as session:
      result = list(session.run(cypher,params))
      if with_connectivity:
        nodes = []
        for record in result:
          node = self.__record2psobj(record[0])
          node.update_with_value(CONNECTIVITY,record[1])
          nodes.append(node)
      else:
        nodes = [self.__record2psobj(record[0]) for record in result]

      if with_childs:
        childs = []
        for node in nodes:
          self.load_children(node,with_connectivity=with_connectivity)
          childs += node[CHILDS]
        nodes += childs
      return set(nodes)
    

  def select_drugs(self,only_from:list[PSObject]=[]):
    cypher,params = Cypher.select_drugs(only_from)
    with self.session() as session:
      result = list(session.run(cypher,params))
      return [self.__record2psobj(record[0]) for record in result]
    


         
    
################ LOAD INTO NEO4J ####################### LOAD INTO NEO4J ########################
  @staticmethod
  def __get_node_labels(node: PSObject):
      lbls = ','.join([k + ':\"' + v[0] + '\"' for k, v in node.items() if k in ENT_PROP_Neo4j])
      return lbls


  @staticmethod
  def _find_and_return_node(tx:tx, node:PSObject):
      query = (
          'MATCH (n:{type}) WHERE n.URN = \"{urn}\"'
          'RETURN n.Name AS Name, n.URN as urn, labels(n) as ObjTypeName'
          # labels(n) returns list
      )
      query = query.format(type=node.objtype(), urn=node.urn())
      neo4j_result = tx.run(query)
      return {(record["Name"], record["urn"], record['ObjTypeName'][0]) for record in neo4j_result}
  

  def __create_node(self, tx:tx, n: PSObject):
      node_found = self._find_and_return_node(tx, n)
      if len(node_found) > 0:
          return {(node[0], node[1], node[2][0]) for node in node_found}
      else:
          node_type = n['ObjTypeName'][0]
          query = ('CREATE (n:' + node_type + '{' + self.__get_node_labels(n) + '})'
                                                                                'RETURN n.Name as Name, n.URN as urn, labels(n) as ObjTypeName'
                    )  # labels(n) returns list
          neo4j_result = tx.run(query)
          try:
              return {(record["Name"], record["urn"], record['ObjTypeName'][0]) for record in neo4j_result}
          # Capture any errors along with the query and data for traceability
          except ServiceUnavailable as exception:
              logging.error(f"{query} raised an error:\n{exception}")
              raise


  def load_nodes_1by1(self, resnet:ResnetGraph):
      nodes = resnet.nodes(data=True)
      with self.driver.session() as session:
          # Write transactions allow the driver to handle retries and transient errors
          for i, d in nodes:
              neo4j_result = session.execute_write(self.__create_node, d)
              for record in neo4j_result:
                  print(f'Neo4j got node: \"{record[0]}\" of type {record[2]} with URN={record[1]}')
  
  
  def __create_nodes(self, tx:tx, nodes:list):
      create_node_count = 0
      for n in nodes:
          node_found = self._find_and_return_node(tx, n)
          if len(node_found) > 0:
              continue
          else:
              node_type = n['ObjTypeName'][0]
              query = ('CREATE (n:' + node_type + '{' + self.__get_node_labels(n) + '})'
                          'RETURN n.Name as Name, n.URN as urn, labels(n) as ObjTypeName'
                      )  # labels(n) returns list
              tx.run(query)
              create_node_count += 1
      return create_node_count

  
  def load_node_list(self, nodes:list):
      with self.driver.session() as session:
          # Write transactions allow the driver to handle retries and transient errors
          create_node_count = session.execute_write(self.__create_nodes, nodes)
      return create_node_count


  def load_nodes(self,resnet:ResnetGraph,max_wokers=5):
      start = time.time()
      nodes = resnet._get_nodes()
      with ThreadPoolExecutor(max_workers=max_wokers, thread_name_prefix='loading nodes') as executor:
          chunk_len = int(len(nodes)/max_wokers)
          chunks = [nodes[x:x+chunk_len] for x in range(0, len(nodes), chunk_len)]
          futures = list()
          for chunk in chunks:
              futures.append(executor.submit(self.load_node_list,chunk))
          
          create_node_count = 0
          for f in futures:
              create_node_count += f.result()

      print('%d nodes were loaded in %s' % (create_node_count,execution_time(start)))
      print('%d nodes were found in the database' % (len(nodes) - int(create_node_count)))


  @staticmethod
  def __get_rel_labels(rel: PSRelation):
      lbls = ','.join([k.replace(':', ' ') + ':\"' + v[0] + '\"' for k, v in rel.items() if k in REL_PROP_Neo4j])
      lbls = lbls + ',RefCount:' + str(rel.count_refs())
      lbls = lbls + ',AbstractCount:' + str(rel.count_refs(count_abstracts=True))
      return lbls


  def __create_rel_query(self, node1: PSObject, node2: PSObject, relation: PSRelation):
      relationType = relation['ObjTypeName'][0]
      node1urn = str(node1['URN'][0])
      node2urn = str(node2['URN'][0])
      node1type = node1['ObjTypeName'][0]
      node2type = node2['ObjTypeName'][0]

      return ('MATCH'
              '(a:' + node1type + '),'
                                  '(b:' + node2type + ')'
                                                      'WHERE a.URN = \"' + node1urn + '\" AND b.URN = \"' + node2urn + '\"'
                                                                                                                        'CREATE (a)-[r:' + relationType + ' {' + self.__get_rel_labels(
          relation) + '}]->(b)'
                      'RETURN a.Name as rName, b.Name as tName, type(r) as rel_type'
              )


  def __create_relation(self, tx:tx, node1: PSObject, node2: PSObject, relation: PSRelation):
      query = self.__create_rel_query(node1, node2, relation)
      neo4j_result = tx.run(query)
      try:
          return {(record['rName'], record['rel_type'], record['tName']) for record in neo4j_result}
      # Capture any errors along with the query and data for traceability
      except ServiceUnavailable as exception:
          logging.error(f"{query} raised an error:\n{exception}")
          raise


  def load_relations_1by1(self, resnet:ResnetGraph):
      with self.driver.session() as session:
          rel_counter = 0
          for regulatorID, targetID, rel in resnet.edges.data('relation'):
              neo4j_result = session.execute_write(self.__create_relation, resnet.nodes[regulatorID],
                                                  resnet.nodes[targetID], rel)
              for record in neo4j_result:
                  print(f"Created {record[1]}: {record[0]} -> {record[2]} relation")
                  rel_counter += 1
                  if rel_counter%10000 == 0:
                      print(f'\n\nImported {rel_counter} relations\n\n')


  def __create_relations(self, tx:tx, rel_tuples:list,resnet:ResnetGraph):
      '''
      Input
      -----
      rel_tuples = [(node1_id, node2_id, PSRelation)]
      '''
      create_rel_count = 0
      for t in rel_tuples:
          query = self.__create_rel_query(resnet.nodes[t[0]], resnet.nodes[t[1]], t[2])
          tx.run(query)
          create_rel_count += 1
      return create_rel_count


  def load_relation_list(self, edge_tuples:list,resnet:ResnetGraph):
      '''
      Input
      -----
      rel_tuples = [(node1_id, node2_id, PSRelation)]
      '''
      with self.session() as session:
          created_relation_count = session.execute_write(self.__create_relations,edge_tuples,resnet)
      return int(created_relation_count)
          

  def load_relations_multithread(self,resnet:ResnetGraph, max_workers=100):
      start = time.time()
      max_edge_in_split = int(resnet.number_of_edges()/max_workers)
      multithreads = resnet.split(max_edge_in_split)

      create_relation_count = 0
      for multithread in multithreads:
          with ThreadPoolExecutor(max_workers=len(multithread), thread_name_prefix='loading relations') as executor:
              futures = list()
              for thread in multithread:
                  tread_edges = thread.edges()
                  futures.append(executor.submit(self.load_relation_list,tread_edges,resnet))
              
              for f in futures:
                  create_relation_count += f.result()

      print('%d relation were loaded in %s' % (create_relation_count,execution_time(start)))


  def load_graph2neo4j(self, resnet:ResnetGraph):
      resnet_size = resnet.number_of_edges()
      print('Importing Resnet with %d edges into local Neo4j' % resnet_size)
      import_start = time.time()
      resnet.load_references()
      if resnet_size > 50000:
          self.load_nodes(resnet)
          self.load_relations_multithread(resnet)
      else:
          self.load_nodes_1by1(resnet)
          self.load_relations_1by1(resnet)

      print("Graph with %d nodes and %d edges was imported into Neo4j in %s ---" % 
          (resnet.number_of_nodes(), resnet_size, execution_time(import_start)))

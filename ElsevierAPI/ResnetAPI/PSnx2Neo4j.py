from .ResnetGraph import ResnetGraph, PSObject, PSRelation
from ..utils import execution_time
from .NetworkxObjects import PS_REFIID_TYPES,OBJECT_TYPE,NONDIRECTIONAL_RELTYPES
import logging,neo4j, time
from neo4j import GraphDatabase
from neo4j import ManagedTransaction as tx
from neo4j.exceptions import ServiceUnavailable
from concurrent.futures import ThreadPoolExecutor
nondirectional_reltype = list(map(str.upper,NONDIRECTIONAL_RELTYPES))

REL_PROP_Neo4j = ['Name', 'Effect', 'Mechanism', 'Source', 'TextRef']
ENT_PROP_Neo4j = ['URN', 'Name', 'Description']
REL_PROPs = list(PS_REFIID_TYPES) + REL_PROP_Neo4j


class nx2neo4j(GraphDatabase):
  def __init__(self, **kwargs):
      '''
      required kwargs: uriNeo4j, userNeo4j, password, database
      '''
      self.NODE_REGISTRY = set()
      self.uri = kwargs['uriNeo4j']
      self.user =  kwargs['userNeo4j']
      self.password = kwargs['password']
      self.database = kwargs['database']
      self.__driver__ = super().driver(self.uri, auth=(self.user, self.password))
    

  def session(self):
    return self.__driver__.session(database=self.database)


  def close(self):
      # Don't forget to close the driver connection when you are finished with it
      self.driver.close()


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
  

  def record2psrel(self, triple:neo4j.Record)->PSRelation:
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


  def load_network(self,cypher:str):
    psrels = []
    def cypher2resnet(cypher:str,tx:tx):
      try:
        neo4j_result = tx.run(cypher)
        for record in neo4j_result:
          psrel = self.record2psrel(record)
          psrels.append(psrel)
        return ResnetGraph.from_rels(psrels)
      except Exception as e:
        print(f"Error during network loading transaction: {e}")
        raise

    with self.session() as session:
      loaded_graph_data = session.execute_read(lambda tx_inner: cypher2resnet(cypher, tx_inner))
      print(f"Network loaded successfully: {loaded_graph_data}")

    return loaded_graph_data


#   def find_node(self, node: PSObject):
#       with self.driver.session() as session:
#           neo4j_result = session.execute_read(self._find_and_return_node, node)
#           for record in neo4j_result:
#               print(f"Found node: {record}")


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
      with self.driver.session() as session:
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

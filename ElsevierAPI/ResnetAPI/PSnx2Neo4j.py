from ElsevierAPI.ResnetAPI.ZeepToNetworkx import PSObject, PSRelation, PSNetworx
from ElsevierAPI.ResnetAPI.NetworkxObjects import REF_ID_TYPES
import logging
from neo4j import GraphDatabase
from neo4j.exceptions import ServiceUnavailable

REL_PROP_Neo4j = ['Name', 'Effect', 'Mechanism', 'Source', 'TextRef']
ENT_PROP_Neo4j = ['URN', 'Name', 'Description']
REL_PROPs = list(REF_ID_TYPES) + REL_PROP_Neo4j


class nx2neo4j:
    def __init__(self, uriNeo4j: str, userNeo4j: str, password: str):
        self.NODE_REGISTRY = set()
        self.uri = uriNeo4j
        self.user = userNeo4j
        self.password = password
        self.driver = GraphDatabase.driver(self.uri, auth=(self.user, self.password))

    def close(self):
        # Don't forget to close the driver connection when you are finished with it
        self.driver.close()

    @staticmethod
    def __get_node_labels(node: PSObject):
        lbls = ','.join([k + ':\"' + v[0] + '\"' for k, v in node.items() if k in ENT_PROP_Neo4j])
        return lbls

    @staticmethod
    def _find_and_return_node(tx, node: PSObject):
        query = (
            'MATCH (n:{type}) WHERE n.URN = \"{urn}\"'
            'RETURN n.Name AS Name, n.URN as urn, labels(n) as ObjTypeName'
            # labels(n) returns list
        )
        query = query.format(type=node['ObjTypeName'][0], urn=node['URN'][0])
        result = tx.run(query)
        return {(record["Name"], record["urn"], record['ObjTypeName'][0]) for record in result}

    def find_node(self, node: PSObject):
        with self.driver.session() as session:
            result = session.read_transaction(self._find_and_return_node, node)
            for record in result:
                print("Found node: {record}".format(record=record))

    def __create_node(self, tx, n: PSObject):
        node_found = self._find_and_return_node(tx, n)
        if len(node_found) > 0:
            return {(node[0], node[1], node[2][0]) for node in node_found}
        else:
            node_type = n['ObjTypeName'][0]
            query = ('CREATE (n:' + node_type + '{' + self.__get_node_labels(n) + '})'
                                                                                  'RETURN n.Name as Name, n.URN as urn, labels(n) as ObjTypeName'
                     )  # labels(n) returns list
            result = tx.run(query)
            try:
                return {(record["Name"], record["urn"], record['ObjTypeName'][0]) for record in result}
            # Capture any errors along with the query and data for traceability
            except ServiceUnavailable as exception:
                logging.error("{query} raised an error: \n {exception}".format(query=query, exception=exception))
                raise

    def load_nodes(self, resnet: PSNetworx):
        nodes = resnet.Graph.nodes(data=True)
        with self.driver.session() as session:
            # Write transactions allow the driver to handle retries and transient errors
            for i, d in nodes:
                result = session.write_transaction(self.__create_node, d)
                for record in result:
                    print("Neo4j got node: \"{n}\" of type {t} with URN={u}".format(n=record[0], t=record[2],
                                                                                    u=record[1]))

    @staticmethod
    def __get_rel_labels(rel: PSRelation):
        lbls = ','.join([k.replace(':', ' ') + ':\"' + v[0] + '\"' for k, v in rel.items() if k in REL_PROP_Neo4j])
        lbls = lbls + ',RefCount:' + str(rel.get_reference_count())
        lbls = lbls + ',AbstractCount:' + str(rel.get_reference_count(count_abstracts=True))
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

    def __create_relation(self, tx, node1: PSObject, node2: PSObject, relation: PSRelation):
        query = self.__create_rel_query(node1, node2, relation)
        result = tx.run(query)
        try:
            return {(record['rName'], record['rel_type'], record['tName']) for record in result}
        # Capture any errors along with the query and data for traceability
        except ServiceUnavailable as exception:
            logging.error("{query} raised an error: \n {exception}".format(
                query=query, exception=exception))
            raise

    def load_relations(self, resnet: PSNetworx):
        with self.driver.session() as session:
            for regulatorID, targetID, rel in resnet.Graph.edges.data('relation'):
                result = session.write_transaction(self.__create_relation, resnet.Graph.nodes[regulatorID],
                                                   resnet.Graph.nodes[targetID], rel)
                for record in result:
                    print("Created {r}: {p1} - {p2} relation".format(r=record[1], p1=record[0], p2=record[2]))

    def load_graph2neo4j(self, resnet: PSNetworx):
        print('Importing Resnet data into local Neo4j')
        #app = nx2neo4j(uriNeo4j, userNeo4j, pswdNeo4j)
        self.load_nodes(resnet)
        resnet.Graph.count_references()
        self.load_relations(resnet)

from ElsevierAPI.ResnetAPI.ResnetGraph import ResnetGraph
from ElsevierAPI.ResnetAPI.PSnx2Neo4j import nx2neo4j
from ElsevierAPI import execution_time
import time

cache_path = 'ElsevierAPI/ResnetAPI/__pscache__/'
cache_file = cache_path+'protein_expression_network'+'.rnef'
nx2neo = nx2neo4j(uriNeo4j='bolt://localhost:11003', userNeo4j='neo4j', password='password')

if __name__ == "__main__":
    expression_network = ResnetGraph.fromRNEF(cache_file)

    # importing ResnetGraph into local Neo4j instance:
    import_start = time.time()
    nx2neo.load_graph2neo4j(expression_network)
    print("Graph was imported into Neo4j in %s ---" % execution_time(import_start))
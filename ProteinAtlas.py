import time
from ElsevierAPI import open_api_session
from ElsevierAPI.ResnetAPI.NetworkxObjects import REF_ID_TYPES
import ElsevierAPI.ResnetAPI.PathwayStudioGOQL as GOQL

global_start = time.time()

REL_PROPs = ['Name', 'Effect', 'Mechanism', 'ChangeType']# add here relation properties to retrieve
# if properties from NetworkxObjects.REF_ID_TYPES or NetworkxObjects.REF_PROPS are added to REL_PROPs then:
# output size may increase dramatically because it will contain one reference per row.
ENT_PROPs = ['Name', 'Description', 'Cell Localization'] 
ps_api = ps_api = open_api_session()
 
ps_api.PageSize = 10000
ps_api.add_rel_props(REL_PROPs+REF_ID_TYPES)
ps_api.add_ent_props(ENT_PROPs)

# this dump file will list all proteins in the database with connectivity >0:
ps_api.add_dump_file('Proteins from database.tsv', replace_main_dump=True)
print('Fetching all proteins from the database')
ProteinsOnyGraph = ps_api.process_oql("Select Entity WHERE objectType = Protein AND Connectivity > 0 AND Name LIKE 'A%'", flush_dump=True)


ps_api.add_dump_file("Protein neighbors dump.tsv", replace_main_dump=True)  # dump file accumulates all data in one big file
out_dir = 'csv'
counter = 0
for node_id, psObj in ProteinsOnyGraph.nodes(data=True):
    protein_name = psObj['Name'][0]
    counter += 1
    print('Finding neighbors for \"%s\", node #%d from %d total' %
          (protein_name, counter, ProteinsOnyGraph.number_of_nodes()))
    
    oql_query = GOQL.expand_entity([node_id], SearchByProperties=['id'])
    ProteinNeighborsGraph = ps_api.process_oql(oql_query)
    protein_neighbors_file = out_dir + '/' + protein_name + '_neighbors.csv'
    reference_pandas = ps_api.to_pandas(ProteinNeighborsGraph)
    reference_pandas.to_csv(protein_neighbors_file, index=False)
    ps_api.Graph.clear()  # need to release memory when performing large dumps

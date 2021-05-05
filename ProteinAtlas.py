import time
global_start = time.time()

REL_PROPs = ['Name','Effect','Mechanism','ChangeType'] #add here relation properties to retreive
#if properties from NetworkxObjects.REF_ID_TYPES or NetworkxObjects.REF_PROPS are added to REL_PROPs
# output size may increase dramaticaly because it will contain one reference per row.   
ENT_PROPs = ['Name','Description','Cell Localization'] #add here node properties to retreive

from ElsevierAPI import networx as PSnx
import ElsevierAPI.ResnetAPI.ResnetAPISession as ssn
sn = ssn.APISession('Select Entity WHERE objectType = Protein AND Connectivity > 0', PSnx)
sn.PageSize = 10000
sn.AddRelProps(REL_PROPs)
sn.AddEntProps(ENT_PROPs)

#this dump file will list all proteins in the database with connectivity >0:
sn.AddDumpFile('Proteins from database.tsv',replace_main_dump=True) 
print('Fetching all proteins from the database')
ProteinsOnyGraph = sn.ProcessOQL(flash_dump=True)

import ElsevierAPI.ResnetAPI.PathwayStudioGOQL as GOQL
sn.AddDumpFile("Protein neighbors dump.tsv",replace_main_dump=True) #dump file accumulates all data in one big file
out_dir = 'csv'
counter = 0
for node_id,psObj in ProteinsOnyGraph.nodes(data=True):
    protein_name = psObj['Name'][0]
    sn.ReplaceGOQL(GOQL.ExpandEntity([node_id],SearchByProperties=['id'],ExpandWithRelationTypes=[],ExpandToNeighborTypes=[],direction=''))
    counter += 1
    print('Finding neighbors for \"%s\", node #%d from %d total' % (protein_name,counter,ProteinsOnyGraph.number_of_nodes()))    
    ProteinNeighborsGraph = sn.ProcessOQL()
    protein_neighbors_file = out_dir+'/'+protein_name + '_neighbors.csv'
    PSnx.PrintReferenceView(protein_neighbors_file,REL_PROPs,ENT_PROPs,ProteinNeighborsGraph,col_sep=',')
    sn.Graph.clear()#need to release memory when performing large dumps 


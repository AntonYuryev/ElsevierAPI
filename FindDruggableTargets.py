import time
import networkx as nx
import ElsevierAPI.ResnetAPI.PathwayStudioGOQL as GOQL
import ElsevierAPI.ResnetAPI.ResnetAPISession as ssn
from ElsevierAPI import networx as PSnx 


global_start = time.time()
#SearchEntitiesBy = ['Ileocolitis']
SearchEntitiesBy = ['Inflammatory Bowel Disease']
InputDiseaseNames = ','.join(SearchEntitiesBy)

#specify files used in your API.ssn.DiseaseNetwork.AddGraph to dump graph data in tab-delimited format:
myDir = ''#'D:\\Python\\PS_API'
foutDiseaseSNPs = myDir+"\\Gene variants linked to "+InputDiseaseNames+'.tsv'
foutDiseaseProteins = myDir+"\\Genes with SNPs linked to "+InputDiseaseNames+'.tsv'
foutDrugsForDiseaseProteins = myDir+"\\Druggable targets for "+InputDiseaseNames+'.tsv'

OQLquery = GOQL.ExpandEntity(PropertyValues=SearchEntitiesBy,SearchByProperties=['Name'],ExpandWithRelationTypes=[],ExpandToNeighborTypes=['GeneticVariant'])
sn = ssn.APISession(OQLquery, PSnx)
print("Finding GeneticVariants linked to %s" % InputDiseaseNames)
sn.AddDumpFile(foutDiseaseSNPs,replace_main_dump=True)
sn.ProcessOQL()

SNPIds = set(sn.GetGraphEntityIds(['GeneticVariant']))
sn.GOQLquery = GOQL.ExpandEntity(PropertyValues=SNPIds,SearchByProperties=['id'],ExpandWithRelationTypes=['GeneticChange'],ExpandToNeighborTypes=['Protein'])
print("Finding Proteins containing GeneticVariants linked to %s" % InputDiseaseNames)
sn.AddDumpFile(foutDiseaseProteins,replace_main_dump=True)
sn.ProcessOQL(flash_dump=True)
 

foutDiseasePPI = myDir+"\\PPIs between genes linked to "+InputDiseaseNames+'.tsv'
PPIgraph = sn.GetPPIgraph(foutDiseasePPI)

#calculating centrality
degree_cent = nx.degree_centrality(PPIgraph)
sorted_centrality = sorted(degree_cent.items(), key=lambda kv: kv[1], reverse=True)
PPIids = set([x for x in PPIgraph.nodes()])
IdToNames = sn.GetProperties(PPIids,'Name')
sorted_centrality_byName = list()
for t in sorted_centrality:
    idx = t[0]
    #if idx in IdToNames.keys():
    name = IdToNames[idx][0]
    sorted_centrality_byName.append((name,t[1])) 
    
print(sorted_centrality_byName)


DiseaseProteins = set(sn.GetGraphEntityIds(['Protein']))
sn.GOQLquery = GOQL.GetDrugs(ForTargetsIDlist=DiseaseProteins)
print("Finding Drugs for Proteins containing GeneticVariants linked to %s" % InputDiseaseNames)
sn.AddDumpFile(foutDrugsForDiseaseProteins,replace_main_dump=True)
start_time = time.time()
sn.ProcessOQL(flash_dump=True)                                                                                                                                
DrugCount =    set([x for x,y in sn.Graph.nodes(data=True) if ((sn.Graph.out_degree(x)>0) & (y['ObjTypeName'][0] in ['Small Molecule', 'SmallMol']))])
FoundTargets = set([x for x,y in sn.Graph.nodes(data=True) if ((sn.Graph.in_degree(x)>1) & (y['ObjTypeName'][0] in ['Protein']))])
execution_time = sn.ExecutionTime(start_time)
print("%d drugs for %d proteins linked to %s were retreived in %s ---" % (len(DrugCount),len(FoundTargets), InputDiseaseNames, execution_time))


#find RMC compounds for nondruggable targets
ProteinNoDrugs = DiseaseProteins.difference(FoundTargets)
if len(ProteinNoDrugs) > 0:
    start_time = time.time()
    print('Searching for RMC compounds binding to %d proteins that do not bind any drugs' % len(ProteinNoDrugs))
    sn.GOQLquery = GOQL.GetReaxysSubstances(ForTargetsIDlist=ProteinNoDrugs)
    sn.ProcessOQL()
    execution_time = sn.ExecutionTime(start_time)
    DrugCompoundCount = set([x for x,y in sn.Graph.nodes(data=True) if ((sn.Graph.out_degree(x)>0) & (y['ObjTypeName'][0] in ['Small Molecule', 'SmallMol']))])
    NewTargetCount = set([x for x,y in sn.Graph.nodes(data=True) if ((sn.Graph.in_degree(x)>1) & (y['ObjTypeName'][0] in ['Protein']))])

    RMCcompoundCount = DrugCompoundCount.difference(DrugCount)
    RMCCoumpundTargetCount = NewTargetCount.difference(FoundTargets)
    print("%d lead compounds for %d undruggable proteins linked to %s were retreived in %s ---" % (len(RMCcompoundCount), len(RMCCoumpundTargetCount), InputDiseaseNames, execution_time))
else:
    print('All proteins are druggable. No need to find RMC compounds at this iteration')


#All relations from all queries are in sn.Graph
#To dump sn.Graph in other formats please study https://networkx.org/documentation/stable/reference/readwrite
execution_time = sn.ExecutionTime(global_start)
print("Entire program ran for %s ---" % execution_time)

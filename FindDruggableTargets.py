import time
import networkx as nx
from ElsevierAPI.ResnetAPI.PathwayStudioGOQL import OQL,len
from ElsevierAPI import open_api_session
from ElsevierAPI.ResnetAPI.PSnx2Neo4j import nx2neo4j
from ElsevierAPI.ResnetAPI.PSnx2Neo4j import REL_PROPs, ENT_PROP_Neo4j

ps_api = open_api_session()

global_start = time.time()
SearchEntitiesBy = ['Ileocolitis']
#SearchEntitiesBy = ['Inflammatory Bowel Disease']
# SearchEntitiesBy = ['Friedreich Ataxia']
InputDiseaseNames = ','.join(SearchEntitiesBy)

# specify files used in ps_api.DiseaseNetwork.AddGraph to dump graph data in tab-delimited format:
myDir = ''  # 'D:\\Python\\PS_API\\'
foutDiseaseSNPs = myDir + "Gene variants linked to " + InputDiseaseNames + '.tsv'
foutDiseaseProteins = myDir + "Genes with SNPs linked to " + InputDiseaseNames + '.tsv'
foutDrugsForDiseaseProteins = myDir + "Druggable targets for " + InputDiseaseNames + '.tsv'

ps_api.add_rel_props(REL_PROPs)
ps_api.add_ent_props(ENT_PROP_Neo4j)

print("Finding GeneticVariants linked to %s" % InputDiseaseNames)
ps_api.add_dump_file(foutDiseaseSNPs, replace_main_dump=True)
ps_api.process_oql(
    OQL.expand_entity(PropertyValues=SearchEntitiesBy, SearchByProperties=['Name', 'Alias'], expand_by_rel_types=[],
                       expand2neighbors=['GeneticVariant']))

SNPIds = ps_api.Graph.get_node_ids(['GeneticVariant'])
print("Finding Proteins containing GeneticVariants linked to %s" % InputDiseaseNames)
ps_api.add_dump_file(foutDiseaseProteins, replace_main_dump=True)
ps_api.process_oql(
    OQL.expand_entity(PropertyValues=SNPIds, SearchByProperties=['id'], expand_by_rel_types=['GeneticChange'],
                       expand2neighbors=['Protein']), flush_dump=True)

foutDiseasePPI = myDir + "\\PPIs between genes linked to " + InputDiseaseNames + '.tsv'
PPIgraph = ps_api._get_ppi_graph(foutDiseasePPI)
# calculating centrality

degree_cent = nx.degree_centrality(PPIgraph)
closness = PPIgraph.closeness()
sorted_centrality = sorted(degree_cent.items(), key=lambda kv: kv[1], reverse=True)
PPIids = set(PPIgraph.nodes())
IdToNames = ps_api.Graph.get_properties('Name',PPIids)
sorted_centrality_byName = list()
for t in sorted_centrality:
    idx = t[0]
    name = IdToNames[idx][0]
    sorted_centrality_byName.append((name, t[1]))

print(sorted_centrality_byName)

DiseaseProteins = set(ps_api.Graph.get_node_ids(['Protein']))
print("Finding Drugs for Proteins containing GeneticVariants linked to %s" % InputDiseaseNames)
ps_api.add_dump_file(foutDrugsForDiseaseProteins, replace_main_dump=True)
start_time = time.time()
ps_api.process_oql(OQL.get_drugs(for_targets_with_ids=list(DiseaseProteins)), flush_dump=True)
DrugCount = set([x for x, y in ps_api.Graph.nodes(data=True) if
                 ((ps_api.Graph.out_degree(x) > 0) & (y['ObjTypeName'][0] in ['Small Molecule', 'SmallMol']))])
FoundTargets = set(
    [x for x, y in ps_api.Graph.nodes(data=True) if ((ps_api.Graph.in_degree(x) > 1) & (y['ObjTypeName'][0] in ['Protein']))])
execution_time = ps_api.execution_time(start_time)
print("%d drugs for %d proteins linked to %s were retrieved in %s ---" %
      (len(DrugCount), len(FoundTargets), InputDiseaseNames, execution_time))

# find RMC compounds for non-druggable targets
ProteinNoDrugs = DiseaseProteins.difference(FoundTargets)
if len(ProteinNoDrugs) > 0:
    start_time = time.time()
    print('Searching for RMC compounds binding to %d proteins that do not bind any drugs' % len(ProteinNoDrugs))
    ps_api.process_oql(OQL.get_reaxys_substances(ForTargetsIDlist=list(ProteinNoDrugs)))
    execution_time = ps_api.execution_time(start_time)
    DrugCompoundCount = set([x for x, y in ps_api.Graph.nodes(data=True) if
                             ((ps_api.Graph.out_degree(x) > 0) & (y['ObjTypeName'][0] in ['Small Molecule', 'SmallMol']))])
    NewTargetCount = set([x for x, y in ps_api.Graph.nodes(data=True) if
                          ((ps_api.Graph.in_degree(x) > 1) & (y['ObjTypeName'][0] in ['Protein']))])

    RMCcompoundCount = DrugCompoundCount.difference(DrugCount)
    RMCCoumpundTargetCount = NewTargetCount.difference(FoundTargets)
    print("%d lead compounds for %d undruggable proteins linked to %s were retrieved in %s ---" %
          (len(RMCcompoundCount), len(RMCCoumpundTargetCount), InputDiseaseNames, execution_time))
else:
    print('All proteins are druggable. No need to find RMC compounds at this iteration')

# All relations from all queries are in ps_api.Graph
# To dump ps_api.Graph in other formats please study https://networkx.org/documentation/stable/reference/readwrite
execution_time = ps_api.execution_time(global_start)
print("Entire program ran for %s ---" % execution_time)

# importing entire disease network into local Neo4j instance:
import_start = time.time()
nx2neo = nx2neo4j(uriNeo4j='bolt://localhost:11003', userNeo4j='neo4j', password='')
nx2neo.load_graph2neo4j(ps_api)
print("Graph was imported into Neo4j in %s ---" % ps_api.execution_time(import_start))
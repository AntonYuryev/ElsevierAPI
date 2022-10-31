from ElsevierAPI.ResnetAPI.ResnetGraph import ResnetGraph
import time
from ElsevierAPI import open_api_session
from ElsevierAPI.ResnetAPI.PathwayStudioGOQL import OQL 
from ElsevierAPI.ResnetAPI.NetworkxObjects import BIBLIO_PROPS,PS_ID_TYPES
from  ElsevierAPI.ResnetAPI.rnef2sbgn import make_file_name
import xml.etree.ElementTree as et
from xml.dom import minidom

working_dir = 'D:/Python/CTs/'
input_drug_list = 'ListOfDrugs.txt'
my_drug_list = working_dir + input_drug_list
my_drug_list = None
ps_api = open_api_session()
ps_api.DumpFiles.clear()
ps_api.add_rel_props(BIBLIO_PROPS+PS_ID_TYPES+['Effect','Mechanism'])
ps_api.PageSize = 1000

def retreive_clinical_trials(drugs_name_file=None):
    request_name = 'Find all completed clinical trials'
    oql_query = 'SELECT Relation WHERE objectType = ClinicalTrial AND TrialStatus = Completed AND NeighborOf (SELECT Entity WHERE objectType = SmallMol)'
    if isinstance(drugs_name_file, str):      
        with open(drugs_name_file) as f:
            drug_names = [line.rstrip('\n') for line in f]
        print('Read %s with %d drug names' %(drugs_name_file,len(drug_names)))

        search_by_prop, drug_names_str = OQL.get_search_strings(['Name','Alias'],drug_names)
        oql_query =  oql_query + f' AND NeighborOf (SELECT Entity WHERE (Name,Alias)= ({drug_names_str}))'
        request_name = f'Find clinical trials for drugs in {drugs_name_file}'

    return ps_api.process_oql(oql_query,request_name, debug=False)

search_antagonists = True
start = time.time()
completed_trials = retreive_clinical_trials(my_drug_list)
drug_ids = set(completed_trials.get_entity_ids(['SmallMol']))
print('Found %d clinical trials with %d drugs' % (completed_trials.number_of_edges(), len(drug_ids)))

effect = 'negative' if search_antagonists else 'positive'
oql_query = r'SELECT Relation WHERE NeighborOf (SELECT Entity WHERE id = ({{ids}})) AND objectType = DirectRegulation AND Effect = {effect}'
oql_query  = oql_query.format(effect=effect)
drug_targets = ps_api.iterate_oql(f'{oql_query}', drug_ids)
drug_with_targets = set(drug_targets.get_entity_ids(['SmallMol']))
print('Found %d drugs with targets and clinical trials' % len(drug_with_targets))

drug_counter = 0
out_dir = working_dir+'DrugTargets/'
output_all_targets = False
start = time.time()
for drug in drug_with_targets:
    drug_counter += 1
    target_ids = set(drug_targets.get_neighbors({drug}))
    indication_ids = completed_trials.get_neighbors({drug})

    effect = 'positive' if search_antagonists else 'negative'
    targets_indications = ps_api.connect_nodes(target_ids,indication_ids,with_effect=[effect])
    if targets_indications.size() == 0: continue

    drug_target_indication = ResnetGraph()
    drug_target_indication.add_graph(completed_trials.get_neighbors_graph({drug}))
    drug_target_indication.add_graph(drug_targets.get_neighbors_graph({drug}))
    drug_target_indication.add_graph(targets_indications)
    drug_target_indication.load_references()
    if not output_all_targets:
        targets_without_indications = set([x for x, y in drug_target_indication.nodes(data=True) if
                        ((drug_target_indication.degree(x) == 1) & 
                        (y['ObjTypeName'][0] in ['Protein','Complex','FunctionalClass']))])
        drug_target_indication.remove_nodes_from(targets_without_indications)

    drug_name = drug_target_indication.get_properties('Name',{drug})[drug][0]
    pathway_name = drug_name +' targets in clinical trials'
    pathway_rnef = et.fromstring(ps_api.to_rnef(drug_target_indication))
    pathway_rnef.set('name', pathway_name)
    pathway_rnef.set('type','Pathway')
    batch_xml = et.Element('batch')
    batch_xml.insert(0,pathway_rnef)
    
    TargetCount = set([x for x, y in drug_target_indication.nodes(data=True) if
                        ((drug_target_indication.degree(x) > 1) & 
                        (y['ObjTypeName'][0] in ['Protein','Complex','FunctionalClass']))])

    print('Found %d targets for %s in %d indications. %d in %d drugs processed in %s' % 
        (len(TargetCount),drug_name,len(indication_ids),drug_counter, len(drug_with_targets),ps_api.execution_time(start)))   
    pathway_rnef_str = et.tostring(batch_xml,encoding='utf-8',xml_declaration=True).decode("utf-8")
    pathway_rnef_str = minidom.parseString(pathway_rnef_str).toprettyxml(indent='   ')
    fout_name = make_file_name(pathway_name)
    with open(out_dir+fout_name+'.rnef', mode='w', encoding='utf-8') as rnefout:
        rnefout.write(pathway_rnef_str)

    ps_api.Graph.clear_resnetgraph() # need to release memory when performing large dumps

print('Pathways with drug targets are in %s into %s directory' % (ps_api.execution_time(start),out_dir))
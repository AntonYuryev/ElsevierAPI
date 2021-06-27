from ElsevierAPI.ResnetAPI.ResnetGraph import ResnetGraph
import time
from ElsevierAPI import APIconfig
from ElsevierAPI.ResnetAPI.ResnetAPISession import APISession
from ElsevierAPI.ResnetAPI.NetworkxObjects import REF_PROPS,REF_ID_TYPES
import networkx as nx
import xml.etree.ElementTree as et

search_antagonists = True
start = time.time()
ps_api = APISession(APIconfig['ResnetURL'],APIconfig['PSuserName'],APIconfig['PSpassword'])
ps_api.DumpFiles.clear()
ps_api.add_rel_props(REF_PROPS+REF_ID_TYPES+['Effect'])
ps_api.PageSize = 1000
oql_query = 'SELECT Relation WHERE objectType = ClinicalTrial AND TrialStatus = Completed AND NeighborOf (SELECT Entity WHERE objectType = SmallMol)'
completed_trials = ps_api.process_oql(oql_query, debug=False)
drug_ids = set(completed_trials.get_entity_ids(['SmallMol']))
print('Found %d clinical trials with %d drugs' % (completed_trials.number_of_edges(), len(drug_ids)))

effect = 'negative' if search_antagonists else 'positive'
oql_query = r'SELECT Relation WHERE NeighborOf (SELECT Entity WHERE id = ({{ids}})) AND objectType = DirectRegulation AND Effect = {effect}'
oql_query  = oql_query.format(effect=effect)
drug_targets = ps_api.iterate_oql(f'{oql_query}', drug_ids)
drug_with_targets = set(drug_targets.get_entity_ids(['SmallMol']))
print('Found %d drugs with targets and clinical trials' % len(drug_with_targets))

drug_counter = 0
from xml.dom import minidom
fout_name = 'drug_targets_for_indications'
with open(fout_name+'.rnef', mode='w', encoding='utf-8') as rnefout:
    rnefout.write('<?xml version="1.0" ?>\n<batch>\n')
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
        drug_target_indication.count_references()

        drug_name = drug_target_indication.get_properties({drug},'Name')[drug][0]
        pathway_name = drug_name +' targets in clinical trials'
        rnef_xml = et.fromstring(ps_api.to_rnef(drug_target_indication))
        rnef_xml.set('name', pathway_name)
        rnef_xml.set('type','Pathway')
        
        TargetCount = set([x for x, y in drug_target_indication.nodes(data=True) if
                            ((drug_target_indication.degree(x) > 1) & (y['ObjTypeName'][0] in ['Protein','Complex','FunctionalClass']))])

        print('Found %d targets for %s in %d indications. %d in %d drugs processed in %s' % 
            (len(TargetCount),drug_name,len(indication_ids),drug_counter, len(drug_with_targets),ps_api.execution_time(start)))   
        #batch_xml.append(rnef_xml)
        pathway_rnef = et.tostring(rnef_xml,encoding='utf-8').decode("utf-8")
        pathway_rnef = minidom.parseString(pathway_rnef).toprettyxml(indent='   ')
        rnefout.write(pathway_rnef)
        ps_api.Graph.clear() # need to release memory when performing large dumps

    rnefout.write('</batch>')

print('Your data was downloaded in %s into %s file' % (ps_api.execution_time(start),fout_name))
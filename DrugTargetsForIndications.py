import time
from ElsevierAPI import APIconfig
from ElsevierAPI.ResnetAPI.ResnetAPISession import APISession
from ElsevierAPI.ResnetAPI.NetworkxObjects import REF_PROPS,REF_ID_TYPES
import networkx as nx
import xml.etree.ElementTree as et

start = time.time()
ps_api = APISession(APIconfig['ResnetURL'],APIconfig['PSuserName'],APIconfig['PSpassword'])
ps_api.add_rel_props(REF_PROPS+REF_ID_TYPES)
ps_api.PageSize = 1000
all_completed_trials = ps_api.process_oql('SELECT Relation WHERE objectType = ClinicalTrial AND TrialStatus = Completed AND NeighborOf (SELECT Entity WHERE objectType = SmallMol)')
batch_xml = et.Element('batch')

resnet_count = 0
for drug, indication, ClinicalTrial in all_completed_trials.edges.data('relation'):
    goql_query = 'SELECT Relation WHERE NeighborOf (SELECT Entity WHERE id = {drug}) AND objectType = DirectRegulation AND Effect = negative AND NeighborOf (SELECT Entity WHERE Connected by (SELECT Relation WHERE objectType = (Regulation,QuantitativeChange) AND Effect = positive)) to (SELECT Entity WHERE id = {indication}))'
    links_to_targets = ps_api.process_oql(goql_query.format(drug=drug,indication=indication))
    if links_to_targets is None: # relaxing Effect constraints
        goql_query = 'SELECT Relation WHERE NeighborOf (SELECT Entity WHERE id = {drug}) AND objectType = DirectRegulation AND NeighborOf (SELECT Entity WHERE Connected by (SELECT Relation WHERE objectType = (Regulation,QuantitativeChange)) to (SELECT Entity WHERE id = {indication}))'
        links_to_targets = ps_api.process_oql(goql_query.format(drug=drug,indication=indication))

    if links_to_targets.size() > 0:
        target_ids = links_to_targets.get_entity_ids(['Protein'])
        if len(target_ids) > 0:
            target_disease_links = ps_api.connect_entities(target_ids,['id'],['Protein'],[indication],['id'],['Disease'],ps_api.relProps)
            if target_disease_links.size() > 0:
                all_together = nx.compose(target_disease_links,links_to_targets)
                all_together.add_edge(drug, indication, relation=ClinicalTrial)
                all_together.count_references()

                drug_name = all_together.get_properties({drug},'Name')[drug][0]
                indication_name = all_together.get_properties({indication},'Name')[indication][0]
                pathway_name = drug_name +' targets in '+indication_name

                rnef_xml = et.fromstring(ps_api.to_rnef(all_together))
                rnef_xml.set('name', pathway_name)
                print('Found %d targets for %s in %s' % (len(target_ids),drug_name,indication_name))   
                batch_xml.append(rnef_xml)


from xml.dom import minidom
pathway_rnef = et.tostring(batch_xml,encoding='utf-8',xml_declaration=True).decode("utf-8")
pathway_rnef = minidom.parseString(pathway_rnef).toprettyxml(indent='   ')
fout_name = 'drug_targets_for_indications'
with open(fout_name+'.rnef', mode='w', encoding='utf-8') as f2: f2.write(pathway_rnef)
print('Your data was downloaded in %s into %s file' % (ps_api.execution_time(start),fout_name))
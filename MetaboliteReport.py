from ElsevierAPI import open_api_session
import pandas as pd
import os
import ElsevierAPI.ResnetAPI.PathwayStudioGOQL as OQL
from ElsevierAPI.ResnetAPI.NetworkxObjects import PROTEIN_TYPES
import time


COMMON_METABOLITES={'H2O','PPi', 'ATP','ADP','AMP','Pi','GDP','GTP','NADP+','NADPH+','NAD+','NADH+','acceptor','reduced acceptor',
                    'oxidized acceptor','oxygen','CoA'}


start_time = time.time()
excel_file_name = 'my_metabolites.xlsx'
excel_file_name = 'D:/Python/MDACC/210908_FattyAcidData_LisaMullany.xlsx'
input_excel = pd.read_excel(excel_file_name)
#metabolites names or aliases must be in the first column in Excel file.
input_metabolite_names = []
[input_metabolite_names.append(x) for x in input_excel[input_excel.columns[0]] if x not in input_metabolite_names] # making alias list unique
#input_metabolite_names = input_metabolite_names[0:3]

# ps_api retreives data from the database and loads it into APISession.Graph derived from Networkx:MultiDiGraph
api_config = 'path2apiconfig.json'
api_config = 'D:/Python/ENTELLECT_API/ElsevierAPI/APIconfigMDACC.json'
ps_api = open_api_session(api_config) # specify here path to your APIconfig file. Defaults to ./ElsevierAPI/APIconfig.json
ps_api.add_ent_props(['Alias']) # need to retreive aliases from the database in case input metabolites are found by Alias
ps_api.PageSize = 10000

# dump file contains references for all relations retreived from database
# do not use dump file unless you need to include references into report:
ps_api.DumpFiles.clear()

# retreive all ChemicalReaction linked to metabolites in excel_file_name as ResnetGraph from the database:
my_goql_query = OQL.expand_entity(input_metabolite_names,['Name','Alias'], expand_by_rel_types=['ChemicalReaction'])
request_name ='Retrieve metabolic reactions graph'
reactions_graph = ps_api.process_oql(my_goql_query,request_name)

input_name2objs, objid2input_names = reactions_graph.get_prop2obj_dic('Name', input_metabolite_names)
aliasinput_2objs, objid2input_alias = reactions_graph.get_prop2obj_dic('Alias', input_metabolite_names)
objid2input_names.update(objid2input_alias)
metabolite_ids = list(objid2input_names.keys())
# objid2input_names = {obj_id:[input_names]} - allows for duplicates when mapping by name+alias


# find enzymes linked to ChemicalReactions and retreive their ontology children (proteins)
enzymes = reactions_graph.get_objects(PROTEIN_TYPES)
enzymes_ids = [x['Id'][0] for x in enzymes]
ps_api._get_obj_ids_by_props(enzymes_ids,['Id'],PROTEIN_TYPES get_childs=True) #loading ontology children for enzymes
#retrieved children are stored in ps_api.ID2Children structure

unique_rows= set() # duplicate row filter. Relations can be duplicated due to diffrent Owner,Effect,Mechanism
temp_report_file = 'temp_report.txt'
with open(temp_report_file, 'w', encoding='utf-8') as f:
    f.write('Metabolite input name'+'\t'+'Metabolite name in Database'+'\t'+'direction'+'\t'+'Substrate or Product'+'\t'+'gene'+'\n')
    for regulatorID, targetID, reaction in reactions_graph.edges.data('relation'):
        if regulatorID in metabolite_ids:
            target = reactions_graph._get_node(targetID)
            if target['ObjTypeName'][0] == 'SmallMol' and target['Name'][0] not in COMMON_METABOLITES:
                enzyme_objects = reactions_graph.find_regulators(for_relation=reaction,filter_by=PROTEIN_TYPES)# enzymes are always regulators in ChemicalReaction
                gene_names = ps_api.get_children_props(for_psobjs=enzyme_objects,prop_name='Name')
                gene_names = ','.join(gene_names)
                
                input_metabolite_obj = reactions_graph._get_node(regulatorID)
                metabolite_input_name = objid2input_names[regulatorID]
                metabolite_input_name = ','.join(metabolite_input_name)
                metabolite_db_name = input_metabolite_obj['Name'][0]
                direction = '-->'
                other_metabolite = target['Name'][0]
            else: continue
        elif targetID in metabolite_ids:
            regulator = reactions_graph._get_node(regulatorID)
            if regulator['ObjTypeName'][0] == 'SmallMol' and regulator['Name'][0] not in COMMON_METABOLITES:
                enzyme_objects = reactions_graph.find_regulators(for_relation=reaction,filter_by=PROTEIN_TYPES)
                gene_names = ps_api.get_children_props(enzyme_objects,'Name')
                gene_names = ','.join(gene_names)

                target = reactions_graph._get_node(targetID)
                metabolite_input_name = objid2input_names[targetID]
                metabolite_input_name = ','.join(metabolite_input_name)
                metabolite_db_name = target['Name'][0]
                direction = '<--'
                other_metabolite = regulator['Name'][0]
            else: continue
        else: continue

        tup = (metabolite_db_name,direction,other_metabolite,gene_names)
        if tup not in unique_rows:
            f.write(metabolite_input_name+'\t'+metabolite_db_name+'\t'+direction+'\t'+other_metabolite+'\t'+gene_names+'\n')
            unique_rows.add(tup)

report = pd.read_csv(temp_report_file,sep='\t')
report.sort_values(['Metabolite input name','Substrate or Product'], inplace=True)
report.to_csv('Metabolite report.txt',index=False,sep='\t')
print('Report is generated in %s' % ps_api.execution_time(start_time))
os.remove(temp_report_file)

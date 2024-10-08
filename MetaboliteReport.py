from ElsevierAPI import open_api_session
import pandas as pd
import os
from ElsevierAPI.ResnetAPI.PathwayStudioGOQL import OQL
from ElsevierAPI.ResnetAPI.NetworkxObjects import PROTEIN_TYPES
import time
from ElsevierAPI.ResnetAPI.ResnetGraph import ResnetGraph


COMMON_METABOLITES={'H2O','PPi', 'ATP','ADP','AMP','Pi','GDP','GTP','NADP+','NADPH+','NAD+','NADH+','acceptor','reduced acceptor',
                    'oxidized acceptor','oxygen','CoA'}


start_time = time.time()
excel_file_name = 'my_metabolites.xlsx'
excel_file_name = 'D:/Python/MDACC/211109_LipidData_MinhNguyen.xlsx'
input_excel = pd.read_excel(excel_file_name)
metabolite_column = 1
#metabolites names or aliases must be in the first column in Excel file.
input_metabolite_names = []
[input_metabolite_names.append(x) for x in input_excel[input_excel.columns[metabolite_column]] if x not in input_metabolite_names] # making alias list unique
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
step = 1000
for i in range(0,len(input_metabolite_names),step):
    name_list = input_metabolite_names[i:i+step]
    my_goql_query = OQL.expand_entity(name_list,['Name','Alias'], expand_by_rel_types=['ChemicalReaction'])
    request_name ='Retrieve metabolic reactions graph for {count} metabolites'.format(count=len(name_list))
    ps_api.process_oql(my_goql_query,request_name)

reactions_graph = ps_api.Graph
input_name2objs, objid2input_names = reactions_graph.get_prop2obj_dic('Name', input_metabolite_names)
aliasinput_2objs, objid2input_alias = reactions_graph.get_prop2obj_dic('Alias', input_metabolite_names)
objid2input_names.update(objid2input_alias)
metabolite_ids = list(objid2input_names.keys())
# objid2input_names = {obj_id:[input_names]} - allows for duplicates when mapping by name+alias


# find enzymes linked to ChemicalReactions and retreive their ontology children (proteins)
enzymes = reactions_graph.psobjs_with(only_with_values=PROTEIN_TYPES)
enzymes_ids = [x['Id'][0] for x in enzymes]
ps_api._props2psobj(enzymes_ids,['Id'], get_childs=True, only_obj_types=PROTEIN_TYPES) #loading ontology children for enzymes
#retrieved children are stored in ps_api.ID2Children structure

unique_rows= set() # duplicate row filter. Relations can be duplicated due to diffrent Owner,Effect,Mechanism
temp_report_file = 'temp_report.txt'
with open(temp_report_file, 'w', encoding='utf-8') as f:
    f.write('Metabolite input name'+'\t'+'Metabolite name in Database'+'\t'+'direction'+'\t'+'Substrate or Product'+'\t'+'gene'+'\n')
    for regulatorID, targetID, reaction in reactions_graph.edges.data('relation'):
        if regulatorID in metabolite_ids:
            target = reactions_graph._psobj(targetID)
            if target['ObjTypeName'][0] == 'SmallMol' and target['Name'][0] not in COMMON_METABOLITES:
                enzyme_objects, products = reactions_graph.find_nodes(for_relation=reaction,filter_by=PROTEIN_TYPES)# enzymes are always regulators in ChemicalReaction
                gene_names = ps_api.get_children_props(for_psobjs=enzyme_objects,prop_name='Name')
                gene_names = ','.join(gene_names)
                
                input_metabolite_obj = reactions_graph._psobj(regulatorID)
                metabolite_input_name = objid2input_names[regulatorID]
                metabolite_input_name = ','.join(metabolite_input_name)
                metabolite_db_name = input_metabolite_obj['Name'][0]
                direction = '-->'
                other_metabolite = target['Name'][0]
            else: continue
        elif targetID in metabolite_ids:
            regulator = reactions_graph._psobj(regulatorID)
            if regulator['ObjTypeName'][0] == 'SmallMol' and regulator['Name'][0] not in COMMON_METABOLITES:
                enzyme_objects, products = reactions_graph.find_nodes(for_relation=reaction,filter_by=PROTEIN_TYPES)
                gene_names = ps_api.get_children_props(enzyme_objects,'Name')
                gene_names = ','.join(gene_names)

                target = reactions_graph._psobj(targetID)
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
metabolites_with_reactions = set(report['Metabolite input name'])
print('Metabolite Report is generated in %s' % ps_api.execution_time(start_time))

print('Generating report on missing metabolites')
noreaction_metabolites = list(set(input_metabolite_names).difference(metabolites_with_reactions))
mapped_noreaction_graph = ResnetGraph()
for i in range(0,len(noreaction_metabolites),step):
    name_list = noreaction_metabolites[i:i+step]
    name_list_str = OQL.join_with_quotes(',',name_list)
    oql_query = 'SELECT Entity WHERE (Name,Alias) = ({names})'.format(names=name_list_str)
    request_name ='Retrieve metabolites for {count} input names'.format(count=len(name_list))
    mapped_noreaction_graph.add_graph(ps_api.process_oql(my_goql_query,request_name))

all_input_name2objs, all_objid2input_names = mapped_noreaction_graph.get_prop2obj_dic('Name', input_metabolite_names)
all_aliasinput_2objs, all_objid2input_alias = mapped_noreaction_graph.get_prop2obj_dic('Alias', input_metabolite_names)
all_input_name2objs.update(all_aliasinput_2objs)
with open("mapped metabolites without reactions.txt", 'w', encoding='utf-8') as f:
    f.write('Metabolite input name\tMetabolite name in Database\tMapped metabolite URNs\n')
    for input_name, psobjects in all_input_name2objs.items():
        mapped_urns = [x['URN'][0] for x in psobjects]
        mapped_names = [x['Name'][0] for x in psobjects]
        f.write(input_name+'\t'+','.join(mapped_names)+'\t'+','.join(mapped_urns) +'\n')

umappped_metabolites = set(noreaction_metabolites).difference(set(all_input_name2objs.keys()))
with open("unmapped metabolites.txt", 'w', encoding='utf-8') as f:
    for m in umappped_metabolites:
        f.write(m +'\n')

os.remove(temp_report_file)
print('Report is generated in %s' % ps_api.execution_time(start_time))

from ElsevierAPI import open_api_session
import pandas as pd
import ElsevierAPI.ResnetAPI.PathwayStudioGOQL as OQL
from ElsevierAPI.ResnetAPI.NetworkxObjects import PROTEIN_TYPES
import time

COMMON_METABOLITES={'H2O','PPi', 'ATP','ADP','AMP','Pi','GDP','GTP','NADP+','NADPH+','NAD+','NADH+','acceptor','reduced acceptor',
                    'oxidized acceptor','oxygen','CoA'}


start_time = time.time()
excel_file_name = 'my_metabolites.xlsx'
input_metabolite_names = pd.read_excel(excel_file_name)
input_metabolite_names = list(set(input_metabolite_names[input_metabolite_names.columns[0]]))
#input_metabolite_names = input_metabolite_names[0:3]

# ps_api retreives data from the database and loads it into APISession.Graph derived from Networkx:MultiDiGraph 
ps_api = open_api_session()#specify here path to your APIconfig file. 
#ps_api = open_api_session()
ps_api.add_ent_props(['Alias'])
ps_api.PageSize = 10000
ps_api.DumpFiles.clear() # no dump file is necessary

# retreiving all ChemicalReaction that have metabolites in excel_file_name as ResnetGraph:
my_goql_query = OQL.expand_entity(input_metabolite_names,['Name','Alias'],expand_by_rel_types=['ChemicalReaction'])
request_name ='Find metabolic reactions_graph'
reactions_graph = ps_api.process_oql(my_goql_query,request_name)

# finding database names for metabolites since some of them can be found by Alias:
metabolites_objects  = reactions_graph.get_objects(input_metabolite_names,['Name','Alias'])
metabolite_names = [x['Name'][0] for x in metabolites_objects]

# find enzymes linked to ChemicalReactions and retreive their ontology children (proteins)
enzymes = reactions_graph.get_objects(['FunctionalClass','Complex'])
enzymes_ids = [x['Id'][0] for x in enzymes]
ps_api._get_obj_ids_by_props(enzymes_ids,['Id'],PROTEIN_TYPES) #loading ontology children for enzymes
chemical_reactions_graph = reactions_graph.get_relations(['ChemicalReaction'])

unique_rows= set() # duplicate row filter. Relations can be duplicated due to diffrent Owner,Effect,Mechanism
report = pd.DataFrame(columns=['Metabolite','direction','Substrate or Product','gene'])
# iterating through  Chemicalreactions to find relevant metabolite pairs:
for reaction in chemical_reactions_graph:
    reaction_enzymes = reactions_graph.find_regulators(for_relation=reaction,filter_by=PROTEIN_TYPES)
    enz_names = str()

    if reaction_enzymes:
        enz_names_list = [x['Name'][0] for x in reaction_enzymes]

        enz_ids = [x['Id'][0] for x in reaction_enzymes]
        protein_ids = ps_api.get_children(enz_ids) # finds a
        if protein_ids:
            protein_names = [x['Name'][0] for nodeid,x in ps_api.Graph.nodes(data=True) if x['Id'][0] in protein_ids]
            enz_names = ','.join(enz_names_list+protein_names)
        else:
            enz_names = ','.join(enz_names_list)
        
    regulator_targets_tuples = reaction.get_regulators_targets() #finds all combinations of regulatorID-targetID pairs
    for regulatorID, targetID in regulator_targets_tuples:
        regulator = reactions_graph._get_node(regulatorID)
        if regulator['ObjTypeName'][0] != 'SmallMol': continue #only metabolite pairs are considred for reprot
        target = reactions_graph._get_node(targetID)
        if target['ObjTypeName'][0] != 'SmallMol': continue #only metabolite pairs are considred for reprot
        
        reg_name = regulator['Name'][0]
        targ_name = target['Name'][0]
        #will write input metabolites always into the first column
        if reg_name in metabolite_names and targ_name not in COMMON_METABOLITES:
            column_0 = reg_name
            column_1 = '-->'
            column_2 = targ_name
        elif targ_name in metabolite_names and reg_name not in COMMON_METABOLITES:
            column_0 = targ_name
            column_1 = '<--'
            column_2 = reg_name
        else: continue

        tup = (column_0,column_1,column_2,enz_names)
        if tup not in unique_rows:
            report = report.append({'Metabolite':column_0,'direction':column_1,'Substrate or Product':column_2,'gene':enz_names},ignore_index=True)
        
        unique_rows.add(tup)

report.sort_values(['Metabolite','Substrate or Product'], inplace=True)
report.to_csv('Metabolite report.txt',index=False,sep='\t')
print('Report is generated in %s' % ps_api.execution_time(start_time))
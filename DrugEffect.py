from ElsevierAPI import load_api_config
from ElsevierAPI.ResnetAPI.TargetIndications import RANK_SUGGESTED_INDICATIONS, PREDICT_RANK_INDICATIONS
from ElsevierAPI.ResnetAPI.RepurposeDrug import RepurposeDrug, INHIBIT,ACTIVATE
from ElsevierAPI.ResnetAPI.RepurposeDrug import AGONIST, ANTAGONIST
from ElsevierAPI.pandas.panda_tricks import ExcelWriter
import time
from datetime import datetime
from contextlib import redirect_stdout


def do_the_job(dcp:RepurposeDrug):
    global_start = time.time()
    targets = ['CNR1','CNR2','GPR55','GPR119','GPR18']

    dcp.params['indication_types'] = ['Disease']
    dcp.params['drug_effect'] = INHIBIT
    report_path = dcp.report_path()
    target_report = ExcelWriter(report_path, engine='xlsxwriter')
    raw_data_path = dcp.report_path('_raw_data.xlsx')
    raw_data_cache = ExcelWriter(raw_data_path, engine='xlsxwriter')
    for target_name in targets:
        dcp.make_report([target_name])
        dcp.add2writer(target_report)
        dcp.addraw2writer(raw_data_cache,dcp._worksheet_prefix())
        dcp.Graph.remove_edges4psobjs(target_name) #for ps_bibliography specificity

    other_effects_graph_disease = dcp.other_effects()
    other_indications_dis = other_effects_graph_disease.snippets2df(df_name='Possbl.Diseases')
    dcp.add2report(other_indications_dis)
    dcp.add2writer(target_report,dcp._worksheet_prefix(),['Possbl.Diseases'])
    target_report.save()
    raw_data_cache.save()

    # to find biological processes inhibited by input compound:
    dcp.params['indication_types'] = ['CellProcess']
    dcp.params['drug_effect'] = INHIBIT
    report_path = dcp.report_path()
    target_report = ExcelWriter(report_path, engine='xlsxwriter')
    raw_data_path = dcp.report_path('_raw_data.xlsx')
    raw_data_cache = ExcelWriter(raw_data_path, engine='xlsxwriter')
    for target_name in targets:
        dcp.make_report([target_name])
        dcp.add2writer(target_report,dcp._worksheet_prefix())
        dcp.addraw2writer(raw_data_cache,dcp._worksheet_prefix())
        dcp.Graph.remove_edges4psobjs(target_name) #for ps_bibliography specificity
    target_report.save()
    raw_data_cache.save()

    # to find biological processes activated by input compound:
    dcp.params['indication_types'] = ['CellProcess']
    dcp.params['drug_effect'] = ACTIVATE
    report_path = dcp.report_path()
    target_report = ExcelWriter(report_path, engine='xlsxwriter')
    raw_data_path = dcp.report_path('_raw_data.xlsx')
    raw_data_cache = ExcelWriter(raw_data_path, engine='xlsxwriter')
    for target_name in targets:
        dcp.make_report([target_name])
        dcp.add2writer(target_report,dcp._worksheet_prefix())
        dcp.addraw2writer(raw_data_cache,dcp._worksheet_prefix())
        dcp.Graph.remove_edges4psobjs(target_name)

    other_effects_graph_cellproc = dcp.other_effects()
    other_indications_cellproc = other_effects_graph_cellproc.snippets2df(df_name='Possbl.CellProcess')
    dcp.add2report(other_indications_cellproc)
    dcp.add2writer(target_report,dcp._worksheet_prefix(),['Possbl.CellProcess'])
    target_report.save()
    raw_data_cache.save()

    dcp.params['drug_effect'] = ACTIVATE
    dcp.params['indication_types'] = ['Disease']
    report_path = dcp.report_path()
    target_report = ExcelWriter(report_path, engine='xlsxwriter')
    raw_data_path = dcp.report_path('_raw_data.xlsx')
    raw_data_cache = ExcelWriter(raw_data_path, engine='xlsxwriter')
    for target_name in targets:
        dcp.make_report([target_name])
        dcp.add2writer(target_report,dcp._worksheet_prefix())
        dcp.addraw2writer(raw_data_cache,dcp._worksheet_prefix())
        dcp.Graph.remove_edges4psobjs(target_name)
    target_report.save()
    raw_data_cache.save()

    print('Report was generated in %s' % dcp.execution_time(global_start))


if __name__ == "__main__":
        instructions = """
        'input_compound' - required name of the drug for repurposing
        'similars'  - optional list of similar drugs that have same mechanism of action (i.e. same target_name(s))
        'indication_types' - required combination of Disease or CellProcess or Virus or Pathogen 
        'drug_effect' - required drug effect on indication: 
                INHIBIT  for Disease to find drig indications,
                ACTIVATE for Disease to find toxicities and side-effects
                INHIBIT  for CellProcess to find cell processes inhibited by input_compound or similars
                ACTIVATE for CellProcess to find cell processes activated by input_compound or similars
        'target_names' - optional. script uses all targets in 'target_names' to find or rank indications. 
                To find indications linked only to single target_name 'target_names' must contain only one target_name name
        'mode_of_action' - required if 'target_names' is specified. specifies effect of input_compound on target_name(s): AGONIST or ANTAGONIST
        'target_type' - optional. Does not influence ranking and used for GOQL query optimization
        'to_inhibit' - internal parameter calculated by set_drug(). Specifies drug effect on targets in 'target_names' 
        'partner_names' - optional. explicit list of endogenous ligands for drug targets. 
                If 'partner_names' is not specified script attempts to find endogenous ligands with object type = Protein. 
                Use 'partner_names' when endogenous ligands for targets in 'target_names' are metabolites
        'partner_class' - required if 'partner_names' is specified
                'partner_class' is used for report headers and messaging. It can be any string e.g. 'Metabolite ligands'
        'strict_mode' - required. must be RANK_SUGGESTED_INDICATIONS or PREDICT_RANK_INDICATIONS:
                RANK_SUGGESTED_INDICATIONS - ranks only indications that exist in knowledge graph for 'input_compound' and 'similars',
                i.e. indications suggested for 'input_compound' or 'smilars' in the literature
                if 'target_names' is specified only indications that exist for both 'input_compound' and 'target_names' are ranked
                PREDICT_RANK_INDICATIONS additionaly ranks indication predicted for 'input_compound'. Predictions are made from:
                    - indications suggested for drugs in 'similars'
                    - idications linked to drug targets from 'target_names'
        'pathway_name_must_include_target' - specifies how to construct downstream signaling pathways for targets from 'target_names'
                if False all curated pathways containg at least one target_name from 'target_names' will be used. 
                set 'pathway_name_must_include_target' to True for targets contained in many curated pathways to increase speed and specificity
        'data_dir' - output directory for report. 
                Report is Excel wokbook containing:
                    - reference count supporting various ways (=scores) each indication is linked to input_compound
                    - no more than five most relevant references linking input_compound and indication
                    - ranking of each indication
                    - high level parent ontology category from MAP2ONTOLOGY for each indication
                    - possible indications that could not be ranked due to missing effect sign in a link to 'input_compound'
     """

        parameters = {
                #'input_compound' : 'cannabidiol',#'2-arachidonoylglycerol', #'anandamide', # , '9-tetrahydrocannabinol'
                #'similars' : ['cannabidivarin', 'Cannabidiolic acid', 'Cannabielsoin'],
                'input_compound' : 'tetrahydrocannabinol', 
                'similars' : ['delta 8-THC','9-tetrahydrocannabinol', 'THC-C4','tetrahydrocannabinolic acid', '11-hydroxy-delta 9-tetrahydrocannabinol'],
                
                'drug_effect': INHIBIT, # for indications
                #'drug_effect': ACTIVATE, # for toxicities

                #'mode_of_action': ANTAGONIST, # for CBD
                'mode_of_action': AGONIST, # for THC

                'indication_types': ['Disease'], #['CellProcess','Virus']
                'target_names':[],
                'target_type':'Protein',
                'to_inhibit':True, # parameter will be recalculated by RepurposeDrug::set_drug()
                'partner_names':['endocannabinoid','anandamide','2-arachidonoylglycerol','oleoylethanolamide','virodhamine',
                                'N-oleoyldopamine','palmitoylethanolamide','N-arachidonoyl dopamine','N-arachidonoylglycine',
                                'eicosapentaenoyl ethanolamide','noladin ether'],
                'partner_class':'Metabolite ligand', # use it only if partner_names not empty
                'strict_mode':RANK_SUGGESTED_INDICATIONS, # PREDICT_RANK_INDICATIONS #
                'pathway_name_must_include_target':True,
                'data_dir':'D:/Python/PMI/',
                'debug':False
                }
        
        APIconfig = load_api_config()
        dcp = RepurposeDrug(APIconfig,**parameters)
        print('Script was started at %s' % datetime.now())
        if parameters['debug']:
                do_the_job(dcp)
        else:
                log_path = dcp.report_path('.log')
                print('Runtime messages will be written to "%s"' % log_path)
                with open(log_path, 'w') as fout:
                        with redirect_stdout(fout):
                                do_the_job(dcp)
        print('Script finished at %s' % datetime.now())
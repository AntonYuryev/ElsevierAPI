from ElsevierAPI import load_api_config
from ElsevierAPI.ResnetAPI.TargetIndications import RANK_SUGGESTED_INDICATIONS, PREDICT_RANK_INDICATIONS
from ElsevierAPI.ResnetAPI.RepurposeDrug import RepurposeDrug, INHIBIT,ACTIVATE
from ElsevierAPI.ResnetAPI.RepurposeDrug import AGONIST, ANTAGONIST
from ElsevierAPI.pandas.panda_tricks import ExcelWriter
import time


DISEASE_MAP2ONTOLOGY = [
                'psychiatric disorder',
                'neurological disorder',
                'musculoskeletal disorder',
                'digestive system disease',
                'endocrine disorder',
                'urogenital disease',
                'skin and connective tissue disease',
                'cardiovascular disease',
                'head and neck disease',
                'respiratory disease',
                'genetic and familial disorder',
                'sclerosis',
                'metabolic disorder',
                'neoplasm',
                'inflammation',
                'degeneration',
                'infection',
                'complications',
                'toxicity',
                'wounds and injuries',
                'pain',
                'nutritional and metabolic diseases',
                'body weight disorder',
                'motor dysfunction',
                'neurocognitive disorder',
                'immunopathology',
                'nausea and vomiting',
                'infertility',
                'sleep disorder',
                'appetite disturbance'
                ]

CELLPROCESS_MAP2ONTOLOGY =[
                            'behavior',
                            'cell proliferation',
                            'apoptosis',
                            'memory',
                            'appetite',
                            'cell death',
                            'locomotion',
                            'cognition',
                            'synaptic transmission',
                            'neurotransmission',
                            'vasodilation',
                            'perception of pain',
                            'long-term synaptic depression',
                            'synaptic plasticity',
                            'learning',
                            'working memory',
                            'cell differentiation',
                            'neuronal activity',
                            'sleep',
                            'long-term synaptic potentiation'
                        ]


if __name__ == "__main__":
    instructions = """
        'input_compound' - required name of the drug for repurposing
        'similars'  - optional list of similar drugs that have same mechanism of action (i.e. same target(s))
        'indication_types' - required combination of Disease or CellProcess or Virus or Pathogen 
        'drug_effect' - required drug effect on indication: 
                INHIBIT  for Disease to find drig indications,
                ACTIVATE for Disease to find toxicities and side-effects
                INHIBIT  for CellProcess to find cell processes inhibited by input_compound or similars
                ACTIVATE for CellProcess to find cell processes activated by input_compound or similars
        'target_names' - optional. script uses all targets in 'target_names' to find or rank indications. 
                To find indications linked only to single target 'target_names' must contain only one target name
        'mode_of_action' - required if 'target_names' is specified. specifies effect of input_compound on target(s): AGONIST or ANTAGONIST
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
                if False all curated pathways containg at least one target from 'target_names' will be used. 
                set 'pathway_name_must_include_target' to True for targets contained in many curated pathways to increase speed and specificity
        'data_dir' - output directory for report. 
                Report is Excel wokbook containing:
                    - reference count supporting various ways (=scores) each indication is linked to input_compound
                    - no more than five most relevant references linking input_compound and indication
                    - ranking of each indication
                    - high level parent ontology category from MAP2ONTOLOGY for each indication
                    - possible indications that could not be ranked due to missing effect sign in a link to 'input_compound'
     """

    global_start = time.time()
    APIconfig = load_api_config()
    parameters = {
                #'input_compound' : 'cannabidiol',#'2-arachidonoylglycerol', #'anandamide', # , '9-tetrahydrocannabinol'
                #'similars' : ['cannabidivarin', 'Cannabidiolic acid', 'Cannabielsoin'],
                'input_compound' : 'tetrahydrocannabinol', 
                'similars' : ['delta 8-THC', '9-tetrahydrocannabinol', 'THC-C4','tetrahydrocannabinolic acid', '11-hydroxy-delta 9-tetrahydrocannabinol'],
                'indication_types': ['Disease'], #['CellProcess','Virus']
                'drug_effect': INHIBIT,
                'mode_of_action': ANTAGONIST,
                'target_names':[],
                'target_type':'Protein',
                'to_inhibit':True, # parameter will be recalculated by RepurposeDrug::set_drug()
                'partner_names':['endocannabinoid','anandamide','2-arachidonoylglycerol','oleoylethanolamide','virodhamine',
                                'N-oleoyldopamine','palmitoylethanolamide','N-arachidonoyl dopamine','N-arachidonoylglycine',
                                'eicosapentaenoyl ethanolamide','noladin ether'],
                'partner_class':'Metabolite ligand', # use it only if partner_names not empty
                'strict_mode':RANK_SUGGESTED_INDICATIONS, # PREDICT_RANK_INDICATIONS #
                'pathway_name_must_include_target':True,
                'data_dir':'D:/Python/PMI/'
                }

    targets = ['CNR1','CNR2','GPR55','GPR119','GPR18']
    report_name = parameters['data_dir']+parameters['input_compound']
    target_report = ExcelWriter(report_name+' effects,indications.xlsx', engine='xlsxwriter')
    raw_data_cache = ExcelWriter(report_name+'_raw_data.xlsx', engine='xlsxwriter')
    
    dcp = RepurposeDrug(APIconfig,parameters)
    dcp.load_ontology(DISEASE_MAP2ONTOLOGY+CELLPROCESS_MAP2ONTOLOGY)
    dcp.param['mode_of_action'] = AGONIST #for THC; ANTAGONIST #for CBD; 

    # to find disease indications
    dcp.param['drug_effect'] = INHIBIT
    dcp.param['indication_types'] = ['Disease']
    dcp.load_drug_indications()
    for target in targets:
         dcp.load_target_indications([target])
         if dcp.perform_semantic_search():
            dcp.add2writer(target_report,dcp._worksheet_prefix())
            dcp.addraw2writer(raw_data_cache,dcp._worksheet_prefix())
            dcp._clear_()

    # to find biological processes inhibited by input compound
    dcp.param['indication_types'] = ['CellProcess']
    dcp.load_drug_indications()
    for target in targets:
        dcp.load_target_indications([target])
        if dcp.perform_semantic_search():
            dcp.add2writer(target_report,dcp._worksheet_prefix())
            dcp.addraw2writer(raw_data_cache,dcp._worksheet_prefix())
            dcp._clear_()

    # to find biological processes activated by input compound
    dcp.param['drug_effect'] = ACTIVATE
    dcp.load_drug_indications()
    for target in targets:
        dcp.load_target_indications([target])
        if dcp.perform_semantic_search():
            dcp.add2writer(target_report,dcp._worksheet_prefix())
            dcp.addraw2writer(raw_data_cache,dcp._worksheet_prefix())
            dcp._clear_()

    other_processes = dcp.rn2pd(dcp.other_effects(),dcp.param['input_compound'])
    other_processes.to_excel(target_report, sheet_name='Possbl.CellProcess', index=False)

    dcp.param['indication_types'] = ['Disease']
    other_indications = dcp.rn2pd(dcp.other_effects(),dcp.param['input_compound'])
    other_indications.to_excel(target_report, sheet_name='Possbl.Diseases', index=False)

    target_report.save()
    raw_data_cache.save()
    print('Report was generated in %s' % dcp.execution_time(global_start))

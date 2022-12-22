from ElsevierAPI import load_api_config
from ElsevierAPI.pandas.panda_tricks import pd
from ElsevierAPI.ResnetAPI.TargetIndications import Indications4targets


if __name__ == "__main__":
    instructions = """
        'indication_types' - required combination of Disease or CellProcess or Virus or Pathogen 
        'target_names' - required. script uses all targets in 'target_names' to find or rank indications. 
                To find indications linked only to single target 'target_names' must contain only one target name
        'target_type' - optional. Does not influence ranking and used for GOQL query optimization
        'to_inhibit' - required. Specifies how targets fro 'target_names' are modulated (by a drug). 
                    if True target should be inhibited by a drug
                    if False target should be activated by a drug
        'partner_names' - optional. explicit list of endogenous ligands for drug targets. 
                If 'partner_names' is not specified script attempts to find endogenous ligands with object type = Protein. 
                Use 'partner_names' when endogenous ligands for targets in 'target_names' are metabolites
        'partner_class' - required if 'partner_names' is specified
                'partner_class' is used for report headers and messaging. It can be any string e.g. 'Metabolite ligands'
        'strict_mode' - required. must be RANK_SUGGESTED_INDICATIONS or PREDICT_RANK_INDICATIONS:
                RANK_SUGGESTED_INDICATIONS - ranks only indications that exist in the knowledge graph for input targets,
                i.e. indications suggeted for input targets in the literature
                PREDICT_RANK_INDICATIONS additionaly ranks indication predicted for input targets. Predictions are made from:
                    - indications that have input targets as clinical biomarkers
        'pathway_name_must_include_target' - specifies how to construct downstream signaling pathways for targets from 'target_names'
                if False all curated pathways containg at least one target from 'target_names' will be used. 
                set 'pathway_name_must_include_target' to True for targets contained in many curated pathways to increase speed and specificity
        'data_dir' - output directory for report. 
                Report is Excel wokbook containing:
                    - reference count supporting various ways (=scores) each indication is linked to targets in 'target_names'
                    - no more than five most relevant references linking targets in 'target_names' and indication
                    - ranking of each indication
                    - high level parent ontology category from MAP2ONTOLOGY for each indication
                    - possible indications that could not be ranked due to missing effect sign in a link to targets in 'target_names'
     """

    APIconfig = load_api_config()
    parameters = {'partner_names':['anandamide','endocannabinoid','2-arachidonoylglycerol'],
    # if partner_names is empty script will try finding Lgands for Receptor targets and Receptors for Ligand targets
                 'partner_class':'Metabolite ligand', # use it only if partner_names not empty
                 'indication_types': ['CellProcess'], #['Disease','Virus']
                 'target_names':['GPR18'],
                 'target_type':'Protein',
                 'to_inhibit':True,
                 'pathway_name_must_include_target':True,
                 'strict_mode':True,
                 'data_dir':'D:/Python/PMI/'
    }

    targets = ['GPR18','GPR119','GPR55','CNR2','CNR1']
    indications = ['CellProcess','Disease']

    path = parameters['data_dir']

    for target in targets:
        parameters['target_names'] = [target]
        target_report = pd.ExcelWriter(path+target+'.xlsx', engine='xlsxwriter')
        raw_data_cache = pd.ExcelWriter(path+target+'raw_data.xlsx', engine='xlsxwriter')

        for indication in indications:
            parameters['indication_types'] = [indication]
            parameters['to_inhibit'] = True
            rd = Indications4targets(APIconfig,parameters)
            rd.make_report()
            rd.columns2drop += [rd.__resnet_name__, rd.__mapped_by__]
            rd.add2writer(target_report)
            rd.addraw2writer(raw_data_cache)

            parameters['to_inhibit'] = False
            rd = Indications4targets(APIconfig,parameters)
            rd.make_report()
            rd.columns2drop += [rd.__resnet_name__, rd.__mapped_by__]
            rd.add2writer(target_report,rd)
            rd.addraw2writer(raw_data_cache,rd)

        target_report.save()
        raw_data_cache.save()


            

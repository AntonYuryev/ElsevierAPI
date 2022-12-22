from ElsevierAPI import load_api_config
import pandas as pd
from ElsevierAPI.ResnetAPI.TargetIndications import Indications4targets
from ElsevierAPI.ResnetAPI.TargetIndications import RANK_SUGGESTED_INDICATIONS,PREDICT_RANK_INDICATIONS
from contextlib import redirect_stdout

if __name__ == "__main__":
    instructions = """
        'indication_types' - required combination of Disease or CellProcess or Virus or Pathogen 
        'target_names' - required. script uses all targets in 'target_names' to find or rank indications. 
                To find indications linked only to single target 'target_names' must contain only one target name
        'target_type' - optional. Does not influence ranking and used for GOQL query optimization
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
    parameters = {
                'target_names':['CSF3R'],
                'target_type':'Protein',
                'partner_names':['CSF3'],
    # if partner_names is empty script will try finding Lgands for Receptor targets and Receptors for Ligand targets
                 'partner_class':'Ligand', # use it only if partner_names not empty
                 'indication_types': ['Disease'],
                 #'to_inhibit':True,
                 'pathway_name_must_include_target':True,
                 'strict_mode': PREDICT_RANK_INDICATIONS, #RANK_SUGGESTED_INDICATIONS, #
                 'data_dir':'D:/Python/CHL/',
                 'add_bibliography' : True
                }

rd = Indications4targets(APIconfig,parameters)
report_name = rd._get_report_name()
target_report = pd.ExcelWriter(report_name+'.xlsx', engine='xlsxwriter')
raw_data_cache = pd.ExcelWriter(report_name+'_raw_data.xlsx', engine='xlsxwriter')

log_name = 'Indications for {}.log'.format(','.join(rd.param['target_names']))
with open(rd.param['data_dir']+log_name, 'w') as fout:
        with redirect_stdout(fout):
                rd.make_report()
                rd.columns2drop += [rd.__resnet_name__, rd.__mapped_by__]
                rd.add2writer(target_report)
                rd.addraw2writer(raw_data_cache)
                target_report.save()
                raw_data_cache.save()


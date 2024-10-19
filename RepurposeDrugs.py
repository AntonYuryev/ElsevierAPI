import time, json
from pandas import ExcelWriter
from ENTELLECT_API.ElsevierAPI.utils import Tee, execution_time
from ENTELLECT_API.ElsevierAPI.ResnetAPI.RepurposeDrug import RepurposeDrug


def do_the_job(dcp:RepurposeDrug):
        dcp.make_report()

        report_path = dcp.report_path()
        report = ExcelWriter(report_path, engine='xlsxwriter')
        dcp.add2writer(report)
        report.close()

        raw_report_path = dcp.report_path('_raw_data','.xlsx')
        raw_data_cache = ExcelWriter(raw_report_path, engine='xlsxwriter')
        dcp.addraw2writer(raw_data_cache)
        raw_data_cache.close()

        print(f'Report was generated in {execution_time(global_start)}')
        print(f'Report is in "{report_path}" file')
        dcp.clear()


if __name__ == "__main__":
        instructions = """
        'input_compound' - required name of the drug for repurposing
        'similars'  - optional list of similar drugs that have same mechanism of action (i.e. same target(s))
        'indication_types' - required. Any combination of [Disease, CellProcess, Virus, Pathogen]
        'drug_effect' - required drug effect on indication: 
                INHIBIT  for Disease to find drug indications,
                ACTIVATE for Disease to find drug toxicities or side-effects
                INHIBIT  for [CellProcess, Virus, Pathogen] to find concepts inhibited by input_compound or similars
                ACTIVATE for [CellProcess, Virus, Pathogen] to find concepts activated by input_compound or similars
        'targets' - optional. dict{mode_of_action:[target_names]}. script uses all targets in 'target_names' to find or rank indications. 
                To find indications linked only to single target 'target_names' must contain only one target name
                if 'target_names' is empty script will attempt to find targets in Reaxys and then in Pathway Studio if no targets found in Reaxys
                'mode_of_action' - specifies effect of input_compound on target(s): AGONIST or ANTAGONIST
        'partners' - optional. dict{target_name:partner_class:[partner_names]}. explicit list of endogenous ligands for drug targets. 
                If 'partner_names' is not specified script attempts to find endogenous ligands with object type = Protein. 
                Use 'partner_names' when endogenous ligands for targets in 'target_names' are metabolites
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
        parameters_list = list(json.load(open(drug_parameters, "r")))

        for parameters in parameters_list:
                assert(isinstance(parameters,dict))
                dcp = RepurposeDrug(**parameters)
                if parameters.pop('skip',False): continue
                else:
                        with Tee(dcp.report_path(extension='.log')):
                                do_the_job(dcp)


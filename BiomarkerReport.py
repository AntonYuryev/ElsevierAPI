from ElsevierAPI import load_api_config
from ElsevierAPI.ResnetAPI.Biomarkers import BiomarkerReport,SOLUBLE,GENETIC,QUANTITATIVE
import time
from contextlib import redirect_stdout

DISEASE_NAME = 'Pulmonary Hypertension'# 'atrial fibrillation'# 'Myocardial Infarction'# 'cardiovascular disease'#'diabetes mellitus'   #'fibrosis'

if __name__ == "__main__":
    params = {
                'disease': DISEASE_NAME,
                'print_references':False,
                'biomarker_type': GENETIC,# SOLUBLE,#  QUANTITATIVE,# 
                'journal_filter_fname':'', 
                'add_specificity4diseases': ['cardiovascular disease'],
                # there is no scientific rational to calculate biomarker specificity for genetic biomarkers
                'data_dir' : 'D:/Python/PMI/'+DISEASE_NAME+'/',
                'debug':False,
                'print_rdf':False
                }

    start_time = time.time()
    api_config = str()
    APIconfig = load_api_config(api_config)
    
    bm = BiomarkerReport(APIconfig,params)
    bm.flush_dump()

    if bm.params['debug']:
        bm.print_report()
    else:
        log_path = bm.report_path('.log')
        with open(log_path, 'w') as fout:
            with redirect_stdout(fout):
                bm.print_report()
        
    print('Biomarkers for %s was found in %s' % (bm.params['disease'], bm.execution_time(start_time)))

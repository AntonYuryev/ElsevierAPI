from ElsevierAPI import load_api_config
from ElsevierAPI.ResnetAPI.Drugs4Disease import Drugs4Targets
from contextlib import redirect_stdout
import pandas as pd
from datetime import datetime

def do_the_job(dt:Drugs4Targets):
    dt.make_report()
    report_name = dt._get_report_name()
    target_report = pd.ExcelWriter(dt.param['data_dir']+report_name+'.xlsx', engine='xlsxwriter')
    raw_data_cache = pd.ExcelWriter(dt.param['data_dir']+report_name+'_raw_data.xlsx', engine='xlsxwriter')

    dt.add2writer(target_report)
    dt.addraw2writer(raw_data_cache)

    raw_data_cache.save()
    target_report.save()
    print('Report is in %s file' % report_name)

if __name__ == "__main__":
    parameters = {
                'disease': ['Carney Complex'] ,
                'symptoms':['primary pigmented nodular adrenocortical disease','pituitary adenoma',
                'fibroadenoma',"Cushing syndrome",'multiple endocrine neoplasia syndrome','Neurilemmoma','myxoma',
                'Lentigo','Melanosis'],
                'processes':[],
                'clinical_parameters':[],
                'pathway_name_must_include_target':True,
                'strict_mode':False,
                'data_dir':'D:/Python/MATRIX/', # needs slash at the end
                'add_bibliography' : True,
                'target_types' : ['Protein','FunctionalClass','Complex'], #['Protein'], #,
                'pathway_folders':[],
                'pathways': [],
                'debug': False,
                'consistency_correction4target_rank':True
                }

print('Script was started at %s' % datetime.now())
dt = Drugs4Targets(load_api_config(), params=parameters)
if dt.param['debug']:
    do_the_job(dt)
else:
    log_name = f'Drugs for {dt._disease2str()}.log'
    with open(dt.param['data_dir']+log_name, 'w') as fout:
        with redirect_stdout(fout):
            do_the_job(dt)

print('Script finished at %s' % datetime.now())

from ElsevierAPI import load_api_config,execution_time
from ElsevierAPI.ResnetAPI.SNEA import SNEA,test_run
from contextlib import redirect_stdout
import time

def do_the_job(snea:SNEA):
    start = time.time()
    snea.expression_regulators()
    snea.activity_df()
    snea.make_drugs_df()
    snea.report()
    print('Job was done in %s' % execution_time(start))


if __name__ == "__main__":
    APIconfig_file = ''
    debug = False

    if debug:
        test_run(load_api_config(APIconfig_file),fast=True)
    else:
        experiment_name = 'PNOC003vsGSE120046'
        snea = SNEA(load_api_config(APIconfig_file),experiment_name,find_drugs=True)
        snea.set_dir('D:/Python/PBTA/')
        log_name = snea.report_path('log')
        print('Runtime messages are written to "%s"'%log_name)
        with open(log_name, 'w') as fout:
            with redirect_stdout(fout):
                do_the_job(snea)



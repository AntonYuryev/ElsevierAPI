
from ElsevierAPI import load_api_config,execution_time
from ElsevierAPI.ResnetAPI.SNEA import SNEA,test_run
from contextlib import redirect_stdout
import time

parameters = {
        'sample_ids' : [1], #for DESeq2 output that has log2FC in the first column and pValue in the second
        'has_pvalue' : True,
        'find_drugs' : True,
        'data_dir' : 'D:/Python/NEO7/AssafAman',
        'no_mess' : False,
        'cache_dir' : 'D:/Python/ENTELLECT_API/ElsevierAPI/ResnetAPI/__pscache__/',
        'limitDrugs2expression_regulators':False,
        'max_diffexp_pval' : 0.05
    }


def do_the_job(snea:SNEA):
    start = time.time()
    snea.expression_regulators()
    snea.activity_df()
    snea.make_drugs_df()
    snea.report()
    print('Job was done in %s' % execution_time(start))


if __name__ == "__main__":
    APIconfig_file = ''
    debug = True
    use_all_samples_without_outliers=False

    if debug:
        fast = not use_all_samples_without_outliers
        test_run(load_api_config(APIconfig_file),fast=fast)
    else:
        experiment_name = 'PNOC003vsGSE120046'
        snea = SNEA(experiment_name,load_api_config(APIconfig_file),find_drugs=True)
        snea.set_dir('D:/Python/PBTA/')
        log_path = snea.report_path('log')
        print('Runtime messages are written to "%s"'%log_path)
        with open(log_path, 'a') as fout:
            # "append" mode is necessary to avoid overriding messages from SNEA initialization
            with redirect_stdout(fout):
                do_the_job(snea)



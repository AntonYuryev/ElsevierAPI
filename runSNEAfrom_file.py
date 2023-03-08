from ElsevierAPI.ResnetAPI.SNEA import SNEA
from contextlib import redirect_stdout
import time


parameters = {
        'sample_ids' : [1], #for DESeq2 output that has log2FC in the first column and pValue in the second
        'has_pvalue' : True,
        'find_drugs' : True,
        'data_dir' : 'D:/Python/NEO7/',
        'no_mess' : False
    }


def do_the_job(snea:SNEA):
    start = time.time()
    snea.expression_regulators()
    snea.activity_df()
    snea.make_drugs_df()
    snea.report()
    print('Job was done in %s' % SNEA.execution_time(start))


if __name__ == "__main__":
    path2experiment = parameters['data_dir']+'GTExvsNEO7US206MAGLB015.DESeq2.txt'
    snea = SNEA.from_files(path2experiment,**parameters)
    if parameters['no_mess']:
        log_name = snea.report_path('log')
        print('Runtime messages are written to "%s"'%log_name)
        with open(log_name, 'a') as fout:
            with redirect_stdout(fout):
                do_the_job(snea)
    else:
        do_the_job(snea)
    




from ElsevierAPI import load_api_config
from ElsevierAPI.ResnetAPI.SNEA import SNEA

#exp = Experiment('PNOC003 transcription regulators activity',load_api_config() )
APIcofig = 'D:/Python/ENTELLECT_API/ElsevierAPI/APIconfigTeva.json'
#APIcofig = ''
snea = SNEA(load_api_config(APIcofig),'Drug pipeline testing') #'PNOC003vsGSE120046') #
snea.expression_regulators()
snea.report()
#snea.test_run(fast=False)

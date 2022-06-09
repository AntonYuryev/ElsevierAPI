
from ElsevierAPI import load_api_config
from ElsevierAPI.ResnetAPI.SNEA import SNEA

APIcofig = ''
snea = SNEA(load_api_config(APIcofig),'Drug pipeline testing')
snea.expression_regulators()
snea.report()
#snea.test_run(fast=False)

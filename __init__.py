import ElsevierAPI.ResnetAPI.PathwayStudioZeepAPI as PSzeepAPI
import ElsevierAPI.ResnetAPI.ZeepToNetworkx as ZeeptoNX
import json

con_file = open("ElsevierAPI/APIconfig.json")#file with your API keys and API URLs
APIconfig = json.load(con_file)
con_file.close()

DBcaller = PSzeepAPI.DataModel(APIconfig['ResnetURL'], APIconfig['PSuserName'], APIconfig['PSpassword'])
networx = ZeeptoNX.PSNetworx(DBcaller)

def ExecutionTime(execution_start):
    import time
    from datetime import timedelta
    return "{}".format(str(timedelta(seconds=time.time()-execution_start)))


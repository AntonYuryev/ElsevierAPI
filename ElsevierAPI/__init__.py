from ElsevierAPI.ResnetAPI.ResnetAPISession import APISession
import json


con_file = open("ElsevierAPI/APIconfig.json")  # file with your API keys and API URLs
APIconfig = json.load(con_file)
con_file.close()

# import ps_api session into your code to start retreival and processing. 
# To create your own API session copy constructor below into your code
ps_api = APISession(APIconfig['ResnetURL'], APIconfig['PSuserName'], APIconfig['PSpassword'])

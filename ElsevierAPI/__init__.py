import json
from ElsevierAPI.ResnetAPI.ResnetAPISession import APISession

#config_dir = 'D:/Python/ENTELLECT_API/ElsevierAPI/'
#config_dir = './ElsevierAPI/.misc/'
config_dir = './ElsevierAPI/'
apiconfig = 'APIconfig.json'
default_apiconfig = config_dir + apiconfig

def load_api_config(api_config_file=None):# file with your API keys and API URLs
    if not isinstance(api_config_file,str):
        print('No API config file was specified\nWil use default %s instead'% default_apiconfig)
        api_config_file = default_apiconfig
    try:
        return json.load(open(api_config_file))
    except FileNotFoundError:
        print("Cannot find API config file: %s" % api_config_file)
        if api_config_file != default_apiconfig:
            return json.load(open(default_apiconfig))
        else:
            print('No working API server was specified!!! Goodbye')

def open_api_session(api_config_file=None) -> APISession:
    APIconfig = load_api_config(api_config_file)
    return APISession(APIconfig['ResnetURL'], APIconfig['PSuserName'], APIconfig['PSpassword'])

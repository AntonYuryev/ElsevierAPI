import json
from .ResnetAPI.ResnetAPISession import APISession
from .ResnetAPI.NetworkxObjects import PSObject, PSRelation
from .ResnetAPI.ResnetGraph import ResnetGraph
from datetime import timedelta
import time

config_dir = './ElsevierAPI/'
apiconfig = 'APIconfig.json'
default_apiconfig = config_dir + apiconfig

def load_api_config(api_config_file=''):# file with your API keys and API URLs
    if not api_config_file:
        print('No API config file was specified\nWill use default %s instead'% default_apiconfig)
        api_config_file = default_apiconfig
        
    try:
        return dict(json.load(open(api_config_file)))
    except FileNotFoundError:
        print("Cannot find API config file: %s" % api_config_file)
        if api_config_file != default_apiconfig:
            print('Cannot open %s config file\nWill use default %s instead'% (api_config_file, default_apiconfig))
            return dict(json.load(open(default_apiconfig)))
        else:
            print('No working API server was specified!!! Goodbye')
            return None


def open_api_session(api_config_file='',what2retrieve=1) -> APISession:
    APIconfig = load_api_config(api_config_file)
    return APISession(APIconfig,what2retrieve=what2retrieve)


def execution_time(execution_start):
    return "{}".format(str(timedelta(seconds=time.time() - execution_start)))

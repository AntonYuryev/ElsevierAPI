import urllib.request
import urllib.parse
import json
import urllib.error as http_error
from time import sleep

class Pharmapendium:
    url = 'https://api.elsevier.com/pharma/'
    page_size = 100 #controls number of records downloaded in one get request

    def __init__(self,APIconfig:dict,add_param=dict()):
        self.headers = {'X-ELS-APIKey':APIconfig['ELSapikey'],'X-ELS-Insttoken':APIconfig['insttoken']}
        self.params = dict()
        self.params.update(add_param)

    def _get_param_str(self):
        return urllib.parse.urlencode(self.params,doseq=True)

    def _url_request(self):
        return self.url+self._get_param_str()
    
    
    def __get_results(self):
        try:
            req = urllib.request.Request(url=self._url_request(), headers=self.headers)
            response = urllib.request.urlopen(req)
        except http_error.HTTPError:
            sleep(30)
            try:
                req = urllib.request.Request(url=self._url_request(), headers=self.headers)
                response = urllib.request.urlopen(req)
            except http_error.HTTPError:
                raise http_error.HTTPError

        pp_view = response.read()
        result = json.loads(pp_view.decode('utf-8'))
        sleep(0.34)
        return result


    def _add_param(self,to_add:dict()):
        self.params.update(to_add)


    def get_results(self):
        result = self.__get_results()
        result_count = result['data']['countTotal']

        all_items = list()
        for page in range(0,result_count,self.page_size):
            self._add_param({'limitation.firstRow':page+self.page_size})
            result = self.__get_results()
            items = result['data']['items']
            all_items += items

        return all_items


class SafetyPP(Pharmapendium):
    def __init__(self,search_type:str, APIconfig:dict):
        self.url = super().url+'safety/'+search_type+'?'
        super().__init__(APIconfig)

    def GetTopEffectCategory(self,taxName):
        add_params = {'taxonomy': 'Effects', 'query':taxName}
        self._add_param(add_params)
        result = self.__get_results()
        if result > 0:
            return result['children'][0]['data']['name']
        else: 
            return ''


class DrugIndications(Pharmapendium):
    def __init__(self,search_type:str, APIconfig:dict):
        self.url = super().url+'drugsindications/'+search_type+'?'
        super().__init__(APIconfig)

    def drugs4indication(self,indications:list):
        self._add_param({'indications':indications})
        result = self.__get_results()
        if result:
            return result['data']['items']
        else: 
            return ''


class DrugActivity(Pharmapendium):
    def __init__(self,search_type:str, APIconfig:dict):
        self.url = super().url+'activity/'+search_type+'?'
        super().__init__(APIconfig)

    def targets4(self,drugs:list):
        self._add_param({'drugs':drugs})
        result = self.__get_results()
        if result:
            return result['data']['items']
        else: 
            return ''
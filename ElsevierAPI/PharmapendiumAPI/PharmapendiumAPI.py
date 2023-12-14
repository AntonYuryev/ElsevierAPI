import urllib.request
import urllib.parse
import json
import urllib.error as http_error
from time import sleep

DEFAULT_APICONFIG = 'D:/Python/ENTELLECT_API/ElsevierAPI/APIconfig.json'
ALPHABET = [chr(i) for i in range(97, 123)] + [str(i) for i in range(10)] + ['-']# + ['(','-','+',' ']

def load_api_config(api_config_file=''):# file with your API keys and API URLs
    if not api_config_file:
        print('No API config file was specified\nWill use default %s instead'% DEFAULT_APICONFIG)
        api_config_file = DEFAULT_APICONFIG
    try:
        return dict(json.load(open(api_config_file,'r')))
    except FileNotFoundError:
        print("Cannot find API config file: %s" % api_config_file)
        if api_config_file != DEFAULT_APICONFIG:
            print('Cannot open %s config file\nWill use default %s instead'% (api_config_file, DEFAULT_APICONFIG))
            return dict(json.load(open(DEFAULT_APICONFIG,'r')))
        else:
            print('No working API server was specified!!! Goodbye')
            return dict()
        

class Pharmapendium:
    url = 'https://api.elsevier.com/pharma/'
    page_size = 100 #controls number of records downloaded in one get request
    cache_dir = 'ElsevierAPI/ResnetAPI/__ppcache__/'

    def __init__(self,APIconfig='',add_param=dict()):
        my_apiconfig = load_api_config(APIconfig)
        self.headers = {'X-ELS-APIKey':my_apiconfig['ELSapikey'],'X-ELS-Insttoken':my_apiconfig['insttoken']}
        self.params = dict()
        self.params.update(add_param)
        self.api_source = ''


    def _add_param(self,to_add:dict={}):
        self.params.update(to_add)


    def _get_param_str(self):
        return urllib.parse.urlencode(self.params,doseq=True)
    

    def _url_request(self,search_type:str, params=dict()):
        my_url = self.url +search_type
        if params:
            my_url += '?'+self._get_param_str()
        else:
            my_url += '/'
        return my_url
    
    
    def _get_results(self, search_type:str, params=dict()):
        my_url = self._url_request(search_type, params)
        try:
            req = urllib.request.Request(my_url, headers=self.headers)
            response = urllib.request.urlopen(req)
        except http_error.HTTPError:
            sleep(30)
            try:
                req = urllib.request.Request(my_url, headers=self.headers)
                response = urllib.request.urlopen(req)
            except http_error.HTTPError as error:
                raise error
            
        pp_view = response.read()
        result = json.loads(pp_view.decode('utf-8'))
        sleep(0.34)
        return result


    def search_results(self,search_type:str, params=dict()):
        result = self._get_results(search_type,params)
        result_count = result['data']['countTotal']

        all_items = list()
        for page in range(0,result_count,self.page_size):
            self._add_param({'limitation.firstRow':page})
            result = self._get_results(search_type,params)
            items = list(result['data']['items'])
            all_items += items

        return all_items
    

    def _taxonomies(self):
        return self._get_results('listTaxonomies')
    

    def _listDataFields(self):
        return self._get_results('listDataFields')
    

    def _listFacets(self):
        return self._get_results('listFacets')
    

    def fetch_all(self,prefix:str):
        self._add_param({'prefix':prefix})
        if len(prefix) > 4: return [] # break for common words like 'Acid', 'acetate'
        childs = self._get_results('suggest',self.params)
        child_collector = list()
        if len(childs) >= 20:
            for l in ALPHABET:
                childs = self.fetch_all(prefix+l)
                child_collector += childs
        else:
            return childs
        return child_collector


    def tax_childs(self,taxonomy:str):
        cache_name = f'{taxonomy} taxonomy 4 PP {self.api_source} module.json'
        cache_path = self.cache_dir + cache_name
        try:
            with open(cache_path,'r',encoding='utf=8') as f:
                return json.load(f)
        except FileNotFoundError:
            children = list()
            self._add_param({'taxonomy':taxonomy})
            chars = ALPHABET + ['(','-','+']
            for l in chars:
                children += self.fetch_all(l)
            
            children = list(set(children))
            with open(cache_path,'w',encoding='utf-8') as f:
                json.dump(children,f,indent=2)

            return children


class SafetyPP(Pharmapendium):
    def __init__(self,APIconfig=''):
        super().__init__(APIconfig)
        self.api_source = 'safety'
        self.url = super().url+self.api_source+'/'
        

    def GetTopEffectCategory(self,taxName):
        add_params = {'taxonomy': 'Effects', 'query':taxName}
        self._add_param(add_params)
        result = self._get_results('lookupFuzzy',self.params)
        if result:
            return result['children'][0]['data']['name']
        else: 
            return ''


class DrugIndications(Pharmapendium):
    def __init__(self,APIconfig=''):
        super().__init__(APIconfig)
        self.api_source = 'drugsindications'
        self.url = super().url+self.api_source+'/'
        

    def drugs4indication(self,indications:list):
        self._add_param({'indications':indications})
        return self.search_results('search',self.params)


class DrugActivity(Pharmapendium):
    def __init__(self,APIconfig=''):
        super().__init__(APIconfig)
        self.api_source = 'activity'
        self.url = super().url+self.api_source+'/'
        self.page_size = 500
        

    def targets4(self,drugs:list):
        self._add_param({'drugs':drugs})
        return self.search_results('search',self.params)


    def drugs(self):
        drugs = list()
        cache_name = f'Drugs 4 PP {self.api_source} module.json'
        cache_path = self.cache_dir + cache_name
        try:
            with open(cache_path,'r',encoding='utf=8') as f:
                drugs =  json.load(f)
                print(f'Read {len(drugs)} drugs from "{cache_name}"')
                return drugs
        except FileNotFoundError:
            drug_names = self.tax_childs('Drugs')
            step_size = 100
            for step in range(0,len(drug_names),step_size):
                self._add_param({'drugs':','.join(drug_names[step:step+step_size])})
                result = self.search_results('search',self.params)
                drugs += result

            with open(cache_path,'w',encoding='utf=8') as f:
                json.dump(drugs,f)
            
            print(f'Dumped {len(drugs)} drugs to "{cache_name}"')
            return drugs


class PPDoc(Pharmapendium):
    def __init__(self,search_type:str, APIconfig=''):
        self.url = super().url+'documents/'+search_type+'?'
        super().__init__(APIconfig)

    
    def get_doc(self,citation_id:str):
        self._add_param({'Id':citation_id})
        return self.search_results('seach')



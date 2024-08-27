import urllib.request,urllib.parse,json,time
from collections import defaultdict
from ..pandas.panda_tricks import df
from ..ResnetAPI.NetworkxObjects import Reference,PUBYEAR
from .. import execution_time
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
    cache_dir = 'ElsevierAPI/PharmapendiumAPI/__ppcache__/'

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


    def retreive_taxonomy(self,taxonomy:str):
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
        
    
    def drugs_by_name(self,drugs:list[str]):
        self._add_param({'drugs':','.join(drugs)})
        return self.search_results('search',self.params)


    def drugs(self)->list[dict]:
        cache_name = f'Drugs 4 PP {self.api_source} module.json'
        cache_path = self.cache_dir + cache_name
        try:
            with open(cache_path,'r',encoding='utf=8') as f:
                drug_records =  list(json.load(f))
                print(f'Read {len(drug_records)} drugs from "{cache_name}"')
                return drug_records
        except FileNotFoundError:
            drug_records = list()
            drug_names = self.retreive_taxonomy('Drugs')
            step_size = 100
            for step in range(0,len(drug_names),step_size):
                result = self.drugs_by_name(drug_names[step:step+step_size])
                drug_records += result
                print(f'Downloaded {len(result)} records for {step_size} drug names')

            drug_records.sort()
            with open(cache_path,'w',encoding='utf=8') as f:
                json.dump(drug_records,f)
            
            print(f'Dumped {len(drug_records)} drugs to "{cache_name}"')
            return drug_records


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
        
    
    def ToxReport(self,drug_names:list[str]):
        '''
        Return
        ------
        df with columns: 'Smiles','Drug','Dose','Dose type','Route','Tox Category','#Ref','References'
        '''
        SMILES = 'Smiles'
        DRUG = 'Drug'
        DOSE = 'Dose'
        DOSETYPE = 'Dose type'
        ROUTE = 'Route'
        TOXICITY = 'Toxicity'
        TOXCATEGORY = 'Tox Category'
        REFCOUNT = '#Ref'
        REFERENCES = 'References'

        col_names = [SMILES,DRUG,DOSE,DOSETYPE,ROUTE,TOXICITY,TOXCATEGORY,REFCOUNT,REFERENCES]
        tox_df = df(columns=col_names)
        tox_df.index.name = "PP name\ttoxicity"
        self._add_param({'taxonomy':'Drugs'})

        start_time = time.time()
        print (f'Will find drug safety documents for {len(drug_names)} drugs')
        for i,drugPPname in enumerate(drug_names):
            self._add_param({'drugs':drugPPname})
            drug_safety_docs = self.search_results('search',self.params)
           # if drugPPname == 'Vorinostat':
           #     print('')
        
            DrugToxicities = defaultdict(list)
            [DrugToxicities[doc['effect']].append(doc) for doc in drug_safety_docs]

            toxCount = 0
            for toxicity, references in DrugToxicities.items():
                toxCount += 1
                print(f'\'{toxicity}\' - {toxCount} from {len(DrugToxicities)} toxicities for \
\"{drugPPname}\" ({i+1} of {len(drug_names)}) was reported in {len(references)} documents')
                
                pandaIndex = drugPPname+'\t'+toxicity
                tox_df.at[pandaIndex,SMILES] = references[0].get('smiles','')
                tox_df.at[pandaIndex,DRUG] = drugPPname
                tox_df.at[pandaIndex,DOSE] = references[0].get('dose','')
                tox_df.at[pandaIndex,DOSETYPE] = references[0].get('doseType','')
                tox_df.at[pandaIndex,ROUTE] = references[0].get('route','')
                tox_df.at[pandaIndex,TOXICITY] = toxicity
                tox_df.at[pandaIndex,TOXCATEGORY] = self.GetTopEffectCategory(toxicity)

                refIndex = dict() # {ref_identifiere:Reference}
                for ref in references:
                    assert(isinstance(ref,dict))
                    document = ref['document']
                    assert(isinstance(document,dict))

                    try:docName = document['name']
                    except KeyError:
                        try: docName = document['article']
                        except KeyError: docName = document['journal']
                    
                    docSource = document['sourceShort']
                    refIdentifier = docSource+':'+docName
                    try: PPRef = refIndex[refIdentifier]
                    except KeyError:
                        PPRef = Reference('Title',refIdentifier)
                        refIndex[refIdentifier] = PPRef

                    assert(isinstance(PPRef,Reference))

                    PubYear = str(document.get('year','historic'))
                    PPRef.update_with_value(PUBYEAR, PubYear)
                    
                    dose = ref.get('dose','')
                    doseType = ref.get('doseType','')
                    route = ref.get('route','')
                    organism=ref['specie']
                    PPRef.update_with_value(route, doseType + ' in ' + organism + ' ' + dose)

                addToPandas = set()
                drug_refs = list(refIndex.values())
                [addToPandas.update([r.to_str(['Title'], col_sep='-',biblio_props=[PUBYEAR],other_props=[route])]) for r in drug_refs]

                tox_df.at[pandaIndex,REFCOUNT] = len(addToPandas)
                reflist = '|'.join(list(addToPandas))
                tox_df.at[pandaIndex,REFERENCES] = reflist

        toxicity_counts = dict(tox_df[TOXICITY].value_counts())
        tox_df = tox_df.merge_dict(toxicity_counts,'Toxicity count','Toxicity')
        tox_df = df.from_pd(tox_df.sort_values(by=['Toxicity count','Toxicity'],ascending=False))
        print(f'Finished finding toxicities in Pharmapendium in {execution_time(start_time)}')
        return tox_df


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


    def drugs4targets(self,target_names:list[str],only_primary=True):
        all_drugs = self.drugs()
        target2drug = defaultdict(list)
        if only_primary:
            [target2drug[d['target']].append(d) for d in all_drugs if d.get('isPrimaryTarget','') == 'Primary']
        else:
            [target2drug[d['target']].append(d) for d in all_drugs] 
        drug_records4targets = {t:drugs for t,drugs in target2drug.items() if t in target_names}
        return drug_records4targets

        

class PPDoc(Pharmapendium):
    def __init__(self,search_type:str, APIconfig=''):
        self.url = super().url+'documents/'+search_type+'?'
        super().__init__(APIconfig)

    
    def get_doc(self,citation_id:str):
        self._add_param({'Id':citation_id})
        return self.search_results('seach')


class FAERS(Pharmapendium):
    def __init__(self,APIconfig=''):
        super().__init__(APIconfig)
        self.api_source = 'faers'
        self.url = super().url+self.api_source+'/'
        self.page_size = 500


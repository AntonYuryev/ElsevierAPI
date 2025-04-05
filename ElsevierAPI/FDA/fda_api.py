import urllib.request, urllib.parse, requests,os,json
import regex as re
import urllib.error as http_error
from time import sleep
from concurrent.futures import ThreadPoolExecutor,as_completed
from ..utils import list2chunks

DRUGNAME_FIELDS = ['openfda.brand_name', 'openfda.generic_name','openfda.substance_name']
DOSAGE_FIELD = 'dosage_and_administration'
CACHE_DIR = os.path.join(os.getcwd(),'ENTELLECT_API/ElsevierAPI/FDA/drug_labels_cache')

def get_field(label:dict,keys:list):
    try:
        return get_field(label[keys[0]], keys[1:]) if keys else label[0]
    except KeyError:
        return ''

def replace_in(s:str,any_pattern_in:list,with_str=''):
    new_str = s
    for pattern in any_pattern_in:
        new_str = re.sub(pattern, lambda x: with_str, new_str)

    return new_str

class FDA:
    page_size = 100 # Controls number of records downloaded in one get request.
    api_documentation_url = 'https://open.fda.gov/apis/drug/drugsfda/explore-the-api-with-an-interactive-chart'
    DATA_DIR = 'D:/Python/Drug labels/'


    def __init__(self):
        self.url = 'https://api.fda.gov/drug/label.json?search='
        self.drug2labels = dict()
        

    def __get_results(self, field:str, query:str,limit=100):
      retries = 3
      for attempt in range(retries):
        try:
          quoted = urllib.parse.quote(query)
          url_request = self.url+field+':'+quoted+'&limit='+str(limit)
          return requests.get(url_request)
        except (http_error.HTTPError,requests.exceptions.RequestException ) as e:
          print(f'An error occurred: {e}')
          if attempt < retries - 1:
            print(f'Pausing for 30 seconds due to HTTPError, attempt {attempt + 1} of {retries}')
            sleep(30)
        return dict()
      
    @staticmethod
    def __cache_name(drug_name:str):
      return drug_name.lower()


    def __load_drug_label_cache(self,drugs:list[str]):
      drug2labels = dict()
      for drug in drugs:
        drug_name_in_cache = self.__cache_name(drug)
        cache_file = os.path.join(CACHE_DIR,drug_name_in_cache+'.json')
        if os.path.exists(cache_file):
          with open(cache_file,'r') as f:
            drug2labels[drug_name_in_cache] = json.load(f)
      self.drug2labels.update(drug2labels)
      return drug2labels
    

    def results(self,field:str, query:str):
      response = self.__get_results(field,query)
      if not response:
        return []
      
      try:
        res = response.json()
      except ValueError:
        print('Error: Response is not valid JSON')
        return []
      
      all_results = []
      if 'error' not in res:
        all_results = res['results']
        while response.links.get('next'):
          url_request = response.links['next']['url']
          response = requests.get(url_request)
          try:
            res_ = response.json()
          except ValueError:
            print('Error: Response is not valid JSON')
            break
          if 'error' not in res_:
            all_results += res_['results']
      return all_results
    

    def drug_labels(self,drug:str)->list[dict]:
      def find_label(drug:str):
        for field in DRUGNAME_FIELDS:
          labels = self.results(field,drug)
          if labels:
            return labels
        return []
      
      drug_name = self.__cache_name(drug)
      if drug_name in self.drug2labels:
        return self.drug2labels[drug_name]
      else:
        if ' ' in drug:
          quoted_drug =  f'"{drug}"'
          labels = find_label(quoted_drug)
        else:
          labels = find_label(drug)
          if not labels:
            labels = find_label(drug+'*')

        self.drug2labels[drug_name] = labels
        cache_file_path = os.path.join(CACHE_DIR, drug_name + '.json')
        os.makedirs(os.path.dirname(cache_file_path), exist_ok=True)
        with open(cache_file_path, 'w') as f:
          json.dump(labels, f,indent=2)
        return labels


    child_dose_pattern = re.compile(r"\b(child|pediat)\b", re.IGNORECASE)
    def __child_doses(self, drug:str):
      drug2dose = dict()
      labels = self.drug_labels(drug)
      for label in labels:
        dosage = str(label.get(DOSAGE_FIELD,''))
        sentences = dosage.split('. ')
        for sentence in sentences:
          if self.child_dose_pattern.search(sentence):
            brand = get_field(label,['openfda','brand_name'])
            sentence = sentence.strip(' []\'".')
            drug2dose[brand]= sentence
            break
      return drug2dose


    def child_doses(self,drugs:list[str])->dict[str,str]:
      self.__load_drug_label_cache(drugs)
      drug2childdose = dict()
      for drug in drugs:
        brands2dose = self.__child_doses(drug)
        if brands2dose:
          drug2childdose[drug] = '. '.join(brands2dose.values())
      return drug2childdose
    

    def child_doses_mt(self,drugs:list[str],max_workers=5):
      drug2childdose = dict()
      print(f'Finding children dosage for {len(drugs)} drugs in {max_workers} threads')
      chunks = list2chunks(drugs, max_workers)
      with ThreadPoolExecutor(max_workers=max_workers, thread_name_prefix='Find child doses') as e:
        futures = [e.submit(self.child_doses, chunk) for chunk in chunks]
        [drug2childdose.update(f.result()) for f in as_completed(futures)]
      
      return drug2childdose
    

    @staticmethod
    def parse_drug_label(label:dict):
        try:
            purposes = label["purpose"]
            purposes = ';'.join(purposes).lower()
            purposes = replace_in(purposes,['uses?:? ','purposes?:? '], '')
        except KeyError:
            purposes = str()
        
        brand = get_field(label,['openfda','brand_name'])
        generic = get_field(label,['openfda','generic_name'])
        substance = get_field(label,['openfda','substance_name'])
        ndc = get_field(label,['openfda','product_ndc'])
        moa = get_field(label,['openfda','pharm_class_moa'])
        if not moa:
            moa = get_field(label,['mechanism_of_action'])
            if moa:
                moa += '[Clin.Pharm]'
        act_ingr = get_field(label,['active_ingredient'])
        return [ndc,brand,generic,substance,act_ingr,moa,purposes]


    def drugs4indication(self,to_file:str, for_indication='pain'):
      labels = self.results('purpose',for_indication)
      with open(to_file, 'w', encoding='utf-8') as f:
        f.write('NDC\tBrand name\tGeneric name\tSubstance\tAct.Ingridient\tMOA\tIndication\n')
        for label in labels:
          purposes = label["purpose"]
          for purpose in purposes:
            if str(purpose).lower().find(for_indication) >= 0:
              row = self.parse_drug_label(label)
              f.write('\t'.join(row)+'\n')


    def write_drugs2table(self,drug_labels:list,to_file:str):
      with open(to_file, 'w', encoding='utf-8') as f:
        f.write('NDC\tBrand name\tGeneric name\tSubstance\tAct.Ingridient\tMOA\tIndication\n')
        for label in drug_labels:
          row = self.parse_drug_label(label)
          f.write('\t'.join(row)+'\n')



if __name__ == "__main__":
    fda_api = FDA()
    #generics = fda_api.results('_exists_','openfda.generic_name')
    #fda_api.write_drugs2table(generics,fda_api.DATA_DIR+'Generics.tsv')

    drugs = ['Melatonin','Cocaine Hydrochloride']
    drug2childdose = dict()
    for drug in drugs:
      brands2dose = fda_api.__child_doses(drug)
      drug2childdose[drug] = brands2dose
      continue

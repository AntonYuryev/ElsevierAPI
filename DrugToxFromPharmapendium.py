import urllib.request
import urllib.parse
import json
import time
import ElsevierAPI.ResnetAPI.PathwayStudioGOQL as OQL
from ElsevierAPI.ResnetAPI.NetworkxObjects import Reference
from ElsevierAPI import load_api_config
from ElsevierAPI.ResnetAPI.ResnetAPISession import APISession
import pandas as pd

APIconfig = load_api_config()
ps_api = APISession(APIconfig['ResnetURL'], APIconfig['PSuserName'], APIconfig['PSpassword'])
fileIn = 'Drugs for Regulators in 4 patients.txt'
InDir = 'D:\\Python\\PBTA\\PNOC003\\4 patients analysis\\'
with open(InDir+fileIn) as f:
    drugs = [line.rstrip('\n') for line in f]

print('Finding drugs in %s in Resnet' %(fileIn))
OQLquery = OQL.get_entities_by_props(drugs, ['Name', 'Alias'], only_object_types=['Small Molecule'])
ps_api.add_ent_props(['Name','PharmaPendium ID'])
resnet_drugs = ps_api.process_oql(OQLquery,'Find all drugs')
print ('Found %d drugs in Resnet' % len(resnet_drugs))

#removing duplicates wiht no PharmaPendium ID
resnet2pharmapendium_map = dict()
for i, drug in resnet_drugs.nodes(data=True):
    try:
        resnet2pharmapendium_map[str(drug['Name'][0]).lower()] = drug['PharmaPendium ID'][0]
    except KeyError: continue

all_drugs = list(resnet_drugs.nodes(data=True))
for i, drug in all_drugs:
    if drug['Name'][0].lower() in resnet2pharmapendium_map.keys() and 'PharmaPendium ID' not in drug.keys():
        resnet_drugs.remove_node(i)
print ('%d drugs left after deduplication' % resnet_drugs.number_of_nodes())


print('Beginning Pharmapendium response download with urllib...')
ELSapiKey = APIconfig['ELSapikey']#Obtain from https://dev.elsevier.com
token = APIconfig['insttoken'] #Obtain from mailto:integrationsupport@elsevier.com 

PP_URL = 'https://api.elsevier.com/pharma/'
PPmodule ='safety/'
RequestType = 'search'

taxonomy = 'Effects'
headers = {'X-ELS-APIKey':ELSapiKey,'X-ELS-Insttoken':token}
baseURL =  PP_URL+PPmodule+RequestType +'?'
PageLimit = 500 #controls number of records downloaded in one request. Cannot exceed 500

LazyDictTox = dict()
def GetTopTaxCategory(taxName):
    try: return LazyDictTox[taxName]
    except KeyError:
        params = {'taxonomy': taxonomy, 'query':taxName}
        paramtURL = urllib.parse.urlencode(params)
        baseURL2 = PP_URL+PPmodule+'lookupFuzzy?'
        req = urllib.request.Request(url=baseURL2+paramtURL,headers=headers)
        with urllib.request.urlopen(req) as response:
            the_page = response.read()
            result = json.loads(the_page.decode('utf-8'))
            if len(result) > 0:
                topCategory = result['children'][0]['data']['name']
                LazyDictTox[taxName]= topCategory
                return topCategory
            else:
                LazyDictTox[taxName] = ''
                return ''


fileOut2 = InDir + fileIn[:len(fileIn)-4]+'_unmapped.txt'
col_names = ['Smiles','Resnet name','Tox Category','#Ref','References']
ToxPandas = pd.DataFrame(columns=col_names)
ToxPandas.index.name = "PP name\ttoxicity"

start_time = time.time()
print ('Will find toxicities for %d drugs found in Resnet' % resnet_drugs.number_of_nodes())
for i, drug in resnet_drugs.nodes(data=True):
    drugPSname = drug['Name'][0]
    try: drugPPname = drug['PharmaPendium ID'][0]
    except KeyError: drugPPname = drugPSname

    params = {'drugs':drugPPname}
    paramtURL = urllib.parse.urlencode(params)
    req = urllib.request.Request(url=baseURL+paramtURL,headers=headers)
    response = urllib.request.urlopen(req)
    the_page = response.read()
    result = json.loads(the_page.decode('utf-8'))
    docCount = result['data']['countTotal']

    if docCount == 0: 
        print('cannot find %s in Pharmapendium' % drugPSname)
        continue

    print('Found %d documents for %s with Pharmapendium ID %s' % (docCount,drugPSname, drugPPname))
    DrugToxicities = dict()

    for page in range(1,docCount,PageLimit):
        for item in result['data']['items']:
            toxicity=item['effect']
            try: DrugToxicities[toxicity].append(item)
            except KeyError:
                DrugToxicities[toxicity] = [item]

        if page+PageLimit < docCount:
            params = {'drugs':drugPPname,'limitation.firstRow':page+PageLimit}
            paramtURL = urllib.parse.urlencode(params)
            req = urllib.request.Request(url=baseURL+paramtURL,headers=headers)
            response = urllib.request.urlopen(req)
            the_page = response.read()
            result = json.loads(the_page.decode('utf-8'))
        else:
            toxCount = 0
            for toxicity, references in DrugToxicities.items():
                toxCount += 1
                print('\'%s\' - %d from %d toxicities for \"%s\" was reported in %d documents' % (toxicity,toxCount,len(DrugToxicities),drugPSname,len(references)))
                pandaIndex = drugPPname+'\t'+toxicity
                
                try: smiles = references[0]['smiles']
                except KeyError: smiles=''

                toxTax = GetTopTaxCategory(toxicity)

                ToxPandas.at[pandaIndex,col_names[0]] = smiles
                ToxPandas.at[pandaIndex,col_names[1]] = drugPSname
                ToxPandas.at[pandaIndex,col_names[2]] = toxTax

                refIndex = dict()
                for ref in references:
                    document = ref['document']

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

                    try:PubYear = str(document['year'])
                    except KeyError: PubYear = 'historic'
                    PPRef.set_property('PubYear', PubYear)
                    
                    try:dose = ref['dose']
                    except KeyError: dose = ''

                    try: doseType=ref['doseType']
                    except KeyError: doseType=''
                    
                    try:route=ref['route']
                    except KeyError: route=''
                    
                    organism=ref['specie']
                    PPRef.update_with_value(route, doseType + ' in ' + organism + ' ' + dose)

                addToPandas = set()
                for ref in refIndex.values():
                    addToPandas.update([ref.to_str('Title', sep='-')])

                ToxPandas.at[pandaIndex,col_names[3]] = len(addToPandas)
                reflist = '|'.join(list(addToPandas))
                ToxPandas.at[pandaIndex,col_names[4]] = reflist  
            break
        
ToxPandas.to_csv(InDir+fileIn[:len(fileIn)-4]+'_PPtaxonomy.txt',sep='\t')
print('Finished finding toxicities in Pharmapendium in %s' % ps_api.execution_time(start_time))

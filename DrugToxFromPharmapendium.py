import urllib.request
import urllib.parse
import json
import time
import ElsevierAPI
import ElsevierAPI.ResnetAPI.PathwayStudioGOQL as OQL
import ElsevierAPI.ResnetAPI.ZeepToNetworkx as znx
from ElsevierAPI import networx as PSnx
import pandas as pd

def OpenFile(fname):
    open(fname, "w", encoding='utf-8').close()
    return open(fname, "a", encoding='utf-8')


print('Beginning Pharmapendium response download with urllib...')
con_file = open("ElsevierAPI/config.json")#file with your APIkeys
config = json.load(con_file)
con_file.close()


fileIn = 'Drugs for Regulators in 4 patients.txt'
InDir = 'D:\\Python\\PBTA\\PNOC003\\'
with open(InDir+fileIn) as f:
    drugs = [line.rstrip('\n') for line in f]


OQLquery = OQL.GetEntitiesByProps(drugs, ['Name','Alias'], OnlyObjectTypes=['Small Molecule'])
PSdrugs = PSnx.GetPSObjects(OQLquery, RetreiveProperties=['Name','PharmaPendium ID'])
print ('Found %d drugs in Resnet' % len(PSdrugs))
#removing duplicates
PStoPPNames = dict()
for drug in PSdrugs:
    try:
        PStoPPNames[(drug['Name'][0]).lower()] = drug['PharmaPendium ID'][0]
    except KeyError: continue

for drug in PSdrugs:
    if drug['Name'][0].lower() in PStoPPNames.keys() and 'PharmaPendium ID' not in drug.keys():
        PSdrugs.remove(drug)
print ('%d drugs left after deduplication' % len(PSdrugs))

ELSapiKey = config['ELSapikey']#Obtain from https://dev.elsevier.com
token = config['insttoken'] #Obtain from mailto:integrationsupport@elsevier.com 

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

PPtoPSName = dict()
start_time = time.time()
print ('Will find toxicities for %d drugs found in Resnet' % len(PSdrugs))
for drug in PSdrugs:
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

    print('Will parse %d documents for %s with Pharmapendium ID %s' % (docCount,drugPSname, drugPPname))
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
                        PPRef = znx.Reference('Title',refIdentifier)
                        refIndex[refIdentifier] = PPRef

                    try:PubYear = str(document['year'])
                    except KeyError: PubYear = 'historic'
                    PPRef.AddSingleProperty('PubYear',PubYear)
                    
                    try:dose = ref['dose']
                    except KeyError: dose = ''

                    try: doseType=ref['doseType']
                    except KeyError: doseType=''
                    
                    try:route=ref['route']
                    except KeyError: route=''
                    
                    organism=ref['specie']
                    PPRef.AddUniqueProperty(route,doseType+' in '+organism+' '+dose)

                addToPandas = set()
                for ref in refIndex.values():
                    addToPandas.update([ref.ToString('Title', sep='-')])

                ToxPandas.at[pandaIndex,col_names[3]] = len(addToPandas)
                reflist = '|'.join(list(addToPandas))
                ToxPandas.at[pandaIndex,col_names[4]] = reflist
                
            break
        
ToxPandas.to_csv(InDir+fileIn[:len(fileIn)-4]+'_PPtaxonomy.txt',sep='\t')
print('Finished finding toxicities in Pharmapendium in %s' % ElsevierAPI.ExecutionTime(start_time))


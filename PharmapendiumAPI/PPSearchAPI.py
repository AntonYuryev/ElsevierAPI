import urllib.request
import urllib.parse
import json
from ElsevierAPI import APIconfig

print('Beginning Pharmapendium response download with urllib...')
ELSapiKey = APIconfig['ELSapikey']#Obtain from https://dev.elsevier.com
token = APIconfig['insttoken'] #Obtain from mailto:integrationsupport@elsevier.com 

PP_URL = 'https://api.elsevier.com/pharma/'
PPmodule ='safety/'
SearchType = 'lookupFuzzy'
taxonomy = 'Effects'
PageSize = 100 #controls number of records downloaded in one get request
headers = {'X-ELS-APIKey':ELSapiKey,'X-ELS-Insttoken':token}
baseURL =  PP_URL+PPmodule+SearchType +'?'

def GetTopTaxCategory(taxName):
    params = {'taxonomy': taxonomy, 'query':taxName}
    paramtURL = urllib.parse.urlencode(params)
    req = urllib.request.Request(url=baseURL+paramtURL, headers=headers)
    with urllib.request.urlopen(req) as response:
        the_page = response.read()
        result = json.loads(the_page.decode('utf-8'))
        if len(result) > 0:
            return result['children'][0]['data']['name']
        else: 
            return ''


fileIn = 'PBTA_PPtoxicities.txt'
InDir = 'D:\\Python\\PBTA\\'
with open(InDir+fileIn) as f:
    TaxNames = [line.rstrip('\n') for line in f]

fileOut = InDir + fileIn[:len(fileIn)-4]+'_PPtaxonomy.txt'
fileOut2 = InDir + fileIn[:len(fileIn)-4]+'_unmapped.txt'
f2 = open(fileOut2, 'w', encoding='utf-8') 

with open(fileOut, 'w', encoding='utf-8') as f:
    for taxName in TaxNames:
        TopTaxCategory = GetTopTaxCategory(taxName)
        new_taxName = ''
        if len(TopTaxCategory) == 0:
            if taxName[:len(taxName)-4] == 'emia':
                new_taxName = taxName[:len(taxName)-4] + 'aemia'
                TopTaxCategory = GetTopTaxCategory(new_taxName)
        if len(TopTaxCategory) == 0:
            if taxName[len(taxName)-1] == 's':
                new_taxName = taxName[:len(taxName)-1]
                TopTaxCategory = GetTopTaxCategory(new_taxName)
        if len(TopTaxCategory) == 0:
            if 'edema' in taxName:
                new_taxName = taxName.replace('edema', 'oedema')
                TopTaxCategory = GetTopTaxCategory(new_taxName)
        if len(TopTaxCategory) == 0:
            if 'ea' in taxName:
                new_taxName = taxName.replace('ea', 'oea')
                TopTaxCategory = GetTopTaxCategory(new_taxName)

        if len(TopTaxCategory) > 0:
            f.write(taxName+'\t'+ TopTaxCategory+'\t'+new_taxName+'\n')
        else:
            f2.write(taxName+'\n')
            print('Cannot find top taxonomy category for %s' % taxName)

f2.close()
print('Finished taxonomy mapping!')

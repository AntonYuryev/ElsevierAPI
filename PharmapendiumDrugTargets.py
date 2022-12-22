from ElsevierAPI import load_api_config
from ElsevierAPI.PharmapendiumAPI.PharmapendiumAPI import DrugActivity,DrugIndications
from ElsevierAPI.ETM_API.medscan import MedScan

APIconfig = load_api_config()
pp_api_di = DrugIndications('search', APIconfig)
pp_api_da = DrugActivity('search', APIconfig)
medscan = MedScan(APIconfig['MedscanLicense'])

def get_drugs(for_indication:str, with_aliases:list):
    drugs = pp_api_di.drugs4indication([for_indication])
    if drugs: return drugs
    for alias in with_aliases:
        drugs = pp_api_di.drugs4indication([alias])
        if drugs: return drugs

    return dict()


def get_rows(drug):
    items = pp_api_da.targets4(drug)
    rows = set()
    for item in items:
        target = item['target']
        if target.lower() in ['unreported'] : continue
        medscan_concept = medscan.find_concepts(target)
        try:
            id2name = dict(list(medscan_concept.values())[0][0])
            medscan_id,medscan_objname = list(id2name.items())[0]
        except KeyError:
            try:
                id2name = dict(list(medscan_concept.values())[0][3000000])
                medscan_id,medscan_objname = list(id2name.items())[0]
            except KeyError:
                try:
                    id2name = dict(list(medscan_concept.values())[0][12000000])
                    medscan_id,medscan_objname = list(id2name.items())[0]
                except KeyError:
                    medscan_id = ''
                    medscan_objname = ''

        primary = item['isPrimaryTarget']
        is_antagonist = item['agonistAntagonist']
        source = item['source']
        tup = tuple([drug,is_antagonist,target,medscan_objname,medscan_id,primary,source])
        rows.add(tup)
    return rows
        

if __name__ == "__main__":
    InDir = 'D:/Python/PMI/'
    drugs = open(InDir+'Analgesics.txt','r').readlines()
    fileOut = InDir + 'Targets4analgesics.txt'

    with open(fileOut, 'w', encoding='utf-8') as f:
        f.write('Drug\tActAs\tPP Target\tMedscan Name\tMedscan ID\tIsPrimary\tSource\n')
        for drug in drugs:
            drug_rows = get_rows(drug.strip())
            for row in drug_rows:
                row_str = '\t'.join(list(row))
                f.write(row_str+'\n')

    print('Finished drug-target retrieval')

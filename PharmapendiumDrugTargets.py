import json
from ElsevierAPI import load_api_config
from ElsevierAPI.PharmapendiumAPI.PharmapendiumAPI import DrugActivity,DrugIndications
from ElsevierAPI.EmbaseAPI.EmbaseSearchAPI import EmbaseAPI,Reference,PUBYEAR,JOURNAL
from ElsevierAPI.ResnetAPI.ResnetAPISession import APISession
from ElsevierAPI.ResnetAPI.ResnetGraph import ResnetGraph,PSRelation,REGULATORS,TARGETS, PSObject,PROTEIN_TYPES
from ElsevierAPI.ETM_API.references import REFERENCE_PROPS,ReferenceEncoder,ReferenceDecoder,SerializableRef
from ElsevierAPI.ETM_API.medscan import MedScan
import urllib.error as http_error
#from requests.exceptions import InvalidURL
from urllib.parse import quote, unquote

def has_valid_url_characters(input_string):
    # Quote the string and then unquote it to check if it changes
    quoted_string = quote(input_string, safe="/:")
    unquoted_string = unquote(quoted_string)
    return input_string == unquoted_string


APIconfig = load_api_config()
pp_api_di = DrugIndications()
pp_api_da = DrugActivity()
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


def make_ref(id_type:str,id:str,pp_drug:dict):
    ref = Reference(id_type,id)
    try:
        ref[PUBYEAR] = [pp_drug['document']['year']]
    except KeyError: pass
    try:
        ref[JOURNAL] = [pp_drug['document']['journal']]
    except KeyError: 
        ref[JOURNAL] = [pp_drug['document']['feature']]
    try:
        ref['Volume'] = [pp_drug['document']['volume']]
    except KeyError: pass
    try:
        ref['Pages'] = [pp_drug['document']['pages']]
    except KeyError: pass


    try:
        ref['TargetType'] = [pp_drug.get('isPrimaryTarget')]
    except KeyError:
        pass
    try:
        ref['Organism'] = [pp_drug['specie']]
    except KeyError: pass

    return ref


def drug_targets():
    drugs = pp_api_da.drugs()
    embase = EmbaseAPI(APIconfig)
    ps_api = APISession()
    drug_target_refs = list()
    drug_names = set()
    target_names = set()
    citationid2ref = dict()
    cache_path = pp_api_da.cache_dir+'pp_drug_target_refs.json'
    try:
        with open(cache_path,'r',encoding='utf-8') as f:
            drug_target_dicts = json.load(f)
        
        drug_target_refs = list()
        for dct in drug_target_dicts:
            my_dict = dct['reference']
            Identifiers = my_dict.pop('Identifiers')
            snippets = my_dict.pop('snippets',dict())
            addresses = my_dict.pop('addresses',dict())

            ref = Reference.from_iddict(Identifiers)
            ref.update(my_dict)
            ref.Identifiers.update(Identifiers)
            ref.snippets.update(snippets)
            ref.addresses.update(addresses)
            drug_target_refs.append(ref)

        drug_names = {ref['Regulator'][0] for ref in drug_target_refs}
        target_names = {ref['Target'][0] for ref in drug_target_refs}
    except FileNotFoundError:
        for i,drug in enumerate(drugs):
            if isinstance(drug,dict):
                drug_name = drug['drug']
                target_name = drug['target']
                citation_id = drug['document']['citation']
                #if has_valid_url_characters(citation_id):
                try:
                    ref = citationid2ref[citation_id]
                    #x = 100
                except KeyError:
                    try:
                        ppdoc = embase.lui2doc(citation_id)
                        ref = EmbaseAPI.article2ref(ppdoc)
                        if isinstance(ref,Reference):
                            try:
                                ref['TargetType'] = [drug.get('isPrimaryTarget')]
                            except KeyError:
                                pass
                            try:
                                ref['Organism'] = [drug['specie']]
                            except KeyError: pass
                    except (http_error.HTTPError):
                        try:
                            ppdoc = embase.pmid2doc(citation_id)
                            ref = EmbaseAPI.article2ref(ppdoc)
                            if isinstance(ref,Reference):
                                try:
                                    ref['TargetType'] = [drug.get('isPrimaryTarget')]
                                except KeyError:
                                    pass
                                try:
                                    ref['Organism'] = [drug['specie']]
                                except KeyError: pass
                        except http_error.HTTPError:
                            try:
                                doi = drug['document']['doi']
                                ref = make_ref('DOI',doi,drug)
                            except KeyError:
                                ref = make_ref('Pharmapendium Document ID',citation_id,drug)

                if isinstance(ref,Reference):
                    ref['Regulator'] = [drug_name]
                    ref['Target'] = [target_name]
                    drug_target_refs.append(SerializableRef(ref))
                    citationid2ref[citation_id] = ref
                    drug_names.add(drug_name)
                    target_names.add(target_name)

        with open(cache_path,'w',encoding='utf-8') as f:
            json.dump(drug_target_refs,f,indent=2,cls=ReferenceEncoder)

    drug_name2psobjs,_ = ps_api.map_props2objs(list(drug_names),["Name","Alias"],['SmallMol'])
    target_name2psobjs,_ = ps_api.map_props2objs(list(target_names),['Name','Alias'],PROTEIN_TYPES)

    PPgraph = ResnetGraph()
    unmapped_drugs = set()
    unmapped_targets = set()
    for i,ref in enumerate(drug_target_refs):
        drug_name = ref['Regulator'][0]
        target_name = ref['Target'][0]
        try:
            regulators = drug_name2psobjs[drug_name]
        except KeyError:
            urn_name = quote(drug_name, safe=':/')
            regulators = [PSObject({'URN':['urn:pp-drug:'+urn_name],'Name':[drug_name]})]
            unmapped_drugs.update(regulators)

        for regulator in regulators:
            regulator['PharmaPendium ID'] = [drug_name]
            try:
                targets = target_name2psobjs[target_name]
            except KeyError:
                urn_name = quote(target_name, safe=':/')
                targets = [PSObject({'URN':['urn:pp-target:'+urn_name],'Name':[target_name]})]
                unmapped_targets.update(targets)

            for target in targets:
                rel = PSRelation({'ObjTypeName':['DirectRegulation']})
                rel.Nodes[REGULATORS] = [(regulator.uid(),0,0)]
                rel.Nodes[TARGETS] = [(target.uid(),1,0)]
                rel._add_refs([ref])
                rel.make_urn([regulator.urn()],[target.urn()])
                PPgraph.add_triple(regulator,target,rel)

    with open(pp_api_da.cache_dir+'unmapped_drugs.json','w',encoding='utf-8') as f:
            json.dump(list(unmapped_drugs),f,indent=2)

    with open(pp_api_da.cache_dir+'unmapped_targets.json','w',encoding='utf-8') as f:
            json.dump(list(unmapped_targets),f,indent=2)

    return PPgraph
        

if __name__ == "__main__":
    drug2target = drug_targets()
    rnef_entprops = ['PharmaPendium ID','Name']
    rnef_relprops= REFERENCE_PROPS + ['Pharmapendium Document ID','TargetType']
    drug2target.dump2rnef('Pharmapendium drug2target.rnef',rnef_entprops,rnef_relprops)
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

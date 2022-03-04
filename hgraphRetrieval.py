import pandas as pd
from ElsevierAPI.ETM_API.references import Reference
from ElsevierAPI.hGraphAPI.wso2processor.hgraphAPIQueryProcessor import HGraphAPIQueryProcessor
from ElsevierAPI.hGraphAPI.cm2access import QPEAccess
from concurrent.futures import ThreadPoolExecutor
import time
from ElsevierAPI import load_api_config
from ElsevierAPI.ResnetAPI.ResnetAPISession import APISession
from ElsevierAPI.ResnetAPI.ResnetGraph import ResnetGraph
from ElsevierAPI.ResnetAPI.NetworkxObjects import PSObject, PSRelation
from ElsevierAPI.ETM_API.etm import ETM
from ElsevierAPI.ETM_API.references import PS_BIBLIO_PROPS,SENTENCE_PROPS,PS_ID_TYPES
import xlsxwriter
from urllib import parse


IMUI_CACHE_CON = {}
IMUI_CACHE_QPE = {}
IMUI_INFO = {}
RELATIONS_OF_INTEREST = ['has diagnostic procedure', 'has clinical finding', 'has risk factor']
FOLDER_NAME = "Data/hGragh/"
RELTYPE2ENTYPE_MAP = {'has diagnostic procedure':'MedicalProcedure', 'has clinical finding':'Disease', 'has risk factor':'RiskFactor'}
DATA_DIR = 'D:/Python/Quest/report tables/'


SEMinPG = '#Practice guidelines'
COOCinPG = 'coocur in practice guidelines sentence'
SEMnoPG = '#Pubs with no practice guidelines'
COOCnoPG = 'coocur in no practice guidelines sentence'

def get_concept_id_hgraph(qstring):
    '''
    First approach to convert strings to IMUIs is using 
    the H-Graph REST API Concept Search Endpoint.
    Using a Cache to avoid repeated querying
    '''
    if qstring in IMUI_CACHE_CON:
        return IMUI_CACHE_CON[qstring]
    resp = hgraphApiProcessor.concept_search(qstring)
    result = resp['result']
    if len(result) > 0:
        concept = [{'imui': result[0]['imui'], 'label': result[0]['label'], 
                   'groupLabel': result[0]['groupLabel']}]
    else:
        concept = [{}]
    IMUI_CACHE_CON[qstring] = concept
    return concept

def get_concept_id_qpe(qstring):
    '''
    Second approach to convert strings to IMUIs is using 
    the H-Graph QPE API.
    Using a Cache to avoid repeated querying
    '''
    if qstring in IMUI_CACHE_QPE:
        return IMUI_CACHE_QPE[qstring]
    tags = QPEAccess.annotateString(qstring)
    concepts = []
    for k in tags:
        concepts.append({"imui": int(k.conceptID), "label": k.medicalName})
    IMUI_CACHE_QPE[qstring] = concepts
    return concepts

def get_concept_info(concept_imui):
    '''
    Get the concept information for selected concept
    Using a Cache to avoid repeated querying
    '''
    if concept_imui in IMUI_INFO:
        concept_info = IMUI_INFO[concept_imui]
    else:
        try:
            concept_info = hgraphApiProcessor.concept_info(concept_imui)
            IMUI_INFO[concept_imui] = concept_info
        except:
            concept_info = {}
    return concept_info

def get_subject_select_relations(concept_imui):
    #This function only retrieves the "object" concept of relations in H-Graph, where 
    #the relation type is selected under RELATIONS_OF_INTEREST
    concept_info = get_concept_info(concept_imui)
    relations = {}
    for k in RELATIONS_OF_INTEREST:
        relations[k] = []
        if 'subjectOfRelation' in concept_info:
            for m in concept_info['subjectOfRelation']:
                if m['label'] == k:
                    relations[k].append({'imui': m['object']['imui'], 'label': m['object']['label']})
    return relations
    
def process_execute_thread(argument):
    # Method to use for each thread
    try:
        (query_list, b, f) = argument 
        resp = f(query_list[b])
    except:
        resp = []
    return resp

def process_parallel(query_list:list, threads=3, f=get_concept_id_hgraph):
    #Function to query any APIs in parallel
    start_time = time.time()
    print("Started at time", start_time)
    ex = ThreadPoolExecutor(max_workers=threads)
    args = ((query_list, b, f) for b in range(len(query_list)))
    results = ex.map(process_execute_thread, args)
    real_results = list(results)
    print("Time taken", time.time() - start_time)
    if real_results:
        return real_results
    else: return []

def get_field(x, stype, field):
    try:
        return concept_matches[x][stype][0][field]
    except:
        return ''


def improve_query(name:str):
    n = name.replace(' measurement', '')
    n = n.replace(' testing', '')
    n = n.replace(' test', '')
    n = n.replace(' procedure', '')
    n = n.replace('serum/plasma ', '')
    n = n.replace(' mutation', '')
    n = n.replace(' gene', '')
    n = n.replace(' function', '')
    n = n.replace('play ', '')
    n = n.replace('  ', ' ')
    if n == 'poor': return 'low income'
    if n == 'FH': 
        return 'fumarate hydratase'
    if n == 'CA 125': 
        return 'MUC16'
    if n == 'female': return ''
    n = n.replace('/', ' ')

    return n.strip().lower()


def get_references_from_ETM(hgraph_name, disease_name, etm_operator:str, add_param:dict, annotate_rel_with:dict, hitcount_only=False):
    #hgraph_name = 'hysteroscopy'
    #disease_name = 'endometrial carcinoma'
    operator2evidence = {'sent':'Coocurence', 'rel':'Semantic'}
    sop2pubtype={'200~':'practice guidelines','dfb~':'excluding practice guidelines' }
    annotator = annotate_rel_with
    annotator.update({'Evidence':[operator2evidence[etm_operator]]})
    pub_type = sop2pubtype[add_param['so_p']]

    string4etm = improve_query(hgraph_name)
    if not string4etm: return []
    etm_dump_dir = 'D:/Python/ENTELLECT_API/Data/ETM/hGraph/'
    etm_stats_dir = 'D:/Python/ENTELLECT_API/Data/ETM/'

    etm_query = etm_operator + '({'+string4etm+'} AND {'+disease_name+'})'
    search_name = operator2evidence[etm_operator]+' search in ({t}) for ({n})-({d})'.format(t=pub_type,n=string4etm, d=disease_name)
    etm_miner = ETM(etm_query, search_name, APIconfig,etm_dump_dir,etm_stats_dir,add_param)
    if hitcount_only:
        return etm_miner.hit_count
    else:
        etm_miner.load_from_json(use_cache=True)
        updated_refs = list()
        for r in etm_miner.references:
            r.update(annotator)
            updated_refs.append(r)
    
        return updated_refs

def print_refs(fname:str, ref_counter:dict):
    sorted_refs = {k:v for k,v in sorted(ref_counter.items(),key=lambda x: len(x[1]),reverse=True)}
    with open(fname, 'w', encoding='utf-8') as f:
        f.write('Citation index\tCitation\tPMID or DOI\n')
        for identifier, refs in sorted_refs.items():
            biblio_str = refs[0].get_biblio_str()
            count = str(len(refs))
            f.write(str(count)+'\t'+biblio_str+'\t'+identifier+'\n')


def graph2excel(excel_fname,hgraph:ResnetGraph,diseases1st_urns:list):
    workbook = xlsxwriter.Workbook(excel_fname)
    pg_worksheet = workbook.add_worksheet('PracticeGuidelines')
    lit_worksheet = workbook.add_worksheet('ScientificLiterature')
    stat_columns = [SEMinPG,COOCinPG,SEMnoPG,COOCnoPG]
    header = ['Disease','Relation','hGraph concept'] + stat_columns + ['Relevant PMIDs or DOIs']
    header = ['Disease','Lab test'] + stat_columns + ['Relevant PMIDs or DOIs']
    for col in range(0,len(header)):
        pg_worksheet.write_string(0,col,header[col])
        lit_worksheet.write_string(0,col,header[col])

    pgw_row_counter = 0
    lit_row_counter = 0
    stat_columns = [SEMinPG,COOCinPG,SEMnoPG,COOCnoPG]
    col2operator = {SEMinPG:'rel',SEMnoPG:'rel',COOCinPG:'sent',COOCnoPG:'sent'}
    col2filter = {SEMinPG:'200~',COOCinPG:'200~',SEMnoPG:'dfb~',COOCnoPG:'dfb~'}
    baseURL = 'https://demo.elseviertextmining.com/advanced?'

    pg_ref_counter = dict() # {ID:[refs]}
    lit_ref_counter = dict() # {ID:[refs]}
    
    for node1id, node2id, rel in hgraph.edges.data('relation'):
        node1 = hgraph._get_node(node1id)
        node2 = hgraph._get_node(node2id)

        if node1['URN'][0] in diseases1st_urns:
            disease_node = node1
            neghbor_node = node2
        else:
            disease_node = node2
            neghbor_node = node1

        if int(rel[SEMinPG][0]) + int(rel[COOCinPG][0]) > 0:
            worksheet = pg_worksheet
            pgw_row_counter += 1
            row_counter = pgw_row_counter
            ref_counter = pg_ref_counter
        else:
            worksheet = lit_worksheet
            lit_row_counter +=1
            row_counter = lit_row_counter
            ref_counter = lit_ref_counter

        worksheet.write_string(row_counter, 0,disease_node['Name'][0])
        worksheet.write_string(row_counter, 1,rel['ObjTypeName'][0])
        worksheet.write_string(row_counter, 2,neghbor_node['Name'][0])
    
        string4etm = improve_query(neghbor_node['Name'][0])

        row_refs = list()
        for i in range(0,len(stat_columns)):
            col_name = stat_columns[i]
            col = i+3
            refcount = rel[col_name][0]

            # finding correct url to ETM
            if refcount > 0:
                etm_operator = col2operator[col_name]
                etm_query = etm_operator + '({'+string4etm+'} AND {'+disease_name+'})'
                filtr = col2filter[col_name]
                params = {'query':etm_query, 'mode':'AdvancedSearch', 'searchTarget':'full_index', 'so_p':filtr}
                data = parse.urlencode(params, quote_via=parse.quote)
                url = baseURL+data
            else:
                url = 'https://viz.graph.hmelsevier.com/concept/'+disease_node['hGraph ID'][0]

            formula = '=HYPERLINK(\"'+url+'",\"{refcount}\")'.format(refcount=str(refcount))
            #worksheet.write_formula(row_counter, col, formula)
            worksheet.write_string(row_counter, col ,str(refcount))

        sorted_refs = rel.sort_references('Relevance rank')
        max_range = 5 if len(sorted_refs) > 5 else len(sorted_refs)
        
        best_ref = list()
        for i in range(0,max_range):
            ref = sorted_refs[i]
            for id_type in PS_ID_TYPES:
                try:
                    identifier = ref.Identifiers[id_type]
                    best_ref.append(identifier)
                    try:
                        ref_counter[id_type+':'+identifier].append(ref)
                    except KeyError:
                        ref_counter[id_type+':'+identifier] = [ref]
                    break
                except KeyError:
                    continue
        best_refs_str = ','.join(best_ref)
        worksheet.write_string(row_counter, 3+len(stat_columns),best_refs_str)

        row_counter += 1
    workbook.close()

    print_refs(DATA_DIR+'hGraph Practice Guidlines References.txt',pg_ref_counter)
    print_refs(DATA_DIR+'hGraph Literature References.txt',lit_ref_counter)


if __name__ == "__main__":
    # Establish an H-GraphAPIProcessor object
    APIconfig = load_api_config()
    hgraphApiProcessor = HGraphAPIQueryProcessor(APIconfig['hGraphLic'])
    ps_api = APISession(APIconfig['ResnetURL'], APIconfig['PSuserName'], APIconfig['PSpassword'])
    prop2print = PS_BIBLIO_PROPS|set(SENTENCE_PROPS)
    rel_table_props = [SEMinPG,COOCinPG,SEMnoPG,COOCnoPG,'ObjTypeName']
    ent_table_props = ['Name']
    global_start = time.time()

    # Read the diseases input file
    disease_names_file = 'D:/Python/ENTELLECT_API/Data/hGraph/Focus uterine cancers from PS.txt'
    diseases_pandas = pd.read_csv(disease_names_file, keep_default_na=False, sep='\t')
    new_rel_id = 0

    concept_matches = {}
    row_counter = 0
    disease_urns = list(diseases_pandas['URN'])
    disease_urns_str = '\''+'\',\''.join(disease_urns)+'\''
    oql_query = 'SELECT Entity WHERE URN = ({urns})'.format(urns=disease_urns_str)
    ps_api.add_ent_props(['Alias','Connectivity'])
    disease_graph = ps_api.process_oql(oql_query, 'Find diseases in Resnet')
    ps_api.entProps = ['Name']
    diseases  = [PSObject(n) for i,n in disease_graph.nodes(data=True)]
    diseases.sort(key=lambda x: int(x['Connectivity'][0]))
    unmapped_log = open('unmapped hGraph concepts.txt', 'a')

    for disease in diseases:
        disease_node = PSObject(disease)
        row_counter += 1
        disease_name = disease_node['Name'][0]
        try:
            disease_aliases = list(disease_node['Alias'])
        except KeyError: disease_aliases = []
        disease_aliases.append(disease_name)
        all_names = list(set(map(lambda x: str(x).lower(),disease_aliases)))

        hgraph_name = disease_name+'_hgraph'
        hgraph = ResnetGraph()
        hgraph.name = hgraph_name   

        concept_list = process_parallel(all_names) # list of h-graph results
        #[{imui:ID, label:disease_name,group_lable:disease}]
        concept_list_qpe = []
        #concept_list_qpe = process_parallel(all_names, f=get_concept_id_qpe)
        all_matches = [x for x in concept_list+concept_list_qpe if x]
        imuis = set()
        for hgraph_match in all_matches:
            try:
                hgraph_match_label = str(hgraph_match[0]['label'])
                if hgraph_match_label.lower() in all_names:
                    try: imui = hgraph_match[0]['imui']
                    except KeyError: continue
                    imuis.add(str(imui))
            except KeyError:
                continue

        if not imuis: 
            continue
        imuis = list(imuis)
        disease_node.update_with_list('hGraph ID', imuis)

        all_concept_rels = process_parallel(imuis, f=get_subject_select_relations)
        all_concept_rels_clean = [x for x in all_concept_rels if x]

        unique_neighbors = set()
        if all_concept_rels_clean:
            hgraph.add1node(disease_node)
            neighbor_labels = set()
            for each_imui_relations in all_concept_rels:
                if isinstance(each_imui_relations,dict):
                    for rel_type, neighbors in each_imui_relations.items():
                        relation_type = rel_type
                        if rel_type == 'has risk factor':
                                for n in neighbors:
                                    if str(n['label']).find('mutation') > -1: 
                                        neighbor_labels.add(n['label'])
                                        unique_neighbors.add(('GeneticChange',n['imui'],n['label']))
                        else:
                            for n in neighbors:
                                neighbor_labels.add(n['label'])
                                unique_neighbors.add((rel_type,n['imui'],n['label']))

            improved_labels = list(map(lambda x: improve_query(x), neighbor_labels))
            label2objs, objid2labels = ps_api.map_props2objs(list(improved_labels),['Name','Alias'],case_insensitive=True)
            unique_neighbors = list(unique_neighbors)
            unique_neighbors.sort(key=lambda x: x[1])

            was_used = set()
            for neigh_tuple in unique_neighbors:
                relation_type = neigh_tuple[0]
                imui_neighbor = neigh_tuple[1]
                neighbor_name = neigh_tuple[2]
                try:
                    ps_neighbors = label2objs[improve_query(neighbor_name)]
                except KeyError:
                    neighbor_node = PSObject(dict())
                    neighbor_node.append_property('Id',imui_neighbor)
                    neighbor_node.append_property('Name',neighbor_name)
                    obj_type = RELTYPE2ENTYPE_MAP[rel_type]
                    neighbor_node.append_property('ObjTypeName',obj_type)
                    neighbor_node.append_property('URN','urn:agi-hgraph:'+str(imui_neighbor))
                    ps_neighbors = [neighbor_node]
                    unmapped_log.write(neighbor_name+'\t'+str(imui_neighbor)+'\n')

                for neighbor in ps_neighbors:
                    #different labels may find the same PSObject in Resnet
                    if neighbor in was_used: continue
                    else: was_used.add(neighbor)

                    neighbor.append_property('hGraph ID',imui_neighbor)
                    rel_obj_type = str(relation_type).replace(' ','_')
                    ps_rel = PSRelation({'ObjTypeName':[rel_obj_type],
                                        'Id':[new_rel_id],
                                        'URN' :'urn:agi-FunctionalAssociation:'+str(new_rel_id)
                                        })
                    
                    textref = 'info:hgraph:'+str(new_rel_id)
                    ref = Reference.from_textref(textref)
                    ref.add_sentence_prop(textref,'Sentence','Imported from Elsevier hGraph')
                    ref['Source'] = ['hGraph']
                    ps_rel.add_references([ref])

                    guidlines_refs = set(get_references_from_ETM(neighbor['Name'][0],disease_node['Name'][0],'rel', 
                                add_param = {'so_p':'200~'}, annotate_rel_with={'PubTypes':['Practice guidelines']}))
                    ps_rel[SEMinPG] = [len(guidlines_refs)] #hit_count uses 1 based index

                    if guidlines_refs:
                        hitcount_only = True 
                        ps_rel.add_references(guidlines_refs)
                        ps_rel[COOCinPG] = [get_references_from_ETM(neighbor['Name'][0],disease_node['Name'][0],'sent', 
                                add_param = {'so_p':'200~'}, annotate_rel_with={'PubTypes':['Practice guidelines']},hitcount_only=True)]
                        ps_rel[SEMnoPG] = [get_references_from_ETM(neighbor['Name'][0],disease_node['Name'][0],'rel', 
                                add_param = {'so_p':'dfb~'}, annotate_rel_with={'PubTypes':['Excludes practice guidelines']},hitcount_only=True)]
                        ps_rel[COOCnoPG] = [get_references_from_ETM(neighbor['Name'][0],disease_node['Name'][0],'sent', 
                                add_param = {'so_p':'dfb~'}, annotate_rel_with={'PubTypes':['Excludes practice guidelines']},hitcount_only=True)]    
                        ps_rel['RelationNumberOfReferences'] = [len(guidlines_refs)]                    
                    else:
                        guidlines_refs.update(get_references_from_ETM(neighbor['Name'][0],disease_node['Name'][0],'sent', 
                                    add_param = {'so_p':'200~'}, annotate_rel_with={'PubTypes':['Practice guidelines']}))
                        ps_rel[COOCinPG] = [len(guidlines_refs)]
                        ps_rel.add_references(guidlines_refs)
                        other_refs = set(get_references_from_ETM(neighbor['Name'][0],disease_node['Name'][0],'rel', 
                                    add_param = {'so_p':'dfb~'}, annotate_rel_with={'PubTypes':['Excludes practice guidelines']}))
                        ps_rel[SEMnoPG] = [len(other_refs)]
                        if other_refs:
                            ps_rel.add_references(other_refs)
                            ps_rel[COOCnoPG] = [get_references_from_ETM(neighbor['Name'][0],disease_node['Name'][0],'sent', 
                                    add_param = {'so_p':'dfb~'}, annotate_rel_with={'PubTypes':['Excludes practice guidelines']},hitcount_only=True)]
                            ps_rel['RelationNumberOfReferences'] = [len(other_refs)]  
                        else:
                            other_refs.update(get_references_from_ETM(neighbor['Name'][0],disease_node['Name'][0],'sent', 
                                    add_param = {'so_p':'dfb~'}, annotate_rel_with={'PubTypes':['Excludes practice guidelines']}))
                            ps_rel[COOCnoPG] = [len(other_refs)]
                            ps_rel['RelationNumberOfReferences'] = [len(other_refs)] 
                            ps_rel.add_references(other_refs)                      
                    
                    #all_refs = guidlines_refs|other_refs
                    #ps_rel['RelationNumberOfReferences'] = [len(all_refs)]

                    ps_rel.Nodes['Regulators'] = [(disease_node['Id'][0], '0', '')]
                    ps_rel.Nodes['Regulators'].append((neighbor['Id'][0], '0', ''))

                    hgraph.add1node(neighbor)
                    hgraph.add_edge(disease_node['Id'][0],neighbor['Id'][0],relation=ps_rel, weight=float(ps_rel['RelationNumberOfReferences'][0]))
                    
                    new_rel_id += 1
        else:
            continue

    unmapped_log.close()
    
    hgraph_rnef_str = hgraph.to_rnef(prop2print)
    hgraph_rnef_str = ps_api.pretty_xml(hgraph_rnef_str,no_declaration=True)
    rnef_file =  open(DATA_DIR+hgraph_name+'.rnef','w',encoding='utf-8')
    rnef_file.write('<batch>\n'+hgraph_rnef_str+'\n</batch>')
    rnef_file.close()

    hgraph.print_references(DATA_DIR+disease_name+'_hgraph_stats.tsv',rel_table_props,ent_table_props,single_rel_row=True)
    
    excel_fname = DATA_DIR+'hgraph_stats.xlsx'
    graph2excel(excel_fname,hgraph, disease_urns)

    print('Job was done in %s' % ps_api.execution_time(global_start))

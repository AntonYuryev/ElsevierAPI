from ElsevierAPI import open_api_session
import pandas as pd
from ElsevierAPI.ResnetAPI.PathwayStudioGOQL import OQL
import xlsxwriter
import urllib.parse
import networkx as nx
from  ElsevierAPI.ResnetAPI.Resnet2rdf import ResnetGraph, ResnetRDF

DATA_DIR = 'D:/Python/Quest/report tables/'
api_config = 'path2apiconfig.json'
api_config = 'D:/Python/ENTELLECT_API/ElsevierAPI/APIconfigMDACC.json'
ps_api = open_api_session(api_config) # specify here path to your APIconfig file. Defaults to ./ElsevierAPI/APIconfig.json
ps_api.add_rel_props(['PMID'])
ps_api.PageSize = 1000
disease_names_file = 'D:/Python/ENTELLECT_API/Data/hGraph/Focus uterine cancers from PS.txt'
diseases_pandas = pd.read_csv(disease_names_file, keep_default_na=False, sep='\t')
disease_urns = list(diseases_pandas['URN'])
oql_query = OQL.expand_entity(disease_urns,['URN'], expand2neighbors=['GeneticVariant'])
disease_graph = ps_api.process_oql(oql_query, 'Find GVs linked to diseases')
add_hyperlinks = False

gv_ids = disease_graph.get_entity_ids(['GeneticVariant'])
prot2gvs_graph = ResnetGraph()
for i in range(0, len(gv_ids),1000):
    chunk = gv_ids[i: i+1000]
    oql_query = OQL.expand_entity(PropertyValues=chunk, SearchByProperties=['id'], 
                        expand_by_rel_types=['GeneticChange'],expand2neighbors=['Protein'])
    prot2gvs_graph.add_graph(ps_api.process_oql(oql_query, 'Find genes containing GVs linked to diseases'))

#annotate GVs with property 'Gene name' before printing
gvid2genes = dict()
for gv_id, protein_id, rel in prot2gvs_graph.edges.data('relation'):
    protein_node = prot2gvs_graph._get_node(protein_id)
    gv_node = prot2gvs_graph._get_node(gv_id)
    gene_name = protein_node['Name'][0]
    try:
        gvid2genes[gv_id].append(gene_name)
    except KeyError:
        gvid2genes[gv_id] = [gene_name]
    

#report_pandas = pd.DataFrame(columns=['Disease','GV','#references'])
table2sort = pd.DataFrame()
workbook = xlsxwriter.Workbook(DATA_DIR+'GV report.xlsx')
worksheet = workbook.add_worksheet('GVs2disease')
header = ['Disease','Gene name','Variant','#references']
for col in range(0,len(header)):
    worksheet.write_string(0,col,header[col])

row_counter = 1
#disease_graph4rdf = ResnetGraph()
for node1id, node2id, rel in disease_graph.edges.data('relation'):
    node1 = disease_graph._get_node(node1id)
    node2 = disease_graph._get_node(node2id)
    if node1['ObjTypeName'][0]=='Disease':
        disease_node = node1
        gv_node = node2
    else:
        disease_node = node2
        gv_node = node1

    worksheet.write_string(row_counter, 0,disease_node['Name'][0])
    
    gene_names = list()
    try:
        gene_names = gvid2genes[gv_node['Id'][0]]
    except KeyError: pass
    #rel.update_with_list('Gene',gene_names)
    #rel.load_references()
    #disease_graph4rdf.add_triple(node1,node2,rel)
    nx.set_node_attributes(disease_graph, {gv_node['Id'][0]:{'Gene':gene_names}})
    gene_names  = ','.join(gene_names)
    worksheet.write_string(row_counter, 1,gene_names)
    
    if add_hyperlinks:
        formula = '=HYPERLINK("https://www.ncbi.nlm.nih.gov/snp/{gv}",\"{gv}\")'.format(gv=gv_node['Name'][0])
        worksheet.write_formula(row_counter, 2, formula)
    else:
        worksheet.write_string(row_counter, 2, gv_node['Name'][0])

    rel_pmids = set()
    refcount = str(rel['RelationNumberOfReferences'][0])

    for prop_set in rel.PropSetToProps.values():
        for k,v in prop_set.items():
            if k == 'PMID':
                rel_pmids.update(v)

    pmids_str = ','.join(list(rel_pmids)[:10])
    if add_hyperlinks:
        base_url = 'https://pubmed.ncbi.nlm.nih.gov/?'
        params = {'term':pmids_str}
        data = urllib.parse.urlencode(params, quote_via=urllib.parse.quote)
        formula = '=HYPERLINK("'+base_url+data+'",\"{refcount}\")'.format(refcount=refcount)
        ref_column_name = '#References'
        ref_col_value = refcount
        worksheet.write_formula(row_counter, 3, formula)
    else:
        worksheet.write_string(row_counter, 3, pmids_str)
        ref_column_name = 'References'
        ref_col_value = pmids_str

    table2sort.at[row_counter,'Disease'] = disease_node['Name'][0]
    table2sort.at[row_counter,'Gene name'] = gene_names
    table2sort.at[row_counter,'GV'] = str(gv_node['Name'][0])
    table2sort.at[row_counter,ref_column_name] = ref_col_value

    row_counter += 1

workbook.close()
ResnetRDF.fromResnetGraph(disease_graph).to_json(DATA_DIR+'GVs2disease.jsonld')
table2sort = table2sort.sort_values(by=['Gene name','Disease'])
table2sort.to_csv(DATA_DIR+'GVs2disease.tsv',sep='\t', index=False)

        

from ElsevierAPI import open_api_session
import pandas as pd
from ElsevierAPI.ResnetAPI.PathwayStudioGOQL import OQL
import xlsxwriter
import urllib.parse
from ElsevierAPI.ResnetAPI.ResnetGraph import ResnetGraph

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
workbook = xlsxwriter.Workbook(DATA_DIR+'GV report.xlsx')
worksheet = workbook.add_worksheet('GVs2disease')
header = ['Disease','Gene name','GV','#references']
for col in range(0,len(header)):
    worksheet.write_string(0,col,header[col])

row_counter = 1
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
    worksheet.write_string(row_counter, 1,','.join(gene_names))

    formula = '=HYPERLINK("https://www.ncbi.nlm.nih.gov/snp/{gv}",\"{gv}\")'.format(gv=gv_node['Name'][0])
    worksheet.write_formula(row_counter, 2, formula)

    rel_pmids = set()
    refcount = str(rel['RelationNumberOfReferences'][0])
    for prop_set in rel.PropSetToProps.values():
        for k,v in prop_set.items():
            if k == 'PMID':
                rel_pmids.update(v)

    rel_pmids = list(rel_pmids)[:20]
    base_url = 'https://pubmed.ncbi.nlm.nih.gov/?'
    params = {'term':','.join(rel_pmids)}
    data = urllib.parse.urlencode(params, quote_via=urllib.parse.quote)
    formula = '=HYPERLINK("'+base_url+data+'",\"{refcount}\")'.format(refcount=refcount)
    worksheet.write_formula(row_counter, 3, formula)

    row_counter += 1
 
workbook.close()




        

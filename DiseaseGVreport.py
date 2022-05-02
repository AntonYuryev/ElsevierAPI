from ElsevierAPI import open_api_session
import pandas as pd
from ElsevierAPI.ResnetAPI.PathwayStudioGOQL import OQL
import xlsxwriter
import urllib.parse
import networkx as nx
from  ElsevierAPI.ResnetAPI.Resnet2rdf import ResnetGraph,ResnetRDF,PSObject, PSRelation
from ElsevierAPI.ETM_API.references import PUBYEAR

DATA_DIR = 'D:/Python/Quest/report tables/'
api_config = str()
ps_api = open_api_session(api_config) # specify here path to your APIconfig file. Defaults to ./ElsevierAPI/APIconfig.json
ps_api.add_rel_props(['PMID',PUBYEAR])
ps_api.PageSize = 1000


def read_disease_urns(file_name:str):
    # one column in file must be called "URN"
    diseases_pandas = pd.read_csv(file_name, keep_default_na=False, sep='\t')
    return list(diseases_pandas['URN'])


def get_disease_childs(disease_name:str):
    ontology_branch_graph = ps_api.child_graph([disease_name],['Name','Alias'])
    branch_nodes = [PSObject(node) for i,node in ontology_branch_graph.nodes(data=True)]
    return [x['URN'][0] for x in branch_nodes]


def map_gvid2genes(disease2gvs:ResnetGraph):
    gv_ids = disease2gvs.get_entity_ids(['GeneticVariant'])
    prot2gvs_graph = ResnetGraph()
    print ('Finding genes for %d genetic variants' % len(gv_ids))
    number_of_iterations = int(len(gv_ids)/1000)+1
    for i in range(0, len(gv_ids),1000):
        chunk = gv_ids[i: i+1000]
        oql_query = OQL.expand_entity(PropertyValues=chunk, SearchByProperties=['id'], 
                            expand_by_rel_types=['GeneticChange'],expand2neighbors=['Protein'])
        request_name = '{iter} iteration out of {iternum} to find genes linked to GVs'.format(iter=str(int(i/1000)+1),iternum=number_of_iterations)
        prot2gvs_graph.add_graph(ps_api.process_oql(oql_query,request_name))

    # making gvid2genes for subsequent annotation
    gvid2genes = dict()
    for gv_id, protein_id, rel in prot2gvs_graph.edges.data('relation'):
        protein_node = prot2gvs_graph._get_node(protein_id)
        gene_name = protein_node['Name'][0]
        try:
            gvid2genes[gv_id].append(gene_name)
        except KeyError:
            gvid2genes[gv_id] = [gene_name]

    [nx.set_node_attributes(disease2gvs, {gvid:{'Gene':gene_names}}) for gvid, gene_names in gvid2genes.items()]


if __name__ == "__main__":
    disease_name = 'diabetes mellitus'
    disease_dir = DATA_DIR+disease_name+'/'
    #disease_urns = read_disease_urns('D:/Python/ENTELLECT_API/Data/hGraph/Focus uterine cancers from PS.txt')
    disease_urns = get_disease_childs(disease_name)

    oql_query = OQL.expand_entity(disease_urns,['URN'], expand2neighbors=['GeneticVariant'])
    disease2gvs = ps_api.process_oql(oql_query, 'Find GVs linked to diseases')
    add_hyperlinks = False
    map_gvid2genes(disease2gvs)

    table2sort = pd.DataFrame(columns=['Disease','Gene name','GV','#References','References'])
    workbook = xlsxwriter.Workbook(disease_dir+disease_name+' GV report.xlsx')
    worksheet = workbook.add_worksheet('GVs2disease')
    header = ['Disease','Gene name','Variant','#references']
    for col in range(0,len(header)):
        worksheet.write_string(0,col,header[col])

    row_counter = 1
    for node1id, node2id, rel in disease2gvs.edges.data('relation'):
        node1 = disease2gvs._get_node(node1id)
        node2 = disease2gvs._get_node(node2id)
        assert isinstance(rel,PSRelation)
        if node1['ObjTypeName'][0]=='Disease':
            disease_node = node1
            gv_node = node2
        else:
            disease_node = node2
            gv_node = node1

        worksheet.write_string(row_counter,0,disease_node['Name'][0]) # column #1
        
        try:
            gene_names = ','.join(gv_node['Gene'])
        except KeyError: 
            gene_names = ''

        worksheet.write_string(row_counter, 1,gene_names) # column #2
        
        if add_hyperlinks:
            formula = '=HYPERLINK("https://www.ncbi.nlm.nih.gov/snp/{gv}",\"{gv}\")'.format(gv=gv_node['Name'][0])
            worksheet.write_formula(row_counter, 2, formula) # column #3
        else:
            worksheet.write_string(row_counter, 2, gv_node['Name'][0]) # column #3


        refcount = rel.get_reference_count()
        last_refs = rel.recent_refs(5)
        rel_pmids = [x.get_doc_id()[1] for x in last_refs]
        ids_str = ','.join(rel_pmids)
        
        if add_hyperlinks:
            base_url = 'https://pubmed.ncbi.nlm.nih.gov/?'
            params = {'term':ids_str}
            data = urllib.parse.urlencode(params, quote_via=urllib.parse.quote)
            formula = '=HYPERLINK("'+base_url+data+'",\"{refcount}\")'.format(refcount=refcount)
            ref_column_name = '#References'
            ref_col_value = refcount
            worksheet.write_formula(row_counter, 3, formula) # column #4
        else:
            worksheet.write_string(row_counter, 3, ids_str) # column #4
            ref_column_name = 'References'
            ref_col_value = ids_str


        table2sort.loc[row_counter] = [disease_node['Name'][0],gene_names,str(gv_node['Name'][0]),int(rel['RelationNumberOfReferences'][0]),ref_col_value]
     #   table2sort.at[row_counter,'Disease'] = disease_node['Name'][0]
     #   table2sort.at[row_counter,'Gene name'] = gene_names
     #   table2sort.at[row_counter,'GV'] = str(gv_node['Name'][0])
     #   table2sort.at[row_counter,'#References'] = refcount
     #   table2sort.at[row_counter,ref_column_name] = ref_col_value
        row_counter += 1

    workbook.close()

    ResnetRDF.fromResnetGraph(disease2gvs).to_json(disease_dir+disease_name+' GVs.jsonld')
    #table2sort = table2sort.sort_values(by=['Gene name','Disease'])
    table2sort.sort_values(by=['#References'],ascending=False, inplace=True)
    table2sort.to_csv(disease_dir+disease_name+' GVs.tsv',sep='\t', index=False)


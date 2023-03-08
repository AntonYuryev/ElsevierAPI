from ElsevierAPI import open_api_session
from ElsevierAPI.pandas.panda_tricks import df,ExcelWriter,NUMBER_OF_REFERENCE
from ElsevierAPI.ResnetAPI.PathwayStudioGOQL import OQL
import networkx as nx
from  ElsevierAPI.ResnetAPI.Resnet2rdf import ResnetGraph,ResnetRDF,PSObject,PSRelation
from ElsevierAPI.ETM_API.references import PUBYEAR
from ElsevierAPI.NCBI.pubmed import pubmed_hyperlink, minor_allele
from ElsevierAPI.ETM_API.etm import ETM_REFS_COLUMN

api_config = str()
ps_api = open_api_session(api_config)
ps_api.add_rel_props(['PMID',PUBYEAR])
ps_api.PageSize = 1000


def get_disease_childs(disease_name:str):
    ontology_branch_graph = ps_api.child_graph([disease_name],['Name','Alias'])
    branch_nodes = [PSObject(node) for i,node in ontology_branch_graph.nodes(data=True)]
    return [x.urn() for x in branch_nodes]


def map_gvid2genes(disease2gvs:ResnetGraph):
    gv_ids = disease2gvs.dbids4nodes(['GeneticVariant'])
    gvid2genes = ps_api.gv2gene(gv_ids)
    [nx.set_node_attributes(disease2gvs, {gvid:{'Gene':gene_names}}) for gvid, gene_names in gvid2genes.items()]


DATA_DIR = 'D:/Python/PMI/'

if __name__ == "__main__":
    disease_name = 'Pulmonary Hypertension' 
    disease_dir = DATA_DIR+disease_name+'/'
    disease_urns = get_disease_childs(disease_name)

    oql_query = OQL.expand_entity(disease_urns,['URN'], expand2neighbors=['GeneticVariant'])
    disease2gvs = ps_api.process_oql(oql_query, 'Find GVs linked to diseases')
    map_gvid2genes(disease2gvs)

    header = ['Disease','Gene name','Variant',NUMBER_OF_REFERENCE,ETM_REFS_COLUMN,'Minor allele','Allele frequency']
    table2sort = df(columns=header)

    rs_ids = [o['Name'][0] for i,o in disease2gvs.nodes(data=True) if o['ObjTypeName'][0] == 'GeneticVariant']
    snpid2allele2freq = minor_allele(rs_ids)

    row_counter = 1
    for node1id, node2id, rel in disease2gvs.edges.data('relation'):
        node1 = disease2gvs._psobj(node1id)
        node2 = disease2gvs._psobj(node2id)
        assert isinstance(rel,PSRelation)
        if node1['ObjTypeName'][0]=='Disease':
            disease_node = node1
            gv_node = node2
        else:
            disease_node = node2
            gv_node = node1

        try:
            minoralleles = dict(snpid2allele2freq[gv_node['Name'][0]])
        except KeyError: continue

        annotation_str = [allele+'-'+str(freq) for allele,freq in minoralleles.items()]
        nx.set_node_attributes(disease2gvs, {gv_node['Id'][0]:{'Minor allele frequencies':annotation_str}})

        for allele, freq in minoralleles.items():
            try:
                gene_names = ','.join(gv_node['Gene'])
            except KeyError: 
                gene_names = ''

            rs_name = '=HYPERLINK("https://www.ncbi.nlm.nih.gov/snp/{gv}",\"{gv}\")'.format(gv=gv_node['Name'][0])

            refcount = rel.get_reference_count()
            last_refs = rel._get_refs(ref_limit=5)
            rel_pmids = [x.get_doc_id()[1] for x in last_refs]
            
            ref_links = pubmed_hyperlink(rel_pmids,refcount)

            table2sort.loc[row_counter] = [disease_node['Name'][0],gene_names,rs_name,refcount,ref_links,allele,freq]
    
            row_counter += 1

    ResnetRDF.fromResnetGraph(disease2gvs).to_json(disease_dir+disease_name+' GVs.jsonld')
    table2sort.sort_values(by=[NUMBER_OF_REFERENCE,'Allele frequency'],ascending=[False,True],inplace=True)
    table2sort.make_header_vertical()
    table2sort.set_hyperlink_color(['Variant',ETM_REFS_COLUMN])
    table2sort.add_column_format(ETM_REFS_COLUMN,'align','center')
    table2sort.add_column_format(NUMBER_OF_REFERENCE,'align','center')
    table2sort.add_column_format('Disease','width',50)
    workbook = ExcelWriter(disease_dir + disease_name + ' GV report.xlsx', engine='xlsxwriter')
    table2sort.df2excel(workbook,'GVs2disease')
    workbook.save()

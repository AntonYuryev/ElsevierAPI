from ENTELLECT_API.ElsevierAPI.pandas.panda_tricks import df,ExcelWriter,NUMBER_OF_REFERENCE
from ENTELLECT_API.ElsevierAPI.ResnetAPI.PathwayStudioGOQL import OQL
import networkx as nx
from ENTELLECT_API.ElsevierAPI.ResnetAPI.Resnet2rdf import ResnetGraph,ResnetRDF,PSObject,PSRelation
from ENTELLECT_API.ElsevierAPI.ResnetAPI.NetworkxObjects import OBJECT_TYPE
from ENTELLECT_API.ElsevierAPI.ETM_API.references import PUBYEAR
from ENTELLECT_API.ElsevierAPI.NCBI.pubmed import pubmed_hyperlink
from ENTELLECT_API.ElsevierAPI.NCBI.dbsnp import minor_allele
from ENTELLECT_API.ElsevierAPI.ETM_API.RefStats import REFCOUNT_COLUMN
from ENTELLECT_API.ElsevierAPI.ResnetAPI.ResnetAPISession import APISession


REFLIMIT = 5

ps_api = APISession()
ps_api.add_rel_props(['PMID',PUBYEAR])
ps_api.PageSize = 1000


def get_disease_childs(disease_name:str):
    ontology_branch_graph = ps_api.child_graph([disease_name],['Name','Alias'])
    branch_nodes = [PSObject(node) for i,node in ontology_branch_graph.nodes(data=True)]
    return [x.urn() for x in branch_nodes]


def map_gvid2genes(disease2gvs:ResnetGraph):
    GVs = disease2gvs._psobjs_with('GeneticVariant',OBJECT_TYPE)
    gvid2genes = ps_api.gv2gene(GVs)
    [nx.set_node_attributes(disease2gvs, {gvid:{'Gene':gene_names}}) for gvid, gene_names in gvid2genes.items()]


DATA_DIR = 'D:/Python/ThermoFisher/'

if __name__ == "__main__":
    disease_name = 'acute rejection' 
    disease_dir = DATA_DIR
  #  disease_dir += disease_name + '/'
    disease_urns = get_disease_childs(disease_name)

    oql_query = OQL.expand_entity(disease_urns,['URN'], expand2neighbors=['GeneticVariant'])
    disease2gvs = ps_api.process_oql(oql_query, 'Find GVs linked to diseases')
    ps_api.annotate_gv_with_genes(disease2gvs)

    header = ['Disease','Gene name','Variant',NUMBER_OF_REFERENCE,REFCOUNT_COLUMN,'Minor allele','Allele frequency']
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

            refcount = rel.count_refs()
            last_refs = rel.refs(REFLIMIT)
            rel_pmids = [x.get_doc_id()[1] for x in last_refs]
            
            ref_links = pubmed_hyperlink(rel_pmids,refcount)

            table2sort.loc[row_counter] = [disease_node['Name'][0],gene_names,rs_name,refcount,ref_links,allele,freq]
    
            row_counter += 1

    ResnetRDF.fromResnetGraph(disease2gvs).to_json(disease_dir+disease_name+' GVs.jsonld')
    table2sort.sort_values(by=[NUMBER_OF_REFERENCE,'Allele frequency'],ascending=[False,True],inplace=True)
    table2sort.make_header_vertical()
    table2sort.set_hyperlink_color(['Variant',REFCOUNT_COLUMN])
    table2sort.add_column_format(REFCOUNT_COLUMN,'align','center')
    table2sort.add_column_format(NUMBER_OF_REFERENCE,'align','center')
    table2sort.add_column_format('Disease','width',50)
    workbook = ExcelWriter(disease_dir + disease_name + ' GV report.xlsx', engine='xlsxwriter')
    table2sort.df2excel(workbook,'GVs2disease')
    workbook.close()

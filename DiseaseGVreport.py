from ElsevierAPI import open_api_session
import pandas as pd
from ElsevierAPI.ResnetAPI.PathwayStudioGOQL import OQL
import xlsxwriter
import urllib.parse
import networkx as nx
from  ElsevierAPI.ResnetAPI.Resnet2rdf import ResnetGraph,ResnetRDF,PSObject,PSRelation
from ElsevierAPI.ETM_API.references import PUBYEAR
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
from ElsevierAPI.NCBI.pubmed import pubmed_hyperlink


api_config = str()
#api_config = 'D:/Python/ENTELLECT_API/ElsevierAPI/APIconfigMDACC.json'
ps_api = open_api_session(api_config)
ps_api.add_rel_props(['PMID',PUBYEAR])
ps_api.PageSize = 1000


def read_disease_urns(file_name:str):
    # one column in file must be called "URN"
    diseases_pandas = pd.read_csv(file_name, keep_default_na=False, sep='\t')
    return list(diseases_pandas['URN'])


def get_disease_childs(disease_name:str):
    ontology_branch_graph = ps_api.child_graph([disease_name],['Name','Alias'])
    branch_nodes = [PSObject(node) for i,node in ontology_branch_graph.nodes(data=True)]
    return [x.urn() for x in branch_nodes]


def map_gvid2genes(disease2gvs:ResnetGraph):
    gv_ids = disease2gvs.get_node_ids(['GeneticVariant'])
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


def minor_allele(rsids:list):
    stepSize = 200
    params = {'db':'snp'}
    snp2allele2freq = dict() # {snp_id:{allele:(allele_count, pop_size)}}
    for i in range(0, len(rsids), stepSize):
        ids = ','.join(s[2:] for s in rsids[i:i+stepSize])

        baseURL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
        params.update({'id':ids})
        req = urllib.request.Request(url=baseURL+'efetch.fcgi?'+urllib.parse.urlencode(params))
        
        response = urllib.request.urlopen(req).read()
        snps = ET.fromstring('<documents>'+response.decode()+'</documents>')
        for snp_record in snps.findall('./DocumentSummary'):
            alleles = dict() # {allele:(allele_count, pop_size)}
            snp_id = snp_record.find('./SNP_ID')
            if type(snp_id) == type(None): continue
            snp_id = snp_record.find('./SNP_ID').text
            for maf in snp_record.findall('./GLOBAL_MAFS/MAF'):
                study = maf.find('STUDY').text
                allele_freq_pop = maf.find('FREQ').text
                eq_pos = allele_freq_pop.find('=')
                slash_pos = allele_freq_pop.find('/', eq_pos)
                allele = allele_freq_pop[:eq_pos]
                population_size = int(allele_freq_pop[slash_pos+1:])
                freq = float(allele_freq_pop[eq_pos+1:slash_pos])
                try:
                    allele_count, pop_size = alleles[allele]
                    alleles[allele] = (allele_count+freq*population_size, pop_size+population_size)
                except KeyError:
                    alleles[allele] = (freq*population_size,population_size)

            alele_freqs = {a:f[0]/f[1] for a,f in alleles.items() if f[1]>0}
            if alele_freqs:
                snp2allele2freq['rs'+snp_id] = alele_freqs

    return snp2allele2freq

DATA_DIR = 'D:/Python/PMI/'

if __name__ == "__main__":
    disease_name = 'Myocardial Infarction'#'Atrial Fibrillation'#'Pulmonary Hypertension' #'diabetes mellitus' #'uterine cancer'# 
    disease_dir = DATA_DIR+disease_name+'/'
    #disease_urns = read_disease_urns('D:/Python/ENTELLECT_API/Data/hGraph/Focus uterine cancers from PS.txt')
    disease_urns = get_disease_childs(disease_name)

    oql_query = OQL.expand_entity(disease_urns,['URN'], expand2neighbors=['GeneticVariant'])
    disease2gvs = ps_api.process_oql(oql_query, 'Find GVs linked to diseases')
    add_hyperlinks = False
    map_gvid2genes(disease2gvs)

    table2sort = pd.DataFrame(columns=['Disease','Gene name','GV','#References','References','Minor allele','Allele frequency'])
    workbook = xlsxwriter.Workbook(disease_dir+disease_name+' GV report.xlsx')
    worksheet = workbook.add_worksheet('GVs2disease')
    header = ['Disease','Gene name','Variant','#references']
    for col in range(0,len(header)):
        worksheet.write_string(0,col,header[col])

    rs_ids = [o['Name'][0] for i,o in disease2gvs.nodes(data=True) if o['ObjTypeName'][0] == 'GeneticVariant']
    snpid2allele2freq = minor_allele(rs_ids)

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

        try:
            minoralleles = snpid2allele2freq[gv_node['Name'][0]]
        except KeyError: continue
        annotation_str = [allele+'-'+str(freq) for allele,freq in minoralleles.items()]
        nx.set_node_attributes(disease2gvs, {gv_node['Id'][0]:{'Minor allele frequencies':annotation_str}})

        for allele, freq in minoralleles.items():
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
            last_refs = rel._get_refs(ref_limit=5)
            rel_pmids = [x.get_doc_id()[1] for x in last_refs]
            
            if add_hyperlinks:
                formula = pubmed_hyperlink(rel_pmids,refcount)
                ref_col_value = refcount
                worksheet.write_formula(row_counter, 3, formula) # column #4
            else:
                ref_col_value = ','.join(rel_pmids)
                worksheet.write_string(row_counter, 3, ref_col_value) # column #4
                
            worksheet.write_string(row_counter, 3, allele) # column #5
            worksheet.write_string(row_counter, 3, str(freq)) # column #6

            table2sort.loc[row_counter] = [disease_node['Name'][0],gene_names,str(gv_node['Name'][0]),
                                        int(rel['RelationNumberOfReferences'][0]),ref_col_value,allele,freq]
    
            row_counter += 1

    workbook.close()

    ResnetRDF.fromResnetGraph(disease2gvs).to_json(disease_dir+disease_name+' GVs.jsonld')
    #table2sort = table2sort.sort_values(by=['Gene name','Disease'])
    table2sort.sort_values(by=['#References','Allele frequency'],ascending=[False,True], inplace=True)
    table2sort.to_csv(disease_dir+disease_name+' GVs.tsv',sep='\t', index=False, float_format='%g')
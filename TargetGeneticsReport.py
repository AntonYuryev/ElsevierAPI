from ElsevierAPI import load_api_config
import pandas as pd
from ElsevierAPI.ETM_API.references import PUBYEAR
from ElsevierAPI.ResnetAPI.PathwayStudioGOQL import OQL
from ElsevierAPI.ResnetAPI.NetworkxObjects import GENETICVARIANT,FUNC_ASSOC,PS_ID_TYPES
from ElsevierAPI.NCBI.dbsnp import minor_allele,dbsnp_hyperlink
from ElsevierAPI.NCBI.pubmed import pubmed_hyperlink


from ElsevierAPI.ResnetAPI.ResnetAPISession import APISession
from ElsevierAPI.ResnetAPI.ResnetGraph import ResnetGraph
DATA_DIR = 'D:/Python/PMI/'

class TargetGenetics(APISession):
    pass
    def __init__(self, APIconfig):
        super().__init__(APIconfig['ResnetURL'], APIconfig['PSuserName'], APIconfig['PSpassword'])
        self.flush_dump_files()
        self.PageSize = 1000
        self.snp2allele2freq = dict()
        self.gv2disease_graph = ResnetGraph()
        self.gv_ids = list()
        
        self.add_rel_props(PS_ID_TYPES+PUBYEAR)

    def load_graph(self,target_name:str):
        oql_query = OQL.expand_entity(PropertyValues=[target_name], SearchByProperties=['Name','Alias'], expand_by_rel_types=['GeneticChange'],
                        expand2neighbors=['Disease'])
        request_name = 'Find diseases genetically linked to {}'.format(target_name)
        self.process_oql(oql_query,request_name)

        oql_query = OQL.expand_entity(PropertyValues=[target_name], SearchByProperties=['Name','Alias'], expand_by_rel_types=['GeneticChange'],
                        expand2neighbors=[GENETICVARIANT])
        request_name = 'Find GVs for {}'.format(target_name)
        target2GV_graph = self.process_oql(oql_query,request_name)


        self.gv_ids = target2GV_graph.get_entity_ids([GENETICVARIANT])
        oql_query = OQL.expand_entity(PropertyValues=self.gv_ids, SearchByProperties=['id'], expand_by_rel_types=[FUNC_ASSOC],
                        expand2neighbors=['Disease'])

        request_name = 'Find diseases for {} GVs'.format(str(len(self.gv_ids)))
        self.gv2disease_graph = self.process_oql(oql_query,request_name)

        rs_ids = target2GV_graph.get_properties(self.gv_ids,'Name')
        rs_ids = [v[0] for v in rs_ids.values()]
        self.snp2allele2freq = minor_allele(rs_ids)


    def make_disease_pd(self):
        diseases_pd = pd.DataFrame(columns=['Disease','#Reference','#SNVs','refcount'])
        rownum = 0
        for disease in self.Graph.get_objects(['Disease']):
            disease_name = disease['Name'][0]
            disease_neighborhood = self.Graph.get_neighbors_graph(set(disease['Id']))
            references = list(disease_neighborhood.load_references())
            references.sort(lambda x: int(x[PUBYEAR][0]), reverse=True)
            # references can be from GeneticChange or FunctionalAssociation
            refcount = len(references)
            refs_with_pmids = [ref for ref in references if 'PMID' in ref.Identifiers.keys()]
            if refs_with_pmids:
                pubmed_link = pubmed_hyperlink([x.Identifiers['PMID'] for x in refs_with_pmids],refcount)
            else:
                pubmed_link = ','.join([x._identifiers_str() for x in references])

            GVs = disease_neighborhood.get_objects(['GeneticVariant'])
            snp_count = len(GVs)
            snp_names = disease_neighborhood.get_properties(self.gv_ids,'Name')
            rs_ids = [v[0] for v in snp_names.values()]
            snp_links = dbsnp_hyperlink(rs_ids) if rs_ids else '0'

            diseases_pd.loc[rownum] = [disease_name,pubmed_link,snp_links,refcount]
            rownum += 1

        diseases_pd.sort_values(by=['refcount'],ascending=False, inplace=True)
        diseases_pd.drop(columns=['refcount'],inplace=True)
        return diseases_pd


    def make_snv_pd(self):
        gv_pd = pd.DataFrame(columns=['SNV','MAF','Disease','#Reference','refcount'])
        rownum = 0
        for gv in self.gv2disease_graph.get_objects(['GeneticVariant']):
            gv_id = gv['Id'][0]
            gv_diseases = self.gv2disease_graph.get_neighbors_graph(set(gv['Id']))
            snp_name = gv['Name'][0]
            snp_link = dbsnp_hyperlink([snp_name],as_count=False)
            try:
                maf = self.snp2allele2freq[snp_name]
            except KeyError:
                maf = 'Minor allele frequency (MAF) is not found in dbSNP'
            for regulatorID, targetID, rel in gv_diseases.edges.data('relation'):
                disease_id = regulatorID if regulatorID != gv_id else targetID
                disease_obj = gv_diseases._get_node(disease_id)

                disease_name = disease_obj['Name'][0]

                refcnt = rel.get_reference_count()
                references = rel._get_refs()
                pmids = [x.Identifiers['PMID'] for x in references if 'PMID' in x.Identifiers.keys()]
                if pmids:
                    refcount_links = pubmed_hyperlink(pmids,refcnt)
                else:
                    refcount_links = ','.join([x._identifiers_str() for x in references])

                gv_pd.loc[rownum] = [snp_link,maf,disease_name,refcount_links,refcnt]
                rownum += 1

        gv_pd.sort_values(by=['refcount'],ascending=False, inplace=True)
        gv_pd.drop(columns=['refcount'],inplace=True)
        return gv_pd



if __name__ == "__main__":

    targets = ['CNR1','CNR2','GPR18','GPR55','GPR119']
    xslx_file = DATA_DIR+','.join(targets)+'_genetics.xlsx'
    writer = pd.ExcelWriter(xslx_file, engine='xlsxwriter')

    for target_name in targets:
        ps_api = TargetGenetics(load_api_config())
        ps_api.load_graph(target_name)
        diseases_pd = ps_api.make_disease_pd()
        gv_pd = ps_api.make_snv_pd()
        
        diseases_pd.to_excel(writer, sheet_name=target_name+'_diseases', index=False)
        gv_pd.to_excel(writer, sheet_name=target_name+'_SNVs', index=False)
        
    writer.save()




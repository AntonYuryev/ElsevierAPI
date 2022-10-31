from ElsevierAPI import load_api_config
from ElsevierAPI.ResnetAPI.SemanticSearch import COUNTS, SemanticSearch,df,BIBLIO_PROPERTIES
from ElsevierAPI.ResnetAPI.ResnetAPISession import SNIPPET_PROPERTIES
from ElsevierAPI.ResnetAPI.PathwayStudioGOQL import OQL
from ElsevierAPI.ResnetAPI.NetworkxObjects import PSObject,REFCOUNT
from ElsevierAPI.ETM_API.references import JOURNAL_PROPS,JOURNAL
from  ElsevierAPI.ResnetAPI.Resnet2rdf import ResnetGraph, ResnetRDF
import numpy as np
import time
#from ElsevierAPI.ETM_API.etm import ETMstat


QUANTITATIVE_BIOMARKER_RELS = ['Biomarker','StateChange','QuantitativeChange','CellExpression']
GENETIC_BIOMARKER_RELS = ['GeneticChange']
DISEASE_TYPES = ['Disease','ClinicalParameter']
QUANTITATIVE_BIOMARKERS = ['Protein','Complex','FunctionalClass','Cell']
GENETIC_BIOMARKERS = ['Protein']
SPECIFICITY = 'Promiscuity'
SOLUBLE_BIOMARKERS_TISSUES = ['blood','plasma','serum']

#Biomarker Types:
QUANTITATIVE = 0
GENETIC = 1
SOLUBLE = 2

#pandas names
WEIGHTED = 'weighted'
SNIPPETS = 'snippets'
RANKED_COUNTS = 'ranked_counts'

class BiomarkerReport(SemanticSearch):
    pass
    journal_filter = dict()

    @staticmethod
    def read_journal_list(fname:str, default_config=True):
        journal_filter = dict()
        header = [JOURNAL,JOURNAL,'ISSN','ESSN'] #default configuration
        with open(fname,'r') as f:
            line = f.readline().strip()
            if not default_config:
                header = line.split('\t') #header must contain one of: Journal,ISSN,ESSN
            line = f.readline().strip()
            while line:
                row = line.split('\t')
                for i in range (0, len(header)):
                    col_name = header[i]
                    try:
                        journal_filter[col_name].update([row[i]])
                    except KeyError: journal_filter[col_name] = set([row[i]])
                line = f.readline().strip()

        return journal_filter
    

    def __init__(self, APIconfig,params:dict):
        what2retrieve = SNIPPET_PROPERTIES if params['print_references'] else BIBLIO_PROPERTIES
        super().__init__(APIconfig,what2retrieve)
        self.add_params(params)
        if self.params['biomarker_type'] == SOLUBLE: 
            self.add_rel_props(['Tissue'])
        if self.params['journal_filter_fname']:
            self.journal_filter = self.read_journal_list(self.params['journal_filter_fname'])
            self.journal_filter = {k:list(v) for k,v in self.journal_filter.items()}
            self.add_rel_props(list(JOURNAL_PROPS))
            self.columns2drop += [self.__resnet_name__, self.__mapped_by__]


    def find_diseases(self, disease_name:str):
        mapping_prop = 'Name'
        prop2obj = self.map_prop2entities([disease_name],mapping_prop,DISEASE_TYPES,get_childs=True)
        if not prop2obj:
            mapping_prop = 'Alias'
            prop2obj = self.map_prop2entities([disease_name],mapping_prop,DISEASE_TYPES,get_childs=True)
        
        if prop2obj: # prop2obj = {propValue:{id:psObj}}
            mapped_diseases = list(prop2obj[disease_name].values())
            self.Disease = PSObject(mapped_diseases[0])
        else:
            print('No diseases were found for %s' % disease_name)
            return

        disease_ids = self.Disease[self.__child_ids__] + self.Disease['Id']
        oql_query = 'SELECT Entity WHERE id = ({ids})'.format(ids=','.join(map(str,disease_ids)))
        request_name = 'Find children for {disease}'.format(disease=self.Disease['Name'][0])
        disease_graph = self.process_oql(oql_query, request_name)
        self.diseases = disease_graph.get_objects(SearchValues=DISEASE_TYPES)


    def load_graph(self, disease_ids:list):
        if self.params['biomarker_type'] == GENETIC:
            biomarker_types = GENETIC_BIOMARKERS 
            biomarker_rel_types = GENETIC_BIOMARKER_RELS
        else:
            biomarker_types = QUANTITATIVE_BIOMARKERS
            biomarker_rel_types = QUANTITATIVE_BIOMARKER_RELS
        
        if self.params['biomarker_type'] ==  SOLUBLE:
            oql_query = 'Select Relation WHERE NeighborOf (SELECT Entity WHERE id = (' + ','.join(map(str,disease_ids)) + '))'
            oql_query += " AND objectType = (" + ','.join(biomarker_rel_types) + ')'
            oql_query += ' AND NeighborOf (SELECT Entity WHERE objectType = (' + ','.join(biomarker_types) + ') AND "Cell Localization" = Secreted)'
            oql_query += ' AND Tissue = ({blood_tissues})'.format(blood_tissues=','.join(SOLUBLE_BIOMARKERS_TISSUES))
            query_name = 'Find soluble {types} biomarkers'.format(types=','.join(biomarker_types))
        else:
            oql_query = OQL.expand_entity_by_id(disease_ids,biomarker_rel_types,biomarker_types)
            query_name = 'Find {types} biomarkers'.format(types=','.join(biomarker_types))

        self.process_oql(oql_query, query_name)

        if self.params['biomarker_type'] != GENETIC:
            select_metabolites = '(SELECT Entity WHERE objectType=SmallMol AND Class=\'Endogenous compound\')'
            oql_query = 'SELECT Relation WHERE NeighborOf (SELECT Entity WHERE id = ({ids})) AND NeighborOf {metabolites} AND objectType = ({biomarkers})'
            query_name = 'Find metabolite biomakers'
            if self.params['biomarker_type'] == SOLUBLE:
                oql_query = oql_query + ' AND Tissue = ({blood_tissues})'.format(blood_tissues=','.join(SOLUBLE_BIOMARKERS_TISSUES))
                query_name = 'Find metabolite biomakers in blood tissues'
            oql_query = oql_query.format(ids=','.join(map(str,disease_ids)),metabolites=select_metabolites, biomarkers=','.join(biomarker_rel_types))
            self.process_oql(oql_query, query_name)

        biomarker_ids = [y['Id'][0] for x,y in self.Graph.nodes(data=True) if y['ObjTypeName'][0] not in DISEASE_TYPES]
        disease_ids =  [y['Id'][0] for x,y in self.Graph.nodes(data=True) if y['ObjTypeName'][0] in DISEASE_TYPES]
        oql_query = OQL.connect_ids(biomarker_ids,disease_ids,['Regulation'])
        self.process_oql(oql_query, 'Find Regulation between biomarkers and diseases')

        self.Graph.load_references()

        if self.params['biomarker_type'] == SOLUBLE:
            self.Graph = self.Graph.filter_references({'Tissue':SOLUBLE_BIOMARKERS_TISSUES},QUANTITATIVE_BIOMARKER_RELS)

        if self.journal_filter:
            self.Graph = self.Graph.filter_references(self.journal_filter)
            # journal_filter = {prop_name:[values]}


    def init_semantic_search (self, reset_pandas = False):
        BiomarkerNames =  [y['Name'][0] for x,y in self.Graph.nodes(data=True) if y['ObjTypeName'][0] not in DISEASE_TYPES]
        BiomarkerScores = df()
        BiomarkerScores['Name'] = np.array(BiomarkerNames)
        print('Will score %d biomarkers linked to %s' % (len(BiomarkerScores),self.Disease['Name'][0]))
        self.RefCountPandas = self.load_pandas(BiomarkerScores,prop_names_in_header=True) #maps entities in ['Name'] column


    def __input_disease_column(self):
        return self._col_name_prefix+self.Disease['Name'][0]


    def semantic_search(self):
        disease_name2ids = dict()
        for disease in self.diseases:
            dis_name = disease['Name'][0]
            disease_ids = list(self._get_obj_ids_by_props(disease['Id'],["Id"],only_obj_types=DISEASE_TYPES)) #finds children for disease
            disease_name2ids[dis_name] = disease_ids

        if self.params['biomarker_type'] == GENETIC:
            biomarker_rel_types = QUANTITATIVE_BIOMARKER_RELS + GENETIC_BIOMARKER_RELS
        else:
            biomarker_rel_types = QUANTITATIVE_BIOMARKER_RELS + ['Regulation']

        for disease in self.diseases:
            dis_name = disease['Name'][0]
            self.set_how2connect(biomarker_rel_types,[],'')
            linked_entities_count = self.link2RefCountPandas(dis_name,disease_name2ids[dis_name])
            print('%d biomarkers are linked to %s' % (linked_entities_count,dis_name))

        #adding biomarker type column
        self.name2objs, objid2names = self.Graph.get_prop2obj_dic('Name', list(bm.RefCountPandas['Name']))
        self.RefCountPandas['Type'] = self.RefCountPandas['Name'].apply(lambda x: self.name2objs[x][0]['ObjTypeName'][0])

        return self.make_count_df()
        

    def biomarker_specificity(self,biomarker_ids:tuple, disease_ids=[]):
        #biomarker specificity can be calculated only for quantitative biomarkers
        rel_type_str = ','.join(QUANTITATIVE_BIOMARKER_RELS)
        biomarker_ids_str = ','.join(map(str,list(biomarker_ids)))
        if disease_ids:
            oql_query = 'SELECT Entity WHERE id = ({diseases}) AND Connected by (SELECT Relation WHERE objectType = ({rel_types})) to (SELECT Entity WHERE id = ({biomarkers}))'
            disease_ids_str = ','.join(map(str,list(disease_ids)))
            oql_query = oql_query.format(diseases=disease_ids_str, rel_types=rel_type_str, biomarkers=biomarker_ids_str)
        else:
            oql_query = 'SELECT Entity WHERE Connected by (SELECT Relation WHERE objectType = ({rel_types})) to (SELECT Entity WHERE id = ({biomarkers}))'
            oql_query = oql_query.format(rel_types=rel_type_str, id=biomarker_ids_str)
        
        oql_query = oql_query.format(diseases=disease_ids_str, rel_types=rel_type_str, biomarkers=biomarker_ids_str)
        disease_count = self.get_result_size(oql_query)
        return disease_count


    def add_specificity(self, to_df:df):
        biomarker_count = len(to_df)
        message = f'Finding biomarker specificity for {str(biomarker_count)} biomarkers'
        disease_ids = list()
        self.specificity_name = SPECIFICITY
        diseases4specificity = ','.join(self.params['add_specificity4diseases'])
        message = message + f' for {diseases4specificity}'
        limit2disease_graph = self.child_graph(self.params['add_specificity4diseases'],['Name','Alias'])
        disease_ids = list(limit2disease_graph)
        self.specificity_name += ' for '+diseases4specificity

        print('Finding biomarker specificity for %d biomarkers' % biomarker_count)

        start_specificity = time.time()
        to_df[self.specificity_name] = to_df[self.__temp_id_col__].apply(lambda x: self.biomarker_specificity(x,disease_ids))
        print ('Specificity for %d biomarkers was found in %s' %(biomarker_count, self.execution_time(start_specificity)))

 
    def weighted_counts(self):
        print('Begin weighted refcount calculation')
        self.weight_prop = 'ObjTypeName'
        self.weight_dict = {'Biomarker':1.0,'StateChange':0.5,'QuantitativeChange':0.5,'GeneticChange':0.5,'CellExpression':0.5,'Regulation':0.25}
        self.__print_refs__ = False
        refcount_cols = [col for col in self.RefCountPandas.columns if self._col_name_prefix in col]
        for col in refcount_cols:
            self.RefCountPandas[col] = self.RefCountPandas[col].astype(float)
        
        return self.semantic_search()

    def make_report(self):
        counts_df = self.raw_data['counts']
        # decide what to add to the report
        try:
            weighted_df = df(self.raw_data[WEIGHTED])
            weighted_df.drop(columns=[self.__mapped_by__,'Resnet name'],inplace=True)
            input_disease_col = self.__input_disease_column()
            weighted_df['rank'] = weighted_df[input_disease_col]
            weighted_df.sort_values(by=['rank'],ascending=False, inplace=True)
            weighted_df.insert(2, 'rank', weighted_df.pop('rank'))
            weighted_df.insert(3, input_disease_col, weighted_df.pop(input_disease_col))
            report_df = self._merge_counts2norm(counts_df,weighted_df)
            report_df._name_ = 'ranked_counts'
            print ('Created %s table' % report_df._name_)
        except KeyError:
            report_df = counts_df

        if self.params['add_specificity4diseases']:
            #self.disease_ontology_categoty = add_specificity4
            self.add_specificity()
            normalized_score_col = 'Normalized score'
            report_df[normalized_score_col] = report_df[input_disease_col]/report_df[self.specificity_name]
            report_df.sort_values(by=[normalized_score_col],ascending=False, inplace=True)
            report_df.insert(2, normalized_score_col, report_df.pop(normalized_score_col))
            report_df.insert(3, input_disease_col, report_df.pop(input_disease_col))
            report_df.insert(4, self.specificity_name, report_df.pop(self.specificity_name)) 
                       
        self.add2report(report_df)
        self.add_ps_bibliography()

        if self.params['print_references']:
            ref_df = self.Graph.snippets2df()
            self.add2report(ref_df)

        
    def biomarker_disease_scores(self, max_etm_row=100)->ResnetGraph:
        biomarker2disease = df()
        row_counter = 0
        ref_count_columns = [col for col in self.RefCountPandas.columns if self._col_name_prefix in col]
        for row in self.RefCountPandas.index:
            biomarker = self.RefCountPandas.loc[row][self.RefCountPandas.columns[0]]
            for col_name in ref_count_columns:
                disease = col_name[len(self._col_name_prefix):]
                refcount = self.RefCountPandas.loc[row][col_name]
                # 'Rank is also available from self.report_pandas[0]['Rank']
                if refcount > 0.0000001:
                    biomarker2disease.at[row_counter,'Biomarker'] = biomarker
                    biomarker2disease.at[row_counter,'Disease'] = disease
                    biomarker2disease.at[row_counter,'Rank'] = float(refcount)

                    row_counter += 1

        biomarker2disease.sort_values(by=['Rank'],ascending=False, inplace=True,ignore_index=True)
        rel_props = {'ObjTypeName':['Biomarker']}
        if self.params['biomarker_type'] == GENETIC:
            add2etm_query = ['terms for genetic variations']
            rel_props.update({'BiomarkerType':['Genetic']})
        elif self.params['biomarker_type'] == SOLUBLE:
            rel_props.update({'BiomarkerType':['Blood']})
            add2etm_query = ['blood']

        #etm_counter = ETMstat(self.APIconfig)
        bm2dis_graph = ResnetGraph() # to convert into RDF
        if self.journal_filter:
            biomarker2disease = self.Graph.add_recent_refs(biomarker2disease,'Biomarker','Disease')
            #biomarker2disease.dropna(subset=['Recent PMIDs','Recent DOIs'],inplace=True, how='all')
        elif max_etm_row:
            biomarker2disease=self.etm_counter.add_etm42columns(biomarker2disease,'Biomarker','Disease',add2etm_query)
            self.etm2df()

            """
            for i in range(0,max_etm_row): 
                biomarker = biomarker2disease.loc[i]['Biomarker']
                disease = biomarker2disease.loc[i]['Disease']
                biomarker_objs = self.Graph.get_obj_by_prop(biomarker)
                disease_objs = self.Graph.get_obj_by_prop(disease)
                rel_props[REFCOUNT] = [int(etm_hit_count)]
                [bm2dis_graph.add_triple(d,b,rel_props,etm_refs) for b in biomarker_objs for d in disease_objs] 
            """

            '''
            for i in range(0,max_etm_row):
                biomarker = biomarker2disease.loc[i]['Biomarker']
                disease = biomarker2disease.loc[i]['Disease']
                search_terms = [biomarker,disease] + add2etm_query
                etm_hit_count,ref_ids,etm_refs = self.etm_counter.relevant_articles(search_terms)
                if ref_ids:
                    biomarker2disease.at[i,'#References'] = int(etm_hit_count)
                    try:
                        biomarker2disease.at[i,'Top PMIDs'] = ','.join(ref_ids['PMID'])
                    except KeyError:
                        biomarker2disease.at[i,'Top PMIDs'] = ''
                    try:
                        biomarker2disease.at[i,'Top DOIs'] = ','.join(ref_ids['DOI'])
                    except KeyError:
                        biomarker2disease.at[i,'Top DOIs'] = ''

                ''' 
            
        biomarker2disease['Type'] = biomarker2disease['Biomarker'].apply(lambda x: self.name2objs[x][0]['ObjTypeName'][0])
        biomarker2disease._name_ = 'biomarker-disease'
        self.add2report(biomarker2disease)
        return bm2dis_graph


    def report_prefix(self):
        biomarker_descr = ' quantitative' 
        if self.params['biomarker_type'] == GENETIC: biomarker_descr = ' genetic'
        if self.params['biomarker_type'] == SOLUBLE: biomarker_descr = ' soluble'

        f_prefix = self.Disease['Name'][0]+biomarker_descr
        return f_prefix+' from selected journals' if self.journal_filter else f_prefix


    def print_report(self):
        report_fname = self.params['data_dir']+self.report_prefix()+' biomarkers.xlsx'
        super().print_report(report_fname)

    def print_rawdata(self):
        raw_datafile = self.params['data_dir']+self.report_prefix()+' biomarkers raw data.xlsx'
        super().print_rawdata(raw_datafile)



DISEASE_NAME = 'Myocardial Infarction'#'atrial fibrillation'#'Pulmonary Hypertension'#'diabetes mellitus'   #'fibrosis'
DATA_DIR = 'D:/Python/PMI/'+DISEASE_NAME+'/'
#DATA_DIR = 'D:/Python/PMI/' 

if __name__ == "__main__":
    params = { 
                'disease': DISEASE_NAME,
                'print_references':True,
                'biomarker_type': SOLUBLE, #GENETIC, #
                'journal_filter_fname':'', #'D:/Python/Quest/raw/'+DISEASE_NAME+'/'+DISEASE_NAME+' high-quality journals.txt'
                'add_specificity4diseases':[], #['respiratory disease']
                # there is no scientific rational to calculate biomarker specificity for genetic biomarkers
                'data_dir' : 'D:/Python/PMI/'+DISEASE_NAME+'/'
                }

    start_time = time.time()
    api_config = str()
    #api_config = 'D:/Python/ENTELLECT_API/ElsevierAPI/APIconfigTeva.json'
    APIconfig = load_api_config(api_config)
    
    bm = BiomarkerReport(APIconfig,params)
    bm.flush_dump()

    bm.find_diseases(DISEASE_NAME) #'Uterine neoplasm'
    
    disease_ids = [y['Id'][0] for x,y in bm.Graph.nodes(data=True) if y['ObjTypeName'][0] in DISEASE_TYPES]
    bm.load_graph(disease_ids)
    bm.init_semantic_search()
    counts_df = bm.semantic_search()
    counts_df._name_ = COUNTS
    bm.add2raw(counts_df)
    
    bm.drop_refcount_columns()
    weighted_df = bm.weighted_counts()
    weighted_df._name_= WEIGHTED
    bm.add2raw(weighted_df)

    bm.make_report()
    bm.biomarker_disease_scores()
    #bm2dis_graph = bm.biomarker_disease_scores()
    #jsonld_fname = bm.params['data_dir']+bm.data_dir+bm.report_prefix()+' biomarker-disease graph.jsonld'
    #ResnetRDF.fromResnetGraph(bm2dis_graph).to_json(jsonld_fname)

    bm.print_report()
    bm.print_rawdata()
    print('Biomarkers for %s was found in %s' % (bm.Disease['Name'][0], bm.execution_time(start_time)))

    

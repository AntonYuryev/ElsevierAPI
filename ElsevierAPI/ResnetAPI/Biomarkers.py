from .SemanticSearch import SemanticSearch,df,COUNTS,RANK
from ..ETM_API.etm import REFCOUNT_COLUMN, ETMstat
from .ResnetAPISession import SNIPPET_PROPERTIES,BIBLIO_PROPERTIES
from .PathwayStudioGOQL import OQL
from .NetworkxObjects import PSObject,REFCOUNT,CHILDS
from ..ETM_API.references import JOURNAL_PROPS,JOURNAL
from  .Resnet2rdf import ResnetGraph,ResnetRDF
import time

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
        params['what2retrieve'] = SNIPPET_PROPERTIES if params['print_snippets'] else BIBLIO_PROPERTIES
        super().__init__(APIconfig,**params)
        if self.params['biomarker_type'] == SOLUBLE: 
            self.add_rel_props(['Tissue'])
        if self.params['journal_filter_fname']:
            self.journal_filter = self.read_journal_list(self.params['journal_filter_fname'])
            self.journal_filter = {k:list(v) for k,v in self.journal_filter.items()}
            self.add_rel_props(list(JOURNAL_PROPS))
            self.columns2drop += [self.__resnet_name__, self.__mapped_by__]

        self.columns2drop += [self.__mapped_by__, self.__resnet_name__]


    def find_diseases(self, disease_name:str):
        my_diseases = self._props2psobj([disease_name],['Name'])
        if not my_diseases:
             my_diseases = self._props2psobj([disease_name],['Alias'])

        if my_diseases:
            self.Disease = my_diseases[0]
        else:
            print('No diseases were found for %s' % disease_name)
            return

        self.diseases = self.Disease[CHILDS] + [self.Disease]


    def load_graph(self, disease_ids:list):
        if self.params['biomarker_type'] == GENETIC:
            biomarker_types = GENETIC_BIOMARKERS 
            self.biomarker_rel_types = GENETIC_BIOMARKER_RELS
        else:
            biomarker_types = QUANTITATIVE_BIOMARKERS
            self.biomarker_rel_types = QUANTITATIVE_BIOMARKER_RELS
        
        if self.params['biomarker_type'] ==  SOLUBLE:
            oql_query = 'Select Relation WHERE NeighborOf (SELECT Entity WHERE id = (' + ','.join(map(str,disease_ids)) + '))'
            oql_query += " AND objectType = (" + ','.join(self.biomarker_rel_types) + ')'
            oql_query += ' AND NeighborOf (SELECT Entity WHERE objectType = (' + ','.join(biomarker_types) + ') AND "Cell Localization" = Secreted)'
            oql_query += f' AND Tissue = ({",".join(SOLUBLE_BIOMARKERS_TISSUES)})'
            query_name = f'Find soluble {",".join(biomarker_types)} biomarkers'
        else:
            oql_query = OQL.expand_entity_by_id(disease_ids,self.biomarker_rel_types,biomarker_types)
            query_name = 'Find {types} biomarkers'.format(types=','.join(biomarker_types))

        self.process_oql(oql_query, query_name)

        if self.params['biomarker_type'] != GENETIC:
            select_metabolites = '(SELECT Entity WHERE objectType=SmallMol AND Class=\'Endogenous compound\')'
            oql_query = 'SELECT Relation WHERE NeighborOf (SELECT Entity WHERE id = ({ids})) AND NeighborOf {metabolites} AND objectType = ({biomarkers})'
            query_name = 'Find metabolite biomakers'
            if self.params['biomarker_type'] == SOLUBLE:
                oql_query = oql_query + ' AND Tissue = ({blood_tissues})'.format(blood_tissues=','.join(SOLUBLE_BIOMARKERS_TISSUES))
                query_name = 'Find metabolite biomakers in blood tissues'
            oql_query = oql_query.format(ids=','.join(map(str,disease_ids)),metabolites=select_metabolites, biomarkers=','.join(self.biomarker_rel_types))
            self.process_oql(oql_query, query_name)
        else:
            oql_query = OQL.expand_entity(disease_ids,['Id'], expand2neighbors=['GeneticVariant'])
            disease2gvs = self.process_oql(oql_query, 'Find GVs linked to diseases')
            gv_ids = disease2gvs.dbids4nodes(['GeneticVariant'])
            gvid2genes = self.gv2gene(gv_ids)
            genes_with_gvs = {name for names in [names for names in gvid2genes.values()] for name in names}
            print('Found %d GeneticVariants linked to input disease in %d genes' % (len(gvid2genes),len(genes_with_gvs)))

        biomarker_ids = list()
        disease_ids = list()
        graph_nodes = self.Graph._get_nodes()
        [(biomarker_ids,disease_ids)[o.objtype() in DISEASE_TYPES].append(o.dbid()) for o in graph_nodes]

        oql_query = OQL.connect_ids(biomarker_ids,disease_ids,['Regulation'])
        self.process_oql(oql_query, 'Find Regulation between biomarkers and diseases')

        self.Graph.load_references()

        if self.params['biomarker_type'] == SOLUBLE:
            self.Graph = self.Graph.filter_references({'Tissue':SOLUBLE_BIOMARKERS_TISSUES},QUANTITATIVE_BIOMARKER_RELS)

        if self.journal_filter:
            self.Graph = self.Graph.filter_references(self.journal_filter)
            # journal_filter = {prop_name:[values]}


    def init_semantic_search (self, reset_pandas = False):
        not_biomarker_types = DISEASE_TYPES+['GeneticVariant']
        biomarkers = [PSObject(y) for x,y in self.Graph.nodes(data=True) if y['ObjTypeName'][0] not in not_biomarker_types]
        self.RefCountPandas = self.load_df(biomarkers,max_child_count=10)


    def __input_disease_column(self):
        return self._col_name_prefix+self.Disease.name()


    def semantic_search(self):
        for disease in self.diseases:
            self.load_children4([disease])

        if self.params['biomarker_type'] == GENETIC:
            biomarker_rel_types = QUANTITATIVE_BIOMARKER_RELS + GENETIC_BIOMARKER_RELS
        else:
            biomarker_rel_types = QUANTITATIVE_BIOMARKER_RELS + ['Regulation']

        disease_count = 0
        for disease in self.diseases:
            dis_name = disease['Name'][0]
            disease_and_childs = disease.childs()+[disease]
            self.set_how2connect(biomarker_rel_types,[],'')
            linked_entities_count = self.link2RefCountPandas(dis_name,disease_and_childs)
            disease_count += 1
            print('%d biomarkers are linked to %s' % (linked_entities_count,dis_name))
            print('%d concepts out of %d were linked' % (disease_count,len(self.diseases)))

        #adding biomarker type column:
        self.name2objs, objid2names = self.Graph.get_prop2obj_dic('Name', list(self.RefCountPandas['Name']))
        self.RefCountPandas['Type'] = self.RefCountPandas['Name'].apply(lambda x: self.name2objs[x][0]['ObjTypeName'][0])

        return self.make_count_df()
        

    def biomarker_specificity(self,biomarker_ids:tuple, disease_ids=[]):
        '''
            # biomarker specificity can be calculated only for quantitative biomarkers
        '''
        rel_type_str = ','.join(self.biomarker_rel_types)
        biomarker_ids_str = ','.join(map(str,list(biomarker_ids)))
        if disease_ids:
            disease_ids_str = ','.join(map(str,list(disease_ids)))
            oql_query = f'SELECT Entity WHERE id = ({disease_ids_str}) AND Connected by (SELECT Relation WHERE \
                objectType = ({rel_type_str})) to (SELECT Entity WHERE id = ({biomarker_ids_str}))'
        else:
            oql_query = f'SELECT Entity WHERE Connected by (SELECT Relation WHERE objectType = ({rel_type_str})) \
                to (SELECT Entity WHERE id = ({biomarker_ids_str}))'
        
        disease_count = self.get_result_size(oql_query)
        return disease_count


    def add_specificity(self, to_df:df):
        if self.params['biomarker_type'] == GENETIC:
            print('!!!!Biomarker specificity cannot be calculated for genetic biomarkers!!!!')
            return
        biomarker_count = len(to_df)
        message = f'Finding biomarker specificity for {str(biomarker_count)} biomarkers'
        #disease_ids = list()
        diseases4specificity = ','.join(self.params['add_specificity4diseases'])
        message = message + f' for {diseases4specificity}'
        limit2disease_graph = self.child_graph(self.params['add_specificity4diseases'],['Name','Alias'])
        disease_dbids = list(limit2disease_graph.dbid2uid().keys())
        specificity_colname = SPECIFICITY+' for '+ diseases4specificity

        print('Finding biomarker specificity for %d biomarkers' % biomarker_count)

        start_specificity = time.time()
        to_df[specificity_colname] = to_df[self.__temp_id_col__].apply(lambda x: self.biomarker_specificity(x,disease_dbids))
        print ('Specificity for %d biomarkers was found in %s' %(biomarker_count, self.execution_time(start_specificity)))
        return specificity_colname
    
 
    def weighted_counts(self):
        self.drop_refcount_columns()
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
        input_disease_col = self.__input_disease_column()
        # decide what to add to the report
        try:
            weighted_df = self.raw_data[WEIGHTED]
            #weighted_df.drop(columns=[self.__mapped_by__,'Resnet name'],inplace=True)
            weighted_df[RANK] = weighted_df[input_disease_col]
            weighted_df.sort_values(by=[RANK],ascending=False, inplace=True)
            weighted_df.insert(2, RANK, weighted_df.pop(RANK))
            weighted_df.insert(3, input_disease_col, weighted_df.pop(input_disease_col))
            report_df = self._merge_counts2norm(counts_df,weighted_df)
        except KeyError:
            report_df = counts_df

        if self.params['add_specificity4diseases']:
            if self.params['biomarker_type'] != GENETIC:
                specificity_colname = self.add_specificity(report_df)
                normalized_score_col = input_disease_col+' normalized by '+specificity_colname.lower()
                report_df[normalized_score_col] = report_df[input_disease_col]/report_df[specificity_colname]
                report_df.sort_values(by=[normalized_score_col],ascending=False, inplace=True)
                report_df.insert(2, normalized_score_col, report_df.pop(normalized_score_col))
                report_df.insert(3, input_disease_col, report_df.pop(input_disease_col))
                report_df.insert(4, specificity_colname, report_df.pop(specificity_colname)) 
                report_df.add_column_format(normalized_score_col,'align','center')
            else:
                print('!!!!Biomarker specificity cannot be calculated for genetic biomarkers!!!!')
        
        report_df._name_ = RANKED_COUNTS
        self.add2report(report_df)
        print ('Created %s table' % report_df._name_)
        self.add_graph_bibliography()

        if self.params['print_snippets']:
            ref_df = self.Graph.snippets2df()
            self.add2report(ref_df)

        
    def biomarker_disease_scores(self, max_etm_row=100)->ResnetGraph:
        '''
        Adds
        ----
        biomarker2disease to report\n
        self.bm2dis_graph if self.params['print_rdf'] = True
        '''
        biomarker2disease = df()
        row_counter = 0
        ref_count_columns = [col for col in self.RefCountPandas.columns if self._col_name_prefix in col]
        for row in self.RefCountPandas.index:
            biomarker = self.RefCountPandas.loc[row][self.RefCountPandas.columns[0]]
            for col_name in ref_count_columns:
                disease = col_name[len(self._col_name_prefix):]
                refcount = self.RefCountPandas.loc[row][col_name]
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
        else:
            add2etm_query = list()

        #etm_counter = ETMstat(self.APIconfig)
        
        if self.journal_filter:
            biomarker2disease = self.Graph.add_recent_refs(biomarker2disease,'Biomarker','Disease')
            #biomarker2disease.dropna(subset=['Recent PMIDs','Recent DOIs'],inplace=True, how='all')
        elif max_etm_row:
            biomarker2disease=self.etm_counter.refs42columns(biomarker2disease,'Biomarker','Disease',ETMstat.basic_query,add2etm_query)
            self.add_bibliography()

        if self.params['print_rdf']:
            self.bm2dis_graph = ResnetGraph() # to convert into RDF
            for i in range(0,max_etm_row): 
                biomarker_name = biomarker2disease.loc[i]['Biomarker']
                disease_name = biomarker2disease.loc[i]['Disease']
                biomarker_objs = self.Graph._psobjs_with(biomarker_name)
                disease_objs = self.Graph._psobjs_with(disease_name)
                rel_props[REFCOUNT] = biomarker2disease.loc[i][REFCOUNT_COLUMN]
                [self.bm2dis_graph.add_triple(d,b,rel_props) for b in biomarker_objs for d in disease_objs] 
            
        biomarker2disease['Type'] = biomarker2disease['Biomarker'].apply(lambda x: self.name2objs[x][0]['ObjTypeName'][0])
        biomarker2disease = df(biomarker2disease.sort_values(RANK,ascending=False))
        biomarker2disease._name_ = 'biomarker-disease'
        biomarker2disease.add_column_format(RANK,'align','center')
        biomarker2disease.add_column_format(REFCOUNT_COLUMN,'align','center')
        biomarker2disease.add_column_format('Disease','width',60)
        biomarker2disease.set_hyperlink_color([REFCOUNT_COLUMN])
        self.add2report(biomarker2disease)


    def report_name(self,extension=''):
        biomarker_descr = ' quantitative' 
        if self.params['biomarker_type'] == GENETIC: biomarker_descr = ' genetic'
        if self.params['biomarker_type'] == SOLUBLE: biomarker_descr = ' soluble'

        f_prefix = self.params['disease']+biomarker_descr+' biomarkers'
        return f_prefix+' from selected journals'+extension if self.journal_filter else f_prefix+extension

    
    def report_path(self,extension=''):
        return self.data_dir+self.report_name(extension)


    def print_report(self):
        self.find_diseases(self.params['disease'])
        disease_ids = self.Graph.dbids4nodes(DISEASE_TYPES)
        self.load_graph(disease_ids)
        self.init_semantic_search()
        counts_df = self.semantic_search()
        counts_df._name_ = COUNTS
        self.add2raw(counts_df)
        
        weighted_df = self.weighted_counts()
        weighted_df._name_= WEIGHTED
        self.add2raw(weighted_df)

        self.make_report()
        self.biomarker_disease_scores()

        if self.params['print_rdf']:
            jsonld_fname = self.report_path(' biomarker-disease graph.jsonld')
            ResnetRDF.fromResnetGraph(self.bm2dis_graph).to_json(jsonld_fname)

        super().print_rawdata(self.report_name(' raw data.xlsx'))
        super().print_report(self.report_name('.xlsx'))
        
from .ResnetGraph import PSObject
from .SemanticSearch import SemanticSearch,df,ResnetGraph
from .ResnetAPISession import BIBLIO_PROPERTIES
from ElsevierAPI import execution_time
import time

SYMPTOMS_INDICATIONS = 'Indications'
DIAGNOSIS_WS = 'Diagnosis'
SYMPTOMS_COUNT = '# symptoms'

class Diagnosis(SemanticSearch):
    pass
    max_threads4ontology = 4
  
    def __init__(self, *args, **kwargs):
        '''
        Input
        -----
        APIconfig - args[0]
        '''
        my_kwargs = {'symptoms_fname':'',
                     'what2retrieve': BIBLIO_PROPERTIES
                     
                }
        
        my_kwargs.update(kwargs)
        super().__init__(*args,**my_kwargs)
        self.columns2drop += [self.__resnet_name__,self.__mapped_by__]
        self.max_ontology_parent = int(my_kwargs.pop('max_ontology_parent',11))


    def load_symptoms(self)->set[PSObject]:
        with open(self.params['symptoms_fname'],'r',encoding='utf-8') as f:
            symptom_names = f.readlines()
        symptom_names = set(map(lambda x: x.strip(), symptom_names))

        oql_query = 'SELECT Entity WHERE (Name,Alias) = ({props})'

        symptoms_g = self.iterate_oql(oql_query,symptom_names,request_name='Retrieve symptoms from database')
        self.symptoms = set(symptoms_g._get_nodes())
        return self.symptoms
 

    def symptom_names(self):
        """
        Returns
        -------
        entity names that can be used both in PS and ETM
        """
        #input names are used for finding references in ETM.
        # RepurposeDrug has it own 'target_names'.
        return [x['Name'][0] for x in self.symptoms]

        
    def _get_report_name(self):
        fname = self.params['symptoms_fname']
        fname = fname[fname.rfind('/')+1:]
        fname = fname[:fname.rfind('.')]
        return str(self.data_dir+ 'Diagnosis' +'4'+ fname)


    def load_indications(self):
        self.get_symptoms = f'SELECT Entity WHERE id = ({ResnetGraph.dbids(list(self.symptoms))})'
        get_indications = f'SELECT Relation WHERE objectType = FunctionalAssociation AND NeighborOf ({self.get_symptoms})'
        #get_indications = f'SELECT Entity Connected by ({get_rels}) to ({get_symptoms})'
        indication_graph = self.process_oql(get_indications,'find indications')
        assert(isinstance(indication_graph,ResnetGraph))
        # need
        self.indications = indication_graph.regulators(['Disease'],self.params['min_symptoms4indication'])
        self.indications = self.indications.difference(self.symptoms)
        return self.indications


    def init_semantic_search(self):
        '''
        Loads
        -----
        DF4ANTAGONISTS and DF4AGONISTS df to raw_data
        '''
        print('\n\nInitializing semantic search')
        symptoms = self.load_symptoms()
        indications = self.load_indications()

        if symptoms:
            self.RefCountPandas = self.load_df(list(indications),
                                            max_child_count=self.max_ontology_parent,
                                            max_threads=self.max_threads4ontology)

           # self.RefCountPandas = self.expand_entities(self.RefCountPandas,['GeneticVariant'])
            self.RefCountPandas._name_ = SYMPTOMS_INDICATIONS
            return True
        else:
            fname = self.params['symptoms_fname']
            print (f'No indications found for symptoms in {fname}')
            return False
        

    def semscore4symptoms(self):
        colname = 'Symptoms'
        how2connect = self.set_how2connect(['FunctionalAssociation'],[],'',['Regulation'])
        self.link2RefCountPandas(colname,list(self.symptoms),how2connect)
        return
    

    def scoreGVs(self):
        my_session = self._clone_session()
        reltypes = ['FunctionalAssociation','Regulation']
        commonGVs_g = my_session.common_neighbors(['GeneticVariant'],list(self._all_dbids()),reltypes,'',
                                               ResnetGraph.dbids(list(self.symptoms)),reltypes,'')
        gvlinkcounter = 0
        for i in self.RefCountPandas.index:
            target_dbids = list(self.RefCountPandas.at[i,self.__temp_id_col__])
            row_indications = self.Graph.psobj_with_dbids(set(target_dbids))
            indicationGVs = commonGVs_g.get_neighbors(set(row_indications))
                
            GVscore = 0         
            self.__colnameGV__ = 'Genetic Variants for symptoms'
            if indicationGVs:
                GVnames = set(ResnetGraph.names(list(indicationGVs)))
                self.RefCountPandas.at[i,self.__colnameGV__] = ';'.join(GVnames)
                gvlinkcounter += 1
                gv_disease_subgraph = commonGVs_g.neighborhood(indicationGVs)
                GVscore = len(gv_disease_subgraph.load_references())
            
            self.RefCountPandas.at[i,self._col_name_prefix+'Genetic Variants'] = GVscore

        GVcount = len(commonGVs_g.psobjs_with(only_with_values=['GeneticVariant']))
        print(f'Found {gvlinkcounter} indications linked to {GVcount} GVs')


    def add_symptoms(self):
        for i in self.RefCountPandas.index:
            row_dbids = list(self.RefCountPandas.at[i,self.__temp_id_col__])
            row_indications = self.Graph.psobj_with_dbids(set(row_dbids))
            row_symptoms = list(self.Graph.get_neighbors(set(row_indications),list(self.symptoms)))
            row_symptoms_names = ResnetGraph.names(row_symptoms)
            row_symptoms_names.sort()
            self.RefCountPandas.at[i,'Symptoms'] = ','.join(row_symptoms_names)
            self.RefCountPandas.at[i,SYMPTOMS_COUNT] = len(row_symptoms_names)


    def refcount_columns(self,counts_df:df=df(),column_prefix=''):
        to_return = super().refcount_columns(counts_df,column_prefix)
        if SYMPTOMS_COUNT in counts_df.columns.tolist():
            to_return = [SYMPTOMS_COUNT] + to_return
        return to_return


    def make_report(self):
        start_time = time.time()
        self.flush_dump_files()
        self.load_symptoms()
        self.load_indications()
        
        if self.init_semantic_search():
            self.semscore4symptoms()
            self.add_symptoms()
            self.scoreGVs()
            self.add2raw(self.make_count_df(self.RefCountPandas,SYMPTOMS_INDICATIONS))
            self.normalize(SYMPTOMS_INDICATIONS,DIAGNOSIS_WS)
            self.add_etm_refs(DIAGNOSIS_WS,self.symptom_names())
            self.add_ps_bibliography()
            self.add_etm_bibliography()
        print(f'{self._get_report_name()} diagnosis is done in {execution_time(start_time)}')

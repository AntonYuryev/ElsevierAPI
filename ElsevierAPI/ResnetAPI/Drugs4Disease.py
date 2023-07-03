from .DiseaseTargets import DiseaseTargets,ResnetGraph
from .DiseaseTargets import ANTAGONIST_TARGETS_WS,AGONIST_TARGETS_WS,REFERENCE_IDENTIFIERS,UNKNOWN_TARGETS_WS
from .SemanticSearch import RANK
from .ResnetAPISession import APISession,NO_REL_PROPERTIES,BIBLIO_PROPERTIES
from ..pandas.panda_tricks import df,PSObject
from numpy import NaN,nan_to_num
import time
import networkx as nx
from concurrent.futures import ThreadPoolExecutor, as_completed
from .DrugTargetConfidence import DrugTargetConsistency

PATHWAY_REGULATOR_SCORE = 'Disease model regulator score'
DRUG2TARGET_REGULATOR_SCORE = 'Regulator score'
PHARMAPENDIUM_ID = 'Marketed Drugs'


class Drugs4Targets(DiseaseTargets):
    pass
    def __init__(self,*args,**kwargs):
        """
        Input
        -----
        APIconfig - args[0]
        what2retrieve - default BIBLIO_PROPERTIES which is required to load self.set_target_disease_state()
        connect2server - default True
        """
        APIconfig = args[0]
        my_kwargs = {
                'disease':[],
                'what2retrieve':BIBLIO_PROPERTIES,
                'strict_mode':True,
                'data_dir':'',
                'add_bibliography' : False,
                'strict_mode' : False,
                'target_types' : ['Protein'],
                'pathway_folders':[],
                'pathways': [],
                'consistency_correction4target_rank':False
                }
        my_kwargs.update(kwargs)

        super().__init__(APIconfig,**my_kwargs)
        self.clone_all_sessions = my_kwargs.pop('clone',True) # must be different from clone that is used by SemanticSearch

        dt_consist_kwargs = dict(my_kwargs)
        dt_consist_kwargs['what2retrieve'] = REFERENCE_IDENTIFIERS
        self.dt_consist = DrugTargetConsistency(self.APIconfig,**dt_consist_kwargs)
        self.target_uid2rank = dict() # {target_dbid:rank}
        self.column_ranks = dict() # {column rank:column name} rank defines the column position in Drugs df and consequently its weight
        self.direct_target2drugs = PSObject() # used for annotation of ANTAGONIST_TARGETS_WS,AGONIST_TARGETS_WS with drugs
        self.indirect_target2drugs = PSObject() # used for annotation of ANTAGONIST_TARGETS_WS,AGONIST_TARGETS_WS with drugs
        self.add_targets2drugs_ws = True


    @classmethod
    def from_files(*args,**kwargs):
        fake_config = {'ResnetURL':'','PSuserName':'','PSpassword':'','ETMURL':'','ETMapikey':''}
        my_kwargs = dict(kwargs)
        my_kwargs['connect2server'] = False
        dt = Drugs4Targets(fake_config,**my_kwargs)
        dt.clone_all_sessions = False # must be different from clone that is used by SemanticSearch

        dt_consist_kwargs = dict(my_kwargs)
        dt_consist_kwargs['what2retrieve'] = NO_REL_PROPERTIES
        dt.dt_consist = DrugTargetConsistency(fake_config,**dt_consist_kwargs)
        return dt

############################ service functions: read/write, load from database, cache######################################            
    def report_name(self):
        rep_pred = 'suggested' if self.params['strict_mode'] else 'predicted'
        return f'{self._disease2str()} {rep_pred} targets,drugs'


    def report_path(self,extension=''):
        if extension:
            ext = extension if extension.find('.')>-1 else '.'+extension
        else:
            ext = ''
        return self.data_dir+self.report_name()+ext


    def __dt_rank(self,drug:PSObject,targets:list,with_effect:str,correct_by_consistency:bool):
        '''
        Input
        -----
        targets = [PSObject]

        Returns
        ------
        target_uid2corrected_rank = {target_uid:rank},\n
        if "correct_by_consistency" = True rank from "self.target_uid2rank" is corrected by "self.dt_consist" 
        '''
        targets_uids = ResnetGraph.uids(targets)
        target_uid2corrected_rank = {k:v for k,v in self.target_uid2rank.items() if k in targets_uids}
        if correct_by_consistency:
            for target in targets:
                correction = self.dt_consist.consistency_coefficient(drug,target,with_effect)         
                target_uid2corrected_rank[target.uid()] *= correction

        return target_uid2corrected_rank
    

    def __add_rank2(self,drug_df:df,from_drug_graph:ResnetGraph):
        """
Adds columns DRUG2TARGET_REGULATOR_SCORE to drug_df.  Used for SNEA
        """
        rank2add2drug = dict()
        for node in from_drug_graph._get_nodes():
            if node['ObjTypeName'][0] == 'SmallMol':
                direct_target_uids, indirect_target_uids = from_drug_graph.direct_indirect_targets(
                node.uid(),
                direct_reltypes=['DirectRegulation'],
                indirect_reltypes=['Regulation','MolTransport','Expression']
                )

                if direct_target_uids or indirect_target_uids:
                    rank2add2drug[node.name()] = node[DRUG2TARGET_REGULATOR_SCORE]

        df_copy = df.copy_df(drug_df)
        if DRUG2TARGET_REGULATOR_SCORE not in list(df_copy.columns):
        # intitalizing columns with default values if df was not visited:
            df_copy[DRUG2TARGET_REGULATOR_SCORE] = NaN

        for idx in df_copy.index:
            drug_name = df_copy.at[idx,'Name']
            try:
                rank2add = rank2add2drug[drug_name]
                new_rank = nan_to_num(df_copy.at[idx,DRUG2TARGET_REGULATOR_SCORE]) + rank2add
                df_copy.at[idx,DRUG2TARGET_REGULATOR_SCORE] = new_rank
            except KeyError:
                continue

        df_copy._name_ = drug_df._name_
        return df_copy


    def __add_rank_targets2(self,drug_df:df,from_drug_graph:ResnetGraph):
        """
Adds following columns to drug_df:\n
DRUG2TARGET_REGULATOR_SCORE,\n'Directly inhibited targets',\n'Indirectly inhibited targets',\n'Directly activated targets',\n'Indirectly activated targets'
        """
        columns2add2drug = dict()
        for node_id, node in from_drug_graph.nodes(data=True):
            if node['ObjTypeName'][0] == 'SmallMol': # node is dict here
                direct_target_ids, indirect_target_ids = from_drug_graph.direct_indirect_targets(node_id,
                direct_reltypes=['DirectRegulation'],
                indirect_reltypes=['Regulation','MolTransport','Expression'])

                if direct_target_ids:
                    direct_target_names = [name[0] for t,name in from_drug_graph.nodes(data='Name') if t in direct_target_ids]
                    direct_target_names_str = ','.join(direct_target_names)
                    [self.direct_target2drugs.update_with_value(n,node['Name'][0]) for n in direct_target_names]
                else:
                    direct_target_names_str = ''

                if indirect_target_ids:
                    indirect_target_names = [name[0] for t,name in from_drug_graph.nodes(data='Name') if t in indirect_target_ids]
                    indirect_target_names_str = ','.join(indirect_target_names)
                    [self.indirect_target2drugs.update_with_value(n,node['Name'][0]) for n in indirect_target_names]
                else:
                    indirect_target_names_str = ''

                if direct_target_names_str or indirect_target_names_str:
                    columns2add2drug[node['Name'][0]] = [node[DRUG2TARGET_REGULATOR_SCORE],node['Drug Class'][0],direct_target_names_str,indirect_target_names_str]

        df_copy = df.copy_df(drug_df)
        if DRUG2TARGET_REGULATOR_SCORE not in list(df_copy.columns):
        # intitalizing columns with default values if df was not visited:
            df_copy[DRUG2TARGET_REGULATOR_SCORE] = NaN
            df_copy['Directly inhibited targets'] = NaN
            df_copy['Directly activated targets'] = NaN
            df_copy['Indirectly inhibited targets'] = NaN
            df_copy['Indirectly activated targets'] = NaN

        
        def __column_name__(drug_class:str,direct_targets:bool):
            first_word = 'Directly ' if direct_targets else 'Indirectly '
            second_word = 'inhibited ' if drug_class=='antagonist' else 'activated '
            return first_word+second_word+'targets'

        for idx in df_copy.index:
            drug_name = df_copy.at[idx,'Name']
            try:
                columns2add = columns2add2drug[drug_name]
                my_rank = nan_to_num(df_copy.at[idx,DRUG2TARGET_REGULATOR_SCORE])
                add_rank = columns2add[0]
                new_rank = my_rank+add_rank
                df_copy.at[idx,DRUG2TARGET_REGULATOR_SCORE] = new_rank
                drug_class = columns2add[1]
                df_copy.loc[idx,__column_name__(drug_class,True)] = columns2add[2]
                df_copy.loc[idx,__column_name__(drug_class,False)] = columns2add[3]
            except KeyError:
                continue

        df_copy._name_ = drug_df._name_
        return df_copy


    def addrank2(self,to_drug_df:df,from_drug_graph:ResnetGraph):
        if self.add_targets2drugs_ws:
            return self.__add_rank_targets2(to_drug_df,from_drug_graph)
        else:
            return self.__add_rank2(to_drug_df,from_drug_graph)
        

    def subtract_antigraph(self,my_df:df, antigraph:ResnetGraph):
        drug_objs = antigraph.psobjs_with(only_with_values=['SmallMol'])
        drug2rank = {n['Name'][0]:n[DRUG2TARGET_REGULATOR_SCORE] for n in drug_objs}
        df_copy = df.copy_df(my_df)

        for idx in df_copy.index:
            drug_name = df_copy.at[idx,'Name']
            try:
                rank2subtract = drug2rank[drug_name]
                my_rank = df_copy.at[idx,DRUG2TARGET_REGULATOR_SCORE]
                new_rank = my_rank-rank2subtract
                df_copy.at[idx,DRUG2TARGET_REGULATOR_SCORE] = new_rank
            except KeyError:
                continue

        return df_copy


    def __rank_drugs4(self,targets:list,with_effect:str,correct_by_consistency:bool):
        """
        Input
        -----
        targets - [PSObject]
        correct4consistency must be supplied explicitly to avoid consitency correction for toxicities in antigraph
        
        Returns
        -------
        ResnetGraph with drugs annotated with DRUG2TARGET_REGULATOR_SCORE and drug class (agonist/antagonist)
        """
        direct_effect_graph = self.drugs2targets.neighborhood(targets,[],['DirectRegulation'],[with_effect])
        indirect_effect_graph = self.drugs2targets.neighborhood(targets,[],['Regulation','Expression','MolTransport'],[with_effect])
        drug_class = 'agonist' if with_effect == 'positive' else 'antagonist'

        # initializing drug ranks
        drug2rank = dict()
        for drug in direct_effect_graph.psobjs_with(only_with_values=['SmallMol']):
            target_ranks4drug = self.__dt_rank(drug,targets,with_effect,correct_by_consistency)
            drug2rank[drug.uid()] = 2 * direct_effect_graph.rank_regulator(drug,target_ranks4drug)

        # adding indirect targets to ranking with smaller weight
        for drug in indirect_effect_graph.psobjs_with(only_with_values=['SmallMol']):
            target_ranks4drug = self.__dt_rank(drug,targets,with_effect,correct_by_consistency)
            try: my_rank = drug2rank[drug.uid()]
            except KeyError: my_rank = 0.0
            drug2rank[drug.uid()] = my_rank + indirect_effect_graph.rank_regulator(drug,target_ranks4drug)

        all_reg_graph = direct_effect_graph.copy()
        all_reg_graph = all_reg_graph.compose(indirect_effect_graph)
        nx.set_node_attributes(all_reg_graph,drug2rank,DRUG2TARGET_REGULATOR_SCORE)

        all_drugs = all_reg_graph._psobjs_with('SmallMol','ObjTypeName')
        druguid2class = {d.uid():[drug_class] for d in all_drugs}
        nx.set_node_attributes(all_reg_graph,druguid2class,'Drug Class')
        
        return all_reg_graph


    def load_target_ranks(self,from_ranked_targets_df:df,for_antagonists=True):
        '''
        Input
        -----
        from_ranked_targets_df must have RANK column

        Loads
        -----
        self.targets4agonists\nself.targets4antagonists
        '''
        my_targets = self._targets()
        name2target = {n.name():n for n in my_targets}
        target_name_column = self.__resnet_name__ if self.__resnet_name__ in from_ranked_targets_df.columns else 'Name'
        
        df_targets = set()
        missing_targets_counter = 0
        for idx in from_ranked_targets_df.index:
            target_name = from_ranked_targets_df.at[idx,target_name_column]
            target_rank = from_ranked_targets_df.at[idx,RANK]
            try:
                target = name2target[target_name]
                df_targets.add(target)
                target_uid = target.uid()
                self.target_uid2rank[target_uid] = float(target_rank)
                children = self.Graph._get_node(target_uid).childs() # only self.Graph is guaranteed to have CHILDS annotation
                for child in children:
                    # boosting drug ranking by adding target children to ranking network with the rank of the parent
                    if child.name() not in name2target:
                        self.target_uid2rank[child.uid()] = float(target_rank) 
                        df_targets.add(child)
            except KeyError:
                print(f'Node named "{target_name}"  is not found in self.Graph') 
                # to skip WEIGHTS row in df
                missing_targets_counter += 1
                continue
        print(f'{missing_targets_counter} targets out of {len(my_targets)} were not found in self.Graph')

        if for_antagonists:
            self.targets4antagonists = list(df_targets)
        else:
            self.targets4agonists = list(df_targets)


    def init_drug_df(self, drugs:list):
        '''
        Input
        -----
        drugs = [PSObject], used to load drugs from SNEA samples
        
        Returns
        -------
        new self.RefCountPandas
        '''
        forbidden_drugs =['DMSO', 'glucose','nicotine']
        print('Loading drug ranking worksheet')
        drug2pharmapendium_id = dict()
        clean_drugs = list(drugs)
        # bad fix for PS ontology problem:
        for drug in drugs:
            drug_name = drug.name()
            if 'antigen' in drug_name:
                clean_drugs.remove(drug)
            if 'vaccine' in drug_name: 
                clean_drugs.remove(drug)
            if drug_name in forbidden_drugs: 
                clean_drugs.remove(drug)
            try:
                drug2pharmapendium_id[drug_name] = drug['PharmaPendium ID'][0]
            except KeyError:
               continue

        drug_df = self.load_df(clean_drugs,max_child_count=None,max_threads=50) 
        # drug children for semantic refcount are loaded by "link2disease_concepts" function
        self.RefCountPandas = drug_df.merge_dict(drug2pharmapendium_id,new_col=PHARMAPENDIUM_ID,map2column='Name')
        print('Initialized "Drugs" worksheet with %d drugs from database for ranking' % len(self.RefCountPandas))


    def regulatory_rank(self):
        '''
        Returns
        -------
        subset of self.RefCountPandas with drugs having positive regulatory score\n
        returned df has new column DRUG2TARGET_REGULATOR_SCORE
        '''
        need_correction = True if self.params['consistency_correction4target_rank'] else False    
        antagonist_graph = self.__rank_drugs4(self.targets4antagonists,with_effect='negative',correct_by_consistency=need_correction)
        agonist_graph=self.__rank_drugs4(self.targets4agonists,with_effect='positive',correct_by_consistency=need_correction) 

        antagonist_antigraph = self.__rank_drugs4(self.targets4antagonists,with_effect='positive',correct_by_consistency=False)
        agonist_antigraph = self.__rank_drugs4(self.targets4agonists,with_effect='negative',correct_by_consistency=False)

        new_df = self.addrank2(self.RefCountPandas,antagonist_graph)
        new_df = self.addrank2(new_df,agonist_graph)
        new_df = self.subtract_antigraph(new_df,antagonist_antigraph)
        new_df = self.subtract_antigraph(new_df,agonist_antigraph)

        new_df = new_df.greater_than(0,DRUG2TARGET_REGULATOR_SCORE)
        new_df._name_ = f'Drugs for {self.input_names()}'
        print('Found %d drugs for %d antagonist targets and %d agonist targets for "%s"' % 
            (len(new_df),len(self.targets4antagonists),len(self.targets4agonists),self.input_names()),flush=True)
        return new_df
        

    def link2disease_concepts(self,in_drugdf:df):
        if not self.input_diseases: return in_drugdf # for SNEA.make_drugs_df()

        if self.__temp_id_col__ not in in_drugdf.columns:
        # case when drug_df was loaded from cache RNEF file:
            drug_df = self.add_temp_id(in_drugdf,max_children_count=0) 
            # max_children_count must be zero to merge with results of in_drugdf.etm_refs2df()
        else:
            drug_df = df.copy_df(in_drugdf)

        print('Linking %d drugs with %s by ClinicalTrial' % (len(drug_df),self._disease2str()),flush=True)
        colname = self._disease2str()+' clinical trials'
        how2connect = self.set_how2connect(['ClinicalTrial'],[],'')
        linked_row_count,linked_entities,drug_df = self.link2concept(
            colname,self.input_diseases,drug_df,how2connect,REFERENCE_IDENTIFIERS)
        print('%d drugs on clinical trials for %s' % (linked_row_count,self._disease2str()),flush=True)
        self.column_ranks[5] = list(drug_df.columns)[-1]

        print('Linking %d drugs with %s by Regulation' % (len(drug_df),self._disease2str()),flush=True)
        colname = 'regulation of '+ self._disease2str()
        how2connect = self.set_how2connect(['Regulation'],[],'')
        linked_row_count,linked_entities,drug_df = self.link2concept(
            colname,self.input_diseases,drug_df,how2connect,REFERENCE_IDENTIFIERS)
        print('%d drugs regulating %s' % (linked_row_count,self._disease2str()),flush=True)
        self.column_ranks[4] = list(drug_df.columns)[-1]

        print('Linking %d drugs with %d Symptoms linked to %s' % (len(drug_df),len(self.input_symptoms),self._disease2str()),flush=True)
        colname = 'symptoms for '+ self._disease2str()
        how2connect = self.set_how2connect(['Regulation'],[],'')
        linked_row_count,linked_entities,drug_df = self.link2concept(
            colname,self.input_symptoms,drug_df,how2connect,REFERENCE_IDENTIFIERS)
        print('%d drugs linked to symptoms for %s' % (linked_row_count,self._disease2str()),flush=True)
        self.column_ranks[1] = list(drug_df.columns)[-1]

        print('Linking %d drugs with %d ClinicalParameters linked to %s' % (len(drug_df),len(self.input_clinpars),self._disease2str()),flush=True)
        colname = 'Clinical parameters for '+ self._disease2str()
        how2connect = self.set_how2connect(['Regulation'],[],'')
        linked_row_count,linked_entities,drug_df = self.link2concept(
            colname,self.input_clinpars,drug_df,how2connect,REFERENCE_IDENTIFIERS)
        print('%d drugs linked to clnical parameters for %s' % (linked_row_count,self._disease2str()),flush=True)
        self.column_ranks[2] = list(drug_df.columns)[-1]

        print('Linking %d drugs with %d CellProcess linked to %s' % (len(drug_df),len(self.input_cellprocs),self._disease2str()),flush=True)
        colname = 'Cell processess affected by '+ self._disease2str()
        how2connect = self.set_how2connect(['Regulation'],[],'')
        linked_row_count,linked_entities,drug_df = self.link2concept(
            colname,self.input_cellprocs,drug_df,how2connect,REFERENCE_IDENTIFIERS)
        print('%d drugs linked to cell processes affected in %s' % (linked_row_count,self._disease2str()),flush=True)
        self.column_ranks[3] = list(drug_df.columns)[-1]
        return drug_df


    def score_drugs(self,normalize=True):
        '''
        Returns
        -------
        df with ranked drugs. df._name_ ='Drugs' if normalize, else 'rawDrugs'\n
        SNEA "normalize" settings is False 
        '''
        start_time = time.time()
        self.load_target_ranks(self.report_pandas[AGONIST_TARGETS_WS],for_antagonists=False)
        #if targets are in both worksheets will assume that it needs antagonist:
        self.load_target_ranks(self.report_pandas[ANTAGONIST_TARGETS_WS]) 
        ranked_df = self.regulatory_rank()
        self.column_ranks = {0:DRUG2TARGET_REGULATOR_SCORE}

        if 'disease' in self.params.keys(): # case of drug repositioning using input disease network and targets from knowledge graph
            if self.params['add_bibliography']:
                with ThreadPoolExecutor(max_workers=30,thread_name_prefix='ScoreDrugs') as sd:
                    drugs_etm_thread = sd.submit(self.etm_refs2df,ranked_df,self.input_names()) # adds etm refs to drug_df
                    concepts_thread = sd.submit(self.link2disease_concepts,ranked_df) # links drugs to disease-related concepts
                    full_drug_df = concepts_thread.result()
                    drugs_df_with_etmrefs = drugs_etm_thread.result()
                    #merging results of two threads:
                    etm_ref_colname = self.etm_ref_column_name('Name',self.input_names())
                    full_drug_df = full_drug_df.merge_df(drugs_df_with_etmrefs,how='left',on='Name',columns=[etm_ref_colname])
            else:
                full_drug_df = self.link2disease_concepts(ranked_df)
        else:
            full_drug_df = ranked_df # case of finding drugs from SNEA results

        raw_drug_df = self.make_count_df(full_drug_df,'rawDrugs')
        self.add2raw(raw_drug_df)

        if normalize:
            ranks = sorted(self.column_ranks)
            columns2norm = [self.column_ranks[k] for k in ranks]
            ranked_drugs_df,normalized_raw_df = self.normalize('rawDrugs','Drugs','Name',columns2norm)

            #reordering drug df columns:
            drugdf_columns = list(ranked_drugs_df.columns)
            drugdf_columns.remove(PHARMAPENDIUM_ID)
            drugdf_columns = [PHARMAPENDIUM_ID]+drugdf_columns
            ranked_drugs_df = ranked_drugs_df.reorder(drugdf_columns)
            
            ranked_drugs_df.tab_format['tab_color'] = 'green'
            self.report_pandas['Drugs'] = ranked_drugs_df

            print("Drug ranking was done in %s" % self.execution_time(start_time)[0], flush=True)
            print('Normalized worksheet named "Drugs" was added to report')
            return ranked_drugs_df
        else:
            print('Ranked drugs are in worksheet "rawDrugs" in raw data')
            return raw_drug_df


    def load_drugs(self,limit2drugs:set={}):
        '''
        Input
        -----
        self.__targets__ - [PSObject] that has to be filled prior to using this function
        limit2drugs - {PSObject}, used by SNEA make_drugs_df

        Output
        ------
        [PSObject] - list of drugs used by self.dt_consist.load_confidence_dict() and self.init_drug_df()
        '''
        self.drugs2targets = self.dt_consist.load_drug_graph(self._targets(),limit2drugs)
        drugs_linked2targets = self.drugs2targets._psobjs_with('SmallMol','ObjTypeName')
        print(f'Found {len(drugs_linked2targets)} drugs linked to {len(self._targets())} targets for ranking')

        drugs_withno_targets = []
        if hasattr(self,'drugs_linked2disease'):
            # case of drug repurposing for disease model in database
            drugs_withno_targets = set(self.drugs_linked2disease).difference(drugs_linked2targets)
            if drugs_withno_targets:
                my_api_session = APISession(self.APIconfig,NO_REL_PROPERTIES)
                my_api_session.entProps = ['Name','PharmaPendium ID']
                oql_query = 'SELECT Entity WHERE URN = ({props})'
                oql_query += " AND InOntology (SELECT Annotation WHERE Ontology='Pathway Studio Ontology' AND Relationship='is-a') under (SELECT OntologicalNode WHERE Name = drugs)"
                r_n = f'Find drugs linked to disease that have no targets linked to disease'
                drugs_withno_targets_urns = ResnetGraph.urns(drugs_withno_targets)
                drugs4disease = self.iterate_oql(oql_query,drugs_withno_targets_urns,request_name=r_n)
                drugs_withno_targets = drugs4disease._get_nodes()
                
                if drugs_withno_targets:
                    print('Found additional %d drugs linked to "%s" but not linked to its targets' % 
                        (len(drugs_withno_targets),self._disease2str()))
        
        return drugs_linked2targets+list(drugs_withno_targets)
        

    def init_load_score(self):
        my_drugs = self.load_drugs()
        if self.params['consistency_correction4target_rank']:
            with ThreadPoolExecutor(max_workers=100, thread_name_prefix='DrugsConsitency') as executor:
                future1 = executor.submit(self.dt_consist.load_confidence_dict,self._targets())
                future2 = executor.submit(self.init_drug_df,my_drugs)
                future1.result()
                future2.result()
 
            with ThreadPoolExecutor(max_workers=10, thread_name_prefix='ScoreDrugsSaveCache') as executor:
                f2 = executor.submit(self.score_drugs)
                f1 = executor.submit(self.dt_consist.save_network)
                f1.result()
                f2.result()
        else:
            self.init_drug_df(my_drugs)
            self.score_drugs()


    def find_rank_targets(self):
        super().make_report()
        self.report_pandas[ANTAGONIST_TARGETS_WS] = self.remove_high_level_entities(self.report_pandas[ANTAGONIST_TARGETS_WS])
        self.report_pandas[AGONIST_TARGETS_WS] = self.remove_high_level_entities(self.report_pandas[AGONIST_TARGETS_WS])
        print('Found %d targets for antagonists, %d targets for agonitsts, %d targets with uknown disease state'%
        (len(self.report_pandas[ANTAGONIST_TARGETS_WS])-1, len(self.report_pandas[AGONIST_TARGETS_WS])-1,len(self.report_pandas[UNKNOWN_TARGETS_WS])-1))
        print('%d targets in total were ranked' % len(self._targets()))
        return self._targets()


    def make_report(self):
        start_time = time.time()
        futures = list()
        with ThreadPoolExecutor(max_workers=30, thread_name_prefix='RankTargets:0,LoadDrugs:1') as executor:
            futures.append(executor.submit(self.find_rank_targets)) # creating target ranking worksheets
            futures.append(executor.submit(self.dt_consist.load_cache))
            [future.result() for future in as_completed(futures)]

        if self.params['add_bibliography']:
            with ThreadPoolExecutor(max_workers=20, thread_name_prefix='RankDrugs:0,ETM4targets:1') as executor:
                future1 = executor.submit(self.init_load_score) #finds and scores drugs
                future2 = executor.submit(self.add_etm_bibliography) # adds etm refs to target ranking worksheets
                future1.result()
                future2.result()
        else:
            self.init_load_score()
        
        self.report_pandas[ANTAGONIST_TARGETS_WS] = self.report_pandas[ANTAGONIST_TARGETS_WS].merge_psobject(self.direct_target2drugs,'Directly Inhibited by','Name',values21cell=True)
        self.report_pandas[ANTAGONIST_TARGETS_WS] = self.report_pandas[ANTAGONIST_TARGETS_WS].merge_psobject(self.indirect_target2drugs,'Indirectly Inhibited by','Name',values21cell=True)
        
        self.report_pandas[AGONIST_TARGETS_WS] =  self.report_pandas[AGONIST_TARGETS_WS].merge_psobject(self.direct_target2drugs,'Directly Activated by','Name',values21cell=True)
        self.report_pandas[AGONIST_TARGETS_WS] = self.report_pandas[AGONIST_TARGETS_WS].merge_psobject(self.indirect_target2drugs,'Indirectly Activated by','Name',values21cell=True)
       # self.drug2target_network.rnef2file(self.path2cache(self.cache_name),self.cache_ent_props,self.cache_rel_props)
        print('Drug repurposing for %s was done in %s' % (self._disease2str(), self.execution_time(start_time)))
       

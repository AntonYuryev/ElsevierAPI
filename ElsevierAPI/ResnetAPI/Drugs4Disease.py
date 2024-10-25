from .DiseaseTargets import DiseaseTargets,execution_time
from .DiseaseTargets import ANTAGONIST_TARGETS_WS,AGONIST_TARGETS_WS,REFERENCE_IDENTIFIERS,UNKNOWN_TARGETS_WS
from .ResnetGraph import ResnetGraph,PSObject
from .SemanticSearch import RANK
from .ResnetAPISession import APISession,NO_REL_PROPERTIES,BIBLIO_PROPERTIES,ALL_CHILDS
from ..pandas.panda_tricks import df,np
from numpy import nan_to_num
import time
import networkx as nx
from concurrent.futures import ThreadPoolExecutor
from .DrugTargetConfidence import DrugTargetConsistency

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
                'consistency_correction4target_rank':False,
                "max_childs" : 11,
                'drug_groups':[]
                }
        my_kwargs.update(kwargs)

        super().__init__(*args,**my_kwargs)
        self.clone_all_sessions = my_kwargs.pop('clone',True) # must be different from clone that is used by SemanticSearch

        #dt_consist_kwargs = dict(**my_kwargs)
        self.dt_consist = DrugTargetConsistency(self.APIconfig,**my_kwargs)
        self.target_uid2rank = dict() # {target_dbid:rank}
        self.column_ranks = dict() # {column rank:column name} rank defines the column position in Drugs df and consequently its weight
        self.direct_target2drugs = PSObject() # used for annotation of ANTAGONIST_TARGETS_WS,AGONIST_TARGETS_WS with drugs
        self.indirect_target2drugs = PSObject() # used for annotation of ANTAGONIST_TARGETS_WS,AGONIST_TARGETS_WS with drugs
        self.add_targets2drugs_ws = True


    @classmethod
    def from_files(cls,**kwargs):
        fake_config = {'ResnetURL':'','PSuserName':'','PSpassword':'','ETMURL':'','ETMapikey':''}
        my_kwargs = dict(kwargs)
        my_kwargs['connect2server'] = False
        dt = cls(fake_config,**my_kwargs)
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


    def __rank(self,targets:list,_4drug:PSObject,with_effect:str,correct_by_consistency:bool)->dict[int,float]:
        '''
        Input
        -----
        targets = [PSObject]

        Returns
        ------
        target_uid2corrected_rank = {target_uid:rank},\n
        if "correct_by_consistency" rank from "self.target_uid2rank" corrected by "self.dt_consist" 
        '''
        targets_uids = ResnetGraph.uids(targets)
        target_uid2corrected_rank = {k:v for k,v in self.target_uid2rank.items() if k in targets_uids}
        meanpX = 6.7 # this value was calculated as average pX of all DirectRegulation from Reaxys
        for target in targets:
            target_uid = target.uid()
            dt_rels = self.drugs2targets._psrels4(_4drug.uid(),target_uid)
            if dt_rels: # dt_rels may be empty because not every target from df is linked to drug
                if dt_rels[0].isprimarytarget():
                    pX = float(dt_rels[0].get_prop('Affinity',0,meanpX))
                    pXcorrection = 1.5 + pX/12.0
                    target_uid2corrected_rank[target_uid] *= pXcorrection

                if correct_by_consistency:
                    consistency_correction = self.dt_consist.consistency_coefficient(_4drug,target,with_effect)
                    target_uid2corrected_rank[target_uid] *= consistency_correction

        return target_uid2corrected_rank
    

    def __add_rank2(self,drug_df:df,from_drug_graph:ResnetGraph):
        """
Adds columns DRUG2TARGET_REGULATOR_SCORE to drug_df.  Used for SNEA
        """
        rank2add2drug = dict()
        for node in from_drug_graph._get_nodes():
            if node.objtype() == 'SmallMol':
                direct_target_uids, indirect_target_uids = from_drug_graph.direct_indirect_targets(node.uid())

                if direct_target_uids or indirect_target_uids:
                    rank2add2drug[node.name()] = float(node.get_prop(DRUG2TARGET_REGULATOR_SCORE))

        df_copy = drug_df.dfcopy()
        if DRUG2TARGET_REGULATOR_SCORE not in list(df_copy.columns):
        # intitalizing columns with default values if df was not visited:
            df_copy[DRUG2TARGET_REGULATOR_SCORE] = np.nan

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


    def __add_rank_targets2(self,drug_df:df,drug2targets:ResnetGraph):
        """
Adds following columns to drug_df:\n
DRUG2TARGET_REGULATOR_SCORE,\n'Directly inhibited targets',\n'Indirectly inhibited targets',\n'Directly activated targets',\n'Indirectly activated targets'
        """
        columns2add2drug = dict() # collecting stats:
        for  drug in drug2targets._get_nodes():
            if drug.objtype() == 'SmallMol': # node is dict here
                direct_targets, indirect_targets = drug2targets.direct_indirect_targets(drug.uid())

                direct_target_names = list()  #direct_targets are sorted by pX
                for target_id, px, consistency in direct_targets:
                    direct_target = drug2targets._get_node(target_id)
                    target_name = direct_target.name()
                    self.direct_target2drugs.update_with_value(target_name,drug.name())
                    target_name += ' ('
                    if float(px) > 0.0:
                        target_name += f'pX={px};'
                    target_name += f'cohesion={consistency})'
                    direct_target_names.append(target_name)

                direct_target_names_str = ','.join(direct_target_names)
                #[self.direct_target2drugs.update_with_value(n,drug.name()) for n in direct_target_names]

                indirect_target_names = list()
                for target_id, px, consistency in indirect_targets:
                    indirect_target = drug2targets._get_node(target_id)
                    target_name = indirect_target.name()
                    self.indirect_target2drugs.update_with_value(target_name,drug.name())
                    target_name += ' ('
                    if float(px) > 0.0:
                        target_name += f'pX={px};'
                    target_name += f'cohesion={consistency})'
                    indirect_target_names.append(target_name)

                indirect_target_names_str = ','.join(indirect_target_names)
               # [self.indirect_target2drugs.update_with_value(n,drug.name()) for n in indirect_target_names]

                if direct_target_names_str or indirect_target_names_str:
                    drug_score = float(drug.get_prop(DRUG2TARGET_REGULATOR_SCORE))
                    drug_class = str(drug.get_prop('Drug Class'))
                    columns2add2drug[drug.name()] = [drug_score,drug_class,direct_target_names_str,indirect_target_names_str]

        # making dataframe
        df_copy = drug_df.dfcopy()
        if DRUG2TARGET_REGULATOR_SCORE not in list(df_copy.columns):
        # intitalizing columns with default values if df was not visited:
            df_copy[DRUG2TARGET_REGULATOR_SCORE] = np.nan
            df_copy['Directly inhibited targets'] = ''
            df_copy['Directly activated targets'] = ''
            df_copy['Indirectly inhibited targets'] = ''
            df_copy['Indirectly activated targets'] = ''

        
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
        drug2rank = {d.name():float(d.get_prop(DRUG2TARGET_REGULATOR_SCORE)) for d in drug_objs}
        df_copy = my_df.dfcopy()

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
        my_dtG = self.drugs2targets.neighborhood(set(targets),[],[],[with_effect])
        my_drugs = my_dtG._psobjs_with('SmallMol','ObjTypeName')
        
        # initializing drug ranks
        drug2rank = {d.uid():[my_dtG.rank_regulator(d,self.__rank(targets,d,with_effect,correct_by_consistency))] for d in my_drugs}
        nx.set_node_attributes(my_dtG,drug2rank,DRUG2TARGET_REGULATOR_SCORE)

        drug_class = 'agonist' if with_effect == 'positive' else 'antagonist'
        druguid2class = {d.uid():[drug_class] for d in my_drugs}
        nx.set_node_attributes(my_dtG,druguid2class,'Drug Class')
        
        return my_dtG


    def load_target_ranks(self,from_ranked_targets_df:df,for_antagonists=True):
        '''
        Input
        -----
        from_ranked_targets_df must have RANK column

        Loads
        -----
        self.targets4agonists\nself.targets4antagonists\nself.target_uid2rank
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


    def init_drug_df(self, drugs:list[PSObject]):
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
            
            pp_id = str(drug.get_prop('PharmaPendium ID'))
            if pp_id:
                drug2pharmapendium_id[drug_name] = pp_id

        count = len(clean_drugs)
        clean_drugs = [d for d in clean_drugs if d not in self.disease_inducers]
        clean_drugs = [d for d in clean_drugs if "'" not in d.urn()] # drugs with "'" in URN are not searchable in PS
        removed_count = count-len(clean_drugs)
        print(f'{removed_count} drugs were removed because they are known to induce {self._disease2str()}')
        drug_df = self.load_df(clean_drugs,max_childs=ALL_CHILDS,max_threads=50) 
        # drug children for semantic refcount are loaded by "link2disease_concepts" function
        self.RefCountPandas = drug_df.merge_dict(drug2pharmapendium_id,new_col=PHARMAPENDIUM_ID,map2column='Name')
        print('Initialized "Drugs" worksheet with %d drugs from database for ranking' % len(self.RefCountPandas))
        return


    def regulatory_rank(self):
        '''
        Returns
        -------
        subset of self.RefCountPandas with drugs having positive regulatory score\n
        returned df has new column DRUG2TARGET_REGULATOR_SCORE
        '''
        need_correction = True if self.params['consistency_correction4target_rank'] else False    
        antagonist_graph = self.__rank_drugs4(self.targets4antagonists,with_effect='negative',correct_by_consistency=need_correction)
        agonist_graph = self.__rank_drugs4(self.targets4agonists,with_effect='positive',correct_by_consistency=need_correction) 

        antagonist_antigraph = self.__rank_drugs4(self.targets4antagonists,with_effect='positive',correct_by_consistency=False)
        agonist_antigraph = self.__rank_drugs4(self.targets4agonists,with_effect='negative',correct_by_consistency=False)

        new_ranked_df = self.addrank2(self.RefCountPandas,antagonist_graph)
        new_ranked_df = self.addrank2(new_ranked_df,agonist_graph)
        new_ranked_df = self.subtract_antigraph(new_ranked_df,antagonist_antigraph)
        new_ranked_df = self.subtract_antigraph(new_ranked_df,agonist_antigraph)

        predicted_drugs = new_ranked_df.greater_than(0,DRUG2TARGET_REGULATOR_SCORE)['Name'].to_list()
        known_drugs = ResnetGraph.names(self.drugs_linked2disease)
        drugs2score = known_drugs + predicted_drugs
        new_ranked_df = new_ranked_df.filter_by(drugs2score,'Name')
        new_ranked_df[DRUG2TARGET_REGULATOR_SCORE] = new_ranked_df[DRUG2TARGET_REGULATOR_SCORE].fillna(0.0)

        disnames = ','.join(self.input_names())
        new_ranked_df._name_ = f'Drugs for {disnames}'
        new_ranked_df = new_ranked_df.sortrows(by=[DRUG2TARGET_REGULATOR_SCORE])
        print('Found %d drugs for %d antagonist targets and %d agonist targets for "%s"' % 
            (len(new_ranked_df),len(self.targets4antagonists),len(self.targets4agonists),self.input_names()),flush=True)
        return new_ranked_df
        

    def link2disease_concepts(self,in_drugdf:df):
        if not self.input_diseases: return in_drugdf # for SNEA.make_drugs_df()
        if len(in_drugdf) < 2:
            print('Drugs worksheet is empty !!!!')
            return in_drugdf

        if self.__temp_id_col__ not in in_drugdf.columns:
        # case when drug_df was loaded from cache RNEF file:
            drug_df = self.add_temp_id(in_drugdf,max_childs=ALL_CHILDS,max_threads=50) 
            # max_childs_count must be zero (=ALL_CHILDS) to avoid removal of drugs with children 
            # and to merge with results of SemanticSearch.bibliography() future
        else:
            drug_df = in_drugdf.dfcopy()

        print('\n\nLinking %d drugs with %s by ClinicalTrial' % (len(drug_df),self._disease2str()),flush=True)
        colname = self._disease2str()+' clinical trials'
        kwargs = {'connect_by_rels':['ClinicalTrial']}
        how2connect = self.set_how2connect(**kwargs)
        linked_row_count,linked_entities,drug_df = self.link2concept(
            colname,self.input_diseases,drug_df,how2connect)
        print('%d drugs on clinical trials for %s' % (linked_row_count,self._disease2str()),flush=True)
        self.column_ranks[5] = list(drug_df.columns)[-4]

        print('\n\nLinking %d drugs with %s by Regulation' % (len(drug_df),self._disease2str()),flush=True)
        colname = 'regulation of '+ self._disease2str()
        kwargs = {'connect_by_rels':['Regulation']}
        how2connect = self.set_how2connect(**kwargs)
        linked_row_count,linked_entities,drug_df = self.link2concept(
            colname,self.input_diseases,drug_df,how2connect)
        print('%d drugs regulating %s' % (linked_row_count,self._disease2str()),flush=True)
        self.column_ranks[4] = list(drug_df.columns)[-4]

        print('\n\nLinking %d drugs with %d Symptoms linked to %s' % (len(drug_df),len(self.input_symptoms),self._disease2str()),flush=True)
        colname = 'symptoms for '+ self._disease2str()
        kwargs = {'connect_by_rels':['Regulation']}
        how2connect = self.set_how2connect(**kwargs)
        linked_row_count,linked_entities,drug_df = self.link2concept(
            colname,self.input_symptoms,drug_df,how2connect,REFERENCE_IDENTIFIERS)
        print('%d drugs linked to symptoms for %s' % (linked_row_count,self._disease2str()),flush=True)
        self.column_ranks[1] = list(drug_df.columns)[-4]

        print('\n\nLinking %d drugs with %d ClinicalParameters linked to %s' % (len(drug_df),len(self.input_clinpars),self._disease2str()),flush=True)
        colname = 'Clinical parameters for '+ self._disease2str()
        kwargs = {'connect_by_rels':['Regulation']}
        how2connect = self.set_how2connect(**kwargs)
        linked_row_count,linked_entities,drug_df = self.link2concept(
            colname,self.input_clinpars,drug_df,how2connect,REFERENCE_IDENTIFIERS)
        print('%d drugs linked to clnical parameters for %s' % (linked_row_count,self._disease2str()),flush=True)
        self.column_ranks[2] = list(drug_df.columns)[-4]

        print('\n\nLinking %d drugs with %d CellProcess linked to %s' % (len(drug_df),len(self.input_cellprocs),self._disease2str()),flush=True)
        colname = 'Cell processess affected by '+ self._disease2str()
        kwargs = {'connect_by_rels':['Regulation']}
        how2connect = self.set_how2connect(**kwargs)
        linked_row_count,linked_entities,drug_df = self.link2concept(
            colname,self.input_cellprocs,drug_df,how2connect,REFERENCE_IDENTIFIERS)
        print('%d drugs linked to cell processes affected in %s' % (linked_row_count,self._disease2str()),flush=True)
        self.column_ranks[3] = list(drug_df.columns)[-4]
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
        self.load_target_ranks(self.report_pandas[ANTAGONIST_TARGETS_WS],for_antagonists=True) 
        ranked_df = self.regulatory_rank()
        self.column_ranks = {0:DRUG2TARGET_REGULATOR_SCORE}

        if 'disease' in self.params.keys(): # case of drug repositioning using input disease network and targets from knowledge graph
            if self.params['add_bibliography']:
                with ThreadPoolExecutor(max_workers=30,thread_name_prefix='ScoreDrugs') as sd:
                    futures = dict()
                    futures['bibliography'] = sd.submit(self.bibliography,ranked_df,self.input_names(),'Name',[],len(ranked_df)) # adds etm refs to drug_df
                    futures['link2disease_concepts'] = sd.submit(self.link2disease_concepts,ranked_df) # links drugs to disease-related concepts

                    try:
                        full_drug_df = futures['link2disease_concepts'].result()
                    except Exception as m:
                        print(f'"link2disease_concepts" finised with exception {m}' )
                        raise m
                    try:
                        drugs_df_with_etmrefs = futures['bibliography'].result()
                    except Exception as m:
                        print(f'"bibliography" finished with exception {m}' )
                        raise m
                    #merging results of two threads:
                    etm_ref_colname = self.etm_refcount_colname('Name',self.input_names())
                    full_drug_df = full_drug_df.merge_df(drugs_df_with_etmrefs,how='left',on='Name',columns=[etm_ref_colname])
            else:
                full_drug_df = self.link2disease_concepts(ranked_df)
        else:
            full_drug_df = ranked_df # case of finding drugs from SNEA results

        raw_drug_df = self.make_count_df(full_drug_df,'rawDrugs')
        self.add2raw(raw_drug_df)

        if normalize:
            ranks = sorted(self.column_ranks) # sorts keys in self.column_ranks = {rank:colname}
            columns2norm = [self.column_ranks[k] for k in ranks]
            ranked_drugs_df,_ = self.normalize('rawDrugs','Drugs','Name',columns2norm)

            #reordering drug df columns:
            drugdf_columns = list(ranked_drugs_df.columns)
            drugdf_columns.remove(PHARMAPENDIUM_ID)
            drugdf_columns = [PHARMAPENDIUM_ID]+drugdf_columns
            ranked_drugs_df = ranked_drugs_df.dfcopy(drugdf_columns)
            
            ranked_drugs_df.tab_format['tab_color'] = 'green'
            ranked_drugs_df._name_ = 'Drugs'
            self.add2report(ranked_drugs_df)

            print("Drug ranking was done in %s" % execution_time(start_time), flush=True)
            print('Normalized worksheet named "Drugs" was added to report')
        else:
            print('Ranked drugs are in worksheet "rawDrugs" in raw data')
            
        return raw_drug_df,ranked_drugs_df


    def load_drugs(self,limit2drugs=set())->list[PSObject]:
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
        drugs_linked2targets = set(self.drugs2targets._psobjs_with('SmallMol','ObjTypeName'))
        drugs_linked2targets = drugs_linked2targets.difference(self._targets()) # to remove metabolite targets
        print(f'Found {len(drugs_linked2targets)} drugs linked to {len(self._targets())} targets for ranking')

        drugs_withno_targets = []
        if hasattr(self,'drugs_linked2disease'):
            # adding drugs that are known to inhibit disease but have no links to disease targets
            drugs_withno_targets = list(set(self.drugs_linked2disease).difference(drugs_linked2targets))
            if drugs_withno_targets:
                my_api_session = APISession(self.APIconfig,NO_REL_PROPERTIES)
                my_api_session.entProps = ['Name','PharmaPendium ID']
                oql_query = 'SELECT Entity WHERE URN = ({props})'
                oql_query += " AND InOntology (SELECT Annotation WHERE Ontology='Pathway Studio Ontology' AND Relationship='is-a') under (SELECT OntologicalNode WHERE Name = drugs)"
                r_n = f'Find drugs linked to disease that have no targets linked to disease'
                drugs_withno_targets_urns = set(ResnetGraph.urns(drugs_withno_targets))
                drugs4disease = self.iterate_oql(oql_query,drugs_withno_targets_urns,request_name=r_n)
                drugs_withno_targets = drugs4disease._get_nodes()
                
                if drugs_withno_targets:
                    print('Found additional %d drugs linked to "%s" but not linked to its targets' % 
                        (len(drugs_withno_targets),self._disease2str()))
        
        return list(drugs_linked2targets) + drugs_withno_targets
        

    def init_load_score(self):
        my_drugs = self.load_drugs()
        if self.params['consistency_correction4target_rank']:
            load_fromdb = self.params.get('recalculate_dtconsistency',False)
            with ThreadPoolExecutor(max_workers=10, thread_name_prefix='Drugs_Consistency') as e1:
                futures = dict()
                futures['init_drug_df'] = e1.submit(self.init_drug_df,my_drugs)
                futures['load_confidence_dict'] = e1.submit(self.dt_consist.load_confidence_dict,self._targets(),{},load_fromdb)
                
                for func_name, future in futures.items():
                    try:
                        future.result()
                    except Exception as e:
                        print(f"{func_name} failed with exception:({e})")
                        raise(e)
                e1.shutdown()

            with ThreadPoolExecutor(max_workers=10, thread_name_prefix='ScoreDrugs_SaveCache') as e2:
                futures = dict()
                futures['score_drugs'] = e2.submit(self.score_drugs)
                futures['dt_consist.save_network'] = e2.submit(self.dt_consist.save_network)

                for func_name,future in futures.items():
                    try:
                        future.result()
                    except Exception as m:
                        print(f'{func_name} finished with exception {m}')
                        raise m
            e2.shutdown()

            self.dt_consist.clear()
        else:
            self.init_drug_df(my_drugs)
            self.score_drugs()


    def find_rank_targets(self):
        super().make_report()
        print('Found %d targets for antagonists, %d targets for agonitsts, %d targets with uknown disease state'%
        (len(self.report_pandas[ANTAGONIST_TARGETS_WS])-1, len(self.report_pandas[AGONIST_TARGETS_WS])-1,len(self.report_pandas[UNKNOWN_TARGETS_WS])-1))
        print('%d targets in total were ranked' % len(self._targets()))
        return self._targets()


    def make_report(self):
        start_time = time.time()
        self.find_rank_targets() # creating target ranking worksheets
        self.refs2report(ANTAGONIST_TARGETS_WS,self.input_names())
        future_dic = dict()
        with ThreadPoolExecutor(max_workers=2) as e:
            future_dic['init_load_score'] = e.submit(self.init_load_score) #finds and scores drugs
            future_dic['refs2report4'+ANTAGONIST_TARGETS_WS] = e.submit(self.refs2report,ANTAGONIST_TARGETS_WS,self.input_names())
            future_dic['refs2report4'+AGONIST_TARGETS_WS] = e.submit(self.refs2report,AGONIST_TARGETS_WS,self.input_names())
            future_dic['refs2report4'+UNKNOWN_TARGETS_WS] = e.submit(self.refs2report,UNKNOWN_TARGETS_WS,self.input_names())

        for name,future in future_dic.items():
            try:
               future.result()
            except Exception as x:
                print(f'{name} thread has finished with error {x}') 
                raise(x)
# this step cannot be multithreaded because previous steps mutate self.RefStats.ref_counter:        
        self.add_tm_bibliography_df()

        self.report_pandas[ANTAGONIST_TARGETS_WS] = self.report_pandas[ANTAGONIST_TARGETS_WS].merge_psobject(self.direct_target2drugs,'Directly Inhibited by','Name',values21cell=True)
        self.report_pandas[ANTAGONIST_TARGETS_WS] = self.report_pandas[ANTAGONIST_TARGETS_WS].merge_psobject(self.indirect_target2drugs,'Indirectly Inhibited by','Name',values21cell=True)
        
        self.report_pandas[AGONIST_TARGETS_WS] =  self.report_pandas[AGONIST_TARGETS_WS].merge_psobject(self.direct_target2drugs,'Directly Activated by','Name',values21cell=True)
        self.report_pandas[AGONIST_TARGETS_WS] = self.report_pandas[AGONIST_TARGETS_WS].merge_psobject(self.indirect_target2drugs,'Indirectly Activated by','Name',values21cell=True)

        try:
            drug_groups = self.params['drug_groups']
            if drug_groups:
                self.add_groups(self.report_pandas['Drugs'],drug_groups)
        except KeyError:
            pass
        print('Drug repurposing for %s was done in %s' % (self._disease2str(), execution_time(start_time)))

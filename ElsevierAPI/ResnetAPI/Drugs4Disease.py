from .DiseaseTargets import DiseaseTargets,execution_time
from .DiseaseTargets import ANTAGONIST_TARGETS_WS,AGONIST_TARGETS_WS,REFERENCE_IDENTIFIERS,UNKNOWN_TARGETS_WS
from .ResnetGraph import ResnetGraph,PSObject,OBJECT_TYPE
from .SemanticSearch import RANK
from .ResnetAPISession import APISession,OQL,NO_REL_PROPERTIES,BIBLIO_PROPERTIES,ALL_CHILDS
from ..pandas.panda_tricks import df,np
from numpy import nan_to_num
import time
import networkx as nx
from .DrugTargetConfidence import DrugTargetConsistency
from ..utils import run_tasks

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

        self.target_uid2rank = dict() # {target_dbid:rank}
        self.direct_target2drugs = PSObject() # used for annotation of ANTAGONIST_TARGETS_WS,AGONIST_TARGETS_WS with drugs
        self.indirect_target2drugs = PSObject() # used for annotation of ANTAGONIST_TARGETS_WS,AGONIST_TARGETS_WS with drugs
        self.add_targets2drugs_ws = True


    def load_dt(self,**kwargs):
        self.dt_consist = DrugTargetConsistency(self.APIconfig,**kwargs)
        return self.dt_consist

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
        input:
          from_ranked_targets_df must have RANK column
        output:
          self.targets4agonists\n
          self.targets4antagonists\n
          self.target_uid2rank
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

        reltype2targets = dict(self.params.get("add_inhibitors4",dict()))
        bad_drugs = set()
        all_targets2inhibit = set()
        rt = list()
        for reltype,target_names in reltype2targets.items():
          targets2expand = [n for n in self.__targets__ if n.name() in target_names]
          targets2expand_dbids = ResnetGraph.dbids(targets2expand)
          reltypes = [reltype] if reltype else []
          rt.append(reltype)
          oql_query = OQL.expand_entity(targets2expand_dbids,['Id'],reltypes,['SmallMol'],'upstream')
          oql_query += ' AND Effect = positive'
          expanded_targets_g = self.process_oql(oql_query,f'Find drugs activating {target_names} by {reltype}')
          bad_drugs.update(expanded_targets_g._psobjs_with('SmallMol',OBJECT_TYPE))
          all_targets2inhibit.update(target_names)
        
        if bad_drugs:
          before_count = len(clean_drugs)
          clean_drugs = [d for d in clean_drugs if d not in bad_drugs]
          print(f'Excluded {before_count-len(clean_drugs)} drugs that activate {all_targets2inhibit} by {rt}')
           
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
        new_ranked_df.col2rank[DRUG2TARGET_REGULATOR_SCORE] = 1
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
          before_mapping = len(drug_df)
          print(f'{before_mapping - len(drug_df)} rows were deleted because database identifier cannot be found for drug names')
          # max_childs_count must be zero (=ALL_CHILDS) to avoid removal of drugs with children 
          # and to merge with results of SemanticSearch.bibliography() future
      else:
          drug_df = in_drugdf.dfcopy()

      kwargs = {'connect_by_rels':['Regulation'],
                'with_effects':['positive'],
                'boost_with_reltypes':['FunctionalAssociation','Regulation'],
                'clone2retrieve' : REFERENCE_IDENTIFIERS,
                'init_refstat' : False,
                'column_name': 'Cell processess to activate'}
      drug_df = self.score_concept('processes2activate',drug_df,**kwargs)[2]

      kwargs['with_effects'] = ['negative']
      kwargs['column_name'] = 'Cell processess to inhibit'
      drug_df = self.score_concept('processes2inhibit',drug_df,**kwargs)[2]
        
      kwargs['with_effects'] = ['unknown']
      kwargs['column_name'] = 'Cell processess affected by '
      drug_df = self.score_concept('processes',drug_df,**kwargs)[2]

      kwargs['with_effects'] = ['positive']
      kwargs['column_name'] = 'Cells to activate'
      drug_df = self.score_concept('cells2activate',drug_df,**kwargs)[2]

      kwargs.pop('with_effects')
      kwargs['column_name'] = 'Clinical parameters for '+ self._disease2str()
      drug_df = self.score_concept('clinical_parameters',drug_df,**kwargs)[2]

      kwargs['column_name'] = 'symptoms for '+ self._disease2str()
      drug_df = self.score_concept('symptoms',drug_df,**kwargs)[2]

      kwargs['column_name'] = 'regulation of diseases similar to '+ self._disease2str()
      drug_df = self.score_concept('similar_diseases',drug_df,**kwargs)[2]

      kwargs['column_name'] = 'regulation of '+ self._disease2str()
      drug_df = self.score_concepts(self.input_diseases,drug_df,**kwargs)[2]

      kwargs['column_name'] = self._disease2str()+' clinical trials'
      drug_df = self.score_concepts(self.input_diseases,drug_df,**kwargs)[2]

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

        if 'disease' in self.params.keys(): # case of drug repositioning using input disease network and targets from knowledge graph
          if self.params['add_bibliography']:
            tasks = [(self.bibliography,(ranked_df,self.names4tmsearch(),'Name',[],len(ranked_df)))]
            tasks.append((self.link2disease_concepts,(ranked_df,)))
            results = run_tasks(tasks)
            full_drug_df = results['link2disease_concepts']
            drugs_df_with_tmrefs = results['bibliography']
            etm_ref_colname = self.tm_refcount_colname('Name',self.names4tmsearch())
            full_drug_df = full_drug_df.merge_df(drugs_df_with_tmrefs,how='left',on='Name',columns=[etm_ref_colname])
          else:
            full_drug_df = self.link2disease_concepts(ranked_df)
        else:
            full_drug_df = ranked_df # case of finding drugs from SNEA results

        raw_drug_df = self.make_count_df(full_drug_df,'rawDrugs')
        self.add2raw(raw_drug_df)

        if normalize:
          ranked_drugs_df,_ = self.normalize('rawDrugs','Drugs','Name')

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
        #self.dt_consist = self.load_dt() now runs in parallele
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
        tasks = [(self.init_drug_df,(my_drugs,))] # need comma after my_drugs to make it a tuple
        tasks.append((self.dt_consist.load_confidence_dict,(self._targets(),{},load_fromdb)))
        run_tasks(tasks)
        run_tasks([(self.score_drugs,()),(self.dt_consist.save_network,())])
        self.dt_consist.clear()
        return
      else:
        self.init_drug_df(my_drugs)
        self.score_drugs()
        return


    def find_rank_targets(self):
        super().make_report()
        print('Found %d targets for antagonists, %d targets for agonitsts, %d targets with uknown disease state'%
        (len(self.report_pandas[ANTAGONIST_TARGETS_WS])-1, len(self.report_pandas[AGONIST_TARGETS_WS])-1,len(self.report_pandas[UNKNOWN_TARGETS_WS])-1))
        print('%d targets in total were ranked' % len(self._targets()))
        return self._targets()


    def annotate_report(self):
      self.report_pandas[ANTAGONIST_TARGETS_WS] = self.report_pandas[ANTAGONIST_TARGETS_WS].merge_psobject(self.direct_target2drugs,'Directly Inhibited by','Name',values21cell=True)
      self.report_pandas[ANTAGONIST_TARGETS_WS] = self.report_pandas[ANTAGONIST_TARGETS_WS].merge_psobject(self.indirect_target2drugs,'Indirectly Inhibited by','Name',values21cell=True)
      
      self.report_pandas[AGONIST_TARGETS_WS] =  self.report_pandas[AGONIST_TARGETS_WS].merge_psobject(self.direct_target2drugs,'Directly Activated by','Name',values21cell=True)
      self.report_pandas[AGONIST_TARGETS_WS] = self.report_pandas[AGONIST_TARGETS_WS].merge_psobject(self.indirect_target2drugs,'Indirectly Activated by','Name',values21cell=True)
 
      drug_groups = self.params.get('drug_groups',[])
      if drug_groups:
          self.add_groups(self.report_pandas['Drugs'],drug_groups)


    def make_report(self):
      start_time = time.time()
      tasks = [(self.load_dt,())]
      tasks.append((self.find_rank_targets,()))
      run_tasks(tasks)
      
      if self.params['add_bibliography']:
        tasks = [(self.init_load_score,())]
        tasks.append((self.refs2report,(ANTAGONIST_TARGETS_WS,self.names4tmsearch())))
        tasks.append((self.refs2report,(AGONIST_TARGETS_WS,self.names4tmsearch())))
        tasks.append((self.refs2report,(UNKNOWN_TARGETS_WS,self.names4tmsearch())))
        run_tasks(tasks)
        self.add_tm_bibliography_df()
      else:
         self.init_load_score()
      
      self.annotate_report()
      print('Drug repurposing for %s was done in %s' % (self._disease2str(), execution_time(start_time)))

from .DiseaseTargets import DiseaseTargets,ResnetGraph,EFFECT
from .DiseaseTargets import ANTAGONIST_TARGETS_WS,AGONIST_TARGETS_WS,REFERENCE_IDENTIFIERS
from .SemanticSearch import RANK,BIBLIO_PROPERTIES
from .ResnetAPISession import APISession,NO_REL_PROPERTIES,OQL,PS_REFIID_TYPES
from ..pandas.panda_tricks import df,PSObject
from numpy import NaN,nan_to_num
import time
import networkx as nx
from threading import Thread
from concurrent.futures import ThreadPoolExecutor
from threading import Thread

PATHWAY_REGULATOR_SCORE = 'Disease model regulator score'
DRUG2TARGET_REGULATOR_SCORE = 'Regulator score'
PHARMAPENDIUM_ID = 'FDA-approved name'

class Drugs4Targets(DiseaseTargets):
    pass
    drug2target_consistency = dict() #{tuple(drug_urn,target_urn):consistency coefficient}
    target_ranks = dict() # {target_id:rank}
    column_ranks = dict() # {column rank:column name} rank defines the column position in Drugs df and consequently its weight
    direct_target2drugs = PSObject()
    indirect_target2drugs = PSObject()
    drug_ids = set() #{int}
    drugs = set() # {PSObject}
    add_targets2drugs_ws = True

    def __init__(self, APIconfig,input_parameters=dict(),clone=True,what2retrieve=BIBLIO_PROPERTIES):
        """
        self.set_target_disease_state() needs references
        """
        super().__init__(APIconfig,input_parameters,what2retrieve)
        self.clone_all_sessions = clone
        if not input_parameters:
            self.add_params({
                'disease':[],
                'strict_mode':True,
                'data_dir':'',
                'add_bibliography' : False,
                'strict_mode' : False,
                'target_types' : ['Protein'],
                'pathway_folders':[],
                'pathways': [],
                'consistency_correction4target_rank':False
                })

            
    def _get_report_name(self):
        rep_pred = 'suggested' if self.params['strict_mode'] else 'suggested,predicted'
        return f'Ranked {rep_pred} Targets and Drugs repositioned for {self._disease2str()}'


    def report_path(self,extension=''):
        return self.data_dir+self._get_report_name()+extension


    def correct_target_ranks_by_confidence(self,drug_urn:str,targetid2urn:dict):
        '''
        Returns
        ------
        corrected_ranks = {target_id:rank}, where target_ids are from targetid2urn\n
        and rank from self.target_ranks was corrected by coeffcients for drug_urn from self.drug2target_consistency 
        '''
        corrected_ranks = {k:v for k,v in self.target_ranks.items() if k in targetid2urn.keys()}
        if self.drug2target_consistency:
            for target_id, rank in corrected_ranks.items():
                target_urn = targetid2urn[target_id]
                try:
                    correction = self.drug2target_consistency[tuple([drug_urn,target_urn])]
                    corrected_ranks[target_id] = correction*rank
                except KeyError:
                    continue

        return corrected_ranks
    

    def __add_rank2(self,drug_df:df,from_drug_graph:ResnetGraph):
        """
Adds following columns to drug_df:\n
DRUG2TARGET_REGULATOR_SCORE,'Directly inhibited targets','Indirectly inhibited targets','Directly activated targets','Indirectly activated targets'
        """
        rank2add2drug = dict()
        for node_id, node in from_drug_graph.nodes(data=True):
            if node['ObjTypeName'][0] == 'SmallMol':
                direct_target_ids, indirect_target_ids = from_drug_graph.direct_indirect_targets(node_id,
                direct_reltypes=['DirectRegulation'],
                indirect_reltypes=['Regulation','MolTransport','Expression']
                )

                if direct_target_ids or indirect_target_ids:
                    rank2add2drug[node['Name'][0]] = node[DRUG2TARGET_REGULATOR_SCORE]

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
DRUG2TARGET_REGULATOR_SCORE,'Directly inhibited targets','Indirectly inhibited targets','Directly activated targets','Indirectly activated targets'
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


    def subtract_antigraph(self,my_df:df, antigraph:ResnetGraph):
        drug_objs = antigraph.get_objects(only_with_values=['SmallMol'])
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


    def addrank2(self,drug_df:df,from_drug_graph:ResnetGraph):
        if self.add_targets2drugs_ws:
            return self.__add_rank_targets2(drug_df,from_drug_graph)
        else:
            return self.__add_rank2(drug_df,from_drug_graph)


    def drug_target_confidence(self,drugs2targets_graph:ResnetGraph,for_agonists:bool):
        '''
        Loads
        -----
        self.drug2target_consistency = {(drug.urn(),target.urn()):1 + consistency_coefficient}\nwhere\n
        consistency_coefficient = (consistency_counter-inconsistency_counter)/\n(consistency_counter+inconsistency_counter)
        '''
        #drug_class = 'agonist' if drug_is_agonist else 'antagonist'
        target_ids = drugs2targets_graph.get_node_ids(self.params['target_types'])
        drugs_ids = drugs2targets_graph.get_node_ids(['SmallMol'])
        print('Calculating drug-target consistency coefficients for %d drugs and %d targets' % (len(drugs_ids),len(target_ids)),flush=True)
        
        if self.clone_all_sessions:
            drug_target_confidence_session = APISession(self.APIconfig,NO_REL_PROPERTIES)
        else:
            drug_target_confidence_session = self

        # drug_target_confidence_session need PS_REFIID_TYPES to execute effect_stats
        drug_target_confidence_session.add_rel_props([EFFECT]+PS_REFIID_TYPES) 
        # do not specify Effect values in common_neighbors because it kills the query
        drug2target2disease_graph = drug_target_confidence_session.common_neighbors_with_effect(
            ['Disease','Virus'],drugs_ids,['Regulation'],'',target_ids,['Regulation','QuantitativeChange'],'')
        
        drug2target2disease_graph.add_graph(drugs2targets_graph)
        drugid2urn = drugs2targets_graph.node_id2urn(['SmallMol'])

        drug2target_consistency_weights = dict()
        for drug_id in drugs_ids:
            drug_urn = drugid2urn[drug_id]
            drug_neighbors_ids = drug2target2disease_graph.neighbors(drug_id)
            neighbors_obj = drug2target2disease_graph._get_nodes(drug_neighbors_ids)

            disease_objs = [n for n in neighbors_obj if n.objtype() in ['Disease','Virus']]
            target_obj = [n for n in neighbors_obj if n.objtype() in self.params['target_types']]
 
            consistency_counter = 0
            inconsistency_counter = 0
            for target in target_obj:
                for disease in disease_objs:      
                    target2disease,pos_refs,neg_refs = drug2target2disease_graph.effect_stats(disease['Id'],target['Id'])
                    if len(pos_refs) > len(neg_refs):
                        target_effect = 1
                    elif len(pos_refs) < len(neg_refs):
                        target_effect = -1
                    else:
                        continue

                    if for_agonists:
                        drug_predicted_effect = target_effect
                    else:
                        drug_predicted_effect = -target_effect
                
                    drug2disease, pos_refs, neg_refs = drug2target2disease_graph.effect_stats(disease['Id'],[drug_id])
                    if len(pos_refs) > len(neg_refs):
                        drug_known_effect = 1
                    elif len(pos_refs) < len(neg_refs):
                        drug_known_effect = -1
                    else:
                        continue

                    if drug_predicted_effect == drug_known_effect:
                        consistency_counter += 1
                    else:
                        inconsistency_counter += 1

            consistency_sum = consistency_counter+inconsistency_counter
            if consistency_sum:
                consistency_coefficient = float(consistency_counter-inconsistency_counter)/(consistency_sum)
            else:
                consistency_coefficient = 0.0

            drug2target_consistency_weights[(drug_urn,target.urn())] = 1 + consistency_coefficient
        
        return self.drug2target_consistency.update(drug2target_consistency_weights)


    '''
       def addrank2_OLD(self,drug_df:df, from_drug_graph:ResnetGraph):
        """
        Adds
        ----
        columns: DRUG2TARGET_REGULATOR_SCORE,'Effect on targets','Direct targets','Indirect targets'
        """
        drug2rows = dict()
        for node_id, node in from_drug_graph.nodes(data=True):
            if node['ObjTypeName'][0] == 'SmallMol':
                direct_target_ids = from_drug_graph.downstream_targets(node_id,['DirectRegulation'])
                if direct_target_ids:
                    direct_target_names = [name[0] for t,name in from_drug_graph.nodes(data='Name') if t in direct_target_ids]
                    direct_target_names_str = ','.join(direct_target_names)
                    [self.direct_target2drugs.update_with_value(n,node['Name'][0]) for n in direct_target_names]
                else:
                    direct_target_names_str = ''
                indirect_target_ids = from_drug_graph.downstream_targets(node_id,['Regulation','MolTransport','Expression'])
                if indirect_target_ids:
                    indirect_target_names = [name[0] for t,name in from_drug_graph.nodes(data='Name') if t in indirect_target_ids]
                    indirect_target_names_str = ','.join(indirect_target_names)
                    [self.indirect_target2drugs.update_with_value(n,node['Name'][0]) for n in indirect_target_names]
                else:
                    indirect_target_names_str = ''
                    
                if direct_target_names_str or indirect_target_names_str:
                    drug2rows[node['Name'][0]] = [node[DRUG2TARGET_REGULATOR_SCORE],node['Drug Class'][0],direct_target_names_str,indirect_target_names_str]
        df_copy = df.copy_df(drug_df)
        if DRUG2TARGET_REGULATOR_SCORE not in list(df_copy.columns):
            df_copy[DRUG2TARGET_REGULATOR_SCORE] = np.NaN
            df_copy['Effect on targets'] = NaN
            df_copy['Direct targets'] = NaN
            df_copy['Indirect targets'] = NaN
        annotate_df = df(df_copy[df_copy[DRUG2TARGET_REGULATOR_SCORE].isnull()])
        # saving df that already has values from previous graphs:
        has_scores_df = df(df_copy[df_copy[DRUG2TARGET_REGULATOR_SCORE].notnull()])
        for idx in annotate_df.index:
            drug_name = annotate_df.at[idx,'Name']
            try:
                rows2add = drug2rows[drug_name]
                annotate_df.loc[idx,DRUG2TARGET_REGULATOR_SCORE] = rows2add[0]
                annotate_df.loc[idx,'Effect on targets'] = rows2add[1]
                annotate_df.loc[idx,'Direct targets'] = rows2add[2]
                annotate_df.loc[idx,'Indirect targets'] = rows2add[3]
                drug2rows.pop(drug_name)
            except KeyError:
                continue
        if drug2rows: # there are drug that both agonists and antagonists
            duplicate_rows = df(columns=has_scores_df.columns)
            for idx in has_scores_df.index:
                drug_name = has_scores_df.at[idx,'Name']
                try:
                    rows2add = drug2rows[drug_name]
                    my_rank = has_scores_df.at[idx,DRUG2TARGET_REGULATOR_SCORE]
                    add_rank = rows2add[0]
                    new_rank = my_rank+add_rank
                    has_scores_df.at[idx,DRUG2TARGET_REGULATOR_SCORE] = new_rank
                    last_row = len(duplicate_rows)
                    old_row = has_scores_df.loc[idx].to_list()
                    duplicate_rows.loc[last_row] = old_row
                    #duplicate_rows.loc[last_row,DRUG2TARGET_REGULATOR_SCORE] = new_rank
                    duplicate_rows.loc[last_row,'Effect on targets'] = rows2add[1]
                    duplicate_rows.loc[last_row,'Direct targets'] = rows2add[2]
                    duplicate_rows.loc[last_row,'Indirect targets'] = rows2add[3]
                except KeyError:
                    continue
            df2return = df(pd.concat([annotate_df,has_scores_df,duplicate_rows],ignore_index=True))
        else:
            df2return = df(pd.concat([annotate_df,has_scores_df],ignore_index=True))
        df2return._name_ = drug_df._name_
        return df2return
    def filter4drugs(self,compound_ids:list):
        """
        Returns
        -------
        set(PSObjects) with compound_ids that also belong to "drug" PS Ontology and has "Pharmapendium ID"
        """
        if self.clone:
            my_api_session = APISession(self.APIconfig,NO_REL_PROPERTIES)
            my_api_session.add_ent_props(["PharmaPendium ID"])
        else:
            my_api_session = self
        oql_query = "SELECT Entity WHERE id = ({ids}) AND InOntology (SELECT Annotation WHERE Ontology='Pathway Studio Ontology' \
            AND Relationship='is-a') under (SELECT OntologicalNode WHERE Name = drugs)"
        r_n = f'Find drugs linked to {self._disease2str()}'
        drugs_graph = my_api_session.iterate_oql(oql_query,compound_ids,use_cache=False,request_name=r_n)
        if self.drugs_must_have_props:
            drugs_graph = drugs_graph.subgraph_by_nodeprops(self.drugs_must_have_props)
        
        drugs = drugs_graph.get_objects(only_with_values=['SmallMol'])
        print ('Found %d drugs among %d compounds' % (len(drugs),len(compound_ids)))
        return set(drugs)
    '''


    '''
    def find_rank_drugs4(self,target_ids:list, with_effect:str,correct4consistency=True):
        """
        DEPRICATED
        ---------
        Input
        -----
        correct4consistency must be supplied explicitly to avoid consitency correction for toxicities in antigraph
        Returns
        -------
        ResnetGraph with drugs annotated with DRUG2TARGET_REGULATOR_SCORE and drug class (agonist/antagonist)
        """
        drug_is_agonist = False if with_effect == 'negative' else True
        if self.clone:
            my_api_session = APISession(self.APIconfig,NO_REL_PROPERTIES)
            my_api_session.add_rel_props([EFFECT])
            my_api_session.Graph = self.Graph.copy() #should speed up drug-target retrieval
        else:
            my_api_session = self
        drug_ids = {d.id() for d in self.drugs} #needs self.filter4drugs to setup self.drugs
        select_drugs = 'SELECT Entity WHERE id = ({ids1})'
        oql_query = f'SELECT Relation WHERE Effect = {with_effect} AND NeighborOf downstream ({select_drugs})'
        oql_query = oql_query + ' AND NeighborOf upstream (SELECT Entity WHERE id = ({ids2}))'
        
        r_n = f'Find drugs {with_effect}ly regulating {str(len(target_ids))} ranked targets'
        drug_graph = my_api_session.iterate_oql2(oql_query,drug_ids,target_ids,request_name=r_n)
        if not drug_graph: return ResnetGraph()
            
        direct_reg_graph = drug_graph.subgraph_by_relprops(['DirectRegulation'])
        indirect_reg_graph = drug_graph.subgraph_by_relprops(['Regulation','Expression','MolTransport'])
        if correct4consistency:
        # correct4consistency must be supplied explicitly to avoid consitency correction for toxicities in antigraph
            drug_class = 'Agonist' if drug_is_agonist else 'Antagonist'
            t1 = Thread(target=self.drug_target_confidence, args=(direct_reg_graph,drug_is_agonist),name=drug_class+'DirectTargetsConfidence')
            t1.start()
            t2 = Thread(target=self.drug_target_confidence, args=(indirect_reg_graph,drug_is_agonist),name=drug_class+'IndirectTargetsConfidence')
            t2.start()
            t1.join()
            t2.join()
        
        # initializing drug ranks
        drug2rank = dict()
        for drug in direct_reg_graph.get_objects(['SmallMol']):
            target_ranks4drug = self.correct_ranks_by_confidence(drug.urn(),direct_reg_graph.node_id2urn(self.params['target_types']))
            drug2rank[drug.id()] = 2 * direct_reg_graph.rank_regulator(drug.id(),target_ranks4drug)
        # adding indirect targets to ranking
        for drug in indirect_reg_graph.get_objects(['SmallMol']):
            target_ranks4drug = self.correct_ranks_by_confidence(drug.urn(),indirect_reg_graph.node_id2urn(self.params['target_types']))
            try:
                my_rank = drug2rank[drug.id()]
            except KeyError:
                my_rank = 0.0
            drug2rank[drug.id()] = my_rank + indirect_reg_graph.rank_regulator(drug.id(),target_ranks4drug)
        all_reg_graph = direct_reg_graph.copy()
        all_reg_graph.add_graph(indirect_reg_graph)
        
        drug_class = 'agonist' if with_effect == 'positive' else 'antagonist'
        drug_id2urn = all_reg_graph.node_id2urn(['SmallMol'])
        drug_urn2class = {urn:[drug_class] for urn in drug_id2urn.values()}
        all_reg_graph.set_node_annotation(drug_urn2class,'Drug Class')
        nx.set_node_attributes(all_reg_graph,drug2rank,DRUG2TARGET_REGULATOR_SCORE)
        return all_reg_graph
    '''


    def rank_drugs4(self,target_ids:list, with_effect:str,correct4consistency=True):
        """
        Input
        -----
        correct4consistency must be supplied explicitly to avoid consitency correction for toxicities in antigraph
        
        Returns
        -------
        ResnetGraph with drugs annotated with DRUG2TARGET_REGULATOR_SCORE and drug class (agonist/antagonist)
        """
        direct_effect_graph = self.drug2targets_graph.get_neighbors_graph(target_ids,[],['DirectRegulation'],[with_effect])
        indirect_effect_graph = self.drug2targets_graph.get_neighbors_graph(target_ids,[],['Regulation','Expression','MolTransport'],[with_effect])

        drug_class = 'agonist' if with_effect == 'positive' else 'antagonist'
        if correct4consistency:
        # correct4consistency must be supplied explicitly to avoid consistency correction for toxicities in antigraph
            ranking_agonists = False if with_effect == 'negative' else True
            #drug_class = 'Agonist' if ranking_agonists else 'Antagonist'
            t1 = Thread(target=self.drug_target_confidence, args=(direct_effect_graph,ranking_agonists),name=drug_class+'DirectTargetsConfidence')
            t1.start()
            t2 = Thread(target=self.drug_target_confidence, args=(indirect_effect_graph,ranking_agonists),name=drug_class+'IndirectTargetsConfidence')
            t2.start()
            t1.join()
            t2.join()
        
        # initializing drug ranks
        drug2rank = dict()
        direct_targetid2urn = direct_effect_graph.node_id2urn(self.params['target_types'])
        for drug in direct_effect_graph.get_objects(only_with_values=['SmallMol']):
            target_ranks4drug = self.correct_target_ranks_by_confidence(drug.urn(),direct_targetid2urn)
            drug2rank[drug.id()] = 2 * direct_effect_graph.rank_regulator(drug.id(),target_ranks4drug)

        # adding indirect targets to ranking with smaller weight
        indirect_targetid2urn = indirect_effect_graph.node_id2urn(self.params['target_types'])
        for drug in indirect_effect_graph.get_objects(only_with_values=['SmallMol']):
            target_ranks4drug = self.correct_target_ranks_by_confidence(drug.urn(),indirect_targetid2urn)
            try: my_rank = drug2rank[drug.id()]
            except KeyError: my_rank = 0.0
            drug2rank[drug.id()] = my_rank + indirect_effect_graph.rank_regulator(drug.id(),target_ranks4drug)

        all_reg_graph = direct_effect_graph.copy()
        all_reg_graph.add_graph(indirect_effect_graph)
        
        drug_id2urn = all_reg_graph.node_id2urn(['SmallMol'])
        drug_urn2class = {urn:[drug_class] for urn in drug_id2urn.values()}
        all_reg_graph.set_node_annotation(drug_urn2class,'Drug Class')
        nx.set_node_attributes(all_reg_graph,drug2rank,DRUG2TARGET_REGULATOR_SCORE)

        return all_reg_graph


    def load_target_ranks(self,from_ranked_targets_df:df,for_antagonists=True):
        '''
        Input
        -----
        from_ranked_targets_df must have self.__temp_id_col__ and RANK columns
        Loads
        -----
        self.target_ranks,self.need_antagonists_ids,self.need_agonists_ids
        '''
        for idx in from_ranked_targets_df.index:
            target_ids = from_ranked_targets_df.at[idx,self.__temp_id_col__]
            target_rank = from_ranked_targets_df.at[idx,RANK]
            for target_id in target_ids:
                self.target_ranks[target_id] = float(target_rank)

        if for_antagonists:
            self.targets4antagonists_ids = self._all_ids(from_ranked_targets_df)
        else:
            self.targets4agonists_ids = self._all_ids(from_ranked_targets_df)


    def load_drug_graph(self,for_target_ids:list,only_drugs_with_ids:list=[]):
        '''
        Input
        -----
        self.targets4agonists_ids, self.targets4antagonists_ids
        Output
        ------
        self.drug2targets_graph used by self.regulatory_rank()\n
        self.drugs used by self.init_drug_df()
        '''
        if hasattr(self, "drug2targets_graph"): return # for SNEA to avoid multiple self.drug2targets_graph loading
        print('Loading drug-target graph')
        start_time = time.time()
        
        my_api_session = APISession(self.APIconfig,NO_REL_PROPERTIES) if self.clone_all_sessions else self
        my_api_session.add_rel_props([EFFECT])
        my_api_session.entProps = ['Name','PharmaPendium ID']

        request_name = f'Find drugs for {len(for_target_ids)} targets'
        # use_cache must be False to ensure retrieval of "PharmaPendium ID" if self.clone_all_sessions = False:
        if only_drugs_with_ids:
            get_targets = 'SELECT Entity WHERE id = ({ids2})'
            select_drugs = 'SELECT Entity WHERE id =({ids1})'
            OQLquery = f'SELECT Relation WHERE NeighborOf downstream ({select_drugs}) AND NeighborOf upstream ({get_targets})'
            self.drug2targets_graph = my_api_session.iterate_oql2(OQLquery,only_drugs_with_ids,for_target_ids,use_cache=False,request_name=request_name)

            select_drug_binding = f'SELECT Relation WHERE objectType = Binding AND NeighborOf ({select_drugs}) AND NeighborOf ({get_targets})'
            request_name = f'Find drugs binding {len(for_target_ids)} targets'
            drug_binding_graph = my_api_session.iterate_oql2(select_drug_binding,only_drugs_with_ids,for_target_ids,use_cache=False,request_name=request_name)
        else:
            get_targets = 'SELECT Entity WHERE id = ({ids})'
            select_drugs = OQL.select_drugs()
            OQLquery = f'SELECT Relation WHERE NeighborOf downstream ({select_drugs}) AND NeighborOf upstream ({get_targets})'
            self.drug2targets_graph = my_api_session.iterate_oql(OQLquery,for_target_ids,use_cache=False,request_name=request_name)

            select_drug_binding = f'SELECT Relation WHERE objectType = Binding AND NeighborOf ({select_drugs}) AND NeighborOf ({get_targets})'
            request_name = f'Find drugs binding {len(for_target_ids)} targets'
            drug_binding_graph = my_api_session.iterate_oql(select_drug_binding,for_target_ids,use_cache=False,request_name=request_name)
            
        self.drug2targets_graph.add_graph(drug_binding_graph)

        reltype_ranking = ['DirectRegulation','Binding','Expression','MolTransport','Regulation']
        self.drug2targets_graph =  self.drug2targets_graph.make_simple(reltype_ranking)
        self.drug2targets_graph = self.drug2targets_graph.subgraph_by_relprops(['positive','negative'],[EFFECT])
        

        drugs4targets_ids = self.drug2targets_graph.get_node_ids(['SmallMol'])
        self.drugs.update(self.drug2targets_graph._get_nodes(drugs4targets_ids))
        
        # self.drug_ids is used for target ranking by parent DiseaseTargets class
        # self.drugs is used by Drugs4Targets to init drug_df with drug names
        # self.drugs holds different set of drugs from self.drug_ids
        if self.drug_ids:
            oql_query = "SELECT Entity WHERE id = ({ids}) AND InOntology (SELECT Annotation WHERE Ontology='Pathway Studio Ontology' \
            AND Relationship='is-a') under (SELECT OntologicalNode WHERE Name = drugs)"
            r_n = f'Find drugs linked to {self._disease2str()}'
            drugs_withno_targets = self.drug_ids.difference(drugs4targets_ids)
            alldrugs4diseases = my_api_session.iterate_oql(oql_query,drugs_withno_targets,request_name=r_n)
            #drugs4disease = alldrugs4diseases.subgraph_by_nodeprops(["PharmaPendium ID"]) # algorithm fetches all drugs now
            len_before_add = len(self.drugs)
            drugs4disease_with_no_targets = alldrugs4diseases._get_nodes()
            self.drugs.update(drugs4disease_with_no_targets)
            added_count = len(self.drugs)-len_before_add
            print('Found additional %d drugs without targets linked to %s' % (added_count,self._disease2str()))
        
        print(f'Found {str(len(drugs4targets_ids))} drugs linked to {str(len(for_target_ids))} targets for ranking')
        print('Drug-target graph was loaded in %s' % self.execution_time(start_time))
        my_api_session.entProps.remove("PharmaPendium ID")
        return


    def init_drug_df(self, drugs:list):
        '''
        Input
        -----
        drugs = [PSObject]
        
        Returns
        -------
        new self.RefCountPandas
        '''
        print('Loading drug ranking worksheet')
        drug_names = list()
        drug_name2ppid = dict()
        for drug in drugs:
            drug_name = drug.name()
            if 'antigen' in drug_name: continue # bad fix for PS ontology problem
            if 'vaccine' in drug_name: continue
            if drug_name == 'DMSO': continue
            drug_names.append(drug_name)
            try:
                drug_name2ppid[drug_name] = drug['PharmaPendium ID'][0]
            except KeyError:
               continue

        drug_names_df = df.from_dict({'Name':list(drug_names)})
        self.entProps = ['Name']
        drug_df = self.load_pandas(drug_names_df,map2type=['SmallMol'],max_children_count=11)
        self.RefCountPandas = drug_df.merge_dict(drug_name2ppid,new_col=PHARMAPENDIUM_ID,map2column='Name')


    def regulatory_rank(self):
        '''
        Returns
        -------
        subset of self.RefCountPandas with drugs having positive regulatory score\n
        returned df has new column DRUG2TARGET_REGULATOR_SCORE
        '''
        antagonist_antigraph = self.rank_drugs4(self.targets4antagonists_ids,'positive', correct4consistency=False)
        agonist_antigraph = self.rank_drugs4(self.targets4agonists_ids,'negative',correct4consistency=False)
        
        with ThreadPoolExecutor(max_workers=4, thread_name_prefix='DrugsRegulatorRanking') as executor:
            future1 = executor.submit(self.rank_drugs4,self.targets4antagonists_ids,'negative',self.params['consistency_correction4target_rank'])
            future3 = executor.submit(self.rank_drugs4,self.targets4agonists_ids,'positive',self.params['consistency_correction4target_rank']) 
            antagonist_graph = future1.result()           
            agonist_graph = future3.result()

        new_df = self.addrank2(self.RefCountPandas,antagonist_graph)
        new_df = self.addrank2(new_df,agonist_graph)
        new_df = self.subtract_antigraph(new_df,antagonist_antigraph)
        new_df = self.subtract_antigraph(new_df,agonist_antigraph)

        new_df = new_df.greater_than(0,DRUG2TARGET_REGULATOR_SCORE)
        new_df._name_ = f'Drugs for {self.input_names()}'
        print('Found %d drugs for %d antagonist targets and %d agonist targets for "%s"' % 
            (len(new_df),len(self.targets4antagonists_ids),len(self.targets4agonists_ids),self.input_names()),flush=True)
        return new_df
        

    def link2disease_concepts(self,in_df:df):
        if not self.disease_ids: return in_df # for SNEA.make_drugs_df()
        colname = self._disease2str()+' clinical trials'
        self.set_how2connect(['ClinicalTrial'],[],'',how2clone=REFERENCE_IDENTIFIERS)
        linked_row_count,linked_ids,drug_df = self.link2concept(colname,self.disease_ids,in_df)
        print('%d drugs on clinical trials for %s' % (linked_row_count,self._disease2str()),flush=True)
        self.column_ranks[5] = list(drug_df.columns)[-1]

        print('Linking %d drugs with %s' % (len(self.drugs),self._disease2str()),flush=True)
        colname = 'regulation of '+ self._disease2str()
        self.set_how2connect(['Regulation'],[],'',how2clone=REFERENCE_IDENTIFIERS)
        linked_row_count,linked_ids,drug_df = self.link2concept(colname,self.disease_ids,drug_df)
        print('%d drugs regulating %s' % (linked_row_count,self._disease2str()),flush=True)
        self.column_ranks[4] = list(drug_df.columns)[-1]

        print('Linking %d drugs with %d Symptoms linked to %s' % (len(self.drugs),len(self.symptoms_ids),self._disease2str()),flush=True)
        colname = 'symptoms for '+ self._disease2str()
        self.set_how2connect(['Regulation'],[],'',how2clone=REFERENCE_IDENTIFIERS)
        linked_row_count,linked_ids,drug_df = self.link2concept(colname,self.symptoms_ids,drug_df)
        print('%d drugs linked to symptoms for %s' % (linked_row_count,self._disease2str()),flush=True)
        self.column_ranks[1] = list(drug_df.columns)[-1]

        print('Linking %d drugs with %d ClinicalParameters linked to %s' % (len(self.drugs),len(self.clinpar_ids),self._disease2str()),flush=True)
        colname = 'Clinical parameters for '+ self._disease2str()
        self.set_how2connect(['Regulation'],[],'',how2clone=REFERENCE_IDENTIFIERS)
        linked_row_count,linked_ids,drug_df = self.link2concept(colname,self.clinpar_ids,drug_df)
        print('%d drugs linked to clnical parameters for %s' % (linked_row_count,self._disease2str()),flush=True)
        self.column_ranks[2] = list(drug_df.columns)[-1]

        print('Linking %d drugs with %d CellProcess linked to %s' % (len(self.drugs),len(self.cellproc_ids),self._disease2str()),flush=True)
        colname = 'Cell processess affected by '+ self._disease2str()
        self.set_how2connect(['Regulation'],[],'',how2clone=REFERENCE_IDENTIFIERS)
        linked_row_count,linked_ids,drug_df = self.link2concept(colname,self.cellproc_ids,drug_df)
        print('%d drugs linked to cell processes affected in %s' % (linked_row_count,self._disease2str()),flush=True)
        self.column_ranks[3] = list(drug_df.columns)[-1]

        return drug_df


    def score_drugs(self,normalize=True):
        # making list of drugs
        start_time = time.time()
        self.load_target_ranks(self.report_pandas[AGONIST_TARGETS_WS],for_antagonists=False)
        #if targets are in both worksheets will assume that it needs antagonist:
        self.load_target_ranks(self.report_pandas[ANTAGONIST_TARGETS_WS]) 
        ranked_df = self.regulatory_rank()
        self.column_ranks = {0:DRUG2TARGET_REGULATOR_SCORE}

        if self.params['add_bibliography']:
            with ThreadPoolExecutor(max_workers=1,thread_name_prefix='ScoreDrugs') as sd:
                drugs_etm_thread = sd.submit(self.etm_refs2df,ranked_df,self.input_names())
                concepts_thread = sd.submit(self.link2disease_concepts,ranked_df)
                full_drug_df = concepts_thread.result()
                drugs_df_with_etmrefs = drugs_etm_thread.result()
                #merging results of two threads:
                ref_column_name = self.etm_counter.etm_ref_column_name[-1]
                full_drug_df[ref_column_name] = drugs_df_with_etmrefs[ref_column_name]
                doi_column_name = self.etm_counter.etm_doi_column_name[-1]
                full_drug_df[doi_column_name] = drugs_df_with_etmrefs[doi_column_name]
                full_drug_df.add_format(from_df=drugs_df_with_etmrefs)
        else:
            full_drug_df = self.link2disease_concepts(ranked_df)

        raw_drug_df = self.make_count_df(full_drug_df,'rawDrugs')
        self.add2raw(raw_drug_df)

        if normalize:
            ranks = sorted(self.column_ranks)
            columns2norm = [self.column_ranks[k] for k in ranks]
            self.normalize('rawDrugs','Drugs','Name',columns2norm)

            #moving RANK column to second position
            drugdf_columns = list(self.report_pandas['Drugs'].columns)
            drugdf_columns.remove(RANK)
            drugdf_columns.insert(1,RANK)
            drugdf_columns.remove(PHARMAPENDIUM_ID)
            drugdf_columns = [PHARMAPENDIUM_ID]+drugdf_columns
            self.report_pandas['Drugs'] =  self.report_pandas['Drugs'].reorder(drugdf_columns)
            
            print("Drug ranking was done in %s" % self.execution_time(start_time), flush=True)
            print('Normalized worksheet named "Drugs" was added to report')
            return self.report_pandas['Drugs']
        else:
            print('Ranked drugs are in worksheet "rawDrugs" in raw data')
            return self.raw_data['rawDrugs']


    def init_load_score(self,for_target_ids:list):
        self.load_drug_graph(for_target_ids)
        self.init_drug_df(self.drugs)
        self.score_drugs()


    def make_report(self):
        start_time = time.time()
        super().make_report()
        self.report_pandas[ANTAGONIST_TARGETS_WS] = self.remove_high_level_entities(self.report_pandas[ANTAGONIST_TARGETS_WS])
        self.report_pandas[AGONIST_TARGETS_WS] = self.remove_high_level_entities(self.report_pandas[AGONIST_TARGETS_WS])
        all_target_ids = list(self._all_ids(self.report_pandas[ANTAGONIST_TARGETS_WS])|self._all_ids(self.report_pandas[AGONIST_TARGETS_WS]))
        
        if self.params['add_bibliography']:
            with ThreadPoolExecutor(max_workers=4, thread_name_prefix='DrugsRanking') as executor:
                future1 = executor.submit(self.init_load_score,all_target_ids)
                future3 = executor.submit(self.add_etm_bibliography) 
                future1.result()           
                future3.result()
        else:
            self.init_load_score(all_target_ids)
        
        self.report_pandas[ANTAGONIST_TARGETS_WS] = self.report_pandas[ANTAGONIST_TARGETS_WS].merge_psobject(self.direct_target2drugs,'Directly Inhibited by','Name',values21cell=True)
        self.report_pandas[ANTAGONIST_TARGETS_WS] = self.report_pandas[ANTAGONIST_TARGETS_WS].merge_psobject(self.indirect_target2drugs,'Indirectly Inhibited by','Name',values21cell=True)
        
        self.report_pandas[AGONIST_TARGETS_WS] =  self.report_pandas[AGONIST_TARGETS_WS].merge_psobject(self.direct_target2drugs,'Directly Activated by','Name',values21cell=True)
        self.report_pandas[AGONIST_TARGETS_WS] = self.report_pandas[AGONIST_TARGETS_WS].merge_psobject(self.indirect_target2drugs,'Indirectly Activated by','Name',values21cell=True)
        print('Drug repurposing for %s was done in %s' % (self._disease2str(), self.execution_time(start_time)))
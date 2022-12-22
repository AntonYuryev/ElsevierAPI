from .DiseaseTargets import DiseaseTargets,ResnetGraph,EFFECT
from .DiseaseTargets import ANTAGONIST_TARGETS_WS,AGONIST_TARGETS_WS,REFERENCE_IDENTIFIERS
from .SemanticSearch import RANK
from ..ETM_API.etm import ETM_REFS_COLUMN
from .ResnetAPISession import APISession, NO_REL_PROPERTIES
from ..pandas.panda_tricks import df,pd,PSObject
from numpy import NaN
import time
import networkx as nx
from threading import Thread
from concurrent.futures import ThreadPoolExecutor
from threading import Thread

PATHWAY_REGULATOR_SCORE = 'Disease model regulator score'
DRUG2TARGET_REGULATOR_SCORE = 'Targets regulator score'

class Drugs4Targets(DiseaseTargets):
    pass
    drug2target_consistency = dict()
    target_ranks = dict()
    antagonist_graph = ResnetGraph()
    antagonist_antigraph = ResnetGraph()
    agonist_graph = ResnetGraph()
    agonist_antigraph = ResnetGraph()
    column_ranks = dict()
    direct_target2drugs = PSObject()
    indirect_target2drugs = PSObject()

    def __init__(self, APIconfig, params=dict()):
        """
        self.set_target_disease_state() needs references
        """
        super().__init__(APIconfig)
        if params:
            self.param = dict(params)
            self.set_dir(self.param['data_dir'])
        else:
            self.param = {
                'disease':[],
                'strict_mode':True,
                'data_dir':'',
                'add_bibliography' : False,
                'strict_mode' : False,
                'target_types' : ['Protein'],
                'pathway_folders':[],
                'pathways': [],
                'consistency_correction4target_rank':False
                }

            
    def _get_report_name(self):
        return f'Ranked Targets and Drugs repositioned for {self._disease2str()}'


    def report_path(self,extension=''):
        return self.data_dir+self._get_report_name()+extension


    def rank2dict(self,from_ranked_targets_df:df):
        """
        Returns
        -------
        {target_id:target_rank}
        """
        for idx in from_ranked_targets_df.index:
            target_ids = from_ranked_targets_df.at[idx,self.__temp_id_col__]
            target_rank = from_ranked_targets_df.at[idx,RANK]
            for target_id in target_ids:
                self.target_ranks[target_id] = float(target_rank)


    def correct_ranks_by_confidence(self,drug_urn:str, targetid2urn:dict):
        corrected_ranks = dict()
        for target_id, rank in self.target_ranks.items():
            try:
                target_urn = targetid2urn[target_id]
                key = tuple([drug_urn,target_urn])
                try:
                    correction = self.drug2target_consistency[key]
                    corrected_ranks[target_id] = correction*rank
                except KeyError:
                    corrected_ranks[target_id] = rank
            except KeyError:
                continue
            
        return corrected_ranks


    def addrank2(self,drug_df:df, from_drug_graph:ResnetGraph):
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
            df_copy[DRUG2TARGET_REGULATOR_SCORE] = NaN
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


    def subtract_antigraph(self,my_df:df, antigraph:ResnetGraph):
        drug_objs = antigraph.get_objects(['SmallMol'])
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


    def drug_target_confidence(self,drugs2targets_graph:ResnetGraph,drug_is_agonist=False):
        drug_class = 'agonist' if drug_is_agonist else 'antagonist'
        target_ids = drugs2targets_graph.get_node_ids(self.param['target_types'])
        drugs_ids = drugs2targets_graph.get_node_ids(['SmallMol'])
        print('Calculating drug-target consistency coefficients for %d %ss and %d targets' % (len(drugs_ids),drug_class,len(target_ids)),flush=True)
        drug_target_confidence_session = APISession(self.APIconfig,NO_REL_PROPERTIES)
        drug_target_confidence_session.add_rel_props([EFFECT])
        #drug_target_confidence_session.Graph = self.Graph.copy() # should speed up drug-target-disease retreival
        drug2target2disease_graph = drug_target_confidence_session.common_neighbors(['Disease','Virus'],drugs_ids,['Regulation','ClinicalTrial'],[],'',
                                target_ids,['Regulation','QuantitativeChange'],[],'')
        
        drug2target2disease_graph.add_graph(drugs2targets_graph)
        drugid2urn = drugs2targets_graph.node_id2urn(['SmallMol'])

        drug2target_consistency_weights = dict()
        for drug_id in drugs_ids:
            drug_urn = drugid2urn[drug_id]
            drug_neighbors_ids = drug2target2disease_graph.neighbors(drug_id)
            neighbors_obj = drug2target2disease_graph._get_nodes(drug_neighbors_ids)

            disease_objs = [n for n in neighbors_obj if n.objtype() in ['Disease','Virus']]
            target_obj = [n for n in neighbors_obj if n.objtype() in self.param['target_types']]
 
            consistency_counter = 0
            inconsistency_counter = 0
            for disease in disease_objs:
                for target in target_obj:
                    target2disease,pos_refs,neg_refs = drug2target2disease_graph.effect_stats(disease['Id'],target['Id'])
                    if len(pos_refs) > len(neg_refs):
                        target_effect = 1
                    elif len(pos_refs) < len(neg_refs):
                        target_effect = -1
                    else:
                        continue

                    if drug_is_agonist:
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


    def filter4drugs(self,compound_ids:list):
        """
        Returns
        -------
        set(PSObjects) with compound_ids that also belong to "drug" PS Ontology and has "Pharmapendium ID"
        """
        new_api_session = APISession(self.APIconfig,NO_REL_PROPERTIES)
        new_api_session.add_ent_props(["PharmaPendium ID"])
        oql_query = "SELECT Entity WHERE id = ({ids}) AND InOntology (SELECT Annotation WHERE Ontology='Pathway Studio Ontology' \
            AND Relationship='is-a') under (SELECT OntologicalNode WHERE Name = drugs)"
        r_n = f'Find drugs linked to {self._disease2str()}'
        # do not use_relation_cache to retrieve "PharmaPendium ID"
        drugs_graph = new_api_session.iterate_oql(oql_query,compound_ids,use_relation_cache=False,request_name=r_n)
        
        maybe_drugs = drugs_graph.get_objects(['SmallMol'])
        drugs = set()
        for drug in maybe_drugs:
            try:
                pp_ids = drug["PharmaPendium ID"]
                drugs.add(drug)
            except KeyError:
                continue
        
        print ('Found %d drugs among %d compounds' % (len(drugs),len(compound_ids)))
        return drugs

    
    def find_rank_drugs4(self,target_ids:list, with_effect:str,correct4consistency=True):
        """
        Returns
        -------
        ResnetGraph with drugs annotated with DRUG2TARGET_REGULATOR_SCORE and drug class (agonist/antagonist)
        """
        drug_is_agonist = False if with_effect == 'negative' else True
        
        my_api_session = APISession(self.APIconfig,NO_REL_PROPERTIES)
        my_api_session.add_rel_props([EFFECT])
        my_api_session.Graph = self.Graph.copy() #should speed up drug-target retrieval

        select_drugs = 'SELECT Entity WHERE objectType = SmallMol AND NOT ("PharmaPendium ID" = NULL)'
        oql_query = f'SELECT Relation WHERE Effect = {with_effect} AND NeighborOf downstream ({select_drugs})'
        oql_query = oql_query + ' AND NeighborOf upstream (SELECT Entity WHERE id = ({ids}))'
        
        r_n = f'Find drugs {with_effect}ly regulating {str(len(target_ids))} ranked targets'
        drug_graph = my_api_session.iterate_oql(oql_query,target_ids, request_name=r_n)

        compound_ids = drug_graph.get_node_ids(['SmallMol'])
        drugs = self.filter4drugs(compound_ids)
        drug_ids = {d.id() for d in drugs}
        remove_ids = set(compound_ids).difference(drug_ids)
        drug_graph.remove_nodes_from(remove_ids)

        if not drug_graph: return ResnetGraph()
            
        direct_reg_graph = drug_graph.subgraph_by_relprops(['DirectRegulation'])
        indirect_reg_graph = drug_graph.subgraph_by_relprops(['Regulation','Expression','MolTransport'])

        if correct4consistency:
            drug_class = 'Agonist' if drug_is_agonist else 'Antagonist'
            t1 = Thread(target=self.drug_target_confidence, args=(direct_reg_graph,drug_is_agonist),name=drug_class+'DirectTargetsConfidence')
            t1.start()
            t2 = Thread(target=self.drug_target_confidence, args=(indirect_reg_graph,drug_is_agonist),name=drug_class+'IndirectTargetsConfidence')
            t2.start()
            t1.join()
            t2.join()
            
        drug_scores = dict()
        for drug in direct_reg_graph.get_objects(['SmallMol']):
            target_ranks4drug = self.correct_ranks_by_confidence(drug.urn(),direct_reg_graph.node_id2urn(self.param['target_types']))
            drug_scores[drug.id()] = 2 * direct_reg_graph.rank_regulator(drug.id(),target_ranks4drug)

        for drug in indirect_reg_graph.get_objects(['SmallMol']):
            target_ranks4drug = self.correct_ranks_by_confidence(drug.urn(),indirect_reg_graph.node_id2urn(self.param['target_types']))
            try:
                my_rank = drug_scores[drug.id()]
            except KeyError:
                my_rank = 0.0
            drug_scores[drug.id()] = my_rank + indirect_reg_graph.rank_regulator(drug.id(),target_ranks4drug)

        all_reg_graph = direct_reg_graph.copy()
        all_reg_graph.add_graph(indirect_reg_graph)
        
        drug_class = 'agonist' if with_effect == 'positive' else 'antagonist'
        drug_id2urn = all_reg_graph.node_id2urn(['SmallMol'])
        drug_urn2class = {urn:[drug_class] for urn in drug_id2urn.values()}
        all_reg_graph.set_node_annotation(drug_urn2class,'Drug Class')
        nx.set_node_attributes(all_reg_graph,drug_scores,DRUG2TARGET_REGULATOR_SCORE)

        return all_reg_graph


    def init_drug_df(self):
        print('Finding drugs for ranked targets', flush=True)
        self.drugs = self.filter4drugs(self.drug_ids)
        
        need_antagonists_ids = self._all_ids(self.report_pandas[ANTAGONIST_TARGETS_WS])
        self.rank2dict(self.report_pandas[ANTAGONIST_TARGETS_WS])
        need_agonists_ids = self._all_ids(self.report_pandas[AGONIST_TARGETS_WS])
        self.rank2dict(self.report_pandas[AGONIST_TARGETS_WS])
        
        # t1 = Thread(target=self.drug_target_confidence, args=(direct_reg_graph,drug_is_agonist))
        with ThreadPoolExecutor(max_workers=4, thread_name_prefix='DrugsRegulatorRanking') as executor:
            future1 = executor.submit(self.find_rank_drugs4,need_antagonists_ids,'negative',self.param['consistency_correction4target_rank'])
            future2 = executor.submit(self.find_rank_drugs4,need_antagonists_ids,'positive', False,)
            future3 = executor.submit(self.find_rank_drugs4,need_agonists_ids,'positive',self.param['consistency_correction4target_rank'])
            future4 = executor.submit(self.find_rank_drugs4,need_agonists_ids,'negative',False)
        
            antagonist_graph = future1.result()
            antagonist_antigraph = future2.result()
            agonist_graph = future3.result()
            agonist_antigraph = future4.result()

        antagonists_obj = antagonist_graph.get_objects(['SmallMol'])
        self.drugs.update(antagonists_obj)

        agonists_obj = agonist_graph.get_objects(['SmallMol'])
        self.drugs.update(agonists_obj)

        drug_names = {d.name() for d in self.drugs}
        drug_names_df = df.from_dict({'Name':list(drug_names)})

        drug_df = self.load_pandas(drug_names_df,map2type=['SmallMol'])
        
        drug_df = self.addrank2(drug_df,antagonist_graph)
        drug_df = self.addrank2(drug_df,agonist_graph)
        drug_df = self.subtract_antigraph(drug_df,antagonist_antigraph)
        drug_df = self.subtract_antigraph(drug_df,agonist_antigraph)

        drug_df = drug_df.greater_than(0,DRUG2TARGET_REGULATOR_SCORE)
        drug_df._name_ = f'Drugs for {self.input_names()}'
        print('Found %d drugs for %d antagonist targets and %d agonist targets for %s' % 
            (len(drug_df),len(need_antagonists_ids),len(need_agonists_ids),self.input_names()),flush=True)
        return drug_df
        

    def link2disease_concepts(self,in_df:df):
        colname = self._disease2str()+' clinical trials'
        self.set_how2connect(['ClinicalTrial'],[],'',how2clone=REFERENCE_IDENTIFIERS)
        linked_row_count,linked_ids,drug_df = self.link2concept(colname,self.disease_ids,in_df)
        print('%d drugs on clinical trials for %s' % (linked_row_count,self._disease2str()),flush=True)
        self.column_ranks[4] = list(drug_df.columns)[-1]

        print('Linking %d drugs with %s' % (len(self.drugs),self._disease2str()),flush=True)
        colname = 'regulation of '+ self._disease2str()
        self.set_how2connect(['Regulation'],[],'',how2clone=REFERENCE_IDENTIFIERS)
        linked_row_count,linked_ids,drug_df = self.link2concept(colname,self.disease_ids,drug_df)
        print('%d drugs regulating %s' % (linked_row_count,self._disease2str()),flush=True)
        self.column_ranks[5] = list(drug_df.columns)[-1]

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


    def score_drugs(self):
        # making list of drugs
        start_time = time.time()
        drug_df = self.init_drug_df()
        self.column_ranks = {0:DRUG2TARGET_REGULATOR_SCORE}

        if self.param['add_bibliography']:
            with ThreadPoolExecutor(max_workers=1,thread_name_prefix='ScoreDrugs') as sd:
                drugs_etm_thread = sd.submit(self.etm_refs2df,drug_df,self.input_names())
                concepts_thread = sd.submit(self.link2disease_concepts,drug_df)
                new_drug_df = concepts_thread.result()
                drugs_df_with_etmrefs = drugs_etm_thread.result()
                #merging results of two threads:
                ref_column_name = self.etm_counter.etm_ref_column_name[-1]
                new_drug_df[ref_column_name] = drugs_df_with_etmrefs[ref_column_name]
                new_drug_df['DOIs'] = drugs_df_with_etmrefs['DOIs']
        else:
            new_drug_df = self.link2disease_concepts(drug_df)

        drug_df_ = self.make_count_df(new_drug_df,'rawDrugs')
        self.add2raw(drug_df_)

        ranks = sorted(self.column_ranks)
        columns2norm = [self.column_ranks[k] for k in ranks]
        self.normalize('rawDrugs','Drugs','Name',columns2norm)

        print("Drug ranking was done in %s" % self.execution_time(start_time), flush=True)


    def make_report(self):
        start_time = time.time()
        super().make_report()
        if self.param['add_bibliography']:
            etm_thread = Thread(target=self.add_etm_bibliography(),name='ETM4diseaseTargets')
            etm_thread.start()

            score_drugs_thread = Thread(target=self.score_drugs(),name='ScoreDrugs')
            score_drugs_thread.start()

            etm_thread.join()
            score_drugs_thread.join()
        else:
            self.score_drugs()
        
        targets4antagonist_df = df(self.report_pandas[ANTAGONIST_TARGETS_WS])
        targets4antagonist_df = targets4antagonist_df.merge_psobject(self.direct_target2drugs,'Directly Inhibited by','Name',values21cell=True)
        self.report_pandas[ANTAGONIST_TARGETS_WS] = targets4antagonist_df.merge_psobject(self.indirect_target2drugs,'Indirectly Inhibited by','Name',values21cell=True)
        
        targets4atagonist_df = df(self.report_pandas[AGONIST_TARGETS_WS])
        targets4atagonist_df =  targets4atagonist_df.merge_psobject(self.direct_target2drugs,'Directly Activated by','Name',values21cell=True)
        self.report_pandas[AGONIST_TARGETS_WS] = targets4atagonist_df.merge_psobject(self.indirect_target2drugs,'Indirectly Activated by','Name',values21cell=True)
        print('Drug repurposing for %s was done in %s' % (self._disease2str(), self.execution_time(start_time)))
        

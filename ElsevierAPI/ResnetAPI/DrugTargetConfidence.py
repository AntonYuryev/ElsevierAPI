from .ResnetAPISession import APISession,NO_REL_PROPERTIES,OQL,PS_REFIID_TYPES, REFERENCE_IDENTIFIERS
from .ResnetGraph import ResnetGraph,EFFECT,REFCOUNT,PROTEIN_TYPES,PSObject
import shutil
from concurrent.futures import ThreadPoolExecutor
import math

CONSISTENCY = 'Consistency coefficient'

class DrugTargetConsistency(APISession):
    '''
    Class to calcuate
    self.drug_target_confidence = {(drug.uid(),target.uid(),effect):1 + consistency_coefficient},  where\n
    consistency_coefficient = (consistency_counter-inconsistency_counter) / (consistency_counter+inconsistency_counter)   
    '''
    network = ResnetGraph()
    drugs2targets = ResnetGraph()
    cache_name = 'drug2target'
    cache_ent_props = ['Name','PharmaPendium ID']
    cache_rel_props = ['URN',EFFECT,REFCOUNT,CONSISTENCY]
    drug_target_confidence = dict()
    debug = False
    cache_was_modified = False


    def __init__(self,*args,**kwargs):
        '''
        APIconfig - args[0]\nwhat2retrieve - defaults NO_REL_PROPERTIES, other options -\n
        [DATABASE_REFCOUNT_ONLY,REFERENCE_IDENTIFIERS,BIBLIO_PROPERTIES,SNIPPET_PROPERTIES,ONLY_REL_PROPERTIES,ALL_PROPERTIES]
        connect2server - default True, set to False to run script using data in __pscache__ files instead of database
        '''
        #session_kwargs = dict(kwargs)
        super().__init__(*args,**kwargs)
        self.add_ent_props(self.cache_ent_props)
        self.add_rel_props(self.cache_rel_props)


    def load_cache(self):
        get_targets = f'SELECT Entity WHERE objectType = ({PROTEIN_TYPES})' # works without join to string!!!!
        oql_query = f'SELECT Relation WHERE Effect = (positive,negative) AND NeighborOf downstream ({OQL.select_drugs()}) AND NeighborOf upstream ({get_targets})'
        self.network = super().load_cache(self.cache_name,[oql_query],self.cache_ent_props,self.cache_rel_props)
        return self.network
    

    def load_drug_graph(self,for_targets:set,limit2drugs:set={},directly_from_db=False):
        '''
        Input
        -----
        for_targets,limit2drugs - {PSObject}

        Return
        ------
        self.drugs2targets - subgraph of self.network limited to "for_targets" and "limit2drugs"
        '''
        if not self.network:
            if directly_from_db:
                rn = f'Loading drugs for {len(for_targets)} targets from database'
                target_names = {t.name() for t in for_targets}
                if limit2drugs:
                    get_targets = 'SELECT Entity WHERE Name = ({props2})' # works without join to string!!!!
                    drug_names = {d.name() for d in limit2drugs}
                    get_drugs = 'SELECT Entity WHERE Name = ({props1})'
                    oql_query = f'SELECT Relation WHERE Effect = (positive,negative) AND NeighborOf downstream ({get_drugs}) AND NeighborOf upstream ({get_targets})'
                    self.drugs2targets = self.iterate_oql2(oql_query,drug_names,target_names,request_name=rn)
                else:
                    get_targets = 'SELECT Entity WHERE Name = ({props})'
                    get_drugs = OQL.select_drugs()
                    oql_query = f'SELECT Relation WHERE Effect = (positive,negative) AND NeighborOf downstream ({get_drugs}) AND NeighborOf upstream ({get_targets})'
                    self.drugs2targets = self.iterate_oql(oql_query,target_names,request_name=rn)
            else:
                self.load_cache()
        
        assert(self.network)
        if limit2drugs:
            self.drugs2targets = self.network.get_subgraph(limit2drugs,for_targets)
        else:
            self.drugs2targets = self.network.neighborhood(for_targets)

        return self.drugs2targets

    
    def annotate_network(self,drug_target_tup:list):
        '''
        Input
        -----
        drug_target_tup - [tup([PSObject,PSObject])]
        '''
        drugs_need_consistency = {t[0] for t in drug_target_tup}
        drug_names = {d.name() for d in drugs_need_consistency}
        target_names = {t[1].name() for t in drug_target_tup}

        drug2target2disease_graph = self.common_neighbors_with_effect(
                    with_entity_types=['Disease','Virus'], 
                    of_dbids1=list(drug_names),
                    reltypes12=['Regulation'],
                    dir12='',
                    and_dbids3=list(target_names),
                    reltypes23=['Regulation','QuantitativeChange'],
                    dir23='',
                    id_type='Name')
        
        if drug2target2disease_graph:
            drug2target2disease_graph= drug2target2disease_graph.compose(self.drugs2targets)
            # re-annotating self.drug2target_network for caching CONSISTENCY after algorithm run
            update4consistency = dict()
            for drug in drugs_need_consistency:
                drug_consistensies = self.__consistency4(drug,drug2target2disease_graph)
                update4consistency.update(drug_consistensies)
                   
            if update4consistency:
                print(f'Additional consistency coefficients for {len(update4consistency)} drug-target pairs were calculated using data from database')
                self.drug_target_confidence.update(update4consistency)

                for (drug_uid,target_uid,effect),consistency_coefficient in update4consistency.items():
                    drug2target_rels = self.network.find_relations([drug_uid],[target_uid],with_effects=[effect])
                    for rel in drug2target_rels:
                        self.network.set_edge_annotation(drug_uid,target_uid,rel.urn(),CONSISTENCY,[consistency_coefficient])
                self.cache_was_modified = True
        
    
    def load_confidence_dict(self,for_targets:set,limit2drugs:set={}): 
        '''
        Input
        -----
        Nodes in drugs2targets_dbgraph must have dbids

        Return
        -----
        self.drug_target_confidence = {(drug.uid(),target.uid()):1 + consistency_coefficient},  where\n
        consistency_coefficient = (consistency_counter-inconsistency_counter) / (consistency_counter+inconsistency_counter)
        '''
        if not self.drugs2targets:
            self.load_drug_graph(for_targets,limit2drugs)

        # first loading consistency coefficients from cache
        dt_need_consistency = list()

        drug_target_with_consistency = 0
        drugs_have_consistency = set()
        targets_have_consistency = set()
        for d,t,rel in self.drugs2targets.edges.data('relation'):
            drug = self.drugs2targets._psobj(d)
            target = self.drugs2targets._psobj(t)
            try:
                consistensy = float(rel[CONSISTENCY][-1])
                self.drug_target_confidence[(drug.uid(),target.uid(),rel.effect())] = consistensy
                drug_target_with_consistency += 1
                drugs_have_consistency.add(drug)
                targets_have_consistency.add(target)
            except KeyError:
                dt_need_consistency.append((drug,target))
    
        print(f'Found consistency coefficients for {drug_target_with_consistency} drug-target pairs ({len(drugs_have_consistency)} drugs, {len(targets_have_consistency)} targets) in "{self.cache_name}.rnef" file')
        if dt_need_consistency:
            print(f'\n\n{len(dt_need_consistency)} drug-target pairs need consistency calculation')
            dt_need_consistency.sort(key=lambda x:self.drugs2targets.in_degree(x[1].uid()),reverse=True)

           # step = 2000
            def next1000(start:int,dt_need_consistency:list):
                drug_set = set()
                for i in range(start, len(dt_need_consistency)):
                    drug_set.add(dt_need_consistency[i])
                    if len(drug_set) == 1000:
                        return i
                return len(dt_need_consistency)

            i = 0
            while i < len(dt_need_consistency):
                start = i
                i = next1000(start, dt_need_consistency)
                dt = dt_need_consistency[start:i]
                self.annotate_network(dt)
                print(f'\Calculated consistency coefficients for {i} out of {len(dt_need_consistency)} drug-target pairs')

 
    def save_network(self):
        if self.cache_was_modified:
            self.replace_cache(self.cache_name,self.network,self.cache_ent_props,self.cache_rel_props)
        return
        

    def consistency_coefficient(self,for_drug:PSObject,acting_on_target:PSObject, with_effect:str):
        try:
            return self.drug_target_confidence[(for_drug.uid(),acting_on_target.uid(),with_effect)]
        except KeyError:
            return 1
            

    @staticmethod
    def __consistency4(drug:PSObject,in_drug2diseases2target:ResnetGraph):
        '''
        Returns
        -------
        {(drug.uid(),target.uid(),effect):1 + consistency_coefficient},  where\n
        consistency_coefficient = (consistency_counter-inconsistency_counter) / 
                                                (consistency_counter+inconsistency_counter)   
        '''
        drug_neighbors = in_drug2diseases2target.get_neighbors([drug])
        if not drug_neighbors:
            return '',0
        diseases = [n for n in drug_neighbors if n.objtype() in ['Disease','Virus','CellProcess']]
        targets = [n for n in drug_neighbors if n.objtype() in PROTEIN_TYPES]

        drug2target_consistencies = dict()
        for target in targets:
            known_drug2target_effect = in_drug2diseases2target.effect_vote(drug,target)
            if not known_drug2target_effect: 
                continue
            effect_str = 'positive' if known_drug2target_effect == 1 else 'negative'
            # marking that consistency was considered for saving imfo to cache file
            my_tuple = tuple([drug.uid(),target.uid(),effect_str])
            drug2target_consistencies[my_tuple] = 1.0        
            consistency_counter = 0
            inconsistency_counter = 0
            for disease in diseases:
                known_target2disease_effect = in_drug2diseases2target.effect_vote(target,disease,any_direction=True)
                # any_direction must be True because Regulation and QuantitativeChange have opposite directions
                if not known_target2disease_effect: continue
                predicted_drug2disease_effect = known_drug2target_effect*known_target2disease_effect
            
                known_drug2disease_effect = in_drug2diseases2target.effect_vote(drug,disease)
                if not known_drug2disease_effect: continue

                if predicted_drug2disease_effect == known_drug2disease_effect:
                    consistency_counter += 1
                else:
                    inconsistency_counter += 1

            consistency_sum = consistency_counter+inconsistency_counter
            if consistency_sum:
                consistency = float(consistency_counter-inconsistency_counter)/(consistency_sum)
            else:
                consistency = 0.0

            drug2target_consistencies[my_tuple] = 1 + consistency
        return drug2target_consistencies 


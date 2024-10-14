from .ResnetAPIcache import APIcache
from .PathwayStudioGOQL import OQL
from .ResnetGraph import ResnetGraph,EFFECT,REFCOUNT,PROTEIN_TYPES,PSObject,CONSISTENCY,PSRelation
from math import log, sqrt
from collections import defaultdict


class DrugTargetConsistency(APIcache):
    '''
    Class to calcuate
    self.drug_target_confidence = {(drug.uid(),target.uid(),effect):1 + consistency_coefficient},  where\n
    consistency_coefficient = (consistency_counter-inconsistency_counter) / (consistency_counter+inconsistency_counter)   
    '''
    drugs2targets = ResnetGraph()
    cache_name = 'drug2target'
    cache_ent_props = ['Name','PharmaPendium ID']
    cache_rel_props = ['URN',EFFECT,REFCOUNT,CONSISTENCY,'pX','Affinity']
    drug_target_confidence = dict()
    debug = False
    predict_effect4 = {'_4enttypes':['SmallMol'],'_4reltypes':['Binding']}


    def __init__(self,*args,**kwargs):
        '''
        APIconfig - args[0]\nwhat2retrieve - defaults NO_REL_PROPERTIES, other options -\n
        [DATABASE_REFCOUNT_ONLY,REFERENCE_IDENTIFIERS,BIBLIO_PROPERTIES,SNIPPET_PROPERTIES,ONLY_REL_PROPERTIES,ALL_PROPERTIES]
        connec
        t2server - default True, set to False to run script using data in __pscache__ files instead of database
        '''
        get_targets = f'SELECT Entity WHERE objectType = ({PROTEIN_TYPES})' # works without join to string!!!!
        get_drugs = OQL.select_drugs()
        oql_query = f'SELECT Relation WHERE Effect = (positive,negative) AND NeighborOf downstream ({get_drugs}) AND NeighborOf upstream ({get_targets})'
        qt = tuple([oql_query,'download drugs2targets regulatory network'])
        oql_queries = [qt]

        oql_query1 = f'SELECT Relation WHERE objectType = Binding AND NeighborOf ({get_drugs}) AND NeighborOf ({get_targets})'
        qt = tuple([oql_query1,'download drugs2targets binding network'])
        oql_queries.append(qt)

        my_kwargs = {
            'cache_name': self.cache_name, # used for folder in self.data_dir and in RNEF dump name
            'oql_queries':oql_queries, # [(oql_query,request_name),...] list of GOQL queries tuples for fetching data from database
            'ent_props' : self.cache_ent_props, # entity props to fetch from database
            'rel_props' : self.cache_rel_props, # relation props to fetch from database
            'ranks4simplifying' : ['DirectRegulation','Binding','Regulation'], # parameter for making simple graph
            'refprop2rel': {'pX':'Affinity'}, # {reference_propname:relation_propname} - moves refprops to relprops using PSRelation._refprop2rel()
            'refprop_minmax':1, # parameter for PSRelation._refprop2rel(refprop2rel,refprop_minmax)
            #if min_max < 0 assigns single value to relprop that is the minimum of all ref_prop values
            #if min_max > 0 assigns single value to relprop that is the maximum of all ref_prop values
            #min_max works only if ref_prop values are numerical
            'predict_effect4' : self.predict_effect4, # [_4enttypes:list,_4reltypes:list] - parameters for ResnetGraph.predict_effect4()
        }

        my_kwargs.update(kwargs)
        my_kwargs.pop('data_dir','') # this cache uses default __pscache__ directory
        super().__init__(*args,**my_kwargs)


    def d2t_from_db(self,for_targets:set[PSObject],limit2drugs:set[PSObject]={}):
        rn = f'Loading drugs for {len(for_targets)} targets from database'
        target_names = {t.name() for t in for_targets}
        if limit2drugs:
            get_targets = 'SELECT Entity WHERE Name = ({props2})' # works without join to string!!!!
            drug_names = {d.name() for d in limit2drugs}
            get_drugs = 'SELECT Entity WHERE Name = ({props1})'
            oql_query = f'SELECT Relation WHERE Effect = (positive,negative) AND NeighborOf downstream ({get_drugs}) AND NeighborOf upstream ({get_targets})'
            return self.iterate_oql2(oql_query,drug_names,target_names,request_name=rn)
        else:
            get_targets = 'SELECT Entity WHERE Name = ({props})'
            get_drugs = OQL.select_drugs()
            oql_query = f'SELECT Relation WHERE Effect = (positive,negative) AND NeighborOf downstream ({get_drugs}) AND NeighborOf upstream ({get_targets})'
            return self.iterate_oql(oql_query,target_names,request_name=rn)


    def load_drug_graph(self,for_targets:set[PSObject],limit2drugs:set[PSObject]={},directly_from_db=False):
        '''
        output:
            self.drugs2targets - subgraph of self.network limited to "for_targets" and "limit2drugs"
        '''
        if not self.network:
            if directly_from_db:
                self.drugs2targets = self.d2t_from_db(for_targets,limit2drugs)
            else:
                self.network = self._load_cache(cache_name=self.cache_name)
        
        assert(self.network)
        if limit2drugs:
            self.drugs2targets = self.network.get_subgraph(limit2drugs,list(for_targets))
        else:
            self.drugs2targets = self.network.neighborhood(for_targets)
            
        all_targets = {PSObject(y) for x,y in self.drugs2targets.nodes(data=True) if self.drugs2targets.in_degree(x)}
        missing_targets = for_targets.difference(all_targets)
        add2network = self.d2t_from_db(missing_targets,limit2drugs)
        add2network = add2network.subtract(self.drugs2targets)
        if add2network:
            self.add2cache(add2network)
            self.drugs2targets = self.drugs2targets.compose(add2network)

        return self.drugs2targets

    
    def annotate_network(self,drug_target_tup:list,replace_consistencies=False):
        '''
        Input
        -----
        drug_target_tup - [tup([PSObject,PSObject])]

        Annotates
        ---------
        in self.network with CONSISTENCY all drug-targets relations for all drug-target pairs from "drug_target_tup"\n
        there more drug-target combinations than drug-target tuples in drug_target_tup
        '''
        drugs_need_consistency = {tup[0] for tup in drug_target_tup}
        drug_names = {d.name() for d in drugs_need_consistency}
        targets4drugs_need_consistency = {tup[1] for tup in drug_target_tup}
        target_names = {t.name() for t in targets4drugs_need_consistency}

        d2t2dG = self.common_neighbors_with_effect(
                    with_entity_types=['Disease','Virus'], 
                    of_dbids1=list(drug_names),
                    reltypes12=['Regulation'],
                    dir12='',
                    and_dbids3=list(target_names),
                    reltypes23=['Regulation','QuantitativeChange'],
                    dir23='',
                    id_type='Name')
        
        if d2t2dG:# re-annotating self.drug2target_network for caching CONSISTENCY after algorithm run
            d2t_need = self.drugs2targets.get_subgraph(list(drugs_need_consistency),list(targets4drugs_need_consistency))
            d2t2dG = d2t2dG.compose(d2t_need)
            update4consistency = defaultdict(float) # {(drug_uid,target_uid,effect):consistency_coefficient}
            for drug in drugs_need_consistency:
                for (drug_uid, target_uid, effect), consistency_coefficient in self.__consistency4(drug, d2t2dG).items():
                    update4consistency[(drug_uid, target_uid, effect)] = consistency_coefficient
                   
            if update4consistency:
                # update4consistency can be empty if drug-target pair have no common diseases linked with known Effect 
                print(f'Additional consistency coefficients for {len(update4consistency)} drug-target pairs were calculated using data from database')
                self.drug_target_confidence.update(update4consistency)

                for (drug_uid,target_uid,effect),consistency_coefficient in update4consistency.items():
                    drug2target_rels = self.network.find_relations([drug_uid],[target_uid],with_effects=[effect])
                    for rel in drug2target_rels:
                        if replace_consistencies or CONSISTENCY not in rel:
                            self.cache_was_modified = True
                            self.network.set_edge_annotation(drug_uid,target_uid,rel.urn(),CONSISTENCY,[consistency_coefficient])

        return
                
    
    def load_confidence_dict(self,for_targets:set,limit2drugs=set(),load_fromdb=False): 
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
            self.load_drug_graph(for_targets,limit2drugs,load_fromdb)

        # first loading consistency coefficients from cache
        dt_need_consistency = list() #[(PSObject,PSObject)]
        if load_fromdb:
            for d,t,rel in self.drugs2targets.edges.data():
                drug = self.drugs2targets._psobj(d)
                target = self.drugs2targets._psobj(t)
                dt_need_consistency.append((drug,target))
        else:
            dt_with_consistency = 0
            drugs_have_consistency = set()
            targets_have_consistency = set()
            for drug,target,rel in self.drugs2targets.iterate():
                if rel.objtype() != 'Binding': # Binding cannot have consistency
                    if CONSISTENCY in rel:
                        consistensy = float(rel[CONSISTENCY][-1])
                        self.drug_target_confidence[(drug.uid(),target.uid(),rel.effect())] = consistensy
                        dt_with_consistency += 1
                        drugs_have_consistency.add(drug)
                        targets_have_consistency.add(target)
                    else:
                        dt_need_consistency.append((drug,target))
    
            print(f'Found consistency coefficients for {dt_with_consistency} drug-target pairs ({len(drugs_have_consistency)} drugs, {len(targets_have_consistency)} targets) in "{self.cache_name}.rnef" file')
        
        if dt_need_consistency:
            print(f'\n\n{len(dt_need_consistency)} drug-target pairs need consistency calculation')
            dt_need_consistency.sort(key=lambda x:self.drugs2targets.in_degree(x[1].uid()),reverse=True)
        
           # step = 2000
            def next1000(start:int,dt_need_consistency:list):
                drug_counter = set()
                for i in range(start, len(dt_need_consistency)):
                    drug_counter.add(dt_need_consistency[i])
                    if len(drug_counter) == 1000: # when number of unique drugs is 1000 stop and caclulate
                        return i
                return len(dt_need_consistency)

            i = 0
            while i < len(dt_need_consistency):
                start = i
                i = next1000(start, dt_need_consistency)
                dt = dt_need_consistency[start:i]
                self.annotate_network(dt,load_fromdb)
                print(f'Calculated consistency coefficients for {i} out of {len(dt_need_consistency)} drug-target pairs\n')
            return
        else:
            print(f'No drug-target pairs need consistency calculation.  All coefficients were found in {self.cache_name}')
            return


    def consistency_coefficient(self,for_drug:PSObject,acting_on_target:PSObject, with_effect:str):
        try:
            return self.drug_target_confidence[(for_drug.uid(),acting_on_target.uid(),with_effect)]
        except KeyError:
            return 1
    
    '''
    @staticmethod
    def __consist_coeffOLD(consistent_count:int,inconsistent_count:int):
       total_count = consistent_count+inconsistent_count
       return 1.0 + float(consistent_count-inconsistent_count)/(total_count) if total_count else 1.0
    '''

    @staticmethod
    def __consist_coeff(consistent_count:int,inconsistent_count:int):
        '''
        consistency becomes negative when inconsistent_count == consistent_count\n
        if drug-target has only one disease in common consistency = 0.315
        '''
        zero_adj = 0.1
        zscore= log((consistent_count+zero_adj)/(inconsistent_count+zero_adj),10)
        std = sqrt((1.0/(consistent_count+zero_adj)) + (1.0/(inconsistent_count+zero_adj)))
        return zscore/std 


    @staticmethod
    def __consistency4(drug:PSObject,in_d2t2d:ResnetGraph)->dict[tuple[int,int,str],float]:
        '''
        Returns
        -------
        {(drug.uid(),target.uid(),effect):1 + consistency_coefficient},  where\n
        consistency_coefficient = self.__consist_coeff()
        '''
        targets = list()
        diseases = list()
        target_types = PROTEIN_TYPES + ['SmallMol']
        for r,target,rel in in_d2t2d.targets_of(drug):
            target_objtype = target.objtype()
            if target_objtype in ['Disease','Virus','CellProcess']:
                diseases.append(target)
            elif target_objtype in target_types:
                if CONSISTENCY not in rel:
                    targets.append(target)

        if not targets or not diseases:
            return dict()

        d2t_consistencies = defaultdict(lambda: 1.0)
        for target in targets:
            known_drug2target_effect = in_d2t2d.effect_vote(drug,target)
            if known_drug2target_effect: # not unknown
                effect_str = 'positive' if known_drug2target_effect == 1 else 'negative'
                consistency_counter = 0
                inconsistency_counter = 0
                for disease in diseases:
                    known_target2disease_effect = in_d2t2d.effect_vote(target,disease,any_direction=True)
                    # any_direction must be True because Regulation and QuantitativeChange have opposite directions
                    if not known_target2disease_effect: continue
                    predicted_drug2disease_effect = known_drug2target_effect*known_target2disease_effect
                
                    known_drug2disease_effect = in_d2t2d.effect_vote(drug,disease)
                    if not known_drug2disease_effect: continue

                    if predicted_drug2disease_effect == known_drug2disease_effect:
                        consistency_counter += 1
                    else:
                        inconsistency_counter += 1

                consistency_coeff = DrugTargetConsistency.__consist_coeff(consistency_counter,inconsistency_counter)
                d2t_consistencies[(drug.uid(),target.uid(),effect_str)] = round(consistency_coeff,3)
        return d2t_consistencies

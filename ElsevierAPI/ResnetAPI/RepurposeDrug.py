from .TargetIndications import Indications4targets,OQL,UNKEFFECTDF,ANTAGONIST,AGONIST,time,SNIPPET_PROPERTIES
from ..pandas.panda_tricks import df
from .ResnetGraph import ResnetGraph,PSObject,PSRelation,OBJECT_TYPE,PROTEIN_TYPES
from ..ETM_API.references import MEASUREMENT,PATENT_APP_NUM,PS_REFERENCE_PROPS,PS_SENTENCE_PROPS
from ..ReaxysAPI.Reaxys_API import Reaxys_API
from ..utils import  execution_time
from concurrent.futures import ThreadPoolExecutor
from .SemanticSearch import RANK


INHIBIT = -1 # indication must be inhibited by a drug
ACTIVATE = 1 # indication must be activated by a drug
ALL_EFFECTS = 0
INDICATION_COUNT = 'IndicationsCnt'
TOXICITY_COUNT = 'ToxicityCnt'
INDICATIONS = 'Indications'
TOXICITIES = 'Toxicities'

class RepurposeDrug(Indications4targets):
    pass
    def __init__(self, *args, **kwargs):
        '''
        APIconfig - args[0]
        '''
        my_kwargs = {
                'input_compound' : '', 
                'similars' : [],
                'drug_effect': INHIBIT, # INHIBIT for indications, ACTIVATE for toxicitites
                'mode_of_action': ANTAGONIST, # AGONIST
                'indication_types': ['Disease'], #['CellProcess','Virus','Pathogen']
                }
        my_kwargs.update(kwargs)

        super().__init__(*args,**my_kwargs)
        self.drugs = list() # list of drug for repurposing.  Usually len(self.drugs) == 1,
        # self.drugs is list because input drug name may find more than one drug in the database 
        # drug(s) are annoated with 'activates_targets' and 'inhibit_targets' attribute after set_drugs() is executed
        self.drug2targetG = ResnetGraph()
        self.SimilarDrugs = set()
        self.DrugIndications = set()
        self.DrugToxicities = set()

        self.activated_targets = set()
        self.inhibited_targets = set()
        self.activated_partners = set()
        self.inhibited_partners = set()
        self.CellSecretingActivatedTargets = set()
        self.CellSecretingInhibitedTargets = set()

        self.target_indications = set()
        self.target_toxicities = set()
        
        self.DrugDirectAgonists = set()
        self.DrugDirectAntagonists = set()
        self.DrugIndirectAgonists = set()
        self.DrugIndirectAntagonists = set()


    def target_names_str(self):
        """
            Overwrites super().target_names_str() to ensure identical column names in report worksheets
        """
        return 'targets'


    def set_drugs(self):
        '''
        input:
            self.params
        output:
            self.drugs
            self.drug2targetG
            self.activated_targets
            self.inhibited_targets
            self.activated_partners
            self.inhibited_partners
        '''
        self.SELECTdrug = OQL.get_childs([self.params['input_compound']],['Name,Alias'],include_parents=True)
        drug_graph = self.load_graph_from_oql(self.SELECTdrug,[],['Connectivity'],get_links=False,add2self=False)
        self.drugs = drug_graph._get_nodes()
        self.drugs.sort(key=lambda x: int(x['Connectivity'][0]), reverse=True)

        if self.drugs:
          ref_ids = ['PMID','DOI',PATENT_APP_NUM]
          refprops2print = list()
          oql = f"SELECT Relation WHERE NeighborOf downstream (SELECT Entity WHERE id = ({ResnetGraph.dbids(self.drugs)}))"
          oql += ' AND Effect = (positive,negative)'
          if self.params['target_names']:
            self.add_ent_props(['Class'])
            oql += f' AND NeighborOf upstream (SELECT Entity WHERE Name = ({self.params['target_names']})'
            oql += ')'
            my_session = self._clone_session(to_retrieve = SNIPPET_PROPERTIES)
            my_session.add_rel_props(self.relProps)
            self.drug2targetG = my_session.process_oql(oql,'Find drug-targets relations in Pathway Studio')
          else:
            print('Input parameters contain no drug target names. Attempting to find targets in Reaxys')
            self.drug2targetG = self.Reaxys_targets()
            print(f'Found {len(self.drug2targetG)-len(self.drugs)} in Reaxys')
            if self.drug2targetG:
              refprops2print = ['Organism',MEASUREMENT,'pX'] # additional props from Reaxys
            else:
              print('Attempting to find targets in Pathway Studio')
              oql += f' AND objectType = DirectRegulation'
              oql += f' AND NeighborOf upstream (SELECT Entity WHERE objectType = ({PROTEIN_TYPES})'
              my_session = self._clone_session(to_retrieve = SNIPPET_PROPERTIES)
              my_session.add_rel_props(self.relProps)
              self.drug2targetG = my_session.process_oql(oql,'Find targets in Pathway Studio')
              if self.drug2targetG:
                dt_rels = list(self.drug2targetG._psrels())
                if len(dt_rels) > 10: # taking top 10 targets in PS
                    dt_rels.sort(key=lambda x: len(x.refs()),reverse=True)
                    dt_rels = dt_rels[:10]
                    self.drug2targetG = self.drug2targetG.subgraph_by_rels(dt_rels)
              
                [self.drug2targetG[r][t][urn]['relation'].set_affinity() for r,t,urn in self.drug2targetG.edges(keys=True)]
                refprops2print = PS_SENTENCE_PROPS
              else:
                print(f'Targets for {self.drug_names_str()} cannot be found')
        else:
          print('Add drug-target relations into Pathway Studio database')

        urn2regweight = dict()
        urn2tarweight = dict()
        for r,t,rel in self.drug2targetG.iterate():
            if 'Affinity' in rel:
              affinity = float(rel['Affinity'][0])
              node_weight =  1 + affinity/12
              target_urn = t.urn()
              urn2regweight[target_urn] = [node_weight]
              urn2tarweight[target_urn] = [node_weight]
            
        if urn2regweight:
            self.drug2targetG.set_node_annotation(urn2regweight,'regulator weight')
            self.targets_have_weights = True
        if urn2tarweight:
            self.drug2targetG.set_node_annotation(urn2tarweight,'target weight')
            self.targets_have_weights = True

        self.activated_targets = self.__targets('positive')
        self.inhibited_targets = self.__targets('negative')
        
        self.drug2targetG = self.drug2targetG.compose(self.set_partners())
        ref_df = self.drug2targetG.snippets2df('Targets',add_nodetype=True,ref_identifiers=ref_ids,ref_sentence_props=refprops2print)
        ref_df.tab_format['tab_color'] = 'brown'
        self.add2report(ref_df)
        
        all_targets = self.activated_targets|self.inhibited_targets
        self.GVs2DiseaseGraph = self._GVindicationsG(all_targets)
        self.targetGVs = self.GVs2DiseaseGraph._psobjs_with('GeneticVariant','ObjTypeName')
        return


    def __targets(self,effect:str)->set[PSObject]:
        '''
        effect in ['positive','negative']
        '''
        assert(effect in ['positive','negative'])
        dt_rels = list(self.drug2targetG._psrels())
        dt_rels = [x for x in dt_rels if x.effect() == effect]
        my_targets = set()
        for rel in dt_rels:
            for reg_uid,target_uid in rel.get_regulators_targets():
                my_targets.add(self.drug2targetG._get_node(target_uid))
        childs,parents = self.load_children4(my_targets)
        my_targets.update(childs)
        return my_targets
    

    def find_drug_indications(self)->set[PSObject]:
        indications = set()
        indic_str = OQL.join_with_quotes(self.params['indication_types'])

        REQUEST_NAME = 'Find {drug} indications by clinical trials'.format(drug=self.params['input_compound'])
        OQLquery = 'SELECT Relation WHERE objectType = ClinicalTrial AND NeighborOf ({select_drug}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_types}))'
        ClinicalTrialIndictions = self.process_oql(OQLquery.format(select_drug=self.SELECTdrug, indication_types = indic_str),REQUEST_NAME)
        if isinstance(ClinicalTrialIndictions,ResnetGraph):
            indications = set(ClinicalTrialIndictions.psobjs_with(only_with_values=self.params['indication_types']))
            print('Found %d indications in %s clinical trials' % (len(indications), self.params['input_compound']))

        REQUEST_NAME = 'Find {drug} other indications'.format(drug=self.params['input_compound'])
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = negative AND NeighborOf({select_drug}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_types}))'
        oql_query = OQLquery.format(select_drug=self.SELECTdrug,indication_types = indic_str)
        LiteratureIndications = self.process_oql(oql_query, REQUEST_NAME)
        if isinstance(LiteratureIndications,ResnetGraph):
            other_indications = LiteratureIndications.psobjs_with(only_with_values=self.params['indication_types'])
            indications.update(other_indications)
            print('Found %d indications reported in scientific literature for %s' % (len(other_indications), self.params['input_compound']))
    
        return indications
    

    def find_drug_toxicities(self)->set[PSObject]:
        toxicities = set()
        indic_str = OQL.join_with_quotes(self.params['indication_types'])

        REQUEST_NAME = 'Find {drug} toxicties'.format(drug=self.params['input_compound'])
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = positive AND NeighborOf({select_drug}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_types}))'
        oql_query = OQLquery.format(select_drug=self.SELECTdrug,indication_types = indic_str)
        LiteratureIndications = self.process_oql(oql_query, REQUEST_NAME)
        if isinstance(LiteratureIndications,ResnetGraph):
            toxicities = LiteratureIndications.psobjs_with(only_with_values=self.params['indication_types'])
            print('Found %d toxicities reported in scientific literature for %s' % (len(toxicities), self.params['input_compound']))
        return set(toxicities)

    
    def drug_names(self)->list:
        return [x.name() for x in self.drugs]
    

    def drug_names_str(self)->str:
        return ','.join(self.drug_names())


    def find_indications4similars(self)->list[PSObject]:
        if not self.params['similars']: return list()
        indications2return = set()
        indic_str = OQL.join_with_quotes(self.params['indication_types'])
        similar_drugs_str = OQL.join_with_quotes(self.params['similars'])
        select_similar_drugs = f'SELECT Entity WHERE (Name,Alias) = ({similar_drugs_str})'

        REQUEST_NAME = 'Find indications by clinical trials for {similars}'.format(similars=','.join(self.params['similars']))
        OQLquery = 'SELECT Relation WHERE objectType = ClinicalTrial AND NeighborOf ({select_drugs}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_types}))'
        SimilarDrugsClinicalTrials = self.process_oql(OQLquery.format(select_drugs=select_similar_drugs, indication_types=indic_str), REQUEST_NAME)
        if isinstance(SimilarDrugsClinicalTrials,ResnetGraph):
            found_indications = SimilarDrugsClinicalTrials.psobjs_with(only_with_values=self.params['indication_types'])
            indications2return.update(found_indications)
            print('Found %d indications in clinical trials for drugs similar to %s' % (len(found_indications), self.params['input_compound']))

        effect = 'negative' if self.params['drug_effect'] == INHIBIT else 'positive'
        REQUEST_NAME = 'Find all other indications for {similars}'.format(similars=','.join(self.params['similars']))
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = {eff} AND NeighborOf ({select_drugs}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_types}))'
        OQLquery = OQLquery.format(eff=effect,select_drugs=select_similar_drugs,indication_types=indic_str)
        SimilarDrugsIndications = self.process_oql(OQLquery, REQUEST_NAME)
        if isinstance(SimilarDrugsIndications,ResnetGraph):
            self.SimilarDrugs = SimilarDrugsIndications.psobjs_with(only_with_values=['SmallMol'])
            found_indications = SimilarDrugsIndications.psobjs_with(only_with_values=self.params['indication_types'])
            indications2return.update(found_indications)
            print('Found %d indications reported in scientific literature or drugs similar to %s' % (len(found_indications), self.params['input_compound']))

        return list(indications2return)
    

    def find_toxicities4similars(self)->list[PSObject]:
        if not self.params['similars']: return list()
        toxicities = set()
        indic_str = OQL.join_with_quotes(self.params['indication_types'])
        similar_drugs_str = OQL.join_with_quotes(self.params['similars'])
        select_similar_drugs = f'SELECT Entity WHERE (Name,Alias) = ({similar_drugs_str})'

        REQUEST_NAME = 'Find all other indications for {similars}'.format(similars=','.join(self.params['similars']))
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = positive AND NeighborOf ({select_drugs}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_types}))'
        OQLquery = OQLquery.format(select_drugs=select_similar_drugs,indication_types=indic_str)
        SimilarDrugsIndications = self.process_oql(OQLquery, REQUEST_NAME)
        if isinstance(SimilarDrugsIndications,ResnetGraph):
            self.SimilarDrugs = SimilarDrugsIndications.psobjs_with(only_with_values=['SmallMol'])
            toxicities = SimilarDrugsIndications.psobjs_with(only_with_values=self.params['indication_types'])
            print('Found %d indications reported in scientific literature or drugs similar to %s' % (len(toxicities), self.params['input_compound']))

        return list(toxicities)


    def clean_indications(self):
      if self.params['drug_effect'] == ACTIVATE: # == report is for toxicities
# clinical trials report their indications as toxicities 
# to check if patient condition is not worsening after drug treatment
        if 'Disease' in self.params['indication_types']:
            #indications = self.DrugIndications
            before_clean_indications_count = len(self.DrugIndications)
            indication_childs,_ = self.load_children4(list(self.DrugIndications))
            all_indications = self.DrugIndications | indication_childs

            all_drugs = self.drugs+self.SimilarDrugs
            oql_query = 'SELECT Entity WHERE objectType = Disease AND id = ({ids1}) AND Connected by (SELECT Relation WHERE objectType = ClinicalTrial)  \
                to (SELECT Entity WHERE id = ({ids2}))'
            all_indications_dbid = ResnetGraph.dbids(list(all_indications))
            all_drugs_dbids = ResnetGraph.dbids(all_drugs)
            clinical_trial_graph = self.iterate_oql2(oql_query,set(all_indications_dbid),set(all_drugs_dbids))
            clinical_trial_indications = clinical_trial_graph._psobjs_with('Disease','ObjTypeName')
            clintrial_indication_parents,_ = self.Graph.find_parents4(clinical_trial_indications)
            clinical_trial_indications += list(clintrial_indication_parents)

            #removing clinical_trial_indications from the list of indications for target drugs
            self.DrugIndications = self.DrugIndications.difference(clinical_trial_indications)
            after_clean_indications_count = len(self.DrugIndications)
            print('%d toxicities were removed because they are indications for %s in clinical trials' % 
                    (before_clean_indications_count-after_clean_indications_count,self.drug_names()))
    

    def _drug_connect_params(self,direct_modulators:bool,score_antagonists:bool):
        '''
        Return
        ------
        how2connect parameters to score substances having desired drug affect on the targets 
        '''
        if score_antagonists:
        # most common case when targets must be inhibited
            if direct_modulators:
                effect_on_indication = 'positive'
                drug_class = 'direct drug antagonists'         
                concepts = self.DrugDirectAntagonists
            else:
                effect_on_indication = 'positive'
                drug_class = 'indirect drug antagonists'             
                concepts = self.DrugIndirectAntagonists
        else:
        # case if drug are agonists
            if direct_modulators:
                effect_on_indication = 'negative'
                drug_class = 'direct drug agonists'
                concepts = self.DrugDirectAgonists
            else:
                effect_on_indication = 'negative'
                drug_class = 'indirect drug agonists'
                concepts = self.DrugIndirectAgonists

        return effect_on_indication, drug_class, concepts


    def _drug_tox_params(self,direct_modulators:bool,score_antagonists:bool):
        '''
        Return
        ------
        how2connect parameters to score substances synergizing with target action on indication\n
        i.e. molecules having effect opposite to desired drug affect on targets
        '''
        if score_antagonists:
            if direct_modulators:
                effect_on_indication = 'positive'
                drug_class = 'direct drug agonists'
                concepts = self.DrugDirectAgonists
            else:
                effect_on_indication = 'positive'
                drug_class = 'indirect drug agonists'
                concepts = self.DrugIndirectAgonists
        else:
        # case if drug are agonists
            if direct_modulators:
                effect_on_indication = 'negative'
                drug_class = 'direct drug antagonists'             
                concepts = self.DrugDirectAntagonists
            else:
                effect_on_indication = 'negative'
                drug_class = 'indirect drug antagonists'             
                concepts = self.DrugIndirectAntagonists

        return effect_on_indication, drug_class, concepts
    

    def semscore4drugs(self,in_df:df)->df:
        '''
        input:
            in_df._name_ is used to determine how to score
        output:
            df with scores
        '''
        my_df = in_df.dfcopy()
        if my_df._name_ == INDICATION_COUNT:
            concept = self.params['input_compound'] + ' clinical trials'
            kwargs = {'connect_by_rels':['ClinicalTrial']}
            how2connect = self.set_how2connect (**kwargs)
            linked_entities_count, _, my_df = self.link2concept(concept, self.drugs,my_df,how2connect)
            if linked_entities_count:
              self._set_rank(my_df,concept)
            print('%s has clinical trials for %d indications' % (self.params['input_compound'], linked_entities_count))

            concept = 'Inhibited by ' + self.params['input_compound']
            kwargs = {'connect_by_rels':['Regulation'],
                  'with_effects' : ['negative'],
                  'boost_by_reltypes' :['Regulation']
                  }
            how2connect = self.set_how2connect (**kwargs)
            regulated = 'inhibited'
        else:
            concept = 'Activated by ' + self.params['input_compound']
            kwargs = {'connect_by_rels':['Regulation'],
                  'with_effects' : ['positive'],
                  'boost_by_reltypes' :['Regulation']
                  }
            how2connect = self.set_how2connect (**kwargs)
            regulated = 'activated'

        linked_entities_count, _, my_df = self.link2concept(concept, self.drugs,my_df,how2connect)
        if linked_entities_count:
          self._set_rank(my_df,concept)
        print('%d indications %s by %s' % (linked_entities_count, regulated, self.params['input_compound']))
        self.drugcolname = self._col_name_prefix+concept
        
        if hasattr(self,'similar_drugs'):
            if my_df._name_ == INDICATION_COUNT:
                concept = 'similar drugs clinical trials'
                kwargs = {'connect_by_rels':['ClinicalTrial']}
                how2connect = self.set_how2connect (**kwargs)
                linked_entities_count, _, my_df = self.link2concept(concept, self.SimilarDrugs,my_df,how2connect)
                if linked_entities_count:
                  self._set_rank(my_df,concept)
                print('%s has clinical trials for %s indications' % (','.join(self.params['similars']), linked_entities_count))

                concept = 'Inhibited by similar drugs'
                kwargs = {'connect_by_rels':['Regulation'],
                  'with_effects' : ['negative'],
                  'boost_by_reltypes' :['Regulation']
                  }
                how2connect = self.set_how2connect (**kwargs)
            else:
                concept = 'Activated by similar drugs'
                kwargs = {'connect_by_rels':['Regulation'],
                  'with_effects' : ['positive'],
                  'boost_by_reltypes' :['Regulation']
                  }
                how2connect = self.set_how2connect (**kwargs)
            linked_entities_count, _, my_df = self.link2concept(concept,self.SimilarDrugs,my_df,how2connect)
            if linked_entities_count:
                self._set_rank(my_df,concept)
            print('%d indications %s by %s' % (linked_entities_count, regulated, ','.join(self.params['similars'])))

        return my_df
    

    def report_name(self):
        indics = OQL.join_with_quotes(self.params['indication_types'])
        rep_pred =  'suggested ' if self.params['strict_mode'] else 'suggested,predicted '
        if indics == 'Disease':
            if self.params['drug_effect'] == INHIBIT:
                return str(rep_pred + 'indications for '+self.params['input_compound'])
            elif self.params['drug_effect'] == ACTIVATE:
                return str(rep_pred + 'toxicities for '+self.params['input_compound'])
            else:
                print('No drug_effect was specified in paramters')
                return ''
        else:
            regulate = ' activated by ' if self.params['drug_effect'] == ACTIVATE else ' inhibited by '
            return str(rep_pred+indics+regulate+self.params['input_compound'])


    def report_path(self,suffix='',extension='.xlsx'):
        ext = '.'+extension if extension[0] != '.' else extension    
        return str(self.data_dir+self.report_name()+suffix+ext)

    '''
    def add_drug_indication_refs(self):
        print('Printing clinical trials')
        clin_trial_graph = self.Graph.subgraph_by_relprops(['ClinicalTrial'])
        ct_pd = clin_trial_graph.snippets2df()
        ct_pd._name_ ='clin. trials'
        self.add2report(ct_pd)

        print('Printing research articles')
        research_graph = self.Graph.subgraph_by_relprops(['Regulation'])
        research_graph = self.Graph.subgraph_by_relprops([effect_str],['Effect'])
        research_ref_pd = research_graph.snippets2df()

        research_ref_pd._name_ ='articles'
        self.add2report(research_ref_pd)
  

    def __ws_suffix(self,with_targets_names=False):
        indications = ','.join(self.params['indication_types'])
        targets = ','.join(self.params['target_names'])
        ef = 'Act.' if self.params['drug_effect'] == ACTIVATE else 'Inh.'
        return targets[:10]+'-'+ef+indications[:3] if with_targets_names else ef+indications[:3]
    '''

    def indications4drugs(self):
        '''
        output:
            self.drugs
            self.drug2targetG
            self.DrugIndications
        '''
        if not self.DrugIndications:
            self.flush_dump_files()
            start_time = time.time()

            self.set_drugs()
            self.DrugIndications = self.find_drug_indications()
            self.DrugToxicities = self.find_drug_toxicities()
            #return #uncomment for fast downstream testing
            self.DrugIndications.update(self.find_indications4similars())
            self.DrugToxicities.update(self.find_toxicities4similars())

            self.clean_indications()
            print("Drug indications were loaded in %s" % execution_time(start_time))
        else:
            print('Drug indications were already loaded')
        return self.DrugIndications


    def other_effects(self)->df:
        my_session = RepurposeDrug(**self.params)
        my_session.relProps = PS_REFERENCE_PROPS
        print(f'Started retreiving sentences for indications with unknown effects')
        select_indications,indication_str = self.oql4indications()
        
        REQUEST_NAME = f'Find {indication_str} modulated by {self.drug_names_str()} with unknown effect'
        select_indications = f'SELECT Entity WHERE objectType = ({indication_str})'
        OQLquery = f'SELECT Relation WHERE objectType = Regulation AND Effect = unknown AND \
            NeighborOf({self.SELECTdrug}) AND NeighborOf ({select_indications})'
        unknown_effectsG = my_session.process_oql(OQLquery,REQUEST_NAME)
        my_targets = self.activated_targets|self.inhibited_targets
        biomarker_iG = my_session._biomarker_indicationsG(my_targets)

        all_unknownsG = unknown_effectsG.compose(self.GVs2DiseaseGraph).compose(biomarker_iG)
        known_indications = self.DrugIndications|self.DrugToxicities
        all_unknownsG.remove_nodes_from(ResnetGraph.uids(known_indications))
        unknown_indication = all_unknownsG.psobjs_with(only_with_values=my_session.params['indication_types'])
        print(f'Found {len(unknown_indication)} indications with unknown effect')
        if unknown_indication:
            other_indications_df = all_unknownsG.snippets2df(df_name=UNKEFFECTDF)
            self.add2report(other_indications_df)
        return other_indications_df

    
    def semscore4(self,targets:list[PSObject],with_effect_on_indication:str, with_partners:list[PSObject],in_df:df):
      if targets:
        targets_indication_df = super().semscore4(targets,with_effect_on_indication, with_partners,in_df) 
        targets_indication_df.fillna(0)
        return targets_indication_df, dict(self.col2rank)
      else:
          return df(),dict()
        
        
    def make_count_df(self,from_rawdf:df,with_name:str,effect_on_indication:str)->tuple[df,dict[str,int]]:
        drug_indication_df = self.semscore4drugs(from_rawdf) # need to go first to assign rank to columns
        drug_col2rank = dict(self.col2rank)

        assert(from_rawdf._name_ in [INDICATION_COUNT,TOXICITY_COUNT])
        a_targets_indication_df,a_col2rank = self.semscore4(self.activated_targets,effect_on_indication,self.activated_partners,from_rawdf)

        self.col2rank = drug_col2rank
        i_targets_indication_df,i_col2rank = self.semscore4(self.inhibited_targets,effect_on_indication,self.inhibited_partners,from_rawdf)
        
        if a_targets_indication_df.empty:
          if i_targets_indication_df.empty:
              print('Drug has not partners.  Aborting execution')
              return None,None
          else:
            cols2merge = [c for c in i_targets_indication_df.columns if c not in ['ObjType','URN',self.__temp_id_col__]]
            i_targets_indication_df = i_targets_indication_df.dfcopy(cols2merge)
            combined_df = drug_indication_df.merge_df(i_targets_indication_df,on='Name')
            return super().make_count_df(combined_df,with_name=with_name),i_col2rank
        elif i_targets_indication_df.empty:
          cols2merge = [c for c in a_targets_indication_df.columns if c not in ['ObjType','URN',self.__temp_id_col__]]
          a_targets_indication_df = a_targets_indication_df.dfcopy(cols2merge)
          combined_df = drug_indication_df.merge_df(a_targets_indication_df,on='Name')
          return super().make_count_df(combined_df,with_name=with_name),a_col2rank
        else: # case when both dfs were created
          assert(len(a_targets_indication_df) == len(i_targets_indication_df))
          assert(len(a_targets_indication_df.columns) == len(i_targets_indication_df.columns))
          assert(a_targets_indication_df['Name'].to_list() == i_targets_indication_df['Name'].to_list())

          ai_indication_df = df()
          colnames = a_targets_indication_df.columns.to_list()
          i = 0
          while (i < len(colnames)):
              colname = colnames[i]
              if colname.startswith(self._col_name_prefix): # column is weighted_refcount_column
                  ai_indication_df[colname] = a_targets_indication_df[colname]+i_targets_indication_df[colname]
                  refcount_column = colnames[i+1]
                  ai_indication_df[refcount_column] = a_targets_indication_df[refcount_column]+i_targets_indication_df[refcount_column]
                  linked_concept_count_col = colnames[i+2]
                  ai_indication_df[linked_concept_count_col] = a_targets_indication_df[linked_concept_count_col]+i_targets_indication_df[linked_concept_count_col]
                  concept_count_col = colnames[i+3]
                  ai_indication_df[concept_count_col] = a_targets_indication_df[concept_count_col]+i_targets_indication_df[concept_count_col]
                  i += 4
              else:
                  if colname not in ['ObjType','URN',self.__temp_id_col__]:
                      # to avoid columns duplication by merge_df. 'ObjType','URN' come from load_df
                      ai_indication_df[colname] = a_targets_indication_df[colname]
                  i += 1

          combined_df = drug_indication_df.merge_df(ai_indication_df,on='Name')
          assert(self.__temp_id_col__ in combined_df.columns) 
          real_col2rank = {k:max(v,i_col2rank[k]) for k,v in a_col2rank.items()}
          return super().make_count_df(combined_df,with_name=with_name),real_col2rank
    

    def perform_semantic_search(self):
        '''
        input:
            self.raw_data[INDICATION_COUNT]
            self.raw_data[TOXICITY_COUNT]
        output:
            self.raw_data[INDICATIONS]
            self.raw_data[TOXICITIES]
        '''
        start_time = time.time()
        indication_df,ind_col2rank = self.make_count_df(self.raw_data[INDICATION_COUNT],INDICATIONS,effect_on_indication='negative')
        self.add2raw(indication_df)

        self.col2rank.clear()
        toxicity_df,tox_col2rank  = self.make_count_df(self.raw_data[TOXICITY_COUNT],TOXICITIES,effect_on_indication='positive')
        self.add2raw(toxicity_df)

        print("%s repurposing using %s as targets was done in %s" % 
        (self.params['input_compound'],self.target_names_str(),execution_time(start_time)))
        return ind_col2rank, tox_col2rank
        

    def add_graph_bibliography(self,suffix:str, _4diseases:list[PSObject]):
        """
        adds:
            df with PS_BIBLIOGRAPHY-suffix name to self.report_pandas
        """
        input_targets = set(self.drugs)|set(self.SimilarDrugs)|self.activated_targets|self.inhibited_targets|self.activated_partners|self.inhibited_partners
        disease_neighborhood =  self.Graph.get_subgraph(input_targets, _4diseases)
        super().add_graph_bibliography(suffix,add_graph=disease_neighborhood)


    def Reaxys_targets(self,minpX = 6.0):
        '''
        output:
            drug to targets ResnetGraph
        '''
        def map2prop(targets2map:list[dict],prop_name):
            prop_values = set()
            [prop_values.update(t[prop_name]) for t in targets2map if prop_name in t.keys()]
            prop2psobj,_ = self.map_props2objs(list(prop_values),[prop_name],
                                        case_insensitive=True,only_objtypes=PROTEIN_TYPES)

            targets2map2 = [t for t in targets2map if set(t.get(prop_name,[])).isdisjoint(prop2psobj)]
            return prop2psobj, targets2map2

        rxapi = Reaxys_API()
        rxapi.OpenSession(self.APIconfig)
        
        all_targets = set()
        all_rels = list()
        for drug in self.drugs:
            reaxys_reltargets = rxapi.getTargets(drug.name())
            targets2map = [t[2] for t in reaxys_reltargets]
            uprot2psobj, targets2map = map2prop(targets2map,'GenBank ID')
            name2psobj, targets2map = map2prop(targets2map,'Name')
            alias2psobj, targets2map = map2prop(targets2map,'Alias')
            for refs,relprops,target in reaxys_reltargets:
                mapped_psobjs = set([g for k,v in uprot2psobj.items() if k in target.get('GenBank ID',[]) for g in v])
                mapped_psobjs.update([n for k,v in name2psobj.items() if k in target.get('Name',[]) for n in v])
                mapped_psobjs.update([a for k,v in alias2psobj.items() if k in target.get('Alias',[]) for a in v])
                relprops[OBJECT_TYPE] = ['DirectRegulation']
                all_targets.update(mapped_psobjs)
                all_rels += [PSRelation.make_rel(drug,o,relprops,refs) for o in mapped_psobjs]

        drugs2targets = ResnetGraph()
        drugs2targets.name = 'Reaxys drug targets'
        drugs2targets.add_psobjs(self.drugs)
        drugs2targets.add_psobjs(all_targets)


        [drugs2targets.add_rel(rel) for rel in all_rels] # need to add rels on by one here to merge their references, pX and other attributes
        drugs2targets = drugs2targets.make_simple()

        [drugs2targets[r][t][urn]['relation'].set_affinity() for r,t,urn in drugs2targets.edges(keys=True)]
        # filtering out weak affinity targets:
        edges2remove = [(r,t,urn,rel) for r,t,urn,rel in drugs2targets.edges.data('relation',keys=True) if rel._affinity()<minpX]
        [drugs2targets.remove_edge(r,t,urn) for r,t,urn,rel in edges2remove]
        drugs2targets.remove_nodes_by_degree(1)

        drugs2targets = self.child_update(drugs2targets,make_new_rels=True)
        drugs2targets = drugs2targets.make_simple()
        
        return drugs2targets
    

    def init_semantic_search(self):
        '''
        input:
            self.DrugIndications
            self.DrugToxicities
        '''
        was_initiated = False
        indication_df = self.load_df(list(self.DrugIndications),
                                         max_childs=self.max_ontology_parent,
                                         max_threads=self.max_threads4ontology)
        if not indication_df.empty:
            indication_df._name_ = INDICATION_COUNT
            self.add2raw(indication_df)
            print(f'Will score {len(indication_df)} indications for {self.drug_names()}')
            was_initiated = True
        else:
            print(f'No idications found for {self.drug_names()}')

        toxicity_df = self.load_df(list(self.DrugToxicities),
                                         max_childs=self.max_ontology_parent,
                                         max_threads=self.max_threads4ontology)
        if not toxicity_df.empty:
            toxicity_df._name_ = TOXICITY_COUNT
            self.add2raw(toxicity_df)
            print(f'Will score {len(toxicity_df)} toxicities for {self.drug_names()}')
            was_initiated = True
        else:
            print(f'No toxicities found for {self.drug_names()}')
        
        return was_initiated


    def set_partners(self):
      '''
      ouput:
          self.activated_partners
          self.inhibited_partners
          self.find_activated_partners_oql
          self.find_inhibited_partners_oql
      '''
      try:
        partner_params = dict(self.params['partners'])
        self.activated_partners,t2pGa = self.params2partners(partner_params,self.activated_targets)
        self.inhibited_partners,t2pGi = self.params2partners(partner_params,self.inhibited_targets)
        t2pG = t2pGa.compose(t2pGi)
      except KeyError:
        partner_classes = ResnetGraph.classes(self.activated_targets)
        t2pG = ResnetGraph()
        if partner_classes:
          for cl in partner_classes:
            partners, t2pGa = self.find_partners(cl)
            self.activated_partners.update(partners)
            t2pG = t2pG.compose(t2pGa)
          
        partner_classes = ResnetGraph.classes(self.inhibited_targets)
        if partner_classes:
          for cl in partner_classes:
            partners, t2pGi = self.find_partners(cl)
            self.inhibited_partners.update(partners)
            t2pG = t2pG.compose(t2pGi)
      return t2pG


    def modulators_effects(self):
        ''' 
        output:
            self.DrugDirectAgonists
            self.DrugDirectAntagonists
            self.DrugIndirectAgonists
            self.DrugIndirectAntagonists
        update:
            self.DrugIndications
            self.DrugToxicities
        '''
        kwargs2findmodulators = {
                        'of_types':['SmallMol'],
                        'with_effect' : 'positive'
                      }

        if self.activated_targets:
            kwargs2findmodulators['on_targets'] = list(self.activated_targets)

            kwargs2findmodulators['with_effect'] = 'positive'
            kwargs2findmodulators['linked_by'] = ['DirectRegulation']
            self.DrugDirectAgonists = set(self.find_modulators(**kwargs2findmodulators))

            kwargs2findmodulators['with_effect'] = 'negative'
            kwargs2findmodulators['linked_by'] = ['DirectRegulation']
            self.DrugDirectAntagonists = set(self.find_modulators(**kwargs2findmodulators))

            kwargs2findmodulators['with_effect'] = 'negative'
            kwargs2findmodulators['linked_by'] = ['Regulation','Expression','MolTransport']           
            self.DrugIndirectAntagonists = set(self.find_modulators(**kwargs2findmodulators))

            kwargs2findmodulators['with_effect'] = 'positive'
            kwargs2findmodulators['linked_by'] = ['Regulation','Expression','MolTransport']
            self.DrugIndirectAgonists = set(self.find_modulators(**kwargs2findmodulators))
        
        if self.inhibited_targets:
            kwargs2findmodulators['on_targets'] = list(self.inhibited_targets)

            kwargs2findmodulators['with_effect'] = 'positive'
            kwargs2findmodulators['linked_by'] = ['Regulation','Expression','MolTransport']
            self.DrugIndirectAntagonists.update(self.find_modulators(**kwargs2findmodulators))

            kwargs2findmodulators['with_effect'] = 'negative'
            kwargs2findmodulators['linked_by'] = ['Regulation','Expression','MolTransport']
            self.DrugIndirectAgonists.update(self.find_modulators(**kwargs2findmodulators))

            kwargs2findmodulators['with_effect'] = 'negative'
            kwargs2findmodulators['linked_by'] = ['DirectRegulation']
            self.DrugDirectAgonists.update(self.find_modulators(**kwargs2findmodulators))

            kwargs2findmodulators['with_effect'] = 'negative'
            kwargs2findmodulators['linked_by'] = ['Regulation','Expression','MolTransport']
            self.DrugIndirectAgonists = self.find_modulators(**kwargs2findmodulators)

        self.DrugIndirectAgonists = self.DrugIndirectAgonists.difference(self.DrugDirectAgonists)
        self.DrugIndirectAntagonists = self.DrugIndirectAntagonists.difference(self.DrugDirectAntagonists)

        intersection_1 = self.DrugDirectAgonists.intersection(self.DrugDirectAntagonists)
        intersection_2 = self.DrugDirectAgonists.intersection(self.DrugIndirectAntagonists)
        intersection_3 = self.DrugDirectAntagonists.intersection(self.DrugIndirectAgonists)
        intersection_4 = self.DrugIndirectAgonists.intersection(self.DrugIndirectAntagonists)
        modulator_conflicts = intersection_1|intersection_2|intersection_3|intersection_4
        # including partners into modulator_conflicts to remove metabolite ligands from the list of drug modulators:
        modulator_conflicts = modulator_conflicts|self.activated_partners|self.inhibited_partners 
        self.DrugDirectAgonists = self.DrugDirectAgonists.difference(modulator_conflicts)
        self.DrugDirectAntagonists = self.DrugDirectAntagonists.difference(modulator_conflicts)
        self.DrugIndirectAgonists = self.DrugIndirectAgonists.difference(modulator_conflicts)
        self.DrugIndirectAntagonists = self.DrugIndirectAntagonists.difference(modulator_conflicts)
 
        if not self._is_strict():
            indication_count = len(self.DrugIndications)
            add2DrugIndications = set(self.find_indications4(self.DrugDirectAgonists))
            toxicity_count = len(self.DrugToxicities)
            add2DrugToxicities = set(self.find_indications4(self.DrugDirectAntagonists))
            conflicts = add2DrugToxicities.intersection(add2DrugIndications)
            self.DrugIndications.update(add2DrugIndications.difference(conflicts))
            self.DrugToxicities.update(add2DrugToxicities.difference(conflicts))

            new_ind_count = len(self.DrugIndications)-indication_count
            new_tox_count = len(self.DrugToxicities) - toxicity_count
            print(f'Added {new_ind_count} indications and {new_tox_count} toxicities from drug agonists')
        return


    def __resolve_conflicts(self,conflicts:set[PSObject],using_neighbors:list[PSObject],effect2become_indication:str)->set[PSObject]:

        def __resolve(conflicting:list[PSRelation], using_rel_type:str):
            only_rels_with_type = [r for r in conflicting if r.objtype() == using_rel_type]
            best_effect = ''
            if only_rels_with_type:
                only_rels_with_type.sort(key=lambda x: x.count_refs(), reverse=True)
                best_effect = only_rels_with_type[0].effect()
            return best_effect if best_effect in ['positive','negative'] else ''

        unresolved_conflicts = set()
        neighbors_uids = ResnetGraph.uids(using_neighbors)
        for conflict in conflicts:
            conflict_rels = list(self.Graph.find_relations(neighbors_uids,[conflict.uid()]))
            best_effect = __resolve(conflict_rels,'Regulation')
            if best_effect == effect2become_indication:
                self.DrugToxicities.discard(conflict)
            elif best_effect:
                self.DrugIndications.discard(conflict)
            else:
                best_effect = __resolve(conflict_rels,'QuantitativeChange')
                if best_effect == effect2become_indication:
                    self.DrugToxicities.discard(conflict)
                elif best_effect:
                    self.DrugIndications.discard(conflict)
                else:
                    unresolved_conflicts.add(conflict)

        return unresolved_conflicts


    def resolve_conflict_indications(self):
        conflicts = self.DrugIndications.intersection(self.DrugToxicities)
        indication_count = len(self.DrugIndications)+len(self.DrugToxicities)-len(conflicts)
        conflict_count = len(conflicts)
        print(f'Resolving {conflict_count} conflicting out of total {indication_count} indications')
    
        conflicts = self.__resolve_conflicts(conflicts,self.drugs,effect2become_indication='negative')
        conflicts = self.__resolve_conflicts(conflicts,self.activated_targets,effect2become_indication='negative')
        conflicts = self.__resolve_conflicts(conflicts,self.inhibited_targets,effect2become_indication='positive')
        conflicts = self.__resolve_conflicts(conflicts,self.activated_partners,effect2become_indication='negative')
        conflicts = self.__resolve_conflicts(conflicts,self.inhibited_partners,effect2become_indication='positive')
        print(f'{conflict_count - len(conflicts)} out of {conflict_count} were resolved')
        if conflicts:
            print(f'{len(conflicts)} unresolved conflicts:')
            [print(c.name()+'\t'+c.urn()) for c in conflicts]
        return


    def target_effects4(self,df_name:str)->set[PSObject]:
        '''
        input:
            self.activated_targets, self.inhibited_targets
            df_name in [INDICATION_COUNT,TOXICITY_COUNT]
        '''
        assert(df_name in [INDICATION_COUNT,TOXICITY_COUNT])
        start_time = time.time()
        all_effects = set()

        target_effect_on_indication = 'negative' if df_name == INDICATION_COUNT else 'positive'
        effects = self._indications4(self.activated_targets,target_effect_on_indication)
        all_effects.update(effects)
        target_effect_on_indication = 'positive' if df_name == INDICATION_COUNT else  'negative'
        effects = self._indications4(self.inhibited_targets,target_effect_on_indication)
        all_effects.update(effects)

        if not self.params['debug']: # avoids long calculations
            effect_on_indication = 'negative' if df_name == INDICATION_COUNT else 'positive'
            effects = self._indications4partners(self.activated_targets,self.activated_partners,effect_on_indication)
            all_effects.update(effects)
            effect_on_indication = 'positive' if df_name == INDICATION_COUNT else 'negative'
            effects = self._indications4partners(self.inhibited_targets,self.inhibited_partners,effect_on_indication)
            all_effects.update(effects)

            effect_on_indication = 'negative' if df_name == INDICATION_COUNT else 'positive'
            effects = self._indications4cells(self.CellSecretingActivatedTargets,self.activated_targets,effect_on_indication)
            all_effects.update(effects)
            effect_on_indication = 'positive' if df_name == INDICATION_COUNT else 'negative'
            effects = self._indications4cells(self.CellSecretingInhibitedTargets,self.inhibited_targets,effect_on_indication)
            all_effects.update(effects)

        eff = 'indications' if df_name == INDICATION_COUNT else 'toxicities'
        print(f"{len(all_effects)} {eff} for {df_name} were found in {execution_time(start_time)}")
        return all_effects
    

    def make_report(self):
        '''
        input:
            self.params['mode_of_action']
            self.params['target_names']
        '''
        self.indications4drugs()
        print('Loading INDICATIONS')
        self.DrugIndications = self.target_effects4(INDICATION_COUNT)
        print('Loading TOXICITIES')
        self.DrugToxicities = self.target_effects4(TOXICITY_COUNT)
        if not self.params['debug']:
            self.modulators_effects()
            self.get_pathway_components(self.activated_targets|self.inhibited_targets,self.activated_partners|self.inhibited_partners)
        self.resolve_conflict_indications()

        tm_other = ThreadPoolExecutor(thread_name_prefix='TMother')
        other_effects_future = tm_other.submit(self.other_effects)
        if self.params['add_bibliography']:
            all_effects = self.DrugIndications|self.DrugToxicities
            indication_names = [n.name() for n in all_effects]
            df4etm = df.from_dict({'Name':indication_names})
            drug_names = self.drug_names()
            df4etm._name_ = f'Indication4 {drug_names}'
            tm_biblio_future = tm_other.submit(self.bibliography,df4etm,drug_names,'Name',[],len(df4etm))
        

        if self.init_semantic_search():
            col2ranks = self.perform_semantic_search()
            ranked_df_names = list()
            for i,df_name in enumerate([INDICATIONS,TOXICITIES]):
                col2rank = col2ranks[i]
                drop_empty_cols = list(self.params.get('drop_empty_columns_from_report',[]))
                ranked_df4report,_ = self.normalize(df_name,df_name,col2rank=col2rank,drop_empty_columns=drop_empty_cols)
                with ThreadPoolExecutor(max_workers=4, thread_name_prefix='Report annotation') as e:
                    ontology_df_future = e.submit(self.add_ontology_df,ranked_df4report)
                    add_parent_future = e.submit(self.id2paths,ranked_df4report)
                    if df_name == INDICATIONS:
                        suffix = '4ind'
                        diseases = self.DrugIndications
                    else:
                        suffix = '4tox'
                        diseases = self.DrugToxicities
                    e.submit(self.add_graph_bibliography,suffix,diseases)
                
                id2paths = add_parent_future.result()
                ranked_df = ranked_df4report.merge_dict(id2paths,'Ontology parents','Name')
                ontology_df = ontology_df_future.result()
                e.shutdown()

                self.add2report(ranked_df)
                self.add2report(ontology_df)

        self.report_pandas[INDICATIONS].tab_format['tab_color'] = 'blue'
        self.report_pandas[TOXICITIES].tab_format['tab_color'] = 'pink'
        other_effects_future.result()
        if self.params['add_bibliography']:
            tm_refs_df = tm_biblio_future.result()
            tm_ref_colname = self.tm_refcount_colname('Name',self.drug_names_str())
            doi_ref_colname = self.tm_doi_colname('Name',self.drug_names_str())
            for ws in ranked_df_names:
                self.report_pandas[ws] = self.report_pandas[ws].merge_df(tm_refs_df,on='Name',columns=[tm_ref_colname,doi_ref_colname])
                #self.report_pandas[ws] = self.report_pandas[ws].move_cols({tm_ref_colname:3})
            self.add_tm_bibliography_df()
        return
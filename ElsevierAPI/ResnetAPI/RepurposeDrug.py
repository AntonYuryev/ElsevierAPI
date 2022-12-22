from .TargetIndications import Indications4targets,OQL,df, DF4ANTAGONISTS,DF4AGONISTS
from .TargetIndications import RANK_SUGGESTED_INDICATIONS, PREDICT_RANK_INDICATIONS,COUNTS,ONTOLOGY_ANALYSIS
import time
INHIBIT = 0 # indication must be inhibited by a drug
ACTIVATE = 1 # indication must be activated by a drug

ANTAGONIST = 0 # drug inhibits its targets
AGONIST = 1 # drug activates its targets
pct = '%'


class RepurposeDrug(Indications4targets):
    pass
    def __init__(self, APIconfig, params={}):
        super().__init__(APIconfig,params)
        self.PageSize = 1000
        self.similar_drug_ids = list()
        

    def set_drug(self):   
        #self.SELECTdrug = 'SELECT Entity WHERE InOntology (SELECT Annotation WHERE Ontology=\'Pathway Studio Ontology\' AND Relationship=\'is-a\') under (SELECT OntologicalNode WHERE (Name,Alias) = (\'{drug}\')) OR (Name,Alias) = (\'{drug}\')'
        #self.SELECTdrug = self.SELECTdrug.format(drug=self.param['input_compound'] )
        self.SELECTdrug = OQL.get_childs([self.param['input_compound']],['Name,Alias'],include_parents=True)
        
        drug_graph = self.load_graph_from_oql(self.SELECTdrug,[],['Connectivity'],get_links=False,add2self=False)
        self.drugs = drug_graph._get_nodes()
        self.drugs.sort(key=lambda x: int(x['Connectivity'][0]), reverse=True)
        
        if self.param['drug_effect'] == INHIBIT and self.param['mode_of_action'] == ANTAGONIST:
            self.param['to_inhibit'] = True #most often case for drugs and their disease indications
        elif self.param['drug_effect'] == ACTIVATE and self.param['mode_of_action'] == ANTAGONIST:
            self.param['to_inhibit'] = False # can be used only with CellProcess
        elif self.param['drug_effect'] == ACTIVATE and self.param['mode_of_action'] == AGONIST:
            self.param['to_inhibit'] = True # can be used only with CellProcess
        elif self.param['drug_effect'] == INHIBIT and self.param['mode_of_action'] == AGONIST:
            self.param['to_inhibit'] = False #find disease indications for drugs that are agonist of their target(s)


    def __get_effect_str(self):
        if self.param['drug_effect'] == INHIBIT: return 'negative'
        if self.param['drug_effect'] == ACTIVATE: return 'positive'
        return 'unknown'


    def input_names(self):
        return [x['Name'][0] for x in self.drugs]


    def needs_clinical_trial(self):
        if 'CellProcess' in self.param['indication_types'] and self.param['drug_effect'] == ACTIVATE: return True
        if 'Disease' in self.param['indication_types'] and self.param['drug_effect'] == INHIBIT: return True
        return False


    def find_drug_indications(self):
        """
        Returns
        -------
        {drug indications ids}

        Loads
        -----
        self.indications4antagonists if self.param['mode_of_action'] == ANTAGONIST
        self.indications4agonists if self.param['mode_of_action'] == AGONIST
        """
        indications2return = set()
        indic_str = OQL.join_with_quotes(self.param['indication_types'])
        # PART I: Finding all possible indications
        if self.needs_clinical_trial():
            REQUEST_NAME = 'Find {drug} indications by clinical trials'.format(drug=self.param['input_compound'])
            OQLquery = 'SELECT Relation WHERE objectType = ClinicalTrial AND NeighborOf ({select_drug}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_types}))'
            ClinicalTrialIndictions = self.process_oql(OQLquery.format(select_drug=self.SELECTdrug, indication_types = indic_str),REQUEST_NAME)
            found_indications = ClinicalTrialIndictions.get_node_ids(self.param['indication_types'])
            indications2return.update(found_indications)
            print('Found %d indications in %s clinical trials' % (len(found_indications), self.param['input_compound']))

        effect = 'negative' if self.param['drug_effect'] == INHIBIT else 'positive'
        REQUEST_NAME = 'Find {drug} all other indications'.format(drug=self.param['input_compound'])
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = {eff} AND NeighborOf({select_drug}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_types}))'
        oql_query = OQLquery.format(eff=effect, select_drug=self.SELECTdrug,indication_types = indic_str)
        LiteratureIndications = self.process_oql(oql_query, REQUEST_NAME)
        found_indications = LiteratureIndications.get_node_ids(self.param['indication_types'])
        indications2return.update(found_indications)
        print('Found %d indications reported in scientific literature for %s' % (len(found_indications), self.param['input_compound']))
        self.drug_ids = self.Graph.get_node_ids(['SmallMol'],['ObjTypeName'])

      #  self.drug_indications4strictmode = set(indications2return)
        if self.param['mode_of_action'] == ANTAGONIST:
            self.indications4antagonists.update(indications2return)
            self.indications4antagonists_strict.update(indications2return)
        else:
            self.indications4agonists.update(indications2return)
            self.indications4agonists_strict.update(indications2return)

        return list(indications2return)

    
    def input_names(self):
        drugs = self.Graph._get_nodes(self.drug_ids)
        return [x['Name'][0] for x in drugs]


    def find_indications4similars(self):
        if not self.param['similars']: return
        indications2return = set()
        indic_str = OQL.join_with_quotes(self.param['indication_types'])
        similar_drugs_str = OQL.join_with_quotes(self.param['similars'])
        select_similar_drugs = f'SELECT Entity WHERE (Name,Alias) = ({similar_drugs_str})'

        if self.needs_clinical_trial():
            REQUEST_NAME = 'Find indications by clinical trials for {similars}'.format(similars=','.join(self.param['similars']))
            OQLquery = 'SELECT Relation WHERE objectType = ClinicalTrial AND NeighborOf ({select_drugs}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_types}))'
            SimilarDrugsClinicalTrials = self.process_oql(OQLquery.format(select_drugs=select_similar_drugs, indication_types=indic_str), REQUEST_NAME)
            found_indications = SimilarDrugsClinicalTrials.get_node_ids(self.param['indication_types'])
            indications2return.update(found_indications)
            print('Found %d indications in clinical trials for drugs similar to %s' % (len(found_indications), self.param['input_compound']))

        effect = 'negative' if self.param['drug_effect'] == INHIBIT else 'positive'
        REQUEST_NAME = 'Find all other indications for {similars}'.format(similars=','.join(self.param['similars']))
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = {eff} AND NeighborOf ({select_drugs}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_types}))'
        OQLquery = OQLquery.format(eff=effect,select_drugs=select_similar_drugs,indication_types=indic_str)
        SimilarDrugsIndications = self.process_oql(OQLquery, REQUEST_NAME)
        self.similar_drug_ids = SimilarDrugsIndications.get_node_ids(['SmallMol'])
        found_indications = SimilarDrugsIndications.get_node_ids(self.param['indication_types'])
        indications2return.update(found_indications)
        print('Found %d indications reported in scientific literature or drugs similar to %s' % (len(found_indications), self.param['input_compound']))
        #self.drug_indications4strictmode.update(indications2return)

        if self.param['mode_of_action'] == ANTAGONIST:
            self.indications4antagonists.update(indications2return)
            self.indications4antagonists_strict.update(indications2return)
        else:
            self.indications4agonists.update(indications2return)
            self.indications4agonists_strict.update(indications2return)

        return list(indications2return)


    def clean_indications(self):
        if self.param['drug_effect'] == ACTIVATE:
            # many clinical trials report their indications as toxicities testing if patient condition is not worsening after drug treatment
            if 'Disease' in self.param['indication_types']:
                indications_ids = self.indications4agonists|self.indications4antagonists
                before_clean_indications_count = len(indications_ids)
                indication_child_ids = self.get_children(indications_ids)
                all_indications = indications_ids | set(indication_child_ids)

                all_drug_ids = self.drug_ids+self.similar_drug_ids
                oql_query = 'SELECT Entity WHERE objectType = Disease AND id = ({ids1}) AND Connected by (SELECT Relation WHERE objectType = ClinicalTrial)  \
                    to (SELECT Entity WHERE id = ({ids2}))'
                clinical_trial_graph = self.iterate_oql2(oql_query,all_indications,all_drug_ids)
                clinical_trial_indication_ids = clinical_trial_graph.get_node_ids(['Disease'])

                parents_with_clinical_trial_children = {p for p,ch in self.ID2Children.items() if not set(ch).isdisjoint(clinical_trial_indication_ids)}
                clinical_trial_indication_ids += parents_with_clinical_trial_children
                for p in parents_with_clinical_trial_children:
                    clinical_trial_indication_ids += self.ID2Children[p]

                self.indications4agonists = self.indications4agonists.difference(clinical_trial_indication_ids)
                self.indications4antagonists = self.indications4antagonists.difference(clinical_trial_indication_ids)
                self.indications4agonists_strict = self.indications4agonists_strict.difference(clinical_trial_indication_ids)
                self.indications4antagonists_strict = self.indications4antagonists_strict.difference(clinical_trial_indication_ids)
                after_clean_indications_count = len(self.indications4agonists|self.indications4antagonists)
                print('%d toxicities were removed because they are indications for %s in clinical trials' % 
                        (before_clean_indications_count-after_clean_indications_count,self.input_names()))
            

    def score_drug_semantics(self):
        if self.param['mode_of_action'] == AGONIST:
            my_df = self.raw_data[DF4AGONISTS]
        else:
            my_df = self.raw_data[DF4ANTAGONISTS]

        colname = self.param['input_compound'] + ' clinical trials'
        self.set_how2connect(['ClinicalTrial'],[],'')
        linked_entities_count, linked_entity_ids, my_df = self.link2concept(colname, self.drug_ids,my_df)
        print('%s has clinical trials for %d indications' % (self.param['input_compound'], linked_entities_count))

        if self.param['drug_effect'] == INHIBIT:
            colname = 'Inhibited by ' + self.param['input_compound']
            self.set_how2connect(['Regulation'],['negative'],'',['Regulation'])
            regulated = 'inhibited'
        else:
            colname = 'Activated by ' + self.param['input_compound']
            self.set_how2connect(['Regulation'],['positive'],'',['Regulation'])
            regulated = 'activated'

        linked_entities_count, linked_entity_ids, my_df = self.link2concept(colname, self.drug_ids,my_df)
        print('%d indications %s by %s' % (linked_entities_count, regulated, self.param['input_compound']))
        self.drugcolname = self._col_name_prefix+colname
        
        if hasattr(self,'similar_drug_ids'):
            colname = 'similar drugs clinical trials'
            self.set_how2connect(['ClinicalTrial'],[],'')
            linked_entities_count, linked_entity_ids, my_df = self.link2concept(colname, self.similar_drug_ids,my_df)
            print('%s has clinical trials for %s indications' % (','.join(self.param['similars']), linked_entities_count))

            if self.param['drug_effect'] == INHIBIT:
                colname = 'Inhibited by similar drugs'
                self.set_how2connect(['Regulation'],['negative'],'')
            else:
                colname = 'Activated by similar drugs'
                self.set_how2connect(['Regulation'],['positive'],'',['Regulation'])
            linked_entities_count, linked_entity_ids, my_df = self.link2concept(colname,self.similar_drug_ids,my_df)
            
            self.add2raw(my_df)
            print('%d indications %s by %s' % (linked_entities_count, regulated, ','.join(self.param['similars'])))
            return my_df
        

    def report_name(self):
        indics = OQL.join_with_quotes(self.param['indication_types'])
        rep_pred =  'suggested ' if self.param['strict_mode'] else 'suggested,predicted '
        regulate = ' activated by ' if self.param['drug_effect'] == ACTIVATE else ' inhibited by '
        return str(rep_pred+indics+regulate+self.param['input_compound'])


    def report_path(self,extension='.xlsx'):
        return self.data_dir+self.report_name()+extension


    def add_drug_indication_refs(self):
        print('Printing clinical trials')
        clin_trial_graph = self.Graph.subgraph_by_relprops(['ClinicalTrial'])
        ct_pd = self.Graph.snippets2df(clin_trial_graph)
        ct_pd._name_ ='clin. trials'
        self.add2report(ct_pd)

        
        print('Printing research articles')
        research_graph = self.Graph.subgraph_by_relprops(['Regulation'])
        effect_str = self.__get_effect_str()
        research_graph = self.Graph.subgraph_by_relprops([effect_str],['Effect'])
        research_ref_pd = self.Graph.snippets2df(research_graph)

        research_ref_pd._name_ ='articles'
        self.add2report(research_ref_pd)
  

    def _worksheet_prefix(self):
        indications = ','.join(self.param['indication_types'])
        targets = ','.join(self.param['target_names'])
        ef = 'Act.' if self.param['drug_effect'] == ACTIVATE else 'Inh.'
        return ef+indications+'-'+targets


    def load_drug_indications(self):
        self.flush_dump_files()
        start_time = time.time()
        self.set_drug() 
        self.find_drug_indications()
        self.find_indications4similars()
        self.clean_indications()
        print("Drug indications were loaded in %s" % self.execution_time(start_time))


    def other_effects(self):
        drug = self.param['input_compound']
        indication_str = ','.join(self.param['indication_types'])
        REQUEST_NAME = 'Find {inds} modulated by {drug} with unknown effect'.format(inds=indication_str,drug=drug)
        select_indications = 'SELECT Entity WHERE objectType = ({})'.format(indication_str)
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = unknown AND NeighborOf({select_target}) AND NeighborOf ({indications})'
        OQLquery = OQLquery.format(select_target=self.SELECTdrug, indications=select_indications)
        to_return = self.process_oql(OQLquery,REQUEST_NAME)
        return to_return


    def clear(self):
        super().clear()
        self.RefCountPandas = df()
  

    def _get_report_name(self):
        indics = ','.join(self.param['indication_types'])
        effect = ' inhibited by ' if self.param['drug_effect'] == ANTAGONIST else ' activated by '
        rep_pred = 'suggested ' if self.param['strict_mode'] == RANK_SUGGESTED_INDICATIONS else 'suggested,predicted ' 
        return self.param['data_dir']+rep_pred+ indics+'s'+effect+ self.param['input_compound']


    def perform_semantic_search(self):
        start_time = time.time()
        if self.init_semantic_search():
            my_df = self.score_drug_semantics()
            target_effect_on_indication = 'positive' if self.param['to_inhibit'] else 'negative'
            self.semantic_score(my_df._name_,target_effect_on_indication)

            trgts = ','.join(self.param['target_names'])
            print("%s repurposing using %s as targets was done in %s" % 
                (self.param['input_compound'],trgts,self.execution_time(start_time)))
            return my_df._name_
        else:
            return ''


    def make_report(self,for_targets:list):
        self.load_drug_indications()
        self.param['target_names'] = for_targets

        if self.param['strict_mode'] == PREDICT_RANK_INDICATIONS:
            self.load_indications4targets()
        else:
            self.set_targets()
            self.set_partners()
            self.get_pathway_componets()

        count_df_name = self.perform_semantic_search()
        if count_df_name:
            worksheet_prefix = self._worksheet_prefix()
            norm_df_name = worksheet_prefix+'ranked_count'
            try:
                drop_empty_cols = self.param['drop_empty_columns_from_report']
            except KeyError:
                drop_empty_cols = False

            self.normalize(count_df_name,norm_df_name,drop_empty_columns=drop_empty_cols)
            self.add_parent_column(norm_df_name)
            self.add_ontology_df(norm_df_name)
            
            self.add_ps_bibliography()
            self.add_etm_refs(norm_df_name,self.input_names())
            biblio_df_name = self.add_etm_bibliography(worksheet_prefix)
            return norm_df_name, biblio_df_name
        else:
            return '',''

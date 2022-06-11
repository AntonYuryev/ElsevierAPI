from .TargetIndications import TargetIndications,OQL,df
from .TargetIndications import RANK_SUGGESTED_INDICATIONS, PREDICT_RANK_INDICATIONS
import time
INHIBIT = 0 # indication must be inhibited by a drug
ACTIVATE = 1 # indication must be activated by a drug

ANTAGONIST = 0 # drug inhibits its targets
AGONIST = 1 # drug activates its targets
pct = '%'


class RepurposeDrug(TargetIndications):
    pass
    def __init__(self, APIconfig, params={}):
        super().__init__(APIconfig,params)
        self.PageSize = 1000

    def set_drug(self):   
        self.SELECTdrug = 'SELECT Entity WHERE InOntology (SELECT Annotation WHERE Ontology=\'Pathway Studio Ontology\' AND Relationship=\'is-a\') under (SELECT OntologicalNode WHERE (Name,Alias) = (\'{drug}\')) OR (Name,Alias) = (\'{drug}\')'
        self.SELECTdrug = self.SELECTdrug.format(drug=self.param['input_compound'] )

        drug_graph = self.load_graph_from_oql(self.SELECTdrug,entity_props=['Connectivity'],get_links=False,add2self=False)
        drugs = drug_graph._get_nodes()
        drugs.sort(key=lambda x: int(x['Connectivity'][0]), reverse=True)
        self.input_names = [x['Name'][0] for x in drugs]

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


    def needs_clinical_trial(self):
        if 'CellProcess' in self.param['indication_types'] and self.param['drug_effect'] == ACTIVATE: return True
        if 'Disease' in self.param['indication_types'] and self.param['drug_effect'] == INHIBIT: return True
        return False


    def find_drug_indications(self):
        """
        returns set of drug indications ids
        """
        indications2return = set()
        indic_str = OQL.join_with_quotes(self.param['indication_types'])
        # PART I: Finding all possible indications
        if self.needs_clinical_trial():
            REQUEST_NAME = 'Find {drug} indications by clinical trials'.format(drug=self.param['input_compound'])
            OQLquery = 'SELECT Relation WHERE objectType = ClinicalTrial AND NeighborOf ({select_drug}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_types}))'
            ClinicalTrialIndictions = self.process_oql(OQLquery.format(select_drug=self.SELECTdrug, indication_types = indic_str),REQUEST_NAME)
            found_indications = ClinicalTrialIndictions.get_node_ids(self.indication_types)
            indications2return.update(found_indications)
            print('Found %d indications in %s clinical trials' % (len(found_indications), self.param['input_compound']))

        effect = 'negative' if self.param['drug_effect'] == INHIBIT else 'positive'
        REQUEST_NAME = 'Find {drug} all other indications'.format(drug=self.param['input_compound'])
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = {eff} AND NeighborOf({select_drug}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_types}))'
        oql_query = OQLquery.format(eff=effect, select_drug=self.SELECTdrug,indication_types = indic_str)
        LiteratureIndications = self.process_oql(oql_query, REQUEST_NAME)
        found_indications= LiteratureIndications.get_node_ids(self.indication_types)
        indications2return.update(found_indications)
        print('Found %d indications reported in scientific literature for %s' % (len(found_indications), self.param['input_compound']))
        self.drug_ids = self.Graph.get_node_ids(['SmallMol'],['ObjTypeName'])
        self.drug_indications4strictmode = set(indications2return)
        return list(indications2return)


    def find_indications4similars(self):
        if not self.param['similars']: return
        indications2return = set()
        indic_str = OQL.join_with_quotes(self.param['indication_types'])
        select_similar_drugs = 'SELECT Entity WHERE (Name,Alias) = ({drugs})'
        similar_drugs_str = OQL.join_with_quotes(self.param['similars'])
        select_similar_drugs = select_similar_drugs.format(drugs=similar_drugs_str)

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
        self.drug_indications4strictmode.update(indications2return)
        return list(indications2return)


    def score_drug_semantics(self):
        colname = self.param['input_compound'] + ' clinical trials'
        self.set_how2connect(['ClinicalTrial'],[],'')
        linked_entities_count = self.link2concept(colname, self.drug_ids)
        print('%s has clinical trials for %d indications' % (self.param['input_compound'], linked_entities_count))

        if self.param['drug_effect'] == INHIBIT:
            colname = 'Inhibited by ' + self.param['input_compound']
            self.set_how2connect(['Regulation'],['negative'],'')
            regulated = 'inhibited'
        else:
            colname = 'Activated by ' + self.param['input_compound']
            self.set_how2connect(['Regulation'],['positive'],'')
            regulated = 'activated'

        linked_entities_count = self.link2concept(colname, self.drug_ids)
        print('%d indications %s by %s' % (linked_entities_count, regulated, self.param['input_compound']))
        self.drugcolname = self._col_name_prefix+colname

        if hasattr(self,'similar_drug_ids'):
            colname = 'similar drugs clinical trials'
            self.set_how2connect(['ClinicalTrial'],[],'')
            linked_entities_count = self.link2concept(colname, self.similar_drug_ids)
            print('%s has clinical trials for %s indications' % (','.join(self.param['similars']), linked_entities_count))

            if self.param['drug_effect'] == INHIBIT:
                colname = 'Inhibited by similar drugs'
                self.set_how2connect(['Regulation'],['negative'],'')
            else:
                colname = 'Activated by similar drugs'
                self.set_how2connect(['Regulation'],['positive'],'')
    
            linked_entities_count = self.link2concept(colname, self.similar_drug_ids)
            print('%d indications %s by %s' % (linked_entities_count, regulated, ','.join(self.param['similars'])))
        

    def fname_prefix(self):
        indics = OQL.join_with_quotes(self.indication_types)
        rep_pred =  'suggested ' if self.param['strict_mode'] else 'suggested,predicted '
        regulate = ' activated by ' if self.param['drug_effect'] == ACTIVATE else ' inhibited by '
        return rep_pred+indics+regulate+self.param['input_compound']+' '


    def print_drug_indictaion_refs(self):
        #fname_out = self.data_dir+self.fname_prefix()+" clinical trials.tsv"
        print('Printing clinical trials')
        ct_pd = self.refs2pd(self.drugcolname, ['ClinicalTrial'],
                            pubid_types=['NCT ID'],biblio_props=['Title','Start'],print_snippets=True)

       # ct_pd.to_excel(self.excel_writer, sheet_name='Clin. Trials', index=False)
        ct_pd.name ='clin. trials'
        self.add2report(ct_pd)

        effect_str = self.__get_effect_str()
        #fname_out = self.data_dir+self.fname_prefix()+" references.tsv"
        print('Printing research articles')
        research_ref_pd = self.refs2pd(self.drugcolname,['Regulation'],[effect_str],
                            pubid_types=['PMID','DOI'],biblio_props=['Title','PubYear'],print_snippets=True)

        research_ref_pd.name ='articles'
        self.add2report(research_ref_pd)
        #research_ref_pd.to_excel(self.excel_writer, sheet_name='Snippets', index=False)

    def set_strictmode_indications(self):
        self.indications4strictmode = self.drug_indications4strictmode.intersection(self.target_indications4strictmode)

    def _worksheet_prefix(self):
        indications = ','.join(self.indication_types)
        targets = ','.join(self.param['target_names'])
        ef = 'Act.' if self.param['drug_effect'] == ACTIVATE else 'Inh.'
        return ef+indications+'-'+targets

    def load_drug_indications(self):
        self.flush_dump_files()
        start_time = time.time()
        self.set_drug() 
        self.find_drug_indications()
        self.find_indications4similars()
        print("Drug indications were loaded in %s" % self.execution_time(start_time))


    def load_target_indications(self,targets:list):
        start_time = time.time()
        self.param['target_names'] = targets
        self.set_targets()

        self.find_target_indications()
        self.indications4chem_modulators()
        self.indications4partners()
        self.get_pathway_componets()

      #  if self.param['strict_mode']: 
      #      self.set_strictmode_indications()
        print("Target indications were loaded in %s" % self.execution_time(start_time))


    def perform_semantic_search(self):
        start_time = time.time()
        if self.init_semantic_search():
            self.score_drug_semantics()
            self.score_semantics()

            self.normalize_counts(bibliography=self.param['add_bibliography'])
            trgts = ','.join(self.param['target_names'])
            print("%s repurposing using %s as targets was done in %s" % 
                (self.param['input_compound'],trgts,self.execution_time(start_time)))
            return True
        else:
            return False


    def make_report(self):
        start_time = time.time()

        self.set_drug() 
        # if partner_names is empty script will try finding Ligands for Receptor targets and Receptors for Ligand targets in the knowledge graph
        self.find_drug_indications()
        self.find_indications4similars()

        self.set_targets()
        # strict mode ranks only indications suggested in the literarure without predicting new indications
        self.flush_dump_files()

        self.find_target_indications()
        self.indications4chem_modulators()
        self.indications4partners()
        self.get_pathway_componets()

        if self.param['strict_mode']: 
            self.set_strictmode_indications()

        if self.init_semantic_search():
            self.score_drug_semantics()
            self.score_semantics()

            self.normalize_counts(bibliography=True)
            trgts = ','.join(self.param['target_names'])
            print("%s repurposing using %s as targets was done in %s" % 
                (self.param['input_compound'],trgts,self.execution_time(start_time)))
            return True
        else:
            return False


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
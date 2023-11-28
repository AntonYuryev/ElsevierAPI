from .TargetIndications import Indications4targets,OQL,df,DF4ANTAGONISTS,DF4AGONISTS,ResnetGraph
from .TargetIndications import ANTAGONIST,AGONIST,BIBLIO_PROPERTIES,PS_SENTENCE_PROPS
from ElsevierAPI import  execution_time
from concurrent.futures import ThreadPoolExecutor
from .SemanticSearch import RANK
import time

INHIBIT = -1 # indication must be inhibited by a drug
ACTIVATE = 1 # indication must be activated by a drug
ALL_EFFECTS = 0

pct = '%'
class RepurposeDrug(Indications4targets):
    pass
    def __init__(self, *args, **kwargs):
        '''
        APIconfig - args[0]
        '''
        APIconfig = args[0] if args else dict()

        my_kwargs = {
                'input_compound' : '', 
                'similars' : [],
                'drug_effect': INHIBIT, # INHIBIT for indications, ACTIVATE for toxicitites
                'mode_of_action': ANTAGONIST,
                'indication_types': ['Disease'], #['CellProcess','Virus']
                }
        my_kwargs.update(kwargs)

        super().__init__(APIconfig,**my_kwargs)
     #   self.PageSize = 1000
        self.similar_drugs = list()
        self.drug_indications = set()
        self.directly_inh = set() 
        self.directly_act = set() 
        self.indirectly_inh = set() 
        self.indirectly_act = set() 
        

    def target2indication_effect(self):
        if self.params['drug_effect'] == INHIBIT and self.params['mode_of_action'] == ANTAGONIST:
            return str('positive') #most often case for drugs and their disease indications
        elif self.params['drug_effect'] == ACTIVATE and self.params['mode_of_action'] == ANTAGONIST:
            return str('negative') # can be used only with CellProcess
        elif self.params['drug_effect'] == ACTIVATE and self.params['mode_of_action'] == AGONIST:
            return str('positive') # can be used only with CellProcess
        elif self.params['drug_effect'] == INHIBIT and self.params['mode_of_action'] == AGONIST:
            return str('negative') #find disease indications for drugs that are agonist of their target(s)
        else:
            return ''



    def set_drug(self):
        self.SELECTdrug = OQL.get_childs([self.params['input_compound']],['Name,Alias'],include_parents=True)
        
        drug_graph = self.load_graph_from_oql(self.SELECTdrug,[],['Connectivity'],get_links=False,add2self=False)
        self.drugs = drug_graph._get_nodes()
        self.drugs.sort(key=lambda x: int(x['Connectivity'][0]), reverse=True)


    def __effect_str(self):
        if self.params['drug_effect'] == INHIBIT: return 'negative'
        if self.params['drug_effect'] == ACTIVATE: return 'positive'
        return 'unknown'


    def needs_clinical_trial(self):
        if 'CellProcess' in self.params['indication_types'] and self.params['drug_effect'] == ACTIVATE: return True
        if 'Disease' in self.params['indication_types'] and self.params['drug_effect'] == INHIBIT: return True
        return False


    def find_drug_indications(self):
        """
        Returns
        -------
        {drug indications ids}

        Loads
        -----
        self.indications4antagonists if self.params['mode_of_action'] == ANTAGONIST
        self.indications4agonists if self.params['mode_of_action'] == AGONIST
        """
        indications2return = set()
        indic_str = OQL.join_with_quotes(self.params['indication_types'])
        # PART I: Finding all possible indications
        if self.needs_clinical_trial():
            REQUEST_NAME = 'Find {drug} indications by clinical trials'.format(drug=self.params['input_compound'])
            OQLquery = 'SELECT Relation WHERE objectType = ClinicalTrial AND NeighborOf ({select_drug}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_types}))'
            ClinicalTrialIndictions = self.process_oql(OQLquery.format(select_drug=self.SELECTdrug, indication_types = indic_str),REQUEST_NAME)
            if isinstance(ClinicalTrialIndictions,ResnetGraph):
                found_indications = ClinicalTrialIndictions.psobjs_with(only_with_values=self.params['indication_types'])
                indications2return.update(found_indications)
                print('Found %d indications in %s clinical trials' % (len(found_indications), self.params['input_compound']))

        effect = 'negative' if self.params['drug_effect'] == INHIBIT else 'positive'
        REQUEST_NAME = 'Find {drug} all other indications'.format(drug=self.params['input_compound'])
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = {eff} AND NeighborOf({select_drug}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_types}))'
        oql_query = OQLquery.format(eff=effect, select_drug=self.SELECTdrug,indication_types = indic_str)
        LiteratureIndications = self.process_oql(oql_query, REQUEST_NAME)
        if isinstance(LiteratureIndications,ResnetGraph):
            found_indications = LiteratureIndications.psobjs_with(only_with_values=self.params['indication_types'])
            indications2return.update(found_indications)
            print('Found %d indications reported in scientific literature for %s' % (len(found_indications), self.params['input_compound']))
    
        self.drug_indications.update(indications2return)
        return list(indications2return)

    
    def input_drug_names(self):
        return [x.name() for x in self.drugs]
    

    def __etm_ref_column_name(self):
        return self.etm_counter._etm_ref_column_name('Name',self.input_drug_names())
    
    def __etm_doi_column_name(self):
        return self.etm_counter._etm_doi_column_name('Name',self.input_drug_names())


    def find_indications4similars(self):
        if not self.params['similars']: return
        indications2return = set()
        indic_str = OQL.join_with_quotes(self.params['indication_types'])
        similar_drugs_str = OQL.join_with_quotes(self.params['similars'])
        select_similar_drugs = f'SELECT Entity WHERE (Name,Alias) = ({similar_drugs_str})'

        if self.needs_clinical_trial():
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
            self.similar_drugs = SimilarDrugsIndications.psobjs_with(only_with_values=['SmallMol'])
            found_indications = SimilarDrugsIndications.psobjs_with(only_with_values=self.params['indication_types'])
            indications2return.update(found_indications)
            print('Found %d indications reported in scientific literature or drugs similar to %s' % (len(found_indications), self.params['input_compound']))

        self.drug_indications.update(indications2return)
        return list(indications2return)


    def clean_indications(self):
        if self.params['drug_effect'] == ACTIVATE: # == report is for toxicities
# clinical trials report their indications as toxicities testing 
# to check if patient condition is not worsening after drug treatment
            if 'Disease' in self.params['indication_types']:
                #indications = self.drug_indications
                before_clean_indications_count = len(self.drug_indications)
                indication_childs, parents_with_children = self.load_children4(list(self.drug_indications))
                all_indications = self.drug_indications | indication_childs

                all_drugs = self.drugs+self.similar_drugs
                oql_query = 'SELECT Entity WHERE objectType = Disease AND id = ({ids1}) AND Connected by (SELECT Relation WHERE objectType = ClinicalTrial)  \
                    to (SELECT Entity WHERE id = ({ids2}))'
                all_indications_dbid = ResnetGraph.dbids(list(all_indications))
                all_drugs_dbids = ResnetGraph.dbids(all_drugs)
                clinical_trial_graph = self.iterate_oql2(oql_query,set(all_indications_dbid),set(all_drugs_dbids))
                clinical_trial_indications = clinical_trial_graph._psobjs_with('Disease','ObjTypeName')
                clintrial_indication_parents, ctuid2parent = self.Graph.find_parents4(clinical_trial_indications)
                clinical_trial_indications += list(clintrial_indication_parents)

                #removing clinical_trial_indications from the list of indications for target drugs
                self.drug_indications = self.drug_indications.difference(clinical_trial_indications)
                after_clean_indications_count = len(self.drug_indications)
                print('%d toxicities were removed because they are indications for %s in clinical trials' % 
                        (before_clean_indications_count-after_clean_indications_count,self.input_drug_names()))
            

    def semscore4drugs(self):
        if self._4agonist():
            indication_df = self.raw_data[DF4AGONISTS]
        else:
            indication_df = self.raw_data[DF4ANTAGONISTS]

        #self.etm_refs2df(indication_df,self.input_drug_names(),'Name',[]) # for speed testing etm
        colname = self.params['input_compound'] + ' clinical trials'
        how2connect = self.set_how2connect (['ClinicalTrial'],[],'')
        linked_entities_count, linked_entity_ids, indication_df = self.link2concept(colname, self.drugs,indication_df,how2connect)
        print('%s has clinical trials for %d indications' % (self.params['input_compound'], linked_entities_count))

        if self.params['drug_effect'] == INHIBIT:
            colname = 'Inhibited by ' + self.params['input_compound']
            how2connect = self.set_how2connect (['Regulation'],['negative'],'',['Regulation'])
            regulated = 'inhibited'
        else:
            colname = 'Activated by ' + self.params['input_compound']
            how2connect = self.set_how2connect (['Regulation'],['positive'],'',['Regulation'])
            regulated = 'activated'

        linked_entities_count, linked_entity_ids, indication_df = self.link2concept(colname, self.drugs,indication_df,how2connect)
        print('%d indications %s by %s' % (linked_entities_count, regulated, self.params['input_compound']))
        self.drugcolname = self._col_name_prefix+colname
        
        if hasattr(self,'similar_drugs'):
            colname = 'similar drugs clinical trials'
            how2connect = self.set_how2connect (['ClinicalTrial'],[],'')
            linked_entities_count, linked_entity_ids, indication_df = self.link2concept(colname, self.similar_drugs,indication_df,how2connect)
            print('%s has clinical trials for %s indications' % (','.join(self.params['similars']), linked_entities_count))

            if self.params['drug_effect'] == INHIBIT:
                colname = 'Inhibited by similar drugs'
                how2connect = self.set_how2connect (['Regulation'],['negative'],'')
            else:
                colname = 'Activated by similar drugs'
                how2connect = self.set_how2connect (['Regulation'],['positive'],'',['Regulation'])
            linked_entities_count, linked_entity_ids, indication_df = self.link2concept(colname,self.similar_drugs,indication_df,how2connect)
            print('%d indications %s by %s' % (linked_entities_count, regulated, ','.join(self.params['similars'])))

        self.add2raw(indication_df)
        return indication_df
        

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


    def add_drug_indication_refs(self):
        print('Printing clinical trials')
        clin_trial_graph = self.Graph.subgraph_by_relprops(['ClinicalTrial'])
        ct_pd = clin_trial_graph.snippets2df()
        ct_pd._name_ ='clin. trials'
        self.add2report(ct_pd)

        
        print('Printing research articles')
        research_graph = self.Graph.subgraph_by_relprops(['Regulation'])
        effect_str = self.__effect_str()
        research_graph = self.Graph.subgraph_by_relprops([effect_str],['Effect'])
        research_ref_pd = research_graph.snippets2df()

        research_ref_pd._name_ ='articles'
        self.add2report(research_ref_pd)
  

    def __ws_suffix(self,with_targets_names=False):
        indications = ','.join(self.params['indication_types'])
        targets = ','.join(self.params['target_names'])
        ef = 'Act.' if self.params['drug_effect'] == ACTIVATE else 'Inh.'
        return targets[:10]+'-'+ef+indications[:3] if with_targets_names else ef+indications[:3]


    def load_drug_indications(self):
        self.flush_dump_files()
        start_time = time.time()
        self.set_drug()
        self.find_drug_indications()
        #return #uncomment for fast downstream testing
        self.find_indications4similars()
        self.clean_indications()
        print("Drug indications were loaded in %s" % execution_time(start_time))


    def other_effects(self):
        drug = self.params['input_compound']
        indication_str = ','.join(self.params['indication_types'])
        REQUEST_NAME = f'Find {indication_str} modulated by {drug} with unknown effect'
        old_rel_prop = list(self.relProps)
        self.add_rel_props(PS_SENTENCE_PROPS)
        select_indications = f'SELECT Entity WHERE objectType = ({indication_str})'
        OQLquery = f'SELECT Relation WHERE objectType = Regulation AND Effect = unknown AND \
            NeighborOf({self.SELECTdrug}) AND NeighborOf ({select_indications})'
        to_return = self.process_oql(OQLquery,REQUEST_NAME)
        self.relProps = old_rel_prop
        return to_return


    def perform_semantic_search(self):
        start_time = time.time()
        if self.init_semantic_search():
            indication_df = self.semscore4drugs()

            # for testing parent path only:
     #       id2child_objs = self.entities(indication_df)
     #       self.ontopaths2(id2child_objs,3)

            target_effect_on_indication = self.target2indication_effect()
            self.semscore4targets(indication_df._name_,target_effect_on_indication)

            trgts = ','.join(self.params['target_names'])
            print("%s repurposing using %s as targets was done in %s" % 
                (self.params['input_compound'],trgts,execution_time(start_time)))
            return indication_df._name_
        else:
            return ''


    def make_report(self,for_targets_with_names:list=[]):
        self.load_drug_indications()
        if for_targets_with_names:
            self.params['target_names'] = for_targets_with_names
        self.load_indications4targets()

        if self._is_strict():
            if self._4agonist():
                self.indications4agonists_strict = self.drug_indications.intersection(self.indications4agonists_strict)
                self.indications4antagonists_strict.clear()
            else:
                self.indications4antagonists_strict = self.drug_indications.intersection(self.indications4antagonists_strict)
                self.indications4agonists_strict.clear()
        else:
            if self._4agonist():
                self.indications4agonists = self.drug_indications|self.indications4agonists
                self.indications4antagonists.clear()
            else:
                self.indications4antagonists = self.drug_indications|self.indications4antagonists
                self.indications4agonists.clear()

        count_df_name = self.perform_semantic_search()
        if count_df_name:
            try:
                drop_empty_cols = self.params['drop_empty_columns_from_report']
            except KeyError:
                drop_empty_cols = []

            ws_suffix = self.__ws_suffix()
            norm_df_name = 'rnkd '+ws_suffix
            self.normalize(count_df_name,norm_df_name,drop_empty_columns=drop_empty_cols)

            with ThreadPoolExecutor(max_workers=3, thread_name_prefix='Report annotation') as e:
                etm_biblio_future = e.submit(self.add_etm_refs, norm_df_name,self.input_drug_names(),'Name',[],False)
                ontology_df_future = e.submit(self.add_ontology_df,norm_df_name)
                add_parent_future = e.submit(self.id2paths,norm_df_name)
                ps_biblio_future = e.submit(self.add_ps_bibliography,ws_suffix)
                
            norm_df = df.copy_df(self.report_pandas[norm_df_name])
            id2paths = add_parent_future.result()
            norm_df = norm_df.merge_dict(id2paths,'Ontology parents','Name')

            ontology_df = ontology_df_future.result()
            ps_biblio_future.result()

            etm_refs_df = etm_biblio_future.result()
            norm_df = norm_df.merge_df(etm_refs_df,on='Name',columns=[self.__etm_ref_column_name(),self.__etm_doi_column_name()])
            
            columns = norm_df.columns.to_list()
            columns.remove('URN')
            columns.append('URN')
            columns.remove(RANK)
            columns.insert(1,RANK)
            norm_df = norm_df.reorder(columns)
            self.add2report(norm_df)

            biblio_df_name = self.add_etm_bibliography(ws_suffix)
            return norm_df_name, biblio_df_name, ontology_df._name_
        else:
            return '','',''

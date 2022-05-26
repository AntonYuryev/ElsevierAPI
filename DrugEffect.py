from ElsevierAPI import load_api_config
from TargetIndications import TargetIndications,OQL,SENTENCE_PROPS
import time
import pandas as pd

INHIBIT = 0 # indication must be inhibited by a drug
ACTIVATE = 1 # indication must be activated by a drug

ANTAGONIST = 0 # drug inhibits its targets
AGONIST = 1 # drug activates its targets

RANK_SUGGESTED_INDICATIONS = True 
# RANK_SUGGESTED_INDICATIONS indications suggested in the lietarure are ranked by amount of supporting evidence
PREDICT_RANK_INDICATIONS = False
# also predict indication using biology of the drug target(s)

pct = '%'
DISEASE_MAP2ONTOLOGY = [
                'psychiatric disorder',
                'neurological disorder',
                'musculoskeletal disorder',
                'digestive system disease',
                'endocrine disorder',
                'urogenital disease',
                'skin and connective tissue disease',
                'cardiovascular disease',
                'head and neck disease',
                'respiratory disease',
                'genetic and familial disorder',
                'sclerosis',
                'metabolic disorder',
                'neoplasm',
                'inflammation',
                'degeneration',
                'infection',
                'complications',
                'toxicity',
                'wounds and injuries',
                'pain',
                'nutritional and metabolic diseases',
                'body weight disorder',
                'motor dysfunction',
                'neurocognitive disorder',
                'immunopathology',
                'nausea and vomiting',
                'infertility',
                'sleep disorder',
                'appetite disturbance'
                ]

CELLPROCESS_MAP2ONTOLOGY =[
                            'behavior',
                            'cell proliferation',
                            'apoptosis',
                            'memory',
                            'appetite',
                            'cell death',
                            'locomotion',
                            'cognition',
                            'synaptic transmission',
                            'neurotransmission',
                            'vasodilation',
                            'perception of pain',
                            'long-term synaptic depression',
                            'synaptic plasticity',
                            'learning',
                            'working memory',
                            'cell differentiation',
                            'neuronal activity',
                            'sleep',
                            'long-term synaptic potentiation'
                        ]

class RepurposeDrug(TargetIndications):
    pass
    def __init__(self, APIconfig, params={}):
        super().__init__(APIconfig,params)
        self.PageSize = 1000

    #VVdrug_name:str, similar_drugs:list, mode=INHIBIT, indication_types=None, type=ANTAGONIST
    def set_drug(self):   
        self.SELECTdrug = 'SELECT Entity WHERE InOntology (SELECT Annotation WHERE Ontology=\'Pathway Studio Ontology\' AND Relationship=\'is-a\') under (SELECT OntologicalNode WHERE (Name,Alias) = (\'{drug}\')) OR (Name,Alias) = (\'{drug}\')'
        self.SELECTdrug = self.SELECTdrug.format(drug=self.param['input_compound'] )

        drug_graph = self.load_graph_from_oql(self.SELECTdrug,entity_props=['Connectivity'],get_links=False,add2self=False)
        drugs = drug_graph._get_nodes()
        drugs.sort(key=lambda x: int(x['Connectivity'][0]), reverse=True)
        self.input_names = [x['Name'][0] for x in drugs]

        #self.required_effect_on_indications = self.param['']
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
            found_indications = ClinicalTrialIndictions.get_entity_ids(self.indication_types)
            indications2return.update(found_indications)
            print('Found %d indications in %s clinical trials' % (len(found_indications), self.param['input_compound']))

        effect = 'negative' if self.param['drug_effect'] == INHIBIT else 'positive'
        REQUEST_NAME = 'Find {drug} all other indications'.format(drug=self.param['input_compound'])
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = {eff} AND NeighborOf({select_drug}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_types}))'
        oql_query = OQLquery.format(eff=effect, select_drug=self.SELECTdrug,indication_types = indic_str)
        LiteratureIndications = self.process_oql(oql_query, REQUEST_NAME)
        found_indications= LiteratureIndications.get_entity_ids(self.indication_types)
        indications2return.update(found_indications)
        print('Found %d indications reported in scientific literature for %s' % (len(found_indications), self.param['input_compound']))
        self.drug_ids = self.Graph.get_entity_ids(['SmallMol'],['ObjTypeName'])
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
            found_indications = SimilarDrugsClinicalTrials.get_entity_ids(self.param['indication_types'])
            indications2return.update(found_indications)
            print('Found %d indications in clinical trials for drugs similar to %s' % (len(found_indications), self.param['input_compound']))

        effect = 'negative' if self.param['drug_effect'] == INHIBIT else 'positive'
        REQUEST_NAME = 'Find all other indications for {similars}'.format(similars=','.join(self.param['similars']))
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = {eff} AND NeighborOf ({select_drugs}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_types}))'
        OQLquery = OQLquery.format(eff=effect,select_drugs=select_similar_drugs,indication_types=indic_str)
        SimilarDrugsIndications = self.process_oql(OQLquery, REQUEST_NAME)
        self.similar_drug_ids = SimilarDrugsIndications.get_entity_ids(['SmallMol'])
        found_indications = SimilarDrugsIndications.get_entity_ids(self.param['indication_types'])
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
        self.report_pandas.append(ct_pd)

        effect_str = self.__get_effect_str()
        #fname_out = self.data_dir+self.fname_prefix()+" references.tsv"
        print('Printing research articles')
        research_ref_pd = self.refs2pd(self.drugcolname,['Regulation'],[effect_str],
                            pubid_types=['PMID','DOI'],biblio_props=['Title','PubYear'],print_snippets=True)

        research_ref_pd.name ='articles'
        self.report_pandas.append(research_ref_pd)
        #research_ref_pd.to_excel(self.excel_writer, sheet_name='Snippets', index=False)

    def set_strictmode_indications(self):
        self.indications4strictmode = self.drug_indications4strictmode.intersection(self.target_indications4strictmode)

    def _worksheet_prefix(self):
        indications = ','.join(self.indication_types)
        targets = ','.join(self.param['target_names'])
        ef = 'Act.' if self.param['drug_effect'] == ACTIVATE else 'Inh.'
        return ef+indications+'-'+targets

    def make_report(self, map2ontology:list):
        global_start = time.time()

        self.set_drug() 
        # if partner_names is empty script will try finding Ligands for Receptor targets and Receptors for Ligand targets in the knowledge graph
        self.set_targets()
        # strict mode ranks only indications suggested in the literarure without predicting new indications
   
        self.flush_dump_files()
        self.find_drug_indications()
        self.find_indications4similars()
        
        self.find_target_indications()
        self.indications4chem_modulators()
        self.indications4partners()
        self.get_pathway_componets()

        if self.param['strict_mode']: 
            self.set_strictmode_indications()

        if self.init_semantic_search():
            self.score_drug_semantics()
            self.score_semantics()
            ontology_map = self.child2parent(map2ontology)
            self.normalize_counts(ontology_map,bibliography=True)
            print("Repurposing %s was done in %s" % (self.param['input_compound'],self.execution_time(global_start)))
            return True
        else:
            return False


    def other_effects(self):
        old_rel_props = self.relProps
        self.add_rel_props(SENTENCE_PROPS+['Title','PubYear'])
        drug = self.param['input_compound']
        REQUEST_NAME = 'Find indications modulated by {} with unknown effect'.format(drug)
        select_indications = 'SELECT Entity WHERE objectType = ({})'.format(','.join(self.param['indication_types']))
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = unknown AND NeighborOf({select_target}) AND NeighborOf ({indications})'
        to_return = self.process_oql(OQLquery.format(select_target=self.SELECTdrug, indications=select_indications),REQUEST_NAME)
        self.relProps = old_rel_props
        return to_return
                                           

if __name__ == "__main__":
    APIconfig = load_api_config()
    parameters = {#  specify here what indications to find and the type of drug: 
                'input_compound' : 'cannabidiol',#'2-arachidonoylglycerol', #'anandamide', # , '9-tetrahydrocannabinol'
                'similars' : ['cannabidivarin', 'Cannabidiolic acid', 'Cannabielsoin'],#['noladin ether'],
                #['oleoylethanolamide','virodhamine','N-oleoyldopamine','palmitoylethanolamide','N-arachidonoyl dopamine','N-arachidonoylglycine','eicosapentaenoyl ethanolamide'],
                #'input_compound' : 'tetrahydrocannabinol', 
                #'similars' : ['delta 8-THC', '9-tetrahydrocannabinol', 'THC-C4','tetrahydrocannabinolic acid', '11-hydroxy-delta 9-tetrahydrocannabinol'],
                # names of similar drugs that have same mechanism of action (i.e. same target):
                'drug_effect': INHIBIT,
                'mode_of_action': ANTAGONIST,
                'partner_class':'Metabolite ligand', # use it only if partner_names not empty
                'indication_types': ['CellProcess'], #['Disease','Virus']
                'target_names':['GPR18'],
                'target_type':'Protein',
                'to_inhibit':True,
                'pathway_name_must_include_target':True,
                'strict_mode':RANK_SUGGESTED_INDICATIONS, #PREDICT_RANK_INDICATIONS#
                'data_dir':'D:/Python/PMI/',
                'partner_names':['endocannabinoid','anandamide','2-arachidonoylglycerol']
                #'oleoylethanolamide','virodhamine','N-oleoyldopamine','palmitoylethanolamide','N-arachidonoyl dopamine','N-arachidonoylglycine','eicosapentaenoyl ethanolamide'
                # 'noladin ether','2-arachidonoylglycerol'
                # specify here endogenous ligands for receptor if known
                }

    targets = ['CNR1','CNR2','GPR55','GPR119','GPR18']
    indications = ['CellProcess','Disease']
    path = parameters['data_dir']
    target_report = pd.ExcelWriter(path+parameters['input_compound']+'.xlsx', engine='xlsxwriter')
    raw_data_cache = pd.ExcelWriter(path+parameters['input_compound']+'raw_data.xlsx', engine='xlsxwriter')

    for target in targets:
        parameters['target_names'] = [target]
        parameters['indication_types'] = ['Disease']
        parameters['drug_effect'] = INHIBIT
        dcp = RepurposeDrug(APIconfig,parameters)
        if dcp.make_report(DISEASE_MAP2ONTOLOGY):
            dcp.add2writer(target_report,dcp._worksheet_prefix())
            dcp.addraw2writer(raw_data_cache,dcp._worksheet_prefix())

        parameters['target_names'] = [target]
        parameters['indication_types'] = ['CellProcess']
        parameters['drug_effect'] = INHIBIT
        dcp = RepurposeDrug(APIconfig,parameters)
        if dcp.make_report(CELLPROCESS_MAP2ONTOLOGY):
            dcp.add2writer(target_report,dcp._worksheet_prefix())
            dcp.addraw2writer(raw_data_cache,dcp._worksheet_prefix())

        parameters['target_names'] = [target]
        parameters['indication_types'] = ['CellProcess']
        parameters['drug_effect'] = ACTIVATE
        dcp = RepurposeDrug(APIconfig,parameters)
        if dcp.make_report(CELLPROCESS_MAP2ONTOLOGY):
            dcp.add2writer(target_report,dcp._worksheet_prefix())
            dcp.addraw2writer(raw_data_cache,dcp._worksheet_prefix())


    dcp = RepurposeDrug(APIconfig,parameters)
    dcp.set_drug()
    parameters['indication_types'] = ['Disease']
    other_indications = dcp.rn2pd(dcp.other_effects(),dcp.param['input_compound'])
    other_indications.to_excel(target_report, sheet_name='Possbl.Diseases', index=False)

    parameters['indication_types'] = ['CellProcess']
    other_indications = dcp.rn2pd(dcp.other_effects(),dcp.param['input_compound'])
    other_indications.to_excel(target_report, sheet_name='Possbl.CellProcess', index=False)

    target_report.save()
    raw_data_cache.save()

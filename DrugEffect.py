from ElsevierAPI import load_api_config
from TargetIndications import TargetIndications
import time

INHIBIT = 0 # indication must be inhibited by a drug
ACTIVATE = 1 # indication must be activated by a drug

ANTAGONIST = 0 # drug inhibits its targets
AGONIST = 1 # drug activates its targets

RANK_SUGGESTED_INDICATIONS = True 
# RANK_SUGGESTED_INDICATIONS indications suggested in the lietarure are ranked by amount of supporting evidence
PREDICT_RANK_INDICATIONS = False
# also predict indication using biology of the drug target(s)

class RepurposeDrug(TargetIndications):
    pass
    def __init__(self, APIconfig):
        super().__init__(APIconfig)
        self.PageSize = 1000
        self.similar_drugs = list()
        # default settings when drug is antagonist of receptor or ligand
        self.required_effect_on_indications = INHIBIT
        self.known_effect_on_targets = ANTAGONIST


    def set_drug(self, drug_name:str, similar_drugs:list, mode=INHIBIT, indication_types=None, type=ANTAGONIST):
        self.drug_name = drug_name
        self.similar_drugs = similar_drugs
        if isinstance(indication_types,list):
            self.indication_types = indication_types
            # defaults to ['Disease','Virus']
        
        self.SELECTdrug = 'SELECT Entity WHERE InOntology (SELECT Annotation WHERE Ontology=\'Pathway Studio Ontology\' AND Relationship=\'is-a\') under (SELECT OntologicalNode WHERE (Name,Alias) = (\'{drug}\')) OR (Name,Alias) = (\'{drug}\')'
        
        # else: self.SELECTdrug = 'SELECT Entity WHERE (Name,Alias) = {drug}'
            
        self.SELECTdrug = self.SELECTdrug.format(drug=drug_name)
        self.required_effect_on_indications = mode
        self.known_effect_on_targets = type
        if self.required_effect_on_indications == INHIBIT and self.known_effect_on_targets == ANTAGONIST:
            self.target_activate_indication = True #most often case for drugs and their disease indications
        elif self.required_effect_on_indications == ACTIVATE and self.known_effect_on_targets == ANTAGONIST:
            self.target_activate_indication = False # can be used only with CellProcess
        elif self.required_effect_on_indications == ACTIVATE and self.known_effect_on_targets == AGONIST:
            self.target_activate_indication = True # can be used only with CellProcess
        elif self.required_effect_on_indications == INHIBIT and self.known_effect_on_targets == AGONIST:
            self.target_activate_indication = False #find disease indications for drugs that are agonist of their target(s)

        # self.target_activate_indication = self.required_effect_on_indications * self.known_effect_on_targets

    def __get_effect_str(self):
        if self.required_effect_on_indications == INHIBIT: return 'negative'
        if self.required_effect_on_indications == ACTIVATE: return 'positive'
        return 'unknown'

    def needs_clinical_trial(self):
        if 'CellProcess' in self.indication_types and self.required_effect_on_indications == ACTIVATE: return True
        if 'Disease' in self.indication_types and self.required_effect_on_indications == INHIBIT: return True
        return False

    def find_drug_indications(self):
        indic_str = ','.join(self.indication_types)
        # PART I: Finding all possible indications
        if self.needs_clinical_trial():
            REQUEST_NAME = 'Find {drug} indications by clinical trials'.format(drug=self.drug_name)
            OQLquery = 'SELECT Relation WHERE objectType = ClinicalTrial AND NeighborOf ({select_drug}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_types}))'
            ClinicalTrialIndictions = self.process_oql(OQLquery.format(select_drug=self.SELECTdrug, indication_types = indic_str),REQUEST_NAME)
            found_indications = ClinicalTrialIndictions.get_entity_ids(self.indication_types)
            print('Found %d indications in %s clinical trials' % (len(found_indications), self.drug_name))

        effect = 'negative' if self.required_effect_on_indications == INHIBIT else 'positive'
        REQUEST_NAME = 'Find {drug} all other indications'.format(drug=self.drug_name)
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = {eff} AND NeighborOf({select_drug}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_types}))'
        oql_query = OQLquery.format(eff=effect, select_drug=self.SELECTdrug,indication_types = indic_str)
        LiteratureIndications = self.process_oql(oql_query, REQUEST_NAME)
        found_indications= LiteratureIndications.get_entity_ids(self.indication_types)
        print('Found %d indications reported in scientific literature for %s' % (len(found_indications), self.drug_name))
        self.drug_ids = self.Graph.get_entity_ids(['SmallMol'],['ObjTypeName'])

    def find_indications4similars(self):
        if not self.similar_drugs: return
        indic_str = ','.join(self.indication_types)
        if self.needs_clinical_trial():
            REQUEST_NAME = 'Find indications by clinical trials for {similars}'.format(similars=self.similar_drugs)
            select_similar_drugs = 'SELECT Entity WHERE (Name,Alias) = ({drugs})'
            select_similar_drugs = select_similar_drugs.format(drugs=self.similar_drugs)
            OQLquery = 'SELECT Relation WHERE objectType = ClinicalTrial AND NeighborOf ({select_drugs}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_types}))'
            SimilarDrugsClinicalTrials = self.process_oql(OQLquery.format(select_drugs=select_similar_drugs, indication_types=indic_str), REQUEST_NAME)
            found_indications = SimilarDrugsClinicalTrials.get_entity_ids(self.indication_types)
            print('Found %d indications in clinical trials for drugs similar to %s' % (len(found_indications), self.drug_name))

        effect = 'negative' if self.required_effect_on_indications == INHIBIT else 'positive'
        REQUEST_NAME = 'Find all other indications for {similars}'.format(similars=self.similar_drugs)
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = {eff} AND NeighborOf ({select_drugs}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_types}))'
        SimilarDrugsIndications = self.process_oql(OQLquery.format(eff=effect,select_drugs=select_similar_drugs,indication_types=indic_str), REQUEST_NAME)
        self.similar_drug_ids = SimilarDrugsIndications.get_entity_ids(['SmallMol'])
        found_indications = SimilarDrugsIndications.get_entity_ids(self.indication_types)
        print('Found %d indications reported in scientific literature or drugs similar to %s' % (len(found_indications), self.drug_name))

    def score_drug_semantics(self):
        colname = self.drug_name + ' cliniclal trials'
        self.set_how2connect(['ClinicalTrial'],[],'')
        linked_entities_count = self.link2concept(colname, self.drug_ids)
        print('%s has clinical trials for %d indications' % (self.drug_name , linked_entities_count))

        if self.required_effect_on_indications == INHIBIT:
            colname = 'Inhibited by ' + self.drug_name
            self.set_how2connect(['Regulation'],['negative'],'')
            regulated = 'inhibited'
        else:
            colname = 'Activated by ' + self.drug_name
            self.set_how2connect(['Regulation'],['positive'],'')
            regulated = 'activated'

        linked_entities_count = self.link2concept(colname, self.drug_ids)
        print('%d indications %s by %s' % (linked_entities_count, regulated, self.drug_name))
        self.drugcolname = self._col_name_prefix+colname

        if hasattr(self,'similar_drug_ids'):
            colname = 'similar drugs cliniclal trials'
            self.set_how2connect(['ClinicalTrial'],[],'')
            linked_entities_count = self.link2concept(colname, self.similar_drug_ids)
            print('%s has clinical trials for %s indications' % (self.similar_drugs, linked_entities_count))

            if self.required_effect_on_indications == INHIBIT:
                colname = 'Inhibited by similar drugs'
                self.set_how2connect(['Regulation'],['negative'],'')
            else:
                colname = 'Activated by similar drugs'
                self.set_how2connect(['Regulation'],['positive'],'')
    
            linked_entities_count = self.link2concept(colname, self.similar_drug_ids)
            print('%d indications %s by %s' % (linked_entities_count, regulated, self.similar_drugs))
        
    def fname_prefix(self):
        indics = ','.join(self.indication_types)
        rep_pred =  'suggested ' if self.strict_mode else 'suggested,predicted '
        regulate = ' activated by ' if self.required_effect_on_indications == ACTIVATE else ' inhibited by '
        return rep_pred+indics+regulate+self.drug_name

    def print_drug_indictaion_refs(self):
        fname_prefix = dcp.fname_prefix()
        dcp.print_references(self.drugcolname,fname_prefix + " clinical trials.tsv", ['ClinicalTrial'])

        effect_str = self.__get_effect_str()
        dcp.print_references(self.drugcolname,fname_prefix + " references.tsv", ['Regulation'],[effect_str])

if __name__ == "__main__":
    global_start = time.time()
    dcp = RepurposeDrug(load_api_config())

    # specify here what indications to find and the type of drug
    similars = [] # names of similar drugs that have same mechanism of action (i.e. same target)
    #dcp.set_drug('cannabinoids', similars, INHIBIT, ['Disease'], AGONIST)
    dcp.set_drug('CBD compounds', similars, INHIBIT, ['Disease'], ANTAGONIST)
    
    # specify here the drug targets and drug mechanism of action
    partner_names = ['anandamide','endocannabinoid','2-arachidonoylglycerol'] # specify here endogenous ligands for receptor if known
    # if partner_names is empty script will try finding Ligands for Receptor targets and Receptors for Ligand targets in the knowledge graph
    partner_class = 'Metabolite ligand' 
    dcp.set_targets(['CNR1'],'Protein',partner_names=partner_names,partner_class=partner_class,
            strict_mode=RANK_SUGGESTED_INDICATIONS)
    # strict mode ranks only indications suggested in the literarure without predicting new indications
    dcp.pathway_name_must_include_target = True
    # must be True if there are a lot of curated pathways with receptor name included into pathway name
    
    dcp.flush_dump_files()
    dcp.find_drug_indications()
    if similars:
        dcp.find_indications4similars()

    dcp.strict_mode = True
    if not dcp.strict_mode:
        dcp.find_target_indications()
        dcp.indications4chem_modulators()
        dcp.indications4partners()
    else:
        dcp.get_GVs()
    
    dcp.get_pathway_componets()

    dcp.init_semantic_search()
    dcp.score_drug_semantics()
    dcp.score_semantics()

    fname_prefix = dcp.fname_prefix()
    dcp.print_ref_count(fname_prefix + " counts.tsv",sep='\t')
    dcp.print_drug_indictaion_refs()

    NormalizedCount = dcp.normalize_counts()
    fout = fname_prefix + ' normalized report.tsv'
    NormalizedCount.to_csv(fout, sep='\t', index=False, float_format='%g')
    print("Repurposing %s was done in %s" % (dcp.drug_name , dcp.execution_time(global_start)))

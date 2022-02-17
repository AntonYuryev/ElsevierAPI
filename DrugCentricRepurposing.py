from ElsevierAPI import load_api_config
from TargetCentricRepurposing import TargetIndications
import time


class RepurposeDrug(TargetIndications):
    pass
    def __init__(self, APIconfig):
        super().__init__(APIconfig)


    def set_drug(self, drug_name:str, similar_drugs:list, include_children=True, indication_types=None):
        self.drug_name = drug_name
        self.similar_drugs = similar_drugs
        if isinstance(indication_types,list):
            self.indication_types = indication_types
            # defaults to ['Disease','Virus']
        if include_children:
            self.SELECTdrug = 'SELECT Entity WHERE InOntology (SELECT Annotation WHERE Ontology=\'Pathway Studio Ontology\' AND Relationship=\'is-a\') under (SELECT OntologicalNode WHERE (Name,Alias) = (\'{drug}\')) OR (Name,Alias) = (\'{drug}\')'
        else:
            self.SELECTdrug = 'SELECT Entity WHERE (Name,Alias) = {drug}'
            
        self.SELECTdrug = self.SELECTdrug.format(drug=drug_name)

    def find_drug_indications(self):
        # PART I: Finding all possible indications
        REQUEST_NAME = 'Find {drug} indications by clinical trials'.format(drug=self.drug_name)
        OQLquery = 'SELECT Relation WHERE objectType = ClinicalTrial AND NeighborOf ({select_drug}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_types}))'
        ClinicalTrialIndictions = self.process_oql(OQLquery.format(select_drug=self.SELECTdrug, indication_types = self.indication_types),REQUEST_NAME)
        found_indications = ClinicalTrialIndictions.get_entity_ids(self.indication_types)
        print('Found %d indications in %s clinical trials' % (len(found_indications), self.drug_name))

        effect = 'negative' if self.target_activate_indication else 'positive'
        REQUEST_NAME = 'Find {drug} all other indications'.format(drug=self.drug_name)
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = {eff} AND NeighborOf({select_drug}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_types})'
        LiteratureIndications = self.process_oql(OQLquery.format(eff=effect, select_drug=self.SELECTdrug,indication_types = self.indication_types), REQUEST_NAME)
        found_indications= LiteratureIndications.get_entity_ids(self.indication_types)
        print('Found %d indications reported in scientific literature for %s' % (len(found_indications), self.drug_name))
        self.drug_ids = self.Graph.get_entity_ids(['SmallMol'],['ObjTypeName'])

    def find_indications4similars(self):
        REQUEST_NAME = 'Find indications by clinical trials for {similars}'.format(similars=self.similar_drugs)
        select_similar_drugs = 'SELECT Entity WHERE (Name,Alias) = ({drugs})'
        select_similar_drugs = select_similar_drugs.format(drugs=self.similar_drugs)
        OQLquery = 'SELECT Relation WHERE objectType = ClinicalTrial AND NeighborOf ({select_drugs}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_types}))'
        SimilarDrugsClinicalTrials = self.process_oql(OQLquery.format(select_drugs=select_similar_drugs, indication_types=self.indication_types), REQUEST_NAME)
        found_indications = SimilarDrugsClinicalTrials.get_entity_ids(self.indication_types)
        print('Found %d indications in clinical trials for drugs similar to %s' % (len(found_indications), self.drug_name))

        effect = 'negative' if self.target_activate_indication else 'positive'
        REQUEST_NAME = 'Find all other indications for {similars}'.format(similars=self.similar_drugs)
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = {eff} AND NeighborOf ({select_drugs}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_types}))'
        SimilarDrugsIndictions = self.process_oql(OQLquery.format(eff=effect,select_drugs=select_similar_drugs,indication_types=self.indication_types), REQUEST_NAME)
        self.similar_drug_ids = self.Graph.get_entity_ids(self.similar_drugs, search_by_properties=['Name','Alias'])
        found_indications = SimilarDrugsIndictions.get_entity_ids(self.indication_types)
        print('Found %d indications reported in scientific literature or drugs similar to %s' % (len(found_indications), self.drug_name))

    def score_drug_semantics(self):
        colname = self.drug_name + ' cliniclal trials'
        self.set_how2connect(['ClinicalTrial'],[],'')
        linked_entities_count = self.link2concept(colname, self.drug_ids)
        print('%s has clinical trials for %d indications' % (self.drug_name , linked_entities_count))

        if self.target_activate_indication:
            colname = 'Inhibited by ' + self.drug_name
            self.set_how2connect(['Regulation'],['negative'],'')
            regulated = 'inhibited'
        else:
            colname = 'Activated by ' + self.drug_name
            self.set_how2connect(['Regulation'],['positive'],'')
            regulated = 'activated'
        
        linked_entities_count = self.link2concept(colname, self.drug_ids)
        print('%d indications %s by %s' % (linked_entities_count, regulated, self.drug_name))

        colname = 'similar drugs cliniclal trials'
        self.set_how2connect(['ClinicalTrial'],[],'')
        linked_entities_count = self.link2concept(colname, self.similar_drug_ids)
        print('%s has clinical trials for %s indications' % (self.similar_drugs, linked_entities_count))

        if self.target_activate_indication:
            colname = 'Inhibited by similar drugs'
            self.set_how2connect(['Regulation'],['negative'],'')
        else:
            colname = 'Activated by similar drugs'
            self.set_how2connect(['Regulation'],['positive'],'')
  
        linked_entities_count = self.link2concept(colname, self.similar_drug_ids)
        print('%d indications %s by %s' % (linked_entities_count, regulated, dcp.similar_drugs))
        
    def fname_prefix(self):
        indics = ','.join(self.indication_types)
        rep_pred =  'suggested,predicted ' if self.target_activate_indication else 'suggested '
        return rep_pred+indics+' for '+ self.drug_name

if __name__ == "__main__":
    global_start = time.time()
    dcp = RepurposeDrug(load_api_config())

    drugname = 'THC compounds'
    targets = ['CNR1','CNR2','GPR18','GPR55','GPR119']
    similars = [] #['Infliximab', 'Etanercept', 'Golimumab', 'Certolizumab', 'Afelimomab'] # for 'Humira'/TNF
    #similars = ['osimertinib', 'gefitinib', 'brigatinib', 'zalutumumab', 'nimotuzumab', 'matuzumab'] # for 'erlotinib'/'EGFR
    partner_names = ['tetrahydrocannabinol','anandamide']
    partner_class = 'Metabolite ligand'

    dcp.set_drug(drugname, similars, include_children=True, indication_types=['CellProcess'])
    dcp.set_targets(targets,'Protein',partner_names=partner_names,partner_class=partner_class,to_inhibit=False,strict_mode=True)
    # strict mode ranks only indications suggested in the lietarure without predicting new indications
    dcp.pathway_name_must_include_target = True
    #dcp.add_ent_props(['Alias'])
    dcp.PageSize = 1000

    dcp.flush_dump_files()
    dcp.find_drug_indications()
    dcp.find_indications4similars()
    dcp.find_target_indications()
    dcp.indications4chem_modulators()
    dcp.indications4partners()
    dcp.get_pathway_componets()

    dcp.init_semantic_search()
    dcp.score_drug_semantics()
    dcp.score_semantics()

    fname_prefix = dcp.fname_prefix()
    dcp.print_ref_count(fname_prefix + " counts.tsv",fname_prefix + " references.tsv",sep='\t')
    NormalizedCount = dcp.normalize_counts()
    fout = fname_prefix + ' normalized report.tsv'
    NormalizedCount.to_csv(fout, sep='\t', index=False, float_format='%g')
    print("Repurposing %s was done in %s" % (dcp.drug_name , dcp.execution_time(global_start)))

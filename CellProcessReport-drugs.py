from ElsevierAPI import load_api_config
from TargetCentricRepurposing import RepurposeDrugs
import time


class RepurposeDrug(RepurposeDrugs):
    pass
    def __init__(self, APIconfig):
        super().__init__(APIconfig)
        #self.PageSize = 500

    def set_drug(self, drug_name:str, similar_drugs:list, include_children=True):
        self.Drug = drug_name
        self.similar_drugs = similar_drugs
        if include_children:
            self.SELECTdrug = 'SELECT Entity WHERE InOntology (SELECT Annotation WHERE Ontology=\'Pathway Studio Ontology\' AND Relationship=\'is-a\') under (SELECT OntologicalNode WHERE (Name,Alias) = (\'{drug}\')) OR (Name,Alias) = (\'{drug}\')'
        else:
            self.SELECTdrug = 'SELECT Entity WHERE (Name,Alias) = {drug}'
            
        self.SELECTdrug = self.SELECTdrug.format(drug=dcp.Drug)
        


if __name__ == "__main__":
    global_start = time.time()
    #similars = ['osimertinib', 'gefitinib', 'brigatinib', 'zalutumumab', 'nimotuzumab', 'matuzumab'] # for 'erlotinib'/'EGFR
    similars = [] #['Infliximab', 'Etanercept', 'Golimumab', 'Certolizumab', 'Afelimomab'] # for 'Humira'/TNF
    dcp = RepurposeDrug(load_api_config())
    dcp.set_drug('CBD compounds', similars, include_children=True)
    dcp.set_targets(['CNR1', 'CNR2', 'GPR18', 'GPR55', 'GPR119'],'Protein', to_inhibit=True, strict_mode=True)
    dcp.add_ent_props(['Alias'])
    dcp.PageSize = 1000

    dcp.flush_dump_files()

    # PART I: Finding all possible indications
    REQUEST_NAME = 'Find {drug} indications by clinical trials'.format(drug=dcp.Drug)
    OQLquery = 'SELECT Relation WHERE objectType = ClinicalTrial AND NeighborOf ({select_drug}) AND NeighborOf (SELECT Entity WHERE objectType = (CellProcess))'
    ClinicalTrialIndictions = dcp.process_oql(OQLquery.format(select_drug=dcp.SELECTdrug),REQUEST_NAME)
    found_diseases = ClinicalTrialIndictions.get_entity_ids(['CellProcess'])
    print('Found %d indications in %s clinical trials' % (len(found_diseases), dcp.Drug))

    REQUEST_NAME = 'Find {drug} all other indications'.format(drug=dcp.Drug)
    OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = negative AND NeighborOf({select_drug}) AND NeighborOf (SELECT Entity WHERE objectType = (CellProcess))'
    LiteratureIndications = dcp.process_oql(OQLquery.format(select_drug=dcp.SELECTdrug), REQUEST_NAME)
    found_diseases= LiteratureIndications.get_entity_ids(['CellProcess'])
    print('Found %d indications reported in scientific literature for %s' % (len(found_diseases), dcp.Drug))
    drug_ids = dcp.Graph.get_entity_ids(['SmallMol'],['ObjTypeName'])

    REQUEST_NAME = 'Find indications by clinical trials for {similars}'.format(similars=dcp.similar_drugs)
    select_similar_drugs = 'SELECT Entity WHERE (Name,Alias) = ({drugs})'
    select_similar_drugs = select_similar_drugs.format(drugs=dcp.similar_drugs)
    OQLquery = 'SELECT Relation WHERE objectType = ClinicalTrial AND NeighborOf ({select_drugs}) AND NeighborOf (SELECT Entity WHERE objectType = (CellProcess))'
    SimilarDrugsClinicalTrials = dcp.process_oql(OQLquery.format(select_drugs=select_similar_drugs), REQUEST_NAME)
    found_diseases= SimilarDrugsClinicalTrials.get_entity_ids(['CellProcess'])
    print('Found %d indications in clinical trials for drugs similar to %s' % (len(found_diseases), dcp.Drug))

    REQUEST_NAME = 'Find all other indications for {similars}'.format(similars=dcp.similar_drugs)
    OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = negative AND NeighborOf ({select_drugs}) AND NeighborOf (SELECT Entity WHERE objectType = (CellProcess))'
    SimilarDrugsIndictions = dcp.process_oql(OQLquery.format(select_drugs=select_similar_drugs), REQUEST_NAME)
    similar_drug_ids = dcp.Graph.get_entity_ids(dcp.similar_drugs, search_by_properties=['Name', 'Alias'])
    found_diseases = SimilarDrugsIndictions.get_entity_ids(['CellProcess'])
    print('Found %d indications reported in scientific literature or drugs similar to %s' % (len(found_diseases), dcp.Drug))

    dcp.find_target_indications()
    dcp.indications4chem_modulators()
    dcp.indications4partners()
    dcp.get_pathway_componets()

    dcp.init_semantic_search()

    colname = dcp.Drug + ' cliniclal trials'
    dcp.set_how2connect(['ClinicalTrial'],[],'')
    linked_entities_count = dcp.link2concept(colname, drug_ids)
    print('%s has clinical trials for %d indications' % (dcp.Drug, linked_entities_count))

    colname = 'Inhibited by ' + dcp.Drug
    dcp.set_how2connect(['Regulation'],['negative'],'')
    linked_entities_count = dcp.link2concept(colname, drug_ids)
    print('%d diseases inhibited by %s' % (linked_entities_count, dcp.Drug))

    colname = 'similar drugs cliniclal trials'
    dcp.set_how2connect(['ClinicalTrial'],[],'')
    linked_entities_count = dcp.link2concept(colname, similar_drug_ids)
    print('%s has clinical trials for %s indications' % (dcp.similar_drugs, linked_entities_count))

    colname = 'Inhibited by similar drugs'
    dcp.set_how2connect(['Regulation'],['negative'],'')
    linked_entities_count = dcp.link2concept(colname, similar_drug_ids)
    print('%d diseases inhibited by %s' % (linked_entities_count, dcp.similar_drugs))

    dcp.score_semantics()
    dcp.print_ref_count(dcp.Drug+" repurposing counts.tsv",dcp.Drug+" repurposing references.tsv",sep='\t')
    NormalizedCount = dcp.normalize_counts()
    fout = dcp.Drug + ' repurposing normalized report.tsv'
    NormalizedCount.to_csv(fout, sep='\t', index=False, float_format='%g')
    print("Repurposing %s was done in %s" % (dcp.Drug, dcp.execution_time(global_start)))
from ElsevierAPI import APIconfig
from TargetCentricRepurposing import RepurposeDrugs
import time


class RepurposeDrug(RepurposeDrugs):
    pass
    def __init__(self, APIconfig):
        super().__init__(APIconfig)
        self.PageSize = 500

    def set_drug(self, drug_name: str, similar_drugs: list):
        self.Drug = drug_name
        self.similar_drugs = similar_drugs

    def pathway_oql(self):
        if self.partners != NotImplemented:
            SELECTpathways = 'SELECT Network WHERE ParentOf ({select_target}) AND ParentOf ({select_partner})'#Name LIKE (\''+percent_sign+'FGF'+percent_sign+'\')'
            SELECTpathways = SELECTpathways.format(select_target=self.find_target_oql,select_partner=self.find_partners_oql)
        else:
            SELECTpathways = 'SELECT Network WHERE ParentOf ({select_target})'#Name LIKE (\''+percent_sign+'FGF'+percent_sign+'\')'
            SELECTpathways = SELECTpathways.format(select_target=self.find_target_oql)
        
        pct = '%'
        drug_target_name = self.Drug_Target['Name'][0]
        SELECTpathways = SELECTpathways + ' AND Name LIKE (\''+pct+drug_target_name+pct+'\')' #additional refinement for popular targets
        SELECTMergedPathway = 'SELECT Relation WHERE objectType = (DirectRegulation,Binding,ProtModification,PromoterBinding,ChemicalReaction) AND MemberOf ({select_networks})'
        return SELECTMergedPathway.format(select_networks=SELECTpathways)

if __name__ == "__main__":
    global_start = time.time()
    #similars = ['osimertinib', 'gefitinib', 'brigatinib', 'zalutumumab', 'nimotuzumab', 'matuzumab'] # for 'erlotinib'/'EGFR
    similars = ['Infliximab', 'Etanercept', 'Golimumab', 'Certolizumab', 'Afelimomab'] # for 'Humira'/TNF
    dcp = RepurposeDrug(APIconfig)
    dcp.set_drug('Humira', similars)
    dcp.set_targets(['TNF'],'Protein')
    dcp.add_ent_props(['Alias'])

    SELECTdrug = 'SELECT Entity WHERE (Name,Alias) = {drug}'
    SELECTdrug = SELECTdrug.format(drug=dcp.Drug)
    dcp.flush_dump_files()

    # PART I: Finding all possible indications
    REQUEST_NAME = 'Find {drug} clinical trials'.format(drug=dcp.Drug)
    OQLquery = 'SELECT Relation WHERE objectType = ClinicalTrial AND NeighborOf ({select_drug}) AND NeighborOf (SELECT Entity WHERE objectType = (Disease))'
    ClinicalTrialIndictions = dcp.process_oql(OQLquery.format(select_drug=SELECTdrug),REQUEST_NAME)
    found_diseases= ClinicalTrialIndictions.get_entity_ids(['Disease'])
    print('Found %d indications in %s clinical trials' %  (len(found_diseases), dcp.Drug))

    drug_id = dcp.Graph.get_entity_ids([dcp.Drug], search_by_properties=['Name', 'Alias'])

    REQUEST_NAME = 'Find {drug} indications'.format(drug=dcp.Drug)
    OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = negative AND NeighborOf({select_drug}) AND NeighborOf (SELECT Entity WHERE objectType = (Disease,Virus))'
    LiteratureIndications = dcp.process_oql(OQLquery.format(select_drug=SELECTdrug), REQUEST_NAME)
    found_diseases= LiteratureIndications.get_entity_ids(['Disease'])
    print('Found %d indications reported in scientific literature for %s' %  (len(found_diseases), dcp.Drug))

     
    similarDrugs = ['Infliximab','Etanercept','Golimumab','Certolizumab','Afelimomab']
    similarDrugsStr = ','.join(similarDrugs)
    REQUEST_NAME = 'Find clinical trials for {similars}'.format(similars=similarDrugsStr)
    select_similar_drugs = 'SELECT Entity WHERE (Name,Alias) = ({drugs})'
    select_similar_drugs = select_similar_drugs.format(drugs=similarDrugsStr)
    OQLquery = 'SELECT Relation WHERE objectType = ClinicalTrial AND NeighborOf ({select_drugs}) AND NeighborOf (SELECT Entity WHERE objectType = (Disease))'
    SimilarDrugsClinicalTrials = dcp.process_oql(OQLquery.format(select_drugs=select_similar_drugs), REQUEST_NAME)
    found_diseases= SimilarDrugsClinicalTrials.get_entity_ids(['Disease'])
    print('Found %d indications in clinical trials for drugs similar to %s' %  (len(found_diseases), dcp.Drug))

    REQUEST_NAME = 'Find indications for {similars}'.format(similars=similarDrugsStr)
    OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = negative AND NeighborOf ({select_drugs}) AND NeighborOf (SELECT Entity WHERE objectType = (Disease))'
    SimilarDrugsIndictions = dcp.process_oql(OQLquery.format(select_drugs=select_similar_drugs), REQUEST_NAME)
    similar_drug_ids = dcp.Graph.get_entity_ids(similarDrugs, search_by_properties=['Name', 'Alias'])
    found_diseases = SimilarDrugsIndictions.get_entity_ids(['Disease'])
    print('Found %d indications reported in scientific literature or drugs similar to %s' %  (len(found_diseases), dcp.Drug))

    dcp.find_target_indications()
    dcp.indications4chem_modulators()
    dcp.indications4partners()
    dcp.get_pathway_componets()
    dcp.init_semantic_search()

    colname = dcp.Drug + ' cliniclal trials'
    dcp.set_how2connect(['ClinicalTrial'],[],'')
    linked_entities_count = dcp.link2concept(colname, drug_id)
    print('%s has clinical trials for %d indications' % (dcp.Drug, linked_entities_count))

    colname = 'Inhibited by ' + dcp.Drug
    dcp.set_how2connect(['Regulation'],['negative'],'')
    linked_entities_count = dcp.link2concept(colname, drug_id)
    print('%d diseases inhibited by %s' % (linked_entities_count, dcp.Drug))

    colname = 'similar drugs cliniclal trials'
    dcp.set_how2connect(['ClinicalTrial'],[],'')
    linked_entities_count = dcp.link2concept(colname, similar_drug_ids)
    print('%s has clinical trials for %s indications' % (similarDrugsStr, linked_entities_count))

    colname = 'Inhibited by similar drugs'
    dcp.set_how2connect(['Regulation'],['negative'],'')
    linked_entities_count = dcp.link2concept(colname, similar_drug_ids)
    print('%d diseases inhibited by %s' % (linked_entities_count, similarDrugsStr))

    dcp.score_semantics()
    dcp.print_ref_count(dcp.Drug+" repurposing counts.tsv",dcp.Drug+" repurposing references.tsv",sep='\t')
    NormalizedCount = dcp.normalize_counts()
    fout = dcp.Drug + ' repurposing normalized report.tsv'
    NormalizedCount.to_csv(fout, sep='\t', index=False, float_format='%g')
    print("Repurposing %s was done in %s" % (dcp.Drug, dcp.execution_time(global_start)))

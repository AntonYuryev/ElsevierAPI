from ElsevierAPI import load_api_config
from ElsevierAPI.ResnetAPI.ResnetGraph import ResnetGraph
from ElsevierAPI.ResnetAPI.SemanticSearch import SemanticSearch
import ElsevierAPI.ResnetAPI.PathwayStudioGOQL as OQL
import pandas as pd
import numpy as np
import time
import networkx as nx

class RepurposeDrugs(SemanticSearch):
    pass
    
    def __init__(self, APIconfig):
        super().__init__(APIconfig)
        self.PageSize = 1000
        self.add_rel_props(['Name','Sentence','PubYear','Title','RelationNumberOfReferences'])
        #self.strict_mode = False
        
    def set_targets(self, target_names:list, target_objtype:str, to_inhibit=True, strict_mode=False):
        self.strict_mode = strict_mode

        prop_names_str, prop_values_str = OQL.get_search_strings(['Name'],target_names)
        self.add_ent_props(['Class'])
        self.find_targets_oql = 'SELECT Entity WHERE ({prop_name}) = ({values}) AND objectType = '+target_objtype 
        self.process_oql(self.find_targets_oql.format(prop_name=prop_names_str,values=prop_values_str))
        
        self.target_ids = self.Graph.get_entity_ids(SearchValues=target_names,search_by_properties=['Name'])
        self.find_targets_oql = 'SELECT Entity WHERE id = ('+ ','.join(map(str,self.target_ids))+')'
        self.Drug_Targets = self.Graph._get_nodes(self.target_ids)

        self.repurpose_antagonist = to_inhibit
        self.partners = self.select_partners()
        if self.partners:
            partners_ids = [p['Id'][0] for p in self.partners]
            self.find_partners_oql = 'SELECT Entity WHERE id = ('+ ','.join(map(str,partners_ids))+')'
        else: self.find_partners_oql = ''


    def _partner_class(self):
        #finding effect for complementary receptor or ligands
        self.target_class = self.Drug_Targets[0]['Class'][0]
        if self.target_class == 'Ligand': self.partner_class = "Receptor"
        elif self.target_class == 'Receptor': self.partner_class = "Ligand"
        else: self.target_class = NotImplemented


    def select_partners(self):
        self._partner_class()
        SELECTpartners = 'SELECT Entity WHERE Class = {partner_class} AND objectType = Protein AND Connected by (SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive) to ({select_target})'
        partners_graph = self.process_oql(SELECTpartners.format(partner_class=self.partner_class, select_target=self.find_targets_oql))
        if self.partner_class == 'Ligand':
            SELECTsecretedpartners = 'SELECT Entity WHERE "Cell Localization" = Secreted AND objectType = Protein AND Connected by (SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive) to ({select_target})'
            secreted_partners = self.process_oql(SELECTsecretedpartners.format(select_target=self.find_targets_oql))
            partners_graph = nx.compose(partners_graph,secreted_partners)

        partners_ids = partners_graph.get_entity_ids(['Protein'])
        target_names = ','.join([x['Name'][0] for x in self.Drug_Targets])
        partners_str = ','.join(n['Name'][0] for n in partners_graph._get_nodes())
        print ('Found %s %ss as partners for %s'%(partners_str, self.partner_class,target_names))
        return self.Graph._get_nodes(partners_ids)


    def find_target_indications(self):
        f_t = self.find_targets_oql
        target_names = [x['Name'][0] for x in self.Drug_Targets]
        t_n = ','.join(target_names)
        effect = 'positive' if self.repurpose_antagonist else 'negative'
        REQUEST_NAME = 'Find indications {effect}ly modulated by {target}'.format(effect=effect, target=t_n)
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = {effect} AND NeighborOf({select_target}) AND NeighborOf (SELECT Entity WHERE objectType = (CellProcess))'      
        ModulatedByTargetNetwork = self.process_oql(OQLquery.format(select_target=f_t, effect = effect),REQUEST_NAME)
        found_disease_ids= ModulatedByTargetNetwork.get_entity_ids(['CellProcess'])
        print('Found %d diseases modulated by target %s' % (len(found_disease_ids), t_n))

        if self.strict_mode:
            select_diseases_oql = 'SELECT Entity WHERE id = ({ids})'.format(ids=','.join(map(str,found_disease_ids)))
        else:
            select_diseases_oql = 'SELECT Entity WHERE objectType = (CellProcess)'


        REQUEST_NAME = 'Find indications {effect}ly regulating {target}'.format(effect=effect, target=t_n)
        OQLquery = 'SELECT Relation WHERE objectType = QuantitativeChange AND Effect = {effect} AND NeighborOf ({select_target}) AND NeighborOf ({select_diseases})'  
        ActivatedInDiseaseNetwork = self.process_oql(OQLquery.format(select_target=f_t, effect = effect,select_diseases=select_diseases_oql),REQUEST_NAME)
        found_diseases= ActivatedInDiseaseNetwork.get_entity_ids(['CellProcess'])
        print('Found %d diseases where target %s is %sly regulated' % (len(found_diseases), t_n, effect))
        
        selectGVs = 'SELECT Entity WHERE objectType = GeneticVariant AND Connected by (SELECT Relation WHERE objectType = GeneticChange) to ({select_target})'
        selectGVs = selectGVs.format(select_target=f_t)
        REQUEST_NAME = 'Find indications linked to {target} genetic variants'.format(target=t_n)
        OQLquery = 'SELECT Relation WHERE NeighborOf({select_gvs}) AND NeighborOf ({select_diseases})'
        self.GVsInDiseaseNetwork = self.process_oql(OQLquery.format(select_gvs=selectGVs,select_diseases=select_diseases_oql),REQUEST_NAME)
        found_diseases= self.GVsInDiseaseNetwork.get_entity_ids(['CellProcess'])
        found_GVs = self.GVsInDiseaseNetwork.get_entity_ids(['GeneticVariant'])
        print('Found %d diseases genetically linked to %d Genetic Variants in %s' % 
             (len(found_diseases), len(found_GVs), t_n))

        REQUEST_NAME = 'Find indications where {target} is biomarker'.format(target=t_n)
        OQLquery = 'SELECT Relation WHERE objectType = Biomarker AND NeighborOf({select_target}) AND NeighborOf ({select_diseases})'
        BiomarkerInDiseaseNetwork = self.process_oql(OQLquery.format(select_target=f_t,select_diseases=select_diseases_oql),REQUEST_NAME)
        found_diseases = BiomarkerInDiseaseNetwork.get_entity_ids(['CellProcess'])
        print('Found %d diseases where target %s is claimed as biomarker' %  (len(found_diseases), t_n))


    def indications4chem_modulators(self, linked_to_target_by=['DirectRegulation']):
        f_t = self.find_targets_oql
        t_n = ','.join([x['Name'][0] for x in self.Drug_Targets])
        effect = 'negative'if self.repurpose_antagonist else 'positive'
        
        REQUEST_NAME = 'Find substances {effect}ly regulating {target}'.format(effect=effect,target=t_n)
        OQLquery = 'SELECT Relation WHERE objectType = {rel_type} AND Effect = {effect} AND NeighborOf ({select_target}) AND NeighborOf (SELECT Entity WHERE objectType = SmallMol AND Connectivity > 1)'   
        TargetInhibitorsNetwork = self.process_oql(OQLquery.format(rel_type=linked_to_target_by, select_target=f_t, effect=effect),REQUEST_NAME)
        self.TargetInhibitorsIDs = TargetInhibitorsNetwork.get_entity_ids(['SmallMol'])
        print('Found %d substances %sly regulating %s' % (len(self.TargetInhibitorsIDs),effect,t_n))
        
        REQUEST_NAME = 'Find indications for substances {effect}ly regulating {target}'.format(effect=effect,target=t_n)
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = {effect}  AND NeighborOf (SELECT Entity WHERE objectType = Disease) AND NeighborOf (SELECT Entity WHERE id = ({drugIDlist}))'
        OQLquery = OQLquery.format(drugIDlist=','.join(map(str, self.TargetInhibitorsIDs)), effect=effect)
        InhibitorsIndicationNetwork = self.process_oql(OQLquery,REQUEST_NAME)

        found_diseases= InhibitorsIndicationNetwork.get_entity_ids(['CellProcess'])
        print('Found %d indications for %d substances %sly regulating %s' %  
             (len(found_diseases), len(self.TargetInhibitorsIDs),effect,t_n))


    def counterindications4chem_antimodulators(self,liked_to_target_by=['DirectRegulation']):
        # liked_to_target_by can be relaxed to (DirectRegulation,Regulation,Expression,MolTransport)
        f_t = self.find_targets_oql
        t_n = ','.join([x['Name'][0] for x in self.Drug_Targets])
        opposite_effect = 'positive' if self.repurpose_antagonist else 'negative'
        REQUEST_NAME = 'Find substances {effect}ly regulating {target}'.format(effect=opposite_effect,target=t_n)
        OQLquery = 'SELECT Relation WHERE objectType = {rel_type} and Effect = {effect} AND NeighborOf ({select_target}) AND NeighborOf (SELECT Entity WHERE objectType = SmallMol AND Connectivity > 1)'     
        TargetAgonistNetwork = self.process_oql(OQLquery.format(rel_type=liked_to_target_by,select_target=f_t, effect=opposite_effect),REQUEST_NAME)
        self.TargetAgonistIDs = list(set([x for x,y in TargetAgonistNetwork.nodes(data=True) if y['ObjTypeName'][0] in ['SmallMol']]))
        print('Found %d substances %sly regulating %s' % (len(self.TargetAgonistIDs),opposite_effect,t_n))
        
        REQUEST_NAME = 'Find counterindications for substances {effect}ly regulating {target}'.format(effect=opposite_effect,target=t_n)
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = {effect} AND NeighborOf (SELECT Entity WHERE objectType = Disease) AND NeighborOf (SELECT Entity WHERE id = ({drugIDlist}))'
        OQLquery = OQLquery.format(drugIDlist=','.join(map(str, self.TargetAgonistIDs)), effect=opposite_effect)
        AgonistToxicitiesnNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        found_diseases= AgonistToxicitiesnNetwork.get_entity_ids(['CellProcess'])
        print('Found %d indications for %d substances %sly regulating %s' %  
             (len(found_diseases), len(self.TargetAgonistIDs),opposite_effect,t_n))

     
    def indications4partners(self):
        t_n = ','.join([x['Name'][0] for x in self.Drug_Targets])
        effect = 'positive' if self.repurpose_antagonist else 'negative'
        REQUEST_NAME = 'Find indications for {partner}s of {targets} '.format(targets=t_n,partner=self.partner_class.lower())
        OQLquery = 'SELECT Relation WHERE objectType = (Regulation,QuantitativeChange) AND Effect = {effect} AND NeighborOf ({select_partners}) AND NeighborOf (SELECT Entity WHERE objectType = Disease)'
        if self.partners:
            OQLquery = OQLquery.format(effect=effect,select_partners=self.find_partners_oql)
            self.PartnerIndicationNetwork = self.process_oql(OQLquery,REQUEST_NAME)
            found_diseases = self.PartnerIndicationNetwork.get_entity_ids(['CellProcess'])
            print('Found %d indications for %d %ss of %s' %  
                 (len(found_diseases), len(self.partners),t_n,self.partner_class.lower()))


    def indications4cells_secreting_target(self):
        t_n = ','.join([x['Name'][0] for x in self.Drug_Targets])
        if self.target_class == 'Ligand':
            REQUEST_NAME = 'Find indications linked to cells secreting the {target}'.format(target=t_n)
            OQLquery = 'SELECT Relation WHERE objectType = (CellExpression,MolTransport) AND NeighborOf ({select_targets}) AND NeighborOf (SELECT Entity WHERE objectType = Cell)'
            cells_make_target = self.process_oql(OQLquery.format(select_targets=self.find_targets_oql),REQUEST_NAME)
            self.ProducingCellsIDs = cells_make_target.get_entity_ids(['CellType'])
            print('Found %d cell types producing %s' % (len(self.ProducingCellsIDs),t_n))

            REQUEST_NAME = 'Find diseases linked to cell secrteing {target}'
            REQUEST_NAME = REQUEST_NAME.format(target = t_n)
            effect = 'positive' if self.repurpose_antagonist else 'negative'
            OQLquery = 'SELECT Relation WHERE Effect = {effect} AND NeighborOf (SELECT Entity WHERE objectType = Disease) AND NeighborOf (SELECT Entity WHERE id = ({cell_ids}))'
            OQLquery = OQLquery.format(effect=effect,cell_ids=','.join(map(str, self.ProducingCellsIDs)))
            CellDiseaseNetwork = self.process_oql(OQLquery,REQUEST_NAME)

            found_diseases= CellDiseaseNetwork.get_entity_ids(['CellProcess'])
            print('Found %d indications for %d cells producing %s' %  
                 (len(found_diseases), len(self.ProducingCellsIDs),t_n))


    def pathway_oql(self, pathway_name_must_include_target=False):
        #set pathway_name_must_include_target to True if targets have a lot of large curated pathways
        if self.partners:
            SELECTpathways = 'SELECT Network WHERE ParentOf ({select_targets}) AND ParentOf ({select_partner})'#Name LIKE (\''+percent_sign+'FGF'+percent_sign+'\')'
            SELECTpathways = SELECTpathways.format(select_targets=self.find_targets_oql,select_partner=self.find_partners_oql)
        else:
            SELECTpathways = 'SELECT Network WHERE ParentOf ({select_targets})'#Name LIKE (\''+percent_sign+'FGF'+percent_sign+'\')'
            SELECTpathways = SELECTpathways.format(select_targets=self.find_targets_oql)
        
        if pathway_name_must_include_target:
            pct = '%'
            oql_queries = dict()
            for target in self.Drug_Targets:
                drug_target_name = target['Name'][0]
                SELECTpathways4target = SELECTpathways + ' AND Name LIKE (\''+pct+drug_target_name+pct+'\')' #additional refinement for popular targets
                SELECTMergedPathway = 'SELECT Relation WHERE objectType = (DirectRegulation,Binding,ProtModification,PromoterBinding,ChemicalReaction) AND MemberOf ({select_networks})'
                oql_queries[drug_target_name] = (SELECTMergedPathway.format(select_networks=SELECTpathways4target))

            return oql_queries
        else:
            target_names = ','.join([x['Name'][0] for x in self.Drug_Targets])
            return {target_names:SELECTpathways}
    

    def get_pathway_componets(self):
        #finding downstream pathway components
        oql_queries = self.pathway_oql() # separate oql_qury for each target
        partners_names = [p['Name'][0] for p in self.partners]
        REQUEST_NAME = 'Find curated pathways containing {targets} and {partners}'
        TargetsRegulome = ResnetGraph()
        self.PathwayComponentsIDs = list()
        target_names_with_regulome = list()
        for target_name, oql_query in oql_queries.items():
            REQUEST_NAME = REQUEST_NAME.format(targets=target_name,partners=','.join(partners_names))
            MergedPathway = self.process_oql(oql_query,REQUEST_NAME)
            if MergedPathway:
                TargetsRegulome.add_graph(MergedPathway.get_regulome(self.target_ids))
                self.PathwayComponentsIDs.append(list(MergedPathway.nodes))
                target_names_with_regulome.append(target_name)

        t_n = ','.join(target_names_with_regulome)
        print ('Found %s regulome with %d components' %(t_n,len(self.PathwayComponentsIDs)))
        return TargetsRegulome


    def init_semantic_search(self):
        print('Initializing semantic search')
        t_n = ','.join([x['Name'][0] for x in self.Drug_Targets])
        DiseaseNames = [y['Name'][0] for x,y in self.Graph.nodes(data=True) if y['ObjTypeName'][0] in ['CellProcess']]
        DiseaseScores = pd.DataFrame()
        DiseaseScores['Name'] = np.array(DiseaseNames)
        print('Will score %d diseases linked to %s' % (len(DiseaseScores),t_n))
        self.load_pandas(DiseaseScores,prop_names_in_header=True,map2type=['CellProcess'])


    def score_GVs(self):
        t_n = ','.join([x['Name'][0] for x in self.Drug_Targets])
        self.__colnameGV__ = t_n+' GVs'
        allGVids = self.Graph.get_entity_ids(['GeneticVariant'])
        gvlinkcounter = 0
        for i in self.RefCountPandas.index:
            row_entity_ids = list(self.RefCountPandas.at[i,self.__temp_id_col__])
            GVneighborsId = self.GVsInDiseaseNetwork.get_neighbors(row_entity_ids, allGVids)
            GVscore = 0
            if len(GVneighborsId) > 0:
                GVnames = set([n[0] for i,n in self.GVsInDiseaseNetwork.nodes.data('Name') if i in GVneighborsId])
                self.RefCountPandas.at[i,self.__colnameGV__] = ';'.join(GVnames)
                gvlinkcounter += 1
                gv_disease_subgraph = self.GVsInDiseaseNetwork.get_subgraph(row_entity_ids,GVneighborsId)
                GVscore = len(gv_disease_subgraph.count_references())
            
            self.RefCountPandas.at[i,self._col_name_prefix+'GVs'] = GVscore

        print('Found %d diseases linked to %d GVs' % (gvlinkcounter, len(allGVids)))


    def score_semantics(self):
        t_n = ','.join([x['Name'][0] for x in self.Drug_Targets])
        #references supporting target modulation of disease
        colname = 'Modulated by '+ t_n
        effect = 'positive' if self.repurpose_antagonist else 'negative'
        if self.strict_mode:
            self.set_how2connect(['Regulation'],[],'')
        else:
            self.set_how2connect(['Regulation'],[effect],'')
        linked_entities_count = self.link2concept(colname,self.target_ids)
        print('%d diseases are %sly regulated by %s' % (linked_entities_count,effect,t_n))
    
        #references where target expression or activity changes in the disease
        colname = t_n + ' is upregulated'
        effect = 'positive' if self.repurpose_antagonist else 'negative'
        if self.strict_mode:
            self.set_how2connect(['Regulation'],[],'')
        else:
            self.set_how2connect(['Regulation'],[effect],'')
        
        linked_entities_count= self.link2concept(colname,self.target_ids)
        print('%d diseases %sly regulate %s' % (linked_entities_count,effect,t_n))

        #references suggested target as a biomarker for indication
        colname = t_n + ' is Biomarker'
        self.set_how2connect(['Biomarker'],[],'')
        linked_entities_count = self.link2concept(colname,self.target_ids)
        print('Found %d diseases where %s is biomarker' % (linked_entities_count,t_n))

        self.score_GVs()

        if hasattr(self,'TargetInhibitorsIDs'):
            #references suggesting target antagonists as treatments for indication 
            colname = t_n + ' antagonists'
            if self.repurpose_antagonist:
                effect = 'negative'
                drug_class = 'antagonists'
            else:
                effect = 'positive'
                drug_class = 'agonists'
            self.set_how2connect(['Regulation'],[effect],'')
            linked_entities_count = self.link2concept(colname,self.TargetInhibitorsIDs)
            print('Found %d diseases are indications for %s %s' % (linked_entities_count,t_n,drug_class))

        #references suggesting target partners as targets for indication
        if hasattr(self, 'PartnerIndicationNetwork'):
            p_cl = self.partner_class
            colname = '{target_name} {partnet_class}s'.format(target_name=t_n,partnet_class=p_cl)
            PartnerIds = self.PartnerIndicationNetwork.get_entity_ids(['Protein'])
            
            effect = 'positive' if self.repurpose_antagonist else 'negative'
            self.set_how2connect(['Regulation'],[effect],'')

            linked_entities_count = self.link2concept(colname,PartnerIds)
            print('Found %d diseases indications for %d %s %ss' % (linked_entities_count,len(PartnerIds),t_n,p_cl))

        #references reporting that cells producing the target linked to indication
        if hasattr(self, 'ProducingCellsIDs'):
            colname = '{target_name} producing cells'.format(target_name=t_n)
            effect = 'positive' if self.repurpose_antagonist else 'negative'
            if self.strict_mode:
                self.set_how2connect(['Regulation'],[],'')
            else:
                self.set_how2connect(['Regulation'],[effect],'')

            linked_entities_count = self.link2concept(colname,self.ProducingCellsIDs)
            print('Found %d diseases linked %d cells producing %s' % (linked_entities_count,len(self.ProducingCellsIDs),t_n))
        
        #references reporting target agonists exacerbating indication or causing indication as adverse events
        if hasattr(self, 'TargetAgonistIDs'):
            colname = t_n + ' agonists'
            if self.repurpose_antagonist:
                effect = 'positive'
                drug_class = 'agonists'
            else:
                effect = 'negative'
                drug_class = 'antagonists'

            print('Finding %s for %s'%(drug_class, t_n))
            self.set_how2connect(['Regulation'],[effect],'')
            linked_entities_count = self.link2concept(colname,self.TargetAgonistIDs)
            print('Found %d diseases are toxicities for %s %s' % (linked_entities_count,t_n,drug_class))
        
        if hasattr(self, 'PathwayComponentsIDs'):
            #references linking target pathway to indication
            colname = t_n + ' pathway components'
            if self.strict_mode: effect = ''
            else: effect = 'positive' if self.repurpose_antagonist else 'negative'
            if self.strict_mode:
                self.set_how2connect(['Regulation'],[],'')
            else:
                self.set_how2connect(['Regulation'],[effect],'')

            linked_entities_count = self.link2concept(colname,self.PathwayComponentsIDs)
            print('%d indications are linked to %s pathway components' % (linked_entities_count,t_n))
        

    def normalize_counts(self):
        NormalizedCount = pd.DataFrame()
        weights = pd.DataFrame()
        refcount_cols = [col for col in self.RefCountPandas.columns if 'RefCount' in col]
        number_of_weights = len(refcount_cols)
        
        NormalizedCount['Name'] = self.RefCountPandas['Name']
        weights.at[0,'Name'] = 'WEIGHTS:'
        weight_index = 0
        for col in refcount_cols:
            col_max = self.RefCountPandas[col].max()
            NormalizedCount[col] = self.RefCountPandas[col]/col_max if col_max > 0 else 0.0
            column_weight = (number_of_weights-weight_index)/number_of_weights
            weights.at[0,col] = column_weight
            weight_index += 1
     
        #calculating cumulative score     
        combined_scores = list()
        for i in NormalizedCount.index:
            scores_row = NormalizedCount.loc[[i]]
            weighted_sum = 0.0
            for col in refcount_cols:
                weighted_sum = weighted_sum + scores_row[col]*weights.at[0,col]
            combined_scores.append(weighted_sum)

        NormalizedCount['Combined score'] = np.array(combined_scores)
        try:
            NormalizedCount[self.__colnameGV__] = self.RefCountPandas[self.__colnameGV__]
        except KeyError:
            pass

        NormalizedCount['#children'] = self.RefCountPandas[self.__temp_id_col__].apply(lambda x: len(x))
        NormalizedCount['Final score'] = NormalizedCount['Combined score']/NormalizedCount['#children']
        NormalizedCount = NormalizedCount.sort_values(by=['Final score'],ascending=False)
        NormalizedCount = NormalizedCount.loc[NormalizedCount['Final score'] > 0.0] # removes rows with all zeros
        try:
            NormalizedCount = NormalizedCount.loc[NormalizedCount[self.__colnameGV__].notna()]
        except KeyError: pass
        
        return weights.append(NormalizedCount)


if __name__ == "__main__":
    start_time = time.time()
    APIconfig = load_api_config()
    rd = RepurposeDrugs(APIconfig)
    #rd.set_targets(['IL15'], 'Protein',to_inhibit=True)
    #rd.set_targets(['FGFR3'], 'Protein',to_inhibit=True)
    #rd.set_targets(['F2RL1'], 'Protein',to_inhibit=True)
    #rd.set_targets(['MRGPRX2'], 'Protein',to_inhibit=True)
    rd.set_targets(['CNR2'], 'Protein', to_inhibit=False, strict_mode=True)
    rd.PageSize = 1000
    rd.flush_dump_files()

    rd.find_target_indications()
    rd.get_pathway_componets()

    if rd.target_class != 'Ligand':
    # if target is Ligand use indications4chem_modulators only if its antibody drugs have relations in Resnet
        rd.indications4chem_modulators() 
        rd.counterindications4chem_antimodulators()

    rd.indications4partners()
    rd.indications4cells_secreting_target() #will work only if target is Ligand
    
    rd.init_semantic_search()
    rd.score_semantics()
    t_n = ','.join([x['Name'][0] for x in rd.Drug_Targets])
    count_file = t_n+" repurposing counts.tsv"
    ref_file = t_n+" repurposing references.tsv"
    rd.print_ref_count(count_file,referencesOut=ref_file,sep='\t')

    NormalizedCount = rd.normalize_counts()
    fout = t_n + ' repurposing normalized report.tsv'
    NormalizedCount.to_csv(fout, sep='\t', index=False,float_format='%g')
    print('Repurposing of %s was done in %s' % (t_n, rd.execution_time(start_time)))
    print ('Ranked indications are in %s' % fout)
    print ('Semantic counts for each indication are in %s' % count_file)
    print('References supporting semantic counts are in %s' % ref_file)
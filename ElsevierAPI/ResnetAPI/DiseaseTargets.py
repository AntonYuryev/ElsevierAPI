from .ResnetGraph import PROTEIN_TYPES,EFFECT,PHYSICAL_INTERACTIONS, ResnetGraph
from .SemanticSearch import SemanticSearch,len,OQL,COUNTS
from .ResnetAPISession import NO_REL_PROPERTIES,ONLY_REL_PROPERTIES,REFERENCE_IDENTIFIERS,BIBLIO_PROPERTIES
from .FolderContent import FolderContent
from ..pandas.panda_tricks import df
import time

RANK_KNOWN_TARGETS = True
# in RANK_KNOWN_TARGETS mode only tarets suggested in the lietarure are ranked by amount of supporting evidence
PREDICT_TARGETS = False
# mode also predicts and ranks targets that are either biomarkers or adjacent possible

# names for raw data dfs
DF4ANTAGONISTS = 'Targets4antag'
DF4AGONISTS = 'Targets4agon'
DF4UNKNOWN = 'OtherTargets'

# names for report worksheets
ANTAGONIST_TARGETS_WS = 'targets4inhibitors'
AGONIST_TARGETS_WS = 'targets4agonists'
UNKNOWN_TARGETS_WS = 'unknown_state_targets'

PATHWAY_REGULATOR_SCORE = 'Disease model regulator score'
DRUG2TARGET_REGULATOR_SCORE = 'Targets regulator score'

class DiseaseTargets(SemanticSearch):
    pass
    def __init__(self,APIconfig,params=dict(),what2retrieve=BIBLIO_PROPERTIES):
        """
        self.set_target_disease_state() needs references
        """
        super().__init__(APIconfig,what2retrieve)
        self.add_rel_props([EFFECT])
        self.entProps = ['Name', 'Class'] #'Class' is for target partners retrieval
        self.symptoms_ids = set()
        self.clinpar_ids = set()
        self.cellproc_ids = set()
        self.disease_ids = list()
        self.targets4strictmode = set()
        if params:
            self.add_params(params)
        else:
            self.add_params({
                'disease':[],
                'strict_mode' : True,
                'data_dir' : '',
                'add_bibliography' : True,
                'strict_mode' : False,
                'target_types' : ['Protein'],
                'pathway_folders' : ['Hypertension Pulmonary'],
                'pathways' : []#['Myocardial Infarction']
                })
        # if strict_mode algorithm only ranks targets suggetsed in the literarure without predicting new targets
        # self.targets4strictmode contains targets linked to Disease by Regulation or targets with GVs linked to Disease
        # if not strict_mode algorithm also predicts and ranks additional targets based on target expression or genetic profiles in disease       
        self.columns2drop += [self.__resnet_name__,self.__mapped_by__,'State in Disease']


    def set_input(self):
        ontology_graph = self.child_graph(self.params['disease'],['Name','Alias'])
        self.diseases = ontology_graph._get_nodes()
        self.disease_ids = list(ontology_graph.nodes())
        self.find_disease_oql = OQL.get_objects(self.disease_ids)

    def input_names(self):
        if hasattr(self,"diseases"):
            prop2values = {'Name':self.params['disease'], 'Alias':self.params['disease']}
            return list(set([x['Name'][0] for x in self.diseases if x.has_value_in(prop2values)]))
        else:
            try:
                return self.params['sample']+' in '+ self.params['experiment']
            except KeyError:
                return ''

    def refcount_columns(self,counts_df=df(),column_prefix=''):
        to_return = super().refcount_columns(counts_df,column_prefix)
        if PATHWAY_REGULATOR_SCORE in counts_df.columns.tolist():
            to_return.append(PATHWAY_REGULATOR_SCORE)
        return to_return

    def _get_report_name(self,extension=''):
        rep_pred = 'suggested ' if self.params['strict_mode'] else 'suggested,predicted ' 
        return str(self.params['data_dir']+rep_pred+ 'targets for '+ ','.join(self.params['disease']))+extension

    def _disease2str(self):
        if hasattr(self,'diseases'):
            if len(self.diseases) > 1:
                return str(len(self.diseases))+' types of '+','.join(self.params['disease'])
            else:
                return ','.join(self.params['disease'])
        else:
            return ','.join(self.params['disease'])

    def __target_ids(self):
        return self.targets4strictmode if self.params['strict_mode'] else self.target_ids

    def __target_types_str(self):
        return ','.join(self.params['target_types'])
        
    def __GVtargets(self):
        disease_names = self._disease2str()

        oql_query = 'SELECT Relation WHERE objectType = FunctionalAssociation \
                    AND NeighborOf ({}) AND NeighborOf (SELECT Entity WHERE objectType = GeneticVariant)'
        oql_query = oql_query.format(self.find_disease_oql)
        request_name = 'Select GeneticVariants for {}'.format(disease_names)
        self.gv2diseases = self.process_oql(oql_query,request_name)
        self.GVids = self.gv2diseases.get_node_ids(['GeneticVariant'])

        oql_query = f'SELECT Relation WHERE objectType = GeneticChange AND \
        NeighborOf (SELECT Entity WHERE objectType = ({self.__target_types_str()})) AND NeighborOf ({OQL.get_objects(self.GVids)})'
        #oql_query = oql_query.format(self.__target_types_str(),OQL.get_objects(self.GVids))
        request_name = 'Find targets with genetic variants linked to {}'.format(disease_names)

        # avoid adding GeneticChange to self.Graph
        new_session = self._clone_(to_retrieve=NO_REL_PROPERTIES)
        gv2target = new_session.process_oql(oql_query,request_name)
        self.gv2diseases.add_graph(gv2target)

        target_withGVs_ids = gv2target.get_node_ids(PROTEIN_TYPES)
        gv2dis_refs = self.gv2diseases.load_references()
        print('Found %d targets with %d GeneticVariants for %s supported by %d references' % 
            (len(target_withGVs_ids), len(self.GVids), disease_names, len(gv2dis_refs)))
        
        self.targets4strictmode.update(target_withGVs_ids)
        return target_withGVs_ids


    def find_targets(self):
        REQUEST_NAME = 'Find targets linked to {}'.format(self._disease2str())
        select_targets = 'SELECT Entity WHERE objectType = ({})'.format(self.__target_types_str())
        OQLquery = 'SELECT Relation WHERE NeighborOf({}) AND NeighborOf ({})'      
        target_disease_graph = self.process_oql(OQLquery.format(select_targets, self.find_disease_oql,REQUEST_NAME))
        self.target_ids = set(target_disease_graph.get_node_ids(self.params['target_types']))
        target_refcount = target_disease_graph.load_references()
        print('Found %d targets linked to %s supported by %d references' 
                    % (len(self.target_ids),self._disease2str(), len(target_refcount)))
    
        self.target_ids.update(self.__GVtargets())

        disease_model_components_ids = set()
        disease_pathways = self.load_pathways()
        if disease_pathways:
            all_pathways = [p for p in disease_pathways.values()]
            disease_model_components = set()
            [disease_model_components.update(p.get_members(self.params['target_types'])) for p in all_pathways]
            disease_model_components_ids = {n.id() for n in disease_model_components}
            self.target_ids.update(disease_model_components_ids)
            before_add = self.Graph.number_of_nodes()
            self.Graph.add_nodes(list(disease_model_components))
            print('Added %d targets from disease model' % (self.Graph.number_of_nodes()-before_add))

        if self.params['strict_mode']:
            disease_regulators = target_disease_graph.subgraph_by_relprops(['Regulation'])
            regulators_ids = disease_regulators.get_node_ids(self.params['target_types'])
            self.targets4strictmode.update(regulators_ids)
            if disease_model_components_ids:
                self.targets4strictmode.update(disease_model_components_ids)


    def find_symptoms(self):
        if self.params['symptoms']:
            REQUEST_NAME = 'Find symptoms linked to {}'.format(self._disease2str())
            OQLquery = 'SELECT Relation WHERE objectType = FunctionalAssociation AND NeighborOf ({}) AND NeighborOf ({})'
            oql2select_symptom = OQL.get_childs(self.params['symptoms'],['Name'],include_parents=True)
            OQLquery = OQLquery.format(oql2select_symptom,self.find_disease_oql)
            cellproc_graph = self.process_oql(OQLquery,REQUEST_NAME)
            self.symptoms_ids = set(cellproc_graph.get_node_ids(['Disease']))
            symptoms_references = cellproc_graph.load_references()
            print('\n\nFound %d symptoms linked to %s supported by %d references' % 
                        (len(self.symptoms_ids),self._disease2str(), len(symptoms_references)))


    def find_clinpar(self):
        if self.params['clinical_parameters']:
            REQUEST_NAME = 'Find clnical parameters linked to {}'.format(self._disease2str())
            OQLquery = 'SELECT Relation WHERE objectType = FunctionalAssociation AND NeighborOf ({}) AND NeighborOf ({})'
            oql2select_clinpar = OQL.get_childs(self.params['clinical_parameters'],['Name'],include_parents=True)
            OQLquery = OQLquery.format(oql2select_clinpar,self.find_disease_oql)
            clinpar_graph = self.process_oql(OQLquery,REQUEST_NAME)
            self.clinpar_ids = set(clinpar_graph.get_node_ids(['ClinicalParameter']))
            cp_references = clinpar_graph.load_references()
            print('\n\nFound %d clinical paramaters linked to %s supported by %d references' % 
                        (len(self.clinpar_ids),self._disease2str(), len(cp_references)))


    def find_cellproc(self):
        if self.params['processes']:
            REQUEST_NAME = 'Find cell processes linked to {}'.format(self._disease2str())
            OQLquery = 'SELECT Relation WHERE objectType = FunctionalAssociation AND NeighborOf ({}) AND NeighborOf ({})'
            oql2select_cellproc = OQL.get_childs(self.params['processes'],['Name'],include_parents=True)
            OQLquery = OQLquery.format(oql2select_cellproc,self.find_disease_oql)
            cellproc_graph = self.process_oql(OQLquery,REQUEST_NAME)
            self.cellproc_ids = set(cellproc_graph.get_node_ids(['CellProcess']))
            cp_references = cellproc_graph.load_references()
            print('\n\nFound %d cell processes linked to %s supported by %d references' % 
                        (len(self.cellproc_ids),self._disease2str(), len(cp_references)))


    def find_drugs(self):      
        REQUEST_NAME = 'Find compounds inhibiting {}'.format(self._disease2str())
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = negative \
            AND NeighborOf (SELECT Entity WHERE objectType = SmallMol) AND NeighborOf ({})'
        OQLquery = OQLquery.format(self.find_disease_oql)
        drugs_graph = self.process_oql(OQLquery,REQUEST_NAME)
        self.drug_ids = set(drugs_graph.get_node_ids(['SmallMol']))
        drugs_references = drugs_graph.load_references()
        print('Found %d compounds inhibiting %s supported by %d references' % 
                    (len(self.drug_ids),self._disease2str(), len(drugs_references)))

        REQUEST_NAME = 'Find clinical trials for {}'.format(self._disease2str())
        OQLquery = 'SELECT Relation WHERE objectType = ClinicalTrial \
            AND NeighborOf (SELECT Entity WHERE objectType = SmallMol) AND NeighborOf ({})'
        OQLquery = OQLquery.format(self.find_disease_oql)
        ct_graph = self.process_oql(OQLquery,REQUEST_NAME)
        ct_drug_ids = ct_graph.get_node_ids(['SmallMol'])
        self.drug_ids.update(ct_drug_ids)
        ct_refs = ct_graph.load_references()
        print('Found %d compounds on %d clinical trilas for %s' % (len(ct_drug_ids),len(ct_refs),self._disease2str()))
       

    def find_inducers(self):
        REQUEST_NAME = 'Find compounds inducing/exacerbating {}'.format(self._disease2str())
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = positive \
            AND NeighborOf (SELECT Entity WHERE objectType = SmallMol) AND NeighborOf ({})'
        OQLquery = OQLquery.format(self.find_disease_oql)
        induction_graph = self.process_oql(OQLquery,REQUEST_NAME)
        self.inducers_ids = set(induction_graph.get_node_ids(['SmallMol']))
        inducers_refs = induction_graph.load_references()
        print('Found %d compounds inducing %s supported by %d references' % 
            (len(self.inducers_ids),self._disease2str(),len(inducers_refs)))

     
    def load_target_partners(self):
        receptor_ids = self.Graph.get_node_ids(['Receptor'],['Class'])
        ligand_ids = self.Graph.get_node_ids(['Ligand'],['Class'])
        #opening new APIsession to avoid adding result graph to self
        if receptor_ids:
            find_receptors = OQL.get_objects(receptor_ids)
            ligands_oql = 'SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive AND RelationNumberOfReferences > 5 AND \
                NeighborOf downstream (SELECT Entity WHERE Class = Ligand) AND NeighborOf upstream ({})'
            ligands_oql = ligands_oql.format(find_receptors)
            new_session = self._clone_(to_retrieve=NO_REL_PROPERTIES)
            request_name = 'Find ligands for {} receptors'.format(str(len(receptor_ids)))
            self.partners = new_session.process_oql(ligands_oql,request_name)
            new_ligands_ids = self.partners.get_node_ids(['Ligand'],['Class'])
            print('Found %d ligands for %d receptors' % (len(new_ligands_ids),len(receptor_ids)))
        else:
            print('No Receptor Class entities were found among %s targets' % self.params['target_types'])

        orphan_ligands = set(ligand_ids).difference(self.partners.nodes())
        if orphan_ligands:
            find_orphan_ligands = OQL.get_objects(list(orphan_ligands))
            receptor_oql = 'SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive AND RelationNumberOfReferences > 5 \
                AND NeighborOf upstream (SELECT Entity WHERE Class = Receptor) AND NeighborOf downstream ({})'
            receptor_oql = receptor_oql.format(find_orphan_ligands)
            request_name = 'Find receptors for {} orphan ligands'.format(str(len(orphan_ligands)))
            new_session = self._clone_(to_retrieve=NO_REL_PROPERTIES)
            orphan_ligand_partners = new_session.process_oql(receptor_oql,request_name)
            self.partners.add_graph(orphan_ligand_partners)
            new_receptors_ids = orphan_ligand_partners.get_node_ids(['Receptor'],['Class'])
            print('Found %d receptors for %d orphan ligands' % (len(new_receptors_ids),len(orphan_ligands)))
        
        
        find_metabolites = OQL.get_childs(['mammal endogenous compounds and their derivatives'],['Name'])
        metabolite_ligands_oql = 'SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive \
            AND NeighborOf downstream ({}) AND NeighborOf upstream ({}) AND RelationNumberOfReferences > 25'
        metabolite_ligands_oql = metabolite_ligands_oql.format(find_metabolites,find_receptors)
        new_session = self._clone_(to_retrieve=NO_REL_PROPERTIES)
        metabolite2receptors_graph = new_session.process_oql(metabolite_ligands_oql,'Find metabolite ligands')
        self.partners.add_graph(metabolite2receptors_graph)
        metabolite_ligand_ids = metabolite2receptors_graph.get_node_ids(['SmallMol'])
        print('Found %d metabolite ligands for %d receptors' % (len(metabolite_ligand_ids),len(receptor_ids)))
        
        
        find_metabolite_products = 'SELECT Relation WHERE objectType = (MolTransport, MolSynthesis,ChemicalReaction) \
            AND Effect = positive AND NeighborOf downstream (SELECT Entity WHERE id =({ids})) AND \
            NeighborOf upstream (SELECT Entity WHERE objectType = SmallMol) AND RelationNumberOfReferences > 50'
        
        req_name = f'Finding downstream metabolite products of targets for {self._disease2str()})'
        new_session = self._clone_(to_retrieve=NO_REL_PROPERTIES)
        metbolite_products_graph = new_session.iterate_oql(find_metabolite_products,self.__target_ids(),request_name=req_name)
        self.partners.add_graph(metbolite_products_graph)
        metabolite_products = metbolite_products_graph.get_node_ids(['SmallMol'])
        print('Found %d metabolites produced by %d targets of %s' % (len(metabolite_products),len(self.target_ids),self._disease2str()))

        # removing promiscous partners
        self.partners.remove_nodes_by_outdegree(max_degree=10,having_values=['SmallMol'])
        self.partners.remove_nodes_by_outdegree(max_degree=10,only_with_prop=['Class'], having_values=['Ligand'])
        self.partners.remove_nodes_by_indegree(max_degree=10,having_values=['SmallMol'])


    def load_pathways(self):
        """
        Loads
        -----
        self.disease_pathways = {CellType:PSPathway}
        PSPathway objects from folders listed in self.params['pathway_folders']
        keeps only pathways listed self.params['pathways'] if latter is not empty 

        Merges
        ------
        pathways with the same CellType

        Returns
        -------
        {CellType:PSPathway}
        where nodes in PSPathway are annotated with 'Closeness' attribute
        """ 
        fc = FolderContent(self.APIconfig,what2retrieve=ONLY_REL_PROPERTIES)
        # references are needed to calculate State Effect for targets in the model
        fc.entProps = ['Name','CellType','Tissue','Organ','Organ System']
        #fc.add_rel_props([EFFECT])
        # effect is needed to calculate State Effect for targets in the model

        disease_pathways = list()
        for folder_name in self.params['pathway_folders']:
            ps_pathways = fc.folder2pspathways(folder_name)
            disease_pathways += ps_pathways

        if self.params['pathways']:
            filtered_pathway_list = list()
            for ps_pathway in disease_pathways:
                if ps_pathway.name() in self.params['pathways']:
                    filtered_pathway_list.append(ps_pathway)
            disease_pathways = filtered_pathway_list

        self.disease_pathways = dict()
        # merging pathways from the same celltype
        for pathway in disease_pathways:
            try:
                cell_types = pathway['CellType']
            except KeyError:
                cell_types = ['disease']

            for cell_type in cell_types:
                try:
                    exist_pathway = self.disease_pathways[cell_type]
                    exist_pathway.merge_pathway(pathway)
                    self.disease_pathways[cell_type] = exist_pathway
                except KeyError:
                    self.disease_pathways[cell_type] = pathway
            
        for cell_type, pathway in self.disease_pathways.items():
            pathway.graph.remove_nodes_by_prop(['CellProcess', 'Disease','Treatment'])

        [p.graph.closeness() for p in self.disease_pathways.values()]
        print('Found %d curated pathways for %s:' %(len(disease_pathways), self._disease2str()))
        [print(p.name()+'\n') for p in disease_pathways]
        return self.disease_pathways


    def init_semantic_search(self):
        print('\n\nInitializing semantic search')
        target_ids = self.__target_ids()
        if not target_ids:
            print ('No targets found for %s' % self._disease2str())
            print('Consider setting "strict_mode" to False')
            return False

        target_names = [y['Name'][0] for x,y in self.Graph.nodes(data=True) if y['Id'][0] in target_ids]
        target2score_df = df.from_dict({'Name':target_names},name=COUNTS)

        print('Will score %d targets linked to %s' % (len(target2score_df),self._disease2str()))
        self.RefCountPandas = self.load_pandas(target2score_df,prop_names_in_header=True,
        map2type=self.params['target_types'],max_children_count=11)
        self.entProps = ['Name'] # Class is no longer needed
        return True
        

    def score_GVs(self):
        print('\n\nScoring targets by number of semantic references linking their Genetic Variants to %s' % self._disease2str())
        if not self.GVids:
            print('%s has no known Genetic Variants' % self._disease2str())
            return

        self.__colnameGV__ = 'Genetic Variants for '+self._disease2str()
        target_gvlinkcounter = 0
        for i in self.RefCountPandas.index:
            target_ids = list(self.RefCountPandas.at[i,self.__temp_id_col__])
            targetGV_ids = self.gv2diseases.get_neighbors(target_ids, only_with_ids=self.GVids)
                
            GVscore = 0
            if targetGV_ids:
                GVnames = set([n[0] for i,n in self.Graph.nodes.data('Name') if i in targetGV_ids])
                self.RefCountPandas.at[i,self.__colnameGV__] = ';'.join(GVnames)
                target_gvlinkcounter += 1
                gv_disease_subgraph = self.gv2diseases.get_subgraph(targetGV_ids,self.disease_ids)
                GVscore = len(gv_disease_subgraph.load_references())
            
            self.RefCountPandas.at[i,self._col_name_prefix+'Genetic Variants'] = GVscore

        print('Found %d targets linked to %d GVs' % (target_gvlinkcounter, len(self.GVids)))


    def __vote4effect(self,target_ids:list):
        between_graph = self.Graph.get_subgraph(target_ids, self.disease_ids)
        positive_refs, negative_refs = between_graph._effect_counts__()

        if len(positive_refs) > len(negative_refs):
            return 'activated'
        elif len(negative_refs) > len(positive_refs):
            return 'repressed'
        else:
            target_partners_ids = set(self.partners.get_neighbors(target_ids))
            between_graph = self.Graph.get_subgraph(target_partners_ids, self.disease_ids)
            positive_refs, negative_refs = between_graph._effect_counts__()
            if len(positive_refs) > len(negative_refs):
                return 'activated'
            elif len(negative_refs) > len(positive_refs):
                return 'repressed'
            else:
                # attempting to deduce effect from GeneticChange
                genetic_change = between_graph.get_relations(['GeneticChange'])
                if genetic_change:
                    return 'repressed'
                else:
                    return 'unknown'


    def set_target_disease_state(self):
        print('\n\nCalculating targets state (activated/repressed) in %s' % self._disease2str())
                       
        for i in self.RefCountPandas.index:
            target_ids = list(self.RefCountPandas.at[i,self.__temp_id_col__])
            self.RefCountPandas.at[i,'State in Disease'] = self.__vote4effect(target_ids)

        model_graph = ResnetGraph()
        for pathway in self.disease_pathways.values():
            model_graph = model_graph.compose(pathway.graph)

        def __disease_state_from_model():
            activ_targets_pd = self.RefCountPandas.loc[(self.RefCountPandas['State in Disease'] == 'activated')]
            inhib_targets_pd = self.RefCountPandas.loc[(self.RefCountPandas['State in Disease'] == 'repressed')]
            activated_target_ids = self._all_ids(activ_targets_pd)
            inhibited_target_ids = self._all_ids(inhib_targets_pd)

            new_regulators_count = 0
            for idx in self.RefCountPandas.index:
                if self.RefCountPandas.at[idx,'State in Disease'] == 'unknown':
                    net_effect = 0
                    for i in self.RefCountPandas.at[idx,self.__temp_id_col__]:
                        activates_activated_targets,inhibits_activated_targets = model_graph.net_regulator_effect(i,activated_target_ids)
                        activates_inhibited_targets,inhibits_inhibited_targets = model_graph.net_regulator_effect(i,inhibited_target_ids)
                        net_effect += len(activates_activated_targets)+len(inhibits_inhibited_targets)-len(inhibits_activated_targets)-len(activates_inhibited_targets)

                    if net_effect > 0:
                        self.RefCountPandas.at[idx,'State in Disease'] = 'activated'
                        new_regulators_count += 1
                    elif net_effect < 0:
                        self.RefCountPandas.at[idx,'State in Disease'] = 'repressed'
                        new_regulators_count += 1

            return new_regulators_count
        
        new_regulators_count = __disease_state_from_model()
        while new_regulators_count:
           new_regulators_count =  __disease_state_from_model()

        return



    def score_partners(self):
        print('\n\nFinding semantic references supporting links between target partners and %s'%
                 self._disease2str())
        self.RefCountPandas['target_ids'] = self.RefCountPandas[self.__temp_id_col__]
        for i in self.RefCountPandas.index:
            target_ids = list(self.RefCountPandas.at[i,'target_ids'])
            target_partners_ids = set(self.partners.get_neighbors(target_ids))
            if target_partners_ids:
                self.RefCountPandas.at[i,self.__temp_id_col__] = tuple(target_partners_ids)
            else:
                self.RefCountPandas.at[i,self.__temp_id_col__] = tuple([0])
            partner_names = [n[0] for i,n in self.partners.nodes(data='Name') if i in set(target_partners_ids)]
            self.RefCountPandas.at[i,'Target partners'] = ';'.join(partner_names)
        
        self.set_how2connect([],[],'')
        colname = 'target partners'
        linked_row_count = self.link2RefCountPandas(colname,self.disease_ids)
        print('%d targets have partners linked to %s' % (linked_row_count,self._disease2str()))

        self.RefCountPandas.drop(columns=self.__temp_id_col__,inplace=True)
        self.RefCountPandas.rename(columns={'target_ids':self.__temp_id_col__}, inplace=True)


    def score_regulators(self):
        print('\n\nScoring regulators by distance to components of disease pathways',flush=True)
        print('Retrieving regulatory network between targets and components of disease pathway')
        disease_pathway_component_ids = set()
        regulation_graph = ResnetGraph()
        [regulation_graph.add_graph(p.graph) for p in self.disease_pathways.values()]
        disease_pathway_component_ids = set(regulation_graph.nodes())

        unconnected_nodes = self.all_entity_ids
        connected_nodes = disease_pathway_component_ids
        for step in range (0,5):
            new_session = self._clone_(to_retrieve=NO_REL_PROPERTIES)
            graph_expansion = new_session.connect_nodes(unconnected_nodes, connected_nodes,
                                                    PHYSICAL_INTERACTIONS, in_direction='>')
            regulation_graph.add_graph(graph_expansion)

            connected_nodes = set(regulation_graph.nodes()).difference(disease_pathway_component_ids)
            connected_nodes = connected_nodes.intersection(unconnected_nodes)
            # only nodes connected at the previous cycle need to be expanded at the next cycle
            unconnected_nodes = set(unconnected_nodes).difference(regulation_graph.nodes())

        for pathway in self.disease_pathways.values():
            closeness_dic = {i:c for i,c in pathway.graph.nodes(data='Closeness')}
            regulation_graph.rank_regulators(closeness_dic,PATHWAY_REGULATOR_SCORE)

        for i in self.RefCountPandas.index:
            target_ids = list(self.RefCountPandas.at[i,self.__temp_id_col__])
            target_regulator_score = 0.0
            for target_id in target_ids:
                try:
                    regulator_score = regulation_graph.nodes[target_id][PATHWAY_REGULATOR_SCORE]
                    target_regulator_score += regulator_score
                except KeyError:
                    continue
            self.RefCountPandas.at[i,PATHWAY_REGULATOR_SCORE] = target_regulator_score


    def score_target_semantics(self):
        t_n = self._disease2str()
        self.score_GVs()

        colname = 'Regulated by '+ t_n
        self.set_how2connect(['Regulation'],[],'')
        linked_row_count = self.link2RefCountPandas(colname,self.disease_ids)
        print('%d targets regulating %s' % (linked_row_count,t_n))

        self.set_target_disease_state()

        colname = 'Genetically linked to '+ t_n
        self.set_how2connect(['GeneticChange'],[],'')
        linked_row_count = self.link2RefCountPandas(colname,self.disease_ids)
        print('%d targets genetically linked to %s' % (linked_row_count,t_n))

        colname = 'Target is Biomarker in '+t_n
        self.set_how2connect(['Biomarker'],[],'')
        linked_row_count = self.link2RefCountPandas(colname,self.disease_ids)
        print('Linked %d indications where %s is biomarker' % (linked_row_count,t_n))

        colname = 'Quantitatively changed in '+ t_n
        self.set_how2connect(['QuantitativeChange'],[],'')
        linked_row_count = self.link2RefCountPandas(colname,self.disease_ids)
        print('%d targets quantitatively changed in %s' % (linked_row_count,t_n))

        self.score_partners()

        if self.symptoms_ids:
            colname = 'symptoms for '+ t_n
            reltypes2connect = ['Regulation','QuantitativeChange','StateChange','Biomarker']
            self.set_how2connect(reltypes2connect,[],'',how2clone=REFERENCE_IDENTIFIERS)
            linked_row_count = self.link2RefCountPandas(colname,self.symptoms_ids)
            print('%d targets linked to symptoms for in %s' % (linked_row_count,t_n))

        if self.clinpar_ids:
            colname = 'Clinical parameters for '+ t_n
            self.set_how2connect(['Regulation'],[],'',how2clone=REFERENCE_IDENTIFIERS)
            linked_row_count = self.link2RefCountPandas(colname,self.clinpar_ids)
            print('%d targets linked to clinical parameters for in %s' % (linked_row_count,t_n))

        self.score_regulators()

        if self.cellproc_ids:
            colname = 'Cell processes affected by '+ t_n
            self.set_how2connect(['Regulation'],[],'',how2clone=REFERENCE_IDENTIFIERS)
            linked_row_count = self.link2RefCountPandas(colname,self.cellproc_ids)
            print('%d targets linked to cell processes for %s' % (linked_row_count,t_n))

        ################## SPLITTING RefCountPandas to score agonists and antagonists differently #####################
        activ_targets_pd = self.RefCountPandas.loc[(self.RefCountPandas['State in Disease'] == 'activated')]
        inhib_targets_pd = self.RefCountPandas.loc[(self.RefCountPandas['State in Disease'] == 'repressed')]
        unk_targets_pd = self.RefCountPandas.loc[(self.RefCountPandas['State in Disease'] == 'unknown')]

        if hasattr(self, 'drug_ids'):
            all_compound_ids = set()
            # drug-target relations are not added to self.Graph to exclude them from bibliography
            colname = 'drugs for '+t_n
            self.set_how2connect(['DirectRegulation'],['negative'],'',['Binding','Expression','Regulation'],REFERENCE_IDENTIFIERS)
            print ('\n\nLinking with Effect negative drugs for "%s" to targets activated in "%s"' % (t_n,t_n))
            was_linked,linked_ids,activ_targets_pd = self.link2concept(colname,self.drug_ids,activ_targets_pd)
            print('Linked %d targets to drugs for "%s"' % (was_linked,t_n))

            self.set_how2connect(['DirectRegulation'],['positive'],'',['Binding','Expression','Regulation'],REFERENCE_IDENTIFIERS)
            print ('\n\nLinking with Effect positive drugs for "%s" to targets inhibiting "%s"' % (t_n,t_n))
            was_linked,linked_ids,inhib_targets_pd = self.link2concept(colname,self.drug_ids,inhib_targets_pd)
            print('Linked %d targets to drugs for "%s"' % (was_linked,t_n))
            all_compound_ids.update(self.drug_ids)

        if hasattr(self, 'inducers_ids'):
            colname = 'inducers of '+t_n
            self.set_how2connect(['DirectRegulation'],['positive'],'',['Binding','Expression','Regulation'],REFERENCE_IDENTIFIERS)
            print ('\n\nLinking with Effect positive compounds inducing "%s" to targets activated in "%s"' % (t_n,t_n))
            was_linked,linked_ids,activ_targets_pd = self.link2concept(colname,self.inducers_ids,activ_targets_pd)
            print('Linked %d targets to inducers of "%s"' % (was_linked,t_n))

            self.set_how2connect(['DirectRegulation'],['negative'],'',['Binding','Expression','Regulation'],REFERENCE_IDENTIFIERS)
            print ('\n\nLinking with effect negative compounds inducing "%s" to targets inhibiting "%s"' % (t_n,t_n))
            was_linked,linked_ids,inhib_targets_pd = self.link2concept(colname,self.inducers_ids,inhib_targets_pd)
            print('Linked %d targets to inducers of "%s"' % (was_linked,t_n))
            all_compound_ids.update(self.inducers_ids)

        if all_compound_ids:
            colname = 'compounds modulating '+t_n
            print ('\n\nLinking compounds modulating "%s" to targets with unknown state in "%s"' % (t_n,t_n))
            self.set_how2connect(['DirectRegulation'],[],'',['Binding','Expression','Regulation'],REFERENCE_IDENTIFIERS)
            was_linked,linked_ids,unk_targets_pd = self.link2concept(colname,all_compound_ids,unk_targets_pd)
            print('Linked %d targets to inducers for %s' % (was_linked,t_n))

        uptarget_pd = self.make_count_df(activ_targets_pd,DF4ANTAGONISTS)
        self.add2raw(uptarget_pd)
        downtarget_pd = self.make_count_df(inhib_targets_pd,DF4AGONISTS)
        self.add2raw(downtarget_pd)
        unktarget_pd = self.make_count_df(unk_targets_pd,DF4UNKNOWN)
        self.add2raw(unktarget_pd)
        return


    def normalize_counts(self):
        self.normalize(DF4ANTAGONISTS,ANTAGONIST_TARGETS_WS)
        self.normalize(DF4AGONISTS,AGONIST_TARGETS_WS)
        self.normalize(DF4UNKNOWN,UNKNOWN_TARGETS_WS)


    def add_etm_bibliography(self):
        '''
        Adds
        ----
        etm references column to ANTAGONIST_TARGETS_WS,AGONIST_TARGETS_WS,UNKNOWN_TARGETS_WS worksheets,\n
        adds ETM_BIBLIOGRAPHY worksheet to report
        '''
        print('Adding ETM bibliography for ranked targets', flush=True)
        self.add_etm_refs(ANTAGONIST_TARGETS_WS,self.input_names())
        self.add_etm_refs(AGONIST_TARGETS_WS,self.input_names())
        self.add_etm_refs(UNKNOWN_TARGETS_WS,self.input_names())
        super().add_etm_bibliography()


    def make_report(self):
        start_time = time.time()
        self.set_input()
        self.flush_dump_files()
        self.find_targets()
        self.load_target_partners()

        self.find_symptoms()
        self.find_clinpar()
        self.find_cellproc()

        self.find_drugs()
        self.find_inducers()
        
        if self.init_semantic_search():
            self.score_target_semantics()

        self.normalize_counts()
        self.add_ps_bibliography()
    
        print('Target ranking for %s was done in %s' % 
        (self._disease2str(), self.execution_time(start_time)))
        

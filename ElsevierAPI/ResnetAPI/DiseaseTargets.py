from .ResnetGraph import PROTEIN_TYPES,EFFECT,PHYSICAL_INTERACTIONS, ResnetGraph, CLOSENESS
from .SemanticSearch import SemanticSearch,len,OQL,TO_RETRIEVE
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
    def __init__(self,*args,**kwargs):
        """
        Input
        -----
        APIconfig - args[0]
        self.set_target_disease_state() needs references
        """
        APIconfig=args[0]
        my_kwargs = {
                'disease':[],
                'what2retrieve':BIBLIO_PROPERTIES,
                'strict_mode' : True,
                'data_dir' : '',
                'add_bibliography' : True,
                'strict_mode' : False,
                'target_types' : ['Protein'],
                'pathway_folders' : ['Hypertension Pulmonary'],
                'pathways' : []
            }
        my_kwargs.update(kwargs)
        super().__init__(APIconfig,**my_kwargs)
        self.add_rel_props([EFFECT])
        self.entProps = ['Name', 'Class'] #'Class' is for target partners retrieval
        self.columns2drop += [self.__resnet_name__,self.__mapped_by__,'State in Disease']
        self.max_threads4ontology = 50

        self.input_diseases = list()
        self.input_symptoms = set()
        self.input_clinpars = set()
        self.input_cellprocs = set()

        self.__targets__ = set()
        self.targets4strictmode = set()
        self.GVs = list()
        self.gv2diseases = ResnetGraph()
        self.partner2targets = ResnetGraph()
        self.drugs_linked2disease = ResnetGraph()
        self.disease_inducers = list()
        self.disease_pathways = dict() # {cell_type:PSPathway}

        # self.targets4strictmode has PSObjects for targets directly Regulating disease or with GVs linked to disease or from disease pathways
        # if strict_mode algorithm only ranks targets suggetsed in the literarure without predicting new targets
        # in not strict_mode algorithm also predicts and ranks additional targets linked to disease by QuantitativeChange or GeneticChange     
       
        # TBD cache networks:
        self.prot2Dis = ResnetGraph() # includes GVs
        self.prot2ClinPar = ResnetGraph() # includes GVs
        self.protein2CelProc = ResnetGraph() # protein-CellProcess regulation network (Regulation) to find partners
        self.ppi  = ResnetGraph() # protein-protein physical interactions (Binding,DirectRegulation,ProteModification) for self.partner2targets
        self.ppmet = ResnetGraph() # protein-metabolite regulation network (MolTransport,MolSynthesis,ChemicalReaction) for self.partner2targets
        self.drug_effects = ResnetGraph() # drug-disease regulation network (Regulation) disease_inducers, drugs_link2disease
        self.drugs2targets = ResnetGraph() # for disease_inducers, drugs_link2disease
    
    def clear(self):
        super().clear()
        

    def set_input(self):
        ontology_graph = self.child_graph(self.params['disease'],['Name','Alias'])
        self.input_diseases = ontology_graph._get_nodes()
        input_diseases_dbids = ResnetGraph.dbids(self.input_diseases)
        self.find_disease_oql = OQL.get_objects(input_diseases_dbids)


    def input_names(self):
        if self.input_diseases:
            prop2values = {'Name':self.params['disease'], 'Alias':self.params['disease']}
            return list(set([x['Name'][0] for x in self.input_diseases if x.has_value_in(prop2values)]))
        else:
            try:# case of SNEA drug_df
                return self.params['sample']+' in '+ self.params['experiment']
            except KeyError:
                return 'Drugs'
            

    def refcount_columns(self,counts_df=df(),column_prefix=''):
        to_return = super().refcount_columns(counts_df,column_prefix)
        if PATHWAY_REGULATOR_SCORE in counts_df.columns.tolist():
            to_return.append(PATHWAY_REGULATOR_SCORE)
        return to_return


    def report_name(self):
        rep_pred = ' suggested ' if self.params['strict_mode'] else ' predicted ' 
        return ','.join(self.params['disease'])+rep_pred+'targets'


    def _disease2str(self):
        # do not put disease name in quotas - it is used to create a file name
        if len(self.input_diseases) > 1:
                return ','.join(self.params['disease'])+' ('+str(len(self.input_diseases))+' types)'
        else:
            return ','.join(self.params['disease'])


    def _targets(self):
        return self.targets4strictmode if self.params['strict_mode'] else self.__targets__
    

    def _targets_dbids(self):
        return ResnetGraph.dbids(self._targets())

    def __target_types_str(self):
        return ','.join(self.params['target_types'])
        
    def __GVtargets(self):
        disease_names = self._disease2str()

        oql_query = f'SELECT Relation WHERE objectType = FunctionalAssociation \
            AND NeighborOf ({self.find_disease_oql}) AND NeighborOf (SELECT Entity WHERE objectType = GeneticVariant)'
        #oql_query = oql_query.format()
        request_name = f'Select GeneticVariants for {disease_names}'
        self.gv2diseases = self.process_oql(oql_query,request_name)
        self.GVs = self.gv2diseases.psobjs_with(only_with_values=['GeneticVariant'])
        GVdbids = ResnetGraph.dbids(self.GVs)

        oql_query = f'SELECT Relation WHERE objectType = GeneticChange AND \
            NeighborOf (SELECT Entity WHERE objectType = ({self.__target_types_str()})) AND NeighborOf ({OQL.get_objects(GVdbids)})'
        request_name = f'Find targets with genetic variants linked to {disease_names}'

        # avoid adding GeneticChange to self.Graph
        new_session = self._clone(to_retrieve=NO_REL_PROPERTIES)
        gv2target = new_session.process_oql(oql_query,request_name)
        self.gv2diseases = self.gv2diseases.compose(gv2target)

        target_withGVs = gv2target.psobjs_with(only_with_values=PROTEIN_TYPES)
        gv2dis_refs = self.gv2diseases.load_references()
        print('Found %d targets with %d GeneticVariants for %s supported by %d references' % 
            (len(target_withGVs), len(self.GVs), disease_names, len(gv2dis_refs)))
        
        return target_withGVs


    def targets_from_cache(self):
        oql = f'SELECT Relation WHERE NeighborOf (SELECT Entity WHERE objectType = ({PROTEIN_TYPES})) AND NeighborOf (SELECT Entity WHERE objectType = Disease)'
        name_dict = {'Name':self.params['disease']}
        rank4simplifying = ['Regulation','Biomarker','StateChange','QuantitativeChange']
        rel_props = self._what2retrieve(REFERENCE_IDENTIFIERS)
        prot2disease_graph = self.load_cache('prot2disease',[oql],['Class'],rel_props,name_dict,rank4simplifying)
        self.Graph = prot2disease_graph.compose(self.Graph)
        return prot2disease_graph


    def targets_from_db(self):
        req_name = f'Find targets linked to {self._disease2str()}'
        select_targets = f'SELECT Entity WHERE objectType = ({self.__target_types_str()})'
        OQLquery = f'SELECT Relation WHERE NeighborOf({select_targets}) AND NeighborOf ({self.find_disease_oql})'      
        target_disease_graph = self.process_oql(OQLquery,req_name)
        return target_disease_graph
        

    def find_targets(self):
        target_disease_graph = self.targets_from_cache() if self.use_cache else self.targets_from_db()
        self.__targets__ = set(target_disease_graph.psobjs_with(only_with_values=self.params['target_types']))
        target_refcount = target_disease_graph.load_references()
        print('Found %d targets linked to %s supported by %d references in database' 
                    % (len(self.__targets__),self._disease2str(), len(target_refcount)))
        
        target_withGVs = self.__GVtargets()
        self.__targets__.update(target_withGVs)
        self.targets4strictmode.update(target_withGVs)
        
       # adding targets from disease model pathways
        disease_pathways = self.load_pathways()
        disease_model_components = set()
        if disease_pathways:
            all_pathways = [p for p in disease_pathways.values()]
            [disease_model_components.update(p.get_members(self.params['target_types'])) for p in all_pathways]
            self.__targets__.update(disease_model_components)
            before_add = self.Graph.number_of_nodes()
            self.Graph.add_psobjs(list(disease_model_components))
            print('Added %d targets from disease model' % (self.Graph.number_of_nodes()-before_add))

        if self.params['strict_mode']:
            regulators2disease_graph = target_disease_graph.subgraph_by_relprops(['Regulation'])
            disease_regulators = regulators2disease_graph.psobjs_with(only_with_values=self.params['target_types'])
            self.targets4strictmode.update(disease_regulators)
            self.targets4strictmode.update(disease_model_components)


    def find_symptoms(self):
        if self.params['symptoms']:
            request_name = '\nLoading symptoms listed in script parameters'
            OQLquery = OQL.get_childs(self.params['symptoms'],['Name'],include_parents=True)
            symptom_graph = self.process_oql(OQLquery,request_name)
            self.input_symptoms = symptom_graph._psobjs_with('Disease','ObjTypeName')
            print('Found %d symptoms in ontology under %d symptoms in parameters' % 
                        (len(self.input_symptoms), len(self.params['symptoms'])))


    def find_clinpar(self):
        if self.params['clinical_parameters']:
            request_name = '\nLoading Clinical Parameters entities listed in script parameters'
            OQLquery = OQL.get_childs(self.params['clinical_parameters'],['Name'],include_parents=True)
            clinpar_graph = self.process_oql(OQLquery,request_name)
            self.input_clinpars = set(clinpar_graph._psobjs_with('ClinicalParameter','ObjTypeName'))
            print('Found %d clinical parameters in ontology under %d "clinical parameters" in parameters' % 
                        (len(self.input_clinpars), len(self.params['clinical_parameters'])))


    def find_cellproc(self):
        if self.params['processes']:
            request_name = '\nLoading Cell Process entities listed in script parameters'
            OQLquery = OQL.get_childs(self.params['processes'],['Name'],include_parents=True)
            cellproc_graph = self.process_oql(OQLquery,request_name)
            self.input_cellprocs = set(cellproc_graph._psobjs_with('CellProcess','ObjTypeName'))
            print('Found %d cell processes in ontology under %d "processes" in parameters' % 
                        (len(self.input_cellprocs), len(self.params['processes'])))


    def find_drugs(self):      
        request_name = 'Find compounds inhibiting {}'.format(self._disease2str())
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = negative \
            AND NeighborOf (SELECT Entity WHERE objectType = SmallMol) AND NeighborOf ({})'
        OQLquery = OQLquery.format(self.find_disease_oql)
        drugs_graph = self.process_oql(OQLquery,request_name)
        self.drugs_linked2disease = set(drugs_graph._psobjs_with('SmallMol','ObjTypeName'))
        drugs_references = drugs_graph.load_references()
        print('Found %d compounds inhibiting %s supported by %d references' % 
                    (len(self.drugs_linked2disease),self._disease2str(), len(drugs_references)))

        request_name = 'Find clinical trials for {}'.format(self._disease2str())
        OQLquery = 'SELECT Relation WHERE objectType = ClinicalTrial \
            AND NeighborOf (SELECT Entity WHERE objectType = SmallMol) AND NeighborOf ({})'
        OQLquery = OQLquery.format(self.find_disease_oql)
        ct_graph = self.process_oql(OQLquery,request_name)
        ct_drugs = ct_graph._psobjs_with('SmallMol','ObjTypeName')
        self.drugs_linked2disease.update(ct_drugs)
        ct_refs = ct_graph.load_references()
        print('Found %d compounds on %d clinical trilas for %s' % (len(ct_drugs),len(ct_refs),self._disease2str()))
       

    def find_inducers(self):
        REQUEST_NAME = 'Find compounds inducing/exacerbating {}'.format(self._disease2str())
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = positive \
            AND NeighborOf (SELECT Entity WHERE objectType = SmallMol) AND NeighborOf ({})'
        OQLquery = OQLquery.format(self.find_disease_oql)
        induction_graph = self.process_oql(OQLquery,REQUEST_NAME)
        self.disease_inducers = induction_graph._psobjs_with('SmallMol','ObjTypeName')
        inducers_refs = induction_graph.load_references()
        print('Found %d compounds inducing %s supported by %d references' % 
            (len(self.disease_inducers),self._disease2str(),len(inducers_refs)))

     
    def load_target_partners(self):
        receptor_dbids = self.Graph.dbids4nodes(['Receptor'],['Class'])
        ref_cutoff = 5
        if receptor_dbids:
            find_receptors = OQL.get_objects(receptor_dbids)
            ligands_oql = f'SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive AND RelationNumberOfReferences > {ref_cutoff} AND \
                NeighborOf downstream (SELECT Entity WHERE Class = Ligand) AND NeighborOf upstream ({find_receptors})'
            #opening new APIsession to avoid adding result graph to self.Graph
            new_session = self._clone(to_retrieve=NO_REL_PROPERTIES)
            request_name = f'Find ligands for {str(len(receptor_dbids))} receptors'
            self.partner2targets = new_session.process_oql(ligands_oql,request_name)
            new_ligands = self.partner2targets._psobjs_with('Ligand','Class')
            print('Found %d additional ligands for %d receptors' % (len(new_ligands),len(receptor_dbids)))
        else:
            find_receptors = ''
            print('No Receptor Class entities were found among %s targets' % self.params['target_types'])

        ligand_dbids = self.Graph.dbids4nodes(['Ligand'],['Class'])
        partners_dbids = ResnetGraph.dbids(self.partner2targets._get_nodes())
        orphan_ligands_dbids = set(ligand_dbids).difference(partners_dbids)
        if orphan_ligands_dbids:
            find_orphan_ligands = OQL.get_objects(list(orphan_ligands_dbids))
            receptor_oql = f'SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive AND RelationNumberOfReferences > {ref_cutoff} \
                AND NeighborOf upstream (SELECT Entity WHERE Class = Receptor) AND NeighborOf downstream ({find_orphan_ligands})'
            request_name = f'Find receptors for {str(len(orphan_ligands_dbids))} orphan ligands'
            #opening new APIsession to avoid adding result graph to self.Graph
            new_session = self._clone(to_retrieve=NO_REL_PROPERTIES)
            orphan_ligand_partners = new_session.process_oql(receptor_oql,request_name)
            self.partner2targets.add_graph(orphan_ligand_partners)
            new_receptors = orphan_ligand_partners._psobjs_with('Receptor','Class')
            print('Found %d receptors for %d orphan ligands' % (len(new_receptors),len(orphan_ligands_dbids)))
        
        if find_receptors:
            find_metabolites = OQL.get_childs(['mammal endogenous compounds and their derivatives'],['Name'])
            metabolite_ligands_oql = f'SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive \
                AND NeighborOf downstream ({find_metabolites}) AND NeighborOf upstream ({find_receptors}) AND RelationNumberOfReferences > 25'
            new_session = self._clone(to_retrieve=NO_REL_PROPERTIES)
            metabolite2receptors_graph = new_session.process_oql(metabolite_ligands_oql,'Find metabolite ligands')
            self.partner2targets.add_graph(metabolite2receptors_graph)
            metabolite_ligands = metabolite2receptors_graph._psobjs_with('SmallMol','ObjTypeName')
            print('Found %d metabolite ligands for %d receptors' % (len(metabolite_ligands),len(receptor_dbids)))
        
        find_metabolite_products = 'SELECT Relation WHERE objectType = (MolTransport, MolSynthesis,ChemicalReaction) \
AND Effect = positive AND NeighborOf downstream (SELECT Entity WHERE id =({ids})) AND \
NeighborOf upstream (SELECT Entity WHERE objectType = SmallMol) AND RelationNumberOfReferences > 50'
        req_name = f'Finding downstream metabolite products of targets for {self._disease2str()}'
        new_session = self._clone(to_retrieve=NO_REL_PROPERTIES)

        metbolite_products_graph = new_session.iterate_oql(find_metabolite_products,self._targets_dbids(),request_name=req_name)
        self.partner2targets.add_graph(metbolite_products_graph)
        metabolite_products = metbolite_products_graph._psobjs_with('SmallMol','ObjTypeName')
        print('Found %d metabolites produced by %d targets of %s' % (len(metabolite_products),len(self.__targets__),self._disease2str()))

        # removing promiscous partners
        self.partner2targets.remove_nodes_by_targets(max_target=10,having_values=['SmallMol'])
        self.partner2targets.remove_nodes_by_targets(max_target=10,only_with_prop=['Class'], having_values=['Ligand'])
        self.partner2targets.remove_nodes_by_regulators(max_regulator=10,having_values=['SmallMol'])


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
        fc.entProps = ['Name','CellType','Tissue','Organ','Organ System']

        disease_pathways = list()
        for folder_name in self.params['pathway_folders']:
            ps_pathways = fc.folder2pspathways(folder_name,with_layout=False)
            disease_pathways += ps_pathways

        if self.params['pathways']:
            filtered_pathway_list = list()
            for ps_pathway in disease_pathways:
                if ps_pathway.name() in self.params['pathways']:
                    filtered_pathway_list.append(ps_pathway)
            disease_pathways = filtered_pathway_list

        print('Found %d curated pathways for %s:' %(len(disease_pathways), self._disease2str()))
        [print(p.name()+'\n') for p in disease_pathways]

        #self.disease_pathways = dict()
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
        return self.disease_pathways


    def init_semantic_search(self):
        print('\n\nInitializing semantic search')
        if not self.__targets__:
            print ('No targets found for %s' % self._disease2str())
            print('Consider setting "strict_mode" to False')
            return False

        targets4ranking = list(self._targets())
        self.RefCountPandas = self.load_df(targets4ranking,max_children_count=11)
        print('Will score %d targets linked to %s' % (len(self.RefCountPandas),self._disease2str()))
        self.entProps = ['Name'] # Class is no longer needed
        return True


    def score_GVs(self):
        print('\n\nScoring targets by number of semantic references linking their Genetic Variants to %s' % self._disease2str())
        if not self.GVs:
            print('%s has no known Genetic Variants' % self._disease2str())
            return

        self.__colnameGV__ = 'Genetic Variants for '+self._disease2str()
        target_gvlinkcounter = 0
        for i in self.RefCountPandas.index:
            target_dbids = list(self.RefCountPandas.at[i,self.__temp_id_col__])
            row_targets = self.Graph.psobj_with_dbids(target_dbids)
            targetGVs = self.gv2diseases.get_neighbors(row_targets, allowed_neigbors=self.GVs)
                
            GVscore = 0
            if targetGVs:
                GVnames = set([n.name() for n in targetGVs])
                self.RefCountPandas.at[i,self.__colnameGV__] = ';'.join(GVnames)
                target_gvlinkcounter += 1
                gv_disease_subgraph = self.gv2diseases.get_subgraph(targetGVs,self.input_diseases)
                GVscore = len(gv_disease_subgraph.load_references())
            
            self.RefCountPandas.at[i,self._col_name_prefix+'Genetic Variants'] = GVscore

        print('Found %d targets linked to %d GVs' % (target_gvlinkcounter, len(self.GVs)))


    def __vote4effect(self,targets:list):
        between_graph = self.Graph.get_subgraph(targets, self.input_diseases)
        positive_refs, negative_refs = between_graph._effect_counts__()

        if len(positive_refs) > len(negative_refs):
            return 'activated'
        elif len(negative_refs) > len(positive_refs):
            return 'repressed'
        else:
            target_partners = set(self.partner2targets.get_neighbors(targets))
            partners2disease_graph = self.Graph.get_subgraph(target_partners, self.input_diseases)
            positive_refs, negative_refs = between_graph._effect_counts__()
            if len(positive_refs) > len(negative_refs):
                return 'activated'
            elif len(negative_refs) > len(positive_refs):
                return 'repressed'
            else:
                # attempting to deduce effect from GeneticChange
                genetic_change = partners2disease_graph.psrels_with(['GeneticChange'])
                if genetic_change:
                    return 'repressed'
                else:
                    return 'unknown'


    def set_target_disease_state(self):
        print('\n\nCalculating targets state (activated/repressed) in %s' % self._disease2str())
                       
        for i in self.RefCountPandas.index:
            target_dbids = list(self.RefCountPandas.at[i,self.__temp_id_col__])
            targets = self.Graph.psobj_with_dbids(target_dbids)
            self.RefCountPandas.at[i,'State in Disease'] = self.__vote4effect(targets)

        # determining disease state of the targets from disease model pathways
        # by finding their regulatory effect on targets linked to disease
        model = ResnetGraph()
        for pathway in self.disease_pathways.values():
            model = model.compose(pathway.graph)

        model_dbid2uid = model.dbid2uid()
        def __disease_state_from_model():
            activ_targets_pd = self.RefCountPandas.loc[(self.RefCountPandas['State in Disease'] == 'activated')]
            inhib_targets_pd = self.RefCountPandas.loc[(self.RefCountPandas['State in Disease'] == 'repressed')]
            activated_targets_dbid2uid = model.dbid2uid(self._all_dbids(activ_targets_pd))
            inhibited_targets_dbid2uid = model.dbid2uid(self._all_dbids(inhib_targets_pd))

            new_regulators_count = 0
            for idx in self.RefCountPandas.index:
                if self.RefCountPandas.at[idx,'State in Disease'] == 'unknown':
                    net_effect = 0
                    for target_dbid in self.RefCountPandas.at[idx,self.__temp_id_col__]:
                        try:
                            target_uid = model_dbid2uid[target_dbid]
                            activates_activated_targets,inhibits_activated_targets = model.net_regulator_effect(target_uid,list(activated_targets_dbid2uid.values()))
                            activates_inhibited_targets,inhibits_inhibited_targets = model.net_regulator_effect(target_uid,list(inhibited_targets_dbid2uid.values()))
                            net_effect += len(activates_activated_targets)+len(inhibits_inhibited_targets)-len(inhibits_activated_targets)-len(activates_inhibited_targets)
                        except KeyError:
                            continue

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
        self.RefCountPandas['target_dbids'] = self.RefCountPandas[self.__temp_id_col__]
        for i in self.RefCountPandas.index:
            target_dbids = list(self.RefCountPandas.at[i,'target_dbids'])
            targets_in_row = self.Graph.psobj_with_dbids(target_dbids)
            # all_targets may include not linked to input_diseases of targets linked that are children of linked targets
            target_partners = set(self.partner2targets.get_neighbors(targets_in_row))
            if target_partners:
                target_partners_dbids = ResnetGraph.dbids(target_partners)
                self.RefCountPandas.at[i,self.__temp_id_col__] = tuple(target_partners_dbids)
            else:
                self.RefCountPandas.at[i,self.__temp_id_col__] = tuple([0]) # fake dbid
            partner_names = [n.name() for n in target_partners]
            self.RefCountPandas.at[i,'Target partners'] = ';'.join(partner_names)
        
        self.set_how2connect([],[],'')
        colname = 'target partners'
        linked_row_count = self.link2RefCountPandas(colname,self.input_diseases)
        print('%d targets have partners linked to %s' % (linked_row_count,self._disease2str()))

        self.RefCountPandas.drop(columns=self.__temp_id_col__,inplace=True)
        self.RefCountPandas.rename(columns={'target_dbids':self.__temp_id_col__}, inplace=True)


    def score_regulators(self):
        print('\n\nScoring regulators by distance to components of disease pathways',flush=True)
        print('Retrieving regulatory network between targets and components of disease pathway')
        regulation_graph = ResnetGraph()
        [regulation_graph.add_graph(p.graph) for p in self.disease_pathways.values()]
        disease_pathway_component_dbids = set(regulation_graph.dbids4nodes())

        connected_nodes_dbids = disease_pathway_component_dbids
        unconnected_nodes_dbids = self.all_entity_dbids
        for step in range (0,5):
            new_session = self._clone(to_retrieve=NO_REL_PROPERTIES)
            graph_expansion = new_session.connect_nodes(unconnected_nodes_dbids, connected_nodes_dbids,
                                                    PHYSICAL_INTERACTIONS, in_direction='>')
            regulation_graph.add_graph(graph_expansion)

            regulation_graph_dbids = regulation_graph.dbids4nodes()
            connected_nodes_dbids = set(regulation_graph_dbids).difference(disease_pathway_component_dbids)
            connected_nodes_dbids = connected_nodes_dbids.intersection(unconnected_nodes_dbids)
            # only nodes connected at the previous cycle need to be expanded at the next cycle
            if not connected_nodes_dbids:
                break
            unconnected_nodes_dbids = set(unconnected_nodes_dbids).difference(regulation_graph_dbids)

        for pathway in self.disease_pathways.values():
            closeness_dic = {i:c[0] for i,c in pathway.graph.nodes(data=CLOSENESS)}
            regulation_graph.rank_regulators(closeness_dic,PATHWAY_REGULATOR_SCORE)

        dbid2uid = regulation_graph.dbid2uid()
        for i in self.RefCountPandas.index:
            target_dbids = list(self.RefCountPandas.at[i,self.__temp_id_col__])
            target_regulator_score = 0.0
            for target_dbid in target_dbids:
                try:
                    target_uid = dbid2uid[target_dbid]
                    regulator_score = regulation_graph.nodes[target_uid][PATHWAY_REGULATOR_SCORE]
                    target_regulator_score += regulator_score
                except KeyError:
                    continue
            self.RefCountPandas.at[i,PATHWAY_REGULATOR_SCORE] = target_regulator_score


    def score_target_semantics(self):
        t_n = self._disease2str()
        self.score_GVs()

        colname = 'Regulate '+ t_n
        self.set_how2connect(['Regulation'],[],'')
        linked_row_count = self.link2RefCountPandas(colname,self.input_diseases)
        print('%d targets regulating %s' % (linked_row_count,t_n))

        self.set_target_disease_state()

        colname = 'Genetically linked to '+ t_n
        self.set_how2connect(['GeneticChange'],[],'')
        linked_row_count = self.link2RefCountPandas(colname,self.input_diseases)
        print('%d targets genetically linked to %s' % (linked_row_count,t_n))

        colname = 'Target is Biomarker in '+t_n
        self.set_how2connect(['Biomarker'],[],'')
        linked_row_count = self.link2RefCountPandas(colname,self.input_diseases)
        print('Linked %d indications where %s is biomarker' % (linked_row_count,t_n))

        colname = 'Quantitatively changed in '+ t_n
        self.set_how2connect(['QuantitativeChange'],[],'')
        linked_row_count = self.link2RefCountPandas(colname,self.input_diseases)
        print('%d targets quantitatively changed in %s' % (linked_row_count,t_n))

        self.score_partners()

        if self.input_symptoms:
            colname = 'symptoms for '+ t_n
            reltypes2connect = ['Regulation','QuantitativeChange','StateChange','Biomarker']
            self.set_how2connect(reltypes2connect,[],'')
            linked_row_count = self.link2RefCountPandas(colname,self.input_symptoms,REFERENCE_IDENTIFIERS)
            print('%d targets linked to symptoms for in %s' % (linked_row_count,t_n))

        if self.input_clinpars:
            colname = 'Clinical parameters for '+ t_n
            self.set_how2connect(['Regulation'],[],'')
            linked_row_count = self.link2RefCountPandas(colname,self.input_clinpars,REFERENCE_IDENTIFIERS)
            print('%d targets linked to clinical parameters for in %s' % (linked_row_count,t_n))

        if self.disease_pathways:
            self.score_regulators()

        if self.input_cellprocs:
            colname = 'Cell processes affected by '+ t_n
            self.set_how2connect(['Regulation'],[],'')
            linked_row_count = self.link2RefCountPandas(colname,self.input_cellprocs,REFERENCE_IDENTIFIERS)
            print('%d targets linked to cell processes for %s' % (linked_row_count,t_n))

        ################## SPLITTING RefCountPandas to score agonists and antagonists differently #####################
        activ_targets_pd = self.RefCountPandas.loc[(self.RefCountPandas['State in Disease'] == 'activated')]
        inhib_targets_pd = self.RefCountPandas.loc[(self.RefCountPandas['State in Disease'] == 'repressed')]
        unk_targets_pd = self.RefCountPandas.loc[(self.RefCountPandas['State in Disease'] == 'unknown')]

        if hasattr(self, 'drugs_linked2disease'):
            all_compounds = set()
            # drug-target relations are not added to self.Graph to exclude them from bibliography
            colname = 'drugs for '+t_n
            self.set_how2connect(['DirectRegulation'],['negative'],'',['Binding','Expression','Regulation'])
            print ('\n\nLinking with Effect negative drugs for "%s" to targets activated in "%s"' % (t_n,t_n))
            was_linked,linked_entities,activ_targets_pd = self.link2concept(
                colname,self.drugs_linked2disease,activ_targets_pd,REFERENCE_IDENTIFIERS)
            print('Linked %d targets to drugs for "%s"' % (was_linked,t_n))
            

            self.set_how2connect(['DirectRegulation'],['positive'],'',['Binding','Expression','Regulation'])
            print ('\n\nLinking with Effect positive drugs for "%s" to targets inhibiting "%s"' % (t_n,t_n))
            was_linked,linked_entities,inhib_targets_pd = self.link2concept(
                colname,self.drugs_linked2disease,inhib_targets_pd,REFERENCE_IDENTIFIERS)
            print('Linked %d targets to drugs for "%s"' % (was_linked,t_n))
            all_compounds.update(self.drugs_linked2disease)

        if hasattr(self, 'disease_inducers'):
            colname = 'inducers of '+t_n
            self.set_how2connect(['DirectRegulation'],['positive'],'',['Binding','Expression','Regulation'])
            print ('\n\nLinking with Effect positive compounds inducing "%s" to targets activated in "%s"' % (t_n,t_n))
            was_linked,linked_entities,activ_targets_pd = self.link2concept(
                colname,self.disease_inducers,activ_targets_pd,REFERENCE_IDENTIFIERS)
            print('Linked %d targets to inducers of "%s"' % (was_linked,t_n))

            self.set_how2connect(['DirectRegulation'],['negative'],'',['Binding','Expression','Regulation'])
            print ('\n\nLinking with effect negative compounds inducing "%s" to targets inhibiting "%s"' % (t_n,t_n))
            was_linked,linked_entities,inhib_targets_pd = self.link2concept(
                colname,self.disease_inducers,inhib_targets_pd,REFERENCE_IDENTIFIERS)
            print('Linked %d targets to inducers of "%s"' % (was_linked,t_n))
            all_compounds.update(self.disease_inducers)

        if all_compounds:
            colname = 'compounds modulating '+t_n
            print ('\n\nLinking compounds modulating "%s" to targets with unknown state in "%s"' % (t_n,t_n))
            self.set_how2connect(['DirectRegulation'],[],'',['Binding','Expression','Regulation'])
            was_linked,linked_entities,unk_targets_pd = self.link2concept(
                colname,all_compounds,unk_targets_pd,REFERENCE_IDENTIFIERS)
            print('Linked %d targets to inducers for %s' % (was_linked,t_n))

        uptarget_df = self.make_count_df(df(activ_targets_pd),DF4ANTAGONISTS)
        self.add2raw(uptarget_df)
        downtarget_df = self.make_count_df(df(inhib_targets_pd),DF4AGONISTS)
        self.add2raw(downtarget_df)
        unktarget_df = self.make_count_df(df(unk_targets_pd),DF4UNKNOWN)
        self.add2raw(unktarget_df)
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


    def add_ps_bibliography(self):
        input_disease_graph = self.Graph.neighborhood(self.input_diseases)
        super().add_ps_bibliography(from_graph=input_disease_graph)


    def add_target_annotation(self,_2column='Class'):
        self.add_entity_annotation(_2column,self.report_pandas[ANTAGONIST_TARGETS_WS])
        self.add_entity_annotation(_2column,self.report_pandas[AGONIST_TARGETS_WS])
        self.add_entity_annotation(_2column,self.report_pandas[UNKNOWN_TARGETS_WS])


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
        self.add_target_annotation()
    
        print('Target ranking for %s was done in %s' % 
        (self._disease2str(), self.execution_time(start_time)))
        

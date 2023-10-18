from .ResnetGraph import PROTEIN_TYPES,EFFECT,PHYSICAL_INTERACTIONS,ResnetGraph,CLOSENESS,CHILDS
from .NetworkxObjects import STATE, ACTIVATED,REPRESSED,UNKNOWN_STATE
from .SemanticSearch import SemanticSearch,len,OQL
from .ResnetAPISession import NO_REL_PROPERTIES,ONLY_REL_PROPERTIES,REFERENCE_IDENTIFIERS,BIBLIO_PROPERTIES,DBID
from .FolderContent import FolderContent
from concurrent.futures import ThreadPoolExecutor
from ..pandas.panda_tricks import df,NaN
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
                'ent_props' : ['Name', 'Class'], #'Class' is for target partners retrieval
                'rel_props' : [EFFECT],
                'strict_mode' : True,
                'data_dir' : '',
                'add_bibliography' : True,
                'strict_mode' : False,
                'target_types' : ['Protein'],
                'pathway_folders' : ['Hypertension Pulmonary'],
                'pathways' : []
            }
        
        entprops = kwargs.pop('ent_props',[])
        relprops = kwargs.pop('rel_props',[])

        my_kwargs.update(kwargs)
        super().__init__(APIconfig,**my_kwargs)
        self.add_ent_props(entprops)
        self.add_rel_props(relprops)

        self.columns2drop += [self.__resnet_name__,self.__mapped_by__,'State in Disease']
        self.max_threads4ontology = 50

        self.input_diseases = list()
        self.input_symptoms = list()
        self.input_clinpars = list()
        self.input_cellprocs = list()

        self.__targets__ = set()
        self.GVs = list()
        self.targets4strictmode = set()
        self.ct_drugs = set()
        self.drugs_linked2disease = set()
        self.disease_inducers = list()
        
        self.gv2diseases = ResnetGraph()
        self.partner2targets = ResnetGraph()
        self.disease_pathways = dict() # {cell_type:PSPathway}

        # self.targets4strictmode has PSObjects for targets directly Regulating disease or with GVs linked to disease or from disease pathways
        # if strict_mode algorithm only ranks targets suggetsed in the literarure without predicting new targets
        # in not strict_mode algorithm also predicts and ranks additional targets linked to disease by QuantitativeChange or GeneticChange     
       
        # TBD cache networks:
        self.prot2Dis = ResnetGraph() # includes GVs
        self.prot2ClinPar = ResnetGraph() # includes GVs
        self.protein2CelProc = ResnetGraph() # protein-CellProcess regulation network (Regulation) to find partners
        self.ppmet = ResnetGraph() # protein-metabolite regulation network (MolTransport,MolSynthesis,ChemicalReaction) for self.partner2targets
        self.drug_effects = ResnetGraph() # drug-disease regulation network (Regulation) disease_inducers, drugs_link2disease
        self.drugs2targets = ResnetGraph() # for disease_inducers, drugs_link2disease
    

    def clear(self):
        super().clear()
        self.gv2diseases.clear_resnetgraph()
        self.partner2targets.clear_resnetgraph()
        self.disease_pathways.clear()
        

    def set_input(self):
        my_session = self._clone_session()
        my_session.add_ent_props(['Connectivity']) # Connectivity is needed for make_disease_df()
        ontology_graph = my_session.child_graph(self.params['disease'],['Name','Alias'])
        if isinstance(ontology_graph, ResnetGraph):
            self.input_diseases = ontology_graph._psobjs_with('Disease','ObjTypeName') 
            # filter by Disease is necessary because some proteins in monogenic disorders may have the same name as disease
            input_diseases_dbids = ResnetGraph.dbids(self.input_diseases)
            self.find_disease_oql = OQL.get_objects(input_diseases_dbids)
        return


    def input_names(self):
        if self.input_diseases:
            prop2values = {'Name':self.params['disease'], 'Alias':self.params['disease']}
            return list(set([x['Name'][0] for x in self.input_diseases if x.has_value_in(prop2values)]))
        else:
            try:# case of SNEA drug_df
                return [self.params['sample']+' in '+ self.params['experiment']]
            except KeyError:
                return ['Drugs']
            

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
        return ResnetGraph.dbids(list(self._targets()))

    def __target_types_str(self):
        tar_types_str = ','.join(self.params['target_types'])
        return tar_types_str
        
    def __GVtargets(self):
        disease_names = self._disease2str()

        oql_query = f'SELECT Relation WHERE objectType = FunctionalAssociation \
            AND NeighborOf ({self.find_disease_oql}) AND NeighborOf (SELECT Entity WHERE objectType = GeneticVariant)'
        #oql_query = oql_query.format()
        request_name = f'Select GeneticVariants for {disease_names}'
        self.gv2diseases = self.process_oql(oql_query,request_name)
        if isinstance(self.gv2diseases,ResnetGraph):
            self.GVs = self.gv2diseases.psobjs_with(only_with_values=['GeneticVariant'])
        
            GVdbids = ResnetGraph.dbids(self.GVs)
            target_withGVs = list()

            if GVdbids:
                oql_query = f'SELECT Relation WHERE objectType = GeneticChange AND \
                    NeighborOf (SELECT Entity WHERE objectType = ({self.__target_types_str()})) AND NeighborOf ({OQL.get_objects(GVdbids)})'
                request_name = f'Find targets with genetic variants linked to {disease_names}'

                # avoid adding GeneticChange to self.Graph
                new_session = self._clone(to_retrieve=NO_REL_PROPERTIES)
                gv2target = new_session.process_oql(oql_query,request_name)
                if isinstance(gv2target, ResnetGraph):
                    self.gv2diseases = self.gv2diseases.compose(gv2target)
                    target_withGVs = gv2target.psobjs_with(only_with_values=PROTEIN_TYPES)
                    gv2dis_refs = self.gv2diseases.load_references()
                    print('Found %d targets with %d GeneticVariants for %s supported by %d references' % 
                        (len(target_withGVs), len(self.GVs), disease_names, len(gv2dis_refs)))
        
                    return target_withGVs
        return ResnetGraph()

 
    def targets_from_db(self):
        req_name = f'Find targets linked to {self._disease2str()}'
        select_targets = f'SELECT Entity WHERE objectType = ({self.__target_types_str()})'
        OQLquery = f'SELECT Relation WHERE NeighborOf ({select_targets}) AND NeighborOf ({self.find_disease_oql})'      
        target_disease_graph = self.process_oql(OQLquery,req_name)
        return target_disease_graph
        

    def find_targets(self):
        target_disease_graph = self.targets_from_db()
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
            self.Graph.add_psobjs(set(disease_model_components))
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
            self.input_clinpars = clinpar_graph._psobjs_with('ClinicalParameter','ObjTypeName')
            print('Found %d clinical parameters in ontology under %d "clinical parameters" in parameters' % 
                        (len(self.input_clinpars), len(self.params['clinical_parameters'])))


    def find_cellproc(self):
        if self.params['processes']:
            request_name = '\nLoading Cell Process entities listed in script parameters'
            OQLquery = OQL.get_childs(self.params['processes'],['Name'],include_parents=True)
            cellproc_graph = self.process_oql(OQLquery,request_name)
            self.input_cellprocs = cellproc_graph._psobjs_with('CellProcess','ObjTypeName')
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
        self.ct_drugs = set(ct_graph._psobjs_with('SmallMol','ObjTypeName'))
        self.drugs_linked2disease.update(self.ct_drugs)
        ct_refs = ct_graph.load_references()
        print('Found %d compounds on %d clinical trilas for %s' % (len(self.ct_drugs),len(ct_refs),self._disease2str()))
       

    def find_inducers(self):
        REQUEST_NAME = 'Find compounds inducing/exacerbating {}'.format(self._disease2str())
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = positive \
            AND NeighborOf (SELECT Entity WHERE objectType = SmallMol) AND NeighborOf ({})'
        OQLquery = OQLquery.format(self.find_disease_oql)
        induction_graph = self.process_oql(OQLquery,REQUEST_NAME)
        self.disease_inducers = induction_graph._psobjs_with('SmallMol','ObjTypeName')
        count = len(self.disease_inducers)
        self.disease_inducers =[x for x in self.disease_inducers if x not in self.ct_drugs]
        remove_count = count - len(self.disease_inducers)
        # to remove indications reported as toxicities in clinical trials
        print(f'{remove_count} {self._disease2str()} inducers were removed because they are linked by clinical trial')
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
            p2t = new_session.process_oql(ligands_oql,request_name)
            if isinstance(p2t,ResnetGraph):
                self.partner2targets = p2t
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
            if isinstance(orphan_ligand_partners,ResnetGraph):
                self.partner2targets.add_graph(orphan_ligand_partners)
                new_receptors = orphan_ligand_partners._psobjs_with('Receptor','Class')
                print('Found %d receptors for %d orphan ligands' % (len(new_receptors),len(orphan_ligands_dbids)))
        
        if find_receptors:
            find_metabolites = OQL.get_childs(['mammal endogenous compounds and their derivatives'],['Name'])
            metabolite_ligands_oql = f'SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive \
                AND NeighborOf downstream ({find_metabolites}) AND NeighborOf upstream ({find_receptors}) AND RelationNumberOfReferences > 25'
            new_session = self._clone(to_retrieve=NO_REL_PROPERTIES)
            metabolite2receptors_graph = new_session.process_oql(metabolite_ligands_oql,'Find metabolite ligands')
            if isinstance(metabolite2receptors_graph, ResnetGraph):
                self.partner2targets.add_graph(metabolite2receptors_graph)
                metabolite_ligands = metabolite2receptors_graph._psobjs_with('SmallMol','ObjTypeName')
                print('Found %d metabolite ligands for %d receptors' % (len(metabolite_ligands),len(receptor_dbids)))
        
        find_metabolite_products = 'SELECT Relation WHERE objectType = (MolTransport, MolSynthesis,ChemicalReaction) \
AND Effect = positive AND NeighborOf downstream (SELECT Entity WHERE id =({ids})) AND \
NeighborOf upstream (SELECT Entity WHERE objectType = SmallMol) AND RelationNumberOfReferences > 50'
        req_name = f'Finding downstream metabolite products of targets for {self._disease2str()}'
        new_session = self._clone(to_retrieve=NO_REL_PROPERTIES)

        metbolite_products_graph = new_session.iterate_oql(find_metabolite_products,set(self._targets_dbids()),request_name=req_name)
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
        self.RefCountPandas = self.load_df(targets4ranking,max_child_count=11,max_threads=10)
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
            row_targets = self.Graph.psobj_with_dbids(set(target_dbids))
            targetGVs = self.gv2diseases.get_neighbors(set(row_targets), allowed_neigbors=self.GVs)
                
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
            return ACTIVATED
        elif len(negative_refs) > len(positive_refs):
            return REPRESSED
        else:
            target_partners = set(self.partner2targets.get_neighbors(set(targets)))
            partners2disease_graph = self.Graph.get_subgraph(list(target_partners), self.input_diseases)
            positive_refs, negative_refs = between_graph._effect_counts__()
            if len(positive_refs) > len(negative_refs):
                return ACTIVATED
            elif len(negative_refs) > len(positive_refs):
                return REPRESSED
            else:
                # attempting to deduce effect from GeneticChange
                genetic_change = partners2disease_graph.psrels_with(['GeneticChange'])
                if genetic_change:
                    return REPRESSED
                else:
                    return UNKNOWN_STATE


    def make_disease_network(self):
        print(f'Creating physical interaction network between targets of {self.input_names()} \
              to calculate target closeness for "score_regulators" function')
        my_session = self._clone_session(to_retrieve=NO_REL_PROPERTIES)
        targets_dbid = ResnetGraph.dbids(list(self.__targets__))
        disease_network = my_session.get_ppi(set(targets_dbid))
        disease_network.name = f'{self._disease2str()} PPPI network'
        return disease_network.make_simple(['DirectRegulation','ProtModification','Binding']) 


    def set_target_disease_state(self):
        print('\n\nCalculating targets state (activated/repressed) in %s' % self._disease2str())
        targets = set()
        for i in self.RefCountPandas.index:
            target_dbids = list(self.RefCountPandas.at[i,self.__temp_id_col__])
            row_targets = self.Graph.psobj_with_dbids(set(target_dbids))
            state = self.__vote4effect(row_targets)
            self.RefCountPandas.at[i,'State in Disease'] = int(state)
            if state != UNKNOWN_STATE:
                [t.set_state(state) for t in row_targets]
                targets.update(row_targets)

        self.disease_model = ResnetGraph()
        for pathway in self.disease_pathways.values():
            self.disease_model = self.disease_model.compose(pathway.graph)

        if self.disease_model and 'propogate_target_state_in_model' in self.params.keys():
            uid2state = dict()
            [uid2state.update(self.disease_model.propagate_state(t)) for t in targets if t.uid() in self.disease_model.nodes()]
            for i in self.RefCountPandas.index:
                if self.RefCountPandas.at[i,'State in Disease'] == UNKNOWN_STATE:
                    target_dbids = list(self.RefCountPandas.at[i,self.__temp_id_col__])
                    row_targets = self.Graph.psobj_with_dbids(set(target_dbids))
                    row_state = UNKNOWN_STATE
                    for t in row_targets:
                        try:
                            t_state = uid2state[t.uid()]
                            row_state += t_state
                        except KeyError:
                            continue
                    self.RefCountPandas.at[i,'State in Disease'] = row_state

            # determining disease state of the targets from disease model pathways
            # by finding their regulatory effect on targets linked to disease
            model_dbid2uid = self.disease_model.dbid2uid()
            def __disease_state_from_model():
                activ_targets_pd = self.RefCountPandas.loc[(self.RefCountPandas['State in Disease'] == ACTIVATED)]
                #activated_targets_dbid2uid = model.dbid2uid(self._all_dbids(activ_targets_pd))
                activated_targets = self.disease_model.psobj_with_dbids(self._all_dbids(activ_targets_pd))
                inhib_targets_pd = self.RefCountPandas.loc[(self.RefCountPandas['State in Disease'] == REPRESSED)]
                inhibited_trgts = self.disease_model.psobj_with_dbids(self._all_dbids(inhib_targets_pd))
                #inhibited_trgts_dbid2uid = model.dbid2uid(self._all_dbids(inhib_targets_pd))

                new_regulators_count = 0
                for idx in self.RefCountPandas.index:
                    if self.RefCountPandas.at[idx,'State in Disease'] == UNKNOWN_STATE:
                        net_effect = 0
                        for target_dbid in self.RefCountPandas.at[idx,self.__temp_id_col__]:
                            try:
                                target_uid = model_dbid2uid[target_dbid]
                                target = self.disease_model._get_node(target_uid)
                                activates_activated_targets,inhibits_activated_targets = self.disease_model.net_regulator_effect(target,activated_targets)
                                activates_inhibited_targets,inhibits_inhibited_targets = self.disease_model.net_regulator_effect(target,inhibited_trgts)
                                net_effect += len(activates_activated_targets)+len(inhibits_inhibited_targets)-len(inhibits_activated_targets)-len(activates_inhibited_targets)
                            except KeyError:
                                continue

                        if net_effect > 0:
                            self.RefCountPandas.at[idx,'State in Disease'] = ACTIVATED
                            new_regulators_count += 1
                        elif net_effect < 0:
                            self.RefCountPandas.at[idx,'State in Disease'] = REPRESSED
                            new_regulators_count += 1
                            
                return new_regulators_count
            
            new_regulators_count = __disease_state_from_model()
            while new_regulators_count:
                new_regulators_count =  __disease_state_from_model()
            return
        return

    def score_partners(self):
        print(f'\n\nFinding semantic references supporting links between target partners and {self._disease2str()}')
        my_df = df.copy_df(self.RefCountPandas)
        partners_dbids = list()
        partner_names = list()
        target_names = list()
        for i in my_df.index:
            target_dbids = list(my_df.at[i,self.__temp_id_col__])
            targets_in_row = self.Graph.psobjs_with([DBID],target_dbids)
            # all_targets may include not linked to input_diseases of targets linked that are children of linked targets
            target_partners = set(self.partner2targets.get_neighbors(set(targets_in_row)))
            names = [n.name() for n in target_partners]
            partner_names.append(','.join(names))
            target_names.append(my_df.at[i,'Name'])
            if target_partners:
                target_partners_dbids = ResnetGraph.dbids(list(target_partners))
                partners_dbids.append(tuple(target_partners_dbids))
            else:
                partners_dbids.append(tuple([0])) # adding fake dbid for link2concept to work

        partners_df = df.from_dict({'Name':target_names,self.__temp_id_col__:partners_dbids,'Target partners':partner_names})
        how2connect = self.set_how2connect([],[],'')
        colname = 'target partners'
        linked_row_count,linked_entity_ids,partners_df = self.link2concept(colname,self.input_diseases,partners_df,how2connect)
        print('%d targets have partners linked to %s' % (linked_row_count,self._disease2str()))
        return df.from_pd(partners_df.drop(columns=[self.__temp_id_col__]))


    def score_partnersOLD(self):
        print('\n\nFinding semantic references supporting links between target partners and %s'%
                 self._disease2str())
        temp_targetdbid_col = 'target_dbids'
        self.RefCountPandas[temp_targetdbid_col] = self.RefCountPandas[self.__temp_id_col__]
        for i in self.RefCountPandas.index:
            target_dbids = list(self.RefCountPandas.at[i,temp_targetdbid_col])
            targets_in_row = self.Graph.psobjs_with([DBID],target_dbids)
            # all_targets may include not linked to input_diseases of targets linked that are children of linked targets
            target_partners = set(self.partner2targets.get_neighbors(set(targets_in_row)))
            if target_partners:
                target_partners_dbids = ResnetGraph.dbids(list(target_partners))
                self.RefCountPandas.at[i,self.__temp_id_col__] = tuple(target_partners_dbids)
            else:
                self.RefCountPandas.at[i,self.__temp_id_col__] = tuple([0]) # fake dbid
            partner_names = [n.name() for n in target_partners]
            self.RefCountPandas.at[i,'Target partners'] = ';'.join(partner_names)
        
        how2connect = self.set_how2connect([],[],'')
        colname = 'target partners'
        linked_row_count = self.link2RefCountPandas(colname,self.input_diseases,how2connect)
        print('%d targets have partners linked to %s' % (linked_row_count,self._disease2str()))

        self.RefCountPandas.drop(columns=self.__temp_id_col__,inplace=True)
        self.RefCountPandas.rename(columns={temp_targetdbid_col:self.__temp_id_col__}, inplace=True)


    def score_regulators(self):
        print('\n\nScoring regulators by distance to components of disease pathways',flush=True)
        print('Retrieving regulatory network between targets and components of disease pathway')
        regulation_graph = self.disease_model if self.disease_model else self.make_disease_network()
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
      #  do not multithread.  self.Graph will leak
      # e = ThreadPoolExecutor(max_workers=1, thread_name_prefix='Target partners score')
      #  partner_score_future = e.submit(self.score_partners)
        t_n = self._disease2str()
        self.score_GVs()

        colname = 'Regulate '+ t_n
        how2connect = self.set_how2connect(['Regulation'],[],'',['FunctionalAssociation'])
        linked_row_count = self.link2RefCountPandas(colname,self.input_diseases,how2connect)
        print('%d targets regulating %s' % (linked_row_count,t_n))

        self.set_target_disease_state()

        colname = 'Genetically linked to '+ t_n
        how2connect = self.set_how2connect(['GeneticChange'],[],'',['FunctionalAssociation'])
        linked_row_count = self.link2RefCountPandas(colname,self.input_diseases,how2connect)
        print('%d targets genetically linked to %s' % (linked_row_count,t_n))

        colname = 'Target is Biomarker in '+t_n
        how2connect = self.set_how2connect(['Biomarker'],[],'',['FunctionalAssociation'])
        linked_row_count = self.link2RefCountPandas(colname,self.input_diseases,how2connect)
        print('Linked %d indications where %s is biomarker' % (linked_row_count,t_n))

        colname = 'Quantitatively changed in '+ t_n
        how2connect = self.set_how2connect(['QuantitativeChange'],[],'',['FunctionalAssociation'])
        linked_row_count = self.link2RefCountPandas(colname,self.input_diseases,how2connect)
        print('%d targets quantitatively changed in %s' % (linked_row_count,t_n))

        if self.input_symptoms:
            colname = 'symptoms for '+ t_n
            reltypes2connect = ['Regulation','QuantitativeChange','StateChange','Biomarker']
            how2connect = self.set_how2connect(reltypes2connect,[],'',['FunctionalAssociation'])
            linked_row_count = self.link2RefCountPandas(colname,self.input_symptoms,how2connect,REFERENCE_IDENTIFIERS)
            print('%d targets linked to symptoms for in %s' % (linked_row_count,t_n))

        if self.input_clinpars:
            colname = 'Clinical parameters for '+ t_n
            how2connect = self.set_how2connect(['Regulation'],[],'',['FunctionalAssociation'])
            linked_row_count = self.link2RefCountPandas(colname,self.input_clinpars,how2connect,REFERENCE_IDENTIFIERS)
            print('%d targets linked to clinical parameters for in %s' % (linked_row_count,t_n))

        if self.input_cellprocs:
            colname = 'Cell processes affected by '+ t_n
            how2connect = self.set_how2connect(['Regulation'],[],'',['FunctionalAssociation'])
            linked_row_count = self.link2RefCountPandas(colname,self.input_cellprocs,how2connect,REFERENCE_IDENTIFIERS)
            print('%d targets linked to cell processes for %s' % (linked_row_count,t_n))

        #partners_df = partner_score_future.result()
        partners_df = self.score_partners()
        self.RefCountPandas = self.RefCountPandas.merge_df(partners_df,on='Name')


        ################## SPLITTING RefCountPandas to score agonists and antagonists differently #####################
        activ_targets_df = df.from_pd(self.RefCountPandas.loc[(self.RefCountPandas['State in Disease'] > UNKNOWN_STATE)])
        inhib_targets_df = df.from_pd(self.RefCountPandas.loc[(self.RefCountPandas['State in Disease'] < UNKNOWN_STATE)])
        unk_targets_df = df.from_pd(self.RefCountPandas.loc[(self.RefCountPandas['State in Disease'] == UNKNOWN_STATE)])

        all_compounds = set()
        if hasattr(self, 'drugs_linked2disease'):
            # drug-target relations are not added to self.Graph to exclude them from bibliography
            colname = 'drugs for '+t_n
            how2connect = self.set_how2connect(['DirectRegulation'],['negative'],'',['Binding','Expression','Regulation'])
            print ('\n\nLinking with Effect negative drugs for "%s" to targets activated in "%s"' % (t_n,t_n))
            was_linked,linked_entities,activ_targets_df = self.link2concept(
                colname,list(self.drugs_linked2disease),activ_targets_df,how2connect,REFERENCE_IDENTIFIERS)
            print('Linked %d targets to drugs for "%s"' % (was_linked,t_n))
        
            how2connect = self.set_how2connect(['DirectRegulation'],['positive'],'',['Binding','Expression','Regulation'])
            print ('\n\nLinking with Effect positive drugs for "%s" to targets inhibiting "%s"' % (t_n,t_n))
            was_linked,linked_entities,inhib_targets_df = self.link2concept(
                colname,list(self.drugs_linked2disease),inhib_targets_df,how2connect,REFERENCE_IDENTIFIERS)
            print('Linked %d targets to drugs for "%s"' % (was_linked,t_n))
            all_compounds.update(self.drugs_linked2disease)

        if hasattr(self, 'disease_inducers'):
            colname = 'inducers of '+t_n
            how2connect = self.set_how2connect(['DirectRegulation'],['positive'],'',['Binding','Expression','Regulation'])
            print ('\n\nLinking with Effect positive compounds inducing "%s" to targets activated in "%s"' % (t_n,t_n))
            was_linked,linked_entities,activ_targets_df = self.link2concept(
                colname,self.disease_inducers,activ_targets_df,how2connect,REFERENCE_IDENTIFIERS)
            print('Linked %d targets to inducers of "%s"' % (was_linked,t_n))

            how2connect = self.set_how2connect(['DirectRegulation'],['negative'],'',['Binding','Expression','Regulation'])
            print ('\n\nLinking with effect negative compounds inducing "%s" to targets inhibiting "%s"' % (t_n,t_n))
            was_linked,linked_entities,inhib_targets_df = self.link2concept(
                colname,self.disease_inducers,inhib_targets_df,how2connect,REFERENCE_IDENTIFIERS)
            print('Linked %d targets to inducers of "%s"' % (was_linked,t_n))
            all_compounds.update(self.disease_inducers)

        if all_compounds:
            colname = 'compounds modulating '+t_n
            print ('\n\nLinking compounds modulating "%s" to targets with unknown state in "%s"' % (t_n,t_n))
            how2connect = self.set_how2connect(['DirectRegulation'],[],'',['Binding','Expression','Regulation'])
            was_linked,linked_entities,unk_targets_df = self.link2concept(
                colname,list(all_compounds),unk_targets_df,how2connect,REFERENCE_IDENTIFIERS)
            print('Linked %d targets to inducers for %s' % (was_linked,t_n))


        uptarget_df = self.make_count_df(activ_targets_df,DF4ANTAGONISTS)
        self.add2raw(uptarget_df)
        downtarget_df = self.make_count_df(inhib_targets_df,DF4AGONISTS)
        self.add2raw(downtarget_df)
        unktarget_df = self.make_count_df(unk_targets_df,DF4UNKNOWN)
        self.add2raw(unktarget_df)
        return


    def normalize_counts(self):
        self.normalize(DF4ANTAGONISTS,ANTAGONIST_TARGETS_WS)
        self.report_pandas[ANTAGONIST_TARGETS_WS].tab_format['tab_color'] = 'red'
        self.normalize(DF4AGONISTS,AGONIST_TARGETS_WS)
        self.report_pandas[AGONIST_TARGETS_WS].tab_format['tab_color'] = 'blue'
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
        input_disease_graph = self.Graph.neighborhood(set(self.input_diseases))
        super().add_ps_bibliography(from_graph=input_disease_graph)


    def add_target_annotation(self,_2column='Class'):
        self.add_entity_annotation(_2column,self.report_pandas[ANTAGONIST_TARGETS_WS])
        self.add_entity_annotation(_2column,self.report_pandas[AGONIST_TARGETS_WS])
        self.add_entity_annotation(_2column,self.report_pandas[UNKNOWN_TARGETS_WS])


    def input_disease_df(self):
        rows = list()
        for d in self.input_diseases:
            d_childs = d.childs()
            child_names = [c.name() for c in d_childs]
            child_names_str = ','.join(child_names)
            rows.append((d.name(),str(d['Connectivity'][0]),child_names_str,d.urn()))

        input_df = df.from_rows(rows,header=['Name','Connectivity','Children','URN'])
        input_df['Connectivity'] = input_df['Connectivity'].astype(int)
        input_df = df.from_pd(input_df.sort_values(by='Connectivity',ascending=False))
        input_df._name_ = 'Disease subtypes'
        return input_df


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
        self.add2report(self.input_disease_df())
    
        self.clear()
        print('Target ranking for %s was done in %s' % 
        (self._disease2str(), self.execution_time(start_time)[0]))

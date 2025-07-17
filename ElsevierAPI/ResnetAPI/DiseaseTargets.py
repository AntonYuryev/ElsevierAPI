from .ResnetGraph import PROTEIN_TYPES,EFFECT,PHYSICAL_INTERACTIONS,ResnetGraph,PSRelation,PSObject
from .PSPathway import PSPathway
from .NetworkxObjects import ACTIVATED,REPRESSED,UNKNOWN_STATE,OBJECT_TYPE,MECHANISM
from .SemanticSearch import SemanticSearch,len,OQL,RANK,execution_time
from .ResnetAPISession import NO_REL_PROPERTIES,ONLY_REL_PROPERTIES,REFERENCE_IDENTIFIERS,BIBLIO_PROPERTIES,DBID
from .FolderContent import FolderContent
from ..pandas.panda_tricks import df
from ..utils import unpack
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
        my_kwargs = {
                'disease':[],
                'what2retrieve':BIBLIO_PROPERTIES,
                'ent_props' : ['Name', 'Class'], #'Class' is for target partners retrieval
                'rel_props' : [EFFECT],
                'data_dir' : '',
                'add_bibliography' : True,
                'strict_mode' : False,
                'target_types' : ['Protein'],
                'pathway_folders' : [],
                'pathways' : [],
                "max_childs" : 11,
                "add_closeness":True, # obsolete
                'propagate_target_state_in_model':True,
                'add_regulators4':dict(),
                'add_targets4':dict(),
                'skip':False
            }
        
        entprops = kwargs.pop('ent_props',[])
        relprops = kwargs.pop('rel_props',[])

        my_kwargs.update(kwargs)
        super().__init__(*args,**my_kwargs)
        self.add_ent_props(entprops)
        self.add_rel_props(relprops)

        self.columns2drop += [self.__resnet_name__,self.__mapped_by__,'State in Disease']
        self.max_threads4ontology = 50

        self.input_diseases = list()
        self.input_symptoms = list()
        self.input_clinpars = list()
        self.input_cellprocs = list()

        self.__targets__ = set()
        self.expanded_targets = dict({}) # {target_neighbor_name:target_name}
        self.target_inhibitorsG= ResnetGraph() 
        # holds proteins, which modulation must inhibit disease targets listed in self.params["add_inhibitors4"] 
        # self.params["add_inhibitors4"] must be a dict {reltype:[target_names]}
        # parameter was introduced to find drugs that activate "Na+/K+-exchanging ATPase" only post-translationally
        # while inhibiting its expression
        self.GVs = list()
        self.targets4strictmode = set()
        self.ct_drugs = set()
        self.drugs_linked2disease = set()
        self.disease_inducers = list()
        
        self.gv2diseases = ResnetGraph()
        self.partner2targets = ResnetGraph()
        self.disease_pathways = dict({}) # {cell_type:PSPathway}

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
        assert (isinstance(self.gv2diseases, ResnetGraph))
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


    def input_disease_names(self)->list[str]:
      if self.input_diseases:
        prop2values = {'Name':self.params['disease'], 'Alias':self.params['disease']}
        diseas_names = list(set([x['Name'][0] for x in self.input_diseases if x.has_value_in(prop2values)]))
        return diseas_names
      else:
        try:# case of SNEA drug_df
          return [self.params['sample']+' in '+ self.params['experiment']]
        except KeyError:
          return ['Drugs']
            

    def names4tmsearch(self):
      search_names = list(self.params['disease'])
      for reltype,targets in self.params.get('add_regulators4',dict()).items():
        for target in targets:# need to preserve the order of target names to support hierachal concept search
          if target not in search_names:
            search_names.append(target)

      for reltype,targets in self.params.get('add_targets4',dict()).items():
        for target in targets:# need to preserve the order of target names to support hierachal concept search
          if target not in search_names:
            search_names.append(target)
            
      return search_names


    def report_name(self):
        rep_pred = ' suggested ' if self.params['strict_mode'] else ' predicted ' 
        return ','.join(self.params['disease'])+rep_pred+'targets'


    def _disease2str(self):
        # do not put disease name in quotas - it is used to create a file name
        if len(self.input_diseases) > 1:
          return ','.join(self.params['disease'])+' ('+str(len(self.input_diseases))+' types)'
        else:
          return ','.join(self.params['disease'])


    def _targets(self)->set[PSObject]:
        return self.targets4strictmode if self.params['strict_mode'] else self.__targets__
    
    def target_types(self):
        return self.params['target_types']

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
        assert isinstance(self.gv2diseases,ResnetGraph)
        self.GVs = self.gv2diseases.psobjs_with(only_with_values=['GeneticVariant'])
    
        GVdbids = ResnetGraph.dbids(self.GVs)
        target_withGVs = list()

        if GVdbids:
            oql_query = f'SELECT Relation WHERE objectType = GeneticChange AND \
                NeighborOf (SELECT Entity WHERE objectType = ({self.__target_types_str()})) AND NeighborOf ({OQL.get_objects(GVdbids)})'
            request_name = f'Find targets with genetic variants linked to {disease_names}'

            # avoid adding GeneticChange to self.Graph
            new_session = self._clone(to_retrieve=NO_REL_PROPERTIES,init_refstat=False)
            gv2target = new_session.process_oql(oql_query,request_name)
            assert isinstance(gv2target, ResnetGraph)
            self.gv2diseases = self.gv2diseases.compose(gv2target)
            target_withGVs = gv2target.psobjs_with(only_with_values=PROTEIN_TYPES)
            gv2dis_refs = self.gv2diseases.load_references()
            print('Found %d targets with %d GeneticVariants for %s supported by %d references' % 
                (len(target_withGVs), len(self.GVs), disease_names, len(gv2dis_refs)))
            return target_withGVs
        else:
            return ResnetGraph()

 
    def targets_from_db(self):
        req_name = f'Find targets linked to {self._disease2str()}'
        select_targets = f'SELECT Entity WHERE objectType = ({self.__target_types_str()})'
        OQLquery = f'SELECT Relation WHERE NeighborOf ({select_targets}) AND NeighborOf ({self.find_disease_oql})'      
        target_disease_graph = self.process_oql(OQLquery,req_name)
        assert(isinstance(target_disease_graph,ResnetGraph))

        return target_disease_graph
        
    
    def metabolite2inhibit(self):
      metabolites2inhibit = self.params.get('metabolites2inhibit',[])
      if metabolites2inhibit:
        metabolites2inhibit = ','.join(self.params['metabolites2inhibit'])
        select_metabolite = f'SELECT Entity WHERE Name = ({metabolites2inhibit})'
        oql = f'SELECT Relation WHERE NeighborOf ({self.find_disease_oql}) AND NeighborOf ({select_metabolite})'
        my_session = self._clone_session()
        metbolite2diseaseG = my_session.process_oql(oql,'Find metabolites2inhibit')
        [metbolite2diseaseG.set_edge_annotation(r.uid(),t.uid(),rel.urn(),EFFECT,['positive']) for r,t,rel in metbolite2diseaseG.iterate()]
        self.Graph.add_graph(metbolite2diseaseG) # to ensure propert targets state in set_target_stae
        metabolites2inhibit_objs = metbolite2diseaseG._psobjs_with('SmallMol',OBJECT_TYPE)
        return metabolites2inhibit_objs
      else:
        return list()
        

    def expand_targets(self,reltype2targets:dict[str,list[str]],dir='upstream'):
      '''
      input:
        reltype2targets = {reltype:[targets]}
        dir = upstream, downstream
      output:
        updated self.expanded_targets
      '''
      expansion = 'regulators' if dir=='upstream' else 'targets'
      for reltype,t_names in reltype2targets.items():
        targets2expand = [n for n in self.__targets__ if n.name() in t_names]
        neighbs = 'regulators' if dir == 'upstream' else 'targets'
        t_dbids = ResnetGraph.dbids(targets2expand)
        oql_query = OQL.expand_entity(t_dbids,['Id'],[reltype],PROTEIN_TYPES,dir) 
        if reltype == 'ProtModification': # ProtModification with Mechanism = cleavage is Expression regulation
          oql_query += ' AND NOT (Mechanism = cleavage)'
        reqname = f'Find {reltype} {neighbs} 4 {t_names}'
        expanded_targetsG = self.process_oql(oql_query,reqname)
        if reltype == 'Expression': # ProtModification with Mechanism = cleavage is Expression regulation with Effect negative
          oql_query = OQL.expand_entity(t_dbids,['Id'],['ProtModification'],PROTEIN_TYPES,dir) 
          oql_query += ' AND Mechanism = cleavage'
          reqname = f'Find {reltype} {neighbs} 4 {t_names}'
          new_session = self._clone_session()
          proteasesG = new_session.process_oql(oql_query,reqname)
          if proteasesG:
            proteasesGrels = proteasesG._psrels()
            for rel in proteasesGrels:
              rel[OBJECT_TYPE] = ['Expression']
              rel[EFFECT] = ['negative']
              rel.urn(refresh=True)
            expanded_targetsG = expanded_targetsG.compose(ResnetGraph.from_rels(proteasesGrels))

        if expanded_targetsG:
          expansion_descr = f'{reltype} {expansion} of {t_names}'
          expanded_targetsG.name = expansion_descr
          rel_rank = ['DirectRegulation','Binding','ProtModification','MolTransport','Regulation','PromoterBinding','Expression']
          # assume that Regulation means post-translational regulation
          expanded_targetsG = expanded_targetsG.make_simple(rel_rank)
          self.__targets__.update(expanded_targetsG._get_nodes())
          for t in targets2expand:
            t_neighbors = expanded_targetsG.get_neighbors({t})
            neigh2target = dict({n.name():t.name() for n in t_neighbors})       
            self.expanded_targets.update(neigh2target)
          print(f'Added {len(neigh2target)} {expansion_descr}')
          self.Graph = self.Graph.compose(expanded_targetsG)


    def find_targets(self):
      target_disease_graph = self.targets_from_db()
      if isinstance(target_disease_graph,ResnetGraph):
        self.__targets__ = set(target_disease_graph.psobjs_with(only_with_values=self.target_types()))
        target_refcount = target_disease_graph.load_references()
        print('Found %d targets linked to %s supported by %d references in database' 
                    % (len(self.__targets__),self._disease2str(), len(target_refcount)))
        
        target_withGVs = self.__GVtargets()
        self.__targets__.update(target_withGVs)
        self.targets4strictmode.update(target_withGVs)

        metabolites2inhibit = self.metabolite2inhibit()
        self.__targets__.update(metabolites2inhibit)
        self.targets4strictmode.update(metabolites2inhibit)
          
      # adding targets from disease model pathways
        disease_pathways = self.load_pathways()
        disease_model_components = set()
        if disease_pathways:
            all_pathways = [p for p in disease_pathways.values()]
            [disease_model_components.update(p.get_members(self.target_types())) for p in all_pathways]
            self.__targets__.update(disease_model_components)
            before_add = self.Graph.number_of_nodes()
            self.Graph.add_psobjs(set(disease_model_components))
            print('Added %d targets from disease model' % (self.Graph.number_of_nodes()-before_add))

        if self.params['strict_mode']:
            regulators2disease_graph = target_disease_graph.subgraph_by_relprops(['Regulation'])
            disease_regulators = regulators2disease_graph.psobjs_with(only_with_values=self.target_types())
            self.targets4strictmode.update(disease_regulators)
            self.targets4strictmode.update(disease_model_components)
        else:
          target2expand_names = dict(self.params.get("add_regulators4",dict()))
          self.expand_targets(target2expand_names,'upstream')
          target2expand_names = dict(self.params.get("add_targets4",dict()))
          self.expand_targets(target2expand_names,'downstream')
         
          reltype2targets = dict(self.params.get("add_inhibitors4",dict()))
          for reltype,target_names in reltype2targets.items():
            targets2expand = [n for n in self.__targets__ if n.name() in target_names]
            targets2expand_dbids = ResnetGraph.dbids(targets2expand)
            reltypes = [reltype] if reltype else []
            oql_query = OQL.expand_entity(targets2expand_dbids,['Id'],reltypes,PROTEIN_TYPES,'upstream')
            oql_query += ' AND Effect = (positive,negative)'
            self.target_inhibitorsG = self.process_oql(oql_query,f'Find inhibitors 4 {target_names}')
            
            if reltype == 'Expression': # ProtModification with Mechanism = cleavage is Expression regulation with Effect negative
              oql_query = OQL.expand_entity(targets2expand_dbids,['Id'],['ProtModification'],PROTEIN_TYPES,'upstream') 
              oql_query += ' AND Mechanism = cleavage'
              new_session = self._clone_session()
              proteasesG = new_session.process_oql(oql_query,f'Find proteases 4 {target_names}')
              if proteasesG:
                proteasesGrels = proteasesG._psrels()
                for rel in proteasesGrels:
                  rel[OBJECT_TYPE] = ['Expression']
                  rel[EFFECT] = ['negative']
                  rel.urn(refresh=True)
                self.target_inhibitorsG = self.target_inhibitorsG.add_graph(ResnetGraph.from_rels(proteasesGrels))
            
            if self.target_inhibitorsG:
              before_add = len(self.__targets__) 
              self.__targets__.update(self.target_inhibitorsG._get_nodes())
              print(f'Added {len(self.__targets__) - before_add} inhibitors of {target_names} {reltype}')
          return
 

    def find_drugs(self):      
        request_name = 'Find compounds inhibiting {}'.format(self._disease2str())
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = negative \
            AND NeighborOf (SELECT Entity WHERE objectType = SmallMol) AND NeighborOf ({})'
        OQLquery = OQLquery.format(self.find_disease_oql)
        drugs_graph = self.process_oql(OQLquery,request_name)
        if isinstance(drugs_graph,ResnetGraph):
            self.drugs_linked2disease = set(drugs_graph._psobjs_with('SmallMol','ObjTypeName'))
            drugs_references = drugs_graph.load_references()
            print('Found %d compounds inhibiting %s supported by %d references' % 
                        (len(self.drugs_linked2disease),self._disease2str(), len(drugs_references)))

        request_name = 'Find clinical trials for {}'.format(self._disease2str())
        OQLquery = 'SELECT Relation WHERE objectType = ClinicalTrial \
            AND NeighborOf (SELECT Entity WHERE objectType = SmallMol) AND NeighborOf ({})'
        OQLquery = OQLquery.format(self.find_disease_oql)
        ct_graph = self.process_oql(OQLquery,request_name)
        if isinstance(ct_graph,ResnetGraph):
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
        if isinstance(induction_graph,ResnetGraph):
            self.disease_inducers = induction_graph._psobjs_with('SmallMol','ObjTypeName')
            count = len(self.disease_inducers)
            self.disease_inducers ={x for x in self.disease_inducers if x not in self.ct_drugs}
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
        new_session = self._clone(to_retrieve=NO_REL_PROPERTIES,init_refstat=False)
        request_name = f'Find ligands for {str(len(receptor_dbids))} receptors'
        p2t = new_session.process_oql(ligands_oql,request_name)
        if isinstance(p2t,ResnetGraph):
            self.partner2targets = p2t
        new_ligands = self.partner2targets._psobjs_with('Ligand','Class')
        print('Found %d additional ligands for %d receptors' % (len(new_ligands),len(receptor_dbids)))
        new_session.close_connection()
      else:
        find_receptors = ''
        print('No Receptor Class entities were found among %s targets' % self.target_types())
        
      ligand_dbids = self.Graph.dbids4nodes(['Ligand'],['Class'])
      partners_dbids = ResnetGraph.dbids(self.partner2targets._get_nodes())
      orphan_ligands_dbids = set(ligand_dbids).difference(partners_dbids)
      if orphan_ligands_dbids:
        find_orphan_ligands = OQL.get_objects(list(orphan_ligands_dbids))
        receptor_oql = f'SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive AND RelationNumberOfReferences > {ref_cutoff} \
            AND NeighborOf upstream (SELECT Entity WHERE Class = Receptor) AND NeighborOf downstream ({find_orphan_ligands})'
        request_name = f'Find receptors for {str(len(orphan_ligands_dbids))} orphan ligands'
        #opening new APIsession to avoid adding result graph to self.Graph
        new_session = self._clone(to_retrieve=NO_REL_PROPERTIES,init_refstat=False)
        orphan_ligand_partners = new_session.process_oql(receptor_oql,request_name)
        if isinstance(orphan_ligand_partners,ResnetGraph):
            self.partner2targets.add_graph(orphan_ligand_partners)
            new_receptors = orphan_ligand_partners._psobjs_with('Receptor','Class')
            print('Found %d receptors for %d orphan ligands' % (len(new_receptors),len(orphan_ligands_dbids)))
        new_session.close_connection()

      if find_receptors:
        find_metabolites = OQL.get_childs(['mammal endogenous compounds and their derivatives'],['Name'])
        metabolite_ligands_oql = f'SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive \
            AND NeighborOf downstream ({find_metabolites}) AND NeighborOf upstream ({find_receptors}) AND RelationNumberOfReferences > 25'
        new_session = self._clone(to_retrieve=NO_REL_PROPERTIES,init_refstat=False)
        metabolite2receptors_graph = new_session.process_oql(metabolite_ligands_oql,'Find metabolite ligands')
        if isinstance(metabolite2receptors_graph, ResnetGraph):
          self.partner2targets.add_graph(metabolite2receptors_graph)
          metabolite_ligands = metabolite2receptors_graph._psobjs_with('SmallMol','ObjTypeName')
          print('Found %d metabolite ligands for %d receptors' % (len(metabolite_ligands),len(receptor_dbids)))
        new_session.close_connection()

      find_metabolite_products = 'SELECT Relation WHERE objectType = (MolTransport, MolSynthesis,ChemicalReaction) \
AND Effect = positive AND NeighborOf downstream (SELECT Entity WHERE id =({ids})) AND \
NeighborOf upstream (SELECT Entity WHERE objectType = SmallMol) AND RelationNumberOfReferences > 50'
      req_name = f'Finding downstream metabolite products of targets for {self._disease2str()}'
      new_session = self._clone(to_retrieve=NO_REL_PROPERTIES,init_refstat=False)
      metbolite_products_graph = new_session.iterate_oql(find_metabolite_products,set(self._targets_dbids()),request_name=req_name)
      self.partner2targets.add_graph(metbolite_products_graph)
      metabolite_products = metbolite_products_graph._psobjs_with('SmallMol','ObjTypeName')
      print('Found %d metabolites produced by %d targets of %s' % (len(metabolite_products),len(self.__targets__),self._disease2str()))
      new_session.close_connection()

      # removing promiscous partners:
      self.partner2targets.remove_nodes_by_targets(max_target=10,having_values=['SmallMol'])
      self.partner2targets.remove_nodes_by_targets(max_target=10,only_with_prop=['Class'], having_values=['Ligand'])
      self.partner2targets.remove_nodes_by_regulators(max_regulator=10,having_values=['SmallMol'])
      return


    def load_pathways(self)->dict[str,PSPathway]:
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

      if self.params.get('pathways',''):
          filtered_pathway_list = list()
          for ps_pathway in disease_pathways:
              assert(isinstance(ps_pathway, PSPathway))
              if ps_pathway.name() in self.params['pathways']:
                  filtered_pathway_list.append(ps_pathway)
          disease_pathways = filtered_pathway_list

      print('Found %d curated pathways for %s:' %(len(disease_pathways), self._disease2str()))
      [print(p.name()+'\n') for p in disease_pathways]

      # merging pathways from the same celltype to build cell type specific models
      for pathway in disease_pathways:
          try:
              cell_types = pathway['CellType']
          except KeyError:
              cell_types = ['disease']

          for cell_type in cell_types:
              try:
                  exist_pathway = self.disease_pathways[cell_type]
                  assert(isinstance(exist_pathway,PSPathway))
                  exist_pathway.merge_pathway(pathway)
                  self.disease_pathways[cell_type] = exist_pathway
              except KeyError:
                  self.disease_pathways[cell_type] = pathway
          
      for pathway in self.disease_pathways.values():
          assert(isinstance(pathway,PSPathway))
          pathway.graph.remove_nodes_by_prop(['CellProcess', 'Disease','Treatment'])

      return self.disease_pathways

        
    def init_semantic_search(self):
        print('\n\nInitializing semantic search')
        if not self.__targets__:
            print ('No targets found for %s' % self._disease2str())
            print('Consider setting "strict_mode" to False')
            return False

        targets4ranking = list(self._targets())
        self.RefCountPandas = self.load_df(targets4ranking,max_childs=self.params['max_childs'],max_threads=10)
        print('Will score %d targets linked to %s' % (len(self.RefCountPandas),self._disease2str()))
        self.entProps = ['Name'] # Class is no longer needed
        return True


    def score_GVs(self):
        print('\n\nScoring targets by number of semantic references linking their Genetic Variants to %s' % self._disease2str())
        if not self.GVs:
            print('%s has no known Genetic Variants' % self._disease2str())
            return

        concept_name = 'GVs for '+self._disease2str()
        self.__colnameGV__ = concept_name
        refcount_column = self._refcount_colname(concept_name)
        weighted_refcount_column = self._weighted_refcount_colname(concept_name) 
        linked_count_column = self._linkedconcepts_colname(concept_name)
        concept_size_column = self._concept_size_colname(concept_name)
        self.RefCountPandas.insert(len(self.RefCountPandas.columns),weighted_refcount_column,[float(0)]*len(self.RefCountPandas))
        self.RefCountPandas.insert(len(self.RefCountPandas.columns),refcount_column,[0]*len(self.RefCountPandas))
        self.RefCountPandas.insert(len(self.RefCountPandas.columns),linked_count_column,[0]*len(self.RefCountPandas))
        self.RefCountPandas.insert(len(self.RefCountPandas.columns),concept_size_column,[len(self.GVs)]*len(self.RefCountPandas))

        target_gvlinkcounter = 0
        assert (isinstance(self.gv2diseases, ResnetGraph))
        for i in self.RefCountPandas.index:
            target_dbids = list(self.RefCountPandas.at[i,self.__temp_id_col__])
            row_targets = self.Graph.psobj_with_dbids(set(target_dbids))
            targetGVs = self.gv2diseases.get_neighbors(set(row_targets), allowed_neigbors=self.GVs)
                
            if targetGVs:
                GVnames = set([n.name() for n in targetGVs])
                self.RefCountPandas.at[i,self.__colnameGV__] = ';'.join(GVnames)
                target_gvlinkcounter += 1
                gv_disease_subgraph = self.gv2diseases.get_subgraph(list(targetGVs),self.input_diseases)
                GVrefcount = len(gv_disease_subgraph.load_references())
            else:
                GVrefcount = 0
            
            GVscore = GVrefcount*(1.0 + len(targetGVs)/len(self.GVs))
            self.RefCountPandas.at[i,weighted_refcount_column] = GVscore
            self.RefCountPandas.at[i,refcount_column] = GVrefcount
            self.RefCountPandas.at[i,linked_count_column] = len(targetGVs)

        if target_gvlinkcounter:
          self._set_rank(self.RefCountPandas,concept_name)
        print('Found %d targets linked to %d GVs' % (target_gvlinkcounter, len(self.GVs)))


    @staticmethod
    def __get_state(fromG:ResnetGraph):
      positive_refs, negative_refs = fromG._effect_counts__()
      if len(positive_refs) > len(negative_refs):
        return ACTIVATED
      elif len(negative_refs) > len(positive_refs):
        return REPRESSED
      else:
        return UNKNOWN_STATE


    def __vote4effect(self,targets:list[PSObject]):
        target2disease = self.Graph.get_subgraph(targets, self.input_diseases)
        state = self.__get_state(target2disease)
        if state != UNKNOWN_STATE:
          return state
        
        # attempting to deduce effect from target_partners
        target_partners = list(self.partner2targets.get_neighbors(set(targets)))
        partners2disease_graph = self.Graph.get_subgraph(target_partners, self.input_diseases)
        state = self.__get_state(partners2disease_graph)
        if state != UNKNOWN_STATE:
          return state
      
        # attempting to deduce effect from GeneticChange
        if partners2disease_graph.psrels_with(['GeneticChange']) or partners2disease_graph.psrels_with(['GeneticChange']):
          return REPRESSED
        else:
          return UNKNOWN_STATE


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

        if self.disease_model and 'propagate_target_state_in_model' in self.params.keys():
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
                activated_targets = self.disease_model.psobj_with_dbids(self._all_dbids(activ_targets_pd))
                inhib_targets_pd = self.RefCountPandas.loc[(self.RefCountPandas['State in Disease'] == REPRESSED)]
                inhibited_trgts = self.disease_model.psobj_with_dbids(self._all_dbids(inhib_targets_pd))

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
    
        if self.expanded_targets:
          for idx in self.RefCountPandas.index:
            if self.RefCountPandas.at[idx,'State in Disease'] == UNKNOWN_STATE:
              target_name = self.RefCountPandas.at[idx,'Name']
              row_targets = [n for n in self.__targets__ if n.name() == target_name]
              try:
                expanded_name = self.expanded_targets[target_name]
                expanded_state = int(self.RefCountPandas.loc[(self.RefCountPandas['Name'] == expanded_name)]['State in Disease'].iloc[0])
                expanded_target = [n for n in self.__targets__ if n.name() == expanded_name]
                my_subgraph = self.Graph.neighborhood(set(expanded_target),row_targets)
                my_subgraph.clone_node(expanded_target[0],self.input_diseases[0],expanded_state,'Regulation')
                self.Graph = my_subgraph.compose(self.Graph)
                state = self.__vote4effect(row_targets)
                self.RefCountPandas.at[idx,'State in Disease'] = int(state)
              except KeyError:
                continue

        # setting disease state for proteins which modulation must inhibit disease targets listed in self.params["add_inhibitors4"]
          if self.target_inhibitorsG:
            inhibitors = self.target_inhibitorsG.regulators()
            inhibitor_names = ResnetGraph.names(inhibitors)
            for idx in self.RefCountPandas.index:
              if self.RefCountPandas.at[idx,'State in Disease'] == UNKNOWN_STATE:
                drugtarget_name = self.RefCountPandas.at[idx,'Name']
                if drugtarget_name in inhibitor_names:
                  row_inhibitors = [x for x in inhibitors if x.name() == drugtarget_name]
                  my_subgraph = self.target_inhibitorsG.neighborhood(row_inhibitors,in_direction='>')
                  if my_subgraph:
                    self.RefCountPandas.at[idx,'State in Disease'] = self.__get_state(my_subgraph)
        return


    def score_partners(self):
        print(f'\n\nFinding semantic references supporting links between target partners and {self._disease2str()}')
        my_df = self.RefCountPandas.dfcopy()
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
        how2connect = self.set_how2connect(**dict())
        colname = 'target partners'
        linked_row_count,_,partners_df = self.link2concept(colname,self.input_diseases,partners_df,how2connect)
        print('%d targets have partners linked to %s' % (linked_row_count,self._disease2str()))
        if linked_row_count:
          partners_df = df.from_pd(partners_df.drop(columns=[self.__temp_id_col__]))
          self._set_rank(self.RefCountPandas,colname)
          return self.RefCountPandas.merge_df(partners_df,on='Name')
        else:
          return self.RefCountPandas


    def make_disease_network(self):
        print(f'Creating physical interaction network between targets of {self.input_disease_names()} \
              to calculate target closeness for "score_regulators" function')
        my_session = self._clone_session(to_retrieve=NO_REL_PROPERTIES,init_refstat=False)
        disease_network = my_session.get_ppi(self.__targets__, self.params.get('ppiRNEFs',[]))
        disease_network.name = f'{self._disease2str()} PPPI network'
        return disease_network.make_simple(['DirectRegulation','ProtModification','Binding']) 
    

    def DiseaseNetworkRegulationScore(self):
        '''
        Add
        ----
        PATHWAY_REGULATOR_SCORE column to self.RefCountPandas
        '''
        print('\n\nScoring regulators by distance to components of disease pathways',flush=True)
        print('Retrieving regulatory network between targets and components of disease pathway')
        regulation_graph = self.disease_model if self.disease_model else self.make_disease_network()
        uid2closeness = regulation_graph.closeness() # will use closeness in disease network to initialize regulatory ranking
        
        disease_pathway_component_dbids = set(regulation_graph.dbids4nodes())
        connected_nodes_dbids = disease_pathway_component_dbids
        unconnected_nodes_dbids = self.all_entity_dbids

        for step in range (0,5): # 5 degree of separation
            new_session = self._clone(to_retrieve=NO_REL_PROPERTIES,init_refstat=False) # new_session.Graph does not contain regulation_graph 
            new_session.Graph = regulation_graph.compose(new_session.Graph) # to speadup connect_nodes

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

        if self.disease_pathways:
            for pathway in self.disease_pathways.values():
                assert(isinstance(pathway,PSPathway))
                uid2clos = pathway.graph.closeness()
                regulation_graph.rank_regulators(uid2clos,PATHWAY_REGULATOR_SCORE)
        else:
            regulation_graph.rank_regulators(uid2closeness,PATHWAY_REGULATOR_SCORE)

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

        my_rank = self.RefCountPandas.max_colrank()
        self.RefCountPandas.col2rank[PATHWAY_REGULATOR_SCORE] = my_rank + 1
        return


    def combine_processes(self):
      processes2inhibit = self.params.get('processes2inhibit',dict())
      if isinstance(processes2inhibit,list):
        processes2inhibit = {p:[1.0] for p in processes2inhibit}
      processes2activate = self.params.get('processes2activate',dict())
      if isinstance(processes2activate,list):
        processes2activate = {p:[1.0] for p in processes2activate}

      processes = self.params.get('processes',dict())
      if isinstance(processes,list):
        processes = {p:[1.0] for p in processes}
      processes.update(processes2inhibit)
      processes.update(processes2activate)
      self.params['processes'] = processes


    def score_target_semantics(self):
      #  do not multithread.  self.Graph will leak
      disease_str = self._disease2str()
      self.score_GVs()
      self.set_target_disease_state()

      kwargs = {'connect_by_rels':['Regulation'],
                'boost_with_reltypes':['FunctionalAssociation'],
                'column_name':'Regulate '+ disease_str
                }
      self.RefCountPandas = self.score_concepts(self.input_diseases,**kwargs)[2]

      kwargs = {'connect_by_rels':['GeneticChange'],
                'boost_with_reltypes':['FunctionalAssociation'],
                'column_name':'Genetically linked to '+ disease_str}
      self.RefCountPandas = self.score_concepts(self.input_diseases,**kwargs)[2]

      kwargs = {'connect_by_rels':['Biomarker'],
                'boost_with_reltypes':['FunctionalAssociation'],
                'column_name':'Target is Biomarker in '+disease_str}
      self.RefCountPandas = self.score_concepts(self.input_diseases,**kwargs)[2]

      kwargs = {'connect_by_rels':['QuantitativeChange'],
                'boost_with_reltypes':['FunctionalAssociation'],
                'column_name':'Quantitatively changed in '+ disease_str}
      self.RefCountPandas = self.score_concepts(self.input_diseases,**kwargs)[2]

      kwargs = {'connect_by_rels':['Regulation','QuantitativeChange','StateChange','Biomarker'],
                'boost_with_reltypes':['FunctionalAssociation','CellExpresion'],
                'column_name':'symptoms for '+ disease_str,
                'clone2retrieve' : REFERENCE_IDENTIFIERS}
      self.RefCountPandas = self.score_concept('symptoms',**kwargs)[2]

      kwargs = {'connect_by_rels':['Regulation'],
                'boost_with_reltypes':['FunctionalAssociation'],
                'column_name':'Clinical parameters for '+ disease_str,
                'clone2retrieve' : REFERENCE_IDENTIFIERS}
      self.RefCountPandas = self.score_concept('clinical_parameters',**kwargs)[2]

      self.combine_processes()

      kwargs = {'connect_by_rels':['Regulation'],
                'boost_with_reltypes':['FunctionalAssociation'],
                'column_name':'Cell processes affected by '+ disease_str,
                'clone2retrieve' : REFERENCE_IDENTIFIERS}
      self.RefCountPandas = self.score_concept('processes',**kwargs)[2]

      self.RefCountPandas = self.score_partners()
      self.DiseaseNetworkRegulationScore()
      ################## SPLITTING RefCountPandas to score agonists and antagonists differently #####################
      activ_targets_df = df.from_pd(self.RefCountPandas.loc[(self.RefCountPandas['State in Disease'] > UNKNOWN_STATE)])
      activ_targets_df.copy_format(self.RefCountPandas)
      inhib_targets_df = df.from_pd(self.RefCountPandas.loc[(self.RefCountPandas['State in Disease'] < UNKNOWN_STATE)])
      inhib_targets_df.copy_format(self.RefCountPandas)
      unk_targets_df = df.from_pd(self.RefCountPandas.loc[(self.RefCountPandas['State in Disease'] == UNKNOWN_STATE)])
      unk_targets_df.copy_format(self.RefCountPandas)
      print(f'Created 3 worksheets with: {len(activ_targets_df)} activated targets, {len(inhib_targets_df)} inhibited_targets, {len(unk_targets_df)} unknown state targets')
      # drug-target relations are not added to self.Graph to exclude them from bibliography
      kwargs = {'connect_by_rels':['DirectRegulation'],
                'with_effects':['negative'],
                'boost_with_reltypes':['Binding','Expression','Regulation'],
                'column_name':'drugs for '+disease_str,
                'clone2retrieve' : REFERENCE_IDENTIFIERS}
      linked_rows,_,activ_targets_df = self.score_concepts(self.drugs_linked2disease,activ_targets_df,**kwargs)

      kwargs['with_effects'] = ['positive']
      if linked_rows:
        kwargs['column_rank'] = activ_targets_df.max_colrank()
      kwargs['column_name'] = 'inducers of '+disease_str
      lr,_,activ_targets_df = self.score_concepts(self.disease_inducers,activ_targets_df,**kwargs)
      linked_rows += lr

      linked_rows,_,inhib_targets_df = self.score_concepts(self.drugs_linked2disease,inhib_targets_df,**kwargs)
      if linked_rows:
        kwargs['column_rank'] = inhib_targets_df.max_colrank()
      kwargs['with_effects'] = ['negative']
      inhib_targets_df = self.score_concepts(self.disease_inducers,inhib_targets_df,**kwargs)[2]

      all_compounds = self.disease_inducers|self.drugs_linked2disease

      if all_compounds:
        print (f'\n\nLinking compounds modulating "{disease_str}" to targets with unknown state in "{self._disease2str()}"')
        kwargs = {'connect_by_rels':['DirectRegulation'],
                  'boost_with_reltypes':['Binding','Expression','Regulation'],
                  'column_name':'compounds modulating '+disease_str,
                  'clone2retrieve' : REFERENCE_IDENTIFIERS,
                  'column_rank': unk_targets_df.max_colrank()}
        unk_targets_df = self.score_concepts(self.disease_inducers,unk_targets_df,**kwargs)[2]

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
        

    def add_bibliography4targets(self):
      '''
      Adds
      ----
      etm references column to ANTAGONIST_TARGETS_WS,AGONIST_TARGETS_WS,UNKNOWN_TARGETS_WS worksheets,\n
      adds ETM_BIBLIOGRAPHY worksheet to report
      '''
      print('Adding ETM bibliography for ranked targets', flush=True)
      #etm_refcount_colname = self.etm_refcount_colname('Name',self.input_disease_names())
      self.refs2report(ANTAGONIST_TARGETS_WS,self.names4tmsearch())
      self.refs2report(AGONIST_TARGETS_WS,self.names4tmsearch())
      self.refs2report(UNKNOWN_TARGETS_WS,self.names4tmsearch())
      self.add_tm_bibliography_df()


    def add_graph_bibliography(self):
        input_disease_graph = self.Graph.neighborhood(set(self.input_diseases))
        super().add_graph_bibliography(from_graph=input_disease_graph)


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

        if len(rows) > 1:
          input_df = df.from_rows(rows,header=['Name','Connectivity','Children','URN'])
          input_df['Connectivity'] = input_df['Connectivity'].astype(int)
          input_df = input_df.sortrows(by='Connectivity')
          input_df._name_ = 'Disease subtypes'
          return input_df
        else: return df()
    

    def paralog_adjustment(self):
      my_apisession = self._clone_session()
      my_apisession.entProps = ['Name']
      my_apisession.relProps = ['Similarity']
      my_targets = self._targets()
      for ws_name in [ANTAGONIST_TARGETS_WS, AGONIST_TARGETS_WS]:
        try:
          ws = self.report_pandas[ws_name]
          assert(isinstance(ws,df))
          name2score = {name:float(score) for name,score in zip(ws['Name'],ws[RANK]) if score} # WEIGHTS Name is skipped
          ws_target_names = set(ws['Name'].to_list())
          ws_targets = [x for x in my_targets if x.name() in ws_target_names]
          ws_paralog_graph = my_apisession.get_network(set(ResnetGraph.dbids(ws_targets)),['Paralog'])
          assert(isinstance(ws_paralog_graph,ResnetGraph))
          paralog_counter = 0
          for p1,p2,rel in ws_paralog_graph.edges.data('relation'):
              assert(isinstance(rel,PSRelation))
              similarity = float(rel.PropSetToProps[0]['Similarity'][0])
              if similarity:
                paralog1 = ws_paralog_graph._get_node(p1)
                paralog2 = ws_paralog_graph._get_node(p2)
                p1name = paralog1.name()
                p2name = paralog2.name()
                p1score = float(name2score[p1name])
                p2score = name2score[p2name]
                p1score_adj = p1score + similarity*p2score
                p2score_adj = p2score + similarity*p1score
                ws.loc[ws['Name'] == p1name, RANK] = p1score_adj
                ws.loc[ws['Name'] == p2name, RANK] = p2score_adj
                paralog_counter += 2
              else:
                print('Paralog relation has no similarity score!')
          self.report_pandas[ws_name] = ws.sortrows(by=[RANK,'Name'],ascending=[False,True],skip_rows=1)
          print(f'Adjusted scores for {paralog_counter} paralogs in {ws_name} worksheet')
        except KeyError:
            continue


    def make_report(self,add_tm_bibliography=False):
        start_time = time.time()
        self.set_input()
        self.find_targets()
        self.load_target_partners()

        self.find_drugs()
        self.find_inducers()
        
        if self.init_semantic_search():
            self.score_target_semantics()

        self.normalize_counts()
        self.paralog_adjustment()
        self.add_graph_bibliography()
        if add_tm_bibliography:
          self.add_bibliography4targets()
        self.add_target_annotation()

        disease_df = self.input_disease_df()
        if not disease_df.empty:
          self.add2report(disease_df)
    
        self.clear()
        print('Target ranking for %s was done in %s' % 
        (self._disease2str(), execution_time(start_time)))

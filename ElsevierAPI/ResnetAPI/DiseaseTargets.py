from .ResnetGraph import PROTEIN_TYPES,PHYSICAL_INTERACTIONS,ResnetGraph,unpack
from .PSPathway import PSPathway
from .FolderContent import FolderContent
from .NetworkxObjects import ACTIVATED,REPRESSED,UNKNOWN_STATE,OBJECT_TYPE,CONNECTIVITY,EFFECT,MECHANISM,REFCOUNT,PSRelation,PSObject
from .SemanticSearch import SemanticSearch,OQL,RANK,execution_time
from .ResnetAPISession import NO_REL_PROPERTIES,ONLY_REL_PROPERTIES,REFERENCE_IDENTIFIERS,BIBLIO_PROPERTIES,DBID
from ..Embio.PSnx2Neo4j import Cypher
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


class DiseaseTargets(SemanticSearch):
  pass
  def __init__(self,*args,**kwargs):
      """
      input:
        APIconfig - args[0]
        self.set_target_disease_state() needs references
      """
      my_kwargs = {
              #'disease':[],
              'what2retrieve':BIBLIO_PROPERTIES,
              'ent_props' : ['Name', 'Class'], #'Class' is for target partners retrieval
              'rel_props' : [EFFECT],
              #'data_dir' : '',
              'add_bibliography' : True,
              'strict_mode' : False,
              'target_types' : ['Protein'],
              #'pathway_folders' : [],
              #'pathways' : [],
              "max_childs" : 11,
              "add_closeness":True, # obsolete
              'propagate_target_state_in_model':True,
              #'add_regulators4':dict(),
              #'add_targets4':dict(),
              #'add_inhibitors4':dict()
              'skip':False,
              #'ontology_file' : ''
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
      self.GVs = list()
      self.targets4strictmode = set()
      self.ct_drugs = set() # compounds from clinical trials for input disease
      self.disease_inhibitors = set() # SmallMol inhibiting input disease
      self.disease_inducers = list() # SmallMol inducing or excerbating input disease
      
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

  
  def input_disease(self,propName:str):
    '''
    output:
      list of first values for each disease in self.input_diseases. Use it to return disease Name, URN, Connectivity etc..
    '''
    return [p.get_prop(propName) for p in self.input_diseases]
          

  def names4tmsearch(self):
    search_names = list(self.params['disease'].values())
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
      return ','.join(self.params['disease'].keys())+rep_pred+'targets'


  def _disease2str(self):
      # do not put disease name in quotas - it is used to create a file name
      disease_abbreviations = ','.join(self.params['disease'].keys())
      if len(self.input_diseases) > 1:
        disease_abbreviations += ' ('+str(len(self.input_diseases))+' types)'
      return disease_abbreviations


  def _targets(self)->set[PSObject]:
      return self.targets4strictmode if self.params['strict_mode'] else self.__targets__
  
  def target_types(self):
      return self.params['target_types']

  def _targets_dbids(self):
      return ResnetGraph.dbids(list(self._targets()))

  def __target_types_str(self):
      '''
        output: 
          ','.join(self.params['target_types'])
      '''
      return ','.join(self.params['target_types'])
  
  # STEP 1 in make_report:
  def set_input(self):
    disease_names = list(self.params['disease'].values())
    if self.useNeo4j():
      self.input_diseases = list(self.neo4j.get_nodes(objtype='Disease',propName='Name',
                            propVals=disease_names,with_childs=True,with_connectivity=True))
    else:
      my_session = self._clone_session()
      my_session.add_ent_props([CONNECTIVITY]) # Connectivity is needed for make_disease_df()
      ontology_graph = my_session.child_graph(disease_names,['Name','Alias'])
      if isinstance(ontology_graph, ResnetGraph):
        self.input_diseases = ontology_graph._psobjs_with('Disease','ObjTypeName')
        # filter by Disease is necessary because some proteins in monogenic disorders may have the same name as disease
        input_diseases_dbids = ResnetGraph.dbids(self.input_diseases)
        self.find_disease_oql = OQL.get_objects(input_diseases_dbids)
    return

  # STEP 2.1 in make_report()->find_targets()
  def targets_from_db(self): 
    if self.useNeo4j():
      disease_urns = ResnetGraph.urns(self.input_diseases)
      target_disease_graph = self.neo4j._neighborhood_(disease_urns,'URN',self.target_types())
      self.Graph = self.Graph.compose(target_disease_graph)
    else:
      req_name = f'Find targets linked to {self._disease2str()}'
      select_targets = f'SELECT Entity WHERE objectType = ({self.__target_types_str()})'
      OQLquery = f'SELECT Relation WHERE NeighborOf ({select_targets}) AND NeighborOf ({self.find_disease_oql})'      
      target_disease_graph = self.process_oql(OQLquery,req_name)
      assert(isinstance(target_disease_graph,ResnetGraph))

    return target_disease_graph
  
  # STEP 2.2 in make_report()->find_targets()
  def __GVtargets(self)->list[PSObject]:
    disease_names = self._disease2str()
    if self.useNeo4j():
      disease_urns = ResnetGraph.urns(self.input_diseases)
      self.gv2diseases = self.neo4j._connect_(['Disease'],disease_urns,'URN',['GeneticVariant'],[],'')
      self.Graph = self.Graph.compose(self.gv2diseases)
    else:
      oql_query = f'SELECT Relation WHERE objectType = FunctionalAssociation \
          AND NeighborOf ({self.find_disease_oql}) AND NeighborOf (SELECT Entity WHERE objectType = GeneticVariant)'
      request_name = f'Select GeneticVariants for {disease_names}'
      self.gv2diseases = self.process_oql(oql_query,request_name)
      assert isinstance(self.gv2diseases,ResnetGraph)

    self.GVs = self.gv2diseases.psobjs_with(only_with_values=['GeneticVariant'])
    target_withGVs = []
    if self.GVs:
      if self.useNeo4j():
        gv_urns = ResnetGraph.urns(self.GVs)
        cypher, params = Cypher.connect(['GeneticVariant'],gv_urns,'URN',self.target_types(),[],'',{OBJECT_TYPE:['GeneticChange']})
        gv2target = self.neo4j.fetch_graph(cypher,params)
      else:
        GVdbids = ResnetGraph.dbids(self.GVs)
        oql_query = f'SELECT Relation WHERE objectType = GeneticChange AND \
            NeighborOf (SELECT Entity WHERE objectType = ({self.__target_types_str()})) AND NeighborOf ({OQL.get_objects(GVdbids)})'
        request_name = f'Find targets with genetic variants linked to {disease_names}'
        # avoid adding GeneticChange to self.Graph
        new_session = self._clone(to_retrieve=NO_REL_PROPERTIES,init_refstat=False)
        gv2target = new_session.process_oql(oql_query,request_name)
        assert isinstance(gv2target, ResnetGraph)

      self.gv2diseases = self.gv2diseases.compose(gv2target)
      target_withGVs = gv2target.psobjs_with(only_with_values=PROTEIN_TYPES)
      gv2dis_refs = self.load_references(self.gv2diseases)

      print('Found %d targets with %d GeneticVariants for %s supported by %d references' % 
          (len(target_withGVs), len(self.GVs), disease_names, len(gv2dis_refs)))
        
    return target_withGVs

  # STEP #2.3 in make_report()->find_targets()
  def metabolite2inhibit(self):
    metabolites2inhibit = self.params.get('metabolites2inhibit',[])
    if metabolites2inhibit:
      if self.useNeo4j():
        disease_urns = self.input_disease('URN')
        metbolite2diseaseG = self.neo4j._connect_(['SmallMol'],self.params['metabolites2inhibit'],'Name',
                                              ['Disease'],disease_urns,'URN')
        [metbolite2diseaseG.set_edge_annotation(r.uid(),t.uid(),rel.urn(),EFFECT,['positive']) for r,t,rel in metbolite2diseaseG.iterate()]
        met_targets = self.neo4j.get_nodes('SmallMol', 'Name',metabolites2inhibit)
        newrels = []
        for met in met_targets:
          if met.uid() not in metbolite2diseaseG:
            for disease in self.input_diseases:
              newrels.append(PSRelation.make_rel(met,disease,{OBJECT_TYPE:'Regulation', EFFECT:'positive'}))
        metbolite2diseaseG = metbolite2diseaseG.compose(ResnetGraph.from_rels(newrels))
      else:
        metabolites2inhibit = OQL.join_with_quotes(self.params['metabolites2inhibit'])
        select_metabolite = f'SELECT Entity WHERE Name = ({metabolites2inhibit})'
        oql = f'SELECT Relation WHERE NeighborOf ({self.find_disease_oql}) AND NeighborOf ({select_metabolite})'
        my_session = self._clone_session()
        metbolite2diseaseG = my_session.process_oql(oql,'Find metabolites2inhibit')
        [metbolite2diseaseG.set_edge_annotation(r.uid(),t.uid(),rel.urn(),EFFECT,['positive']) for r,t,rel in metbolite2diseaseG.iterate()]

    
      if metbolite2diseaseG: # force EFFECT since metabolite must be inhibited
        # to ensure proper targets state in set_target_state:
        self.Graph.add_graph(metbolite2diseaseG) 
        metabolites2inhibit_objs = metbolite2diseaseG._psobjs_with('SmallMol',OBJECT_TYPE)
        return metabolites2inhibit_objs
      else:
        print(f'Cannot find metabolites {metabolites2inhibit} to inhibit for {self._disease2str()}')
        return []
    else:
      return []
    

  def cleavage2expression(self,target_names:list[str])->ResnetGraph:
    relProps = {OBJECT_TYPE:['ProtModification'],MECHANISM:['cleavage']}
    reqname = f'Find ProtModification targets 4 {target_names}'
    if self.useNeo4j():
      cypher,param = Cypher.expand(target_names,'Name',PROTEIN_TYPES,relProps,'upstream')
      proteasesG = self.neo4j.fetch_graph(cypher,param,reqname)
    else:
      oql_query = OQL.expand_entity(target_names,['Name'],relProps,PROTEIN_TYPES,'upstream')
      new_session = self._clone_session()
      proteasesG = new_session.process_oql(oql_query,reqname)
    
    if proteasesG:
      proteasesGrels = proteasesG._psrels()
      for rel in proteasesGrels:
        rel[OBJECT_TYPE] = ['Expression']
        rel[EFFECT] = ['negative']
        rel.urn(refresh=True)
      return ResnetGraph.from_rels(proteasesGrels)
    else:
      print(f'cleavage2expression for {target_names} did not return any data')
      return ResnetGraph()
    

  # STEP 2.5 in make_report()->find_targets()
  def expand_targets(self,reltype2targets:dict[str,list[str]],by_relProps:dict[str,list[str|int|float]]={},dir='upstream'):
    '''
    input:
      reltype2targets = {reltype:[targets]}
      by_relProps = {reltype:[propValue1,propValue2,...]},
      dir = upstream, downstream,
    output:
      updated self.expanded_targets
    '''
    neighbs = 'regulators' if dir=='upstream' else 'targets'
    for reltype,target_names in reltype2targets.items():
      relProps = {OBJECT_TYPE:[reltype]}
      relProps.update(by_relProps)
      reqname = f'Find {reltype} {neighbs} 4 {target_names}'
      if self.useNeo4j():
        if  reltype != 'ProtModification':
          expanded_targetsG = self.neo4j._neighborhood_(target_names,'Name',PROTEIN_TYPES,relProps,dir=dir)
        else:
          cypher,param = Cypher.expand(target_names,'Name',PROTEIN_TYPES,relProps,dir)
          return_pos = cypher.rfind('RETURN')      
          cypher = cypher[:return_pos] + "AND r.Mechanism IS NULL OR r.Mechanism <> 'cleavage'\n"+ cypher[return_pos:]
          expanded_targetsG = self.neo4j.fetch_graph(cypher,param,reqname)
      else:
        oql_query = OQL.expand_entity(target_names,['Name'],relProps,PROTEIN_TYPES,dir)
        if reltype == 'ProtModification': # ProtModification with Mechanism = cleavage is Expression regulation
          oql_query += f' AND NOT ({MECHANISM} = cleavage)'
        expanded_targetsG = self.process_oql(oql_query,reqname)

      if reltype == 'Expression': # ProtModification with Mechanism = cleavage is Expression regulation with Effect negative
        expanded_targetsG = expanded_targetsG.compose(self.cleavage2expression(target_names))

      if expanded_targetsG:
        expanded_targets = [n for n in self.__targets__ if n.name() in target_names]
        expansion_descr = f'{reltype} {neighbs} of {target_names}'
        expanded_targetsG.name = expansion_descr
        rel_rank = ['DirectRegulation','Binding','ProtModification','MolTransport','Regulation','PromoterBinding','Expression']
        # assume that Regulation means post-translational regulation
        expanded_targetsG = expanded_targetsG.make_simple(rel_rank)
        self.__targets__.update(expanded_targetsG._get_nodes())
        for t in expanded_targets:
          t_neighbors = expanded_targetsG.get_neighbors({t})
          neigh2target = dict({n.name():t.name() for n in t_neighbors})       
          self.expanded_targets.update(neigh2target)
        print(f'Added {len(neigh2target)} {expansion_descr}')
        self.Graph = self.Graph.compose(expanded_targetsG)
    return

  # STEP 2 in make_report()
  def find_targets(self):
    target_disease_graph = self.targets_from_db() # STEP 2.1

    if isinstance(target_disease_graph,ResnetGraph):
      self.__targets__ = set(target_disease_graph.psobjs_with(only_with_values=self.target_types()))
      target_refs = self.load_references(target_disease_graph)
      print('Found %d targets linked to %s supported by %d references in database' 
                  % (len(self.__targets__),self._disease2str(), len(target_refs)))
      
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
        self.expand_targets(self.params.get("add_regulators4",dict()),dir='upstream')
        self.expand_targets(self.params.get("add_targets4",dict()),dir='downstream')
        targets4inhibitors = self.params.get("add_inhibitors4",dict())
        self.expand_targets(targets4inhibitors,{EFFECT:['positive','negative']},dir='downstream')
        return


  def find_disease_inhibitors(self): #Step #4.1 in make_report()->find_drugs()
    if self.useNeo4j():
      disease_urns = ResnetGraph.urns(self.input_diseases)
      disease_inhibitorsG = self.neo4j._connect_(['SmallMol'],[],'',
                                          ['Disease'],disease_urns,'URN',
                                        {OBJECT_TYPE:['Regulation'],EFFECT:['negative']})
      self.Graph = self.Graph.compose(disease_inhibitorsG)
    else:
      request_name = 'Find compounds inhibiting {}'.format(self._disease2str())
      OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = negative \
          AND NeighborOf (SELECT Entity WHERE objectType = SmallMol) AND NeighborOf ({})'
      OQLquery = OQLquery.format(self.find_disease_oql)
      disease_inhibitorsG = self.process_oql(OQLquery,request_name)
      
    drugs_references = self.load_references(disease_inhibitorsG)
    print('Found %d compounds inhibiting %s supported by %d references' % 
                (len(disease_inhibitorsG),self._disease2str(), len(drugs_references)))
    return set(disease_inhibitorsG._psobjs_with('SmallMol','ObjTypeName'))


  def drugs_from_clinical_trials(self): #Step #4.2 in make_report()->find_drugs()
    if self.useNeo4j():
      disease_urns = ResnetGraph.urns(self.input_diseases)
      self.ct_graph = self.neo4j._connect_(['SmallMol'],[],'', ['Disease'],disease_urns,'URN',
                      {OBJECT_TYPE:['ClinicalTrial']},
                      dir=True
                      )
      self.Graph = self.Graph.compose(self.ct_graph)
    else:
      request_name = 'Find clinical trials for {}'.format(self._disease2str())
      OQLquery = 'SELECT Relation WHERE objectType = ClinicalTrial \
          AND NeighborOf (SELECT Entity WHERE objectType = SmallMol) AND NeighborOf ({})'
      OQLquery = OQLquery.format(self.find_disease_oql)
      self.ct_graph = self.process_oql(OQLquery,request_name)

    ct_refs = self.load_references(self.ct_graph)
    print('Found %d compounds on %d clinical trilas for %s' % (len(self.ct_drugs),len(ct_refs),self._disease2str()))
    return set(self.ct_graph._psobjs_with('SmallMol',OBJECT_TYPE))

  def find_drugs(self):  #Step #4 in make_report()
      self.disease_inhibitors = self.find_disease_inhibitors()
      self.disease_inhibitors.update(self.drugs_from_clinical_trials())


  def find_inducers(self): #Step 5 in make_report()
    if self.useNeo4j():
      disease_urns = ResnetGraph.urns(self.input_diseases)
      induction_graph = self.neo4j._connect_(['SmallMol'],[],'',
                                          ['Disease'],disease_urns,'URN',
                                        {OBJECT_TYPE:['Regulation'],EFFECT:['positive']})
      self.Graph = self.Graph.compose(induction_graph)
    else:
      REQUEST_NAME = 'Find compounds inducing/exacerbating {}'.format(self._disease2str())
      OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = positive \
          AND NeighborOf (SELECT Entity WHERE objectType = SmallMol) AND NeighborOf ({})'
      OQLquery = OQLquery.format(self.find_disease_oql)
      induction_graph = self.process_oql(OQLquery,REQUEST_NAME)
    
    inducers_refs = self.load_references(induction_graph)
    self.disease_inducers = induction_graph._psobjs_with('SmallMol','ObjTypeName')
    before_count = len(self.disease_inducers)
    # to remove indications reported as toxicities in clinical trials:
    self.disease_inducers ={x for x in self.disease_inducers if x not in self.ct_drugs}
    remove_count = before_count - len(self.disease_inducers)
    print(f'{remove_count} {self._disease2str()} inducers were removed because they are linked by clinical trial')
    print('Found %d compounds inducing %s supported by %d references' % 
        (len(self.disease_inducers),self._disease2str(),len(inducers_refs)))
    return


  # STEP #3.1 in make_report()->load_target_partners()
  def ligands4receptors(self):
    all_receptors = self.Graph._psobjs_with('Receptor','Class')
    added_receptors = self.partner2targets._psobjs_with('Ligand','Class')
    orphan_receptors = set(all_receptors).difference(set(added_receptors))
    if orphan_receptors: # find ligands for receptors among targets
      if self.useNeo4j():
        # finding protein ligands:
        receptor_urns = ResnetGraph.urns(orphan_receptors)
        ligands4receptorsG = self.neo4j._connect_(PROTEIN_TYPES,['Ligand'],'Class',
                                                  PROTEIN_TYPES,receptor_urns,'URN',
                      {OBJECT_TYPE:['DirectRegulation'],EFFECT:['positive'],REFCOUNT:'> 5'},
                      dir=True)
        
        # finding metabolite ligands:
        cypher, params = Cypher.match_psobjs(orphan_receptors,'t')
        cypher += Cypher.match_childs('mammal endogenous compounds and their derivative','m')
        cypher += 'MATCH (m)-[r]->(t)\n'
        cypher = Cypher.add_relProps(cypher,{OBJECT_TYPE:'DirectRegulation',EFFECT:'positive',REFCOUNT:'> 25'})
        cypher += 'RETURN m,r,t\n'
        metabolite2receptorsG = self.neo4j.fetch_graph(cypher,params)
      else:
        ref_cutoff = 5
        find_receptors_oql = OQL.get_objects(ResnetGraph.dbids(orphan_receptors))
        ligands_oql = f'SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive AND RelationNumberOfReferences > {ref_cutoff} AND \
            NeighborOf downstream (SELECT Entity WHERE Class = Ligand) AND NeighborOf upstream ({find_receptors_oql})'
        #opening new APIsession to avoid adding result graph to self.Graph
        new_session = self._clone(to_retrieve=NO_REL_PROPERTIES,init_refstat=False)
        request_name = f'Find ligands for {str(len(orphan_receptors))} receptors'
        ligands4receptorsG = new_session.process_oql(ligands_oql,request_name)
        
        find_metabolites = OQL.get_childs(['mammal endogenous compounds and their derivatives'],['Name'])
        metabolite_ligands_oql = f'SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive \
          AND NeighborOf downstream ({find_metabolites}) AND NeighborOf upstream ({find_receptors_oql}) AND RelationNumberOfReferences > 25'
        new_session = self._clone(to_retrieve=NO_REL_PROPERTIES,init_refstat=False)
        metabolite2receptorsG = new_session.process_oql(metabolite_ligands_oql,'Find metabolite ligands')
        new_session.close_connection()

      if ligands4receptorsG:
        new_ligands = ligands4receptorsG._psobjs_with('Ligand','Class')
        print('Found %d additional ligands for %d receptors' % (len(new_ligands),len(orphan_receptors)))

      if metabolite2receptorsG:
        metabolite_ligands = metabolite2receptorsG._psobjs_with('SmallMol','ObjTypeName')
        print('Found %d metabolite ligands for %d receptors' % (len(metabolite_ligands),len(orphan_receptors)))
      return ligands4receptorsG.compose(metabolite2receptorsG)
    else:
      print('No Receptor Class entities were found among %s targets' % self.target_types())
      return ResnetGraph()
      

  # STEP #3.2 in make_report()->load_target_partners()
  def receptors4ligands(self):
    ligands = self.Graph._psobjs_with('Ligand','Class')
    added_ligands = self.partner2targets._psobjs_with('Ligand','Class')
    orphan_ligands = set(ligands).difference(set(added_ligands))
    if orphan_ligands: # find receptors for ligands among targets
      if self.useNeo4j():
        ligand_urns = ResnetGraph.urns(ligands)
        orphan_ligand_partners = self.neo4j._connect_(PROTEIN_TYPES,ligand_urns,'URN',
                                                  PROTEIN_TYPES,['Receptor'],'Class',         
                      {OBJECT_TYPE:['DirectRegulation'],EFFECT:['positive'],REFCOUNT:'> 5'},
                      dir=True)
      else:
        orphan_ligands_dbids = ResnetGraph.dbids(orphan_ligands)
        find_orphan_ligands = OQL.get_objects(list(orphan_ligands_dbids))
        receptor_oql = f'SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive AND RelationNumberOfReferences > 5 \
            AND NeighborOf upstream (SELECT Entity WHERE Class = Receptor) AND NeighborOf downstream ({find_orphan_ligands})'
        request_name = f'Find receptors for {str(len(orphan_ligands_dbids))} orphan ligands'
        #opening new APIsession to avoid adding result graph to self.Graph
        new_session = self._clone(to_retrieve=NO_REL_PROPERTIES,init_refstat=False)
        orphan_ligand_partners = new_session.process_oql(receptor_oql,request_name)
        new_session.close_connection()

      new_receptors = orphan_ligand_partners._psobjs_with('Receptor','Class')
      print('Found %d receptors for %d orphan ligands' % (len(new_receptors),len(orphan_ligands_dbids)))
      return orphan_ligand_partners
      
    # STEP #3.3 in make_report()->load_target_partners()
  def metabolite_products(self):
    if self.useNeo4j():
      # finding metabolite products:
      target_urns = ResnetGraph.urns(self._targets())
      relProps = {OBJECT_TYPE:['MolTransport','MolSynthesis','ChemicalReaction'],EFFECT:['positive'],REFCOUNT:'> 50'}
      metbolite_productsG = self.neo4j._connect_(PROTEIN_TYPES,target_urns,'URN',
                                      ['SmallMol'],[],'',relProps,dir=True)
    else:
      find_metabolite_products = 'SELECT Relation WHERE objectType = (MolTransport, MolSynthesis,ChemicalReaction) \
        AND Effect = positive AND NeighborOf downstream (SELECT Entity WHERE id =({ids})) AND \
        NeighborOf upstream (SELECT Entity WHERE objectType = SmallMol) AND RelationNumberOfReferences > 50'
      req_name = f'Finding downstream metabolite products of targets for {self._disease2str()}'
      new_session = self._clone(to_retrieve=NO_REL_PROPERTIES,init_refstat=False)
      metbolite_productsG = new_session.iterate_oql(find_metabolite_products,set(self._targets_dbids()),request_name=req_name)
      #self.partner2targets.add_graph(metbolite_productsG)
      new_session.close_connection()
    metabolite_products = metbolite_productsG._psobjs_with('SmallMol','ObjTypeName')
    print('Found %d metabolites produced by %d targets of %s' % (len(metabolite_products),len(self.__targets__),self._disease2str()))
    return metbolite_productsG


  # STEP #3 in make_report()
  def load_target_partners(self):
    self.partner2targets.add_graph(self.ligands4receptors())
    self.partner2targets.add_graph(self.receptors4ligands())
    self.partner2targets.add_graph(self.metabolite_products())

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
    for folder_name in self.params.get('pathway_folders',[]):
      ps_pathways = fc.folder2pspathways(folder_name,with_layout=False)
      disease_pathways += ps_pathways

    if self.params.get('pathways',[]):
      pathway_names_lower = list(map(str.lower,self.params['pathways']))
      filtered_pathway_list = list()
      for ps_pathway in disease_pathways:
        assert(isinstance(ps_pathway, PSPathway))
        if ps_pathway.name().lower() in pathway_names_lower:
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

      
  def init_semantic_search(self): # STEP 6 in make_report()
    print('\n\nInitializing semantic search')
    if  self.__targets__:
      targets4ranking = list(self._targets())
      if self.useNeo4j():
        self.RefCountPandas = self.load_df_neo4j(targets4ranking,max_childs=self.params['max_childs'])
      else:
        self.RefCountPandas = self.load_df(targets4ranking,max_childs=self.params['max_childs'],max_threads=10)
        self.entProps = ['Name'] # Class is no longer needed
      self.RefCountPandas.add_entities(targets4ranking)
      self.RefCountPandas._name_ = 'DiseaseTargets'
      print(f'Will score {len(self.RefCountPandas)} targets linked to {self._disease2str()}')       
      return True
    else:
      print ('No targets found for %s' % self._disease2str())
      print('Consider setting "strict_mode" to False')
      return False


  def score_GVs(self): # STEP 7a in make_report()->score_target_semantics()
      print(f'\n\nScoring targets by number of semantic references linking their Genetic Variants to {self._disease2str()}')
      if not self.GVs:
        print(f'{self._disease2str()} has no known Genetic Variants')
        return

      concept_name = f'GVs for {self._disease2str()}'
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
      id_type = 'URN' if self.useNeo4j() else DBID
      for i in self.RefCountPandas.index:
        target_ids = list(self.RefCountPandas.at[i,self.__temp_id_col__])
        row_targets = self.Graph.psobj_with_ids(set(target_ids),id_type)
        targetGVs = self.gv2diseases.get_neighbors(set(row_targets), allowed_neigbors=self.GVs)
            
        if targetGVs:
          GVnames = set([n.name() for n in targetGVs])
          self.RefCountPandas.at[i,self.__colnameGV__] = ';'.join(GVnames)
          target_gvlinkcounter += 1
          gv_disease_subgraph = self.gv2diseases.get_subgraph(list(targetGVs),self.input_diseases)
          GVrefcount = len(self.load_references(gv_disease_subgraph))
        else:
          GVrefcount = 0
        
        GVscore = GVrefcount*(1.0 + len(targetGVs)/len(self.GVs))
        self.RefCountPandas.at[i,weighted_refcount_column] = GVscore
        self.RefCountPandas.at[i,refcount_column] = GVrefcount
        self.RefCountPandas.at[i,linked_count_column] = len(targetGVs)

      if target_gvlinkcounter:
        self.set_rank4(concept_name,self.RefCountPandas)
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
        # attempting to deduce effect from GeneticChange:
      elif partners2disease_graph.psrels_with(['GeneticChange']):
        return REPRESSED
      else:
        return UNKNOWN_STATE


  def set_target_disease_state(self): # STEP 7b in make_report()->score_target_semantics()
    print('\n\nCalculating targets state (activated/repressed) in %s' % self._disease2str())
    targets = set()
    self.RefCountPandas['State in Disease'] = UNKNOWN_STATE
    id_type = self.idtype()
    for i in self.RefCountPandas.index:
        target_dbids = list(self.RefCountPandas.at[i,self.__temp_id_col__])
        row_targets = self.Graph.psobj_with_ids(set(target_dbids),id_type)
        state = self.__vote4effect(row_targets)
        self.RefCountPandas.at[i,'State in Disease'] = int(state)
        if state != UNKNOWN_STATE:
          [t.set_state(state) for t in row_targets]
          targets.update(row_targets)

    self.disease_model = ResnetGraph()
    for graph in (p.graph for p in self.disease_pathways.values()):
      self.disease_model = self.disease_model.compose(graph)

    if self.disease_model and 'propagate_target_state_in_model' in self.params:
      uid2state = dict()
      [uid2state.update(self.disease_model.propagate_state(t)) for t in targets if t.uid() in self.disease_model.nodes()]
      #id_type = 'URN' if self.useNeo4j() else DBID

      def __disease_state_from_model():
        activ_targets_pd = self.RefCountPandas.loc[(self.RefCountPandas['State in Disease'] == ACTIVATED)]
        activ_targets_ids = unpack(activ_targets_pd[self.__temp_id_col__].to_list())
        inhib_targets_pd = self.RefCountPandas.loc[(self.RefCountPandas['State in Disease'] == REPRESSED)]
        inhib_targets_ids = unpack(inhib_targets_pd[self.__temp_id_col__].to_list())

        activated_targets = self.disease_model.psobj_with_ids(activ_targets_ids,id_type)
        inhibited_trgts = self.disease_model.psobj_with_ids(inhib_targets_ids,id_type)

        new_regulators_count = 0
        for idx in self.RefCountPandas.index:
          if self.RefCountPandas.at[idx,'State in Disease'] == UNKNOWN_STATE:
            net_effect = 0
            for target_id in self.RefCountPandas.at[idx,self.__temp_id_col__]:
              if target_id in model_id2uid:
                target_uid = model_id2uid[target_id]
                target = self.disease_model._get_node(target_uid)
                activates_activated_targets,inhibits_activated_targets = self.disease_model.net_regulator_effect(target,activated_targets)
                activates_inhibited_targets,inhibits_inhibited_targets = self.disease_model.net_regulator_effect(target,inhibited_trgts)
                net_effect += len(activates_activated_targets)+len(inhibits_inhibited_targets)-len(inhibits_activated_targets)-len(activates_inhibited_targets)

            target_state = ACTIVATED if net_effect > 0 else REPRESSED if net_effect < 0 else UNKNOWN_STATE
            self.RefCountPandas.at[idx,'State in Disease'] = target_state
            new_regulators_count += 1             
        return new_regulators_count
      
      for i in self.RefCountPandas.index:
        if self.RefCountPandas.at[i,'State in Disease'] == UNKNOWN_STATE:
          target_dbids = list(self.RefCountPandas.at[i,self.__temp_id_col__])
          row_targets = self.Graph.psobj_with_ids(set(target_dbids), id_type)
          row_state = UNKNOWN_STATE
          for t in row_targets:
            tuid = t.uid()
            if tuid in uid2state: 
              row_state += uid2state[t.uid()]
          self.RefCountPandas.at[i,'State in Disease'] = row_state

        # determining disease state of the targets from disease model pathways
        # by finding their regulatory effect on targets linked to disease
        model_id2uid = self.disease_model.dbid2uid(id_type=self.idtype())          
        while __disease_state_from_model(): pass

    if self.expanded_targets:
      for idx in self.RefCountPandas.index:
        if self.RefCountPandas.at[idx,'State in Disease'] == UNKNOWN_STATE:
          target_name = self.RefCountPandas.at[idx,'Name']
          if target_name in self.expanded_targets:      
            expanded_name = self.expanded_targets[target_name]
            expanded_state = int(self.RefCountPandas.loc[(self.RefCountPandas['Name'] == expanded_name)]['State in Disease'].iloc[0])
            expanded_target = [n for n in self.__targets__ if n.name() == expanded_name]
            row_targets = [n for n in self.__targets__ if n.name() == target_name]
            my_subgraph = self.Graph.neighborhood(set(expanded_target),row_targets)
            my_subgraph.clone_node(expanded_target[0],self.input_diseases[0],expanded_state,'Regulation')
            self.Graph = my_subgraph.compose(self.Graph)
            state = self.__vote4effect(row_targets)
            self.RefCountPandas.at[idx,'State in Disease'] = int(state)
    return


  def score_partners(self):
    print(f'\n\nFinding semantic references supporting links between target partners and {self._disease2str()}')
    my_df = self.RefCountPandas.dfcopy()
    id_type = self.idtype()
    partners_ids = []
    partner_names = []
    target_names = []
    partners = set()
    for i in my_df.index:
      target_ids = list(my_df.at[i,self.__temp_id_col__])
      row_targets = self.Graph.psobj_with_ids(set(target_ids),id_type)
      row_partners = set(self.partner2targets.get_neighbors(set(row_targets)))
      partner_names.append(','.join([n.name() for n in row_partners]))
      target_names.append(str(my_df.at[i,'Name']))

      if row_partners:
        partners_ids.append(tuple(ResnetGraph.ids(row_partners,id_type)))
        partners.update(row_partners)
      else:
        partners_ids.append(tuple([0])) # adding fake id for link2concept to work
      
    partners_df = df.from_dict({'Name':target_names,self.__temp_id_col__:partners_ids,'Target partners':partner_names})
    partners_df.add_entities(partners)

    how2connect = self.set_how2connect(**dict())
    concept_name = 'target partners'
    connectionG,partners_df = self.link2concept(concept_name,self.input_diseases,partners_df,how2connect)
    linked_partners_count = len(connectionG._get_nodes(my_df.uids()))
    print(f'{len(partners_df)} targets have {linked_partners_count} partners linked to {self._disease2str()}')

    if linked_partners_count:
      partners_df = df.from_pd(partners_df.drop(columns=[self.__temp_id_col__]))
      updated_RefCountPandas = self.RefCountPandas.merge_df(partners_df,on='Name')
      self.set_rank4(concept_name,updated_RefCountPandas)
      return updated_RefCountPandas
    else:
      return self.RefCountPandas


  def make_disease_network(self): # STEP 8.1 in make_report()->DiseaseNetworkRegulationScore()
    print(f'Creating physical interaction network between targets of {self.input_disease('Name')} to calculate target closeness for "score_regulators" function')
    if self.useNeo4j():
      disease_network = self.neo4j.get_ppi(self.__targets__, minref=2,with_references=False)
    else:
      my_session = self._clone_session(to_retrieve=NO_REL_PROPERTIES,init_refstat=False)
      disease_network = my_session.get_ppi(self.__targets__, self.params.get('ppiRNEFs',[]), minref=2)
    
    disease_network.name = f'{self._disease2str()} PPPI network'
    return disease_network.make_simple(['DirectRegulation','ProtModification','Binding']) 


  def DiseaseNetworkRegulationScore(self): #` STEP 7e in make_report()->score_target_semantics()
      '''
      adds:
        PATHWAY_REGULATOR_SCORE column to self.RefCountPandas
      '''
      print('\n\nScoring regulators by distance to components of disease pathways',flush=True)
      print('Retrieving regulatory network between targets and components of disease pathway')
      regulation_graph = self.disease_model if self.disease_model else self.make_disease_network()
      uid2closeness = regulation_graph.closeness() # will use closeness in disease network to initialize regulatory ranking
      
      disease_pathway_components = set(regulation_graph._get_nodes())
      connected_nodes = disease_pathway_components
      unconnected_nodes = self.RefCountPandas.entities()

      for _ in range (0,5): # 5 degree of separation
        if self.useNeo4j():
          graph_expansion = self.connect_entities_neo4j(unconnected_nodes, connected_nodes,
                                                PHYSICAL_INTERACTIONS, in_direction='>')
        else:
          new_session = self._clone(to_retrieve=NO_REL_PROPERTIES,init_refstat=False) # new_session.Graph does not contain regulation_graph 
          new_session.Graph = regulation_graph.compose(new_session.Graph) # to speadup connect_nodes
          unconnected_nodes_ids = ResnetGraph.dbids(unconnected_nodes)
          connected_nodes_ids = ResnetGraph.dbids(connected_nodes)
          graph_expansion = new_session.connect_nodes(unconnected_nodes_ids, connected_nodes_ids,
                                                PHYSICAL_INTERACTIONS, in_direction='>')
        
        regulation_graph = regulation_graph.compose(graph_expansion)
        regulation_graph_nodes = set(regulation_graph._get_nodes())
        connected_nodes = regulation_graph_nodes - disease_pathway_components
        connected_nodes = connected_nodes.intersection(unconnected_nodes)
        # only nodes connected at the previous cycle need to be expanded at the next cycle
        if not connected_nodes:
          break
        unconnected_nodes = unconnected_nodes - regulation_graph_nodes

      print('Ranking targets by regulation score')
      if self.disease_pathways:
        for pathway in self.disease_pathways.values():
          assert(isinstance(pathway,PSPathway))
          uid2clos = pathway.graph.closeness()
          #regulation_graph.rank_regulatorsOLD(uid2clos,PATHWAY_REGULATOR_SCORE)
          regulation_graph.rank_regulators(uid2clos,PATHWAY_REGULATOR_SCORE)
      else:
        #regulation_graph.rank_regulatorsOLD(uid2closeness,PATHWAY_REGULATOR_SCORE)
        regulation_graph.rank_regulators(uid2closeness,PATHWAY_REGULATOR_SCORE)

      dbid2uid = regulation_graph.dbid2uid(id_type=self.idtype())
      for i in self.RefCountPandas.index:
          target_ids = list(self.RefCountPandas.at[i,self.__temp_id_col__])
          target_regulator_score = 0.0
          for target_id in target_ids:
              try:
                  target_uid = dbid2uid[target_id]
                  regulator_score = regulation_graph.nodes[target_uid][PATHWAY_REGULATOR_SCORE]
                  target_regulator_score += regulator_score
              except KeyError:
                  continue
          self.RefCountPandas.at[i,PATHWAY_REGULATOR_SCORE] = target_regulator_score

      self.RefCountPandas.set_rank(PATHWAY_REGULATOR_SCORE)
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


  def score_target_semantics(self): # STEP 7 in make_report()
    ###### do not multithread.  self.Graph will leak !!!! #####
    disease_str = self._disease2str()
    self.score_GVs() # STEP 7a
    self.set_target_disease_state() # STEP 7b

    kwargs = {'connect_by_rels':['Regulation'],
              'boost_with_reltypes':['FunctionalAssociation'],
              'concept_name':'Regulate '+ disease_str
              }
    _,self.RefCountPandas = self.score_concepts(self.input_diseases,**kwargs)

    kwargs = {'connect_by_rels':['GeneticChange'],
              'boost_with_reltypes':['FunctionalAssociation'],
              'concept_name':'Genetically linked to '+ disease_str}
    _,self.RefCountPandas = self.score_concepts(self.input_diseases,**kwargs)

    kwargs = {'connect_by_rels':['Biomarker'],
              'boost_with_reltypes':['FunctionalAssociation'],
              'concept_name':'Target is Biomarker in '+disease_str}
    _,self.RefCountPandas = self.score_concepts(self.input_diseases,**kwargs)

    kwargs = {'connect_by_rels':['QuantitativeChange'],
              'boost_with_reltypes':['FunctionalAssociation'],
              'concept_name':'Quantitatively changed in '+ disease_str}
    _,self.RefCountPandas = self.score_concepts(self.input_diseases,**kwargs)

    kwargs = {'connect_by_rels':['Regulation','QuantitativeChange','StateChange','Biomarker'],
              'boost_with_reltypes':['FunctionalAssociation','CellExpresion'],
              'concept_name':'symptoms for '+ disease_str,
              'clone2retrieve' : REFERENCE_IDENTIFIERS}
    _,self.RefCountPandas,_ = self.score_concept('symptoms',**kwargs)

    kwargs = {'connect_by_rels':['Regulation'],
              'boost_with_reltypes':['FunctionalAssociation'],
              'concept_name':'Clinical parameters for '+ disease_str,
              'clone2retrieve' : REFERENCE_IDENTIFIERS}
    _,self.RefCountPandas,_ = self.score_concept('clinical_parameters',**kwargs)

    self.combine_processes()

    kwargs = {'connect_by_rels':['Regulation'],
              'boost_with_reltypes':['FunctionalAssociation'],
              'concept_name':'Cell processes affected by '+ disease_str,
              'clone2retrieve' : REFERENCE_IDENTIFIERS}
    _,self.RefCountPandas,_ = self.score_concept('processes',**kwargs)

    self.RefCountPandas = self.score_partners() #STEP 7d
    self.DiseaseNetworkRegulationScore() # STEP 7e

    ################## SPLITTING RefCountPandas to score agonists and antagonists differently #####################
    activ_targets_df = df.from_pd(self.RefCountPandas.loc[(self.RefCountPandas['State in Disease'] > UNKNOWN_STATE)])
    activ_targets_df.copy_format(self.RefCountPandas)
    activ_targets_df.add_entities(self.RefCountPandas.entities(include_children=False))
    inhib_targets_df = df.from_pd(self.RefCountPandas.loc[(self.RefCountPandas['State in Disease'] < UNKNOWN_STATE)])
    inhib_targets_df.copy_format(self.RefCountPandas)
    inhib_targets_df.add_entities(self.RefCountPandas.entities(include_children=False))
    unk_targets_df = df.from_pd(self.RefCountPandas.loc[(self.RefCountPandas['State in Disease'] == UNKNOWN_STATE)])
    unk_targets_df.copy_format(self.RefCountPandas)
    unk_targets_df.add_entities(self.RefCountPandas.entities(include_children=False))
    print(f'Created 3 worksheets with: {len(activ_targets_df)} activated targets, {len(inhib_targets_df)} inhibited_targets, {len(unk_targets_df)} unknown state targets')
    
    kwargs = {'connect_by_rels':['DirectRegulation'],
              'with_effects':['negative'], 
              'boost_with_reltypes':['Binding','Expression','Regulation'],
              'concept_name':'drugs for '+disease_str,
              'clone2retrieve' : REFERENCE_IDENTIFIERS}
    connectionG,activ_targets_df = self.score_concepts(self.disease_inhibitors,activ_targets_df,**kwargs)
    activ_targets_uids = activ_targets_df.uids()
    linked_row_count = len(connectionG._get_nodes(activ_targets_uids))
    kwargs['with_effects'] = ['positive']
    if linked_row_count:
      kwargs['column_rank'] = activ_targets_df.max_colrank()
    kwargs['concept_name'] = 'inducers of '+disease_str
    connectionG,activ_targets_df = self.score_concepts(self.disease_inducers,activ_targets_df,**kwargs) 
    linked_row_count += len(connectionG._get_nodes(activ_targets_uids))

    # linking compounds inhibiting disease to targets inhibited in disease by effect positive:
    # measures synergy between disease drugs and inhibited targets
    kwargs['concept_name'] = 'drugs for '+disease_str
    connectionG,inhib_targets_df = self.score_concepts(self.disease_inhibitors,inhib_targets_df,**kwargs)
    linked_row_count = len(connectionG._get_nodes(inhib_targets_df.uids()))
    # linking compounds inducing disease to targets inhibited in disease by effect negative:
    # measures synergy between disease activators and inhibited targets
    if linked_row_count:
      kwargs['column_rank'] = inhib_targets_df.max_colrank()
    kwargs['concept_name'] = 'inducers for ' + disease_str
    kwargs['with_effects'] = ['negative']
    _,inhib_targets_df = self.score_concepts(self.disease_inducers,inhib_targets_df,**kwargs)

    all_compounds = self.disease_inducers|self.disease_inhibitors
    if all_compounds:
      print (f'\n\nLinking compounds modulating "{disease_str}" to targets with unknown state in "{self._disease2str()}"')
      kwargs = {'connect_by_rels':['DirectRegulation'],
                'boost_with_reltypes':['Binding','Expression','Regulation'],
                'concept_name':'compounds modulating '+disease_str,
                'clone2retrieve' : REFERENCE_IDENTIFIERS,
                'column_rank': unk_targets_df.max_colrank()}
      _,unk_targets_df = self.score_concepts(self.disease_inducers,unk_targets_df,**kwargs)

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
          rows.append((d.name(),str(d[CONNECTIVITY][0]),child_names_str,d.urn()))

      if len(rows) > 1:
        input_df = df.from_rows(rows,header=['Name',CONNECTIVITY,'Children','URN'])
        input_df[CONNECTIVITY] = input_df[CONNECTIVITY].astype(int)
        input_df = input_df.sortrows(by=CONNECTIVITY)
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
    self.set_input() #Step 1
    self.find_targets() #Step 2
    self.load_target_partners() #Step 3

    self.find_drugs() #Step 4
    self.find_inducers() #Step 5
    
    if self.init_semantic_search(): #Step 6
      self.score_target_semantics() #Step 7

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
    return

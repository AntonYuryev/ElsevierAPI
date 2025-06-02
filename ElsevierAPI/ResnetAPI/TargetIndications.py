from .SemanticSearch import SemanticSearch,len,df,pd,TOTAL_REFCOUNT,RELEVANT_CONCEPTS
from .ResnetAPISession import BIBLIO_PROPERTIES,NO_REL_PROPERTIES,SNIPPET_PROPERTIES,REFERENCE_IDENTIFIERS
from .ResnetGraph import REFCOUNT, PSObject, ResnetGraph, OBJECT_TYPE, EFFECT
from ..ETM_API.references import PS_SENTENCE_PROPS,PS_BIBLIO_PROPS,SENTENCE
from .PathwayStudioGOQL import OQL
from .FolderContent import FolderContent,PSPathway
from concurrent.futures import ThreadPoolExecutor,as_completed
from ..utils import execution_time
import time,os


RANK_SUGGESTED_INDICATIONS = True 
# in RANK_SUGGESTED_INDICATIONS mode only indications suggested in the lietarure are ranked by amount of supporting evidence
PREDICT_RANK_INDICATIONS = False
# mode also predicts and ranks indications from diseases having input target as biomarker

RAWDF4ANTAGONISTS = 'antagonistInd'
RAWDF4AGONISTS = 'agonistInd'

ANTAGONISTSDF = 'indications4antagonists'
AGONISTSDF = 'indications4agonists'
UNKEFFECTDF = 'possibilities'
TOXICITIES = 'TargetToxicities'

ANTAGONIST = -1 # drug inhibits its targets
AGONIST = 1 # drug activates its targets
ANY_MOA = 0

class Indications4targets(SemanticSearch):
    pass
    max_threads4ontology = 50
    max_cell_indications = 1000 # to limit number indications for popular ligands like TNF
  
    def __init__(self, *args, **kwargs):
        '''
        input:
          APIconfig - args[0], if empty uses default APIconfig
        '''
        my_kwargs = {'partner_names':[],
    # if partner_names is empty script will try finding Ligands for Receptor targets and Receptors for Ligand targets
                 'partner_class':'', # use it only if partner_names not empty
                 'indication_types': ['Disease','Virus','Pathogen'], #['Disease','Virus','CellProcess']
                 'target_names':[],
                 'pathway_name_must_include_target':True,
    # if 'pathway_name_must_include_target' True only pathways depicting target signaling are considered for regulome construction
    # if 'pathway_name_must_include_target' False all pathways containing both Targets and Partners will be considered
                 'strict_mode':RANK_SUGGESTED_INDICATIONS,
    # if strict mode algorithm only ranks indications suggetsed in the literarure without predicting new indications
    # if not in strict mode algorithm also predicts and rank additional indications based on target expression or genetic profiles in disease
                 'data_dir':'',
                 'add_bibliography' : True,
                 'what2retrieve':BIBLIO_PROPERTIES, # need biblio props to create PSbibliography
                 'mode_of_action':ANY_MOA,
                 'max_ontology_parent': 10
                }
        my_kwargs.update(kwargs)

        super().__init__(*args,**my_kwargs)
        self.add_rel_props([EFFECT])
        self.columns2drop += [self.__resnet_name__,self.__mapped_by__]
        # 4 threads perform a bit faster than 8 threads 
        # 288 disease parents out of 288 were processed in 0:19:51.732684 by 8 threads
        # 288 disease parents out of 288 were processed in 0:18:52.715937 by 4 threads

        # sets of PSObject:
        self.__targets__ = set()
        self.__partners__ = set()
        self.__targets__secretingCells = set()
        self.__GVs__ = list()

        self.__indications4antagonists__ = set()
        self.__indications4agonists__ = set()
        self.__unknown_effect_indications__ = set()

        self.__IndirectAgonists = set() # name reflects action on target
        self.__DirectAgonists = set() # activate input targets
        self.__IndirectAntagonists = set() # inhibit input targets
        self.__DirectAntagonists = set()

        self.targets_have_weights = False
        self.__coocG__ = ResnetGraph() # to add coocurence count
        self.ws_prefix = ''

############################## UTILITIES ############################ UTILITIES #############################   
    def input_target_names_str(self):
      name_str = ','.join(self.params['target_names'][:2])
      if len(self.params['target_names']) > 2:
          name_str += '...'         
      return name_str
    

    def target_names(self):
      """
      Returns
      -------
      entity names that can be used both in PS and ETM
      """
      #input names are used for finding references in ETM.
      # RepurposeDrug has it own 'target_names'.
      return [x['Name'][0] for x in self.__targets__]


    def target_names_str(self):
        """
        output:
            comma-separted string with self.target names as they appear in database.
            return may be different from self.params['target_names']
        """
        return ','.join([t.name() for t in self.__targets__]) if len(self.__targets__) < 3 else 'targets'
    
    def ws_names(self):
      '''
      output:
        Indication_ws = 'ind4'+ self.ws_prefix[0:21]+antag/agon
        Toxicities_ws = 'tox4'+ self.ws_prefix[0:21]+antag/agon
      '''
      ind_prefix = 'ind4'
      tox_prefix = 'tox4'
      trunc_name = self.ws_prefix[0:21]
      sep = '-' if trunc_name[-1].isalnum() else ''
      suffix = 'antag' if self.__moa() == ANTAGONIST else 'agon'
      return ind_prefix + trunc_name+sep+suffix, tox_prefix +trunc_name+sep+suffix
       

    def __known_effect_indications(self):
        known_effect_indication = self.__indications4agonists__ | self.__indications4antagonists__
        return known_effect_indication
    

    def __moa(self):
        return self.params['mode_of_action']
    

    def clear_indications(self):
      self.__targets__.clear()
      self.__indications4agonists__.clear()
      self.__indications4antagonists__.clear()
      self.__unknown_effect_indications__.clear()
      return


    def _my_indications(self):
        """
        Returns
        -------
        all found indications depending on self.params['strict_mode']
        """
        return self.__indications4antagonists__ | self.__indications4agonists__


    def _is_strict(self):
        return self.params['strict_mode']
    

    def oql4indications(self):
        '''
        output:
          GOQL query, indication_types_str
        '''
        indication_types_str = ','.join(self.params['indication_types'])
        return f'SELECT Entity WHERE objectType = ({indication_types_str})', indication_types_str
  

    def __load_targets(self)->set[PSObject]:
        assert(self.params['target_names'])
        i_t_n = self.params['target_names']

        self.add_ent_props(['Class'])
        oql = OQL.get_childs(i_t_n,['Name','Alias'],include_parents=True)
        targetsG = self.process_oql(oql,f'Find {i_t_n}')
        self.entProps.remove('Class')
        return set(targetsG._get_nodes())
    

    def _load_children4targets(self,targets:list[PSObject]):
        my_target_names = ResnetGraph.names(targets)
        oql = OQL.get_childs(my_target_names,['Name','Alias'],include_parents=True)
        twc = self.process_oql(oql,'Loading children4targets') # twc - target_with_children
        return twc._get_nodes()


    def _set_targets(self):
        '''
        input:
            self.params['target_names']
        output:
            self.__targets__
        '''
        self.__targets__ = self.__load_targets()
        if self.__targets__:
            self.__targets__.update(self._load_children4targets(self.__targets__))
            self.oql4targets = OQL.get_objects(ResnetGraph.dbids(self.__targets__))
        else:
            print(f'No targets were found for {self.params['target_names']}')


    @staticmethod
    def _partner_class(target_class:str):
        if target_class == 'Ligand': return "Receptor"
        elif target_class == 'Receptor': return "Ligand"
        else: return ''


    def find_partners(self,_4targets:list[PSObject])->tuple[set[PSObject],ResnetGraph]:
      """
      Finds
      -----
      Ligands if targets are Receptor, Receptors if targets are Ligand linked to target(s) with Effect=positive\n
      dose not work for any other target class
      output:
          self.partnets
      """
      allowed_target_clases = ['Ligand','Receptor']
      partners_graph = ResnetGraph()
      my_session = self._clone_session(to_retrieve = SNIPPET_PROPERTIES)
      for clas in allowed_target_clases:
        class_targets = [o for o in _4targets if o.get_prop('Class') == clas]
        if class_targets:
            partner_class = self._partner_class(clas)
            select_targets = OQL.get_objects(ResnetGraph.dbids(class_targets))
            target_names = OQL.get_objects(ResnetGraph.names(class_targets))

            REQUEST_NAME = f'Find {partner_class}s for {target_names}'
            #SELECTpartners = f'SELECT Entity WHERE Class = {partner_class} AND Connected by (SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive) to ({select_targets})'
            SELECTpartners = f'SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive AND NeigborOf (SELECT Entity WHERE Class = {partner_class}) AND NeigborOf ({select_targets})'
            partners_graph = partners_graph.compose(my_session.process_oql(SELECTpartners,REQUEST_NAME))
            assert(isinstance(partners_graph,ResnetGraph))
            if partner_class == 'Ligand':
                # request to find additional secreted molecules not annotated with Class=Ligand
                #SELECTsecretedpartners = 'SELECT Entity WHERE "Cell Localization" = Secreted AND objectType = Protein AND Connected by (SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive) to ({select_target})'
                SELECTsecretedpartners = f'SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive AND NeigborOf (SELECT Entity WHERE "Cell Localization" = Secreted AND objectType = Protein) AND NeigborOf ({select_targets})'
                secreted_partners = my_session.process_oql(SELECTsecretedpartners,f'Find secreted partners for {target_names}')
                assert isinstance(secreted_partners,ResnetGraph)
                partners_graph = secreted_partners.compose(partners_graph)

      partners = set(partners_graph._get_nodes()).difference(_4targets)
      partners_names = ','.join(ResnetGraph.names(partners))
      print(f'Found {partners_names} {clas}s as partners for {target_names}')
      return partners, partners_graph
    

    def params2partners(self, partner_params:dict,_4targets:list[PSObject])->tuple[set[PSObject],ResnetGraph]:
      '''
      case when self.params has 'partner_names' and 'partner_class' explicitly\n
      use it when receptor partners are metabolites and therefore cannot be found by select_partners()\n
      input:
          partner_params - {target_name:partner_class:[partner_names]}
      '''
      if not _4targets: return set(),ResnetGraph()
      _4target_names = {n.name():n for n in _4targets}
      partner_names = dict() #{name:(class,weight)}
      t2pG = ResnetGraph()
      my_session = self._clone_session(to_retrieve = SNIPPET_PROPERTIES)
      for target_name, c2p in partner_params.items():
        if target_name in _4target_names:
          assert(isinstance(c2p,dict))
          for p_class, p_names in c2p.items():
              for p_name in p_names:
                partner_weight = _4target_names[target_name].get_prop('target weight',if_missing_return=0)
                partner_names[str(p_name).lower()] = (p_class,partner_weight)
                target_partner_oql = f"SELECT Relation WHERE NeighborOf (SELECT Entity WHERE Name = '{target_name}') AND  NeighborOf (SELECT Entity WHERE Name = '{p_name}')"
                req_name = f'Finding partners4{target_name}'
                t2pG =  t2pG.compose(my_session.process_oql(target_partner_oql,req_name))

      nodes = t2pG._get_nodes()
      partners = set()
      for p in nodes:
          try:
            p_class, p_weight = partner_names[p.name().lower()]
            p['Class'] = [p_class]
            p['target weight'] = [p_weight]
            p['regulator weight'] = [p_weight]
            partners.add(p)
          except KeyError:
              continue

      return partners,t2pG


    def set_partners(self,_4targets:list[PSObject]):
        """
        input:
            self.params['partners'] = {target_name:{partner_class:[partner_names]}}
            Assumes partners are linked to Targets with Effect=positive
        output:
            self.__partners__
            self.partner_class - is used for column names only
            self.find_partners_oql
        """
        try:
            partner_params = dict(self.params['partners'])
            self.__partners__,t2pG = self.params2partners(partner_params,_4targets)
        except KeyError:
            self.partner_class = ResnetGraph.classes(self.__targets__)
            if self.partner_class:
              self.__partners__,t2pG = self.find_partners(_4targets)
        return t2pG
    

    def _partner_class_str(self):
        clases = {ResnetGraph.classes(self.__partners__)}
        return ','.join(clases)
   

    def report_path(self, extension='.xlsx'):
      indics = ','.join(self.params['indication_types'])
      rep_pred = 'suggested ' if self.params['strict_mode'] else 'suggested,predicted '
      if self.__moa() == ANTAGONIST:
        mode = ' inhibition'
      elif self.__moa() == AGONIST:
        mode = ' activation'
      else:
        mode = ''
      
      fname = rep_pred+ indics+'4'+ self.input_target_names_str()+mode+extension
      return os.path.join(self.data_dir, fname)


    def __find_cells_secreting(self,ligands:list[PSObject])->list[PSObject]:
        t_n = ','.join(ResnetGraph.names(ligands))
        oql4ligands = OQL.get_objects(ResnetGraph.dbids(ligands))
        REQUEST_NAME = f'Find cells secreting {t_n}'
        OQLquery = 'SELECT Relation WHERE objectType = (CellExpression,MolTransport) AND NeighborOf ({select_targets}) AND NeighborOf (SELECT Entity WHERE objectType = Cell)'
        secreting_cellsG = self.process_oql(OQLquery.format(select_targets=oql4ligands),REQUEST_NAME)
        if isinstance(secreting_cellsG,ResnetGraph):
            return secreting_cellsG._psobjs_with('CellType','ObjTypeName')
        else:
            print(f'No secreting cells found for {t_n}')
            return []
    

    def find_modulators(self,of_types:list[str],with_effect:str,on_targets:list[PSObject],linked_by:list[str],
                        with_min_refcount=1,drugs_only=False)->list[PSObject]:
        t_n = ','.join(ResnetGraph.names(on_targets))
        f_t = OQL.get_objects(ResnetGraph.dbids(on_targets))
        reltypes = ','.join(linked_by)
        modulator_types =  ','.join(of_types)

        REQUEST_NAME = f'Find substances {with_effect}ly regulating {t_n} by {reltypes}'
        OQLquery = f'SELECT Relation WHERE objectType = ({reltypes}) AND Effect = {with_effect} AND NeighborOf upstream ({f_t})'
        
        if drugs_only:
            nb = f' AND NeighborOf downstream ({OQL.select_drugs()} AND Connectivity > 1)'
        else:
            nb = f' AND NeighborOf downstream (SELECT Entity WHERE objectType = {modulator_types} AND Connectivity > 1)'
        OQLquery += nb

        if with_min_refcount > 1:
             OQLquery += ' AND '+REFCOUNT+' >= '+str(with_min_refcount)

        modulators = list()
        TargetModulatorsNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        if isinstance(TargetModulatorsNetwork,ResnetGraph):
            if self.targets_have_weights:
                TargetModulatorsNetwork.add_weights2neighbors(on_targets,'target weight','regulator weight','<')
            modulators = TargetModulatorsNetwork._psobjs_with('SmallMol','ObjTypeName')
            print(f'Found {len(modulators)} {modulator_types}(s) {with_effect}ly regulating {t_n} by {reltypes} with minimum reference count = {with_min_refcount}')
      
        return set(modulators)
    
    
    def modulators_effects(self):
      '''
      output:
          self.__DirectAgonists
          self.__DirectAntagonists
          self.__IndirectAgonists
          self.__IndirectAntagonists
      updates:
          self.__indications4antagonists__
          self.__indications4agonists__
      '''
      kwargs2findmodulators = {
                      'of_types':['SmallMol'],
                      'with_effect' : 'negative',
                      'on_targets' : self.__targets__
                    }

      #if self.__moa() in [ANTAGONIST,ANY_MOA]:
      kwargs2findmodulators['linked_by'] = ['DirectRegulation']
      self.__DirectAntagonists = self.find_modulators(**kwargs2findmodulators)
      kwargs2findmodulators['linked_by'] = ['Regulation','Expression','MolTransport']
      kwargs2findmodulators['drugs_only'] = True
      self.__IndirectAntagonists = self.find_modulators(**kwargs2findmodulators)
        
      #if self.__moa() in [AGONIST,ANY_MOA]:
      kwargs2findmodulators['with_effect'] = 'positive'
      kwargs2findmodulators['linked_by'] = ['DirectRegulation']
      self.__DirectAgonists = self.find_modulators(**kwargs2findmodulators)
      kwargs2findmodulators['linked_by'] = ['Regulation','Expression','MolTransport']
      kwargs2findmodulators['drugs_only'] = True
      self.__IndirectAgonists = self.find_modulators(**kwargs2findmodulators)

      if not self._is_strict():
        if self.__DirectAntagonists:
          self.__indications4antagonists__.update(self.find_indications4(self.__DirectAntagonists))
          self.__indications4agonists__.update(self.find_toxicities4(self.__DirectAntagonists))
        if self.__DirectAgonists:
          self.__indications4agonists__.update(self.find_indications4(self.__DirectAgonists))
          self.__indications4antagonists__.update(self.find_toxicities4(self.__DirectAgonists))
      return


    def pathway_oql4(self,targets:list[PSObject])->dict[tuple[str,str,int],str]:
        '''
        Return
        ------
        {(target_name,partner_name(),partner_dbid):pathway_oql}
        '''
        #set pathway_name_must_include_target to True if targets have a lot of large curated pathways
        target_oqls = dict()
        pct = '%'
        merged_pathways = 'SELECT Relation WHERE objectType = (DirectRegulation,Binding,ProtModification,PromoterBinding,ChemicalReaction) AND MemberOf ({select_networks})'

        for target in targets:
            target_id = target['Id'][0]
            target_name = target['Name'][0]
            SELECTpathways = 'SELECT Network WHERE ParentOf (SELECT Entity WHERE id = ({i}))'.format(i=str(target_id))
            if self.params['pathway_name_must_include_target']:
                SELECTpathways = SELECTpathways + ' AND Name LIKE (\''+pct+target_name+pct+'\')' #additional refinement for popular targets

            target_oqls[target_name] = SELECTpathways

        oqls = dict() # {tuple[target_name,partner_name,partner_dbid]:str}
        if self.__partners__:
            for partner in self.__partners__:
                for target_name, oql_query in target_oqls.items():
                    find_pathways_query = oql_query + ' AND ParentOf (SELECT Entity WHERE id = ({i}))'.format(i=str(partner.dbid()))
                    oqls[(target_name,partner.name(),partner.dbid())] = merged_pathways.format(select_networks=find_pathways_query)
                    # some ligands have duplicate names - must include dbid into tuple
        else:
            for target_name, oql_query in target_oqls.items():
                oqls[(target_name, '',0)] = merged_pathways.format(select_networks=oql_query)

        return oqls
    

    def load_pathways4(self,targets:list[PSObject]):
        #finding downstream pathway components
        oql_queries = self.pathway_oql4(targets) # separate oql_query for each target
        merged_pathway = ResnetGraph()
        futures = list()
        with ThreadPoolExecutor(max_workers=len(oql_queries), thread_name_prefix='Find curated pathways') as e: 
          for components_tuple, oql_query in oql_queries.items():
            target_name = components_tuple[0]
            partner_name = components_tuple[1]
            request_name = f'Find curated pathways containing {target_name}'
            if partner_name:
              request_name = request_name + ' and ' + partner_name

            futures.append(e.submit(self.process_oql,oql_query,request_name))
          
          for f in as_completed(futures):
              merged_pathway = merged_pathway.compose(f.result())

        if merged_pathway:
            print (f'Found pathway for {self.target_names_str()} with {len(merged_pathway)} components')
            return merged_pathway
        else: 
            print('No curated pathways were found for %s' % self.target_names_str())
            return ResnetGraph()
        

    def  load_pathways_from_folders(self)->ResnetGraph:
        '''
        Requires
        --------
        self.params['pathway_folders'], self.params['pathways']
        '''
        fc = FolderContent(self.APIconfig,what2retrieve=NO_REL_PROPERTIES)
        my_folders = self.params.get('pathway_folders','')
        if not my_folders: 
            raise KeyError

        folder_pathways = list()
        for folder_name in my_folders:
            folder_pathways += fc.folder2pspathways(folder_name,with_layout=False)

        if self.params.get('pathways',''):
            filtered_pathway_list = list()
            for ps_pathway in folder_pathways:
                assert(isinstance(ps_pathway, PSPathway))
                if ps_pathway.name() in self.params['pathways']:
                    filtered_pathway_list.append(ps_pathway)
            folder_pathways = filtered_pathway_list

        print(f'Found {len(folder_pathways)} curated pathways in {my_folders}:')
        [print(p.name()+'\n') for p in folder_pathways]

        # merging pathways 
        merged_pathway = PSPathway()
        [merged_pathway.merge_pathway(p) for p in folder_pathways]
        return merged_pathway.graph


    def get_pathway_components(self,targets:list[PSObject],partners:list[PSObject]):
        try:
            target_pathways = self.load_pathways_from_folders()
        except KeyError:
            target_pathways = self.load_pathways4(targets)
        
        # non-Protein components make link2concept unresponsive:
        target_pathways.remove_nodes_by_prop(['CellProcess','Disease','Treatment'])
        targets_regulome = target_pathways.get_regulome(set(targets))
        self.PathwayComponents = set(targets_regulome._get_nodes())
        self.PathwayComponents = self.PathwayComponents.difference(targets|partners)
        self.Graph.add_psobjs(self.PathwayComponents) # to avoid their retrieval during linking
        print (f'Found regulome for {self.target_names_str()} with {len(self.PathwayComponents)} components')

########################### FIND INDICATIONS ########################## FIND INDICATIONS ###################
    def _indications4(self,targets:list[PSObject],with_effect_on_indications:str)->list[PSObject]:
      '''
      input:
          moa = [AGONIST, ANTAGONIST]
          targets - list[PSObject]
      output:
          indications, indications_strict
      '''
      assert(with_effect_on_indications in ['positive','negative'])
      if not targets: 
        print('No targets are known for the drug')
        return []
      
      indications = set()
      t_n = ResnetGraph.names(targets)
      t_oql = OQL.get_objects(ResnetGraph.dbids(targets))
      i_oql,_ = self.oql4indications()
      my_relprops = self.relProps
      self.add_rel_props(PS_SENTENCE_PROPS)
      REQUEST_NAME = f'Find indications {with_effect_on_indications}ly regulating {t_n}'
      OQLquery = f'SELECT Relation WHERE objectType = Regulation AND Effect = {with_effect_on_indications} AND \
        NeighborOf({t_oql}) AND NeighborOf ({i_oql})' 
      ModulatedByTargetNetwork = self.process_oql(OQLquery,REQUEST_NAME)
      if isinstance(ModulatedByTargetNetwork,ResnetGraph):
        indications = set(ModulatedByTargetNetwork.psobjs_with(only_with_values=self.params['indication_types']))
        print(f'Found {len(indications)} diseases {with_effect_on_indications}ly regulated by {t_n}')

      if not self._is_strict():
        REQUEST_NAME = f'Find indications {with_effect_on_indications}ly modulated by {t_n}'
        OQLquery = f'SELECT Relation WHERE objectType = QuantitativeChange AND Effect = {with_effect_on_indications} AND \
          NeighborOf ({t_oql}) AND NeighborOf ({i_oql})'
        ActivatedInDiseaseNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        if with_effect_on_indications == 'positive':
          REQUEST_NAME = f'Find indications where {t_n} is biomarker'
          OQLquery = f'SELECT Relation WHERE objectType = Biomarker AND NeighborOf ({t_oql}) AND NeighborOf ({i_oql})'
          BiomarkerG = self.process_oql(OQLquery,REQUEST_NAME)
          ActivatedInDiseaseNetwork = ActivatedInDiseaseNetwork.compose(BiomarkerG)

        if isinstance(ActivatedInDiseaseNetwork,ResnetGraph):
          add2indications = ActivatedInDiseaseNetwork.psobjs_with(only_with_values=self.params['indication_types'])
          indications.update(add2indications)
          print(f'Found {len(add2indications)} diseases where {t_n} is {with_effect_on_indications}ly regulated')
      self.relProps = my_relprops
      return indications


    def _indications4targets(self):
      '''
      input:
          if _4targets is empty will use target names from self.__targets__
      output:
          self.__indications4antagonists__
          self.__indications4agonists__
          self.__unknown_effect_indications__
      '''
     # moa = self.params['mode_of_action']
     # if moa in [ANTAGONIST,ANY_MOA]:
      #  effect_on_indication = 'positive'
      self.__indications4antagonists__ = self._indications4(self.__targets__,'positive')
     # if moa in [AGONIST,ANY_MOA]:
      #  effect_on_indication = 'negative'
      self.__indications4agonists__ = self._indications4(self.__targets__,'negative')
      self.___unknown_effect_indications__(self.__known_effect_indications())
      return


    def _indications4partners(self,targets:list[PSObject],partners:list[PSObject],effect_on_indications:str)->list[PSObject]:
        OQLtemplate = 'SELECT Relation WHERE objectType = (Regulation,QuantitativeChange) AND \
Effect = {eff} AND NeighborOf ({partners}) AND NeighborOf ({indications})'
        assert(effect_on_indications in ['positive','negative'])

        #effect = 'positive' if moa == ANTAGONIST else 'negative'
        if not targets or not partners: return []
        t_n = ','.join(ResnetGraph.names(targets))
        p_c = ','.join(ResnetGraph.classes(partners))
        oql4indications,_ = self.oql4indications()
        oql4partners = OQL.get_objects(ResnetGraph.dbids(partners))

        OQLquery = OQLtemplate.format(eff=effect_on_indications,partners=oql4partners,indications=oql4indications)
        REQUEST_NAME = f'Find indications {effect_on_indications}ly regulated by {p_c}s of {t_n}'
        my_relprops = self.relProps
        self.add_rel_props(PS_SENTENCE_PROPS)
        PartnerIndicationNetwork4anatagonists = self.process_oql(OQLquery,request_name=REQUEST_NAME)
        indications = PartnerIndicationNetwork4anatagonists.psobjs_with(only_with_values=self.params['indication_types'])
        print(f'Found {len(indications)} indications for {len(partners)} {t_n} {p_c}')
        self.relProps = my_relprops
        return indications


    def __indications4partners(self):
        """
        input:
            self.__partners__
        Assumes partners are linked to Target with Effect=positive
        """
        t_n = self.target_names_str()
        exist_indication_count = len(self._my_indications())
    #    moa = self.params['mode_of_action']
   #     if moa in [ANTAGONIST,ANY_MOA]:
   #         effect_on_indications = 'positive'
        self.__indications4agonists__.update(self._indications4partners(self.__targets__,self.__partners__,'positive'))
   #     if moa in [AGONIST,ANY_MOA]:
   #         effect_on_indications = 'negative'
        self.__indications4agonists__.update(self._indications4partners(self.__targets__,self.__partners__,'negative'))
  
        new_indication_count = len(self._my_indications()) - exist_indication_count
        print('%d indications for %d %s partnerss were not found by previous searches' %  
                (new_indication_count, len(self.__partners__),t_n))
        return
    

    def find_indications4(self,modulators:list[PSObject])->set[PSObject]:
        if modulators:
          assert(len(modulators) < 1000)
          indications = set()
          get_modulators = OQL.get_objects(ResnetGraph.dbids(modulators))
          indications_type=','.join(self.params['indication_types'])
          REQUEST_NAME = f'Find indications for {len(modulators)} modulators'
          OQLquery = f'SELECT Relation WHERE objectType = Regulation AND Effect = negative AND \
          NeighborOf (SELECT Entity WHERE objectType = ({indications_type})) AND NeighborOf ({get_modulators})'
          InhibitorsIndicationNetwork = self.process_oql(OQLquery,REQUEST_NAME)
          if isinstance(InhibitorsIndicationNetwork,ResnetGraph):
              indications = set(InhibitorsIndicationNetwork.psobjs_with(only_with_values=self.params['indication_types']))
              print(f'Found {len(indications)} indications for {len(modulators)} substances')

          REQUEST_NAME = f'Find clinical trials for {len(modulators)} modulators'
          OQLquery = f'SELECT Relation WHERE objectType = ClinicalTrial AND \
          NeighborOf (SELECT Entity WHERE objectType = ({indications_type})) AND NeighborOf ({get_modulators})'
          ct_g = self.process_oql(OQLquery,REQUEST_NAME)
          ct_indications = ct_g.psobjs_with(only_with_values=self.params['indication_types'])
          if isinstance(ct_g,ResnetGraph):
              indications.update(ct_indications)

          print(f'Found {len(ct_indications)} indications on clinical trials with {len(modulators)} substances')
          return indications
        else: 
            return []


    def find_toxicities4(self,modulators:list[PSObject])->set[PSObject]:
        if modulators:
          assert(len(modulators) < 1000)
          toxicities = set()
          get_modulators = OQL.get_objects(ResnetGraph.dbids(modulators))
          indications_type=','.join(self.params['indication_types'])
          REQUEST_NAME = f'Find toxicities for {len(modulators)} modulators'
          OQLquery = f'SELECT Relation WHERE objectType = Regulation AND Effect = positive AND \
          NeighborOf (SELECT Entity WHERE objectType = ({indications_type})) AND NeighborOf ({get_modulators})'
          InhibitorsIndicationNetwork = self.process_oql(OQLquery,REQUEST_NAME)
          if isinstance(InhibitorsIndicationNetwork,ResnetGraph):
              toxicities = set(InhibitorsIndicationNetwork.psobjs_with(only_with_values=self.params['indication_types']))
              print(f'Found {len(toxicities)} indications for {len(modulators)} substances')

          print(f'Found {len(toxicities)} toxicities on clinical trials with {len(modulators)} substances')
          return toxicities
        else:
            return []


    def __indications4cells(self,secreting:list[PSObject],targets:list[PSObject],with_effect_on_indication:str)->list[PSObject]:
        '''
           effect_on_indication in ['positive','negative']
        '''
        assert(with_effect_on_indication in ['positive','negative'])
        t_n = ResnetGraph.names(targets)
        def best_cell_indications(max_indication_count:int, cell2disease:ResnetGraph):
            '''
            Return
            ------
            Indications that belong to cells having most # of indications (>max_indication_count)\n
            These are the most popular indications for cells producing input targets above max_indication_count threashold
            '''
            disease_uid2indigree = {uid:cell2disease.in_degree(uid) for uid,o in cell2disease.nodes(data=True) if o['ObjTypeName'][0] in self.params['indication_types']}
            disease_uid2indigree = sorted(disease_uid2indigree.items(), key=lambda x: x[1],reverse=True)
            min_indication_indegree = disease_uid2indigree[max_indication_count][1]
            best_indications_uids = [uid for uid,v in disease_uid2indigree if v >= min_indication_indegree]
            return cell2disease._get_nodes(best_indications_uids)

        #effect = 'positive' if moa == ANTAGONIST else 'negative'
        REQUEST_NAME = f'Find indications {with_effect_on_indication}ly linked to cell secreting {t_n}'
        OQLtemplate = 'SELECT Relation WHERE Effect = {effect} AND NeighborOf (SELECT Entity WHERE \
            objectType = ({indication_type})) AND NeighborOf (SELECT Entity WHERE id = ({cell_ids}))'
        
        secreting_cells_dbids = ResnetGraph.dbids(secreting)
        OQLquery = OQLtemplate.format(effect=with_effect_on_indication,cell_ids=OQL.id2str(secreting_cells_dbids),indication_type=','.join(self.params['indication_types']))
        CellDiseaseNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        if isinstance(CellDiseaseNetwork,ResnetGraph):
            if self.targets_have_weights:
                CellDiseaseNetwork.add_weights2neighbors(targets,'target weight','regulator weight')

            indications = CellDiseaseNetwork.psobjs_with(only_with_values=self.params['indication_types'])
            if len(indications) > self.max_cell_indications:
                indications = best_cell_indications(self.max_cell_indications,CellDiseaseNetwork)
        
        return indications


    def _indications4cells(self,secreting:set[PSObject],targets:list[PSObject],effect_on_indications:str)->list[PSObject]:
        '''
        input:
          "secreting" - list of cells secreting "targets"
        ouput:
          indications,
          if "secreting" is empty adds cells to "secreting"
        '''
        if targets:
          assert(effect_on_indications in ['positive','negative'])
          if self.params.get('include_indications4cells_expressing_targets',False):
            my_targets = targets
          else:
            my_targets = [o for o in targets if o.get_prop('Class') == 'Ligand']
          
          if my_targets:
            t_n = ResnetGraph.names(my_targets)
            if not secreting:
              secreting.update(self.__find_cells_secreting(my_targets))

            if secreting:
              my_relProps = self.relProps
              self.add_rel_props(PS_SENTENCE_PROPS)
              indications = self.__indications4cells(secreting,my_targets,effect_on_indications)
              print(f'Found {len(indications)} indications for {len(secreting)} cells producing {t_n}')
              self.relProps = my_relProps
              return indications
            else:
              print(f'No cells producing {t_n} were found')
              return []
          else:
              print('Target list contains no ligands. No producing cell can be found')
              return []
        else:
            print(f'Target list for {effect_on_indications} effect on indications is empty')
            return []
        

    def indications4secreting_cells(self):
        if self._is_strict(): 
            print('Script runs in strict mode.  No indications for cells secreting targets will be used')
            self.__targets__secretingCells = self.__find_cells_secreting(self.__targets__)
            t_n = ','.join(ResnetGraph.names(self.__targets__))
            print('Found %d cell types producing %s' % (len(self.__targets__secretingCells),t_n))
            # loaded here self.__targets__secretingCells is used for target ranking
            return
        else:
          exist_indication_count = len(self._my_indications())
          if self.__moa() in [ANTAGONIST,ANY_MOA]:
              indications,secreting_cells = self._indications4cells(self.__targets__secretingCells,self.__targets__,ANTAGONIST)
              self.__targets__secretingCells.update(secreting_cells)
              self.__indications4antagonists__.update(indications)

          if self.__moa() in [AGONIST,ANY_MOA]:
              indications, secreting_cells = self._indications4cells(self.__targets__secretingCells,self.__targets__,AGONIST)
              self.__targets__secretingCells.update(secreting_cells)
              self.__indications4agonists__.update(indications)
              
          new_indication_count = len(self._my_indications()) - exist_indication_count
          print('%d indications for %d cells producing %s were not found by previous searches' %  
                  (new_indication_count, len(self.__targets__secretingCells),self.target_names()))
          return


    def _resolve_conflict_indications(self):
      def __resolve(conflict:PSObject, all_rels:list, using_rel_type:str):
        only_rels_with_type = [r for r in all_rels if r.objtype() == using_rel_type]
        if only_rels_with_type:
          only_rels_with_type.sort(key=lambda x: x.count_refs(), reverse=True)
          best_effect = only_rels_with_type[0]['Effect'][0]
          if best_effect == 'positive':
            self.__indications4agonists__.discard(conflict)
            return True
          elif best_effect == 'negative':
            self.__indications4antagonists__.discard(conflict)
            return True
          else:
            return False
        return False

      conflict_indications = self.__indications4antagonists__.intersection(self.__indications4agonists__)
      unique_indications = len(self.__indications4antagonists__)+len(self.__indications4agonists__)-len(conflict_indications)
      if conflict_indications:
        print(f'Resolving {len(conflict_indications)} conflicting out of total {unique_indications} indications')
        unresolved_uids = set()
        target_uids = ResnetGraph.uids(self.__targets__)
        partner_uids = ResnetGraph.uids(self.__partners__)
        for conflict in conflict_indications:
            target_rels = list(self.Graph.find_relations(target_uids,[conflict.uid()]))
            if not __resolve(conflict,target_rels,'Regulation'):
                if not __resolve(conflict,target_rels,'QuantitativeChange'):
                    partner_rels = list(self.Graph.find_relations(partner_uids,[conflict.uid()]))
                    if not __resolve(conflict,partner_rels,'Regulation'):
                        if not __resolve(conflict,partner_rels,'QuantitativeChange'):
                            unresolved_uids.add(conflict.uid())
        if unresolved_uids:
          print(f'{len(unresolved_uids)} indications cannot be resolved.\n They will appear in both worksheets for agonists and antagonist:')
          for uid in unresolved_uids:
            try:
              print(self.Graph._psobj(uid).name())
            except KeyError:
              continue
            
      if 'include_indications4ontology_groups' in self.params:
        if not hasattr(self,'mustbe_indications'):
          ontology_group = self.params['include_indications4ontology_groups']
          print(f'Loading indications from {ontology_group} ontology group')
          oql = OQL.get_childs(ontology_group,['Name'])
          mustbe_indications = self.process_oql(oql)._get_nodes()
          childs,self.mustbe_indications = self.load_children4(mustbe_indications,
                                           max_childs=self.max_ontology_parent)
          
          #terminal_childs = [o for o in childs if not o.childs()]
          #self.mustbe_indications = set(self.mustbe_indications+terminal_childs)
        if self.__moa() == ANTAGONIST:
          if self.params.get('add_abstcooccur',False) or self.__coocG__:
            self.__indications4antagonists__.update(self.mustbe_indications)
          len_before = len(self.__indications4agonists__)
          self.__indications4agonists__ = self.__indications4agonists__.difference(self.mustbe_indications)
          print(f'Deleted {len_before-len(self.__indications4agonists__)} toxicities belonging to {self.params['include_indications4ontology_groups']}')
        elif self.__moa() == AGONIST:
          if self.params.get('add_abstcooccur',False) or self.__coocG__:
            self.__indications4agonists__.update(self.mustbe_indications)
          len_before = len(self.__indications4antagonists__)
          self.__indications4antagonists__ = self.__indications4antagonists__.difference(self.mustbe_indications)
          print(f'Deleted {len_before-len(self.__indications4antagonists__)} toxicities belonging to {self.params['include_indications4ontology_groups']}')

      if 'include_toxicity4ontology_groups' in self.params:
        if not hasattr(self,'mustbe_toxicities'):
          oql = OQL.get_childs(self.params['include_toxicity4ontology_groups'],['Name'])
          mustbe_toxicities = self.process_oql(oql)._get_nodes()
          childs,self.mustbe_toxicities = self.load_children4(mustbe_toxicities,
                                           max_childs=self.max_ontology_parent)
          
        if self.__moa() == ANTAGONIST:
          len_before = len(self.__indications4antagonists__)
          self.__indications4antagonists__ = self.__indications4antagonists__.difference(self.mustbe_toxicities)
          print(f'Deleted {len_before-len(self.__indications4antagonists__)} indications belonging to {self.params['include_toxicity4ontology_groups']}')
        elif self.__moa() == AGONIST:
          len_before = len(self.__indications4agonists__)
          self.__indications4agonists__ = self.__indications4agonists__.difference(self.mustbe_toxicities)
          print(f'Deleted {len_before-len(self.__indications4agonists__)} indications belonging to {self.params['include_toxicity4ontology_groups']}')
       
      if 'exclude_ontology_groups' in self.params:
        if not hasattr(self,'mustnotbe'):
          oql = OQL.get_childs(self.params['exclude_ontology_groups'],['Name'],include_parents=True)
          self.mustnotbe = self.process_oql(oql)._get_nodes()
        len_before = len(self.__indications4antagonists__)
        self.__indications4antagonists__ = self.__indications4antagonists__.difference(self.mustnotbe)
        print(f'Deleted {len_before-len(self.__indications4antagonists__)} indications belonging to {self.params['exclude_ontology_groups']}')
        len_before = len(self.__indications4agonists__)
        self.__indications4agonists__ = self.__indications4agonists__.difference(self.mustnotbe)
        print(f'Deleted {len_before-len(self.__indications4agonists__)} indications belonging to {self.params['exclude_ontology_groups']}')
        

    def indications4targets(self):
        '''
        Input
        -----
        target names must be in self.params['target_names']
        moa in [ANTAGONIST, AGONIST, ANY_MOA]
        '''
        #assert(moa in [ANTAGONIST, AGONIST, ANY_MOA])
        start_time = time.time()
        self._set_targets()
        self.set_partners(self.__targets__)
        self._indications4targets()
        self.__indications4partners()

        self._resolve_conflict_indications()
        print("%d indications for %s were found in %s" % 
        (len(self._my_indications()), self.target_names_str(), execution_time(start_time)))

        return self._my_indications()

############### UNKNOWN EFFECT INDICATIONS ######################### UNKNOWN EFFECT INDICATIONS ####################
    def _GVindicationsG(self,_4targets:list[PSObject])->ResnetGraph:
        '''
        output:
            GVsInDiseaseNetwork where GVs are annoated with target weights if self.targets_have_weights
        '''
        t_n = ','.join([x.name() for x in _4targets])
        oql4targets = OQL.get_objects(ResnetGraph.dbids(_4targets))
        #selectGVs = f'SELECT Entity WHERE objectType = GeneticVariant AND Connected by \
         #   (SELECT Relation WHERE objectType = GeneticChange) to ({oql4targets})'
        
        findGVs_oql = f'SELECT Relation WHERE objectType = GeneticChange AND NeighborOf upstream ({oql4targets}) \
AND NeighborOf downstream (SELECT Entity WHERE objectType = GeneticVariant)'
        
        REQUEST_NAME = f'Find genetic variants linked to {t_n}'
        GV2targetsG = self.process_oql(findGVs_oql,REQUEST_NAME)
        allGVs = GV2targetsG._psobjs_with('GeneticVariant',OBJECT_TYPE)
        if allGVs:
          selectGVs = OQL.get_objects(ResnetGraph.dbids(allGVs))
          REQUEST_NAME = f'Find indications linked to {t_n} genetic variants'
          indication_type=','.join(self.params['indication_types'])
          OQLquery = f'SELECT Relation WHERE NeighborOf({selectGVs}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_type}))'
          GVsInDiseaseNetwork = self.process_oql(OQLquery,REQUEST_NAME)
          if isinstance(GVsInDiseaseNetwork,ResnetGraph):
              if self.targets_have_weights:
                  # need to do it in both directions because GV-Disease link is non-directional
                  GVsInDiseaseNetwork = GVsInDiseaseNetwork.compose(GV2targetsG)
                  GVsInDiseaseNetwork.add_weights2neighbors(_4targets,'target weight','regulator weight','<')
                  #GVsInDiseaseNetwork.add_weights2neighbors(_4targets,'regulator weight','target weight','>')
              return GVsInDiseaseNetwork
          else:
             return ResnetGraph()
        else:
            return ResnetGraph()
        

    def GVindications(self)->list[PSObject]:
        self.GVs2DiseaseGraph = self._GVindicationsG(self.__targets__)
        if isinstance(self.GVs2DiseaseGraph,ResnetGraph):
            indications = self.GVs2DiseaseGraph.psobjs_with(only_with_values=self.params['indication_types'])
            self.__GVs__ = set(self.GVs2DiseaseGraph._psobjs_with('GeneticVariant',OBJECT_TYPE))

            t_n = self.target_names_str()
            print('Found %d indications genetically linked to %d Genetic Variants in %s' % 
                (len(indications), len(self.__GVs__), t_n))
            return indications
        else:
            return list()


    def _biomarker_indicationsG(self,_4targets:list[PSObject]):
        t_n = self.target_names_str()
        REQUEST_NAME = f'Find indications where {t_n} is biomarker'
        OQLquery = 'SELECT Relation WHERE objectType = Biomarker AND NeighborOf({select_target}) AND NeighborOf ({indications})'

        f_t = OQL.get_objects(ResnetGraph.dbids(_4targets))
        oql4indications,_ = self.oql4indications()
        OQLquery = OQLquery.format(select_target=f_t,indications=oql4indications)
        BiomarkerInDiseaseNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        if isinstance(BiomarkerInDiseaseNetwork,ResnetGraph):
            return BiomarkerInDiseaseNetwork
        else:
            return ResnetGraph()
         

    def ___unknown_effect_indications__(self,known_indications:list[PSObject]=[]):
        t_n = self.target_names_str()
        gv_indications = set(self.GVindications())
        print(f'{len(gv_indications)} indications linked to {t_n} genetic variations')

        BiomarkerInDiseaseNetwork = self._biomarker_indicationsG(self.__targets__)
        biomarker_indications = set(BiomarkerInDiseaseNetwork.psobjs_with(only_with_values=self.params['indication_types']))
        print('Found %d indications where target %s is claimed as biomarker' %  (len(biomarker_indications),t_n))
        
        indications = gv_indications|biomarker_indications
        indications = indications.difference(known_indications)
        print(f'Found {len(indications)} indications linked with unknown effect to targets')
        self.__unknown_effect_indications__ = indications
        return indications

        
##################  SCORING SCORE SCORING ######################### SCORING SCORE SCORING ####################
    def score_GVs(self, df2score:df):
      if not self.__GVs__: 
          return
      t_n = self.target_names_str()
      concept_name = ' genetic variants'
      self.__colname4GV__ = t_n+concept_name
      refcount_column = self._refcount_colname(concept_name)
      weighted_refcount_column = self._weighted_refcount_colname(concept_name) 
      linked_count_column = self._linkedconcepts_colname(concept_name)
      concept_size_column = self._concept_size_colname(concept_name)

      gvlinkcounter = 0
      if hasattr(self,'GVs2DiseaseGraph'):
          if self.targets_have_weights:
            regurn2weight = {gv.urn():gv.get_prop('regulator weight') for gv in self.__GVs__}
          
          totalGVcount = len(self.__GVs__)
          df2score.insert(len(df2score.columns),weighted_refcount_column,[float(0)]*len(df2score))
          df2score.insert(len(df2score.columns),refcount_column,[0]*len(df2score))
          df2score.insert(len(df2score.columns),linked_count_column,[0]*len(df2score))
          df2score.insert(len(df2score.columns),concept_size_column,[totalGVcount]*len(df2score))
          
          for i in df2score.index:
            row_indication_dbids = set(df2score.at[i,self.__temp_id_col__])
            row_indications = set(self.Graph.psobj_with_dbids(row_indication_dbids))
            GVscore = 0.0
            rowGVs_count = 0
            refcount = 0
            if isinstance(self.GVs2DiseaseGraph,ResnetGraph):
              rowGVs = list(self.GVs2DiseaseGraph.get_neighbors(row_indications,self.__GVs__))
              rowGVs_count = len(rowGVs)
              if len(rowGVs) > 0:
                GVnames = [g.name() for g in rowGVs]
                df2score.at[i,self.__colname4GV__] = ';'.join(GVnames)
                gvlinkcounter += 1
                gv_disease_subgraph = self.GVs2DiseaseGraph.get_subgraph(list(row_indications),rowGVs)
                if self.targets_have_weights:
                  # GV-disease associations are non-directional, therfore GV may be a target or regulator in gv_disease_subgraph
                  gv_disease_subgraph.add_node_weight2ref(regurn2weight,regurn2weight)
                  subgraph_refs = gv_disease_subgraph.load_references()
                  ref_weights = [r.get_weight('nodeweight') for r in subgraph_refs]
                  GVscore = float(sum(ref_weights))
                  refcount =  len(subgraph_refs)
                else:
                  subgraph_refs = gv_disease_subgraph.load_references()
                  GVscore = len(subgraph_refs)
                  refcount =  len(subgraph_refs)
        
                df2score.at[i,weighted_refcount_column] = GVscore
                df2score.at[i,refcount_column] = refcount
                df2score.at[i,linked_count_column] = rowGVs_count
          
          self._set_rank(df2score,concept_name)
          print('Found %d indications linked to %d GVs' % (gvlinkcounter, len(self.__GVs__)))


    def _drug_connect_params(self,direct_modulators:bool,score_antagonists:bool):
        '''
        Return
        ------
        how2connect parameters to score substances having desired drug affect on the targets 
        '''
        if score_antagonists:
        # most common case when targets must be inhibited
            if direct_modulators:
                effect_on_indication = 'negative'
                drug_class = 'direct antagonists'         
                concepts = self.__DirectAntagonists
            else:
                effect_on_indication = 'negative'
                drug_class = 'indirect antagonists'             
                concepts = self.__IndirectAntagonists
        else:
        # case if drug are agonists
            if direct_modulators:
                effect_on_indication = 'positive'
                drug_class = 'direct agonists'
                concepts = self.__DirectAgonists
            else:
                effect_on_indication = 'positive'
                drug_class = 'indirect agonists'
                concepts = self.__IndirectAgonists

        return effect_on_indication, drug_class, concepts


    def _drug_tox_params(self,direct_modulators:bool,score_antagonists:bool):
        '''
        Return
        ------
        how2connect parameters to score substances synergizing with target action on indication\n
        i.e. molecules having effect opposite to desired drug affect on targets
        '''
        if score_antagonists:
        # most common case when targets must be inhibited
            if direct_modulators:
                effect_on_indication = 'positive'
                drug_class = 'direct agonists'
                concepts = self.__DirectAgonists
            else:
                effect_on_indication = 'positive'
                drug_class = 'indirect agonists'
                concepts = self.__IndirectAgonists
        else:
        # case if drug are agonists
            if direct_modulators:
                effect_on_indication = 'negative'
                drug_class = 'direct antagonists'             
                concepts = self.__DirectAntagonists
            else:
                effect_on_indication = 'negative'
                drug_class = 'indirect antagonists'             
                concepts = self.__IndirectAntagonists

        return effect_on_indication, drug_class, concepts


    def semscore4(self,targets:list[PSObject],with_effect_on_indication:str,
                         with_partners:list[PSObject],in_df:df):
        """
        Input
        -----
        target_effect_on_indication: required Effect sign between target and Indication 
        target_effect_on_indication = 'positive' to score antagonists
        target_effect_on_indication = 'negative' to score agonists

        Output
        ------
        adds Refcount score columns "in_worksheet" from self.raw_data
        """
        my_df = in_df.dfcopy()
        booster_reltypes = ['Regulation','Biomarker','GeneticChange','QuantitativeChange','StateChange','FunctionalAssociation']
        
        t_n = self.target_names_str()
        target_in_header = t_n if len(t_n) < 45 else 'targets'
        if target_in_header == 'targets':
             concept = 'Regulated by ' + target_in_header
        elif with_effect_on_indication == 'positive':
            concept = 'Activated by '+target_in_header
        else:
            concept = 'Inhibited by '+target_in_header
        
        kwargs = {'connect_by_rels' : ['Regulation'],
                  'with_effects' : [with_effect_on_indication],
                  'boost_by_reltypes' : booster_reltypes
                  }
        if self.targets_have_weights:
            kwargs.update({'nodeweight_prop': 'regulator weight'})
        how2connect = self.set_how2connect(**kwargs)
        linked_row_count,_,my_df = self.link2concept(concept,targets,my_df,how2connect)
        self._set_rank(my_df,concept)
        print('%d indications are %sly regulated by %s' % 
            (linked_row_count,with_effect_on_indication,t_n))
        self.nodeweight_prop = ''

        self.score_GVs(my_df)

        score4antagonists = True if with_effect_on_indication == 'positive' else False
        link_effect, drug_class, drug_connect_concepts = self._drug_connect_params(direct_modulators=True,score_antagonists=score4antagonists)
        if drug_connect_concepts:
            # references suggesting that known drugs for the target as treatments for indication
            concept = target_in_header+' '+drug_class+' clin. trials'
            concept = target_in_header+' '+drug_class
            kwargs = {'connect_by_rels':['ClinicalTrial']}
            if self.targets_have_weights:
                kwargs.update({'nodeweight_prop': 'regulator weight'})
            how2connect = self.set_how2connect(**kwargs)
            linked_row_count,_,my_df = self.link2concept(concept,drug_connect_concepts,my_df,how2connect)
            self._set_rank(my_df,concept)
            print(f'Linked {linked_row_count} clinical trial indictions for {t_n} {drug_class}')

            concept = target_in_header+' '+drug_class
            kwargs = {'connect_by_rels':['Regulation'],
                  'with_effects' : [link_effect],
                  'boost_by_reltypes' :['Regulation','FunctionalAssociation']
                  }
            how2connect = self.set_how2connect(**kwargs)
            linked_row_count,_,my_df = self.link2concept(concept,drug_connect_concepts,my_df,how2connect)
            self._set_rank(my_df,concept)
            print(f'Linked {linked_row_count} indications for {t_n} {drug_class}')

        #references reporting target agonists exacerbating indication or causing indication as adverse events
        link_effect, drug_class, concepts = self._drug_tox_params(direct_modulators=True,score_antagonists=score4antagonists)
        if concepts:
            concept = target_in_header+' '+drug_class
            kwargs = {'connect_by_rels':['Regulation'],
                  'with_effects' : [link_effect],
                  'boost_by_reltypes' :['Regulation','FunctionalAssociation']
                  }
            how2connect = self.set_how2connect(**kwargs)
            linked_row_count,_,my_df = self.link2concept(concept,concepts,my_df,how2connect)
            max_rank = my_df.max_colrank()
            if not drug_connect_concepts:
              max_rank += 1 # advancing max_rank only if drug_connect_concepts were empty, otherwise they have the same max_rank
            self._set_rank(my_df,concept,max_rank)

            print(f'Linked {linked_row_count} indications as toxicities for {t_n} {drug_class}')

        #references where target expression or activity changes in the indication
        if target_in_header == 'targets':
             indication_types = ','.join(self.params['indication_types'])
             concept = target_in_header + f' changes in {indication_types}'
        elif with_effect_on_indication == 'positive':
            concept = target_in_header+' is upregulated'
        else:
            concept = target_in_header+' is downregulated'

        # indications that could not be connected in "Regulated by" are considered here
        kwargs = {'connect_by_rels':['QuantitativeChange'],
                  'with_effects' : [with_effect_on_indication],
                  'boost_by_reltypes' : ['Biomarker','StateChange','FunctionalAssociation'],
                  }
        if self.targets_have_weights:
            kwargs.update({'nodeweight_prop': 'target weight'})
        how2connect = self.set_how2connect(**kwargs)
        linked_row_count,_,my_df= self.link2concept(concept,targets,my_df,how2connect)
        self._set_rank(my_df,concept)
        print(f'{linked_row_count} indications {with_effect_on_indication}ly regulate {t_n}')
        self.nodeweight_prop = ''

        #references suggesting target partners as targets for indication
        if with_partners:
          concept = f'{target_in_header} partners'
          kwargs = {'connect_by_rels':['Regulation'],
                'with_effects' : [with_effect_on_indication],
                'boost_by_reltypes' : ['Regulation','FunctionalAssociation']
                }
          if self.targets_have_weights:
            kwargs.update({'nodeweight_prop': 'regulator weight'})
          how2connect = self.set_how2connect(**kwargs)
          linked_row_count,_,my_df = self.link2concept(concept,with_partners,my_df,how2connect)           
          self._set_rank(my_df,concept)
          print('Linked %d indications for %d %s partners' % 
              (linked_row_count,len(with_partners),t_n))
          self.nodeweight_prop = ''

        # references reporting that cells producing the target linked to indication  
        # only used if taregts are secretred ligands
        if hasattr(self, '__targets__secretingCells'):
            concept = f'{target_in_header} secreting cells'
            kwargs = {'connect_by_rels':['Regulation'],
                  'with_effects' : [with_effect_on_indication],
                  'boost_by_reltypes' : ['Regulation','FunctionalAssociation']
                  }
            if self.targets_have_weights:
                kwargs.update({'nodeweight_prop': 'regulator weight'})
            how2connect = self.set_how2connect(**kwargs)
            linked_row_count,_,my_df = self.link2concept(concept,self.__targets__secretingCells,my_df,how2connect)           
            self._set_rank(my_df,concept)
            print(f'Liked {linked_row_count} indications linked {len(self.__targets__secretingCells)} cells producing {t_n}')


        # references suggesting that known drugs for the target as treatments for indication
        link_effect, drug_class, drug_connect_concepts = self._drug_connect_params(direct_modulators=False,
                                                                      score_antagonists=score4antagonists)
        if drug_connect_concepts:
            concept = target_in_header+' '+drug_class+' clin. trials'
            kwargs = {'connect_by_rels':['ClinicalTrial']}
            if self.targets_have_weights:
                kwargs.update({'nodeweight_prop': 'regulator weight'})
            # cloning session to avoid adding relations to self.Graph
            new_session = self._clone(to_retrieve=REFERENCE_IDENTIFIERS)
            how2connect = new_session.set_how2connect(**kwargs)
            linked_row_count,_,my_df = new_session.link2concept(concept,drug_connect_concepts,my_df,how2connect)           
            self._set_rank(my_df,concept)
            print(f'Linked {linked_row_count} clinical trial indications for {t_n} {drug_class}')
            new_session.close_connection()

            concept = target_in_header+' '+drug_class
            kwargs = {'connect_by_rels':['Regulation'],
                  'with_effects' : [link_effect],
                  'boost_by_reltypes' : ['Regulation'] # boosting with unknown effect Regulation
                  }
            if self.targets_have_weights:
                kwargs.update({'nodeweight_prop': 'regulator weight'})
            # cloning session to avoid adding relations to self.Graph
            new_session = self._clone(to_retrieve=REFERENCE_IDENTIFIERS)
            how2connect = new_session.set_how2connect(**kwargs)
            linked_row_count,_,my_df = new_session.link2concept(concept,drug_connect_concepts,my_df,how2connect)          
            self._set_rank(my_df,concept)
            print('Linked %d indications for %s %s' % (linked_row_count,t_n,drug_class))
            new_session.close_connection()

        #references reporting target agonists exacerbating indication or causing indication as adverse events
        link_effect, drug_class, concepts = self._drug_tox_params(direct_modulators=False,
                                                                   score_antagonists=score4antagonists)
        if concepts:
            concept = target_in_header+' '+drug_class
            kwargs = {'connect_by_rels':['Regulation'],
                  'with_effects' : [link_effect],
                  'boost_by_reltypes' : ['Regulation'] # boosting with unknown effect Regulation
                  }
            if self.targets_have_weights:
                kwargs.update({'nodeweight_prop': 'regulator weight'})
            # cloning session to avoid adding relations to self.Graph
            new_session = self._clone(to_retrieve=REFERENCE_IDENTIFIERS)
            how2connect = new_session.set_how2connect(**kwargs)
            linked_row_count,_,my_df = new_session.link2concept(concept,concepts,my_df,how2connect)
            max_colrank = my_df.max_colrank()
            if not drug_connect_concepts:
              max_colrank += 1 # advancing max_rank only if drug_connect_concepts were empty, otherwise they have the same max_rank
            self._set_rank(my_df,concept,max_colrank)
            print('Linked %d indications as toxicities for %s %s' % (linked_row_count,t_n,drug_class))
            new_session.close_connection()
        
        #references linking target pathways to indication
        if hasattr(self, 'PathwayComponents'):
            concept = target_in_header + ' pathway components'
            kwargs = {'connect_by_rels':['Regulation'],
                  'with_effects' : [with_effect_on_indication],
                  'boost_by_reltypes' : ['Regulation'],  # boosting with unknown effect Regulation
                  'step' : 50
                  }
            # cloning session to avoid adding relations to self.Graph:
            new_session = self._clone(to_retrieve=REFERENCE_IDENTIFIERS)
            how2connect = new_session.set_how2connect(**kwargs)
            linked_row_count,_,my_df = new_session.link2concept(concept,list(self.PathwayComponents),my_df,how2connect)
            
            self._set_rank(my_df,concept)
            print('Linked %d indications to %s pathway components' % (linked_row_count,t_n))
            new_session.close_connection()
 
        return my_df


################ REPORT MAKING REPORT ######################### REPORT MAKING REPORT ##################
    def init_semantic_search(self):
        '''            
        Loads:
            RAWDF4ANTAGONISTS and/or RAWDF4AGONISTS worksheets to raw_data
        '''
        print('\n\nInitializing semantic search')
        t_n = self.target_names_str()

        if self.__indications4antagonists__:
            antagonist_indication_df = self.load_df(list(self.__indications4antagonists__),
                                         max_childs=self.max_ontology_parent,
                                         max_threads=self.max_threads4ontology)
            antagonist_indication_df._name_ = RAWDF4ANTAGONISTS
            self.add2raw(antagonist_indication_df)
            print(f'Will score {len(antagonist_indication_df)} indications for {t_n} antagonists')

        if self.__indications4agonists__:
            agonist_indication_df = self.load_df(list(self.__indications4agonists__),
                                         max_childs=self.max_ontology_parent,
                                         max_threads=self.max_threads4ontology)
            agonist_indication_df._name_ = RAWDF4AGONISTS
            self.add2raw(agonist_indication_df)
            print(f'Will score {len(agonist_indication_df)} indications for {t_n} agonists')
        
        if self.__indications4antagonists__ or self.__indications4agonists__:
            return True
        else:
            print (f'No indications found for {t_n}')
            if self._is_strict(): print('Try setting strict_mode to False')
            return False


    def perform_semantic_search(self):
       self._perform_semantic_search()
       self.add_ontologies()


    def raw2report(self):
      '''
      output:
        normalized copies of worksheets from raw_data are moved to report_pandas
        raw_data[ws4ind],raw_data[ws4tox], report_pandas[ws4ind],report_pandas[ws4tox]
      '''
      empty_cols2drop = list(self.params.get('drop_empty_columns_from_report',[]))
      # at this point raw_data has RAWDF4ANTAGONISTS and RAWDF4AGONISTS worksheets
      # we need to know which worksheet is for indications and which for toxicities based on self.params['mode_of_action']
      ws4ind, ws4tox = self.ws_names()
      if self.__moa() == ANTAGONIST:
        df4ind = self.raw_data.get(RAWDF4ANTAGONISTS,df())
        df4tox = self.raw_data.get(RAWDF4AGONISTS,df())
      else:
        df4ind = self.raw_data.get(RAWDF4AGONISTS,df())
        df4tox = self.raw_data.get(RAWDF4ANTAGONISTS,df())
         
      if not df4ind.empty:
        count_df = self.make_count_df(df4ind,ws4ind)
        self.add2raw(count_df)
        self.normalize(ws4ind,ws4ind,drop_empty_columns=empty_cols2drop,add_pvalue=False)
        self.report_pandas[ws4ind].tab_format['tab_color'] = 'blue'

      if not df4tox.empty:
        count_df = self.make_count_df(df4tox,ws4tox)
        self.add2raw(count_df)
        self.normalize(ws4tox,ws4tox,drop_empty_columns=empty_cols2drop,add_pvalue=False)
        self.report_pandas[ws4tox].tab_format['tab_color'] = 'red'
      return
    

    def add_target_columns(self):
      for ws in self.ws_names():
        if ws in self.report_pandas.keys():
          self.report_pandas[ws],_t2iG = self.add_target_column(self.report_pandas[ws])
          snippets_ws = ws[0:22]+'_snippets'
          snippets_df = _t2iG.snippets2df(snippets_ws,ref_sentence_props=[SENTENCE])
          self.add2report(snippets_df)
      return


    def _perform_semantic_search(self):
      '''
      output:
        raw_data[ws4ind],raw_data[ws4tox], report_pandas[ws4ind],report_pandas[ws4tox]
      '''
      start = time.time()
      # cannot multithread here yet - self.Graph is mutating.  Need to have self.clone function
      self.raw_data[RAWDF4ANTAGONISTS] = self.semscore4(self.__targets__,'positive',self.__partners__,self.raw_data[RAWDF4ANTAGONISTS])
      self.raw_data[RAWDF4AGONISTS] = self.semscore4(self.__targets__,'negative',self.__partners__,self.raw_data[RAWDF4AGONISTS])
      self.raw2report()
      self.add_target_columns()
      print(f'TargetIndications semantic search is finished in {execution_time(start)}')
      return
    
    def add_ontologies(self):
      #dfnames_map = self.__dfnames_map(self.params['mode_of_action'])
      with ThreadPoolExecutor(max_workers=5, thread_name_prefix='AddAnnot') as b:
        id_path_futures = list()
        for ws in self.ws_names():
          if ws in self.report_pandas.keys():
            ranked_df = self.report_pandas[ws]
            b.submit(self.add_ontology_df,ranked_df)
            id_path_futures.append((ws,b.submit(self.id2paths,ranked_df)))
            
        for ws, future in id_path_futures:
          id2paths = future.result()
          self.report_pandas[ws] = self.report_pandas[ws].merge_dict(id2paths,'Ontology parents','Name')
        b.shutdown()



    def other_effects(self)->ResnetGraph:
        # need to be called after ranking to subtract self.all_entity_ids
        print('Findind indication linked with unknown effect to %s' % self.target_names_str())
        old_rel_props = self.relProps
        self.add_rel_props(PS_SENTENCE_PROPS+list(PS_BIBLIO_PROPS))
        t_n = self.target_names_str()
        
        request_name = f'Find indications modulated by {t_n} with unknown effect'
        oql4indications,_ = self.oql4indications()
        OQLquery = f'SELECT Relation WHERE objectType = (Regulation,QuantitativeChange) AND Effect = unknown AND \
NeighborOf({self.oql4targets}) AND NeighborOf ({oql4indications})'
        to_return = self.process_oql(OQLquery,request_name)
        
        request_name = f'Find indications where {t_n} was suggested as Biomarker'
        OQLquery = f'SELECT Relation WHERE objectType = (Biomarker,StateChange,GeneticChange) AND \
NeighborOf({self.oql4targets}) AND NeighborOf ({oql4indications})'
        if isinstance(to_return,ResnetGraph):
            bm_g = self.process_oql(OQLquery,request_name)
            if isinstance(bm_g,ResnetGraph):
                to_return = bm_g.compose(to_return)
        
        request_name = f'Find indications with genetically linked {t_n} Genetic Variants'
        OQLquery = f'SELECT Relation WHERE objectType = FunctionalAssociation AND NeighborOf ({oql4indications})'
        OQLquery += ' AND NeighborOf (SELECT Entity WHERE id = ({ids}))'
        GVdbids = set(ResnetGraph.dbids(self.__GVs__))
        if isinstance(to_return,ResnetGraph):
            gv_g = self.iterate_oql(OQLquery,GVdbids,request_name=request_name)
            to_return = gv_g.compose(to_return)

        assert(isinstance(to_return,ResnetGraph))
        ranked_indication_ids = self.Graph.dbids4nodes(self.params['indication_types'])
        to_return.remove_nodes_from(ranked_indication_ids) # now delete indications with known effect
        indications = to_return.psobjs_with(only_with_values=self.params['indication_types'])
        print('Found %d new indications linked to %s with unknown effect' % (len(indications),self.target_names_str()))
        
        self.relProps = old_rel_props
        return to_return


    def add_graph_bibliography(self,suffix=''):
      """
      adds:
        df with PS_BIBLIOGRAPHY-suffix name to self.report_pandas
      """
      target_concepts = list(self.__targets__) + list(self.__GVs__)+list(self.__partners__)
      antagonistsG = self.Graph.get_subgraph(target_concepts,self.__indications4antagonists__)
      agonistsG = self.Graph.get_subgraph(target_concepts,self.__indications4agonists__)
      ws_ind, ws_tox = self.ws_names()
      ind_refs = f'EBKGrefs4{ws_ind[:22]}'
      tox_refs = f'EBKGrefs4{ws_tox[:22]}'

      if self.__moa() == ANTAGONIST:
        refs4antagonist_df = antagonistsG.bibliography(ind_refs+suffix)
        refs4agonist_df = agonistsG.bibliography(tox_refs+suffix)
      else:
        refs4antagonist_df = antagonistsG.bibliography(tox_refs+suffix)
        refs4agonist_df = agonistsG.bibliography(ind_refs+suffix)
      
      self.add2report(refs4antagonist_df)
      self.add2report(refs4agonist_df)
      return


    def other_effects2df(self):
        '''
        Adds
        ----
        __unknown_effect_indications__ worksheet to self.report_data with snippets for unknown effect indications
        '''
        other_effects_graph = self.other_effects()
        other_indications = other_effects_graph.snippets2df(df_name=UNKEFFECTDF)
        self.add2report(other_indications)


    def add_target_column(self,_2df:df,indication_col='Name'):
      my_copy = _2df.dfcopy()
      my_copy[TOTAL_REFCOUNT] = pd.Series(dtype='int')
      my_copy[RELEVANT_CONCEPTS] = pd.Series(dtype='str')
      t2iG = ResnetGraph()
      target_concepts = list(self.__targets__) + list(self.__GVs__)+list(self.__partners__)
      for idx in my_copy.index:
        ind_name = str(_2df.at[idx,indication_col])
        indications = self.Graph._psobjs_with(ind_name)
        if indications:
          indications += indications[0].childs()
          target_ref = list()
          all_refs = set()
          for target in target_concepts:
            sub_graph = self.Graph.get_subgraph(indications,[target])
            t2iG = t2iG.compose(sub_graph)
            t2i_refs = sub_graph.load_references()
            if t2i_refs:
              target_ref.append(f'{target.name()} ({str(len(t2i_refs))})')
              all_refs.update(t2i_refs)
          if target_ref:
            my_copy.loc[idx,TOTAL_REFCOUNT] = len(all_refs)
            my_copy.loc[idx,RELEVANT_CONCEPTS] = ','.join(target_ref)
        else:
          assert(ind_name == 'WEIGHTS:')
      
      return my_copy,t2iG


    def make_report(self):
      start_time = time.time()
      self.flush_dump_files()
      all_indications = self.indications4targets()
      
      if self.params['add_bibliography']:
        indication_names = [n.name() for n in all_indications]
        df4etm = df.from_dict({'Name':indication_names})
        target_names = self.target_names()
        df4etm._name_ = f'Targets4 {target_names}'
        etm_other = ThreadPoolExecutor(thread_name_prefix='ETMother')
        etm_future = etm_other.submit(self.RefStats.reflinks,df4etm,'Name',target_names)
        other_effects_future = etm_other.submit(self.other_effects2df)
      
      if self.init_semantic_search():
        self.perform_semantic_search()
        
        if self.params['add_bibliography']:
          names2reflinks = etm_future.result()
          refcount_col = self.tm_refcount_colname('Name',self.target_names())
          for ws in self.ws_names():
            if ws in self.report_pandas.keys():
              self.report_pandas[ws][refcount_col] = self.report_pandas[ws]['Name'].map(names2reflinks)
              self.report_pandas[ws] = self.RefStats.add_reflinks(names2reflinks,refcount_col,self.report_pandas[ws],'Name')
          self.add_tm_bibliography_df()
          self.add_graph_bibliography()
          other_effects_future.result()

        print(f'{self.report_path()} repurposing is done in {execution_time(start_time)}')
      else:
          print('Failed to initialize semantic search for TargetIndications')

            
    def write_report(self,extension='.xlsx'):
      report = pd.ExcelWriter(self.report_path(extension), engine='xlsxwriter')
      ordered_worksheets = list()
      for ws in self.ws_names():
        if ws in self.report_pandas.keys():
          ordered_worksheets.append(ws)
          ordered_worksheets.append('ontology4'+ws)
      ordered_worksheets.append(UNKEFFECTDF)
      other_ws = [ x for x in self.report_pandas.keys() if x not in ordered_worksheets]
      ordered_worksheets += other_ws
      self.add2writer(report,df_names=ordered_worksheets)
      report.close()

      raw_data = pd.ExcelWriter(self.report_path(extension='_raw.xlsx'), engine='xlsxwriter')
      self.addraw2writer(raw_data)
      raw_data.close()
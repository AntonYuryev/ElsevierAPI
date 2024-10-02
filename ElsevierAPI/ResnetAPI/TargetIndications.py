from .SemanticSearch import SemanticSearch,len,df,pd
from .ResnetAPISession import REFERENCE_IDENTIFIERS,BIBLIO_PROPERTIES,NO_REL_PROPERTIES
from .ResnetGraph import REFCOUNT, PSObject, ResnetGraph, OBJECT_TYPE
from ..ETM_API.references import PS_SENTENCE_PROPS,EFFECT,PS_BIBLIO_PROPS
from .PathwayStudioGOQL import OQL
from .FolderContent import FolderContent,PSPathway
from concurrent.futures import ThreadPoolExecutor,as_completed
from ..utils import execution_time
import time,os


RANK_SUGGESTED_INDICATIONS = True 
# in RANK_SUGGESTED_INDICATIONS mode only indications suggested in the lietarure are ranked by amount of supporting evidence
PREDICT_RANK_INDICATIONS = False
# mode also predicts and ranks indications from diseases having input target as biomarker

RAWDF4ANTAGONISTS = 'antagonistInd.'
RAWDF4AGONISTS = 'agonistInd.'

ANTAGONISTSDF = 'indications4antagonists'
AGONISTSDF = 'indications4agonists'
UNKEFFECTDF = 'possibilities'
TOXICITIES = 'TargetToxicities'

ANTAGONIST = -1 # drug inhibits its targets
AGONIST = 1 # drug activates its targets
ANY_MOA = 0

class Indications4targets(SemanticSearch):
    pass
    max_threads4ontology = 10
    max_cell_indications = 1000 # to limit number indications for popular ligands like TNF
  
    def __init__(self, *args, **kwargs):
        '''
        Input
        -----
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
                 'what2retrieve':BIBLIO_PROPERTIES,
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
        self.__targets = set()
        self.__partners = set()
        self.__TargetSecretingCells = set()
        self.targetGVs = list()

        self.__indications4antagonists = set()
        self.__indications4agonists = set()
        self.unknown_effect_indications = set()

        self.__IndirectAgonists = set() # name reflects action on target
        self.__DirectAgonists = set() # activate input targets
        self.__IndirectAntagonists = set() # inhibit input targets
        self.__DirectAntagonists = set()

        self.targets_have_weights = False

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
        return [x['Name'][0] for x in self.__targets]


    def target_names_str(self):
        """
        output:
            comma-separted string with self.target names as they appear in database.
            return may be different from self.params['target_names']
        """
        return ','.join([t.name() for t in self.__targets]) if len(self.__targets) < 3 else 'targets'
    
    
    def __dfnames_map(self):
        '''
        Return
        ------
        {rawDFname:reportDFname}
        '''
        if self.params['mode_of_action'] == ANTAGONIST:
            return {RAWDF4ANTAGONISTS:ANTAGONISTSDF, RAWDF4AGONISTS:TOXICITIES}
        elif self.params['mode_of_action'] == AGONIST:
            return {RAWDF4ANTAGONISTS:TOXICITIES, RAWDF4AGONISTS:AGONISTSDF}
        else:
            assert(self.params['mode_of_action'] == ANY_MOA)
            return {RAWDF4ANTAGONISTS:ANTAGONISTSDF, RAWDF4AGONISTS:AGONISTSDF}


    def __known_effect_indications(self):
        known_effect_indication = self.__indications4agonists | self.__indications4antagonists
        return known_effect_indication
    

    def __moa(self):
        return self.params['mode_of_action']


    def __my_indications(self):
        """
        Returns
        -------
        all found indications depending on self.params['strict_mode']
        """
        return self.__indications4antagonists | self.__indications4agonists


    def _is_strict(self):
        return self.params['strict_mode']
    

    def _oql4indications_type(self):
        indication_types_str = ','.join(self.params['indication_types'])
        return f'SELECT Entity WHERE objectType = ({indication_types_str})', indication_types_str
  

    def __load_targets(self)->set[PSObject]:
        assert(self.params['target_names'])
        i_t_n = self.params['target_names']

        self.add_ent_props(['Class'])
        oql = OQL.get_childs(i_t_n,['Name'],include_parents=True)
        targetsG = self.process_oql(oql,f'Find {i_t_n}')
        self.entProps.remove('Class')
        return set(targetsG._get_nodes())
    

    def _load_children4targets(self,targets:list[PSObject]):
        my_target_names = ResnetGraph.names(targets)
        oql = OQL.get_childs(my_target_names,['Name'],include_parents=True)
        twc = self.process_oql(oql,'Loading children4targets') # twc - target_with_children
        return twc._get_nodes()


    def __set_targets(self):
        '''
        input:
            self.params['target_names']
        output:
            self.__targets
        '''
        self.__targets = self.__load_targets()
        if self.__targets:
            self.__targets.update(self._load_children4targets(self.__targets))
            self.oql4targets = OQL.get_objects(ResnetGraph.dbids(self.__targets))
        else:
            print(f'No targets were found for {self.params['target_names']}')


    @staticmethod
    def _partner_class(target_class:str):
        if target_class == 'Ligand': return "Receptor"
        elif target_class == 'Receptor': return "Ligand"
        else: return ''


    def find_partners(self,_4targets:list[PSObject])->list[PSObject]:
        """
        Finds
        -----
        Ligands if targets are Receptor, Receptors if targets are Ligand linked to target(s) with Effect=positive\n
        dose not work for any other target class
        output:
            self.partnets
        """
        allowed_target_clases = ['Ligand','Receptor']
        all_partners = set()
        for clas in allowed_target_clases:
            class_targets = [o for o in _4targets if o.get_prop('Class') == clas]
            if class_targets:
                partner_class = self._partner_class(clas)
                select_targets = OQL.get_objects(ResnetGraph.dbids(class_targets))
                target_names = OQL.get_objects(ResnetGraph.names(class_targets))

                REQUEST_NAME = f'Find {partner_class}s for {target_names}'
                SELECTpartners = f'SELECT Entity WHERE Class = {partner_class} AND Connected by (SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive) to ({select_targets})'
                partners_graph = self.process_oql(SELECTpartners,REQUEST_NAME)
                if isinstance(partners_graph,ResnetGraph):
                    if partner_class == 'Ligand':
                        # request to find additional secreted molecules not annotated with Class=Ligand
                        SELECTsecretedpartners = 'SELECT Entity WHERE "Cell Localization" = Secreted AND objectType = Protein AND Connected by (SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive) to ({select_target})'
                        secreted_partners = self.process_oql(SELECTsecretedpartners.format(select_target=self.oql4targets),f'Find secreted partners for {target_names}')
                        if isinstance(secreted_partners,ResnetGraph):
                            partners_graph = secreted_partners.compose(partners_graph)

                    partners = partners_graph._get_nodes()
                    partners_names = ','.join(ResnetGraph.names(partners))
                    print(f'Found {partners_names} {clas}s as partners for {target_names}')
                    all_partners.update(partners)
        
        return all_partners
    

    def params2partners(self, partner_params:dict,_4targets:list[PSObject])->list[PSObject]:
        '''
        case when self.params has 'partner_names' and 'partner_class' explicitly\n
        use it when receptor partners are metabolites and therefore cannot be found by select_partners()\n
        input:
            partner_params - {target_name:partner_class:[partner_names]}
        '''
        _4target_names = {n.name():n for n in _4targets}
        partner_names = dict() #{name:(class,weight)}
        for target_name, c2p in partner_params.items():
            if target_name in _4target_names:
                assert(isinstance(c2p,dict))
                for p_class, p_names in c2p.items():
                    for name in p_names:
                        partner_weight = _4target_names[target_name].get_prop('target weight',if_missing_return=0)
                        partner_names[str(name).lower()] = (p_class,partner_weight)

        partners_names_s = OQL.join_with_quotes(list(partner_names.keys()))
        SELECTpartners = f'SELECT Entity WHERE Name = ({partners_names_s})'
        req_name = f'Finding partners4{_4target_names.keys()}'
        partners_graph = self.process_oql(SELECTpartners,req_name)
        if isinstance(partners_graph,ResnetGraph):
            partners = partners_graph._get_nodes()
            for p in partners:
                p_class, p_weight = partner_names[p.name().lower()]
                p['Class'] = [p_class]
                p['target weight'] = [p_weight]
                p['regulator weight'] = [p_weight]
            return partners
        else:
            return []


    def set_partners(self,_4targets:list[PSObject]):
        """
        input:
            self.params['partners'] = {target_name:{partner_class:[partner_names]}}
            Assumes partners are linked to Targets with Effect=positive
        output:
            self.__partners
            self.partner_class - is used for column names only
            self.find_partners_oql
        """
        try:
            partner_params = dict(self.params['partners'])
            self.__partners = self.params2partners(partner_params,_4targets)
        except KeyError:
            self.partner_class = ResnetGraph.classes(self.__targets)
            if self.partner_class:
                self.__partners = self.find_partners()
        return
    

    def _partner_class_str(self):
        clases = {ResnetGraph.classes(self.__partners)}
        return ','.join(clases)

        
    def report_path(self, extension='.xlsx'):
        indics = ','.join(self.params['indication_types'])
        rep_pred = 'suggested ' if self.params['strict_mode'] else 'suggested,predicted '
        if self.params['mode_of_action'] == ANTAGONIST:
            mode = ' inhibition'
        elif self.params['mode_of_action'] == AGONIST:
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
            print('Found %d cell types producing %s' % (len(self.__TargetSecretingCells),t_n))
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
            self.__indications4antagonists
            self.__indications4agonists
        '''
        kwargs2findmodulators = {
                        'of_types':['SmallMol'],
                        'with_effect' : 'negative',
                        'on_targets' : self.__targets
                      }

        if self.__moa() in [ANTAGONIST,ANY_MOA]:
            kwargs2findmodulators['linked_by'] = ['DirectRegulation']
            self.__DirectAntagonists = self.find_modulators(**kwargs2findmodulators)
            kwargs2findmodulators['linked_by'] = ['Regulation','Expression','MolTransport']
            kwargs2findmodulators['drugs_only'] = True
            self.__IndirectAntagonists = self.find_modulators(**kwargs2findmodulators)
            
        if self.__moa() in [AGONIST,ANY_MOA]:
            kwargs2findmodulators['with_effect'] = 'positive'
            kwargs2findmodulators['linked_by'] = ['DirectRegulation']
            self.__DirectAgonists = self.find_modulators(**kwargs2findmodulators)
            kwargs2findmodulators['linked_by'] = ['Regulation','Expression','MolTransport']
            kwargs2findmodulators['drugs_only'] = True
            self.__IndirectAgonists = self.find_modulators(**kwargs2findmodulators)


        if not self._is_strict():
            if self.__DirectAntagonists:
                self.__indications4antagonists.update(self.find_indications4(self.__DirectAntagonists))
                self.__indications4agonists.update(self.find_toxicities4(self.__DirectAntagonists))
            if self.__DirectAgonists:
                self.__indications4agonists.update(self.find_indications4(self.__DirectAgonists))
                self.__indications4antagonists.update(self.find_toxicities4(self.__DirectAgonists))
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
        if self.__partners:
            for partner in self.__partners:
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
        #effect = 'positive' if moa == ANTAGONIST else 'negative'
        t_n = ResnetGraph.names(targets)

        REQUEST_NAME = 'Find indications {effect}ly regulating {target}'.format(effect=with_effect_on_indications, target=t_n)
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = {effect} AND \
            NeighborOf({target}) AND NeighborOf ({indications})' 
            
        f_t = OQL.get_objects(ResnetGraph.dbids(targets))
        f_i,_ = self._oql4indications_type()
        indications = set()

        OQLquery =  OQLquery.format(target=f_t,effect = with_effect_on_indications,indications=f_i)  
        ModulatedByTargetNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        if isinstance(ModulatedByTargetNetwork,ResnetGraph):
            indications = set(ModulatedByTargetNetwork.psobjs_with(only_with_values=self.params['indication_types']))
            print('Found %d diseases %sly regulated by %s' % (len(indications),with_effect_on_indications,t_n))

        if not self._is_strict():
            REQUEST_NAME = f'Find indications {with_effect_on_indications}ly modulated by {t_n}'
            OQLquery = 'SELECT Relation WHERE objectType = QuantitativeChange AND Effect = {effect} AND \
                NeighborOf ({target}) AND NeighborOf ({indications})'
            OQLquery = OQLquery.format(target=f_t,effect = with_effect_on_indications,indications=f_i)
            ActivatedInDiseaseNetwork = self.process_oql(OQLquery,REQUEST_NAME)
            if isinstance(ActivatedInDiseaseNetwork,ResnetGraph):
                add2indications = ActivatedInDiseaseNetwork.psobjs_with(only_with_values=self.params['indication_types'])
                indications.update(add2indications)
                print('Found %d diseases where %s is %sly regulated' % (len(add2indications),t_n,with_effect_on_indications))
        return indications


    def __indications4targets(self):
        '''
        input:
            if _4targets is empty will use target names from self.__targets
        output:
            self.__indications4antagonists
            self.__indications4agonists
            self.unknown_effect_indications
        '''
        moa = self.params['mode_of_action']
        
        if moa in [ANTAGONIST,ANY_MOA]:
            effect_on_indication = 'positive'
            self.__indications4antagonists = self._indications4(self.__targets,effect_on_indication)

        if moa in [AGONIST,ANY_MOA]:
            effect_on_indication = 'negative'
            self.__indications4agonists = self._indications4(self.__targets,effect_on_indication)

        self._unknown_effect_indications(self.__targets,self.__known_effect_indications())
        return


    def _indications4partners(self,targets:list[PSObject],partners:list[PSObject],effect_on_indications:str):
        OQLtemplate = 'SELECT Relation WHERE objectType = (Regulation,QuantitativeChange) AND \
Effect = {eff} AND NeighborOf ({partners}) AND NeighborOf ({indications})'
        assert(effect_on_indications in ['positive','negative'])

        #effect = 'positive' if moa == ANTAGONIST else 'negative'
        t_n = ','.join(ResnetGraph.names(targets))
        p_c = ','.join(ResnetGraph.classes(partners))
        oql4indications,_ = self._oql4indications_type()
        oql4partners = OQL.get_objects(ResnetGraph.dbids(partners))

        OQLquery = OQLtemplate.format(eff=effect_on_indications,partners=oql4partners,indications=oql4indications)
        REQUEST_NAME = f'Find indications {effect_on_indications}ly regulated by {p_c}s of {t_n}'
        PartnerIndicationNetwork4anatagonists = self.process_oql(OQLquery,request_name=REQUEST_NAME)
        indications = PartnerIndicationNetwork4anatagonists.psobjs_with(only_with_values=self.params['indication_types'])
        print(f'Found {len(indications)} indications for {len(partners)} {t_n} {p_c}')
        return indications


    def __indications4partners(self):
        """
        input:
            self.params['mode_of_action']
            self.__partners
        Assumes partners are linked to Target with Effect=positive
        """
        t_n = self.target_names_str()
        exist_indication_count = len(self.__my_indications())
        moa = self.params['mode_of_action']

        if moa in [ANTAGONIST,ANY_MOA]:
            effect_on_indications = 'positive'
            partner_indications = self._indications4partners(self.__targets,self.__partners,effect_on_indications)
            self.__indications4agonists.update(partner_indications)
        if moa in [AGONIST,ANY_MOA]:
            effect_on_indications = 'negative'
            partner_indications = self._indications4partners(self.__targets,self.__partners,effect_on_indications)
            self.__indications4agonists.update(partner_indications)
  
        new_indication_count = len(self.__my_indications()) - exist_indication_count
        print('%d indications for %d %s partnerss were not found by previous searches' %  
                (new_indication_count, len(self.__partners),t_n))
        return
    

    def find_indications4(self,modulators:list[PSObject])->set[PSObject]:
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


    def find_toxicities4(self,modulators:list[PSObject])->set[PSObject]:
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


    def __indications4cells(self,secreting:list[PSObject],targets:list[PSObject],with_effect_on_indication:str):
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
        REQUEST_NAME = f'Find indications {with_effect_on_indication}ely linked to cell secreting {t_n}'
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


    def _indications4cells(self,secreting:set[PSObject],targets:list[PSObject],effect_on_targets:str)->list[PSObject]:
        '''
        ouput:
            indications, 
            if "secreting" is empty adds cells to "secreting"
        '''
        assert(effect_on_targets in ['positive','negative'])
        my_ligands = [o for o in targets if o.get_prop('Class') == 'Ligand']
        if my_ligands:
            if not secreting:
                secreting.update(self.__find_cells_secreting(my_ligands))

            indications = list()
            if secreting:
                indications = self.__indications4cells(secreting,my_ligands,effect_on_targets)
            return indications
        else:
            print('Target list contains no ligands. No secreting cell can be found')
            return []
        

    def indications4secreting_cells(self):
        if self._is_strict(): 
            print('Scipt runs in strict mode.  No indications for cells secreting targets will be used')
            self.__TargetSecretingCells = self.__find_cells_secreting(self.__targets)
            # loaded here self.__TargetSecretingCells is used for target ranking
            return
        
        exist_indication_count = len(self.__my_indications())
        if self.__moa() in [ANTAGONIST,ANY_MOA]:
            indications,secreting_cells = self._indications4cells(self.__TargetSecretingCells,self.__targets,ANTAGONIST)
            self.__TargetSecretingCells.update(secreting_cells)
            self.__indications4antagonists.update(indications)

        if self.__moa() in [AGONIST,ANY_MOA]:
            indications, secreting_cells = self._indications4cells(self.__TargetSecretingCells,self.__targets,AGONIST)
            self.__TargetSecretingCells.update(secreting_cells)
            self.__indications4agonists.update(indications)
            
        new_indication_count = len(self.__my_indications()) - exist_indication_count
        print('%d indications for %d cells producing %s were not found by previous searches' %  
                (new_indication_count, len(self.__TargetSecretingCells),self.target_names()))


    def __resolve_conflict_indications(self):
        def __resolve(conflict:PSObject, all_rels:list, using_rel_type:str):
            only_rels_with_type = [r for r in all_rels if r.objtype() == using_rel_type]
            if only_rels_with_type:
                only_rels_with_type.sort(key=lambda x: x.count_refs(), reverse=True)
                best_effect = only_rels_with_type[0]['Effect'][0]
                if best_effect == 'positive':
                    self.__indications4agonists.discard(conflict)
                    #self.__indications4agonists_strict.discard(conflict)
                    return True
                elif best_effect == 'negative':
                    self.__indications4antagonists.discard(conflict)
                    #self.__indications4antagonists_strict.discard(conflict)
                    return True
                else:
                    return False
            return False

        conflict_indications = self.__indications4antagonists.intersection(self.__indications4agonists)
        unique_indications = len(self.__indications4antagonists)+len(self.__indications4agonists)-len(conflict_indications)
        print(f'Resolving {len(conflict_indications)} conflicting out of total {unique_indications} indications')
        unresolved_ids = set()
        target_uids = ResnetGraph.uids(self.__targets)
        partner_uids = ResnetGraph.uids(self.__partners)
        for conflict in conflict_indications:
            target_rels = list(self.Graph.find_relations(target_uids,[conflict.uid()]))
            if not __resolve(conflict,target_rels,'Regulation'):
                if not __resolve(conflict,target_rels,'QuantitativeChange'):
                    partner_rels = list(self.Graph.find_relations(partner_uids,[conflict.uid()]))
                    if not __resolve(conflict,partner_rels,'Regulation'):
                        if not __resolve(conflict,partner_rels,'QuantitativeChange'):
                            unresolved_ids.add(conflict.uid())

        print('%d indications cannot be resolved.\n They will appear in both worksheets for agonists and antagonist:' % 
                    len(unresolved_ids))
        for uid in unresolved_ids:
            try:
                print(self.Graph._psobj(uid).name())
            except KeyError:
                continue


    def indications4targets(self):
        '''
        Input
        -----
        target names must be in self.params['target_names']
        moa in [ANTAGONIST, AGONIST, ANY_MOA]
        '''
        #assert(moa in [ANTAGONIST, AGONIST, ANY_MOA])
        start_time = time.time()
        self.__set_targets()
        self.set_partners()
        self.__indications4targets()
        self.__indications4partners()

        self.__resolve_conflict_indications()
        print("%d indications for %s were found in %s" % 
        (len(self.__my_indications()), self.target_names_str(), execution_time(start_time)))

        return self.__my_indications()

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
        

    def GVindications(self)->list[PSObject]:
        self.GVs2DiseaseGraph = self._GVindicationsG(self.__targets)
        if isinstance(self.GVs2DiseaseGraph,ResnetGraph):
            indications = self.GVs2DiseaseGraph.psobjs_with(only_with_values=self.params['indication_types'])
            self.targetGVs = self.GVs2DiseaseGraph._psobjs_with('GeneticVariant',OBJECT_TYPE)

            t_n = self.target_names_str()
            print('Found %d indications genetically linked to %d Genetic Variants in %s' % 
                (len(indications), len(self.targetGVs), t_n))
            return indications
        else:
            return list()


    def _biomarker_indicationsG(self,_4targets:list[PSObject]):
        t_n = self.target_names_str()
        REQUEST_NAME = f'Find indications where {t_n} is biomarker'
        OQLquery = 'SELECT Relation WHERE objectType = Biomarker AND NeighborOf({select_target}) AND NeighborOf ({indications})'

        f_t = OQL.get_objects(ResnetGraph.dbids(_4targets))
        oql4indications,_ = self._oql4indications_type()
        OQLquery = OQLquery.format(select_target=f_t,indications=oql4indications)
        BiomarkerInDiseaseNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        if isinstance(BiomarkerInDiseaseNetwork,ResnetGraph):
            return BiomarkerInDiseaseNetwork
        else:
            return ResnetGraph()
         

    def _unknown_effect_indications(self,known_indications:list[PSObject]=[]):
        t_n = self.target_names_str()
        gv_indications = set(self.GVindications())
        print('%d indications linked to %s genetic variations' % (len(indications),t_n))

        BiomarkerInDiseaseNetwork = self._biomarker_indicationsG(self.__targets)
        biomarker_indications = set(BiomarkerInDiseaseNetwork.psobjs_with(only_with_values=self.params['indication_types']))
        print('Found %d indications where target %s is claimed as biomarker' %  (len(biomarker_indications),t_n))
        
        indications = gv_indications|biomarker_indications
        indications = indications.difference(known_indications)
        print(f'Found {len(indications)} indications linked with unknown effect to targets')
        self.unknown_effect_indications = indications
        return indications

        
##################  SCORING SCORE SCORING ######################### SCORING SCORE SCORING ####################
    def score_GVs(self, df2score:df):
        if not self.targetGVs: 
            return
        t_n = self.target_names_str()
        concept_name = 'GVs'
        self.__colname4GV__ = t_n+concept_name
        refcount_column = self._refcount_colname(concept_name)
        weighted_refcount_column = self._weighted_refcount_colname(concept_name) 
        linked_count_column = self._linkedconcepts_colname(concept_name)
        concept_size_column = self._concept_size_colname(concept_name)

        gvlinkcounter = 0
        if hasattr(self,'GVs2DiseaseGraph'):
            if self.targets_have_weights:
                regurn2weight = {gv.urn():gv.get_prop('regulator weight') for gv in self.targetGVs}
            
            totalGVcount = len(self.targetGVs)
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
                    rowGVs = list(self.GVs2DiseaseGraph.get_neighbors(row_indications,self.targetGVs))
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
                
            print('Found %d indications linked to %d GVs' % (gvlinkcounter, len(self.targetGVs)))


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
        my_df = df.copy_df(in_df)
        score4antagonists = True if with_effect_on_indication == 'positive' else False
        if self._is_strict():
            booster_reltypes = ['Regulation','Biomarker','GeneticChange','QuantitativeChange','StateChange','FunctionalAssociation']
        else:
            booster_reltypes = ['Regulation','GeneticChange','FunctionalAssociation']
        
        t_n = self.target_names_str()
        target_in_header = t_n if len(t_n) < 45 else 'targets'
        if target_in_header == 'targets':
             colname = 'Regulated by ' + target_in_header
        elif with_effect_on_indication == 'positive':
            colname = 'Activated by '+target_in_header
        else:
            colname = 'Inhibited by '+target_in_header
        
        kwargs = {'connect_by_rels':['Regulation'],
                  'with_effects' : [with_effect_on_indication],
                  'boost_by_reltypes' :booster_reltypes
                  }
        if self.targets_have_weights:
            kwargs.update({'nodeweight_prop': 'regulator weight'})
        how2connect = self.set_how2connect(**kwargs)
        linked_row_count,_,my_df = self.link2concept(colname,targets,my_df,how2connect)
        print('%d indications are %sly regulated by %s' % 
            (linked_row_count,with_effect_on_indication,t_n))
        self.nodeweight_prop = ''

        self.score_GVs(my_df)

        link_effect, drug_class, concepts = self._drug_connect_params(True,score4antagonists)
        if concepts:
            # references suggesting that known drugs for the target as treatments for indication
            colname = target_in_header+' '+drug_class+' clin. trials'
            colname = target_in_header+' '+drug_class
            kwargs = {'connect_by_rels':['ClinicalTrial']}
            if self.targets_have_weights:
                kwargs.update({'nodeweight_prop': 'regulator weight'})
            how2connect = self.set_how2connect(**kwargs)
            linked_row_count,_,my_df = self.link2concept(colname,concepts,my_df,how2connect)
            print('Linked %d clinical trial indictions for %s %s' % 
                  (linked_row_count,t_n,drug_class))

            colname = target_in_header+' '+drug_class
            kwargs = {'connect_by_rels':['Regulation'],
                  'with_effects' : [link_effect],
                  'boost_by_reltypes' :['Regulation','FunctionalAssociation']
                  }
            how2connect = self.set_how2connect(**kwargs)
            linked_row_count,_,my_df = self.link2concept(colname,concepts,my_df,how2connect)
            print('Linked %d indications for %s %s' % 
                  (linked_row_count,t_n,drug_class))

        #references reporting target agonists exacerbating indication or causing indication as adverse events
        link_effect, drug_class, concepts = self._drug_tox_params(direct_modulators=True,
                                                                   score_antagonists=score4antagonists)
        if concepts:
            colname = target_in_header+' '+drug_class
            kwargs = {'connect_by_rels':['Regulation'],
                  'with_effects' : [link_effect],
                  'boost_by_reltypes' :['Regulation','FunctionalAssociation']
                  }
            how2connect = self.set_how2connect(**kwargs)
            linked_row_count,_,my_df = self.link2concept(colname,concepts,my_df,how2connect)
            print('Linked %d indications as toxicities for %s %s' % 
                  (linked_row_count,t_n,drug_class))

        #references where target expression or activity changes in the indication
        if target_in_header == 'targets':
             indication_types = ','.join(self.params['indication_types'])
             colname = target_in_header + f' changes in {indication_types}'
        elif with_effect_on_indication == 'positive':
            colname = target_in_header+' is upregulated'
        else:
            colname = target_in_header+' is downregulated'


        kwargs = {'connect_by_rels':['QuantitativeChange'],
                  'with_effects' : [with_effect_on_indication],
                  'boost_by_reltypes' : ['Biomarker','StateChange','FunctionalAssociation'],
                  }
        if self.targets_have_weights:
            kwargs.update({'nodeweight_prop': 'target weight'})
        how2connect = self.set_how2connect(**kwargs)
        linked_row_count,_,my_df= self.link2concept(colname,targets,my_df,how2connect)
        print('%d indications %sly regulate %s' % 
            (linked_row_count,with_effect_on_indication,t_n))
        self.nodeweight_prop = ''

        #references suggesting target partners as targets for indication
        if with_partners:
            colname = f'{target_in_header} partners'
            kwargs = {'connect_by_rels':['Regulation'],
                  'with_effects' : [with_effect_on_indication],
                  'boost_by_reltypes' : ['Regulation','FunctionalAssociation']
                  }
            if self.targets_have_weights:
                 kwargs.update({'nodeweight_prop': 'regulator weight'})
            how2connect = self.set_how2connect(**kwargs)
            linked_row_count,_,my_df = self.link2concept(colname,with_partners,my_df,how2connect)
            print('Linked %d indications for %d %s partners' % 
                (linked_row_count,len(with_partners),t_n))
            self.nodeweight_prop = ''

        # references reporting that cells producing the target linked to indication  
        # only used if taregts are secretred ligands
        if hasattr(self, '__TargetSecretingCells'):
            colname = f'{target_in_header} secreting cells'
            kwargs = {'connect_by_rels':['Regulation'],
                  'with_effects' : [with_effect_on_indication],
                  'boost_by_reltypes' : ['Regulation','FunctionalAssociation']
                  }
            if self.targets_have_weights:
                kwargs.update({'nodeweight_prop': 'regulator weight'})
            how2connect = self.set_how2connect(**kwargs)
            linked_row_count,_,my_df = self.link2concept(colname,self.__TargetSecretingCells,my_df,how2connect)
            print('Liked %d indications linked %d cells producing %s' % 
                  (linked_row_count,len(self.__TargetSecretingCells),t_n))

        link_effect, drug_class, concepts = self._drug_connect_params(direct_modulators=False,
                                                                      score_antagonists=score4antagonists)
        if concepts:
            # references suggesting that known drugs for the target as treatments for indication
            colname = target_in_header+' '+drug_class+' clin. trials'
            kwargs = {'connect_by_rels':['ClinicalTrial']}
            if self.targets_have_weights:
                kwargs.update({'nodeweight_prop': 'regulator weight'})
            # cloning session to avoid adding relations to self.Graph
            new_session = self._clone(to_retrieve=REFERENCE_IDENTIFIERS)
            how2connect = new_session.set_how2connect(**kwargs)
            linked_row_count,_,my_df = new_session.link2concept(colname,concepts,my_df,how2connect)
            print('Linked %d clinical trial indications for %s %s' % 
                  (linked_row_count,t_n,drug_class))
            new_session.close_connection()


            colname = target_in_header+' '+drug_class
            kwargs = {'connect_by_rels':['Regulation'],
                  'with_effects' : [link_effect],
                  'boost_by_reltypes' : ['Regulation'] # boosting with unknown effect Regulation
                  }
            if self.targets_have_weights:
                kwargs.update({'nodeweight_prop': 'regulator weight'})
            # cloning session to avoid adding relations to self.Graph
            new_session = self._clone(to_retrieve=REFERENCE_IDENTIFIERS)
            how2connect = new_session.set_how2connect(**kwargs)
            linked_row_count,_,my_df = new_session.link2concept(colname,concepts,my_df,how2connect)
            print('Linked %d indications for %s %s' % (linked_row_count,t_n,drug_class))
            new_session.close_connection()

        #references reporting target agonists exacerbating indication or causing indication as adverse events
        link_effect, drug_class, concepts = self._drug_tox_params(direct_modulators=False,
                                                                   score_antagonists=score4antagonists)
        if concepts:
            colname = target_in_header+' '+drug_class
            kwargs = {'connect_by_rels':['Regulation'],
                  'with_effects' : [link_effect],
                  'boost_by_reltypes' : ['Regulation'] # boosting with unknown effect Regulation
                  }
            if self.targets_have_weights:
                kwargs.update({'nodeweight_prop': 'regulator weight'})
            # cloning session to avoid adding relations to self.Graph
            new_session = self._clone(to_retrieve=REFERENCE_IDENTIFIERS)
            how2connect = new_session.set_how2connect(**kwargs)
            linked_row_count,_,my_df = new_session.link2concept(colname,concepts,my_df,how2connect)
            print('Linked %d indications as toxicities for %s %s' % (linked_row_count,t_n,drug_class))
            new_session.close_connection()
        
        
        if hasattr(self, 'PathwayComponents'):
            #references linking target pathways to indication
            colname = target_in_header + ' pathway components'
            kwargs = {'connect_by_rels':['Regulation'],
                  'with_effects' : [with_effect_on_indication],
                  'step' : 125 # boosting with unknown effect Regulation
                  }
            # cloning session to avoid adding relations to self.Graph
            new_session = self._clone(to_retrieve=REFERENCE_IDENTIFIERS)
            how2connect = new_session.set_how2connect(**kwargs)
            linked_row_count,_,my_df = new_session.link2concept(colname,list(self.PathwayComponents),my_df,how2connect)
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

        if self.__indications4antagonists and self.params['mode_of_action'] in [ANTAGONIST,ANY_MOA]:
            antagonist_indication_df = self.load_df(list(self.__indications4antagonists),
                                         max_child_count=self.max_ontology_parent,
                                         max_threads=self.max_threads4ontology)
            antagonist_indication_df._name_ = RAWDF4ANTAGONISTS
            self.add2raw(antagonist_indication_df)
            print(f'Will score {len(antagonist_indication_df)} indications for {t_n} antagonists')

        if self.__indications4agonists and self.params['mode_of_action'] in [AGONIST,ANY_MOA]:
            agonist_indication_df = self.load_df(list(self.__indications4agonists),
                                         max_child_count=self.max_ontology_parent,
                                         max_threads=self.max_threads4ontology)
            agonist_indication_df._name_ = RAWDF4AGONISTS
            self.add2raw(agonist_indication_df)
            print(f'Will score {len(agonist_indication_df)} indications for {t_n} agonists')
        
        if self.__indications4antagonists or self.__indications4agonists:
            return True
        else:
            print (f'No indications found for {t_n}')
            if self._is_strict(): print('Try setting strict_mode to False')
            return False


    def perform_semantic_search(self):
        '''
        Create:
            worksheets in self.report_data: ANTAGONISTSDF,AGONISTSDF
        '''
        start = time.time()
        # cannot multithread here yet - self.Graph is mutating.  Need to have self.clone function
        self.semscore4(self.__targets,RAWDF4ANTAGONISTS,'positive')
        self.semscore4(self.__targets,RAWDF4AGONISTS,'negative')
        
        empty_cols2drop = list(self.params.get('drop_empty_columns_from_report',[]))
        dfnames_map = self.__dfnames_map()
        for raw_dfname, report_dfname in dfnames_map.items():
            raw_df = self.raw_data[raw_dfname]
            count_df = self.make_count_df(raw_df,raw_df.__name__)
            self.add2raw(count_df)
            self.normalize(raw_dfname,report_dfname,drop_empty_columns=empty_cols2drop)

        with ThreadPoolExecutor(max_workers=5, thread_name_prefix='AddAnnot') as b:
            id_path_futures = list()
            for ws in dfnames_map.values():
                if ws in self.report_pandas.keys():
                    ranked_df = self.report_pandas[ws]
                    b.submit(self.add_ontology_df,ranked_df)
                    id_path_futures.append((worksheet_name,b.submit(self.id2paths,ranked_df)))

            for worksheet_name, future in id_path_futures:
                id2paths = future.result()
                self.report_pandas[worksheet_name] = self.report_pandas[worksheet_name].merge_dict(id2paths,'Ontology parents','Name')
            
            b.shutdown()
        print(f'TargetIndications semantic search is finished in {execution_time(start)}')
        return 


    def other_effects(self)->ResnetGraph:
        # need to be called after ranking to subtract self.all_entity_ids
        print('Findind indication linked with unknown effect to %s' % self.target_names_str())
        old_rel_props = self.relProps
        self.add_rel_props(PS_SENTENCE_PROPS+list(PS_BIBLIO_PROPS))
        t_n = self.target_names_str()
        
        request_name = f'Find indications modulated by {t_n} with unknown effect'
        oql4indications,_ = self._oql4indications_type()
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
        GVdbids = set(ResnetGraph.dbids(self.targetGVs))
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


    def add_graph_bibliography(self,suffix='',add_graph=ResnetGraph()):
        """
        adds:
            df with PS_BIBLIOGRAPHY-suffix name to self.report_pandas
        """
        targets_neighbors = self.Graph.neighborhood(set(self.__targets))
        targets_neighbors = self.Graph.neighborhood(set(self.targetGVs)).compose(targets_neighbors)
        targets_neighbors = self.Graph.neighborhood(set(self.__partners)).compose(targets_neighbors)
        if add_graph:
            targets_neighbors = add_graph.compose(targets_neighbors)
        super().add_graph_bibliography(suffix,from_graph=targets_neighbors)


    def other_effects2df(self):
        '''
        Adds
        ----
        UNKNOWN_EFFECT_INDICATIONS worksheet to self.report_data with snippets for unknown effect indications
        '''
        other_effects_graph = self.other_effects()
        other_indications = other_effects_graph.snippets2df(df_name=UNKEFFECTDF)
        self.add2report(other_indications)


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
            etm_future = etm_other.submit(self.bibliography,df4etm,target_names,'Name',[],len(df4etm))
            other_effects_future = etm_other.submit(self.other_effects2df)
        
        if self.init_semantic_search():
            self.perform_semantic_search()

            if self.params['add_bibliography']:
                indication_etmrefs = etm_future.result()
                etm_ref_colname = self.refcount_column_name('Name',self.target_names())
                doi_ref_colname = self.doi_column_name('Name',self.target_names())
                for ws in self.__dfnames_map().values():
                    if ws in self.report_pandas.keys():
                        self.report_pandas[ws] = self.report_pandas[ws].merge_df(indication_etmrefs,how='left',on='Name',columns=[etm_ref_colname,doi_ref_colname])
                        self.report_pandas[ws] = self.report_pandas[ws].move_cols({etm_ref_colname:2})
                self.add_etm_bibliography()
                self.add_graph_bibliography()
                other_effects_future.result()

            print(f'{self.report_path()} repurposing is done in {execution_time(start_time)}')
        else:
            print('Failed to initialize semantic search for TargetIndications')

            
    def write_report(self):
        report = pd.ExcelWriter(self.report_path(), engine='xlsxwriter')
        ordered_worksheets = list()
        for worksheet_name in self.__dfnames_map().values():
            if worksheet_name in self.report_pandas.keys():
                ordered_worksheets.append(worksheet_name)
                ordered_worksheets.append('ontology4'+worksheet_name)
        ordered_worksheets.append(UNKEFFECTDF)
        
        self.add2writer(report,df_names=ordered_worksheets)
        report.close()


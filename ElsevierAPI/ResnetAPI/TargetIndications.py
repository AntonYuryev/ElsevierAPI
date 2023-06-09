from .SemanticSearch import SemanticSearch,ResnetGraph, len,df
from .ResnetAPISession import REFERENCE_IDENTIFIERS,BIBLIO_PROPERTIES
from .ResnetGraph import REFCOUNT, PSObject
from ..ETM_API.references import PS_SENTENCE_PROPS,EFFECT,PS_BIBLIO_PROPS
from .PathwayStudioGOQL import OQL
from concurrent.futures import ThreadPoolExecutor,as_completed
import time

RANK_SUGGESTED_INDICATIONS = True 
# in RANK_SUGGESTED_INDICATIONS mode only indications suggested in the lietarure are ranked by amount of supporting evidence
PREDICT_RANK_INDICATIONS = False
# mode also predicts and ranks indications from diseases having input target as biomarker

DF4ANTAGONISTS = 'antagonistInd.'
DF4AGONISTS = 'agonistInd.'

ANTAGONISTS_WS = 'indications4antagonists'
AGONISTS_WS = 'indications4agonists'

class Indications4targets(SemanticSearch):
    pass
    indications4antagonists = set()
    indications4agonists = set()
    indications4antagonists_strict = set()
    indications4agonists_strict = set()
    unknown_effect_indications = set()
    unknown_effect_indications_strict = set()
    partners = list() # list of PSObject
    GVids = list()    


    def __init__(self, *args, **kwargs):
        '''
        Input
        -----
        APIconfig - args[0]
        '''
        APIconfig = args[0]
        my_kwargs = {'partner_names':[],
    # if partner_names is empty script will try finding Ligands for Receptor targets and Receptors for Ligand targets
                 'partner_class':'', # use it only if partner_names not empty
                 'indication_types': ['Disease','Virus','Pathogen'], #['Disease','Virus','CellProcess']
                 'target_names':[],
                 'pathway_name_must_include_target':True,
    # if 'pathway_name_must_include_target' True only pathways depicting target signaling are considered for regulome construction
    # if 'pathway_name_must_include_target' False all pathways containing both Targets and Partners will be considered
                 'strict_mode':True,
    # if strict mode algorithm only ranks indications suggetsed in the literarure without predicting new indications
    # if not in strict mode algorithm also predicts and rank additional indications based on target expression or genetic profiles in disease
                 'data_dir':'',
                 'add_bibliography' : True,
                 'what2retrieve':BIBLIO_PROPERTIES
                }
        my_kwargs.update(kwargs)

        super().__init__(APIconfig,**my_kwargs)
        self.add_rel_props([EFFECT])
        self.columns2drop += [self.__resnet_name__,self.__mapped_by__]
        self.max_ontology_parent = 11
        self.max_threads4ontology = 10
        # 4 threads perform a bit faster than 8 threads 
        # 288 disease parents out of 288 were processed in 0:19:51.732684 by 8 threads
        # 288 disease parents out of 288 were processed in 0:18:52.715937 by 4 threads


    def _target_names(self):
        """
        Returns target names as they were entered in configuration.
        return values may be different from self.input_names()
        """
        return ','.join(self.params['target_names'])


    def input_names(self):
        """
        Returns
        -------
        entity names that can be used both in PS and ETM
        """
        #input names are used for finding references in ETM.
        # RepurposeDrug has it own 'input_names'.
        return [x['Name'][0] for x in self.targets]


    def __known_effect_indications(self):
        known_effect_indication = self.indications4agonists | self.indications4antagonists
        return known_effect_indication


    def __indications(self):
        """
        Returns
        -------
        all found indications depending on self.params['strict_mode']
        """
        if self.params['strict_mode']:
            return self.indications4antagonists_strict | self.indications4agonists_strict
        else:
            return self.indications4antagonists | self.indications4agonists


    def _indications4antagonists(self):
        if self.params['strict_mode']:
            return self.indications4antagonists_strict
        else:
            return self.indications4antagonists

    def _indications4agonists(self):
        if self.params['strict_mode']:
            return self.indications4agonists_strict
        else:
            return self.indications4agonists

    def __is_strict(self):
        return self.params['strict_mode']
  

    def set_targets(self):
        try:
            target_objtype_str = str(','.join(self.params['target_objtypes']))
        except KeyError: 
            target_objtype_str = 'Protein'

        prop_names_str, prop_values_str = OQL.get_search_strings(['Name'],self.params['target_names'])
        self.add_ent_props(['Class'])
        self.oql4targets = f'SELECT Entity WHERE ({prop_names_str}) = ({prop_values_str}) AND objectType = ({target_objtype_str})'
        targets_graph = self.process_oql(self.oql4targets)
        
        self.targets = targets_graph._get_nodes()
        target_dbids = [x['Id'][0] for x in self.targets]
        self.oql4targets = OQL.get_objects(target_dbids)
        for t in self.targets:
            try:
                self.target_class = t['Class'][0]
                break
            # assumes all targets have the same class
            except KeyError: continue


    def set_partner_class(self):
        #finding effect for complementary receptor or ligands
        if self.target_class == 'Ligand': 
            self.partner_class = "Receptor"
            return True
        elif self.target_class == 'Receptor': 
            self.partner_class = "Ligand"
            return True
        else: return False


    def find_partners(self):
        """
        Finds
        -----
        Ligands if targets are Receptor, Receptors if targets are Ligand linked to target(s) with Effect=positive\n
        dose not work for any other target class
        """
        if self.set_partner_class():
            SELECTpartners = f'SELECT Entity WHERE Class = {self.partner_class} AND Connected by (SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive) to ({self.oql4targets})'
            partners_graph = self.process_oql(SELECTpartners)
            if self.partner_class == 'Ligand':
                # additional request to find secreted molecules that are not annotated with class Ligand
                SELECTsecretedpartners = 'SELECT Entity WHERE "Cell Localization" = Secreted AND objectType = Protein AND Connected by (SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive) to ({select_target})'
                secreted_partners = self.process_oql(SELECTsecretedpartners.format(select_target=self.oql4targets))
                partners_graph = partners_graph.add_graph(secreted_partners)

            self.partners = partners_graph._get_nodes()
            target_names = self._target_names()
            partners_names = ','.join([n.name() for n in self.partners])
            message = f'Found {partners_names} as partners for {target_names}'
            if self.partner_class:
                message = f'Found {partners_names} {self.partner_class}s as partners for {target_names}'
            print (message)
        
        return self.partners


    def set_partners(self):
        """
        Assumes partners are linked to Targets with Effect=positive
        """
        try:
        # case when self.params has 'partner_names' and 'partner_class' explicitly
        # use it when receptor partners are metabolites and therefore cannot be found by select_partners()
            partner_names = self.params['partner_names']
            try:
                self.partner_class = self.params['partner_class']
            except KeyError:
                print ('"partner_class" parameter is not specified !!!')
                self.partner_class = ''
        except KeyError:
            partner_names = ''

        if partner_names: 
                partners_names_s = OQL.join_with_quotes(partner_names)
                SELECTpartners = f'SELECT Entity WHERE Name = ({partners_names_s})'
                partners_graph = self.process_oql(SELECTpartners)
                self.partners = partners_graph._get_nodes()
                self.find_partners_oql = OQL.get_objects(ResnetGraph.dbids(self.partners))
        else:
            self.find_partners()

        
    def _get_report_name(self):
        indics = ','.join(self.params['indication_types'])
        rep_pred = 'suggested ' if self.params['strict_mode'] else 'suggested,predicted ' 
        return str(self.data_dir+rep_pred+ indics+' for '+ self._target_names())


    def GVindications(self):
        target_names = [x['Name'][0] for x in self.targets]
        t_n = ','.join(target_names)
        #select_targets = OQL.get_objects(self.target_ids)
  
        selectGVs = f'SELECT Entity WHERE objectType = GeneticVariant AND Connected by (SELECT Relation WHERE objectType = GeneticChange) to ({self.oql4targets})'
        #selectGVs = selectGVs.format(select_target=self.oql4targets)
        REQUEST_NAME = f'Find indications linked to {t_n} genetic variants'
        indication_type=','.join(self.params['indication_types'])
        OQLquery = f'SELECT Relation WHERE NeighborOf({selectGVs}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_type}))'
        self.GVsInDiseaseNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        indications = self.GVsInDiseaseNetwork.psobjs_with(only_with_values=self.params['indication_types'])
        self.GVs = self.GVsInDiseaseNetwork._psobjs_with('GeneticVariant','ObjTypeName')
        
        print('Found %d indications genetically linked to %d Genetic Variants in %s' % 
            (len(indications), len(self.GVs), t_n))
        
        return indications


    def __resolve_conflict_indications(self):

        def __resolve(conflict:PSObject, all_rels:list, using_rel_type:str):
            only_rels_with_type = [r for r in all_rels if r.objtype() == using_rel_type]
            if only_rels_with_type:
                only_rels_with_type.sort(key=lambda x: x.get_reference_count(), reverse=True)
                best_effect = only_rels_with_type[0]['Effect'][0]
                if best_effect == 'positive':
                    self.indications4agonists.remove(conflict)
                    return True
                elif best_effect == 'negative':
                    self.indications4antagonists.remove(conflict)
                    return True
                else:
                    return False
            return False

        conflict_indications = self.indications4antagonists.intersection(self.indications4agonists)
        unique_indications = len(self.indications4antagonists)+len(self.indications4agonists)-len(conflict_indications)
        print(f'Resolving {len(conflict_indications)} conflicting out of total {unique_indications} indications')
        unresolved_ids = set()
        target_uids = ResnetGraph.uids(self.targets)
        partner_uids = ResnetGraph.uids(self.partners)
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


    def __oql4indications_type(self):
        indication_types_str = ','.join(self.params['indication_types'])
        return f'SELECT Entity WHERE objectType = ({indication_types_str})'
    

    def find_target_indications(self):
        f_t = self.oql4targets
        t_n = self._target_names()
        # initializing self.indications4antagonists
        effect = 'positive'
        select_indications_by_type = self.__oql4indications_type()
        
        REQUEST_NAME = f'Find indications {effect}ly modulated by {t_n}'
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = {effect} AND \
            NeighborOf({target}) AND NeighborOf ({indications})' 
        OQLquery =  OQLquery.format(target=f_t,effect = effect,indications=select_indications_by_type)  
        ModulatedByTargetNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        self.indications4antagonists = set(ModulatedByTargetNetwork.psobjs_with(only_with_values=self.params['indication_types']))
        self.indications4antagonists_strict = set(self.indications4antagonists)
        print('Found %d diseases %sly regulated by target %s' % (len(self.indications4antagonists),effect,t_n))

        REQUEST_NAME = 'Find indications {effect}ly regulating {target}'.format(effect=effect, target=t_n)
        OQLquery = 'SELECT Relation WHERE objectType = QuantitativeChange AND Effect = {effect} AND \
            NeighborOf ({target}) AND NeighborOf ({indications})'
        OQLquery = OQLquery.format(target=f_t,effect = effect,indications=select_indications_by_type)
        ActivatedInDiseaseNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        add2indications = ActivatedInDiseaseNetwork.psobjs_with(only_with_values=self.params['indication_types'])
        self.indications4antagonists.update(add2indications)
        print('Found %d diseases where %s is %sly regulated' % (len(add2indications),t_n,effect))

        # initializing self.indications4agonists
        effect =  'negative'
        REQUEST_NAME = 'Find indications {effect}ly modulated by {target}'.format(effect=effect, target=t_n)
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = {effect} AND NeighborOf({select_target}) AND NeighborOf ({indications})'      
        OQLquery = OQLquery.format(select_target=f_t,effect=effect,indications=select_indications_by_type)
        ModulatedByTargetNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        self.indications4agonists = set(ModulatedByTargetNetwork.psobjs_with(only_with_values=self.params['indication_types']))
        self.indications4agonists_strict = set(self.indications4agonists)
        print('Found %d diseases %sly regulated by %s' % (len(self.indications4agonists),effect,t_n))

        REQUEST_NAME = 'Find indications {effect}ly regulating {target}'.format(effect=effect, target=t_n)
        OQLquery = 'SELECT Relation WHERE objectType = QuantitativeChange AND Effect = {effect} AND NeighborOf ({select_target}) AND NeighborOf ({indications})'  
        OQLquery = OQLquery.format(select_target=f_t,effect = effect,indications=select_indications_by_type)
        ActivatedInDiseaseNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        add2indications = ActivatedInDiseaseNetwork.psobjs_with(only_with_values=self.params['indication_types'])
        self.indications4agonists.update(add2indications)
        print('Found %d diseases where %s is %sly regulated' % (len(add2indications),t_n,effect))

        # initializing unknown_effect_indications
        gv_indications = self.GVindications()
        self.unknown_effect_indications_strict = set(gv_indications).difference(self.__known_effect_indications())
        print('%d indications linked to %s genetic variations were not found by previous searches' % (len(self.unknown_effect_indications_strict),t_n))

        REQUEST_NAME = 'Find indications where {target} is biomarker'.format(target=t_n)
        OQLquery = 'SELECT Relation WHERE objectType = Biomarker AND NeighborOf({select_target}) AND NeighborOf ({indications})'
        OQLquery = OQLquery.format(select_target=f_t,indications=select_indications_by_type)
        BiomarkerInDiseaseNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        biomarker_indications = BiomarkerInDiseaseNetwork.psobjs_with(only_with_values=self.params['indication_types'])
        print('Found %d indications where target %s is claimed as biomarker' %  (len(biomarker_indications),t_n))
        add2indications = set(biomarker_indications).difference(self.__known_effect_indications())
        print('%d indications having %s as biomarker were not found by previous searches' % (len(add2indications),t_n))
        self.unknown_effect_indications = set(add2indications)

        

    def modulators_effects(self,linked_by:list(),with_effect_on_target:str,find_indications=True,min_refcount=1,drugs_only=False):
        f_t = self.oql4targets
        t_n = self._target_names()
        reltype_str = ','.join(linked_by)

        REQUEST_NAME = f'Find substances {with_effect_on_target}ly regulating {t_n} by {reltype_str}'
        OQLquery = f'SELECT Relation WHERE objectType = ({reltype_str}) AND Effect = {with_effect_on_target} AND \
            NeighborOf upstream ({f_t})'
        if drugs_only:
            nb = f' AND NeighborOf downstream ({OQL.select_drugs()})'
        else:
            nb = ' AND NeighborOf downstream (SELECT Entity WHERE objectType = SmallMol)'

        OQLquery += nb

        if min_refcount > 1:
             OQLquery += ' AND '+REFCOUNT+' >= '+str(min_refcount)

        TargetInhibitorsNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        modulators = TargetInhibitorsNetwork._psobjs_with('SmallMol','ObjTypeName')
        print('Found %d substances %sly regulating %s by %s' % (len(modulators),with_effect_on_target,t_n,reltype_str))
        
        # now find indications for modulators found on previous step
        if modulators:
            if self.__is_strict():
                all_indications = list()
            else:
                indication_type=','.join(self.params['indication_types'])
                modulators_dbids = ResnetGraph.dbids(modulators)
                get_modulators=OQL.get_objects(modulators_dbids)
                effect_type = 'negative' if find_indications else 'positive'

                REQUEST_NAME = f'Find indications for substances {with_effect_on_target}ly regulating {t_n} by {reltype_str}'
                OQLquery = f'SELECT Relation WHERE objectType = Regulation AND Effect = {effect_type} AND \
                    NeighborOf (SELECT Entity WHERE objectType = ({indication_type})) AND NeighborOf ({get_modulators})'
                InhibitorsIndicationNetwork = self.process_oql(OQLquery,REQUEST_NAME)
                indication_by_regulation = InhibitorsIndicationNetwork.psobjs_with(only_with_values=self.params['indication_types'])
                print('Found %d indications for %d substances %sly regulating %s by %s relations' %  
                    (len(indication_by_regulation), len(modulators),with_effect_on_target,t_n,reltype_str))

                REQUEST_NAME = f'Find clinical trials for substances {with_effect_on_target}ly regulating {t_n} by {reltype_str}'
                OQLquery = f'SELECT Relation WHERE objectType = ClinicalTrial AND \
                    NeighborOf (SELECT Entity WHERE objectType = ({indication_type})) AND NeighborOf ({get_modulators})'
                InhibitorsIndicationNetwork.add_graph(self.process_oql(OQLquery,REQUEST_NAME))

                all_indications = InhibitorsIndicationNetwork.psobjs_with(only_with_values=self.params['indication_types'])
                print('Found %d indications on clinical trials with %d substances %sly regulating %s by %s relations' %  
                    (len(all_indications)- len(indication_by_regulation), len(modulators),with_effect_on_target,t_n,reltype_str))

            return modulators, all_indications
        else:
            return [],[]


    def modulators_indications(self,linked_by:list(),with_effect_on_target:str,min_refcount=1,drugs_only=False):
        return self.modulators_effects(linked_by,with_effect_on_target,True,min_refcount,drugs_only)


    def modulators_toxicities(self,linked_by:list(),with_effect_on_target:str,min_refcount=1,drugs_only=False):
        return self.modulators_effects(linked_by,with_effect_on_target,False,min_refcount,drugs_only)


    def indications4chem_modulators(self):
        exist_indication_count = len(self.__indications())
        self.DirectAntagonists, indications = self.modulators_indications(['DirectRegulation'],'negative')
        self.indications4antagonists.update(indications)
        new_indication_count = len(self.__indications()) - exist_indication_count
        print('%d indications for drugs directy inhibiting %s were not found by previous searches' %  
                    (new_indication_count,self._target_names()))
        self.IndirectAgonists, indications = self.modulators_indications(
        ['Regulation','Expression','MolTransport'],'positive',min_refcount=1,drugs_only=True)

        self.DirectAgonists, indications = self.modulators_indications(['DirectRegulation'],'positive')
        self.indications4agonists.update(indications)
        self.IndirectAntagonists, indications = self.modulators_indications(
            ['Regulation','Expression','MolTransport'],'negative',min_refcount=1,drugs_only=True)

 
    def indications4partners(self):
        """
        Assumes partners are linked to Target with Effect=positive
        """
        if not self.partners: return

        t_n = self._target_names()
        indications = self.__indications()
        indications_dbids = ResnetGraph.dbids(indications)
        exist_indication_count = len(indications)
        partners_s = self.partner_class.lower() if self.partner_class else 'partner'
        oql4indications_type = self.__oql4indications_type()

        REQUEST_NAME = f'Find indications positively regulated by {partners_s}s of {t_n}'
        OQLtemplate = 'SELECT Relation WHERE objectType = (Regulation,QuantitativeChange) AND \
                Effect = {eff} AND NeighborOf ({partners}) AND NeighborOf ({indications})'

        OQLquery = OQLtemplate.format(eff='positive',partners=self.find_partners_oql,indications=oql4indications_type)
        PartnerIndicationNetwork4anatagonists = self.iterate_oql(OQLquery,indications_dbids,request_name=REQUEST_NAME)
        indications = PartnerIndicationNetwork4anatagonists.psobjs_with(only_with_values=self.params['indication_types'])
        self.indications4antagonists.update(indications)
        self.indications4antagonists_strict.update(indications)
        print('Found %d indications for %d %s %ss' %  
                (len(indications), len(self.partners),t_n,self.partner_class.lower()))

        REQUEST_NAME = f'Find indications negatively regulated by {self.partner_class.lower()}s of {t_n}'
        OQLquery = OQLtemplate.format(eff='negative',partners=self.find_partners_oql,indications=oql4indications_type)
        PartnerIndicationNetwork4agonists = self.iterate_oql(OQLquery,indications_dbids,request_name=REQUEST_NAME)
        indications = PartnerIndicationNetwork4agonists.psobjs_with(only_with_values=self.params['indication_types'])
        self.indications4agonists.update(indications)
        self.indications4agonists_strict.update(indications)

        new_indication_count = len(self.__indications()) - exist_indication_count
        print('%d indications for %d %s %ss were not found by previous searches' %  
                (new_indication_count, len(self.partners),t_n,self.partner_class.lower()))
        

    def indications4cells_secreted_target(self):
        if not self.target_class == 'Ligand': return

        t_n = self._target_names()
        exist_indication_count = len(self.__indications())

        REQUEST_NAME = f'Find cells secreting {t_n}'
        OQLquery = 'SELECT Relation WHERE objectType = (CellExpression,MolTransport) AND NeighborOf ({select_targets}) AND NeighborOf (SELECT Entity WHERE objectType = Cell)'
        cells_make_target = self.process_oql(OQLquery.format(select_targets=self.oql4targets),REQUEST_NAME)
        self.ProducingCells = cells_make_target._psobjs_with('CellType','ObjTypeName')
        print('Found %d cell types producing %s' % (len(self.ProducingCells),t_n))

        if self.__is_strict(): # no further predictions of indications
            return # found Cells will be used for target ranking
            
        ProducingCells_dbids = ResnetGraph.dbids(self.ProducingCells)
        REQUEST_NAME = f'Find indications positively linked to cell secreting {t_n}'
        OQLtemplate = 'SELECT Relation WHERE Effect = {effect} AND NeighborOf (SELECT Entity WHERE objectType = ({indication_type})) AND NeighborOf (SELECT Entity WHERE id = ({cell_ids}))'
        OQLquery = OQLtemplate.format(effect='positive',cell_ids=OQL.id2str(ProducingCells_dbids),indication_type=','.join(self.params['indication_types']))
        CellDiseaseNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        indications = CellDiseaseNetwork.psobjs_with(only_with_values=self.params['indication_types'])
        self.indications4antagonists.update(indications)

        REQUEST_NAME = f'Find indications negatively linked to cell secreting {t_n}'
        OQLquery = OQLtemplate.format(effect='negative',cell_ids=OQL.id2str(ProducingCells_dbids),indication_type=','.join(self.params['indication_types']))
        CellDiseaseNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        indications = CellDiseaseNetwork.psobjs_with(only_with_values=self.params['indication_types'])
        self.indications4antagonists.update(indications)
            
        new_indication_count = len(self.__indications()) - exist_indication_count
        print('%d indications for %d cells producing %s were not found by previous searches' %  
                (new_indication_count, len(self.ProducingCells),t_n))


    def pathway_oql(self):
        #set pathway_name_must_include_target to True if targets have a lot of large curated pathways
        target_oqls = dict()
        pct = '%'
        merged_pathways = 'SELECT Relation WHERE objectType = (DirectRegulation,Binding,ProtModification,PromoterBinding,ChemicalReaction) AND MemberOf ({select_networks})'

        for target in self.targets:
            target_id = target['Id'][0]
            target_name = target['Name'][0]
            SELECTpathways = 'SELECT Network WHERE ParentOf (SELECT Entity WHERE id = ({i}))'.format(i=str(target_id))
            if self.params['pathway_name_must_include_target']:
                SELECTpathways = SELECTpathways + ' AND Name LIKE (\''+pct+target_name+pct+'\')' #additional refinement for popular targets

            target_oqls[target_name] = SELECTpathways

        if self.partners:
            target_partner_oqls = dict()
            for partner in self.partners:
                for target_name, oql_query in target_oqls.items():
                    find_pathways_query = oql_query + ' AND ParentOf (SELECT Entity WHERE id = ({i}))'.format(i=str(partner.dbid()))
                    target_partner_oqls[(target_name,partner.name())] = merged_pathways.format(select_networks=find_pathways_query)

            return target_partner_oqls
        else:
            tuple_target_oqls = dict()
            for target_name, oql_query in target_oqls.items():
                tuple_target_oqls[(target_name, '')] = merged_pathways.format(select_networks=oql_query)

            return tuple_target_oqls
    

    def get_pathway_componets(self):
        #finding downstream pathway components
        oql_queries = self.pathway_oql() # separate oql_query for each target
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
            targets_regulome = merged_pathway.get_regulome(self.targets)
            self.PathwayComponents = set(targets_regulome._get_nodes())
            print (f'Found regulome with {len(self.PathwayComponents)} components')
            return targets_regulome
        else: 
            print('No curated pathways were found for %s' % self._target_names())
            return ResnetGraph()


    def init_semantic_search(self):
        '''
        Loads DF4ANTAGONISTS and DF4AGONISTS df to raw_data
        '''
        print('\n\nInitializing semantic search')
        t_n = self._target_names()

        indications4antagonists = self._indications4antagonists()
        if indications4antagonists:
            indication_df = self.load_df(list(indications4antagonists),max_children_count=11)
            indication_df._name_ = DF4ANTAGONISTS
            self.add2raw(indication_df)
            print('Will score %d indications for antagonists of %s' % (len(indication_df),t_n))

        indications4agonists = self._indications4agonists()
        if indications4agonists:
            indication_df = self.load_df(list(indications4agonists),max_children_count=11)
            indication_df._name_ = DF4AGONISTS
            self.add2raw(indication_df)
            print('Will score %d indications for agonists of %s' % (len(indication_df),t_n))
        
        if indications4antagonists or indications4agonists:
            return True
        else:
            print ('No indications found for %s' % t_n)
            if self.params['strict_mode']:
                print('Try setting strict_mode to False')
            return False
        

    def score_GVs(self, df2score:df):
        if not self.GVids: return
        t_n = self._target_names()
        self.__colnameGV__ = t_n+' GVs'
        gvlinkcounter = 0
        if hasattr(self,'GVsInDiseaseNetwork'):
            for i in df2score.index:
                target_dbids = list(df2score.at[i,self.__temp_id_col__])
                row_targets = self.Graph.psobj_with_dbids(target_dbids)
                targetGVs = self.GVsInDiseaseNetwork.get_neighbors(row_targets,self.GVs)
                    
                GVscore = 0
                if len(targetGVs) > 0:
                    GVnames = [g.name() for g in targetGVs]
                    df2score.at[i,self.__colnameGV__] = ';'.join(GVnames)
                    gvlinkcounter += 1
                    gv_disease_subgraph = self.GVsInDiseaseNetwork.get_subgraph(row_targets,targetGVs)
                    GVscore = len(gv_disease_subgraph.load_references())
                
                df2score.at[i,self._col_name_prefix+'GVs'] = GVscore
            print('Found %d indications linked to %d GVs' % (gvlinkcounter, len(self.GVids)))


    def __drug_connect_params(self,direct_modulators:bool,score_antagonists:bool):
        # function returns how2connect paramters to score molecules similar to desired drug affect on the targets 
        effect = None
        drug_class = None

        if score_antagonists:
        # most common case when targets must be inhibited
            if direct_modulators:
                    effect = 'negative'
                    drug_class = 'direct antagonists'             
                    concepts = self.DirectAntagonists if hasattr(self,'DirectAntagonists') else []
            else:
                    effect = 'negative'
                    drug_class = 'indirect antagonists'             
                    concepts = self.IndirectAntagonists if hasattr(self,'IndirectAntagonists') else []
        else:
        # case if drug are agonists
            if direct_modulators:
                    effect = 'positive'
                    drug_class = 'direct agonists'
                    concepts = self.DirectAgonists if hasattr(self,'DirectAgonists') else []
            else:
                    effect = 'positive'
                    drug_class = 'indirect agonists'
                    concepts = self.IndirectAgonists if hasattr(self,'IndirectAgonists') else []

        return effect, drug_class, concepts


    def __drug_tox_params(self,direct_modulators:bool,score_antagonists:bool):
        effect = None
        drug_class = None
        concepts = list()
        # function returns how2connect parameters to score molecules synergizing with target action on indication
        # i.e. molecules that have effect opposite to desired drug affect on the targets
        if score_antagonists:
        # most common case when targets must be inhibited
            if direct_modulators:
                    effect = 'positive'
                    drug_class = 'direct agonists'
                    concepts = self.DirectAgonists if hasattr(self,'DirectAgonists') else []
            else:
                    effect = 'positive'
                    drug_class = 'indirect agonists'
                    concepts = self.IndirectAgonists if hasattr(self,'IndirectAgonists') else []
        else:
        # case if drug are agonists
            if direct_modulators:
                    effect = 'negative'
                    drug_class = 'direct antagonists'             
                    concepts = self.DirectAntagonists if hasattr(self,'DirectAntagonists') else []
            else:
                    effect = 'negative'
                    drug_class = 'indirect antagonists'             
                    concepts = self.IndirectAntagonists if hasattr(self,'IndirectAntagonists') else []

        return effect, drug_class, concepts


    def semantic_score(self,in_worksheet:str,target_effect_on_indication:str):
        """
        Input
        -----
        target_effect_on_indication: required Effect sign between target and Indication 
        target_effect_on_indication = 'positive' to score antagonists
        target_effect_on_indication = 'negative' to score agonists

        Output
        ------
        adds score columns to "in_worksheet" from self.raw_data
        """
        t_n = self._target_names()
        target_in_header = t_n if len(t_n) < 45 else 'targets'
        colname = 'Activated by ' if target_effect_on_indication == 'positive' else 'Inhibited by '
        colname += target_in_header

        indication_df = self.raw_data[in_worksheet]
        score4antagonists = True if target_effect_on_indication == 'positive' else False

        if self.__is_strict():
            booster_reltypes = ['Regulation','Biomarker','GeneticChange','QuantitativeChange','StateChange']
        else:
            booster_reltypes = ['Regulation','GeneticChange']
        self.set_how2connect(['Regulation'],[target_effect_on_indication],'',booster_reltypes)
        linked_row_count,linked_ent_ids,indication_df = self.link2concept(colname,self.targets,indication_df)
        print('%d indications are %sly regulated by %s' % (linked_row_count,target_effect_on_indication,t_n))

        self.score_GVs(indication_df)

        link_effect, drug_class, concepts = self.__drug_connect_params(True,score4antagonists)
        if concepts:
            # references suggesting that known drugs for the target as treatments for indication
            colname = target_in_header+' '+drug_class+' clin. trials'
            self.set_how2connect(['ClinicalTrial'],[],'')
            linked_row_count,linked_ent_ids,indication_df = self.link2concept(colname,concepts,indication_df)
            print('Linked %d clinical trial indictions for %s %s' % (linked_row_count,t_n,drug_class))

            colname = target_in_header+' '+drug_class
            self.set_how2connect(['Regulation'],[link_effect],'',['Regulation'])
            linked_row_count,linked_ent_ids,indication_df = self.link2concept(colname,concepts,indication_df)
            print('Linked %d indications for %s %s' % (linked_row_count,t_n,drug_class))

        #references reporting target agonists exacerbating indication or causing indication as adverse events
        link_effect, drug_class, concepts = self.__drug_tox_params(True,score4antagonists)
        if concepts:
            colname = target_in_header+' '+drug_class
            self.set_how2connect(['Regulation'],[link_effect],'',['Regulation'])
            linked_row_count,linked_ent_ids,indication_df = self.link2concept(colname,concepts,indication_df)
            print('Linked %d indications as toxicities for %s %s' % (linked_row_count,t_n,drug_class))

    
        #references where target expression or activity changes in the indication
        colname = ' is upregulated' if target_effect_on_indication == 'positive' else ' is downregulated'
        colname = target_in_header + colname
        self.set_how2connect(['QuantitativeChange'],[target_effect_on_indication],'',['Biomarker','StateChange'])
        linked_row_count,linked_ent_ids,indication_df= self.link2concept(colname,self.targets,indication_df)
        print('%d indications %sly regulate %s' % (linked_row_count,target_effect_on_indication,t_n))

        #references suggesting target partners as targets for indication
        if self.partners:
            p_cl = self.partner_class if self.partner_class else 'partner'
            colname = f'{target_in_header} {p_cl}s'
            self.set_how2connect(['Regulation'],[target_effect_on_indication],'',['Regulation'])
            linked_row_count,linked_ent_ids,indication_df = self.link2concept(colname,self.partners,indication_df)
            print('Linked %d indications for %d %s %ss' % (linked_row_count,len(self.partners),t_n,p_cl))

        # references reporting that cells producing the target linked to indication  
        # only used if taregts are secretred ligands
        if hasattr(self, 'ProducingCells'):
            colname = f'{target_in_header} producing cells'
            self.set_how2connect(['Regulation'],[target_effect_on_indication],'',['Regulation'])
            linked_row_count,linked_ent_ids,indication_df = self.link2concept(colname,self.ProducingCells,indication_df)
            print('Liked %d indications linked %d cells producing %s' % (linked_row_count,len(self.ProducingCells),t_n))

        link_effect, drug_class, concepts = self.__drug_connect_params(False,score4antagonists)
        if concepts:
            # references suggesting that known drugs for the target as treatments for indication
            colname = target_in_header+' '+drug_class+' clin. trials'
            self.set_how2connect(['ClinicalTrial'],[],'')
            new_session = self._clone(to_retrieve=REFERENCE_IDENTIFIERS)
            linked_row_count,linked_ent_ids,indication_df = new_session.link2concept(colname,concepts,indication_df)
            print('Linked %d clinical trial indications for %s %s' % (linked_row_count,t_n,drug_class))
            new_session.close_connection()

            colname = target_in_header+' '+drug_class
            self.set_how2connect(['Regulation'],[link_effect],'',['Regulation'])
            new_session = self._clone(to_retrieve=REFERENCE_IDENTIFIERS)
            linked_row_count,linked_ent_ids,indication_df = new_session.link2concept(colname,concepts,indication_df)
            print('Linked %d indications for %s %s' % (linked_row_count,t_n,drug_class))
            new_session.close_connection()

        #references reporting target agonists exacerbating indication or causing indication as adverse events
        link_effect, drug_class, concepts = self.__drug_tox_params(False,score4antagonists)
        if concepts:
            colname = target_in_header+' '+drug_class
            self.set_how2connect(['Regulation'],[link_effect],'',['Regulation'])
            new_session = self._clone(to_retrieve=REFERENCE_IDENTIFIERS)
            linked_row_count,linked_ent_ids,indication_df = new_session.link2concept(colname,concepts,indication_df)
            print('Linked %d indications as toxicities for %s %s' % (linked_row_count,t_n,drug_class))
            new_session.close_connection()
        
        if hasattr(self, 'PathwayComponents'):
            #references linking target pathway to indication
            colname = target_in_header + ' pathway components'
            self.set_how2connect(['Regulation'],[target_effect_on_indication],'',step=125)
            new_session = self._clone(to_retrieve=REFERENCE_IDENTIFIERS)
            linked_row_count,linked_ent_ids,indication_df = new_session.link2concept(colname,list(self.PathwayComponents),indication_df)
            print('Linked %d indications to %s pathway components' % (linked_row_count,t_n))
            new_session.close_connection()

        counts_df = self.make_count_df(indication_df,with_name=in_worksheet)
        self.add2raw(counts_df)
        return


    def other_effects(self):
        # need to be called after ranking to subtract self.all_entity_ids
        print('Findind indication linked with unknown effect to %s' % self._target_names())
        old_rel_props = self.relProps
        self.add_rel_props(PS_SENTENCE_PROPS+list(PS_BIBLIO_PROPS))
        t_n = self._target_names()
        ranked_indication_ids = self.Graph.dbids4nodes(self.params['indication_types'])
        REQUEST_NAME = 'Find indications modulated by {target} with unknown effect'.format(target=t_n)
        oql4indications_type = self.__oql4indications_type()
        OQLquery = f'SELECT Relation WHERE objectType = (Regulation,QuantitativeChange) AND Effect = unknown AND \
            NeighborOf({self.oql4targets}) AND NeighborOf ({oql4indications_type})'
        to_return = self.process_oql(OQLquery,REQUEST_NAME)
        
        REQUEST_NAME = 'Find indications where {target} was suggested as Biomarker'.format(target=t_n)
        OQLquery = f'SELECT Relation WHERE objectType = (Biomarker,StateChange,GeneticChange) AND \
            NeighborOf({self.oql4targets}) AND NeighborOf ({oql4indications_type})'
        to_return.add_graph(self.process_oql(OQLquery,REQUEST_NAME))
        
        REQUEST_NAME = f'Find indications with genetically linked {t_n} Genetic Variants'
        OQLquery = f'SELECT Relation WHERE objectType = FunctionalAssociation AND NeighborOf ({oql4indications_type})'
        OQLquery += ' AND SELECT Entity WHERE id = ({ids})'
        GVdbids = ResnetGraph.dbids(self.GVs)
        to_return.add_graph(self.iterate_oql(OQLquery,GVdbids,request_name=REQUEST_NAME))

        to_return.remove_nodes_from(ranked_indication_ids) # now delete all indication with known effect
        indications = to_return.psobjs_with(only_with_values=self.params['indication_types'])

        print('Found %d new indications linked to %s with unknown effect' % (len(indications),self._target_names()))
        self.relProps = old_rel_props
        return to_return


    def load_indications4targets(self):
        '''
        Input
        -----
        target names must be in self.params['target_names']
        '''
        start_time = time.time()
        self.set_targets()
        self.set_partners()

        self.find_target_indications()
        #return  #uncomment to enable fast downstream testing
        self.indications4partners()

        if self.target_class == 'Ligand':
            self.indications4cells_secreted_target()
        else:
            self.indications4chem_modulators()
            
        self.get_pathway_componets()
        self.__resolve_conflict_indications()
        print("%d indications for %s were found in %s" % 
              (len(self.__indications()), self._target_names(), self.execution_time(start_time)))


    def make_report(self):
        start_time = time.time()
        self.flush_dump_files()
        self.load_indications4targets()
        
        if self.init_semantic_search():
            self.semantic_score(DF4ANTAGONISTS,target_effect_on_indication='positive')
            self.semantic_score(DF4AGONISTS,target_effect_on_indication='negative')

        self.normalize(DF4ANTAGONISTS,ANTAGONISTS_WS)
        self.normalize(DF4AGONISTS,AGONISTS_WS)
        other_effects_graph = self.other_effects()
        other_indications = other_effects_graph.snippets2df(df_name='possibilities')
        self.add2report(other_indications)

        self.add_ps_bibliography()

        self.add_etm_refs(ANTAGONISTS_WS,self.input_names())
        self.add_etm_refs(AGONISTS_WS,self.input_names())
        self.add_etm_bibliography()
        print('Repurposing of %s was done in %s' % 
              (self._get_report_name(), self.execution_time(start_time)))

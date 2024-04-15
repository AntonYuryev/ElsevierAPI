from .SemanticSearch import SemanticSearch,len,df,pd
from .ResnetAPISession import REFERENCE_IDENTIFIERS,BIBLIO_PROPERTIES
from .ResnetGraph import REFCOUNT, PSObject, ResnetGraph
from ..ETM_API.references import PS_SENTENCE_PROPS,EFFECT,PS_BIBLIO_PROPS
from .PathwayStudioGOQL import OQL
from concurrent.futures import ThreadPoolExecutor,as_completed
from ElsevierAPI import  execution_time
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
    max_ontology_parent = 11
    max_threads4ontology = 4
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
                 'mode_of_action':ANY_MOA
                }
        my_kwargs.update(kwargs)

        super().__init__(*args,**my_kwargs)
        self.add_rel_props([EFFECT])
        self.columns2drop += [self.__resnet_name__,self.__mapped_by__]
        # 4 threads perform a bit faster than 8 threads 
        # 288 disease parents out of 288 were processed in 0:19:51.732684 by 8 threads
        # 288 disease parents out of 288 were processed in 0:18:52.715937 by 4 threads

        self.indications4antagonists = set()
        self.indications4antagonists_strict = set()
        self.indications4agonists = set()
        self.indications4agonists_strict = set()
        self.unknown_effect_indications = set()
        self.unknown_effect_indications_strict = set()
        self.partners = list() # list of PSObject
        self.GVs = list()
        self.target_class = ''


    def _input_targets(self):
        """
        Returns target names as they were entered in configuration.
        return values may be different from self.target_names()
        """
        return ','.join(self.params['target_names']) if len(self.params['target_names']) < 3 else 'targets'
    
    
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


    def input_names(self):
        """
        Returns
        -------
        entity names that can be used both in PS and ETM
        """
        #input names are used for finding references in ETM.
        # RepurposeDrug has it own 'target_names'.
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


    def __indications4antagonists(self)->set[PSObject]:
        return self.indications4antagonists_strict if self.params['strict_mode'] else self.indications4antagonists


    def __indications4agonists(self)->set[PSObject]:    
        return self.indications4agonists_strict if self.params['strict_mode'] else self.indications4agonists


    def __all_indications(self):
        return self.__indications4antagonists()|self.__indications4agonists()
    

    def _is_strict(self):
        return self.params['strict_mode']
  

    def set_targets(self):
        try:
            target_objtype_str = str(','.join(self.params['target_objtypes']))
        except KeyError: 
            target_objtype_str = 'Protein'

        if self.params['target_names']:
            prop_names_str, prop_values_str = OQL.get_search_strings(['Name'],self.params['target_names'])
            self.add_ent_props(['Class'])
            self.oql4targets = f'SELECT Entity WHERE ({prop_names_str}) = ({prop_values_str}) AND objectType = ({target_objtype_str})'
            targets_graph = self.process_oql(self.oql4targets)
            if isinstance(targets_graph,ResnetGraph):
                self.targets = targets_graph._get_nodes()
                target_dbids = [x['Id'][0] for x in self.targets]
                self.oql4targets = OQL.get_objects(target_dbids)
                for t in self.targets:
                    try:
                        self.target_class = t['Class'][0]
                        break
                    # assumes all targets have the same class
                    except KeyError:
                        continue


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
            if isinstance(partners_graph,ResnetGraph):
                if self.partner_class == 'Ligand':
                    # additional request to find secreted molecules that are not annotated with class Ligand
                    SELECTsecretedpartners = 'SELECT Entity WHERE "Cell Localization" = Secreted AND objectType = Protein AND Connected by (SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive) to ({select_target})'
                    secreted_partners = self.process_oql(SELECTsecretedpartners.format(select_target=self.oql4targets))
                    if isinstance(secreted_partners,ResnetGraph):
                        partners_graph = secreted_partners.compose(partners_graph)

                self.partners = partners_graph._get_nodes()
                target_names = self._input_targets()
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
            if isinstance(partners_graph,ResnetGraph):
                self.partners = partners_graph._get_nodes()
        else:
            self.find_partners()

        self.find_partners_oql = OQL.get_objects(ResnetGraph.dbids(self.partners))
        return

        
    def report_path(self, extension='.xlsx'):
        indics = ','.join(self.params['indication_types'])
        rep_pred = 'suggested ' if self.params['strict_mode'] else 'suggested,predicted '
        if self.params['mode_of_action'] == ANTAGONIST:
            mode = ' inhibition'
        elif self.params['mode_of_action'] == AGONIST:
            mode = ' activation'
        else:
            mode = ''
        
        fname = rep_pred+ indics+' for '+ self._input_targets()+mode+extension
        return os.path.join(self.data_dir, fname)


    def GVindications(self)->list[PSObject]:
        target_names = [x['Name'][0] for x in self.targets]
        t_n = ','.join(target_names)
        selectGVs = f'SELECT Entity WHERE objectType = GeneticVariant AND Connected by (SELECT Relation WHERE objectType = GeneticChange) to ({self.oql4targets})'

        REQUEST_NAME = f'Find indications linked to {t_n} genetic variants'
        indication_type=','.join(self.params['indication_types'])
        OQLquery = f'SELECT Relation WHERE NeighborOf({selectGVs}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_type}))'
        self.GVs2DiseaseGraph = self.process_oql(OQLquery,REQUEST_NAME)
        if isinstance(self.GVs2DiseaseGraph,ResnetGraph):
            indications = self.GVs2DiseaseGraph.psobjs_with(only_with_values=self.params['indication_types'])
            self.GVs = self.GVs2DiseaseGraph._psobjs_with('GeneticVariant','ObjTypeName')
            
            print('Found %d indications genetically linked to %d Genetic Variants in %s' % 
                (len(indications), len(self.GVs), t_n))
        
            return indications
        else:
            return list()


    def __resolve_conflict_indications(self):

        def __resolve(conflict:PSObject, all_rels:list, using_rel_type:str):
            only_rels_with_type = [r for r in all_rels if r.objtype() == using_rel_type]
            if only_rels_with_type:
                only_rels_with_type.sort(key=lambda x: x.count_refs(), reverse=True)
                best_effect = only_rels_with_type[0]['Effect'][0]
                if best_effect == 'positive':
                    self.indications4agonists.discard(conflict)
                    self.indications4agonists_strict.discard(conflict)
                    return True
                elif best_effect == 'negative':
                    self.indications4antagonists.discard(conflict)
                    self.indications4antagonists_strict.discard(conflict)
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
    

    def find_target_indications(self,moa=ANY_MOA):
        if not hasattr(self,'oql4targets'): 
            return
        f_t = self.oql4targets
        t_n = self._input_targets()
        select_indications_by_type = self.__oql4indications_type()
        
        if moa in [ANTAGONIST,ANY_MOA]:
            effect = 'positive'
            REQUEST_NAME = f'Find indications {effect}ly modulated by {t_n}'
            OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = {effect} AND \
                NeighborOf({target}) AND NeighborOf ({indications})' 
            OQLquery =  OQLquery.format(target=f_t,effect = effect,indications=select_indications_by_type)  
            ModulatedByTargetNetwork = self.process_oql(OQLquery,REQUEST_NAME)
            if isinstance(ModulatedByTargetNetwork,ResnetGraph):
                self.indications4antagonists = set(ModulatedByTargetNetwork.psobjs_with(only_with_values=self.params['indication_types']))
                self.indications4antagonists_strict = set(self.indications4antagonists)
                print('Found %d diseases %sly regulated by target %s' % (len(self.indications4antagonists),effect,t_n))

            REQUEST_NAME = 'Find indications {effect}ly regulating {target}'.format(effect=effect, target=t_n)
            OQLquery = 'SELECT Relation WHERE objectType = QuantitativeChange AND Effect = {effect} AND \
                NeighborOf ({target}) AND NeighborOf ({indications})'
            OQLquery = OQLquery.format(target=f_t,effect = effect,indications=select_indications_by_type)
            ActivatedInDiseaseNetwork = self.process_oql(OQLquery,REQUEST_NAME)
            if isinstance(ActivatedInDiseaseNetwork,ResnetGraph):
                add2indications = ActivatedInDiseaseNetwork.psobjs_with(only_with_values=self.params['indication_types'])
                self.indications4antagonists.update(add2indications)
                print('Found %d diseases where %s is %sly regulated' % (len(add2indications),t_n,effect))
        
        if moa in [AGONIST,ANY_MOA]:
            effect =  'negative'
            REQUEST_NAME = f'Find indications {effect}ly modulated by {t_n}'
            OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = {effect} AND NeighborOf({select_target}) AND NeighborOf ({indications})'      
            OQLquery = OQLquery.format(select_target=f_t,effect=effect,indications=select_indications_by_type)
            ModulatedByTargetNetwork = self.process_oql(OQLquery,REQUEST_NAME)
            if isinstance(ModulatedByTargetNetwork,ResnetGraph):
                self.indications4agonists = set(ModulatedByTargetNetwork.psobjs_with(only_with_values=self.params['indication_types']))
                self.indications4agonists_strict = set(self.indications4agonists)
                print('Found %d diseases %sly regulated by %s' % (len(self.indications4agonists),effect,t_n))

            REQUEST_NAME = 'Find indications {effect}ly regulating {target}'.format(effect=effect, target=t_n)
            OQLquery = 'SELECT Relation WHERE objectType = QuantitativeChange AND Effect = {effect} AND NeighborOf ({select_target}) AND NeighborOf ({indications})'  
            OQLquery = OQLquery.format(select_target=f_t,effect = effect,indications=select_indications_by_type)
            ActivatedInDiseaseNetwork = self.process_oql(OQLquery,REQUEST_NAME)
            if isinstance(ActivatedInDiseaseNetwork,ResnetGraph):
                add2indications = ActivatedInDiseaseNetwork.psobjs_with(only_with_values=self.params['indication_types'])
                self.indications4agonists.update(add2indications)
                print('Found %d diseases where %s is %sly regulated' % (len(add2indications),t_n,effect))

        # initializing unknown_effect_indications
        gv_indications = self.GVindications()
        self.unknown_effect_indications_strict = set(gv_indications).difference(self.__known_effect_indications())
        print('%d indications linked to %s genetic variations were not found by previous searches' % (len(self.unknown_effect_indications_strict),t_n))

        REQUEST_NAME = f'Find indications where {t_n} is biomarker'
        OQLquery = 'SELECT Relation WHERE objectType = Biomarker AND NeighborOf({select_target}) AND NeighborOf ({indications})'
        OQLquery = OQLquery.format(select_target=f_t,indications=select_indications_by_type)
        BiomarkerInDiseaseNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        if isinstance(BiomarkerInDiseaseNetwork,ResnetGraph):
            biomarker_indications = BiomarkerInDiseaseNetwork.psobjs_with(only_with_values=self.params['indication_types'])
            print('Found %d indications where target %s is claimed as biomarker' %  (len(biomarker_indications),t_n))
            add2indications = set(biomarker_indications).difference(self.__known_effect_indications())
            print('%d indications having %s as biomarker were not found by previous searches' % (len(add2indications),t_n))
            self.unknown_effect_indications = set(add2indications)


    def indications4partners(self,moa=ANY_MOA):
        """
        Assumes partners are linked to Target with Effect=positive
        """
        if not self.partners: return

        t_n = self._input_targets()
        indications = self.__indications()
        indications_dbids = ResnetGraph.dbids(list(indications))
        exist_indication_count = len(indications)
        partners_s = self.partner_class.lower() if self.partner_class else 'partner'
        oql4indications_type = self.__oql4indications_type()

        OQLtemplate = 'SELECT Relation WHERE objectType = (Regulation,QuantitativeChange) AND \
Effect = {eff} AND NeighborOf ({partners}) AND NeighborOf ({indications})'

        if moa in [ANTAGONIST,ANY_MOA]:
            REQUEST_NAME = f'Find indications positively regulated by {partners_s}s of {t_n}'
            OQLquery = OQLtemplate.format(eff='positive',partners=self.find_partners_oql,indications=oql4indications_type)
            PartnerIndicationNetwork4anatagonists = self.iterate_oql(OQLquery,set(indications_dbids),request_name=REQUEST_NAME)
            indications = PartnerIndicationNetwork4anatagonists.psobjs_with(only_with_values=self.params['indication_types'])
            self.indications4antagonists.update(indications)
            self.indications4antagonists_strict.update(indications)
            print('Found %d indications for %d %s %ss' %  
                    (len(indications), len(self.partners),t_n,partners_s))

        if moa in [AGONIST,ANY_MOA]:
            REQUEST_NAME = f'Find indications negatively regulated by {partners_s}s of {t_n}'
            OQLquery = OQLtemplate.format(eff='negative',partners=self.find_partners_oql,indications=oql4indications_type)
            PartnerIndicationNetwork4agonists = self.iterate_oql(OQLquery,set(indications_dbids),request_name=REQUEST_NAME)
            indications = PartnerIndicationNetwork4agonists.psobjs_with(only_with_values=self.params['indication_types'])
            self.indications4agonists.update(indications)
            self.indications4agonists_strict.update(indications)
            print('Found %d indications for %d %s %ss' %  
                    (len(indications), len(self.partners),t_n,partners_s))

        new_indication_count = len(self.__indications()) - exist_indication_count
        print('%d indications for %d %s %ss were not found by previous searches' %  
                (new_indication_count, len(self.partners),t_n,partners_s))
    

    def modulators_effects(self,linked_by:list,with_effect_on_target:str,find_indications=True,min_refcount=1,drugs_only=False):
        f_t = self.oql4targets
        t_n = self._input_targets()
        reltype_str = ','.join(linked_by)

        REQUEST_NAME = f'Find substances {with_effect_on_target}ly regulating {t_n} by {reltype_str}'
        OQLquery = f'SELECT Relation WHERE objectType = ({reltype_str}) AND Effect = {with_effect_on_target} AND \
            NeighborOf upstream ({f_t})'
        if drugs_only:
            nb = f' AND NeighborOf downstream ({OQL.select_drugs()} AND Connectivity > 1)'
        else:
            nb = ' AND NeighborOf downstream (SELECT Entity WHERE objectType = SmallMol AND Connectivity > 1)'

        OQLquery += nb

        if min_refcount > 1:
             OQLquery += ' AND '+REFCOUNT+' >= '+str(min_refcount)

        modulators = list()
        TargetInhibitorsNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        if isinstance(TargetInhibitorsNetwork,ResnetGraph):
            modulators = TargetInhibitorsNetwork._psobjs_with('SmallMol','ObjTypeName')
            print('Found %d substances %sly regulating %s by %s' % (len(modulators),with_effect_on_target,t_n,reltype_str))
        if not modulators: return [],[]

        # now find indications for modulators found on previous step
        all_indications = list()
        if not self._is_strict():
            indication_type=','.join(self.params['indication_types'])
            modulators_dbids = ResnetGraph.dbids(modulators)
            get_modulators=OQL.get_objects(modulators_dbids)
            effect_type = 'negative' if find_indications else 'positive'

            REQUEST_NAME = f'Find indications for substances {with_effect_on_target}ly regulating {t_n} by {reltype_str}'
            OQLquery = f'SELECT Relation WHERE objectType = Regulation AND Effect = {effect_type} AND \
                NeighborOf (SELECT Entity WHERE objectType = ({indication_type})) AND NeighborOf ({get_modulators})'
            InhibitorsIndicationNetwork = self.process_oql(OQLquery,REQUEST_NAME)
            if isinstance(InhibitorsIndicationNetwork,ResnetGraph):
                indication_by_regulation = InhibitorsIndicationNetwork.psobjs_with(only_with_values=self.params['indication_types'])
                print('Found %d indications for %d substances %sly regulating %s by %s relations' %  
                    (len(indication_by_regulation), len(modulators),with_effect_on_target,t_n,reltype_str))

                REQUEST_NAME = f'Find clinical trials for substances {with_effect_on_target}ly regulating {t_n} by {reltype_str}'
                OQLquery = f'SELECT Relation WHERE objectType = ClinicalTrial AND \
NeighborOf (SELECT Entity WHERE objectType = ({indication_type})) AND NeighborOf ({get_modulators})'
                ct_g = self.process_oql(OQLquery,REQUEST_NAME)
                if isinstance(ct_g,ResnetGraph):
                    InhibitorsIndicationNetwork = ct_g.compose(InhibitorsIndicationNetwork)

                all_indications = InhibitorsIndicationNetwork.psobjs_with(only_with_values=self.params['indication_types'])
                print('Found %d indications on clinical trials with %d substances %sly regulating %s by %s relations' %  
                    (len(all_indications)- len(indication_by_regulation), len(modulators),with_effect_on_target,t_n,reltype_str))

        return modulators, all_indications


    def modulators_indications(self,linked_by:list,with_effect_on_target:str,min_refcount=1,drugs_only=False):
        return self.modulators_effects(linked_by,with_effect_on_target,True,min_refcount,drugs_only)


    def modulators_toxicities(self,linked_by:list,with_effect_on_target:str,min_refcount=1,drugs_only=False):
        return self.modulators_effects(linked_by,with_effect_on_target,False,min_refcount,drugs_only)


    def indications4chem_modulators(self,moa=ANY_MOA):
        indication_count = len(self.__indications())

        self.DirectAntagonists, indications = self.modulators_indications(['DirectRegulation'],'negative')
        if moa in [ANTAGONIST,ANY_MOA]:
            self.indications4antagonists.update(indications)
            new_indication_count = len(self.__indications()) - indication_count
            print('%d indications for drugs directy inhibiting %s were not found by previous searches' %  
                        (new_indication_count,self._input_targets()))
            
        self.DirectAgonists, indications = self.modulators_indications(['DirectRegulation'],'positive')
        if moa in [AGONIST,ANY_MOA]:
            self.indications4agonists.update(indications)
            new_indication_count = len(self.__indications()) - indication_count
            print('%d indications for drugs directy activating %s were not found by previous searches' %  
                        (new_indication_count,self._input_targets()))
            
        self.IndirectAgonists,_ = self.modulators_indications(
        ['Regulation','Expression','MolTransport'],'positive',min_refcount=1,drugs_only=True)

        self.IndirectAntagonists,_ = self.modulators_indications(
            ['Regulation','Expression','MolTransport'],'negative',min_refcount=1,drugs_only=True)
 

    def indications4cells_secreted_target(self,moa=ANY_MOA):
        if not self.target_class == 'Ligand': return

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

        t_n = self._input_targets()
        exist_indication_count = len(self.__indications())

        REQUEST_NAME = f'Find cells secreting {t_n}'
        OQLquery = 'SELECT Relation WHERE objectType = (CellExpression,MolTransport) AND NeighborOf ({select_targets}) AND NeighborOf (SELECT Entity WHERE objectType = Cell)'
        cells_make_target = self.process_oql(OQLquery.format(select_targets=self.oql4targets),REQUEST_NAME)
        if isinstance(cells_make_target,ResnetGraph):
            self.ProducingCells = cells_make_target._psobjs_with('CellType','ObjTypeName')
            print('Found %d cell types producing %s' % (len(self.ProducingCells),t_n))

        if self._is_strict(): # no further predictions of indications
            return # found self.ProducingCells will be used for target ranking
            

        OQLtemplate = 'SELECT Relation WHERE Effect = {effect} AND NeighborOf (SELECT Entity WHERE \
            objectType = ({indication_type})) AND NeighborOf (SELECT Entity WHERE id = ({cell_ids}))'
        
        ProducingCells_dbids = ResnetGraph.dbids(self.ProducingCells)
        if moa in [ANTAGONIST,ANY_MOA]:
            REQUEST_NAME = f'Find indications positively linked to cell secreting {t_n}'
            OQLquery = OQLtemplate.format(effect='positive',cell_ids=OQL.id2str(ProducingCells_dbids),indication_type=','.join(self.params['indication_types']))
            CellDiseaseNetwork = self.process_oql(OQLquery,REQUEST_NAME)
            if isinstance(CellDiseaseNetwork,ResnetGraph):
                indications = CellDiseaseNetwork.psobjs_with(only_with_values=self.params['indication_types'])
                if len(indications) > self.max_cell_indications:
                    indications = best_cell_indications(self.max_cell_indications,CellDiseaseNetwork)
                self.indications4antagonists.update(indications)

        if moa in [AGONIST,ANY_MOA]:
            REQUEST_NAME = f'Find indications negatively linked to cell secreting {t_n}'
            OQLquery = OQLtemplate.format(effect='negative',cell_ids=OQL.id2str(ProducingCells_dbids),indication_type=','.join(self.params['indication_types']))
            CellDiseaseNetwork = self.process_oql(OQLquery,REQUEST_NAME)
            if isinstance(CellDiseaseNetwork,ResnetGraph):
                indications = CellDiseaseNetwork.psobjs_with(only_with_values=self.params['indication_types'])
                if len(indications) > self.max_cell_indications:
                    indications = best_cell_indications(self.max_cell_indications,CellDiseaseNetwork)
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

        oqls = dict() # {tuple[target_name,partner_name,partner_dbid]:str}
        if self.partners:
            for partner in self.partners:
                for target_name, oql_query in target_oqls.items():
                    find_pathways_query = oql_query + ' AND ParentOf (SELECT Entity WHERE id = ({i}))'.format(i=str(partner.dbid()))
                    oqls[(target_name,partner.name(),partner.dbid())] = merged_pathways.format(select_networks=find_pathways_query)
                    # some ligands have duplicate names - must include dbid into tuple
        else:
            for target_name, oql_query in target_oqls.items():
                oqls[(target_name, '',0)] = merged_pathways.format(select_networks=oql_query)

        return oqls
    

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
            targets_regulome = merged_pathway.get_regulome(set(self.targets))
            self.PathwayComponents = set(targets_regulome._get_nodes())
            print (f'Found regulome with {len(self.PathwayComponents)} components')
            return targets_regulome
        else: 
            print('No curated pathways were found for %s' % self._input_targets())
            return ResnetGraph()


    def init_semantic_search(self):
        '''
        Loads
        -----
        RAWDF4ANTAGONISTS and RAWDF4AGONISTS df to raw_data
        '''
        print('\n\nInitializing semantic search')
        t_n = self._input_targets()

        indications4antagonists = self.__indications4antagonists()
        if indications4antagonists:
            indication_df = self.load_df(list(indications4antagonists),
                                         max_child_count=self.max_ontology_parent,
                                         max_threads=self.max_threads4ontology)
            indication_df._name_ = RAWDF4ANTAGONISTS
            self.add2raw(indication_df)
            print(f'Will score {len(indication_df)} indications for {t_n} antagonists')

        indications4agonists = self.__indications4agonists()
        if indications4agonists:
            indication_df = self.load_df(list(indications4agonists),
                                         max_child_count=self.max_ontology_parent,
                                         max_threads=self.max_threads4ontology)
            indication_df._name_ = RAWDF4AGONISTS
            self.add2raw(indication_df)
            print(f'Will score {len(indication_df)} indications for {t_n} agonists')
        
        if indications4antagonists or indications4agonists:
            return True
        else:
            print (f'No indications found for {t_n}')
            if self._is_strict(): print('Try setting strict_mode to False')
            return False
        


    def score_GVs(self, df2score:df):
        if not self.GVs: return
        t_n = self._input_targets()
        self.__colnameGV__ = t_n+' GVs'
        gvlinkcounter = 0
        if hasattr(self,'GVs2DiseaseGraph'):
            for i in df2score.index:
                target_dbids = set(df2score.at[i,self.__temp_id_col__])
                row_targets = set(self.Graph.psobj_with_dbids(target_dbids))
                GVscore = 0
                if isinstance(self.GVs2DiseaseGraph,ResnetGraph):
                    targetGVs = list(self.GVs2DiseaseGraph.get_neighbors(row_targets,self.GVs))
                    if len(targetGVs) > 0:
                        GVnames = [g.name() for g in targetGVs]
                        df2score.at[i,self.__colnameGV__] = ';'.join(GVnames)
                        gvlinkcounter += 1
                        gv_disease_subgraph = self.GVs2DiseaseGraph.get_subgraph(list(row_targets),targetGVs)
                        GVscore = len(gv_disease_subgraph.load_references())
                
                df2score.at[i,self._col_name_prefix+'GVs'] = GVscore
            print('Found %d indications linked to %d GVs' % (gvlinkcounter, len(self.GVs)))


    def __drug_connect_params(self,direct_modulators:bool,score_antagonists:bool):
        '''
        Return
        ------
        how2connect paramters to score substances having desired drug affect on the targets 
        '''

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
        '''
        Return
        ------
        how2connect parameters to score substances synergizing with target action on indication\n
        i.e. molecules having effect opposite to desired drug affect on targets
        '''
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


    def semscore4targets(self,in_worksheet:str,target_effect_on_indication:str):
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
        t_n = self._input_targets()
        target_in_header = t_n if len(t_n) < 45 else 'targets'
        colname = 'Activated by ' if target_effect_on_indication == 'positive' else 'Inhibited by '
        colname += target_in_header

        indication_df = self.raw_data[in_worksheet]
        score4antagonists = True if target_effect_on_indication == 'positive' else False

        if self._is_strict():
            booster_reltypes = ['Regulation','Biomarker','GeneticChange','QuantitativeChange','StateChange','FunctionalAssociation']
        else:
            booster_reltypes = ['Regulation','GeneticChange','FunctionalAssociation']
        how2connect = self.set_how2connect(['Regulation'],[target_effect_on_indication],'',booster_reltypes)
        linked_row_count,_,indication_df = self.link2concept(colname,self.targets,indication_df,how2connect)
        print('%d indications are %sly regulated by %s' % 
            (linked_row_count,target_effect_on_indication,t_n))

        self.score_GVs(indication_df)

        link_effect, drug_class, concepts = self.__drug_connect_params(True,score4antagonists)
        if concepts:
            # references suggesting that known drugs for the target as treatments for indication
            colname = target_in_header+' '+drug_class+' clin. trials'
            how2connect = self.set_how2connect(['ClinicalTrial'],[],'')
            linked_row_count,_,indication_df = self.link2concept(colname,concepts,indication_df,how2connect)
            print('Linked %d clinical trial indictions for %s %s' % 
                  (linked_row_count,t_n,drug_class))

            colname = target_in_header+' '+drug_class
            how2connect = self.set_how2connect(['Regulation'],[link_effect],'',['Regulation','FunctionalAssociation'])
            linked_row_count,_,indication_df = self.link2concept(colname,concepts,indication_df,how2connect)
            print('Linked %d indications for %s %s' % 
                  (linked_row_count,t_n,drug_class))

        #references reporting target agonists exacerbating indication or causing indication as adverse events
        link_effect, drug_class, concepts = self.__drug_tox_params(direct_modulators=True,
                                                                   score_antagonists=score4antagonists)
        if concepts:
            colname = target_in_header+' '+drug_class
            how2connect = self.set_how2connect(['Regulation'],[link_effect],'',['Regulation','FunctionalAssociation'])
            linked_row_count,_,indication_df = self.link2concept(colname,concepts,indication_df,how2connect)
            print('Linked %d indications as toxicities for %s %s' % 
                  (linked_row_count,t_n,drug_class))

    
        #references where target expression or activity changes in the indication
        colname = ' is upregulated' if target_effect_on_indication == 'positive' else ' is downregulated'
        colname = target_in_header + colname
        how2connect = self.set_how2connect(['QuantitativeChange'],[target_effect_on_indication],'',['Biomarker','StateChange','FunctionalAssociation'])
        linked_row_count,_,indication_df= self.link2concept(colname,self.targets,indication_df,how2connect)
        print('%d indications %sly regulate %s' % 
            (linked_row_count,target_effect_on_indication,t_n))

        #references suggesting target partners as targets for indication
        if self.partners:
            p_cl = self.partner_class if self.partner_class else 'partner'
            colname = f'{target_in_header} {p_cl}s'
            how2connect = self.set_how2connect(['Regulation'],[target_effect_on_indication],'',['Regulation','FunctionalAssociation'])
            linked_row_count,_,indication_df = self.link2concept(colname,self.partners,indication_df,how2connect)
            print('Linked %d indications for %d %s %ss' % 
                (linked_row_count,len(self.partners),t_n,p_cl))

        # references reporting that cells producing the target linked to indication  
        # only used if taregts are secretred ligands
        if hasattr(self, 'ProducingCells'):
            colname = f'{target_in_header} producing cells'
            how2connect = self.set_how2connect(['Regulation'],[target_effect_on_indication],'',['Regulation','FunctionalAssociation'])
            linked_row_count,_,indication_df = self.link2concept(colname,self.ProducingCells,indication_df,how2connect)
            print('Liked %d indications linked %d cells producing %s' % 
                  (linked_row_count,len(self.ProducingCells),t_n))

        link_effect, drug_class, concepts = self.__drug_connect_params(False,score4antagonists)
        if concepts:
            # references suggesting that known drugs for the target as treatments for indication
            colname = target_in_header+' '+drug_class+' clin. trials'
            # cloning session to avoid adding relations to self.Graph
            new_session = self._clone(to_retrieve=REFERENCE_IDENTIFIERS)
            how2connect = new_session.set_how2connect(['ClinicalTrial'])
            linked_row_count,_,indication_df = new_session.link2concept(colname,concepts,indication_df,how2connect)
            print('Linked %d clinical trial indications for %s %s' % 
                  (linked_row_count,t_n,drug_class))
            new_session.close_connection()

            colname = target_in_header+' '+drug_class
            # cloning session to avoid adding relations to self.Graph
            new_session = self._clone(to_retrieve=REFERENCE_IDENTIFIERS)
            how2connect = new_session.set_how2connect(['Regulation'],[link_effect],'',['Regulation'])
            linked_row_count,_,indication_df = new_session.link2concept(colname,concepts,indication_df,how2connect)
            print('Linked %d indications for %s %s' % (linked_row_count,t_n,drug_class))
            new_session.close_connection()

        #references reporting target agonists exacerbating indication or causing indication as adverse events
        link_effect, drug_class, concepts = self.__drug_tox_params(direct_modulators=False,
                                                                   score_antagonists=score4antagonists)
        if concepts:
            colname = target_in_header+' '+drug_class
            # cloning session to avoid adding relations to self.Graph
            new_session = self._clone(to_retrieve=REFERENCE_IDENTIFIERS)
            how2connect = new_session.set_how2connect(['Regulation'],[link_effect],'',['Regulation'])
            linked_row_count,_,indication_df = new_session.link2concept(colname,concepts,indication_df,how2connect)
            print('Linked %d indications as toxicities for %s %s' % (linked_row_count,t_n,drug_class))
            new_session.close_connection()
        
        if hasattr(self, 'PathwayComponents'):
            #references linking target pathway to indication
            colname = target_in_header + ' pathway components'
            # cloning session to avoid adding relations to self.Graph
            new_session = self._clone(to_retrieve=REFERENCE_IDENTIFIERS)
            how2connect = new_session.set_how2connect(['Regulation'],[target_effect_on_indication],step=125)
            linked_row_count,_,indication_df = new_session.link2concept(colname,list(self.PathwayComponents),indication_df,how2connect)
            print('Linked %d indications to %s pathway components' % (linked_row_count,t_n))
            new_session.close_connection()
 
        if not indication_df.empty:
            counts_df = self.make_count_df(indication_df,with_name=in_worksheet)
            self.add2raw(counts_df)
        return


    def other_effects(self)->ResnetGraph:
        # need to be called after ranking to subtract self.all_entity_ids
        print('Findind indication linked with unknown effect to %s' % self._input_targets())
        old_rel_props = self.relProps
        self.add_rel_props(PS_SENTENCE_PROPS+list(PS_BIBLIO_PROPS))
        t_n = self._input_targets()
        ranked_indication_ids = self.Graph.dbids4nodes(self.params['indication_types'])
        request_name = f'Find indications modulated by {t_n} with unknown effect'
        oql4indications_type = self.__oql4indications_type()
        OQLquery = f'SELECT Relation WHERE objectType = (Regulation,QuantitativeChange) AND Effect = unknown AND \
NeighborOf({self.oql4targets}) AND NeighborOf ({oql4indications_type})'
        to_return = self.process_oql(OQLquery,request_name)
        
        request_name = f'Find indications where {t_n} was suggested as Biomarker'
        OQLquery = f'SELECT Relation WHERE objectType = (Biomarker,StateChange,GeneticChange) AND \
NeighborOf({self.oql4targets}) AND NeighborOf ({oql4indications_type})'
        if isinstance(to_return,ResnetGraph):
            bm_g = self.process_oql(OQLquery,request_name)
            if isinstance(bm_g,ResnetGraph):
                to_return = bm_g.compose(to_return)
        
        request_name = f'Find indications with genetically linked {t_n} Genetic Variants'
        OQLquery = f'SELECT Relation WHERE objectType = FunctionalAssociation AND NeighborOf ({oql4indications_type})'
        OQLquery += ' AND NeighborOf (SELECT Entity WHERE id = ({ids}))'
        GVdbids = set(ResnetGraph.dbids(self.GVs))
        if isinstance(to_return,ResnetGraph):
            gv_g = self.iterate_oql(OQLquery,GVdbids,request_name=request_name)
            to_return = gv_g.compose(to_return)

            to_return.remove_nodes_from(ranked_indication_ids) # now delete all indication with known effect
            indications = to_return.psobjs_with(only_with_values=self.params['indication_types'])

            print('Found %d new indications linked to %s with unknown effect' % (len(indications),self._input_targets()))
        self.relProps = old_rel_props
        assert(isinstance(to_return,ResnetGraph))
        return to_return


    def load_indications4targets(self,moa=ANY_MOA):
        '''
        Input
        -----
        target names must be in self.params['target_names']
        moa in [ANTAGONIST, AGONIST, ANY_MOA]
        '''
        assert(moa in [ANTAGONIST, AGONIST, ANY_MOA])
        start_time = time.time()
        self.set_targets()
        self.set_partners()

        self.find_target_indications(moa)
        #return  #uncomment to enable fast downstream testing
        self.indications4partners(moa)

        if self.target_class == 'Ligand':
            self.indications4cells_secreted_target(moa)
        else:
            self.indications4chem_modulators(moa)

        self.get_pathway_componets()
        self.__resolve_conflict_indications()
        print("%d indications for %s were found in %s" % 
        (len(self.__indications()), self._input_targets(), execution_time(start_time)))

        return self.__all_indications()


    def add_ps_bibliography(self,suffix='',add_graph=ResnetGraph()):
        targets_neighbors = self.Graph.neighborhood(set(self.targets))
        targets_neighbors = self.Graph.neighborhood(set(self.GVs)).compose(targets_neighbors)
        targets_neighbors = self.Graph.neighborhood(set(self.partners)).compose(targets_neighbors)
        if add_graph:
            targets_neighbors = add_graph.compose(targets_neighbors)
        super().add_ps_bibliography(suffix,from_graph=targets_neighbors)


    def other_effects2df(self):
        '''
        Adds
        ----
        UNKNOWN_EFFECT_INDICATIONS worksheet to self.report_data with snippets for unknown effect indications
        '''
        other_effects_graph = self.other_effects()
        other_indications = other_effects_graph.snippets2df(df_name=UNKEFFECTDF)
        self.add2report(other_indications)


    def perform_semantic_search(self):
        '''
        Create
        ------
        worksheets in self.report_data: ANTAGONISTSDF,AGONISTSDF,UNKEFFECTDF
        '''
        start = time.time()
        # cannot multithread here yet - self.Graph is mutating.  Need to have self.clone function
        self.semscore4targets(RAWDF4ANTAGONISTS,'positive')
        self.semscore4targets(RAWDF4AGONISTS,'negative')
        self.other_effects2df()

        empty_cols2drop = list(self.params.get('drop_empty_columns_from_report',[]))
        dfnames_map = self.__dfnames_map()
        for raw_dfname, report_dfname in dfnames_map.items():
            self.normalize(raw_dfname,report_dfname,drop_empty_columns=empty_cols2drop)

        with ThreadPoolExecutor(max_workers=5, thread_name_prefix='AddAnnot') as b:
            id_path_futures = list()
            for report_dfname in dfnames_map.values():
                b.submit(self.add_ontology_df,report_dfname)
                id_path_futures.append((report_dfname,b.submit(self.id2paths,report_dfname)))

            for report_dfname, future in id_path_futures:
                id2paths = future.result()
                self.report_pandas[report_dfname] = self.report_pandas[report_dfname].merge_dict(id2paths,'Ontology parents','Name')
            
            b.shutdown()
        print(f'TargetIndications semantic search is finished in {execution_time(start)}')
        return 


    def make_report(self):
        start_time = time.time()
        self.flush_dump_files()
        all_indications = self.load_indications4targets()
        
        if self.params['add_bibliography']:
            indication_names = [n.name() for n in all_indications]
            df4etm = df.from_dict({'Name':indication_names})
            e = ThreadPoolExecutor(thread_name_prefix='ETMannot')
            etm_future = e.submit(self.etm_refs2df,df4etm,self.input_names(),'Name',[],len(df4etm))
        
        if self.init_semantic_search():
            self.perform_semantic_search()

            if self.params['add_bibliography']:
                indication_etmrefs = etm_future.result()
                etm_ref_colname = self.etm_ref_column_name('Name',self.input_names())
                doi_ref_colname = self.etm_doi_column_name('Name',self.input_names())
                for worksheet_name in self.__dfnames_map().values():
                    self.report_pandas[worksheet_name] = self.report_pandas[worksheet_name].merge_df(indication_etmrefs,how='left',on='Name',columns=[etm_ref_colname,doi_ref_colname])
                self.add_etm_bibliography()
                self.add_ps_bibliography()

            print(f'{self.report_path()} repurposing is done in {execution_time(start_time)}')
        else:
            print('Failed to initialize semantic search for TargetIndications')


    def write_report(self):
        report = pd.ExcelWriter(self.report_path(), engine='xlsxwriter')
        ordered_worksheets = list()
        for report_worksheet_name in self.__dfnames_map().values():
            ordered_worksheets.append(report_worksheet_name)
            ordered_worksheets.append('ontology4'+report_worksheet_name)
        ordered_worksheets.append(UNKEFFECTDF)
        
        self.add2writer(report,df_names=ordered_worksheets)
        report.close()


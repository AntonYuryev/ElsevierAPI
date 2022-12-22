from .SemanticSearch import SemanticSearch,ResnetGraph, len,df,REFERENCE_IDENTIFIERS,COUNTS,ONTOLOGY_ANALYSIS
from .ResnetGraph import REFCOUNT
from ..ETM_API.references import PS_SENTENCE_PROPS, EFFECT,PS_BIBLIO_PROPS
from .PathwayStudioGOQL import OQL
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
    DirectAntagonistsIDs = list()
    DirectAgonistsIDs = list()
    IndirectAgonistsIDs = list()
    IndirectAntagonistsIDs = list()
    GVids = list()


    def __init__(self, APIconfig, params={}):
        super().__init__(APIconfig)
        self.add_rel_props([EFFECT])
        self.partners = list() # list of PSObject
        self.partners_ids = list() #
        self.columns2drop += [self.__resnet_name__,self.__mapped_by__]

        if params: 
            self.param = dict(params)
        else: # default params:
            self.param = {'partner_names':[],
    # if partner_names is empty script will try finding Ligands for Receptor targets and Receptors for Ligand targets
                 'partner_class':'', # use it only if partner_names not empty
                 'indication_types': ['Disease','Virus','Pathogen'], #['Disease','Virus','CellProcess']
                 'target_names':[],
                 'target_type':'',
                 'pathway_name_must_include_target':True,
    # if 'pathway_name_must_include_target' True only pathways depicting target signaling are considered for regulome construction
    # if 'pathway_name_must_include_target' False all pathways containing both Targets and Partners will be considered
                 'strict_mode':True,
        # if strict mode algorithm only ranks indications suggetsed in the literarure without predicting new indications
        # if not strict mode algorithm also predicts and rank additional indications based on target expression or genetic profiles in disease
                 'data_dir':'',
                 'add_bibliography' : True
                }

        self.set_dir(self.param['data_dir'])
    
    
    def _target_names(self):
        """
        Returns target names as they were entered in configuration.
        return values may be different from self.input_names()
        """
        return ','.join(self.param['target_names'])

    def input_names(self):
        """
        Returns
        -------
        entity names that can be used both in PS and ETM
        """
        #input names are used for finding references in ETM.
        # RepurposeDrug has it own 'input_names'.
        return [x['Name'][0] for x in self.Drug_Targets]

    def __known_effect_indications(self):
        known_effect_indication = self.indications4agonists | self.indications4antagonists
        return known_effect_indication

    def _indications(self):
        """
        Returns
        -------
        all found indications depending on self.param['strict_mode']
        """
        if self.param['strict_mode']:
            return self.indications4antagonists_strict | self.indications4agonists_strict
        else:
            return self.indications4antagonists | self.indications4agonists

    def _indications4antagonists(self):
        if self.param['strict_mode']:
            return self.indications4antagonists_strict
        else:
            return self.indications4antagonists

    def _indications4agonists(self):
        if self.param['strict_mode']:
            return self.indications4agonists_strict
        else:
            return self.indications4agonists

    def __is_strict(self):
        return self.param['strict_mode']
  
    def set_targets(self):
        try:
            target_objtype = str(self.param['target_objtype'])
        except KeyError: 
            target_objtype = 'Protein'

        prop_names_str, prop_values_str = OQL.get_search_strings(['Name'],self.param['target_names'])
        self.add_ent_props(['Class'])
        self.oql4targets = 'SELECT Entity WHERE ({prop_name}) = ({values}) AND objectType = '+target_objtype 
        targets_graph = self.process_oql(self.oql4targets.format(prop_name=prop_names_str,values=prop_values_str))
        
        self.Drug_Targets = targets_graph._get_nodes()
        self.target_ids = [x['Id'][0] for x in self.Drug_Targets]

        self.oql4targets = OQL.get_objects(self.target_ids)
        for t in self.Drug_Targets:
            try:
                self.target_class = t['Class'][0]
                break
            # assumes all targets have the same class
            except KeyError: continue


    def set_partner_class(self):
        #finding effect for complementary receptor or ligands
        if self.target_class == 'Ligand': self.partner_class = "Receptor"
        elif self.target_class == 'Receptor': self.partner_class = "Ligand"
        else: self.target_class = NotImplemented


    def find_partners(self):
        """
        Finds
        -----
        Ligands for receptors, Receptors for Ligand\
            linked to target with Effect=positive
        """
        self.set_partner_class()
        SELECTpartners = 'SELECT Entity WHERE Class = {partner_class} AND objectType = Protein AND Connected by (SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive) to ({select_target})'
        partners_graph = self.process_oql(SELECTpartners.format(partner_class=self.partner_class, select_target=self.oql4targets))
        if self.partner_class == 'Ligand':
            SELECTsecretedpartners = 'SELECT Entity WHERE "Cell Localization" = Secreted AND objectType = Protein AND Connected by (SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive) to ({select_target})'
            secreted_partners = self.process_oql(SELECTsecretedpartners.format(select_target=self.oql4targets))
            partners_graph = partners_graph.add_graph(secreted_partners)

        self.partners_ids = partners_graph.get_node_ids(['Protein'])
        target_names = self._target_names()
        partners_str = ','.join(n['Name'][0] for n in partners_graph._get_nodes())
        print ('Found %s %ss as partners for %s'%(partners_str, self.partner_class,target_names))
        return self.Graph._get_nodes(self.partners_ids)


    def set_partners(self):
        """
        Assumes partners are linked to Targets with Effect=positive
        """
        try:
        # case when self.param has 'partner_names' and 'partner_class' explicitly
        # use it when receptor partners are metabolites and therefore cannot be found by select_partners()
            partner_names = self.param['partner_names']
            try:
                self.partner_class = self.param['partner_class']
            except KeyError:
                print ('"partner_class" parameter is not specified !!!')
                self.partner_class = ''
        except KeyError:
            partner_names = ''

        if partner_names: 
                partners_str = OQL.join_with_quotes(partner_names)
                SELECTpartners = 'SELECT Entity WHERE Name = ({partner_names})'.format(partner_names=partners_str)
                partners_graph = self.process_oql(SELECTpartners)
                self.partners = partners_graph._get_nodes()
                self.partners_ids = [p['Id'][0] for p in self.partners]
                self.find_partners_oql = OQL.get_objects(self.partners_ids)
        else:
            self.find_partners()

        
    def _get_report_name(self):
        indics = ','.join(self.param['indication_types'])
        rep_pred = 'suggested ' if self.param['strict_mode'] else 'suggested,predicted ' 
        return str(self.param['data_dir']+rep_pred+ indics+' for '+ self._target_names())


    def GVindications(self):
        target_names = [x['Name'][0] for x in self.Drug_Targets]
        t_n = ','.join(target_names)
        select_targets = OQL.get_objects(self.target_ids)
  
        selectGVs = 'SELECT Entity WHERE objectType = GeneticVariant AND Connected by (SELECT Relation WHERE objectType = GeneticChange) to ({select_target})'
        selectGVs = selectGVs.format(select_target=select_targets)
        REQUEST_NAME = 'Find indications linked to {target} genetic variants'.format(target=t_n)
        OQLquery = 'SELECT Relation WHERE NeighborOf({select_gvs}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_type}))'
        self.GVsInDiseaseNetwork = self.process_oql(OQLquery.format(select_gvs=selectGVs,indication_type=','.join(self.param['indication_types'])),REQUEST_NAME)
        indication_ids = self.GVsInDiseaseNetwork.get_node_ids(self.param['indication_types'])
        self.GVids = self.GVsInDiseaseNetwork.get_node_ids(['GeneticVariant'])
        
        print('Found %d indications genetically linked to %d Genetic Variants in %s' % 
            (len(indication_ids), len(self.GVids), t_n))
        
        return indication_ids


  #  def set_strictmode_indications(self):
        # this function is overidden in RepurposeDrug
  #      self.indications4strictmode = self.target_indications4strictmode


    def __resolve_conflict_indications(self):

        def __resolve(all_rels:list, using_rel_type:str):
            only_rels_with_type = [r for r in all_rels if r.objtype() == using_rel_type]
            if only_rels_with_type:
                only_rels_with_type.sort(key=lambda x: x.get_reference_count(), reverse=True)
                best_effect = only_rels_with_type[0]['Effect'][0]
                if best_effect == 'positive':
                    self.indications4agonists.remove(i)
                    return True
                elif best_effect == 'negative':
                    self.indications4antagonists.remove(i)
                    return True
                else:
                    return False

            return False

        conflict_ids = self.indications4antagonists.intersection(self.indications4agonists)
        unresolved_ids = set()
        for i in conflict_ids:
            target_rels = list(self.Graph.find_relations(self.target_ids,[i]))
            if not __resolve(target_rels,'Regulation'):
                if not __resolve(target_rels,'QuantitativeChange'):
                    partner_rels = list(self.Graph.find_relations(self.partners_ids,[i]))
                    if not __resolve(partner_rels,'Regulation'):
                        if not __resolve(partner_rels,'QuantitativeChange'):
                            unresolved_ids.add(i)

        print('%d indications cannot be resolved.\n They will appear in both worksheets for agonists and antagonist:' % 
                    len(unresolved_ids))
        for i in unresolved_ids:
            try:
                print(self.Graph._get_node(i).name())
            except KeyError:
                continue


    def find_target_indications(self):
        f_t = self.oql4targets
        t_n = self._target_names()
        # finding indications for antagonists
        effect = 'positive'
        select_indications_by_type = 'SELECT Entity WHERE objectType = ({indication_type})'.format(indication_type=','.join(self.param['indication_types']))
        
        REQUEST_NAME = 'Find indications {effect}ly modulated by {target}'.format(effect=effect, target=t_n)
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = {effect} AND \
            NeighborOf({target}) AND NeighborOf ({indications})' 
        OQLquery =  OQLquery.format(target=f_t,effect = effect,indications=select_indications_by_type)  
        ModulatedByTargetNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        indication_ids = ModulatedByTargetNetwork.get_node_ids(self.param['indication_types'])
        self.indications4antagonists.update(indication_ids)
        self.indications4antagonists_strict.update(indication_ids)
        #self.target_indications4strictmode = set(indication_ids)
        print('Found %d diseases %sly regulated by target %s' % (len(indication_ids),effect,t_n))

        REQUEST_NAME = 'Find indications {effect}ly regulating {target}'.format(effect=effect, target=t_n)
        OQLquery = 'SELECT Relation WHERE objectType = QuantitativeChange AND Effect = {effect} AND \
            NeighborOf ({target}) AND NeighborOf ({indications})'
        OQLquery = OQLquery.format(target=f_t,effect = effect,indications=select_indications_by_type)
        ActivatedInDiseaseNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        indication_ids = ActivatedInDiseaseNetwork.get_node_ids(self.param['indication_types'])
        self.indications4antagonists.update(indication_ids)
        print('Found %d diseases where %s is %sly regulated' % (len(indication_ids),t_n,effect))

        # finding indications for agonists
        effect =  'negative'
        REQUEST_NAME = 'Find indications {effect}ly modulated by {target}'.format(effect=effect, target=t_n)
        #oql2select_indications = 'SELECT Entity WHERE objectType = ({indication_type})'.format(indication_type=','.join(self.param['indication_types']))
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = {effect} AND NeighborOf({select_target}) AND NeighborOf ({indications})'      
        OQLquery = OQLquery.format(select_target=f_t,effect=effect,indications=select_indications_by_type)
        ModulatedByTargetNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        indication_ids = ModulatedByTargetNetwork.get_node_ids(self.param['indication_types'])
        self.indications4agonists.update(indication_ids)
        self.indications4agonists_strict.update(indication_ids)
        #self.target_indications4strictmode = set(indication_ids)
        print('Found %d diseases %sly regulated by %s' % (len(indication_ids),effect,t_n))

        REQUEST_NAME = 'Find indications {effect}ly regulating {target}'.format(effect=effect, target=t_n)
        OQLquery = 'SELECT Relation WHERE objectType = QuantitativeChange AND Effect = {effect} AND NeighborOf ({select_target}) AND NeighborOf ({indications})'  
        OQLquery = OQLquery.format(select_target=f_t,effect = effect,indications=select_indications_by_type)
        ActivatedInDiseaseNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        indication_ids = ActivatedInDiseaseNetwork.get_node_ids(self.param['indication_types'])
        self.indications4agonists.update(indication_ids)
        print('Found %d diseases where %s is %sly regulated' % (len(indication_ids),t_n,effect))

        # finding indications with unknown effect
        indication_ids = self.GVindications()
        new_indications = set(indication_ids).difference(self.__known_effect_indications())
        self.unknown_effect_indications_strict.update(new_indications)
        print('%d indications linked to %s genetic variations were not found by previous searches' % (len(new_indications),t_n))

        REQUEST_NAME = 'Find indications where {target} is biomarker'.format(target=t_n)
        OQLquery = 'SELECT Relation WHERE objectType = Biomarker AND NeighborOf({select_target}) AND NeighborOf ({indications})'
        OQLquery = OQLquery.format(select_target=f_t,indications=select_indications_by_type)
        BiomarkerInDiseaseNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        indication_ids = BiomarkerInDiseaseNetwork.get_node_ids(self.param['indication_types'])
        print('Found %d indications where target %s is claimed as biomarker' %  (len(indication_ids),t_n))
        new_indications = set(indication_ids).difference(self.__known_effect_indications())
        print('%d indications having %s as biomarker were not found by previous searches' % (len(new_indications),t_n))
        self.unknown_effect_indications.update(new_indications)

        self.oql2select_indications = 'SELECT Entity WHERE id = ({ids})'
        

    def modulators_effects(self,linked_by:list(),with_effect_on_target:str,find_indications=True,min_refcount=1,drugs_only=False):
        f_t = self.oql4targets
        t_n = self._target_names()
        reltype_str = ','.join(linked_by)

        REQUEST_NAME = 'Find substances {effect}ly regulating {target} by {rt}'.format(effect=with_effect_on_target,target=t_n,rt=reltype_str)
        OQLquery = 'SELECT Relation WHERE objectType = ({rel_type}) AND Effect = {effect} AND \
            NeighborOf upstream ({select_target})'
            
        OQLquery = OQLquery.format(rel_type=reltype_str,select_target=f_t,effect=with_effect_on_target)
        if drugs_only:
            nb = ' AND NeighborOf downstream ({})'.format(OQL.get_childs(['drugs'],['Name']))
        else:
            nb = ' AND NeighborOf downstream (SELECT Entity WHERE objectType = SmallMol)'

        OQLquery += nb

        if min_refcount > 1:
             OQLquery += ' AND '+REFCOUNT+' >= '+str(min_refcount)

        TargetInhibitorsNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        modulators_ids = TargetInhibitorsNetwork.get_node_ids(['SmallMol'])
        print('Found %d substances %sly regulating %s by %s' % (len(modulators_ids),with_effect_on_target,t_n,reltype_str))
        
        # now find indications for modulators found on previous step
        if modulators_ids:
            indication_type=','.join(self.param['indication_types'])
            modulators=OQL.get_objects(modulators_ids)
            effect_type = 'negative' if find_indications else 'positive'

            REQUEST_NAME = 'Find indications for substances {effect}ly regulating {target} by {rt}'.format(effect=with_effect_on_target,target=t_n,rt=reltype_str)
            OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = {eftype} AND \
                NeighborOf (SELECT Entity WHERE objectType = ({inditypes})) AND NeighborOf ({modulators})'
            OQLquery = OQLquery.format(eftype=effect_type, inditypes=indication_type,modulators=modulators)
            InhibitorsIndicationNetwork = self.process_oql(OQLquery,REQUEST_NAME)
            indication_ids = InhibitorsIndicationNetwork.get_node_ids(self.param['indication_types'])
            print('Found %d indications for %d substances %sly regulating %s by %s relations' %  
                (len(indication_ids), len(modulators_ids),with_effect_on_target,t_n,reltype_str))

            REQUEST_NAME = 'Find clinical trials for substances {effect}ly regulating {target} by {rt}'.format(
                                                        effect=with_effect_on_target,target=t_n,rt=reltype_str)
            OQLquery = 'SELECT Relation WHERE objectType = ClinicalTrial AND \
                NeighborOf (SELECT Entity WHERE objectType = ({inditypes})) AND NeighborOf ({modulators})'
            OQLquery = OQLquery.format(inditypes=indication_type,modulators=modulators)
            InhibitorsIndicationNetwork.add_graph(self.process_oql(OQLquery,REQUEST_NAME))

            indication_by_regulation = indication_ids
            indication_ids = InhibitorsIndicationNetwork.get_node_ids(self.param['indication_types'])
            print('Found %d indications on clinical trials with %d substances %sly regulating %s by %s relations' %  
                (len(indication_ids)- len(indication_by_regulation), len(modulators_ids),with_effect_on_target,t_n,reltype_str))

        return modulators_ids, indication_ids


    def modulators_indications(self,linked_by:list(),with_effect_on_target:str,min_refcount=1,drugs_only=False):
        return self.modulators_effects(linked_by,with_effect_on_target,True,min_refcount,drugs_only)


    def modulators_toxicities(self,linked_by:list(),with_effect_on_target:str,min_refcount=1,drugs_only=False):
        return self.modulators_effects(linked_by,with_effect_on_target,False,min_refcount,drugs_only)


    def indications4chem_modulators(self):
        exist_indication_count = len(self._indications())
        self.DirectAntagonistsIDs, indication_ids = self.modulators_indications(['DirectRegulation'],'negative')
        self.indications4antagonists.update(indication_ids)
        new_indication_count = len(self._indications()) - exist_indication_count
        print('%d indications for drugs directy inhibiting %s were not found by previous searches' %  
                    (new_indication_count,self._target_names()))
        self.IndirectAgonistsIDs, indication_ids = self.modulators_indications(
        ['Regulation','Expression','MolTransport'],'positive',min_refcount=1,drugs_only=True)

        self.DirectAgonistsIDs, indication_ids = self.modulators_indications(['DirectRegulation'],'positive')
        self.indications4agonists.update(indication_ids)
        self.IndirectAntagonistsIDs, indication_ids = self.modulators_indications(
            ['Regulation','Expression','MolTransport'],'negative',min_refcount=1,drugs_only=True)

 
    def indications4partners(self):
        """
        Assumes partners are linked to Target with Effect=positive
        """
        if not self.partners: return

        t_n = self._target_names()
        exist_indication_count = len(self._indications())

        REQUEST_NAME = 'Find indications positively regulated by {partner}s of {targets}'.format(targets=t_n,partner=self.partner_class.lower())
        OQLtemplate = 'SELECT Relation WHERE objectType = (Regulation,QuantitativeChange) AND \
                Effect = {eff} AND NeighborOf ({partners}) AND NeighborOf ({indications})'

        OQLquery = OQLtemplate.format(eff='positive',partners=self.find_partners_oql,indications=self.oql2select_indications)
        self.PartnerIndicationNetwork = self.iterate_oql(OQLquery,self._indications(),request_name=REQUEST_NAME)
        indication_ids = self.PartnerIndicationNetwork.get_node_ids(self.param['indication_types'])
        self.indications4antagonists.update(indication_ids)
        print('Found %d indications for %d %s %ss' %  
                (len(indication_ids), len(self.partners),t_n,self.partner_class.lower()))

        REQUEST_NAME = 'Find indications negatively regulated by {partner}s of {targets}'.format(targets=t_n,partner=self.partner_class.lower())
        OQLquery = OQLtemplate.format(eff='negative',partners=self.find_partners_oql,indications=self.oql2select_indications)
        self.PartnerIndicationNetwork = self.iterate_oql(OQLquery,self._indications(),request_name=REQUEST_NAME)
        indication_ids = self.PartnerIndicationNetwork.get_node_ids(self.param['indication_types'])
        self.indications4antagonists.update(indication_ids)

        new_indication_count = len(self._indications()) - exist_indication_count
        print('%d indications for %d %s %ss were not found by previous searches' %  
                (new_indication_count, len(self.partners),t_n,self.partner_class.lower()))
        

    def indications4cells_secreted_target(self):
        if not self.target_class == 'Ligand': return

        t_n = self._target_names()
        exist_indication_count = len(self._indications())

        REQUEST_NAME = 'Find indications linked to cells secreting the {target}'.format(target=t_n)
        OQLquery = 'SELECT Relation WHERE objectType = (CellExpression,MolTransport) AND NeighborOf ({select_targets}) AND NeighborOf (SELECT Entity WHERE objectType = Cell)'
        cells_make_target = self.process_oql(OQLquery.format(select_targets=self.oql4targets),REQUEST_NAME)
        self.ProducingCellsIDs = cells_make_target.get_node_ids(['CellType'])
        print('Found %d cell types producing %s' % (len(self.ProducingCellsIDs),t_n))

        if self.param['strict_mode']: return # no further predictions
        REQUEST_NAME = 'Find indications linked to cell secrteing {target}'
        REQUEST_NAME = REQUEST_NAME.format(target = t_n)

        OQLtemplate = 'SELECT Relation WHERE Effect = {effect} AND NeighborOf (SELECT Entity WHERE objectType = ({indication_type})) AND NeighborOf (SELECT Entity WHERE id = ({cell_ids}))'
        
        OQLquery = OQLtemplate.format(effect='positive',cell_ids=OQL.id2str(self.ProducingCellsIDs),indication_type=','.join(self.param['indication_types']))
        CellDiseaseNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        indication_ids = CellDiseaseNetwork.get_node_ids(self.param['indication_types'])
        self.indications4antagonists.update(indication_ids)

        OQLquery = OQLtemplate.format(effect='negative',cell_ids=OQL.id2str(self.ProducingCellsIDs),indication_type=','.join(self.param['indication_types']))
        CellDiseaseNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        indication_ids = CellDiseaseNetwork.get_node_ids(self.param['indication_types'])
        self.indications4antagonists.update(indication_ids)
            
        new_indication_count = len(self._indications()) - exist_indication_count
        print('%d indications for %d cells producing %s were not found by previous searches' %  
                (new_indication_count, len(self.ProducingCellsIDs),t_n))


    def pathway_oql(self):
        #set pathway_name_must_include_target to True if targets have a lot of large curated pathways
        target_oqls = dict()
        pct = '%'
        merged_pathways = 'SELECT Relation WHERE objectType = (DirectRegulation,Binding,ProtModification,PromoterBinding,ChemicalReaction) AND MemberOf ({select_networks})'

        for target in self.Drug_Targets:
            target_id = target['Id'][0]
            target_name = target['Name'][0]
            SELECTpathways = 'SELECT Network WHERE ParentOf (SELECT Entity WHERE id = ({i}))'.format(i=str(target_id))
            if self.param['pathway_name_must_include_target']:
                SELECTpathways = SELECTpathways + ' AND Name LIKE (\''+pct+target_name+pct+'\')' #additional refinement for popular targets

            target_oqls[target_name] = SELECTpathways

        if self.partners_ids:
            target_partner_oqls = dict()
            for partner in self.partners:
                partner_id = partner['Id'][0]
                partner_name = partner['Name'][0]
                for target_name, oql_query in target_oqls.items():
                    find_pathways_query = oql_query + ' AND ParentOf (SELECT Entity WHERE id = ({i}))'.format(i=str(partner_id))
                    target_partner_oqls[(target_name,partner_name)] = merged_pathways.format(select_networks=find_pathways_query)

            return target_partner_oqls
        else:
            tuple_target_oqls = dict()
            for target_name, oql_query in target_oqls.items():
                tuple_target_oqls[(target_name, '')] = merged_pathways.format(select_networks=oql_query)

            return tuple_target_oqls
    

    def get_pathway_componets(self):
        #finding downstream pathway components
        oql_queries = self.pathway_oql() # separate oql_query for each target
        REQUEST_NAME = 'Find curated pathways containing {targets}'
        merged_pathway = ResnetGraph()
        component_names_with_regulome = set()
        for components_tuple, oql_query in oql_queries.items():
            target_name = components_tuple[0]
            partner_name = components_tuple[1]
            REQUEST_NAME = REQUEST_NAME.format(targets=target_name)
            if partner_name:
                REQUEST_NAME = REQUEST_NAME + ' and ' + partner_name

            pathway_with_target = self.process_oql(oql_query,REQUEST_NAME)
            if pathway_with_target:
                merged_pathway.add_graph(pathway_with_target)
                component_names_with_regulome.add(target_name)
                if partner_name: component_names_with_regulome.add(partner_name)

        if merged_pathway:
            targets_regulome = merged_pathway.get_regulome(self.target_ids)
            self.PathwayComponentsIDs = set()
            self.PathwayComponentsIDs.update(list(targets_regulome.nodes))

            t_n = ','.join(component_names_with_regulome)
            print ('Found %s regulome with %d components' %(t_n,len(self.PathwayComponentsIDs)))
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
        IndicationNames4antagonists = [y[0] for i,y in self.Graph.nodes(data='Name') if i in indications4antagonists]

        if IndicationNames4antagonists:
            indications2load = df.from_dict({'Name':IndicationNames4antagonists})
            indication_df = self.load_pandas(indications2load,prop_names_in_header=True,map2type=self.param['indication_types'])
            indication_df._name_ = DF4ANTAGONISTS
            self.add2raw(indication_df)
            print('Will score %d indications for antagonists of %s' % (len(indication_df),t_n))
        
        indications4agonists = self._indications4agonists()
        IndicationNames4agonists = [y[0] for i,y in self.Graph.nodes(data='Name') if i in indications4agonists]

        if IndicationNames4agonists:
            indications2load = df.from_dict({'Name':IndicationNames4agonists})
            indication_df = self.load_pandas(indications2load,prop_names_in_header=True,map2type=self.param['indication_types'])
            indication_df._name_ = DF4AGONISTS
            self.add2raw(indication_df)
            print('Will score %d indications for agonists of %s' % (len(indication_df),t_n))
        
        if IndicationNames4antagonists or IndicationNames4agonists:
            return True
        else:
            print ('No indications found for %s' % t_n)
            if self.param['strict_mode']:
                print('Try setting strict_mode to False')
            return False
        

    def score_GVs(self, df2score:df):
        if not self.GVids: return
        t_n = self._target_names()
        self.__colnameGV__ = t_n+' GVs'
        gvlinkcounter = 0
        for i in df2score.index:
            row_entity_ids = list(df2score.at[i,self.__temp_id_col__])
            if hasattr(self,'GVsInDiseaseNetwork'):
                GVneighborsId = self.GVsInDiseaseNetwork.get_neighbors(row_entity_ids, self.GVids)
                
            GVscore = 0
            if len(GVneighborsId) > 0:
                GVnames = set([n[0] for i,n in self.Graph.nodes.data('Name') if i in GVneighborsId])
                df2score.at[i,self.__colnameGV__] = ';'.join(GVnames)
                gvlinkcounter += 1
                gv_disease_subgraph = self.GVsInDiseaseNetwork.get_subgraph(row_entity_ids,GVneighborsId)
                GVscore = len(gv_disease_subgraph.load_references())
            
            df2score.at[i,self._col_name_prefix+'GVs'] = GVscore

        print('Found %d indications linked to %d GVs' % (gvlinkcounter, len(self.GVids)))


    def __drug_connect_params(self,direct_modulators:bool,score_antagonists:bool):
        # function returns how2connect paramters to score molecules similar to desired drug affect on the targets 
        effect = None
        drug_class = None
        concept_ids = list()

        if score_antagonists:
        # most common case when targets must be inhibited
            if direct_modulators:
                    effect = 'negative'
                    drug_class = 'direct antagonists'             
                    concept_ids = self.DirectAntagonistsIDs
            else:
                    effect = 'negative'
                    drug_class = 'indirect antagonists'             
                    concept_ids = self.IndirectAntagonistsIDs
        else:
        # case if drug are agonists
            if direct_modulators:
                    effect = 'positive'
                    drug_class = 'direct agonists'
                    concept_ids = self.DirectAgonistsIDs
            else:
                    effect = 'positive'
                    drug_class = 'indirect agonists'
                    concept_ids = self.IndirectAgonistsIDs

        return effect, drug_class, concept_ids


    def __drug_tox_params(self,direct_modulators:bool,score_antagonists:bool):
        effect = None
        drug_class = None
        concept_ids = list()
        # function returns how2connect parameters to score molecules synergizing with target action on indication
        # i.e. molecules that have effect opposite to desired drug affect on the targets
        if score_antagonists:
        # most common case when targets must be inhibited
            if direct_modulators:
                    effect = 'positive'
                    drug_class = 'direct agonists'
                    concept_ids = self.DirectAgonistsIDs
            else:
                    effect = 'positive'
                    drug_class = 'indirect agonists'
                    concept_ids = self.IndirectAgonistsIDs
        else:
        # case if drug are agonists
            if direct_modulators:
                    effect = 'negative'
                    drug_class = 'direct antagonists'             
                    concept_ids = self.DirectAntagonistsIDs
            else:
                    effect = 'negative'
                    drug_class = 'indirect antagonists'             
                    concept_ids = self.IndirectAntagonistsIDs

        return effect, drug_class, concept_ids


    def semantic_score(self,in_worksheet:str,target_effect_on_indication:str):
        """
        Input
        -----
        target_effect_on_indication: required Effect sign between target and Indication 
        target_effect_on_indication = 'positive' to score antagonists
        target_effect_on_indication = 'negative' to score agonists

        Output
        ------
        adds score columns to "in_worksheet" from raw_data
        """
        t_n = self._target_names()   
        colname = 'Activated by ' if target_effect_on_indication == 'positive' else 'Inhibited by '
        colname += t_n
        my_df = self.raw_data[in_worksheet]
        score4antagonists = True if target_effect_on_indication == 'positive' else False

        if self.__is_strict():
            booster_reltypes = ['Regulation','Biomarker','GeneticChange','QuantitativeChange','StateChange']
        else:
            booster_reltypes = ['Regulation','GeneticChange']
        self.set_how2connect(['Regulation'],[target_effect_on_indication],'',booster_reltypes)
        linked_row_count,linked_ent_ids,my_df = self.link2concept(colname,self.target_ids,my_df)
        print('%d indications are %sly regulated by %s' % (linked_row_count,target_effect_on_indication,t_n))

        self.score_GVs(my_df)

        link_effect, drug_class, concept_ids = self.__drug_connect_params(True,score4antagonists)
        if concept_ids:
            # references suggesting that known drugs for the target as treatments for indication
            colname = t_n+' '+drug_class+' clin. trials'
            self.set_how2connect(['ClinicalTrial'],[],'')
            linked_row_count,linked_ent_ids,my_df = self.link2concept(colname,concept_ids,my_df)
            print('Linked %d clinical trial indictions for %s %s' % (linked_row_count,t_n,drug_class))

            colname = t_n+' '+drug_class
            self.set_how2connect(['Regulation'],[link_effect],'',['Regulation'])
            linked_row_count,linked_ent_ids,my_df = self.link2concept(colname,concept_ids,my_df)
            print('Linked %d indications for %s %s' % (linked_row_count,t_n,drug_class))

        #references reporting target agonists exacerbating indication or causing indication as adverse events
        link_effect, drug_class, concept_ids = self.__drug_tox_params(True,score4antagonists)
        if concept_ids:
            colname = t_n+' '+drug_class
            self.set_how2connect(['Regulation'],[link_effect],'',['Regulation'])
            linked_row_count,linked_ent_ids,my_df = self.link2concept(colname,concept_ids,my_df)
            print('Linked %d indications as toxicities for %s %s' % (linked_row_count,t_n,drug_class))

    
        #references where target expression or activity changes in the indication
        colname = ' is upregulated' if target_effect_on_indication == 'positive' else ' is downregulated'
        colname = t_n + colname
        self.set_how2connect(['QuantitativeChange'],[target_effect_on_indication],'',['Biomarker','StateChange'])
        linked_row_count,linked_ent_ids,my_df= self.link2concept(colname,self.target_ids,my_df)
        print('%d indications %sly regulate %s' % (linked_row_count,target_effect_on_indication,t_n))


        #references suggesting target partners as targets for indication
        if hasattr(self, 'PartnerIndicationNetwork'):
            p_cl = self.partner_class
            colname = '{target_name} {partnet_class}s'.format(target_name=t_n,partnet_class=p_cl)
            self.set_how2connect(['Regulation'],[target_effect_on_indication],'',['Regulation'])
            linked_row_count,linked_ent_ids,my_df = self.link2concept(colname,self.partners_ids,my_df)
            print('Linked %d indications for %d %s %ss' % (linked_row_count,len(self.partners_ids),t_n,p_cl))

        #references reporting that cells producing the target linked to indication
        if hasattr(self, 'ProducingCellsIDs'):
            colname = '{target_name} producing cells'.format(target_name=t_n)
            self.set_how2connect(['Regulation'],[target_effect_on_indication],'',['Regulation'])
            linked_row_count,linked_ent_ids,my_df = self.link2concept(colname,self.ProducingCellsIDs,my_df)
            print('Liked %d indications linked %d cells producing %s' % (linked_row_count,len(self.ProducingCellsIDs),t_n))

        link_effect, drug_class, concept_ids = self.__drug_connect_params(False,score4antagonists)
        if concept_ids:
            # references suggesting that known drugs for the target as treatments for indication
            colname = t_n+' '+drug_class+' clin. trials'
            self.set_how2connect(['ClinicalTrial'],[],'')
            self._clone_(to_retrieve=REFERENCE_IDENTIFIERS)
            linked_row_count,linked_ent_ids,my_df = self.link2concept(colname,concept_ids,my_df)
            print('Linked %d clinical trial indications for %s %s' % (linked_row_count,t_n,drug_class))

            colname = t_n+' '+drug_class
            self.set_how2connect(['Regulation'],[link_effect],'',['Regulation'])
            self._clone_(to_retrieve=REFERENCE_IDENTIFIERS)
            linked_row_count,linked_ent_ids,my_df = self.link2concept(colname,concept_ids,my_df)
            print('Linked %d indications for %s %s' % (linked_row_count,t_n,drug_class))

        #references reporting target agonists exacerbating indication or causing indication as adverse events
        link_effect, drug_class, concept_ids = self.__drug_tox_params(False,score4antagonists)
        if concept_ids:
            colname = t_n+' '+drug_class
            self.set_how2connect(['Regulation'],[link_effect],'',['Regulation'])
            self._clone_(to_retrieve=REFERENCE_IDENTIFIERS)
            linked_row_count,linked_ent_ids,my_df = self.link2concept(colname,concept_ids,my_df)
            print('Linked %d indications as toxicities for %s %s' % (linked_row_count,t_n,drug_class))
        
        if hasattr(self, 'PathwayComponentsIDs'):
            #references linking target pathway to indication
            colname = t_n + ' pathway components'
            self.set_how2connect(['Regulation'],[target_effect_on_indication],'')
            self._clone_(to_retrieve=REFERENCE_IDENTIFIERS)
            linked_row_count,linked_ent_ids,my_df = self.link2concept(colname,list(self.PathwayComponentsIDs),my_df)
            print('Linked %d indications to %s pathway components' % (linked_row_count,t_n))

        counts_df = self.make_count_df(my_df,with_name=in_worksheet)
        self.add2raw(counts_df)
        return


    def other_effects(self):
        # need to be called after ranking to subtract self.all_entity_ids
        print('Findind indication linked with unknown effect to %s' % self._target_names())
        old_rel_props = self.relProps
        self.add_rel_props(PS_SENTENCE_PROPS+list(PS_BIBLIO_PROPS))
        t_n = self._target_names()
        ranked_indication_ids = self.Graph.get_node_ids(self.param['indication_types'])
        REQUEST_NAME = 'Find indications modulated by {target} with unknown effect'.format(target=t_n)
        oql2select_indications = 'SELECT Entity WHERE objectType = ({indication_type})'.format(indication_type=','.join(self.param['indication_types']))
        OQLquery = 'SELECT Relation WHERE objectType = (Regulation,QuantitativeChange) AND Effect = unknown AND NeighborOf({targets}) AND NeighborOf ({indications})'
        to_return = self.process_oql(OQLquery.format(targets=self.oql4targets, indications=oql2select_indications),REQUEST_NAME)
        
        REQUEST_NAME = 'Find indications where {target} was suggested as Biomarker'.format(target=t_n)
        OQLquery = 'SELECT Relation WHERE objectType = (Biomarker,StateChange,GeneticChange) AND NeighborOf({targets}) AND NeighborOf ({indications})'
        OQLquery = OQLquery.format(targets=self.oql4targets, indications=oql2select_indications)
        to_return.add_graph(self.process_oql(OQLquery,REQUEST_NAME))
        
        REQUEST_NAME = 'Find indications with genetically linked {target} Genetic Variants'.format(target=t_n)
        OQLquery = 'SELECT Relation WHERE objectType = FunctionalAssociation AND NeighborOf ({indications}) AND NeighborOf({GVs})'
        selectGVs = 'SELECT Entity WHERE id = ({ids})'
        OQLquery = OQLquery.format(indications=oql2select_indications, GVs=selectGVs)
        to_return.add_graph(self.iterate_oql(OQLquery,self.GVids,request_name=REQUEST_NAME))

        to_return.remove_nodes_from(ranked_indication_ids) # now delete all indication with known effect
        indication_ids = to_return.get_node_ids(self.param['indication_types'])
       # new_ids = set(indication_ids).difference(ranked_indication_ids)
       # new_ids = set(new_ids).difference(self.all_entity_ids)

        print('Found %d new indications linked to %s with unknown effect' % (len(indication_ids),self._target_names()))
        self.relProps = old_rel_props
        return to_return


    def load_indications4targets(self):
        '''
        assumes self.param['target_names'] is not empty
        '''
        start_time = time.time()
        self.set_targets()
        self.set_partners()

        self.find_target_indications()
        if self.target_class != 'Ligand':
            self.indications4chem_modulators()
        else:
            self.indications4cells_secreted_target()
        self.indications4partners()
        self.get_pathway_componets()
        self.__resolve_conflict_indications()

        print("Target indications were loaded in %s" % self.execution_time(start_time))

    def make_report(self):
        start_time = time.time()
        self.flush_dump_files()
        self.load_indications4targets()
        '''
        self.set_targets()
        self.set_partners()
        self.find_target_indications()
        self.get_pathway_componets()

        if self.target_class != 'Ligand':
        # if target is Ligand use indications4chem_modulators only if its antibody drugs have relations in Resnet
            self.indications4chem_modulators()
            #self.counterindications4chem_antimodulators()
        else:
            self.indications4cells_secreted_target() #will work only if target is Ligand

        self.indications4partners()
        self.__resolve_conflict_indications()
        '''
        
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
        print('Repurposing of %s was done in %s' % (self._get_report_name(), self.execution_time(start_time)))

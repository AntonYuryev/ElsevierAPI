from .ResnetGraph import ResnetGraph
from .SemanticSearch import SemanticSearch
from .SemanticSearch import PS_ID_TYPES,PS_SENTENCE_PROPS,COUNTS,BIBLIOGRAPHY,CHILDREN_COUNT,RANK
from .PathwayStudioGOQL import OQL
from ..ETM_API.etm import ETMstat
from ..pandas.panda_tricks import df
import time
import numpy as np

RANK_SUGGESTED_INDICATIONS = True 
# in RANK_SUGGESTED_INDICATIONS mode only indications suggested in the lietarure are ranked by amount of supporting evidence
PREDICT_RANK_INDICATIONS = False
# mode also predicts and ranks indications from diseases having input target as biomarker

class TargetIndications(SemanticSearch):
    pass
    def __init__(self, APIconfig, params={}):
        super().__init__(APIconfig)
        self.indication_types = ['Disease','Virus']
        self.pathway_name_must_include_target = True 
        # if True only pathways depicting target signaling are considered for regulome construction
        # if False all pathways containing both Targets and Partners will be considered
        self.target_activate_indication = True
        # use True when looking for Indications of Targets antagonists or for Toxicities of Targets agonists
        # use False when looking for Indications of Targets agonists or for Toxicities of Targets antagonists
        self.strict_mode = False
        # if strict mode algorithm only ranks indications suggetsed in the literarure without predicting new indications
        # if not strict mode algorithm also predicts and rank additional indications based on target expression or genetic profiles in disease
        
        if params: self.param = dict(params)
        else:
            self.param = {'partner_names':[],
    # if partner_names is empty script will try finding Ligands for Receptor targets and Receptors for Ligand targets
                 'partner_class':'', # use it only if partner_names not empty
                 'indication_types': [], #['Disease','Virus','CellProcess']
                 'target_names':[],
                 'target_type':'',
                 'to_inhibit':True,
                 'pathway_name_must_include_target':True,
                 'strict_mode':True,
                 'data_dir':'',
                 'add_bibliography' : True
                }

        self.child2parent = dict() # child2parent = {indication:parent_ontology_category}
        

    def _target_names(self): return ','.join([x['Name'][0] for x in self.Drug_Targets])


    def _partner_class(self):
        #finding effect for complementary receptor or ligands
        if self.target_class == 'Ligand': self.partner_class = "Receptor"
        elif self.target_class == 'Receptor': self.partner_class = "Ligand"
        else: self.target_class = NotImplemented


    def select_partners(self):
        self._partner_class()
        SELECTpartners = 'SELECT Entity WHERE Class = {partner_class} AND objectType = Protein AND Connected by (SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive) to ({select_target})'
        partners_graph = self.process_oql(SELECTpartners.format(partner_class=self.partner_class, select_target=self.find_targets_oql))
        if self.partner_class == 'Ligand':
            SELECTsecretedpartners = 'SELECT Entity WHERE "Cell Localization" = Secreted AND objectType = Protein AND Connected by (SELECT Relation WHERE objectType = DirectRegulation AND Effect = positive) to ({select_target})'
            secreted_partners = self.process_oql(SELECTsecretedpartners.format(select_target=self.find_targets_oql))
            partners_graph = partners_graph.add_graph(secreted_partners)

        self.partners_ids = partners_graph.get_node_ids(['Protein'])
        target_names = self._target_names()
        partners_str = ','.join(n['Name'][0] for n in partners_graph._get_nodes())
        print ('Found %s %ss as partners for %s'%(partners_str, self.partner_class,target_names))
        return self.Graph._get_nodes(self.partners_ids)


    def set_partners(self, partner_names:list, partner_class:str):
        # use it when receptor partners are metabolites and therfore cannot be found by select_partners()
        partners_str = OQL.join_with_quotes(partner_names)
        SELECTpartners = 'SELECT Entity WHERE Name = ({partner_names})'.format(partner_names=partners_str)
        partners_graph = self.process_oql(SELECTpartners)
        self.partners = partners_graph._get_nodes()
        self.partner_class = partner_class
        self.partners_ids = [p['Id'][0] for p in self.partners]
        self.find_partners_oql = 'SELECT Entity WHERE id = ('+ ','.join(map(str,self.partners_ids))+')'

        
    def set_targets(self):
        target_names =  list(self.param['target_names'])
        try:
            target_objtype = str(self.param['target_objtype'])
        except KeyError: target_objtype = 'Protein'
        partner_names=list(self.param['partner_names'])
        partner_class=str(self.param['partner_class'])
        to_inhibit=bool(self.param['to_inhibit'])
        strict_mode=bool(self.param['strict_mode'])
        indication_types=list(self.param['indication_types'])

        self.strict_mode = strict_mode
        if isinstance(to_inhibit, bool):
            self.target_activate_indication = to_inhibit
            # defaults to False if nothing is specified
        if isinstance(indication_types,list):
            self.indication_types = indication_types
        # defaults to [Disease,Virus] if nothing is specified

        prop_names_str, prop_values_str = OQL.get_search_strings(['Name'],target_names)
        self.add_ent_props(['Class'])
        self.find_targets_oql = 'SELECT Entity WHERE ({prop_name}) = ({values}) AND objectType = '+target_objtype 
        targets_graph = self.process_oql(self.find_targets_oql.format(prop_name=prop_names_str,values=prop_values_str))
        
        self.Drug_Targets = targets_graph._get_nodes()
        self.target_ids = [x['Id'][0] for x in self.Drug_Targets]

        if not hasattr(self,'input_names'): #input names are used for finding references in ETM.
        # RepurposeDrug has it own 'input_names'
            self.input_names = [x['Name'][0] for x in self.Drug_Targets]

        self.find_targets_oql = 'SELECT Entity WHERE id = ('+ ','.join(map(str,self.target_ids))+')'
        for t in self.Drug_Targets:
            try:
                self.target_class = t['Class'][0]
                break
            except KeyError: continue
        
        if partner_names:
            self.set_partners(partner_names,partner_class)
        else:
            self.partners = self.select_partners()
            if self.partners:
                self.partners_ids = [p['Id'][0] for p in self.partners]
                self.find_partners_oql = 'SELECT Entity WHERE id = ('+ ','.join(map(str,self.partners_ids))+')'
            else: self.find_partners_oql = ''
        

    def _get_report_name(self):
        indics = ','.join(self.indication_types)
        mode = ' antagonists' if self.target_activate_indication else ' agonists'
        rep_pred = 'suggested ' if self.strict_mode else 'suggested,predicted ' 
        return self.param['data_dir']+rep_pred+ indics+' for '+ self._target_names() + mode


    def GVindications(self):
        target_names = [x['Name'][0] for x in self.Drug_Targets]
        t_n = ','.join(target_names)
        select_targets = 'SELECT Entity WHERE id = ({target_ids})'.format(target_ids = ','.join(map(str,self.target_ids)))
  
        selectGVs = 'SELECT Entity WHERE objectType = GeneticVariant AND Connected by (SELECT Relation WHERE objectType = GeneticChange) to ({select_target})'
        selectGVs = selectGVs.format(select_target=select_targets)
        REQUEST_NAME = 'Find indications linked to {target} genetic variants'.format(target=t_n)
        OQLquery = 'SELECT Relation WHERE NeighborOf({select_gvs}) AND NeighborOf (SELECT Entity WHERE objectType = ({indication_type}))'
        self.GVsInDiseaseNetwork = self.process_oql(OQLquery.format(select_gvs=selectGVs,indication_type=','.join(self.indication_types)),REQUEST_NAME)
        found_indications = self.GVsInDiseaseNetwork.get_node_ids(self.indication_types)
        self.GVids = self.GVsInDiseaseNetwork.get_node_ids(['GeneticVariant'])
        
        print('Found %d indications genetically linked to %d Genetic Variants in %s' % 
            (len(found_indications), len(self.GVids), t_n))
        
        self.target_indications4strictmode.update(found_indications)
        return found_indications


    def set_strictmode_indications(self):
        self.indications4strictmode = self.target_indications4strictmode


    def find_target_indications(self):
        f_t = self.find_targets_oql
        t_n = self._target_names()
        indications2return = set()
        effect = 'positive' if self.target_activate_indication else 'negative'
        REQUEST_NAME = 'Find indications {effect}ly modulated by {target}'.format(effect=effect, target=t_n)
        select_indications = 'SELECT Entity WHERE objectType = ({indication_type})'.format(indication_type=','.join(self.indication_types))
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = {effect} AND NeighborOf({select_target}) AND NeighborOf ({indications})'      
        ModulatedByTargetNetwork = self.process_oql(OQLquery.format(select_target=f_t, effect = effect, indications=select_indications),REQUEST_NAME)
        found_indication_ids = ModulatedByTargetNetwork.get_node_ids(self.indication_types)
        self.target_indications4strictmode = set(found_indication_ids)
        indications2return.update(found_indication_ids)
        
        print('Found %d diseases %sly regulated by target %s' % (len(found_indication_ids),effect,t_n))
        
        REQUEST_NAME = 'Find indications {effect}ly regulating {target}'.format(effect=effect, target=t_n)
        OQLquery = 'SELECT Relation WHERE objectType = QuantitativeChange AND Effect = {effect} AND NeighborOf ({select_target}) AND NeighborOf ({indications})'  
        ActivatedInDiseaseNetwork = self.process_oql(OQLquery.format(select_target=f_t,effect = effect,indications=select_indications),REQUEST_NAME)
        found_indications= ActivatedInDiseaseNetwork.get_node_ids(self.indication_types)
        indications2return.update(found_indications)
        print('Found %d diseases where target %s is %sly regulated' % (len(found_indications),t_n,effect))
        
        REQUEST_NAME = 'Find indications where {target} is biomarker'.format(target=t_n)
        OQLquery = 'SELECT Relation WHERE objectType = Biomarker AND NeighborOf({select_target}) AND NeighborOf ({indications})'
        BiomarkerInDiseaseNetwork = self.process_oql(OQLquery.format(select_target=f_t,indications=select_indications),REQUEST_NAME)
        found_indications = BiomarkerInDiseaseNetwork.get_node_ids(self.indication_types)
        indications2return.update(found_indications)
        print('Found %d indications where target %s is claimed as biomarker' %  (len(found_indications),t_n))

        indications2return.update(self.GVindications())

        if self.strict_mode:
            self.set_strictmode_indications()
            self.select_indications = 'SELECT Entity WHERE id = ({ids})'.format(ids=','.join(map(str,self.indications4strictmode)))
        else:
            self.select_indications = 'SELECT Entity WHERE objectType = ({indication_type})'.format(indication_type=','.join(self.indication_types))

        return list(indications2return)
       

    def indications4chem_modulators(self, linked_to_target_by=['DirectRegulation']):
        f_t = self.find_targets_oql
        t_n = self._target_names()
        effect = 'negative'if self.target_activate_indication else 'positive'
        
        REQUEST_NAME = 'Find substances {effect}ly regulating {target}'.format(effect=effect,target=t_n)
        OQLquery = 'SELECT Relation WHERE objectType = {rel_type} AND Effect = {effect} AND NeighborOf ({select_target}) AND NeighborOf (SELECT Entity WHERE objectType = SmallMol AND Connectivity > 1)'   
        TargetInhibitorsNetwork = self.process_oql(OQLquery.format(rel_type=linked_to_target_by, select_target=f_t, effect=effect),REQUEST_NAME)
        self.TargetInhibitorsIDs = TargetInhibitorsNetwork.get_node_ids(['SmallMol'])
        print('Found %d substances %sly regulating %s' % (len(self.TargetInhibitorsIDs),effect,t_n))
        
        if self.strict_mode: return # no further predictions
        REQUEST_NAME = 'Find indications for substances {effect}ly regulating {target}'.format(effect=effect,target=t_n)
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = {effect}  AND NeighborOf (SELECT Entity WHERE objectType = ({indication_type})) AND NeighborOf (SELECT Entity WHERE id = ({drugIDlist}))'
        OQLquery = OQLquery.format(drugIDlist=','.join(map(str, self.TargetInhibitorsIDs)), effect=effect,indication_type=','.join(self.indication_types))
        InhibitorsIndicationNetwork = self.process_oql(OQLquery,REQUEST_NAME)

        found_indications= InhibitorsIndicationNetwork.get_node_ids(self.indication_types)
        print('Found %d indications for %d substances %sly regulating %s' %  
             (len(found_indications), len(self.TargetInhibitorsIDs),effect,t_n))


    def counterindications4chem_antimodulators(self,liked_to_target_by=['DirectRegulation']):
        # liked_to_target_by can be relaxed to (DirectRegulation,Regulation,Expression,MolTransport)
        f_t = self.find_targets_oql
        t_n = self._target_names()
        opposite_effect = 'positive' if self.target_activate_indication else 'negative'
        REQUEST_NAME = 'Find substances {effect}ly regulating {target}'.format(effect=opposite_effect,target=t_n)
        OQLquery = 'SELECT Relation WHERE objectType = {rel_type} and Effect = {effect} AND NeighborOf ({select_target}) AND NeighborOf (SELECT Entity WHERE objectType = SmallMol AND Connectivity > 1)'     
        TargetAgonistNetwork = self.process_oql(OQLquery.format(rel_type=liked_to_target_by,select_target=f_t, effect=opposite_effect),REQUEST_NAME)
        self.TargetAgonistIDs = list(set([x for x,y in TargetAgonistNetwork.nodes(data=True) if y['ObjTypeName'][0] in ['SmallMol']]))
        print('Found %d substances %sly regulating %s' % (len(self.TargetAgonistIDs),opposite_effect,t_n))
        
        if self.strict_mode: return # no further predictions
        REQUEST_NAME = 'Find counterindications for substances {effect}ly regulating {target}'.format(effect=opposite_effect,target=t_n)
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = {effect} AND NeighborOf (SELECT Entity WHERE objectType = ({indication_type})) AND NeighborOf (SELECT Entity WHERE id = ({drugIDlist}))'
        OQLquery = OQLquery.format(drugIDlist=','.join(map(str, self.TargetAgonistIDs)), effect=opposite_effect,indication_type=','.join(self.indication_types))
        AgonistToxicitiesnNetwork = self.process_oql(OQLquery,REQUEST_NAME)
        found_indications= AgonistToxicitiesnNetwork.get_node_ids(self.indication_types)
        print('Found %d indications for %d substances %sly regulating %s' %  
             (len(found_indications), len(self.TargetAgonistIDs),opposite_effect,t_n))

     
    def indications4partners(self):
        t_n = self._target_names()
        effect = 'positive' if self.target_activate_indication else 'negative'
        REQUEST_NAME = 'Find indications for {partner}s of {targets} '.format(targets=t_n,partner=self.partner_class.lower())
        OQLquery = 'SELECT Relation WHERE objectType = (Regulation,QuantitativeChange) AND Effect = {effect} AND NeighborOf ({select_partners}) AND NeighborOf ({select_indications})'
        if self.partners:
            OQLquery = OQLquery.format(effect=effect,select_partners=self.find_partners_oql,select_indications=self.select_indications)
            self.PartnerIndicationNetwork = self.process_oql(OQLquery,REQUEST_NAME)
            found_indications = self.PartnerIndicationNetwork.get_node_ids(self.indication_types)
            print('Found %d indications for %d %s %ss' %  
                 (len(found_indications), len(self.partners),t_n,self.partner_class.lower()))


    def indications4cells_secreting_target(self):
        t_n = self._target_names()
        if self.target_class == 'Ligand':
            REQUEST_NAME = 'Find indications linked to cells secreting the {target}'.format(target=t_n)
            OQLquery = 'SELECT Relation WHERE objectType = (CellExpression,MolTransport) AND NeighborOf ({select_targets}) AND NeighborOf (SELECT Entity WHERE objectType = Cell)'
            cells_make_target = self.process_oql(OQLquery.format(select_targets=self.find_targets_oql),REQUEST_NAME)
            self.ProducingCellsIDs = cells_make_target.get_node_ids(['CellType'])
            print('Found %d cell types producing %s' % (len(self.ProducingCellsIDs),t_n))

            if self.strict_mode: return # no further predictions
            REQUEST_NAME = 'Find indications linked to cell secrteing {target}'
            REQUEST_NAME = REQUEST_NAME.format(target = t_n)
            effect = 'positive' if self.target_activate_indication else 'negative'
            OQLquery = 'SELECT Relation WHERE Effect = {effect} AND NeighborOf (SELECT Entity WHERE objectType = ({indication_type})) AND NeighborOf (SELECT Entity WHERE id = ({cell_ids}))'
            OQLquery = OQLquery.format(effect=effect,cell_ids=','.join(map(str, self.ProducingCellsIDs)),indication_type=','.join(self.indication_types))
            CellDiseaseNetwork = self.process_oql(OQLquery,REQUEST_NAME)

            found_indications= CellDiseaseNetwork.get_node_ids(self.indication_types)
            print('Found %d indications for %d cells producing %s' %  
                 (len(found_indications), len(self.ProducingCellsIDs),t_n))


    def pathway_oql(self):
        #set pathway_name_must_include_target to True if targets have a lot of large curated pathways
        target_oqls = dict()
        pct = '%'
        merged_pathways = 'SELECT Relation WHERE objectType = (DirectRegulation,Binding,ProtModification,PromoterBinding,ChemicalReaction) AND MemberOf ({select_networks})'

        for target in self.Drug_Targets:
            target_id = target['Id'][0]
            target_name = target['Name'][0]
            SELECTpathways = 'SELECT Network WHERE ParentOf (SELECT Entity WHERE id = ({i}))'.format(i=str(target_id))
            if self.pathway_name_must_include_target:
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
        oql_queries = self.pathway_oql() # separate oql_qury for each target
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
        else: return ResnetGraph()


    def init_semantic_search(self):
        print('\n\nInitializing semantic search')
        t_n = self._target_names()
        if self.strict_mode:
            IndicationNames = [y['Name'][0] for x,y in self.Graph.nodes(data=True) if y['Id'][0] in self.indications4strictmode]
            if not IndicationNames:
                print ('No indications found for %s in strict mode' % t_n)
                return False
        else:
            IndicationNames = [y['Name'][0] for x,y in self.Graph.nodes(data=True) if y['ObjTypeName'][0] in self.indication_types]
        
        if IndicationNames:
            indication2score = df()
            indication2score['Name'] = np.array(IndicationNames)
            print('Will score %d indications linked to %s' % (len(indication2score),t_n))
            self.load_pandas(indication2score,prop_names_in_header=True,map2type=self.indication_types)
            self.add_rel_props(PS_ID_TYPES)
            return True
        else:
            return False
        

    def score_GVs(self):
        if not self.GVids: return
        t_n = self._target_names()
        self.__colnameGV__ = t_n+' GVs'
        gvlinkcounter = 0
        for i in self.RefCountPandas.index:
            row_entity_ids = list(self.RefCountPandas.at[i,self.__temp_id_col__])
            if hasattr(self,'GVsInDiseaseNetwork'):
                GVneighborsId = self.GVsInDiseaseNetwork.get_neighbors(row_entity_ids, self.GVids)
                
            GVscore = 0
            if len(GVneighborsId) > 0:
                GVnames = set([n[0] for i,n in self.Graph.nodes.data('Name') if i in GVneighborsId])
                self.RefCountPandas.at[i,self.__colnameGV__] = ';'.join(GVnames)
                gvlinkcounter += 1
                gv_disease_subgraph = self.GVsInDiseaseNetwork.get_subgraph(row_entity_ids,GVneighborsId)
                GVscore = len(gv_disease_subgraph.load_references())
            
            self.RefCountPandas.at[i,self._col_name_prefix+'GVs'] = GVscore

        print('Found %d indications linked to %d GVs' % (gvlinkcounter, len(self.GVids)))


    def __drug_connect_params(self):
        # function returns how2connect paramters to score molecules similar to desired drug affect on the targets 
        effect = None
        drug_class = None
        concept_ids = None

        if self.target_activate_indication:
        # most common case when targets must be inhibited
            if hasattr(self, 'TargetInhibitorsIDs'):
                effect = 'negative'
                drug_class = 'antagonists'
                concept_ids = self.TargetInhibitorsIDs
        else:
        # case if drug are agonists
            if hasattr(self, 'TargetAgonistIDs'):
                effect = 'positive'
                drug_class = 'agonists'
                concept_ids = self.TargetAgonistIDs

        return effect, drug_class, concept_ids

    def __drug_tox_params(self):
        effect = None
        drug_class = None
        concept_ids = None
        # function returns how2connect parameters to score molecules synergizing with target action on indication
        # i.e. molecules that have effect opposite to desired drug affect on the targets
        if self.target_activate_indication:
        # most common case when targets must be inhibited
            if hasattr(self, 'TargetAgonistIDs'):
                effect = 'positive'
                drug_class = 'agonists'
                concept_ids = self.TargetAgonistIDs
        else:
        # case if drug are agonists
            if hasattr(self, 'TargetInhibitorsIDs'):
                effect = 'negative'
                drug_class = 'antagonists'
                concept_ids = self.TargetInhibitorsIDs

        return effect, drug_class, concept_ids


    def score_semantics(self):
        t_n = self._target_names()
        #references supporting target modulation of indication
        if self.target_activate_indication:
            effect = 'positive'
            colname = 'Activated by '+ t_n
        else:
            effect = 'negative'
            colname = 'Inhibited by '+ t_n

        self.set_how2connect(['Regulation'],[effect],'')
        linked_entities_count = self.link2concept(colname,self.target_ids)
        print('%d indications are %sly regulated by %s' % (linked_entities_count,effect,t_n))
    
        #references where target expression or activity changes in the indication
        if self.target_activate_indication:
            effect = 'positive' 
            colname = t_n + ' is upregulated'
        else:
            effect = 'negative'
            colname = t_n + ' is downregulated'

        self.set_how2connect(['QuantitativeChange'],[effect],'')
        linked_entities_count= self.link2concept(colname,self.target_ids)
        print('%d indications %sly regulate %s' % (linked_entities_count,effect,t_n))

        #references suggested target as a biomarker for indication
        colname = t_n + ' is Biomarker'
        self.set_how2connect(['Biomarker'],[],'')
        linked_entities_count = self.link2concept(colname,self.target_ids)
        print('Linked %d indications where %s is biomarker' % (linked_entities_count,t_n))

        self.score_GVs()

        effect, drug_class, concept_ids = self.__drug_connect_params()
        if isinstance(effect,str):
            # references suggesting that existing drugs for the target as treatments for indication
            colname = t_n+' '+drug_class
            self.set_how2connect(['Regulation'],[effect],'')
            linked_entities_count = self.link2concept(colname,concept_ids)
            print('Linked %d indications for %s %s' % (linked_entities_count,t_n,drug_class))

        #references suggesting target partners as targets for indication
        if hasattr(self, 'PartnerIndicationNetwork'):
            p_cl = self.partner_class
            colname = '{target_name} {partnet_class}s'.format(target_name=t_n,partnet_class=p_cl)         
            effect = 'positive' if self.target_activate_indication else 'negative'
            self.set_how2connect(['Regulation'],[effect],'')
            linked_entities_count = self.link2concept(colname,self.partners_ids)
            print('Linked %d indications for %d %s %ss' % (linked_entities_count,len(self.partners_ids),t_n,p_cl))

        #references reporting that cells producing the target linked to indication
        if hasattr(self, 'ProducingCellsIDs'):
            colname = '{target_name} producing cells'.format(target_name=t_n)
            effect = 'positive' if self.target_activate_indication else 'negative'
            self.set_how2connect(['Regulation'],[effect],'')

            linked_entities_count = self.link2concept(colname,self.ProducingCellsIDs)
            print('Liked %d indications linked %d cells producing %s' % (linked_entities_count,len(self.ProducingCellsIDs),t_n))
        
        #references reporting target agonists exacerbating indication or causing indication as adverse events
        effect, drug_class, concept_ids = self.__drug_tox_params()
        if isinstance(effect, str):
            colname = t_n+' '+drug_class
            self.set_how2connect(['Regulation'],[effect],'')
            linked_entities_count = self.link2concept(colname,concept_ids)
            print('Linked %d indications are toxicities for %s %s' % (linked_entities_count,t_n,drug_class))
        
        if hasattr(self, 'PathwayComponentsIDs'):
            #references linking target pathway to indication
            colname = t_n + ' pathway components'
            effect = 'positive' if self.target_activate_indication else 'negative'
            self.set_how2connect(['Regulation'],[effect],'')

            linked_entities_count = self.link2concept(colname,list(self.PathwayComponentsIDs))
            print('Linked %d indications to %s pathway components' % (linked_entities_count,t_n))

        counts_pd = self.make_count_pd()
        counts_pd.name = COUNTS
        self.add2raw(counts_pd)
        return

    def normalize_counts(self,bibliography=True):
        normalized_count_pd = df()
        weights = df()
        refcount_cols = [col for col in self.RefCountPandas.columns if self._col_name_prefix in col]
        number_of_weights = len(refcount_cols)
        
        normalized_count_pd['Name'] = self.RefCountPandas['Name']
        weights.at[0,'Name'] = 'WEIGHTS:'
        weight_index = 0
        for col in refcount_cols:
            col_max = self.RefCountPandas[col].max()
            normalized_count_pd[col] = self.RefCountPandas[col]/col_max if col_max > 0 else 0.0
            column_weight = (number_of_weights-weight_index)/number_of_weights
            weights.at[0,col] = column_weight
            weight_index += 1
     
        #calculating cumulative score  
        combined_scores = list()
        for i in normalized_count_pd.index:
            scores_row = normalized_count_pd.loc[[i]]
            weighted_sum = 0.0
            for col in refcount_cols:
                weighted_sum = weighted_sum + scores_row[col]*weights.at[0,col]
            combined_scores.append(weighted_sum)

        normalized_count_pd['Combined score'] = np.array(combined_scores)
        if hasattr(self,'__colnameGV__'):
            try:
                normalized_count_pd[self.__colnameGV__] = self.RefCountPandas[self.__colnameGV__]
            except KeyError:
                pass

        normalized_count_pd[CHILDREN_COUNT] = self.RefCountPandas[self.__temp_id_col__].apply(lambda x: len(x))
        normalized_count_pd[RANK] = normalized_count_pd['Combined score']/normalized_count_pd[CHILDREN_COUNT]
        normalized_count_pd = normalized_count_pd.sort_values(by=[RANK],ascending=False)
        normalized_count_pd = normalized_count_pd.loc[normalized_count_pd[RANK] > 0.0] # removes rows with all zeros
        normalized_count_pd['Combined score'] = normalized_count_pd['Combined score'].map(lambda x: '%2.3f' % x)
        normalized_count_pd[RANK] = normalized_count_pd[RANK].map(lambda x: '%2.3f' % x)

        if self.child2parent:
            def get_ontology_parent(indication_name):
                try:
                    parents = self.child2parent[indication_name]
                    return ';'.join(parents)
                except KeyError: return ''
            normalized_count_pd['Ontology group'] = normalized_count_pd['Name'].apply(get_ontology_parent)


        #fprefix = self.fname_prefix()
        if bibliography:
            print('Finding relevant articles for each indication in ETM' )
            max_ref_count = 5
            etm_counter = ETMstat(self.APIconfig,limit=max_ref_count)
            def add_refs(indication:str):
                refid_list = list() # list will be sorted by relevance
                for drug_or_target in self.input_names:
                    search_terms = [drug_or_target,indication]
                    hit_count,ref_ids,etm_refs = etm_counter.relevant_articles(search_terms)
                    for id_type,identifiers in ref_ids.items():
                        for identifier in identifiers:
                            if identifier not in refid_list:
                                refid_list.append(identifier)

                    if len(refid_list) >= max_ref_count: 
                        refid_list = refid_list[:max_ref_count]
                        break
                return ';'.join(refid_list)

            normalized_count_pd['References'] = normalized_count_pd['Name'].apply(add_refs)
            biblio_pd = etm_counter.counter2pd()
            

      # use this filter only if indications linked to GVs are of interest
      #  try: normalized_count_pd = normalized_count_pd.loc[normalized_count_pd[self.__colnameGV__].notna()]
      #  except KeyError: pass

        counts_pd = self.raw_data[COUNTS]
        for i in range(1,len(weights.columns)):
            col = weights.columns[i]
            weight = weights.loc[0,col]
            if weight > 0.0:
                weights.loc[0,col] = '{:,.4f}'.format(weight)

        normalized_count_pd = weights.append(normalized_count_pd,ignore_index=True)
        weigths_header = list(normalized_count_pd.iloc[0].replace(np.nan, ''))
        normalized_count_pd.loc[0] = weigths_header
        normalized_count_pd.name = 'normalized'
        self.add2raw(normalized_count_pd)

        ranked_counts_pd = self.merge_counts2norm(counts_pd,normalized_count_pd)
        ranked_counts_pd.loc[0] = weigths_header
        # merge_counts2norm blanks out weigths_header because counts_pd does not have it

        ranked_counts_pd.name = 'ranked_counts'
        self.add2report(ranked_counts_pd)

        if bibliography:
            biblio_pd.name = BIBLIOGRAPHY
            self.add2report(biblio_pd)
        

    def load_ontology(self, parent_ontology_groups:list):
        """
        loads self.child2parent = {child_name:[ontology_parent_names]}
        """
        self.add_ent_props(['Name'])
        for group_name in parent_ontology_groups:
            childs = self.child_graph([group_name],['Name'],include_parents=False)
            for id, child in self.Graph.nodes(data=True):
                child_name = child['Name'][0]
                try:
                    self.child2parent[child_name].append(group_name)
                except KeyError:
                    self.child2parent[child_name] = [group_name]


    def other_effects(self):
        old_rel_props = self.relProps
        self.add_rel_props(PS_SENTENCE_PROPS+['Title','PubYear'])
        t_n = self._target_names()
        REQUEST_NAME = 'Find indications modulated by {target} with unknown effect'.format(target=t_n)
        select_indications = 'SELECT Entity WHERE objectType = ({indication_type})'.format(indication_type=','.join(self.indication_types))
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = unknown AND NeighborOf({select_target}) AND NeighborOf ({indications})'
        to_return = self.process_oql(OQLquery.format(select_target=self.find_targets_oql, indications=select_indications),REQUEST_NAME)
        self.relProps = old_rel_props
        return to_return


    def rn2pd(self, ModulatedNetwork:ResnetGraph, by_entity:str):
        old_rel_props = self.relProps
        self.add_rel_props(PS_SENTENCE_PROPS+['Title','PubYear'])
        ref_pd = df(columns=['Indication','Entity','#Reference','PMID','DOI','Title','PubYear','Snippet'])
        for n1,n2,rel in ModulatedNetwork.edges.data('relation'):
            node1 = ModulatedNetwork._get_node(n1)
            node2 = ModulatedNetwork._get_node(n2)
            indication = node1 if node1['ObjTypeName'][0] in self.indication_types else node2

            ref_counter = set()
            links_between = self.connect_nodes([n1],[n2])
            references = links_between.load_references()
            ETMstat.count_refs(ref_counter,references)
            #refs come from Resnet without Relevance 
        
            for ref in list(ref_counter):
                ref_list = ref.to_list(id_types=['PMID','DOI'],
                                print_snippets=True,biblio_props=['Title','PubYear']) #other_props=['Citation index']
                row = [indication['Name'][0],by_entity,len(references)]+ref_list
                ref_pd.loc[len(ref_pd.index)] = row
        
        ref_pd.sort_values(by=['#Reference','Indication','PMID','DOI'],ascending=False, inplace=True)
        ref_pd.name = 'possibilities'
        
        self.relProps = old_rel_props
        return ref_pd


    def _worksheet_prefix(self):
        indications = ','.join(self.indication_types)
        act = 'Act.' if self.target_activate_indication else 'Inh.'
        return act+indications+'.'


    def make_report(self):
        start_time = time.time()
        self.data_dir = self.param ['data_dir']
        
        self.set_targets()
        # strict mode ranks only indications suggetsed in the literarure without predicting new indications
        self.pathway_name_must_include_target = self.param['pathway_name_must_include_target']
        self.flush_dump_files()

        self.find_target_indications()
        self.get_pathway_componets()

        if self.target_class != 'Ligand':
        # if target is Ligand use indications4chem_modulators only if its antibody drugs have relations in Resnet
            self.indications4chem_modulators()
            self.counterindications4chem_antimodulators()

        self.indications4partners()
        self.indications4cells_secreting_target() #will work only if target is Ligand
        
        self.init_semantic_search()
        self.score_semantics()

        self.normalize_counts(bibliography=self.param['add_bibliography'])
        print('Repurposing of %s was done in %s' % (self.fname_prefix(), self.execution_time(start_time)))

import time,os,math
from .PathwayStudioGOQL import OQL
from ..ETM_API.references import Reference
from .ResnetGraph import ResnetGraph,PSObject,EFFECT,df
from .ResnetAPISession import APISession,len
from .ResnetAPISession import DO_NOT_CLONE,BELONGS2GROUPS,NO_REL_PROPERTIES,REFERENCE_IDENTIFIERS
from ..ETM_API.etm import ETMstat
from ..utils import execution_time
#from ..ETM_API.scibite import SBSstat
import pandas as pd
import numpy as np

COUNTS = 'counts'
PS_BIBLIOGRAPHY = 'EBKGrefs'
ETM_BIBLIOGRAPHY = 'ETMrefs'
CHILDREN_COUNT = 'Number of ontology children'
RANK = 'PRTS rank' # PRTS = Probability of Technical and Regulatory Success
ONTOLOGY_ANALYSIS = 'ontology'
WEIGHTED = 'weighted '

class SemanticSearch (APISession):
    __refCache__='reference_cache.tsv'
    __cntCache__='counts_cache.tsv'
    __concept_name__ = 'Concept name'
    __mapped_by__ = 'mapped_by'
    __resnet_name__ = 'Resnet name'
    _col_name_prefix = WEIGHTED
    _col_name_prefix2 = "RefCount to "
    __temp_id_col__ = 'entity_IDs'
    __rel_dir__ = '' # relation directions: allowed values: '>', '<',''
    __iter_size__ = 1000 #controls the iteration size in sematic reference count
    __print_refs__ = True
    data_dir = ''
    iteration_step = 1000
    min_etm_relevance = 0.0
    __colname4GV__  = ''

    
    def __init__(self,*args,**kwargs):
        '''
        input:
            APIconfig - args[0]
            kwargs:
            what2retrieve - defaults NO_REL_PROPERTIES
            what2retrieve other options:
            [DATABASE_REFCOUNT_ONLY,REFERENCE_IDENTIFIERS,BIBLIO_PROPERTIES,SNIPPET_PROPERTIES,ONLY_REL_PROPERTIES,ALL_PROPERTIES]
            connect2server - default True, set to False to run script using data in __pscache__ files instead of database
            max_ontology_parent - default: 10
        '''
        my_kwargs = {
                'what2retrieve':REFERENCE_IDENTIFIERS,
                'max_ontology_parent': 10
                }
        my_kwargs.update(kwargs)
        session_kwargs, parameters = self.split_kwargs(my_kwargs)
        super().__init__(*args,**session_kwargs)

        self.params = dict(parameters)

        self.PageSize = 1000
        self.RefCountPandas = df() # stores result of semantic retreival
        self.RefCountPandas._name_ = COUNTS

        self.report_pandas=dict() # stores pandas used to generate report file
        self.raw_data = dict() # stores raw data pandas used for generation pandas in self.report_pandas
        self.__only_map2types__ = list()
        self.__connect_by_rels__= list()
        self.__rel_effect__ = list()
        self.references = dict() # {identifier:Reference} made by self.Graph.citation_index()
        self.etm_counter = ETMstat(self.APIconfig,limit=10)
        #self.etm_counter = SBSstat(self.APIconfig,limit=10)
        self.etm_counter.min_relevance = self.min_etm_relevance
        self.all_entity_dbids = set()
        self.columns2drop = [self.__temp_id_col__] # columns to drop before printing pandas
        self.__boost_with__ = list()
        self.relweight_prop = str()
        self.relweight_dict = dict() # = dict[str,float]
        self.max_ontology_parent = self.params.get('max_ontology_parent',10)
        self.nodeweight_prop = str()
        self._Ontology = list()


    def reset(self):
        score_columns = [col for col in self.RefCountPandas.columns if self._col_name_prefix in col]
        for col in score_columns:
                self.RefCountPandas[col] = 0.0
        print('DataFrame was reset')


    def refcount_columns(self,counts_df=df(),column_prefix=''):
        d = self.RefCountPandas if counts_df.empty else counts_df
        refcount_column_prefix = column_prefix if column_prefix else self._col_name_prefix
        to_return = [col for col in d.columns if refcount_column_prefix in col]
        return to_return


    def _clone(self, **kwargs):
        '''
        Return
        ------
        Copy of SemanticSearch object with copy of self.Graph
        '''
        my_kwargs = dict(kwargs)
        my_kwargs['copy_graph'] = True # need to copy self.Graph to enable using cache in link2concept
        api_session = super()._clone_session(**my_kwargs) 
        new_session = SemanticSearch(api_session.APIconfig,**my_kwargs)
        new_session.entProps = api_session.entProps
        new_session.relProps = api_session.relProps
        new_session.data_dir = self.data_dir
        new_session.__connect_by_rels__ = self.__connect_by_rels__
        new_session.__rel_effect__ = self.__rel_effect__
        new_session.__rel_dir__ = self.__rel_dir__
        new_session.__boost_with__ = self.__boost_with__
        new_session.relweight_prop = self.relweight_prop
        new_session.relweight_dict = self.relweight_dict
        new_session.iteration_step = self.iteration_step
        if self.__rel_effect__:
            new_session.add_rel_props([EFFECT])
        
        new_session.add2self = False 
        # cloning is done to avoid adding irrelevant references to self.Graph to avoid in PS_Bibliography worksheet
        # therefore new_session.add2self is set to false
        return new_session
    

    def __concept(self,colname:str):
        return colname[len(self._col_name_prefix+self._col_name_prefix2):]
    
    def _refcount_colname(self,concept_name:str):
        return self._col_name_prefix2 + concept_name
    
    def _weighted_refcount_colname(self,concept_name:str):
        return self._col_name_prefix+self._col_name_prefix2 + concept_name
    
    def _linkedconcepts_colname(self,concept_name:str):
        return 'Linked '+ concept_name + ' count'
    
    def _concept_size_colname(self,concept_name:str):
        return concept_name + ' count'


    def drop_refcount_columns(self,counts_df=df.from_pd(pd.DataFrame())):
        """
        Removes
        -----
        columns with refcount prefix before redoing semantic search with weighted reference count
        """
        d = self.RefCountPandas if counts_df.empty else counts_df
        ref_count_columns = self.refcount_columns(d)
        d.drop(columns=ref_count_columns, inplace=True) #defaults to deleting rows, column must be specified
        print('%d refcount columns were dropped' % len(ref_count_columns))


    def report_name(self):
        return 'semantic refcount'


    def report_path(self,extension=''):
        if extension:
            ext = extension if extension.find('.')>-1 else '.'+extension
        else:
            ext = ''
        return self.data_dir+self.report_name()+ext


    def _all_dbids(self,df2score=None):
        """
        Returns
        -------
        list of all ids from the tuples in column self.__temp_id_col__].\n
        if df2score is empty returns all ids from self.RefCountPandas
        """
        if df2score is None:
            return self.all_entity_dbids
        else:
            id_set = set()
            try:
                list_of_tuples = list(df2score[self.__temp_id_col__])          
                [id_set.update(list(lst)) for lst in list_of_tuples]
            except KeyError:
                pass
            return id_set


    def entities(self,from_df=df(),id_column='Name',map_by_property='Name'):
        '''
        Return
        ------
        propval2objs = {map_by_property:[PSObject]}
        '''
        my_df = self.RefCountPandas if from_df.empty else from_df
        entity_names = my_df[id_column].to_list()
        propval2objs,_ = self.Graph.get_prop2obj_dic(map_by_property,entity_names)
        return propval2objs


    def expand_entities(self, in_df:df,with_objtypes:list,linked_by=list()):
        my_session = self._clone_session(what2retrieve=NO_REL_PROPERTIES)
        my_session.entProps = []

        oql = 'SELECT Relation WHERE NeighborOf (SELECT Entity WHERE id = ({ids}))'
        oql += f' AND NeighborOf (SELECT Entity WHERE objectType = ({with_objtypes}))'
        if linked_by:
            oql += f' AND objectType = ({linked_by})'

        my_dbids = self._all_dbids(in_df)
        reqname = f'Expanding entities with {with_objtypes} linked by {linked_by}'
        expanded_entity_graph = my_session.iterate_oql(oql,my_dbids,request_name=reqname)

        def __expand_row(cell):
            row_dbids = set(list(cell)) # row_dbids is a tuple
            enities2expand = expanded_entity_graph.psobj_with_dbids(row_dbids)
            enity_neighbors = expanded_entity_graph.get_neighbors(set(enities2expand))
            row_dbids.update(ResnetGraph.dbids(list(enity_neighbors)))
            return tuple(list(row_dbids)), ','.join(ResnetGraph.names(list(enity_neighbors))),len(enity_neighbors)

        my_df_copy = df.copy_df(in_df)
        my_df_copy[[self.__temp_id_col__,'ExpandedBy','Expand Size']] = my_df_copy[self.__temp_id_col__].apply(lambda x: pd.Series(__expand_row(x)))
        return my_df_copy


    def load_pandas(self,from_entity_df:df,prop_names_in_header=True,use_cache=False,
                    map2type:list=[],max_children_count=0,expand2=[],with_rels=[]):
        '''
        Input
        -----
        if prop_names_in_header=True uses column headers as mapping property\n
        iterates through all columns in for mapping
        '''
        self.__only_map2types__ = map2type
        if use_cache:
            try: 
                self.read_cache()
                refcount_df = self.RefCountPandas
            except FileNotFoundError:
                print ('Cannot find cache file %s.  You will have to map entities on input identifiers again!' % self.__cntCache__)
                refcount_df = self.__map_entities(from_entity_df,prop_names_in_header)
        else:
            refcount_df = self.__map_entities(from_entity_df,prop_names_in_header)
        
        if max_children_count:
            refcount_df = self.remove_high_level_entities(refcount_df,max_children_count)

        if expand2:
            refcount_df = self.expand_entities(refcount_df,expand2,with_rels)
        
        self.all_entity_dbids.update(self._all_dbids(refcount_df))
        return refcount_df
    

    def add_temp_id(self,to_df:df,map2column='URN',max_child_count=10,max_threads=10):
        mapping_values = to_df[map2column].to_list()
        db_entities = self._props2psobj(mapping_values,[map2column],get_childs=False)
        [o.update_with_value(self.__mapped_by__,o.name()) for o in db_entities]
        
        children,_ = self.load_children4(db_entities,max_childs=max_child_count,max_threads=max_threads)
        for c in children:
            self.Graph.nodes[c.uid()][self.__mapped_by__] = 'Name'
            self.Graph.nodes[c.uid()][self.__resnet_name__] = c.name()

        graph_psobjects = self.Graph._get_nodes(ResnetGraph.uids(db_entities))
        df_psobjects = [o for o in graph_psobjects if len(o.childs()) <= max_child_count] 
        print(f'{len(graph_psobjects)-len(df_psobjects)} entities with > ontology children {max_child_count}  were removed from further calculations')
        urn2tempids = {o.urn():tuple(o.child_dbids()+[o.dbid()]) for o in df_psobjects}

        new_df = to_df.merge_dict(urn2tempids,self.__temp_id_col__,map2column)
        self.all_entity_dbids.update(self._all_dbids(new_df))
        return new_df
        

    def load_df(self,from_entities:list[PSObject],max_child_count=0,max_threads=10): 
        # do not use self.max_ontology_parent instead of max_child_count to avoid complications for session cloning  
        '''
        Input
        -----
        from_entities - [PSObject]
        
        output:
            refcount_df with columns 'Name',
            if max_child_count > 0 refcount_df will have "self.__temp_id_col__" column
        '''
        [o.update_with_value(self.__mapped_by__,o.name()) for o in from_entities] 
        if max_child_count:
            if from_entities[0].is_from_rnef():# database ids were not loaded
                self.load_dbids4(from_entities)

            # load_children4() adds empty PSObject in CHILDS property if parents has children > max_child_count 
            children,_ = self.load_children4(from_entities,
                                                max_childs=max_child_count,
                                                max_threads=max_threads)
            for c in children:
                self.Graph.nodes[c.uid()][self.__mapped_by__] = 'Name'
                self.Graph.nodes[c.uid()][self.__resnet_name__] = c.name()

            # load_children4() adds empty PSObject in CHILDS property if parents has children > max_child_count 
            graph_psobjects = self.Graph._get_nodes(ResnetGraph.uids(from_entities))
             # child_dbids() will return empty list if o[CHILDS] has empty PSObject in CHILDS property
            df_psobjects = [o for o in graph_psobjects if len(o.childs()) <= max_child_count]
            print(f'{len(graph_psobjects)-len(df_psobjects)} entities with > {max_child_count} ontology children were removed from further calculations')

            # "graph_psobjects" and "from_entities" may not have the same 'Name' !!!!
            names = [o.name() for o in df_psobjects]
            urns = [o.urn() for o in df_psobjects]
            objtypes = [o.objtype() for o in df_psobjects]
            child_dbids = [tuple(o.child_dbids()+[o.dbid()]) for o in df_psobjects]
            assert(len(names) == len(child_dbids))
            refcount_df = df.from_dict({'Name':names,'ObjType':objtypes,'URN':urns,self.__temp_id_col__:child_dbids})
            #refcount_df = self.remove_high_level_entities(refcount_df,max_child_count)
            
            refcount_dbids = self._all_dbids(refcount_df)
            self.all_entity_dbids.update(refcount_dbids)
        else:
            names = [o.name() for o in from_entities]
            urns = [o.urn() for o in from_entities]
            objtypes = [o.objtype() for o in from_entities]
            refcount_df = df.from_dict({'Name':names,'ObjType':objtypes,'URN':urns})
            
        return refcount_df


    def __map_entities(self,EntityPandas:df,prop_names_in_header=False):
        start_mapping_time  = time.time()
        if not prop_names_in_header:
            print('Entity infile does not have header with property names:\nwill use 1st column for mapping as Name and then as Alias')
            EntityPandas.rename(columns={0:'Name'}, inplace=True)
            EntityPandas['Alias'] = EntityPandas['Name']
            self.add_ent_props(['Alias'])
        
        map2types = self.__only_map2types__
        PropName2Prop2psobj = dict()
        RemainToMap = df.copy_df(EntityPandas)
        mapped_count = 0
        for propName in EntityPandas.columns:
            identifiers = list(map(str,list(RemainToMap[propName])))
            prop2dbid2psobj = self.map_prop2entities(identifiers, propName, map2types, get_childs=True)
            PropName2Prop2psobj[propName] = prop2dbid2psobj
            mapped_count = mapped_count + len(prop2dbid2psobj)
            RemainToMap = RemainToMap[~RemainToMap[propName].isin(list(prop2dbid2psobj.keys()))]
        print('%d rows out of %d were mapped to %d entities that have at least one connection in Resnet using %s' % 
             (mapped_count,len(EntityPandas),self.Graph.number_of_nodes(),','.join(EntityPandas.columns)))

        def get_entIDs(col_name,cell_value):
            prop2dbid2psobj = PropName2Prop2psobj[col_name]
            try:
                mapped_objs = list(map(PSObject,prop2dbid2psobj[str(cell_value)].values()))
                resnet_names = ','.join(ResnetGraph.names(mapped_objs))
                all_obj_dbids = set(ResnetGraph.dbids(mapped_objs + ResnetGraph.childs(mapped_objs)))
                return resnet_names,str(cell_value),tuple(all_obj_dbids)
            except KeyError:
                return None,None,None
        
        df_name = EntityPandas._name_ if EntityPandas._name_ else COUNTS
        df2return = df()
        df2return._name_ = df_name
        for propName in EntityPandas.columns:
            MappedEntitiesByProp = df.copy_df(EntityPandas)
            MappedEntitiesByProp = df.apply_and_concat(MappedEntitiesByProp,propName,get_entIDs,
                                        [self.__resnet_name__,self.__mapped_by__,self.__temp_id_col__])
            MappedEntitiesByProp = df.from_pd(MappedEntitiesByProp[MappedEntitiesByProp[self.__temp_id_col__].notnull()])
            df2return = df2return.append_df(MappedEntitiesByProp)

        ex_time = execution_time(start_mapping_time)
        print ('Mapped %d out of %d identifiers to entities in database in %s\n' % 
              (len(df2return), len(EntityPandas.index), ex_time))
        
        return df2return


    def map_prop2entities(self,propValues:list,propName:str,map2types=[],get_childs=False,MinConnectivity=1,max_children=11,max_threads=10)->dict[str,dict[int,PSObject]]:
        """
        Returns
        -------
        prop2psobj = {prop_val:{dbid:PSObject}}
        where PSObject are annotated with self.__mapped_by__,CHILD_DBIDS,CHILD_UIDS properties
        """
        ent_props = self.entProps
        if propName not in ent_props: 
            ent_props.append(propName)

        step = 950 # must be slightly less than 1000 to accomodate names with commas
        iteration_counter = math.ceil(len(propValues) / step)

        print('Will use %d %s identifiers to find entities in %d iterations' % 
             (len(propValues), propName, iteration_counter))

        dbid2entity = dict()
        for i in range(0, len(propValues), step):
            start_time = time.time()
            propval_chunk = propValues[i:i + step]
            query_node = OQL.get_entities_by_props(propval_chunk, [propName], map2types, MinConnectivity)
            zeep_entities = self.get_data(query_node, ent_props, getLinks=False)
            if type(zeep_entities) != type(None):
                dbid2entity_chunk = self._zeep2psobj(zeep_entities)
                dbid2entity.update(dbid2entity_chunk)
                ex_time = execution_time(start_time)
                print("%d in %d iteration found %d entities for %d %s identifiers in %s" % 
                    (i / step + 1, iteration_counter, len(dbid2entity_chunk), len(propval_chunk), propName, ex_time))

       # lazy_child_dict = dict()
        propValues_set = list(map(lambda x: str(x).lower(),propValues))
        def my_prop_vals(psobj:PSObject):
            lower_case_values = set(map(lambda x: str(x).lower(),psobj[propName]))
            intersection = list()
            for propval in lower_case_values:
                try:
                    pval_idx = propValues_set.index(propval)
                    intersection.append(propValues[pval_idx])
                except ValueError:
                    continue
            return intersection

        prop2psobj = dict()
        for psobj in dbid2entity.values():
            prop_values = my_prop_vals(psobj)
            mapped_by_propvalue = propName + ':' + ','.join(prop_values)
            psobj.update_with_value(self.__mapped_by__, mapped_by_propvalue)
            for prop_val in prop_values: 
                try:
                    prop2psobj[prop_val][psobj.dbid()] = psobj
                except KeyError:
                    prop2psobj[prop_val] = {psobj.dbid(): psobj}

        if get_childs:
            children,_ = self.load_children4(list(dbid2entity.values()),
                                                                  max_childs=max_children,
                                                                  max_threads=max_threads)
            dbid2entity.update({child.dbid():child for child in children})
            for c in children:
                self.Graph.nodes[c.uid()][self.__mapped_by__] = 'Name'

        self.Graph.add_nodes_from([(v.uid(), v.items()) for k, v in dbid2entity.items()])
        print('%d out of %d %s identifiers were mapped on entities in the database' % 
             (len(prop2psobj), len(propValues), propName))
        return prop2psobj


    def __refcount_by_dbids(self,node1dbids:list,node2dbids:list,how2connect) -> "ResnetGraph":
        '''
        how2connect - func(node1dbids,node2dbids)\n
        contains instruction how to connect node1dbids with node1dbids using class parameters
        '''
        start_time = time.time()
        graph_connects = how2connect(node1dbids,node2dbids)
        assert(isinstance(graph_connects,ResnetGraph))
        print('%d nodes were linked by %d relations supported by %d references in %s' %
             (len(graph_connects),graph_connects.number_of_edges(),graph_connects.weight(),execution_time(start_time)))
        
        return graph_connects
    

    def connect(self, my_df:df,concepts:list,how2connect) -> ResnetGraph:
        '''
        combines relation from database with relation that exist only in self.Graph and not in database
        such relation can exist in RNEF cache
        '''
        concepts_db_ids = ResnetGraph.dbids(concepts)
        all_df_dbids = list(self._all_dbids(my_df))
        graph_from_db = self.__refcount_by_dbids(concepts_db_ids, all_df_dbids,how2connect)
        all_df_entities = self.Graph.psobj_with_dbids(set(all_df_dbids))
        graph_from_self = self.Graph.get_subgraph(concepts,all_df_entities,self.__connect_by_rels__,self.__rel_effect__,self.__rel_dir__)
        return graph_from_db.compose(graph_from_self)


    def set_how2connect(self,**kwargs):
        '''
        kwargs:
            connect_by_rels - desired relation types for connection. Defaults to []
            with_effects - desired relation effect for connection [positive,negative]. Deafults to [] (any)
            in_dir allowed values: 
                '>' - from entities in row to concept in column
                '<' - from concept in column to entities in row
                '' - any direction (default)
            boost_with_reltypes - if not empty adds to refcount references from other relation types listed in [boost_with_reltypes] regardless Effect sign and direction.  Equivalent to merging relations for refcount if at least one relation specified by other parameters exist between entities and concept
            step = step for iterate_oql() function. Defaults to 500
            nodeweight_prop - indicates the name of the property holding node weightd. Deafaults to '' (no weight property)
        '''

        self.__connect_by_rels__ = kwargs.get('connect_by_rels',[])
        self.__rel_effect__ = kwargs.get('with_effects',[])
        self.__rel_dir__ = kwargs.get('in_dir','')
        self.__boost_with__ = kwargs.get('boost_by_reltypes',[])
        self.iteration_step = kwargs.get('step',500)
        self.nodeweight_prop = kwargs.get('nodeweight_prop','')

        if self.__boost_with__:
            reltypes2load = list(set(self.__connect_by_rels__+self.__boost_with__))
            def connect(node1dbids:list, node2dbids:list):
                return self.connect_nodes(set(node1dbids),set(node2dbids),reltypes2load,step=self.iteration_step)
                # ResnetGraph.relation_exist() checks for self.__connect_by_rels__, self.__rel_effect__
            return connect
        else:
            def connect(node1dbids:list, node2dbids:list):
                return self.connect_nodes(set(node1dbids),set(node2dbids),self.__connect_by_rels__,self.__rel_effect__,self.__rel_dir__,step=self.iteration_step)
            return connect
            

    def __annotate_rels(self, from_graph:ResnetGraph, with_concept_name:str):
        for regulatorID, targetID, rel in from_graph.edges.data('relation'):
            try:
                entity_search_attrs = self.Graph.nodes[regulatorID][self.__mapped_by__]
            except KeyError:
                try:
                    entity_search_attrs = self.Graph.nodes[targetID][self.__mapped_by__]
                except KeyError:
                    continue

            self.Graph.set_edge_annotation(regulatorID, targetID,rel.urn(),self.__mapped_by__, entity_search_attrs)
            self.Graph.set_edge_annotation(regulatorID, targetID,rel.urn(),self.__concept_name__, [with_concept_name])
        return


    def init_concept(self,my_df:df, concept_name:str):
        in_df = df.copy_df(my_df)
        refcount_column = self._refcount_colname(concept_name)
        weighted_refcount_column = self._weighted_refcount_colname(concept_name) 
        linked_count_column = self._linkedconcepts_colname(concept_name)
        concept_size_column = self._concept_size_colname(concept_name)

        in_df.insert(len(in_df.columns),weighted_refcount_column,[float(0)]*len(in_df))
        in_df.insert(len(in_df.columns),refcount_column,[0]*len(in_df))
        in_df.insert(len(in_df.columns),linked_count_column,[0]*len(in_df))
        in_df.insert(len(in_df.columns),concept_size_column,[1]*len(in_df))
        return in_df


    def __link2concept(self,ConceptName:str,concepts:list[PSObject],to_entities:df|pd.DataFrame,how2connect):
        """
        Input
        -----
        to_entities.columns must have self.__temp_id_col__
        concepts - [PSObject]
        how2connect - function with instructions how to connect "concepts","to_entities"

        Returns
        -------
        linked_row_count, linked_entity_ids,
        my_df = with 3 columns added: "Refcount to {ConceptName}", "Linked {ConceptName} count", "{ConceptName} count"
        """
        
        my_df = df.copy_df(to_entities) if isinstance(to_entities,df) else df.from_pd(to_entities)

        concept_size = len(concepts)-1 if concepts else 0
        if (len(concepts) > 501 and len(my_df) > 500):
            print(f'"{ConceptName}" concept has {concept_size} ontology children! Linking may take a while, be patient' )
        else:
            print(f'\nLinking row entities to "{ConceptName}" column which has {concept_size} ontology children')

        refcount_column = self._refcount_colname(ConceptName)
        weighted_refcount_column = self._weighted_refcount_colname(ConceptName) 
        linked_count_column = self._linkedconcepts_colname(ConceptName)
        concept_size_column = self._concept_size_colname(ConceptName)

        concepts_uids = set(ResnetGraph.uids(concepts))
        concept_count = len(concepts_uids)

        if self.nodeweight_prop:
            regurn2weight = {o.urn():o.get_prop('regulator weight') for o in concepts}
            tarurn2weight = {o.urn():o.get_prop('target weight') for o in concepts}

        try:
            my_df.insert(len(my_df.columns),weighted_refcount_column,[float(0)]*len(my_df))
            my_df.insert(len(my_df.columns),refcount_column,[0]*len(my_df))
            my_df.insert(len(my_df.columns),linked_count_column,[0]*len(my_df))
            my_df.insert(len(my_df.columns),concept_size_column,[concept_count]*len(my_df)) # assuming all concepts have at least one member 
            # occurence of linked concepts. Format str(#linked concepts)/str(#total concepts)
        except ValueError:
            print('%s column already exists in dataframe!!!' % refcount_column)
            pass

        linked_row_count = 0
        linked_entities_counter = set()
        start_time  = time.time()
        connection_graph = self.connect(my_df,concepts, how2connect)
   
        if connection_graph.size() > 0:
            self.__annotate_rels(connection_graph, ConceptName)
            ref_sum = set()
            for idx in my_df.index:
                idx_entities = connection_graph.psobj_with_dbids(set(my_df.at[idx,self.__temp_id_col__]))
                if self.__rel_dir__ =='>':
                    row_has_connection = connection_graph.relation_exist(idx_entities,concepts,
                                                self.__connect_by_rels__,self.__rel_effect__,[],False)
                elif self.__rel_dir__ =='<':
                    row_has_connection = connection_graph.relation_exist(concepts,idx_entities,
                                                self.__connect_by_rels__,self.__rel_effect__,[],False)
                else:
                    row_has_connection = connection_graph.relation_exist(concepts,idx_entities,
                                                self.__connect_by_rels__,self.__rel_effect__,[],True)

                if not row_has_connection:
                    continue

                row_subgraph = connection_graph.get_subgraph(idx_entities, concepts)
                connected_concepts = [c for c in row_subgraph.nodes() if row_subgraph.degree(c) and c in concepts_uids]
                my_df.at[idx,linked_count_column] = len(connected_concepts)

                if self.nodeweight_prop:
                    row_subgraph.add_node_weight2ref(regurn2weight,tarurn2weight)
                    
                references = list(row_subgraph.load_references(self.relweight_prop,self.relweight_dict))
                if self.relweight_prop:
                    ref_relweights = [r.get_weight('relweight') for r in references]
                    if self.nodeweight_prop:
                        nodeweights = [r.get_weight('nodeweight') for r in references]
                        ref_weights =  [x + y for x, y in zip(ref_relweights, nodeweights)]
                        row_score = float(sum(ref_weights))
                elif self.nodeweight_prop:
                        ref_weights = [r.get_weight('nodeweight') for r in references]
                        row_score = float(sum(ref_weights))
                else:
                    row_score = len(references)
                
                my_df.at[idx,weighted_refcount_column] = row_score
                my_df.at[idx,refcount_column] = len(references)
                # concept incidence measures the occurence of concepts linked to row entities among all input concepts
               
                if references:
                    ref_sum.update(references)
                    linked_row_count += 1
                    linked_entities_counter.update(idx_entities)

            effecStr = ','.join(self.__rel_effect__) if len(self.__rel_effect__)>0 else 'all'
            relTypeStr = ','.join(self.__connect_by_rels__) if len(self.__connect_by_rels__)>0 else 'all'
            exec_time = execution_time(start_time)
            if linked_row_count:
                print("Concept \"%s\" is linked to %d entities by %s relations of type \"%s\" supported by %d references with effect \"%s\" in %s" %
                    (ConceptName,linked_row_count, connection_graph.number_of_edges(),relTypeStr,len(ref_sum),effecStr,exec_time))
            elif connection_graph:
                print("Concept \"%s\" has no links of type \"%s\" with effect \"%s\" to entities in boosted graph" %
                (ConceptName,relTypeStr,effecStr,))
        else:  # connection_graph is empty
            print("Concept \"%s\" has no links to entities" % (ConceptName))

        return linked_row_count, linked_entities_counter, my_df


    def link2concept(self,to_concept_named:str,concepts:list,to_entities:df|pd.DataFrame,how2connect,clone2retrieve=DO_NOT_CLONE):
        """
        Input
        -----
        to_entities.columns must have self.__temp_id_col__
        concepts - [PSObject]
        how2connect - function with instructions how to connect "concepts","to_entities"

        Returns
        -------
        linked_row_count, linked_entity_ids, my_df = to_entities|ConceptName
        """
        if clone2retrieve:
            my_session = self._clone(to_retrieve=clone2retrieve) # careful to_retrieve MPSVI ClinicalTrial
            linked_row_count,linked_entity_ids,return_df = my_session.__link2concept(
                                            to_concept_named,concepts,to_entities,how2connect)
            my_session.close_connection()
            return linked_row_count,linked_entity_ids,return_df
        else: 
            return self.__link2concept(to_concept_named,concepts,to_entities,how2connect)


    def link2RefCountPandas(self,to_concept_named:str,concepts:list,how2connect=set_how2connect,clone2retrieve=DO_NOT_CLONE):
        '''
        Input
        -----
        concepts - [PSObject]
        wrapper for backward compatibility
        '''
        linked_row_count,linked_entity_ids,self.RefCountPandas = self.link2concept(
                        to_concept_named,concepts,self.RefCountPandas,how2connect,clone2retrieve)

        return linked_row_count


    def flush_dump(self):
        self.flush_dump_files()
        open(self.__cntCache__, 'w').close()
        open(self.__refCache__, 'w').close()


    def add2report(self,table:df):
        table_name = table._name_
        if not table_name:
            table_name = 'Sheet' + str(len(self.report_pandas))
        self.report_pandas[str(table_name)] = table


    def add2raw(self,table:df):
        assert(table._name_)
        self.raw_data[table._name_] = table


    def find_ref(self,ref:Reference):
        id_type, identifier = ref.get_doc_id()
        try:
            return self.references[id_type+':'+identifier]
        except KeyError:
            return dict()


    def add2writer(self,writer:pd.ExcelWriter,ws_prefix='',df_names:list[str]=[]):
        '''
        input:
            self.columns2drop
        '''
        if df_names:
            my_pandas = {k:self.report_pandas[k] for k in df_names if k in self.report_pandas}
        else:
            my_pandas = self.report_pandas

        for ws_name, report_df in my_pandas.items():
            df2print = df.copy_df(report_df)
            report_df_columns = set(df2print.columns)
            my_drop = [c for c in self.columns2drop if c in report_df_columns]
            clean_df = df2print.drop(columns=my_drop)
            clean_df = df.from_pd(clean_df.loc[(clean_df!=0).any(axis=1)]) # removes rows with all zeros
            clean_df.copy_format(report_df)
            if 'bibliography' in ws_name:
                clean_df.make_header_horizontal()
            else:
                clean_df.make_header_vertical()
            sh_name = ws_prefix+'_'+ws_name if ws_prefix else ws_name
            clean_df.df2excel(writer, sheet_name=sh_name[:30])


    def addraw2writer(self, writer:pd.ExcelWriter, ws_prefix=''):
        if hasattr(self,'raw_data'):
            for ws_name, rawdf in self.raw_data.items():
                sh_name = ws_prefix+'_'+ws_name if ws_prefix else ws_name
                raw_df = df.copy_df(rawdf)
                raw_df.make_header_vertical()
                raw_df.df2excel(writer, sheet_name=sh_name[:30])


    def print_report(self, xslx_file:str, ws_prefix=''):
        fout = self.data_dir+xslx_file
        writer = pd.ExcelWriter(fout, engine='xlsxwriter')
        self.add2writer(writer,ws_prefix)
        writer.close()


    def print_rawdata(self, xslx_file:str, ws_prefix=''):
        fout = self.data_dir+xslx_file
        writer = pd.ExcelWriter(fout, engine='xlsxwriter')
        self.addraw2writer(writer,ws_prefix)
        writer.close()


    def read_cache(self):# returns last concept linked in cache
        try:
            self.RefCountPandas = df.from_pd(pd.read_csv(self.__cntCache__,sep='\t',header=0))
            last_col_name = list(self.RefCountPandas.columns)[-1]
            return last_col_name[len(self._col_name_prefix):]
        except FileNotFoundError:
            print('Cache was not found! Will start processing all input files from scratch')
            return FileNotFoundError


    def make_count_df(self,from_rawdf=df(),with_name=COUNTS,sort_by_1stcol=True):
        '''
        Returns
        -------
        df with_name from_df with formatted CHILDREN_COUNT column, soreted by first refcount_column
        '''
        my_df = self.RefCountPandas if from_rawdf.empty else from_rawdf
        count_df = df.copy_df(my_df)
        if self.__temp_id_col__ in count_df.columns:
            def sortlen(x):
                return 0 if x is np.NaN else max(len(list(x)),1)
            count_df[CHILDREN_COUNT] = count_df[self.__temp_id_col__].apply(lambda x: sortlen(x))
            count_df.add_column_format(CHILDREN_COUNT,'align','center')

        if sort_by_1stcol:
            def __col2sort__(my_df:df):
                refcount_columns = self.refcount_columns(my_df)
                if refcount_columns: return refcount_columns[0]
                else:
                    for c in my_df.columns.tolist():
                        if my_df.is_numeric(c): return c

                    return str(my_df.columns[1])

            sort_by = __col2sort__(my_df)
            count_df.sort_values(sort_by,ascending=False,inplace=True)

        if self.__colname4GV__ in count_df.columns:
            count_df.move_cols({self.__colname4GV__:len(count_df.columns)})

        count_df.copy_format(my_df)
        count_df._name_ = with_name
        print (f'Created {count_df._name_} table with {len(count_df)} rows' )
        return count_df
        

    '''
    def __merge_counts2norm(self,rawcounts_df:df,normalized_df:df,entity_id_column='Name',columns2merge=list()):
        """
        Replaces normalized score values in "normalized_df" by original refcounts from "counts_df"
        """
        dict2rename = dict()
        for c in normalized_df.columns.to_list():
            if c.startswith(self._col_name_prefix):
                dict2rename[c] = c[len(self._col_name_prefix)+1:]

        merged_df = df.copy_df(normalized_df,rename2=dict2rename)
        if columns2merge:
            common_columns = columns2merge
        else:
            common_columns = set(merged_df.columns).intersection(set(rawcounts_df.columns))

        for col in common_columns:
            merge_dict = pd.Series(rawcounts_df[col].values,index=rawcounts_df[entity_id_column]).to_dict()
            merged_df[col] = merged_df[entity_id_column].map(merge_dict)
            if pd.api.types.is_integer_dtype(rawcounts_df[col]):
                merged_df[col] = merged_df[col].astype(int)
        return merged_df
    '''

    def normalize(self,raw_df_named:str,to_df_named:str,entity_column='Name',columns2norm:list[str]=[],drop_empty_columns:list[str]=[]):
        """
        input:
            df with name "raw_df_named" must be in self.raw_data
            if "columns2norm" is empty uses list of columns with prefix WEIGHTED
        output:
            rankedf4report with name "to_df_named" in self.report_pandas[to_df_named]
            normdf4raw with name "norm.to_df_named" in self.raw_data[norm.to_df_named]
        """
        try:
            if self.raw_data[raw_df_named].empty:
                print(f'{raw_df_named} worksheet is empty')
                return df(),df()
        except KeyError:
            print(f'No worksheet named "{raw_df_named}" is available in raw_data')
            return df(),df()
        
        my_raw_df = self.raw_data[raw_df_named]
        assert(isinstance(my_raw_df,df))
        counts_df = df(name= my_raw_df._name_)
        rawdf_colnames = my_raw_df.columns.to_list()
        i = 0
        while i < len(rawdf_colnames): # in range() controls i and cannot be used here
            colname = rawdf_colnames[i]
            if colname.startswith(self._col_name_prefix):
                linked_concepts_count_col = rawdf_colnames[i+2]
                max_concept_count = max(my_raw_df[linked_concepts_count_col])
                if max_concept_count >= 1:# ensures that refcount column with all zeros is not copied
                    counts_df[colname] = my_raw_df[colname]
                    refcount_colname = rawdf_colnames[i+1]
                    counts_df[refcount_colname] = my_raw_df[refcount_colname]
                
                    concepts_count_col = rawdf_colnames[i+3]
                    if max_concept_count > 1: # boosts multicomponent concepts
                        counts_df[colname] = my_raw_df[colname]*(1+my_raw_df[linked_concepts_count_col]/my_raw_df[concepts_count_col])
                    
                    concept_name = self.__concept(colname)
                    incidence_colname = concept_name + ' incidence'
                    temp_df = df()
                    temp_df['percentage'] = 100*my_raw_df[linked_concepts_count_col]/my_raw_df[concepts_count_col]
                    temp_df['percentage'] = temp_df['percentage'].round(decimals=2)
                    temp_df['percentage_str'] = temp_df['percentage'].apply(lambda x: f"{x:.2f}")
                    counts_df[incidence_colname] = my_raw_df[linked_concepts_count_col].astype(str)
                    counts_df[incidence_colname] =  counts_df[incidence_colname] + '/' +my_raw_df[concepts_count_col].astype(str)
                    counts_df[incidence_colname] = counts_df[incidence_colname] +'('+temp_df['percentage_str']+')'
                i += 4
            else:
                if not colname.startswith(self._col_name_prefix): 
                    # skipping real refcount columns because theyr are not used for ranking but only used for final report
                    counts_df[colname] = my_raw_df[colname]
                i += 1

        refcount_cols = columns2norm if columns2norm else self.refcount_columns(counts_df)
        refcount_header = [entity_column]+refcount_cols # first adding refcount columns that start with "weighed "
        weights_df = df(columns=refcount_header)
        weights_df.at[0,entity_column] = 'WEIGHTS:'
        # calculating weights for each column in refcount_cols:
        weight_index = 0
        number_of_weights = len(refcount_cols)
        for col in refcount_cols:
            column_weight = (number_of_weights-weight_index)/number_of_weights
            weights_df.at[0,col] = column_weight
            weight_index += 1

        #calculating weighted cumulative score
        normdf4raw = df.copy_df(counts_df,refcount_header)
        normdf4raw = normdf4raw.l2norm(refcount_cols)
        weights = weights_df.loc[0, refcount_cols].values.tolist()
        for i in normdf4raw.index:
            row_scores = normdf4raw.loc[i,refcount_cols].values.tolist()
            assert(len(weights) == len(row_scores))
            weighted_sum = sum(s*w for s,w in zip(row_scores, weights))
            normdf4raw.loc[i,'Combined score'] = weighted_sum

        # moving all other columns including CHILDREN_COUNT to normdf4raw 
        # except count columns that were used to create incidence column
        added_columns = set(normdf4raw.columns)
        for col in list(counts_df.columns):
            if col not in added_columns and ' count' not in col:
                normdf4raw[col] = counts_df[col]

        normdf4raw[RANK] = normdf4raw['Combined score']/normdf4raw[CHILDREN_COUNT]
        normdf4raw = normdf4raw.sort_values(by=[RANK,entity_column],ascending=False)

        normdf4raw = df.from_pd(normdf4raw.loc[normdf4raw[RANK] >= 0.001]) # removes rows with all zeros
        print(f'Removed {len(counts_df)-len(normdf4raw)} rows from normalized worksheet with score=0')
        # prettyfying scores in df:
        normdf4raw['Combined score'] = normdf4raw['Combined score'].map(lambda x: '%2.3f' % x)
        normdf4raw[RANK] = normdf4raw[RANK].map(lambda x: '%2.3f' % x)
            
        if drop_empty_columns:
            normdf4raw = normdf4raw.drop_empty_columns() 

        rename = dict()
        header4rankedf = [entity_column]
        for c in refcount_header:
            if c.startswith(self._col_name_prefix):
                refcount_colname = c[len(self._col_name_prefix):]
                rename[c] = refcount_colname
                header4rankedf.append(refcount_colname)
        #prettyfying weighter_df header:
        weights_df = df.from_pd(weights_df.map(lambda x: f'{x:,.4f}' if isinstance(x,float) else x))
        # renaming weight_df columns from "weighted " to "RefCount to ":
        weights_df = df.copy_df(weights_df,rename2=rename) 

        normdf4raw_cols = normdf4raw.columns.to_list()
        forbidden_cols = set(refcount_header+header4rankedf+[self.__temp_id_col__,'URN',RANK,'ObjType',self.__colname4GV__])
        other_cols4rankedf = [x for x in normdf4raw_cols if x not in forbidden_cols]
        header4rankedf += other_cols4rankedf
        header4rankedf +=  [RANK,'URN','ObjType']
        if self.__colname4GV__ in normdf4raw_cols:
            header4rankedf.append(self.__colname4GV__)
        # at this point: header4rankedf = [entity_column]+refcount_colnames+other_cols4rankedf+[RANK,'URN','ObjType',self.__colname4GV__]

        rankedf4report = df.copy_df(normdf4raw,only_columns=header4rankedf)
        # now finish normdf4raw by adding weights to df and df to raw data:
        normdf4raw = weights_df.append_df(normdf4raw)
        normdf4raw._name_ = 'norm.'+counts_df._name_
        self.add2raw(normdf4raw)

        assert(weights_df.columns.to_list() == rankedf4report.columns.to_list()[:len(weights_df.columns)])
        rankedf4report = weights_df.append_df(rankedf4report)# append will add columns missing in weights_str 
        rankedf4report = rankedf4report.move_cols({RANK:1})

        rankedf4report.copy_format(counts_df)
        rankedf4report.add_column_format(RANK,'align','center')
        rankedf4report.add_column_format('Combined score','align','center') 
        rankedf4report._name_ = to_df_named
        self.add2report(rankedf4report)
        
        return rankedf4report,normdf4raw


    def etm_refcount_colname(self,between_names_in_col,and_concepts):
        return self.etm_counter.refcount_column_name(between_names_in_col,and_concepts)
    

    def doi_column_name(self,between_names_in_col,and_concepts):
        return self.etm_counter.doi_column_name(between_names_in_col,and_concepts)


    def bibliography(self,to_df:df,input_names:list,entity_name_col:str='Name',add2query=[],max_row=300):
        '''
        Returns
        -------
        copy of to_df with added columns REFS_COLUMN,DOIs
        '''
        if isinstance(self.etm_counter,ETMstat):
            TMsoft = 'ETM' 
            my_query_func = ETMstat.basic_query
        else:
            TMsoft = 'SBS' 
     #       my_query_func = SBSstat.basic_query

        print(f'Finding {self.etm_counter._limit()} most relevant articles in {TMsoft} for {len(to_df)} rows \
in "{to_df._name_} worksheet and {input_names}', flush=True)
        
        return self.etm_counter.add_refs(to_df,entity_name_col,input_names,my_query_func,add2query,max_row)


    def __refs2report(self,to_df_named:str,input_names:list,entity_name_col:str='Name',add2query=[],add2report=True,skip1strow=True):
        """
        input:
            df with "to_df_named" must be in self.report_pandas and have column RANK
        output:
            copy of "to_df_named" with added columns:
                etm.ETM_REFS_COLUMN
                etm._etm_doi_column_name()
        """
        my_df = df.copy_df(self.report_pandas[to_df_named])
        if skip1strow:
            weights = my_df.iloc[[0]].copy()
            my_df = df.from_pd(my_df.drop(0, axis=0,inplace=False), self.report_pandas[to_df_named]._name_)
            if my_df.empty: return my_df
            my_df.copy_format(self.report_pandas[to_df_named])
            my_df = self.bibliography(my_df,input_names,entity_name_col,add2query,max_row=len(my_df))

            #rank_counts_pd = rank_counts_df.sort_values(RANK,ascending=False)
            rank_counts_pd = pd.concat([weights, my_df]).reset_index(drop=True)
            my_df = df.from_pd(rank_counts_pd,self.report_pandas[to_df_named]._name_)     
            my_df.copy_format(self.report_pandas[to_df_named])
        else:
            my_df = self.bibliography(my_df,input_names,entity_name_col,add2query,max_row=len(my_df))
            
        if add2report:
            self.add2report(my_df)

        return my_df


    def add_etm_bibliography(self,suffix=''):
        """
        Adds
        ----
        df with ETM_BIBLIOGRAPHY name to self.report_pandas from self.etm_counter.counter2df()
        """
        biblio_df = self.etm_counter.counter2df()
        biblio_df._name_ = ETM_BIBLIOGRAPHY
        if suffix: biblio_df._name_ += '-'+suffix
        biblio_df._name_ = biblio_df._name_[:31]
        self.add2report(biblio_df)
        return biblio_df._name_
  

    def add_graph_bibliography(self,suffix='',from_graph=ResnetGraph(),add2report=True):
        """
        adds:
            df with PS_BIBLIOGRAPHY-suffix name to self.report_pandas
        """
        my_graph = from_graph if from_graph else self.Graph
        ref_df_name = PS_BIBLIOGRAPHY
        if suffix: ref_df_name += '-'+suffix
        ref_df_name = ref_df_name[:31]
        ref_df = my_graph.bibliography(ref_df_name)
        if add2report:
            self.add2report(ref_df)

        return ref_df


    def id2paths(self,_4df:df, for_entities_in_column='Name', 
                 map_by_graph_property='Name',ontology_depth=3)->dict[str,str]:

        name2child_objs = self.entities(_4df,for_entities_in_column,map_by_graph_property)
        return self.ontopaths2(name2child_objs,ontology_depth)


    def add_groups(self,_2df:df,group_names:list,for_entities_in_column='Name',map_by_graph_property='Name'):
        self.add_group_annotation(group_names)
        prop2objs = self.entities(_2df,for_entities_in_column,map_by_graph_property)
        prop2groups = dict()
        for prop, psobjs in prop2objs.items():
            prop_groups = set()
            for psobj in psobjs:
                obj_groups = psobj.propvalues(BELONGS2GROUPS)
                if obj_groups:
                    prop_groups.update(obj_groups)
            prop2groups[prop] = ',\n'.join(prop_groups)

        new_df = df.copy_df(_2df)
        new_df = new_df.merge_dict(prop2groups,'Groups',for_entities_in_column)
        self.add2report(new_df)
        return new_df
    

    def add_entity_annotation(self,_2column:str,in_df:df,from_node_property='', for_entities_in_column='Name',map_by_graph_property='Name'):
        node_property = from_node_property if from_node_property else _2column
        prop2objs = self.entities(in_df,for_entities_in_column,map_by_graph_property)
        prop2values = dict()
        for prop, psobjs in prop2objs.items():
            prop_values = set()
            for o in psobjs:
                try:
                    o_props = o[node_property]
                    prop_values.update(o_props)
                except KeyError:
                    continue
            
            if prop_values:
                prop2values[prop] = ','.join(prop_values)

        new_df = df.copy_df(in_df)
        new_df = new_df.merge_dict(prop2values,_2column,for_entities_in_column)
        self.add2report(new_df)
        return new_df


    def __load_ontology(self):
        disease_ontology = os.path.join(os.getcwd(),'ENTELLECT_API/ElsevierAPI/ResnetAPI/ontology/disease4ontology_analysis.txt')
        disease_ontology_groups = [x.strip() for x in open(disease_ontology, 'r').readlines()]

        query_node = OQL.get_entities_by_props(disease_ontology_groups, ['Name'])
        new_session = self._clone_session() #to avoid mutating self.Graph in parallel calculations
        disease_ontology_groups_graph = new_session.process_oql(query_node)
        if isinstance(disease_ontology_groups_graph,ResnetGraph):
            disease_ontology_groups = disease_ontology_groups_graph._get_nodes()
        else:
            disease_ontology_groups = list()
        _, ontology = new_session.load_children4(disease_ontology_groups,
                                                                max_threads=50)
        new_session.close_connection()
        return ontology


    def add_ontology_df(self,_4df:df,add2report=True):
        '''
        input:
            df with name = for_df_name must be in self.report_pandas
        output:
            Adds worksheet to "report_pandas" containing statistics for ontology groups listed in 'ElsevierAPI/ResnetAPI/ontology/disease4ontology_analysis.txt'
        '''
        if not self._Ontology:
            self._Ontology = self.__load_ontology()

        try:
            scored_disease_names = _4df[self.__resnet_name__].to_list()
        except KeyError:
            scored_disease_names = _4df['Name'].to_list()  # df made by load_df does not have self.__resnet_name__

        all_diseases = self.Graph._psobjs_with('Disease','ObjTypeName')
        scored_diseases = [x for x in all_diseases if x.name() in scored_disease_names]
        grant_total = len(scored_diseases)
        rows = [['All diseases',grant_total,100]]
        all_scored_children_counter = set()
        for parent_disease in self._Ontology:
            parent_ontology = [parent_disease] + parent_disease.childs()
            scored_children_in_parent_ontology = set(parent_ontology).intersection(scored_diseases)
            group_stat = len(scored_children_in_parent_ontology)
            rows.append([parent_disease.name(),group_stat, 100*group_stat/grant_total])
            all_scored_children_counter.update(scored_children_in_parent_ontology)

        stats_colname = 'Number of indications'
        other_count = len(scored_diseases)-len(all_scored_children_counter)
        rows.append(['Other diseases',other_count, 100*other_count/grant_total])
        ontology_df = df.from_rows(rows,['Ontology category',stats_colname,'Percentage'])
        ontology_df._name_ = ONTOLOGY_ANALYSIS+'4'+ _4df._name_ if _4df._name_ else ONTOLOGY_ANALYSIS

        ontology_df.sort_values(by=stats_colname,ascending=False,inplace=True)
        print('Created "Ontology analysis" table')

        if add2report:
            self.add2report(ontology_df)
        return ontology_df


    def remove_high_level_entities(self,from_df:df, max_children_count=11):
        return_df = df.copy_df(from_df)
        if max_children_count:
            return_df = df.from_pd(return_df[return_df.apply(lambda x: len(x[self.__temp_id_col__]) <= max_children_count, axis=1)])
            print('%d entities with > %d ontology children were removed from further calculations' %
            ( (len(from_df)-len(return_df)), max_children_count-1))
            return_df.copy_format(from_df)
            return_df._name_ = from_df._name_
        return return_df


    def clear(self):
        super().clear()
        self.RefCountPandas = df()
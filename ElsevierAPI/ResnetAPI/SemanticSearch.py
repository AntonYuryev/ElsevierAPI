import time
import math
from .PathwayStudioGOQL import OQL
from ..ETM_API.references import Reference
from .ResnetGraph import ResnetGraph,PSObject,EFFECT,df,CHILDS
from .ResnetAPISession import  APISession,len
from .ResnetAPISession import DO_NOT_CLONE,BELONGS2GROUPS,NO_REL_PROPERTIES,CURRENT_SPECS
from ..ETM_API.etm import ETMstat
import pandas as pd
import numpy as np

COUNTS = 'counts'
PS_BIBLIOGRAPHY = 'PSrefs'
ETM_BIBLIOGRAPHY = 'ETMrefs'
CHILDREN_COUNT = 'Number of ontology children'
RANK = 'Score'
ONTOLOGY_ANALYSIS = 'ontology'

class SemanticSearch (APISession):
    __refCache__='reference_cache.tsv'
    __cntCache__='counts_cache.tsv'
    __concept_name__ = 'Concept name'
    __mapped_by__ = 'mapped_by'
    __resnet_name__ = 'Resnet name'
    _col_name_prefix = "RefCount to "
    __temp_id_col__ = 'entity_IDs'
    __rel_dir__ = '' # relation directions: allowed values: '>', '<',''
    __iter_size__ = 1000 #controls the iteration size in sematic reference count
    __print_refs__ = True
    data_dir = ''
    iteration_step = 1000
    __boost_with__ = list()
    weight_prop = str()
    params = dict()
    

    def __init__(self,*args,**kwargs):
        '''
        Input
        -----
        APIconfig - args[0]\nwhat2retrieve - defaults NO_REL_PROPERTIES, other options -\n
        [DATABASE_REFCOUNT_ONLY,REFERENCE_IDENTIFIERS,BIBLIO_PROPERTIES,SNIPPET_PROPERTIES,ONLY_REL_PROPERTIES,ALL_PROPERTIES]
        connect2server - default True, set to False to run script using data in __pscache__ files instead of database
        '''
        APIconfig = args[0]
        session_kwargs, parameters = self.split_kwargs(kwargs)
        super().__init__(APIconfig,**session_kwargs)
        self.params.update(parameters)

        self.PageSize = 1000
        self.RefCountPandas = df(name=COUNTS) # stores result of semantic retreival

        self.report_pandas=dict() # stores pandas used to generate report file
        self.raw_data = dict()
        self.__only_map2types__ = list()
        self.__connect_by_rels__= list()
        self.__rel_effect__ = list()
        self.references = dict() # {identifier:Reference} made by self.Graph.citation_index()
        self.etm_counter = ETMstat(self.APIconfig,limit=5)
        self.etm_counter.min_relevance = 0.0
        self.all_entity_dbids = set()
        self.weight_dict = dict()
        self.columns2drop = [self.__temp_id_col__] # columns to drop before printing pandas
        self.max_ontology_parent = 11


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
        my_kwargs = dict(kwargs)
        api_session = self._clone_session(**my_kwargs)
        new_session = SemanticSearch(api_session.APIconfig,**my_kwargs)
        new_session.entProps = api_session.entProps
        new_session.relProps = api_session.relProps
        new_session.data_dir = self.data_dir
        new_session.__connect_by_rels__ = self.__connect_by_rels__
        new_session.__rel_effect__ = self.__rel_effect__
        new_session.__rel_dir__ = self.__rel_dir__
        new_session.__boost_with__ = self.__boost_with__
        new_session.weight_prop = self.weight_prop
        new_session.weight_dict = self.weight_dict
        new_session.iteration_step = self.iteration_step
        if self.__rel_effect__:
            new_session.add_rel_props([EFFECT])
        
        return new_session


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


    def _all_dbids(self,df2score=df()):
        """
        Returns
        -------
        list of all ids from the tuples in column self.__temp_id_col__].\n
        if df2score is empty returns all ids from self.RefCountPandas
        """
        if df2score.empty:
            return self.all_entity_dbids
        else:
            list_of_tuples = list(df2score[self.__temp_id_col__])
            id_set = set()
            [id_set.update(list(lst)) for lst in list_of_tuples]
            return id_set


    def entities(self,from_df=df(),id_column='Name',map_by_property='Name'):
        '''
        Return
        ------
        propval2objs = {map_by_property:[PSObject]}
        '''
        my_df = self.RefCountPandas if from_df.empty else from_df
        entity_names = my_df[id_column].to_list()
        propval2objs,objid2propval = self.Graph.get_prop2obj_dic(map_by_property,entity_names)
        return propval2objs


    def load_pandas(self,from_entity_df:df,prop_names_in_header=True,use_cache=False,map2type:list=[],max_children_count=0):
        '''
        Input
        -----
        if prop_names_in_header=True uses column headers as mapping property\n
        iterates through all columns in for mapping
        '''
        self.__only_map2types__ = map2type
        if use_cache:
            try: self.read_cache()
            except FileNotFoundError:
                print ('Cannot find cache file %s.  You will have to map entities on input identifiers again!' % self.__cntCache__)
                refcount_df = self.__map_entities(from_entity_df,prop_names_in_header)
        else:
            refcount_df = self.__map_entities(from_entity_df,prop_names_in_header)
        
        if max_children_count:
            refcount_df = self.remove_high_level_entities(refcount_df,max_children_count)

        refcount_dbids = self._all_dbids(refcount_df)
        self.all_entity_dbids.update(refcount_dbids)
        return refcount_df
    

    def add_temp_id(self,to_df:df,map2column='URN',max_children_count=11):
        mapping_values = to_df[map2column].to_list()
        db_entities = self._props2psobj(mapping_values,[map2column],get_childs=False)
        [o.update_with_value(self.__mapped_by__,o.name()) for o in db_entities]
        
        children = self.load_children4(db_entities,max_childs=self.max_ontology_parent)
        for c in children:
            self.Graph.nodes[c.uid()][self.__mapped_by__] = 'Name'
            self.Graph.nodes[c.uid()][self.__resnet_name__] = c.name()

        graph_psobjects = self.Graph._get_nodes(ResnetGraph.uids(db_entities))
        urn2tempids = {o.urn():tuple(o.child_dbids()+[o.dbid()]) for o in graph_psobjects}

        new_df = to_df.merge_dict(urn2tempids,self.__temp_id_col__,map2column)
        new_df.fillna('', inplace=True) # need to replace NaN to enable cutoff by len(self.__temp_id_col__) next
        if max_children_count:
            new_df = self.remove_high_level_entities(new_df,max_children_count)
        self.all_entity_dbids.update(self._all_dbids(new_df))
        return new_df
        

    def load_df(self,from_entities:list,max_children_count=0):
        '''
        Input
        -----
        from_entities - [PSObject]
        
        Return
        ------
        refcount_df with columns 'Name',\n
        if max_children_count > 0 refcount_df will have "self.__temp_id_col__" column
        '''
        [o.update_with_value(self.__mapped_by__,o.name()) for o in from_entities] 
        if max_children_count:
            if from_entities[0].is_from_rnef():# database ids were not loaded
                self.load_dbids4(from_entities)

            children = self.load_children4(from_entities, min_connectivity=1,max_childs=self.max_ontology_parent)
            for c in children:
                self.Graph.nodes[c.uid()][self.__mapped_by__] = 'Name'
                self.Graph.nodes[c.uid()][self.__resnet_name__] = c.name()

            graph_psobjects = self.Graph._get_nodes(ResnetGraph.uids(from_entities))
            # "graph_psobjects" and "from_entities" may not have the same 'Name' !!!!
            names = [o.name() for o in graph_psobjects]
            urns = [o.urn() for o in graph_psobjects]
            child_dbids = [tuple(o.child_dbids()+[o.dbid()]) for o in graph_psobjects]
            assert(len(names) == len(child_dbids))
            refcount_df = df.from_dict({'Name':names,'URN':urns,self.__temp_id_col__:child_dbids})
            refcount_df = self.remove_high_level_entities(refcount_df,max_children_count)
            
            refcount_dbids = self._all_dbids(refcount_df)
            self.all_entity_dbids.update(refcount_dbids)
        else:
            names = [o.name() for o in from_entities]
            urns = [o.urn() for o in from_entities]
            refcount_df = df.from_dict({'Name':names,'URN':urns})
            
        return refcount_df


    def __map_entities(self,EntityPandas:df,prop_names_in_header=False):
        start_mapping_time  = time.time()
        if not prop_names_in_header:
            print('Entity infile does not have header with property names:\nwill use 1st column for mapping as Name and then as Alias')
            EntityPandas.rename(columns={0:'Name'}, inplace=True)
            EntityPandas['Alias'] = EntityPandas['Name']
            self.add_ent_props(['Alias'])
        
        map2types = self.__only_map2types__
        PropName2Prop2EntityID = dict()
        RemainToMap= df.copy_df(EntityPandas)
        mapped_count = 0
        for propName in EntityPandas.columns:
            identifiers = list(map(str,list(RemainToMap[propName])))
            propVal2psobj = self.map_prop2entities(identifiers, propName, map2types, get_childs=True)
            PropName2Prop2EntityID[propName] = propVal2psobj
            mapped_count = mapped_count + len(propVal2psobj)
            RemainToMap = RemainToMap[~RemainToMap[propName].isin(list(propVal2psobj.keys()))]
        print('%d rows out of %d were mapped to %d entities that have at least one connection in Resnet using %s' % 
             (mapped_count,len(EntityPandas),self.Graph.number_of_nodes(),','.join(EntityPandas.columns)))

        def get_entIDs(col_name,cell_value):
            propVal2psobj = PropName2Prop2EntityID[col_name]
            try:
                mappedPSobj = propVal2psobj[str(cell_value)].values()
                obj_names = [o.name() for o in mappedPSobj]
                resnet_names = ','.join(obj_names)

                all_obj_dbids = {o.dbid() for o in mappedPSobj}
                try:
                    childs = [o[CHILDS] for o in mappedPSobj]
                    child_dbids = ResnetGraph.dbids(childs)
                    all_obj_dbids.update(set().union(*child_dbids))
                    return resnet_names,str(cell_value),tuple(all_obj_dbids)
                except KeyError:
                    return resnet_names,str(cell_value),tuple(all_obj_dbids)
            except KeyError:
                return None,None,None
        
        df_name = EntityPandas._name_ if EntityPandas._name_ else COUNTS
        df2return = df(name=df_name)
        for propName in EntityPandas.columns:
            MappedEntitiesByProp = df(EntityPandas)
            MappedEntitiesByProp = df.apply_and_concat(MappedEntitiesByProp,propName,get_entIDs,
                                        [self.__resnet_name__,self.__mapped_by__,self.__temp_id_col__])
            MappedEntitiesByProp = MappedEntitiesByProp[MappedEntitiesByProp[self.__temp_id_col__].notnull()]
            df2return = df2return.append_df(MappedEntitiesByProp)

        ex_time = self.execution_time(start_mapping_time)
        print ('Mapped %d out of %d identifiers to entities in database in %s\n' % 
              (len(df2return), len(EntityPandas.index), ex_time))
        
        return df2return


    def map_prop2entities(self,propValues:list,propName:str,map2types=[],get_childs=False,MinConnectivity=1):
        """
        Returns
        -------
        prop2psobj = {prop_val:{PSObject.dbid():PSObject}}
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
                ex_time = self.execution_time(start_time)
                print("%d in %d iteration found %d entities for %d %s identifiers in %s" % 
                    (i / step + 1, iteration_counter, len(dbid2entity_chunk), len(propval_chunk), propName, ex_time))

       # lazy_child_dict = dict()
        propValues_set = set(propValues)
        def my_prop_vals(psobj:PSObject):
            lower_case_values = set(map(lambda x: str(x).lower(),psobj[propName]))
            return lower_case_values.intersection(propValues_set)

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
            children = self.load_children4(dbid2entity.values(),max_childs=self.max_ontology_parent)
            dbid2entity.update({child.dbid():child for child in children})
            for c in children:
                self.Graph.nodes[c.uid()][self.__mapped_by__] = 'Name'

        self.Graph.add_nodes_from([(k, v.items()) for k, v in dbid2entity.items()])
        print('%d out of %d %s identifiers were mapped on entities in the database' % 
             (len(prop2psobj), len(propValues), propName))
        return prop2psobj


    def __refcount_by_dbids(self, node1dbids:list, node2dbids:list, step=500) -> "ResnetGraph":
        '''
        step - regulates the speed of retrieval:\n
        iteration through node1ids is done by step size\n
        inetration through node2ids is done 2*step size in parallele threads
        '''
        start_time = time.time()
        if self.__boost_with__:
            reltypes2load = list(set(self.__connect_by_rels__+self.__boost_with__))
            cumulative = self.connect_nodes(set(node1dbids),set(node2dbids),reltypes2load,step=step)
        else:
            cumulative = self.connect_nodes(set(node1dbids),set(node2dbids),self.__connect_by_rels__,self.__rel_effect__,self.__rel_dir__,step=step)
        
        print('%d nodes were linked by %d relations supported by %d references in %s' %
             (cumulative.number_of_nodes(), cumulative.number_of_edges(), cumulative.weight(),self.execution_time(start_time)))
        return cumulative


    def set_how2connect(self,connect_by_rels:list,with_effects:list,in_dir:str,boost_by_reltypes=list(),step=500):
        '''
        Input
        -----
        rel_dir allowed values: 
            '>' - from entities in row to concept in column
            '<' - from concept in column to entities in row
            '' - any direction

        boost_with_reltypes - if not empty adds to refcount refrences from other relations listed in boost_with_reltypes\
            regardless Effect sign and direction.  Equivalent to merging relations for refcount if at least one relation\
                specified by other parameters exist between entities and concept
        '''
        self.__connect_by_rels__ = connect_by_rels
        self.__rel_effect__ = with_effects
        self.__rel_dir__ = in_dir
        self.__boost_with__ = boost_by_reltypes
        self.iteration_step = step


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


    def __link2concept(self,ConceptName:str,concepts:list,to_entities:df|pd.DataFrame):
        """
        Input
        -----
        to_entities.columns must have self.__temp_id_col__
        concepts - [PSObject]

        Returns
        -------
        linked_row_count, linked_entity_ids, my_df = to_entities|ConceptName
        """
        my_df = df(to_entities.copy())
        if isinstance(to_entities,df): my_df._name_ = to_entities._name_

        if (len(concepts) > 500 and len(my_df) > 500):
            print('"%s" concept has %d ontology children! Linking may take a while, be patient' % (ConceptName,len(concepts)-1))
        else:
            print('\nLinking row entities to \"%s\" concept column which has %d ontology children' % (ConceptName,len(concepts)-1))
        
        new_column = self._col_name_prefix + ConceptName
        try:
            my_df.insert(len(my_df.columns),new_column,0)
        except ValueError:
            print('%s column already exists in dataframe!!!' % new_column)
            pass

        if self.weight_prop:
            my_df[new_column] = pd.to_numeric(my_df[new_column], downcast="float")

        linked_row_count = 0
        linked_entities_counter = set()
        start_time  = time.time()
        concepts_db_ids = ResnetGraph.dbids(concepts)
        all_entity_dbids = list(self._all_dbids(my_df))
        relations_graph = self.__refcount_by_dbids(concepts_db_ids, all_entity_dbids,self.iteration_step)
   
        if relations_graph.size() > 0:
            self.__annotate_rels(relations_graph, ConceptName)
            ref_sum = set()
            for idx in my_df.index:
                idx_entities = self.Graph.psobj_with_dbids(list(my_df.at[idx,self.__temp_id_col__]))
                if self.__rel_dir__ =='>':
                    has_connection = relations_graph.relation_exist(idx_entities,concepts,
                                                self.__connect_by_rels__,self.__rel_effect__,[],False)
                elif self.__rel_dir__ =='<':
                    has_connection = relations_graph.relation_exist(concepts,idx_entities,
                                                self.__connect_by_rels__,self.__rel_effect__,[],False)
                else:
                    has_connection = relations_graph.relation_exist(concepts,idx_entities,
                                                self.__connect_by_rels__,self.__rel_effect__,[],True)

                if not has_connection:
                    continue

                if self.weight_prop:
                    references = relations_graph.load_references_between(idx_entities, concepts,self.weight_prop,self.weight_dict)
                    ref_weights = [r.get_weight() for r in references]
                    weighted_count = float(sum(ref_weights))
                    my_df.at[idx,new_column] = weighted_count
                else:
                    references = relations_graph.load_references_between(idx_entities, concepts)
                    my_df.at[idx,new_column] = len(references)

                if len(references) > 0:
                    ref_sum.update(references)
                    linked_row_count += 1
                    linked_entities_counter.update(idx_entities)

            effecStr = ','.join(self.__rel_effect__) if len(self.__rel_effect__)>0 else 'all'
            relTypeStr = ','.join(self.__connect_by_rels__) if len(self.__connect_by_rels__)>0 else 'all'
            exec_time = self.execution_time(start_time)
            print("Concept \"%s\" is linked to %d entities by %s relations of type \"%s\" supported by %d references with effect \"%s\" in %s" %
                 (ConceptName,linked_row_count, relations_graph.number_of_edges(),relTypeStr,len(ref_sum),effecStr,exec_time))

        else: print("Concept \"%s\" has no links to entities" % (ConceptName))
        self.__how2clone__ = DO_NOT_CLONE
        return linked_row_count, linked_entities_counter, my_df


    def link2concept(self,to_concept_named:str,concepts:list,to_entities:df|pd.DataFrame,clone2retrieve=DO_NOT_CLONE):
        if clone2retrieve:
            my_session = self._clone(to_retrive=clone2retrieve)
            linked_row_count,linked_entity_ids,return_df = my_session.__link2concept(
                                            to_concept_named,concepts,to_entities)
            my_session.close_connection()
            return linked_row_count,linked_entity_ids,return_df
        else:
            return self.__link2concept(to_concept_named,concepts,to_entities)


    def link2RefCountPandas(self,to_concept_named:str,concepts:list,clone2retrieve=DO_NOT_CLONE):
        '''
        Input
        -----
        concepts - [PSObject]
        wrapper for backward compatibility
        '''
        if clone2retrieve:
            my_session = self._clone(to_retrive=clone2retrieve)
            linked_row_count,linked_entity_ids,self.RefCountPandas = my_session.link2concept(
                                            to_concept_named,concepts,self.RefCountPandas)
            my_session.close_connection()
        else:
            linked_row_count,linked_entity_ids,self.RefCountPandas = self.link2concept(
                                            to_concept_named,concepts,self.RefCountPandas)

        return linked_row_count


    def flush_dump(self):
        self.flush_dump_files()
        open(self.__cntCache__, 'w').close()
        open(self.__refCache__, 'w').close()


    def add2report(self,table:df):
        self.report_pandas[table._name_] = table


    def add2raw(self,table:df):
        self.raw_data[table._name_] = table


    def find_ref(self,ref:Reference):
        id_type, identifier = ref.get_doc_id()
        try:
            return self.references[id_type+':'+identifier]
        except KeyError:
            return dict()


    def add2writer(self,writer:pd.ExcelWriter,ws_prefix='',only_df_with_names=[]):
        if only_df_with_names:
            my_pandas = {n:d for n,d in self.report_pandas.items() if n in only_df_with_names}
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


    @staticmethod
    def _merge_counts2norm(counts_df:df,normalized_df:df,entity_id_column='Name',columns2merge=list()):
        merged_df = df.copy_df(normalized_df)
        common_columns = columns2merge if columns2merge else set(merged_df.columns).intersection(set(counts_df.columns))
        for col in common_columns:
            merge_dict = pd.Series(counts_df[col].values,index=counts_df[entity_id_column]).to_dict()
            merged_df[col] = merged_df[entity_id_column].map(merge_dict)
        return merged_df


    def read_cache(self):# returns last concept linked in cache
        try:
            self.RefCountPandas = pd.read_csv(self.__cntCache__,sep='\t',header=0)
            last_col_name = list(self.RefCountPandas.columns)[-1]
            return last_col_name[len(self._col_name_prefix):]
        except FileNotFoundError:
            print('Cache was not found! Will start processing all input files from scratch')
            return FileNotFoundError


    def make_count_df(self,from_df=df(),with_name=COUNTS):
        '''
        Returns
        -------
        df with_name from_df with formatted CHILDREN_COUNT column, soreted by first refcount_column
        '''
        my_df = self.RefCountPandas if from_df.empty else from_df

        pandas2print = df(my_df)
        if self.__temp_id_col__ in pandas2print.columns:
            pandas2print[CHILDREN_COUNT] = pandas2print[self.__temp_id_col__].apply(lambda x: max(len(x),1))
            pandas2print.add_column_format(CHILDREN_COUNT,'align','center')

        def __col2sort__(my_df:df):
            refcount_columns = self.refcount_columns(my_df)
            if refcount_columns: return refcount_columns[0]
            else:
                for c in my_df.columns.tolist():
                    if my_df.is_numeric(c): return c

                return str(my_df.columns[1])

        sort_by = __col2sort__(my_df)
        pandas2print.sort_values(sort_by,ascending=False,inplace=True)

        pandas2print.copy_format(from_df)
        pandas2print._name_ = with_name
        print ('Created "%s" table' % pandas2print._name_)
        return pandas2print
        

    def normalize(self,raw_df_named:str,to_df_named:str,entity_column='Name',columns2norm=list(),drop_empty_columns=[]):
        """
        Returns
        -------
        df with _name_ 'norm.raw_df_named' added to self.raw_data.\n
        'norm.raw_df' has normalized values from 'raw_df' df and 'Combined score' and 'RANK' columns\n
        df with ._name_='to_df_named' added to self.report_data where normalized values from 'norm.raw_df' are replaced by original counts from raw_df
        """
        counts_df = df.copy_df(self.raw_data[raw_df_named])
        refcount_cols = columns2norm if columns2norm else self.refcount_columns(counts_df)
        
        # removing empty columns
        no_score_cols = list()
        [no_score_cols.append(col) for col in refcount_cols if counts_df[col].max() < 0.000000000000001]
        if no_score_cols:
            [refcount_cols.remove(c) for c in no_score_cols]
            print('%s columns were excluded from %s because they have no scores' % (no_score_cols,counts_df._name_))
            
        header = [entity_column]+refcount_cols
        
        weights = df(columns=header)
        normalized_count_df = df(columns=header)

        weights.at[0,entity_column] = 'WEIGHTS:'
        normalized_count_df[entity_column] = counts_df[entity_column]

        weight_index = 0
        number_of_weights = len(refcount_cols)
        for col in refcount_cols:
            col_max = counts_df[col].max()
            normalized_count_df[col] = counts_df[col]/col_max
            column_weight = (number_of_weights-weight_index)/number_of_weights
            weights.at[0,col] = column_weight
            weight_index += 1

        #calculating cumulative score  
        combined_scores = list()
        for i in normalized_count_df.index:
            scores_row = normalized_count_df.loc[[i]]
            weighted_sum = 0.0
            for col in refcount_cols:
                weighted_sum = weighted_sum + scores_row[col]*weights.at[0,col]
            combined_scores.append(weighted_sum)

        # moving all other columns to normalized_count_df
        added_columns = set(normalized_count_df.columns)
        for col in list(counts_df.columns):
            if col not in added_columns:
                normalized_count_df[col] = counts_df[col]

        normalized_count_df['Combined score'] = np.array(combined_scores)
        normalized_count_df[RANK] = normalized_count_df['Combined score']/normalized_count_df[CHILDREN_COUNT]
        normalized_count_df = normalized_count_df.sort_values(by=[RANK,entity_column],ascending=False)

        # removes rows with all zeros
        normalized_count_df = df.from_pd(normalized_count_df.loc[normalized_count_df[RANK] > 0.0]) 

        # this converts values to pretty string and has to be performed at the end
        normalized_count_df['Combined score'] = normalized_count_df['Combined score'].map(lambda x: '%2.3f' % x)
        normalized_count_df[RANK] = normalized_count_df[RANK].map(lambda x: '%2.3f' % x)

        if drop_empty_columns:
            normalized_count_df = normalized_count_df.drop_empty_columns()

        # prettifying weights row
        first_row = ['WEIGHTS:']
        for c in refcount_cols:
            weight = weights.loc[0,c]
            first_row.append('{:,.4f}'.format(weight))

        # must initialize with list(weights.columns) to avoid automatic dtype copying
        weight_str = df.from_pd(pd.DataFrame(columns=list(weights.columns)))
        weight_str.loc[0] = first_row
        normalized_count_df = df.from_pd(weight_str.append_df(normalized_count_df))
        weigths_header = list(normalized_count_df.iloc[0].replace(np.nan, ''))
        normalized_count_df.loc[0] = weigths_header

        normalized_count_df._name_ = 'norm.'+counts_df._name_
        self.add2raw(normalized_count_df)

        # merge_counts2norm blanks out weigths_header because counts_df does not have it
        ranked_counts_df = self._merge_counts2norm(counts_df,normalized_count_df,columns2merge=columns2norm)
        ranked_counts_df.loc[0] = weigths_header
        
        ranked_counts_df.copy_format(counts_df)
        ranked_counts_df.add_column_format(RANK,'align','center')
        ranked_counts_df.add_column_format('Combined score','align','center')
        ranked_counts_df._name_ = to_df_named
        self.add2report(ranked_counts_df)
        return ranked_counts_df,normalized_count_df


    def etm_ref_column_name(self,between_names_in_col,and_concepts):
        return self.etm_counter._etm_ref_column_name(between_names_in_col,and_concepts)
    
    def etm_doi_column_name(self,between_names_in_col,and_concepts):
        return self.etm_counter._etm_doi_column_name(between_names_in_col,and_concepts)


    def etm_refs2df(self,to_df:df,input_names:list,entity_name_col:str='Name',add2query=[]):
        print('Finding %d most relevant articles in ETM for %d rows in %s and %s' 
                % (self.etm_counter._limit(),len(to_df),to_df._name_,input_names), flush=True)
        return self.etm_counter.add_etm_references(to_df,entity_name_col,input_names,ETMstat.basic_query,add2query)


    def add_etm_refs(self,to_df_named:str,input_names:list,entity_name_col:str='Name',add2query=[],add2report=True):
        """
        Adds
        ----
        columns etm.ETM_REFS_COLUMN, etm._etm_doi_column_name() to df with name "to_df_name" from self.report_pandas
        """
        rank_counts_df = df.copy_df(self.report_pandas[to_df_named])
        rank_counts_df = self.etm_refs2df(rank_counts_df,input_names,entity_name_col,add2query)
        if add2report:
            self.add2report(rank_counts_df)

        return rank_counts_df


    def add_etm_bibliography(self,suffix=''):
        """
        Adds
        ----
        df with ETM_BIBLIOGRAPHY name to self.report_pandas
        """
        biblio_df = self.etm_counter.counter2df()
        biblio_df._name_ = ETM_BIBLIOGRAPHY+'-'+suffix
        biblio_df._name_ = biblio_df._name_[:31]
        self.add2report(biblio_df)
        return biblio_df._name_
  

    def add_ps_bibliography(self,suffix='',from_graph=ResnetGraph(),add2report=True):
        my_graph = from_graph if from_graph else self.Graph
        ref_df_name = PS_BIBLIOGRAPHY+'-'+suffix
        ref_df_name = ref_df_name[:31]
        ref_df = my_graph.refs2df(ref_df_name)
        if add2report:
            self.add2report(ref_df)

        return ref_df


    def id2paths(self,to_report_named:str, for_entities_in_column='Name', map_by_graph_property='Name',ontology_depth=3):
        to_df = self.report_pandas[to_report_named]
        id2child_objs = self.entities(to_df,for_entities_in_column,map_by_graph_property)
        return self.ontopaths2(id2child_objs,ontology_depth)


    def add_parent_column(self,to_report_named:str, _4entities_in_column='Name', map_by_graph_property='Name',ontology_depth=3):
        '''
        DEPRICATED
        '''
        id2paths = self.id2paths(to_report_named,_4entities_in_column,map_by_graph_property,ontology_depth)

        new_df = df.copy_df(self.report_pandas[to_report_named])
        new_df = new_df.merge_dict(id2paths,'Ontology parents',_4entities_in_column)
        #self.add2report(new_df)
        return new_df
    

    def add_groups(self,_2df:df,group_names:list,for_entities_in_column='Name',map_by_graph_property='Name'):
        self.add_group_annotation(group_names)
        prop2objs = self.entities(_2df,for_entities_in_column,map_by_graph_property)
        prop2groups = dict()
        for prop, psobjs in prop2objs.items():
            prop_groups = set()
            for psobj in psobjs:
                obj_groups = psobj.get_prop(BELONGS2GROUPS)
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


    def add_ontology_df(self, for_df_name='', add2report=True):
        '''
        Return
        ------
        Adds worksheet to "report_pandas" containing statistics for ontology groups listed in 'ElsevierAPI/ResnetAPI/ontology/disease4ontology_analysis.txt'
        '''
        my_df = self.report_pandas[for_df_name] if for_df_name else self.RefCountPandas

        disease_ontology = 'ElsevierAPI/ResnetAPI/ontology/disease4ontology_analysis.txt'
        disease_ontology_groups = [x.strip() for x in open(disease_ontology, 'r').readlines()]

        query_node = OQL.get_entities_by_props(disease_ontology_groups, ['Name'])
        new_session = self._clone_session() #to avoid mutating self.Graph in parallel calculations
        disease_ontology_groups = new_session.process_oql(query_node)._get_nodes()
        new_session.load_children4(disease_ontology_groups)

        try:
            scored_disease_names = my_df[self.__resnet_name__].to_list()
        except KeyError:
            # df made by load_df does not have self.__resnet_name__
            scored_disease_names = my_df['Name'].to_list()

        all_diseases = self.Graph._psobjs_with('Disease','ObjTypeName')
        scored_diseases = [x for x in all_diseases if x.name() in scored_disease_names]
        ontology_stats = {'All diseases': len(scored_diseases)}
        all_scored_children_counter = set()
        for parent_disease in disease_ontology_groups:
            parent_ontology = [parent_disease] + parent_disease.childs()
            scored_children_in_parent_ontology = set(parent_ontology).intersection(scored_diseases)
            group_stat = len(scored_children_in_parent_ontology)
            ontology_stats[parent_disease.name()] = group_stat
            all_scored_children_counter.update(scored_children_in_parent_ontology)

        ontology_stats['Other diseases'] = len(scored_diseases)-len(all_scored_children_counter)
        ontology_df = df.from_dict2(ontology_stats,'Ontology category','Number of diseases')
        ontology_df._name_ = ONTOLOGY_ANALYSIS+'4'+ for_df_name if for_df_name else ONTOLOGY_ANALYSIS

        ontology_df.sort_values(by='Number of diseases',ascending=False,inplace=True)
        print('Created "Ontology analysis" table')

        if add2report:
            self.add2report(ontology_df)

        new_session.close_connection()
        return ontology_df


    def remove_high_level_entities(self,from_df:df, max_children_count=11):
        return_df = df.copy_df(from_df)
        if max_children_count:
            return_df = df(return_df[return_df.apply(lambda x: len(x[self.__temp_id_col__]) <= max_children_count, axis=1)])
            print('%d entities with > %d ontology children were removed from further calculations' %
            ( (len(from_df)-len(return_df)), max_children_count-1))
            return_df.copy_format(from_df)
            return_df._name_ = from_df._name_
        return return_df

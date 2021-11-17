import os
import time
import pandas as pd
import networkx as nx
import math
import re
import ElsevierAPI.ResnetAPI.PathwayStudioGOQL as OQL
from ElsevierAPI.ResnetAPI.NetworkxObjects import REF_ID_TYPES
from ElsevierAPI.ResnetAPI.ResnetGraph import ResnetGraph
from ElsevierAPI.ResnetAPI.ResnetAPISession import APISession


class SemanticSearch (APISession): 
    __refCache__='reference_cache.tsv'
    __cntCache__='counts_cache.tsv'
    __concept_name__ = 'Concept name'
    __mapped_by__ = 'mapped_by'
    _col_name_prefix = "RefCount to "
    __child_ids__ = 'Child Ids'
    __temp_id_col__ = 'entity_IDs'
    __connect_by_rels__= list()
    __rel_effect__ = list()
    __rel_dir__ = '' # allowed values: '>', '<',''
    __only_map2types__ = list()
    __iter_size__ = 500 #controls the iteration size in sematic reference count
    RefCountPandas = pd.DataFrame()
    relProps = list(REF_ID_TYPES) # if need_references() add here ['Name','Sentence','PubYear','Title']

    def __init__(self, APIconfig):
        super().__init__(APIconfig['ResnetURL'], APIconfig['PSuserName'], APIconfig['PSpassword'])
        self.PageSize = 500
        self.DumpFiles = []

    def need_references(self):
        return 'Sentence' in self.relProps

    def load_pandas(self, EntityPandas:pd.DataFrame,prop_names_in_header, use_cache=False, map2type = None ):
        self.__only_map2types__ = map2type if isinstance(map2type,list) else []
        if use_cache:
            try: self.read_cache()
            except FileNotFoundError:
                print ('Cannot find cache file %s.  You will have to map entities on input identifiers again!' % self.__cntCache__)
                self.all_entity_ids = self.map_entities(EntityPandas,prop_names_in_header)
        else:
            self.all_entity_ids = self.map_entities(EntityPandas,prop_names_in_header)
            
    def set_how2connect(self, connect_by_rels, rel_effect, rel_dir):
        self.__connect_by_rels__ = connect_by_rels
        self.__rel_effect__ = rel_effect
        self.__rel_dir__ = rel_dir

    def map_entities(self,EntityPandas:pd.DataFrame, prop_names_in_header = False):
        start_mapping_time  = time.time()
        if not prop_names_in_header:
            print('Entity infile does not have header with property names:\nwill use 1st column for mapping as Name and then as Alias')
            EntityPandas.rename(columns={0:'Name'}, inplace=True)
            EntityPandas['Alias'] = EntityPandas['Name']
            self.add_ent_props(['Alias'])
        
        map2types = self.__only_map2types__
        PropName2Prop2EntityID = dict()
        RemainToMap = pd.DataFrame(EntityPandas)
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
                obj_names = [o['Name'] for o in mappedPSobj]
                name_union = set().union(*obj_names)
                resnet_names = ','.join(name_union)

                obj_ids = [o['Id'] for o in mappedPSobj]
                all_obj_ids = set().union(*obj_ids)
                try:
                    child_ids = [o[self.__child_ids__] for o in mappedPSobj]
                    all_obj_ids.update(set().union(*child_ids))
                    return resnet_names,str(cell_value),tuple(all_obj_ids)
                except KeyError:
                    return resnet_names,str(cell_value),tuple(all_obj_ids)
            except KeyError:
                return None,None,None
        
        def apply_and_concat(dataframe:pd.DataFrame, field, func, column_names):
            return pd.concat((dataframe,dataframe[field].apply(lambda cell: pd.Series(func(field,cell),index=column_names))),axis=1)

        for propName in EntityPandas.columns:
            MappedEntitiesByProp = pd.DataFrame(EntityPandas)
            MappedEntitiesByProp = apply_and_concat(MappedEntitiesByProp,propName,get_entIDs,['Resnet name',self.__mapped_by__,self.__temp_id_col__])
            MappedEntitiesByProp = MappedEntitiesByProp[MappedEntitiesByProp[self.__temp_id_col__].notnull()]
            self.RefCountPandas = self.RefCountPandas.append(MappedEntitiesByProp)

        ex_time = self.execution_time(start_mapping_time)
        print ('Mapped %d out of %d identifiers to entities in database in %s\n' % 
              (len(self.RefCountPandas), len(EntityPandas.index), ex_time))

        entity_id_tuples = self.RefCountPandas[self.__temp_id_col__].tolist()
        set_list = [set(x) for x in entity_id_tuples]
        return list(set.union(*set_list))


    def map_prop2entities(self, propValues: list, propName: str,map2types=None,get_childs=False,MinConnectivity=1):
        only_object_types = [] if map2types is None else map2types
        ent_props = self.entProps
        step = 1000
        iteration_counter = math.ceil(len(propValues) / step)

        print('Will use %d %s identifiers to find entities in %d iterations' % 
             (len(propValues), propName, iteration_counter))

        id2entity = dict()
        for i in range(0, len(propValues), step):
            start_time = time.time()
            propval_chunk = propValues[i:i + step]
            query_node = OQL.get_entities_by_props(propval_chunk, [propName], only_object_types, MinConnectivity)
            zeep_entities = self.get_data(query_node, ent_props, getLinks=False)
            if type(zeep_entities) != type(None):
                id2entity_chunk = self._zeep2psobj(zeep_entities)
                id2entity.update(id2entity_chunk)
                ex_time = self.execution_time(start_time)
                print("%d in %d iteration found %d entities for %d %s identifiers in %s" % 
                    (i / step + 1, iteration_counter, len(id2entity_chunk), len(propval_chunk), propName, ex_time))

        lazy_child_dict = dict()
        child_id2psobj = dict()
        has_childs = 0
        prop2psobj = dict()
        for psobj in id2entity.values():
            psobj_id = psobj['Id']
            prop_values = [x for x in psobj[propName] if x in propValues]
            mapped_by_propvalue = propName + ':' + ','.join(prop_values)
            psobj.update_with_value(self.__mapped_by__, mapped_by_propvalue)

            if get_childs:
                lazy_key = tuple(psobj_id)
                try:
                    child_ids = lazy_child_dict[lazy_key]
                except KeyError:
                    query_ontology = OQL.get_childs(psobj_id, ['Id'], only_object_types)
                    zeep_entities = self.get_data(query_ontology, ent_props, getLinks=False)
                    if type(zeep_entities) != type(None):
                        has_childs += 1
                        child_ids2entities = self._zeep2psobj(zeep_entities)
                        for child in child_ids2entities.values():
                            child.update_with_value(self.__mapped_by__, mapped_by_propvalue)
                        child_id2psobj.update(child_ids2entities)
                        child_ids = list(child_ids2entities.keys())
                    else:
                        child_ids = []
                    lazy_child_dict[lazy_key] = child_ids

                psobj.update_with_list(self.__child_ids__, child_ids)

            for prop_val in prop_values:
                try:
                    prop2psobj[prop_val][psobj['Id'][0]] = psobj
                except KeyError:
                    prop2psobj[prop_val] = {psobj['Id'][0]: psobj}

        if get_childs:
            print('Childs for %d entities were found in database' % has_childs)
            id2entity.update(child_id2psobj)

        self.Graph.add_nodes_from([(k, v.items()) for k, v in id2entity.items()])
        print('%d out of %d %s identifiers were mapped on entities in the database' % 
             (len(prop2psobj), len(propValues), propName))
        return prop2psobj

    def semantic_refcount_by_ids(self, node1ids: list, node2ids: list, no_mess=True) -> "ResnetGraph":
        rel_effect = self.__rel_effect__
        connect_by_rel_types = self.__connect_by_rels__
        rel_props = self.relProps
        ent_props = self.entProps
        rel_dir = self.__rel_dir__
        step_size = self.__iter_size__

        iteration_counter = 0
        cumulative = ResnetGraph()
        number_of_iterations = math.ceil(len(node1ids)/step_size) * math.ceil(len(node2ids) / step_size)
        print('Semantic linking will be done in %d iterations' % number_of_iterations)
        if no_mess:
            print('Progress report is suppressed. Linking may take a long time - be patient!')

        start_time = time.time()
        for n1 in range(0, len(node1ids), step_size):
            n1end = min(n1 + step_size, len(node1ids))
            n1ids = node1ids[n1:n1end]
            for n2 in range(0, len(node2ids), step_size):
                n2end = min(n2 + step_size, len(node2ids))
                n2ids = node2ids[n2:n2end]
                oql_query = OQL.connect_ids(n1ids, n2ids, connect_by_rel_types, rel_effect, rel_dir)
                relations = self.load_graph_from_oql(oql_query, relation_props=list(rel_props),
                                                           entity_props=list(ent_props))
                iteration_counter += 1
                if isinstance(relations, ResnetGraph):
                    if not no_mess:
                        print('Iteration %d in %d found %d semantic relations between %d nodes supported by %d references'
                            % (iteration_counter, number_of_iterations, relations.number_of_edges(), relations.number_of_nodes(),
                                relations.size()))
                    cumulative = nx.compose(relations, cumulative)

        print('%d nodes were linked by %d relations supported by %d references in %s' %
             (cumulative.number_of_nodes(), cumulative.number_of_edges(), cumulative.size(),self.execution_time(start_time)))
        return cumulative


    def semantic_refcount(self, node1PropValues: list, node1PropTypes: list, node1ObjTypes: list, node2PropValues: list,
                          node2PropTypes: list, node2obj_types: list):
        node1ids = list(self._get_obj_ids_by_props(node1PropValues, node1PropTypes, only_object_types=node1ObjTypes))
        accumulate_reference = set()
        accumulate_relation = ResnetGraph()
        if len(node1ids) > 0:
            node2ids = list(self._get_obj_ids_by_props(node2PropValues, node2PropTypes, True, node2obj_types))
            if len(node2ids) > 0:
                accumulate_reference, accumulate_relation = self.semantic_refcount_by_ids(node1ids, node2ids)

        return accumulate_reference, accumulate_relation

    def link2concept(self,ConceptName,concept_ids:list,no_mess=True):
        print('\nLinking input entities to \"%s\" concept with %d ontology children' % (ConceptName,len(concept_ids)))
        if (len(concept_ids) > 500 and len(self.RefCountPandas) > 500):
            print('%s concept has %d ontology children! Linking may take a while, be patient' % (ConceptName,len(concept_ids)-1))
        
        new_column = self._col_name_prefix + ConceptName
        self.RefCountPandas.insert(len(self.RefCountPandas.columns),new_column,0)

        effecStr = ','.join(self.__rel_effect__) if len(self.__rel_effect__)>0 else 'all'
        relTypeStr = ','.join(self.__connect_by_rels__) if len(self.__connect_by_rels__)>0 else 'all'

        linked_entities_counter = 0
        start_time  = time.time()
        relations = self.semantic_refcount_by_ids(concept_ids, self.all_entity_ids,no_mess)
        if relations.size() > 0:
            for regulatorID, targetID, rel in relations.edges.data('relation'):
                try:
                    entity_search_attr = ','.join(self.Graph.nodes[regulatorID][self.__mapped_by__])
                except KeyError:
                    entity_search_attr = ','.join(self.Graph.nodes[targetID][self.__mapped_by__])

                rel[self.__mapped_by__] = [entity_search_attr]
                rel[self.__concept_name__] = [ConceptName]
      
            ref_sum = set()
            for idx in self.RefCountPandas.index:
                idx_entity_ids = list(self.RefCountPandas.at[idx,self.__temp_id_col__])
                references = relations.count_references_between(idx_entity_ids, concept_ids)
                self.RefCountPandas.at[idx,new_column] = len(references)
                if len(references) > 0:
                    ref_sum.update(references)
                    linked_entities_counter += 1
            
            exec_time = self.execution_time(start_time)
            self.print_ref_count() #prints cache files for temp storage to handle network interruptions
            print("Concept \"%s\" is linked to %d rows from infile by %s relations of type \"%s\" supported by %d references with effect \"%s\" in %s" %
                 (ConceptName, linked_entities_counter, relations.number_of_edges(), relTypeStr, len(ref_sum), effecStr, exec_time))

        else: print("Concept \"%s\" has no links to entities in infile" % (ConceptName))
        return linked_entities_counter
    

    def print_ref_count(self, refCountsOut='', referencesOut='',**kwargs):
        PandasToPrint = pd.DataFrame(self.RefCountPandas)    
  
        if len(refCountsOut)>0: 
            countFile = refCountsOut
            PandasToPrint = PandasToPrint.drop(self.__temp_id_col__,axis=1)
            PandasToPrint = PandasToPrint.loc[(PandasToPrint!=0).any(axis=1)] # removes rows with all zeros
        else:
            print('Data is saved data to %s for protection or re-use' % self.__cntCache__)
            countFile = self.__cntCache__

        if 'index' not in kwargs:
            kwargs['index'] = False

        if 'sep' not in kwargs:
            kwargs['sep'] = self.csv_delimeter

        PandasToPrint.to_csv(countFile,**kwargs)

        if self.need_references():
            if len(referencesOut) > 0:
                #ref_pandas = self.to_pandas()
                temp_fname = '__temp__.tsv'
                sep = '\t'
                self.Graph.print_references(temp_fname,self.relProps,self.entProps,col_sep=sep)

                search_entity_names = list(self.RefCountPandas['Resnet name'])
                f = open(temp_fname,'r',encoding='utf-8')
                with open(referencesOut, 'w',encoding='utf-8') as fout:
                    header = f.readline()
                    col_names = re.split('\t|\n',header)
                    regulator_column = col_names.index('Regulator:Name')
                    target_column = col_names.index('Target:Name')
                    fout.write('Search Entity'+'\t'+header)

                    for line in f.readlines():
                        columns = re.split('\t|\n',line)
                        search_entity = '' 
                        if columns[regulator_column] in search_entity_names:
                            search_entity = columns[regulator_column]
                        elif columns[target_column] in search_entity_names:
                            search_entity = columns[target_column]
                        fout.write(search_entity+'\t'+line)
                f.close()
                os.remove(temp_fname)
            else:
                self.Graph.print_references(self.__refCache__,relPropNames=self.relProps,entity_prop_names=['Name'],access_mode='a')
                print('Supporting semantic triples are saved to %s for protection or re-use' % self.__refCache__)
 

    def read_cache(self):#returns last concept linked in cache
        try:
            self.RefCountPandas = pd.read_csv(self.__cntCache__,sep='\t',header=0)
            last_col_name = list(self.RefCountPandas.columns)[-1]
            return last_col_name[len(self._col_name_prefix):]
        except FileNotFoundError:
            print('Cache was not found! Will start processing all input files from scratch')
            return FileNotFoundError

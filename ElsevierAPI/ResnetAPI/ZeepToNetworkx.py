import networkx as nx
from  .PathwayStudioGOQL import OQL
from  .PathwayStudioZeepAPI import DataModel
from  .NetworkxObjects import PSObject,PSRelation,len,REGULATORS,TARGETS,EFFECT
from  .ResnetGraph import ResnetGraph
import math
import time
from datetime import timedelta


class PSNetworx(DataModel):
    '''
    PSNetworx is not aware of retreived properties.\n
    Properties for retreival must be passed to all functions explicitly
    '''
    def __init__(self, url, username, password):
        super().__init__(url, username, password)
        self.id2relation = dict()  # {relID:{node_id1,node_id2,PSRelation}} needs to be - Resnet relations may not be binary
        self.Graph = ResnetGraph()
        self.ID2Children = dict()

    @staticmethod
    def _zeep2psobj(zeep_objects):
        id2entity = dict()
        if type(zeep_objects) == type(None):
            return id2entity
        for o in zeep_objects.Objects.ObjectRef:
            ps_obj = PSObject.from_zeep(o)
            obj_id = o.Id
            id2entity[obj_id] = ps_obj

        for prop in zeep_objects.Properties.ObjectProperty:
            obj_id = prop.ObjId
            if obj_id == 72057594037935395:
                print('')
            prop_id = prop.PropId
            prop_name = prop.PropName
            if type(prop.PropValues) != type(None):
                values = prop.PropValues.string
                id2entity[obj_id][prop_id] = values
                id2entity[obj_id][prop_name] = values
                try:
                    prop_display_name = prop.PropDisplayName
                    id2entity[obj_id][prop_display_name] = values
                except AttributeError:
                    continue

        return id2entity

    
    @staticmethod
    def execution_time(execution_start):
        return "{}".format(str(timedelta(seconds=time.time() - execution_start)))


    def _load_graph(self, zeep_relations, zeep_objects, add2self = True):
        new_graph = ResnetGraph()
        # loading entities and their properties
        id2entity = self._zeep2psobj(zeep_objects)
        new_graph.add_nodes(id2entity)

        new_relations = dict()
        if type(zeep_relations) != type(None):
            for rel in zeep_relations.Objects.ObjectRef:
                ps_rel = PSRelation.from_zeep(rel)
                rel_id = rel['Id']
                new_relations[rel_id] = ps_rel

            # loading relations and their properties
            for prop in zeep_relations.Properties.ObjectProperty:
                rel_id = prop['ObjId']
                prop_id = prop['PropId']
                prop_set_id = prop['PropSet']
                prop_name = prop['PropName']
                prop_display_name = prop['PropDisplayName']
                values = prop['PropValues']['string']

                if not self.IdToPropType[prop_id]['IsMultiple']:
                    new_relations[rel_id][prop_id] = values
                    new_relations[rel_id][prop_name] = values
                    new_relations[rel_id][prop_display_name] = values
                elif prop_set_id in new_relations[rel_id].PropSetToProps.keys():
                    new_relations[rel_id].PropSetToProps[prop_set_id][prop_id] = values
                    new_relations[rel_id].PropSetToProps[prop_set_id][prop_name] = values
                    new_relations[rel_id].PropSetToProps[prop_set_id][prop_display_name] = values
                else:
                    new_relations[rel_id].PropSetToProps[prop_set_id] = {prop_id: values}
                    new_relations[rel_id].PropSetToProps[prop_set_id] = {prop_name: values}
                    new_relations[rel_id].PropSetToProps[prop_set_id] = {prop_display_name: values}

            # loading connected entities from Links
            for l in zeep_relations.Links.Link:
                rel_id = l['RelationId']
                direction = l['Dir'] # 0:no arrow, 1:arrow to entity, -1:arrow from entity
                link = (l['EntityId'], direction, l[EFFECT])

                if direction == 1:
                    if len(new_relations[rel_id].Nodes) < 2:
                        new_relations[rel_id].Nodes[TARGETS] = [link]
                    else:
                        new_relations[rel_id].Nodes[TARGETS].append(link)
                else:
                    if not len(new_relations[rel_id].Nodes):
                        new_relations[rel_id].Nodes[REGULATORS] = [link]
                    else:
                        new_relations[rel_id].Nodes[REGULATORS].append(link)

            try:
                new_relations[rel_id].Nodes[TARGETS].sort(key=lambda x:x[0])
                # sorting TARGETS by id for speed
            except KeyError: pass

            try:
                new_relations[rel_id].Nodes[REGULATORS].sort(key=lambda x:x[0])
                # sorting REGULATORS by id for speed
            except KeyError: pass
            
            [new_graph.add_rel(rel) for rel in new_relations.values()]
            
        if add2self:
            self.Graph.add_graph(new_graph)
            self.id2relation.update(new_relations)

        return new_graph


    def load_graph_from_oql(self,oql_query:str,relation_props:list=[],entity_props:list=[],get_links=True,add2self=True):
        '''
        # to retrieve only IDs for speed cache specify relation_props,entity_props as empty lists
        '''
        if get_links:
            zeep_relations = self.get_data(oql_query, list(relation_props), getLinks=True)
            if type(zeep_relations) != type(None):
                obj_id_list = list(set([x['EntityId'] for x in zeep_relations.Links.Link]))
                zeep_objects = self.get_object_properties(obj_id_list, list(entity_props))
                return self._load_graph(zeep_relations, zeep_objects,add2self)
            else: return ResnetGraph()
        else:
            zeep_objects = self.get_data(oql_query, list(entity_props), getLinks=False)
            return self._load_graph(None, zeep_objects,add2self)


    def _obj_id_by_oql(self, oql_query: str):
        zeep_entities = self.get_data(oql_query, retrieve_props=['Name'], getLinks=False)
        if type(zeep_entities) != type(None):
            obj_ids = set([x['Id'] for x in zeep_entities.Objects.ObjectRef])
            return obj_ids
        else:
            return set()


    def _id2children_(self,parent_ids:list):
        '''
        Returns
        -------
        parent2children = {parent_id:[children_ids]}

        Updates
        -------
        self.ID2Children with parent2children
        '''
        parent2children = dict()
        for parent_id in parent_ids:
            query_ontology = OQL.get_childs([parent_id],['id'])
            children_ids = list(self._obj_id_by_oql(query_ontology))
            parent2children[parent_id] = children_ids

        self.ID2Children.update(parent2children)
        return parent2children


    def _get_obj_ids_by_props(self, propValues: list, search_by_properties=[], get_childs=True,
                             only_obj_types=[]):
        '''
        Loads
        -----
        self.ID2Children if get_childs == True  

        Returns
        -------
        {database_ids} for all parents and children combined
        '''
        if not search_by_properties: search_by_properties = ['Name','Alias']
        query_node = OQL.get_entities_by_props(propValues, search_by_properties, only_obj_types)
        parent_ids = self._obj_id_by_oql(query_node)

        if get_childs:
            parent2children = self._id2children_(parent_ids)
            all_ids = parent_ids
            for parent_id, children_ids in parent2children.items():
                all_ids.add(parent_id)
                all_ids.update(children_ids)
            return all_ids
        else:
            return parent_ids


    def get_children(self, parent_ids:list):
        '''
        Returns
        -------
        ids of children of parents with parent_ids

        Loads
        -----
        self.ID2Children
        '''
        child_ids = set()
        for parent_id in parent_ids:
            try:
                children = list(self.ID2Children[parent_id])
                child_ids.update(children)
            except KeyError:
                query_ontology = OQL.get_childs([parent_id],['id'])
                children = list(self._obj_id_by_oql(query_ontology))
                self.ID2Children[parent_id] = children
                child_ids.update(children)

        return list(child_ids)


    def get_children_props(self, for_psobjs:list, prop_name:str):
        """
        finds children of for_psobjs and then returns their annotation by prop_name together with for_psobjs annotation
        """
        psobjs_prop_lists = [x[prop_name] for x in for_psobjs] 
        psobjs_props = [prop for x in psobjs_prop_lists for prop in x]
        children_props = list()
        psobjs_ids = [x['Id'][0] for x in for_psobjs]
        children_ids = self.get_children(psobjs_ids) # finds children in cache for the enzyme
        if children_ids:
            children_prop_lists = [x[prop_name] for nodeid,x in self.Graph.nodes(data=True) if x['Id'][0] in children_ids]
            children_props = [prop for x in children_prop_lists for prop in x]
                  
        return list(set(psobjs_props + children_props))


    def find_drugs(self, for_targets_with_ids: list, REL_PROPS: list, ENTITY_PROPS: list):
        oql_query = OQL.get_drugs(for_targets_with_ids)
        zeep_relations = self.get_data(oql_query, REL_PROPS)
        if type(zeep_relations) != type(None):
            obj_ids = list(set([x['EntityId'] for x in zeep_relations.Links.Link]))
            zeep_objects = self.get_object_properties(obj_ids, ENTITY_PROPS)
            new_ps_relations = self._load_graph(zeep_relations, zeep_objects)
            return new_ps_relations
        else:
            return ResnetGraph()


    def find_reaxys_substances(self, ForTargetsIDlist: list, REL_PROPS: list, ENTITY_PROPS: list):
        oql_query = OQL.get_reaxys_substances(ForTargetsIDlist)
        zeep_relations = self.get_data(oql_query, REL_PROPS)
        if type(zeep_relations) != type(None):
            obj_ids = list(set([x['EntityId'] for x in zeep_relations.Links.Link]))
            zeep_objects = self.get_object_properties(obj_ids, ENTITY_PROPS)
            return self._load_graph(zeep_relations, zeep_objects)
        else:
            return ResnetGraph()


    def connect_entities(self, PropertyValues1: list, SearchByProperties1: list, EntityTypes1: list,
                         PropertyValues2: list, SearchByProperties2: list, EntityTypes2: list,
                         REL_PROPS=None, connect_by_rel_types=None, ENTITY_PROPS=None):

        if not isinstance(connect_by_rel_types,list): connect_by_rel_types = list()
        rel_props = {'Name', 'RelationNumberOfReferences'} 
        if isinstance(REL_PROPS,list): rel_props.update(REL_PROPS)
        ent_props = {'Name'}
        if isinstance(ENTITY_PROPS,list): ent_props.update(ENTITY_PROPS)

        oql_query = OQL.connect_entities(PropertyValues1, SearchByProperties1, EntityTypes1, PropertyValues2,
                                         SearchByProperties2, EntityTypes2, connect_by_rel_types)
        zeep_relations = self.get_data(oql_query, list(rel_props))
        
        if type(zeep_relations) != type(None):
            obj_ids = list(set([x['EntityId'] for x in zeep_relations.Links.Link]))
            zeep_objects = self.get_object_properties(obj_ids, list(ent_props))
            return self._load_graph(zeep_relations, zeep_objects)
        else:
            return ResnetGraph()


    def get_ppi(self, InteractorIdList:set, REL_PROPS:list, ENTITY_PROPS:list):
        splitter = list() #holds lists of ids splits 
        splitter.append(list(InteractorIdList))
        number_of_splits = int(math.log2(len(InteractorIdList)))
        ppi_keeper = ResnetGraph()
        for s in range(1, number_of_splits):
            new_splitter = list()
            half = int(len(splitter[0]) / 2)
            for split in splitter:
                uq_list1 = split[0:half]
                uq_list2 = split[half:]
                oql_query = OQL.get_ppi(set(uq_list1), set(uq_list2))
                zeep_relations = self.get_data(oql_query, REL_PROPS)
                if type(zeep_relations) != type(None):
                    obj_ids = list(set([x['EntityId'] for x in zeep_relations.Links.Link]))
                    zeep_objects = self.get_object_properties(obj_ids, ENTITY_PROPS)
                    new_ps_relations = self._load_graph(zeep_relations, zeep_objects)
                    ppi_keeper.add_graph(new_ps_relations)

                new_splitter.append(uq_list1)
                new_splitter.append(uq_list2)

            splitter = new_splitter
            s += 1

        return ppi_keeper


    def get_network(self, InteractorIdList:set, connect_by_rel_types:list=None, REL_PROPS:list=None, ENTITY_PROPS:list=None):
        splitter = list() #holds lists of ids splits 
        splitter.append(list(InteractorIdList))
        number_of_splits = int(math.log2(len(InteractorIdList)))
        print('Will load network of %d nodes in %d iterations' % (len(InteractorIdList),number_of_splits))
        accumulate_network = ResnetGraph()
        for s in range(1, number_of_splits):
            new_splitter = list()
            half = int(len(splitter[0]) / 2)
            iter_start = time.time()
            ids1 = set()
            ids2 = set()
            for split in splitter:
                uq_list1 = split[:half]
                uq_list2 = split[half:]
                ids1.update(uq_list1)
                ids2.update(uq_list2)
                new_splitter.append(uq_list1)
                new_splitter.append(uq_list2)

            oql_query = 'SELECT Relation WHERE NeighborOf (SELECT Entity WHERE id = ({ids1})) AND NeighborOf (SELECT Entity WHERE id = ({ids2}))'
            oql_query = oql_query.format(ids1=','.join(map(str,ids1)), ids2=','.join(map(str,ids2)))
            if isinstance(connect_by_rel_types,list):
                oql_query = oql_query +  'AND objectType = ({rel_types})'.format(rel_types=','.join(connect_by_rel_types))

            network_iter = self.load_graph_from_oql(oql_query, REL_PROPS,ENTITY_PROPS)
            accumulate_network = nx.compose(accumulate_network, network_iter)

            splitter = new_splitter
            executiontime = self.execution_time(iter_start)
            print('Iteration %d out of %d was completed in %s' % (s,number_of_splits, executiontime))
            s += 1
        return accumulate_network


    def _iterate_oql(self,oql_query:str,id_set:set,REL_PROPS=list(),ENTITY_PROPS=list(),add2self=True,get_links=True):
        # oql_query MUST contain string placeholder called {ids} 
        entire_graph = ResnetGraph()
        id_list = list(id_set)
        step = 1000
        for i in range(0,len(id_list), step):
            ids = id_list[i:i+step]
            oql_query_with_ids = oql_query.format(ids=','.join(map(str,ids)))
            iter_graph = self.load_graph_from_oql(oql_query_with_ids,REL_PROPS,ENTITY_PROPS,add2self=add2self,get_links=get_links)
            entire_graph.add_graph(iter_graph)
        return entire_graph


    def _iterate_oql_s(self,oql_query:str,prop_set:set,REL_PROPS=list(),ENTITY_PROPS=list(),add2self=True,get_links=True):
        '''
            # oql_query MUST contain string placeholder called {props} 
        '''
        entire_graph = ResnetGraph()
        prop_list = list(prop_set)
        step = 1000
        for i in range(0,len(prop_list), step):
            props = prop_list[i:i+step]
            oql_query_with_props = oql_query.format(props=OQL.join_with_quotes(props))
            iter_graph = self.load_graph_from_oql(oql_query_with_props,REL_PROPS,ENTITY_PROPS,add2self=add2self,get_links=get_links)
            entire_graph.add_graph(iter_graph)
        return entire_graph


    def _iterate_oql2(self, oql_query:str, id_set1:set, id_set2:set, REL_PROPS=list(), ENTITY_PROPS=list(),add2self=True,get_links=True):
        # oql_query MUST contain 2 string placeholders called {ids1} and {ids2}
        entire_graph = ResnetGraph()
        id_list1 = list(id_set1)
        id_list2 = list(id_set2)
        step = 1000
        number_of_iterations = math.ceil(len(id_set1)/step) * math.ceil(len(id_set2)/step)
        print('\nConnecting %d with %d entities' % (len(id_set1), len(id_set2)))
        if number_of_iterations > 2:
            print('Query will be executed in %d iterations' % number_of_iterations)

        iteration_counter = 1
        start  = time.time()
        for i1 in range(0,len(id_list1), step):
            ids1 = id_list1[i1:i1+step]
            for i2 in range(0,len(id_list2), step):
                ids2 = id_list2[i2:i2+step]
                oql_query_with_ids = oql_query.format(ids1=','.join(map(str,ids1)),ids2=','.join(map(str,ids2)))
                iter_graph = self.load_graph_from_oql(oql_query_with_ids,REL_PROPS,ENTITY_PROPS,add2self=add2self,get_links=get_links)
                entire_graph.add_graph(iter_graph)
                if number_of_iterations > 2:
                    print('Iteration %d out of %d performed in %s' % 
                        (iteration_counter,number_of_iterations,self.execution_time(start)))
                iteration_counter +=1
        return entire_graph
    
    
    def get_pathway_members(self, pathway_ids: list, search_pathways_by=None, only_entities=None,
                               with_properties=None):
        """
        returns id2psobj {id:PSObject} of entities from pathways found by 'search_pathways_by' or from 'pathway_ids'
        """

        if with_properties is None:
            with_properties = ['Name', 'Alias']
        if only_entities is None:
            only_entities = []

        if search_pathways_by is None:
            oql_query = 'SELECT Entity WHERE MemberOf (SELECT Network WHERE id = (' + ','.join(list(map(str,pathway_ids))) + '))'
        else:
            property_names, values = OQL.get_search_strings(search_pathways_by, pathway_ids)
            oql_query = 'SELECT Entity WHERE MemberOf (SELECT Network WHERE (' + property_names + ') = (' + values + '))'

        if len(only_entities) > 0:
            filter_prop_name, filter_values = OQL.get_search_strings(with_properties, only_entities)
            if len(with_properties) > 1:
                oql_query = oql_query + ' AND (' + filter_prop_name + ') = (' + filter_values + ')'
            else:
                oql_query = oql_query + ' AND ' + filter_prop_name + ' = (' + filter_values + ')'

        zeep_objects = self.get_data(oql_query, retrieve_props=['Name'], getLinks=False)
        return self._zeep2psobj(zeep_objects)


    def find_targets_in_pathways(self, DrugProps: list, DrugSearchPropertyNames: list, PathwayNames: list,
                                 relation_types=None, target_types=None):
        if target_types is None:
            target_types = ['Protein']
        if relation_types is None:
            relation_types = ['DirectRegulation']

        rel_props = ['Name', 'RelationNumberOfReferences', 'DOI', 'PMID', 'Source']
        rel_types = ','.join(relation_types)
        target_types = ','.join(target_types)
        drug_prop_names, drug_prop_values = OQL.get_search_strings(DrugSearchPropertyNames, DrugProps)
        drug_query = 'Select Entity WHERE (' + drug_prop_names + ') = (' + drug_prop_values + ')'
        pathway_query = 'SELECT Network WHERE Name = (' + OQL.join_with_quotes( PathwayNames) + ')'
        oql_query = 'SELECT Relation WHERE objectType = ({RelTypes}) AND NeighborOf upstream (SELECT Entity WHERE MemberOf ({pathways}) AND objectType = ({Target_Types})) AND NeighborOf downstream ({drug})'

        direct_regulations = self.load_graph_from_oql(
            oql_query.format(RelTypes=rel_types, pathways=pathway_query, Target_Types=target_types, drug=drug_query),
            relation_props=rel_props)
        oql_query = 'SELECT Relation WHERE objectType = Binding AND NeighborOf (SELECT Entity WHERE MemberOf ({pathways}) AND objectType = ({Target_Types})) AND NeighborOf ({drug})'
        bindings = self.load_graph_from_oql(
            oql_query.format(pathways=pathway_query, Target_Types=target_types, drug=drug_query),
            relation_props=rel_props)

        if isinstance(direct_regulations, ResnetGraph):
            if isinstance(bindings, ResnetGraph):
                return nx.compose(direct_regulations, bindings)
            else:
                return direct_regulations
        else:
            if isinstance(bindings, ResnetGraph):
                return bindings
            else:
                return ResnetGraph()

    def find_drug_toxicities(self, DrugIds: list, DrugSearchPropertyNames: list, min_ref_count=0,
                             relation_properties=None, entity_properties=None):
        if entity_properties is None:
            entity_properties = ['Name']
        if relation_properties is None:
            relation_properties = ['Name', 'RelationNumberOfReferences']

        drug_prop_names, drug_prop_values = OQL.get_search_strings(DrugSearchPropertyNames, DrugIds)
        drug_query = 'Select Entity WHERE (' + drug_prop_names + ') = (' + drug_prop_values + ')'

        oql_query = 'SELECT Relation WHERE objectType = Regulation AND Effect = positive AND ' \
                    'RelationNumberOfReferences >= {minref} AND NeighborOf upstream (SELECT Entity WHERE objectType = ' \
                    'Disease) AND NeighborOf downstream ({drug})'
        oql_query = oql_query.format(minref=str(min_ref_count), drug=drug_query)
        diseases = self.load_graph_from_oql(oql_query, relation_properties, entity_properties)

        clinical_params = 'SELECT Entity WHERE InOntology (SELECT Annotation WHERE Ontology=\'Pathway Studio Ontology\' AND Relationship=\'is-a\') under (SELECT OntologicalNode WHERE Name=\'disease-related parameters\')'
        oql_query = 'SELECT Relation WHERE objectType = Regulation AND Effect = positive AND RelationNumberOfReferences >= {minref} AND NeighborOf upstream ({ClinicalParameters}) AND NeighborOf downstream ({drug})'
        oql_query = oql_query.format(minref=str(min_ref_count), drug=drug_query, ClinicalParameters=clinical_params)
        clinical_parameters = self.load_graph_from_oql(oql_query, relation_properties, entity_properties)

        if isinstance(diseases, ResnetGraph):
            if isinstance(clinical_parameters, ResnetGraph):
                return nx.compose(diseases, clinical_parameters)
            else:
                return diseases
        else:
            if isinstance(clinical_parameters, ResnetGraph):
                return clinical_parameters
            else:
                return ResnetGraph()


    def pathway_components(self,by_pathway_props:list,in_prop_type:str, 
                relprops2load:list=[], entprops2load:list=[]):
        """
        Returns
        -------
        ResnetGraph containing graphs merged from all pathways found\n
        with by_pathway_props+in_prop_type
        """
         
        rel_query = 'SELECT Relation WHERE MemberOf (SELECT Network WHERE {propName} = \'{pathway}\')'
        subgraph_relation_ids = set()
        for pathway_prop in by_pathway_props:
            loaded_node_ids = set(self.Graph.nodes())
            loaded_relation_ids = set(self.Graph._relations_ids())

            rel_q = rel_query.format(propName=in_prop_type,pathway=pathway_prop)
            pathway_id_only_graph = self.load_graph_from_oql(rel_q,[],[],get_links=True,add2self=False)
            pathway_relation_ids = set(pathway_id_only_graph._relations_ids())
            subgraph_relation_ids.update(pathway_relation_ids)
            new_relation_ids = pathway_relation_ids.difference(loaded_relation_ids)
    
            if new_relation_ids:
                oql_query = OQL.get_relations_by_props(list(new_relation_ids),['id'])
                new_relation_graph = self.load_graph_from_oql(oql_query,relprops2load)
                # new_relation_graph is now fully loaded with desired rel_props
                # still need desired props for new nodes that did not exist in self.Graph
                new_nodes_ids = set(new_relation_graph.nodes()).difference(loaded_node_ids)
                if new_nodes_ids:
                    entity_query = OQL.get_objects(list(new_nodes_ids))
                    self.load_graph_from_oql(entity_query,[],entprops2load,get_links=False)
        
        rels4subgraph = [rel for i,rel in self.id2relation.items() if i in subgraph_relation_ids]
        return_subgraph = self.Graph.subgraph_by_rels(rels4subgraph)
        return return_subgraph


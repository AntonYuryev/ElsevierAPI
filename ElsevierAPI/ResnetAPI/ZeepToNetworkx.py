import networkx as nx
import ElsevierAPI.ResnetAPI.PathwayStudioGOQL as OQL
from ElsevierAPI.ResnetAPI.PathwayStudioZeepAPI import DataModel
from ElsevierAPI.ResnetAPI.NetworkxObjects import PSObject,PSRelation
from ElsevierAPI.ResnetAPI.ResnetGraph import ResnetGraph
import math
import time
from datetime import timedelta

class PSNetworx(DataModel):
    def __init__(self, url, username, password):
        super().__init__(url, username, password)
        self.IDtoRelation = dict()  # {relID:PSRelation} needs to be - Resnet relations may not be binary
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
            prop_id = prop.PropId
            prop_name = prop.PropName
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
    def link_id(t:tuple):
        return t[0]

    def _load_graph(self, zeep_relations, zeep_objects, add2self = True):
        new_graph = ResnetGraph()
        # loading entities and their properties
        id2entity = self._zeep2psobj(zeep_objects)
        new_graph.add_nodes_from([(k, v.items()) for k, v in id2entity.items()])

        if type(zeep_relations) != type(None):
            new_relations = dict()
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
                direction = l['Dir']
                link = (l['EntityId'], direction, l['Effect'])

                if direction == 1:
                    if len(new_relations[rel_id].Nodes) < 2:
                        new_relations[rel_id].Nodes['Targets'] = [link]
                    else:
                        new_relations[rel_id].Nodes['Targets'].append(link)
                else:
                    if len(new_relations[rel_id].Nodes) < 1:
                        new_relations[rel_id].Nodes['Regulators'] = [link]
                    else:
                        new_relations[rel_id].Nodes['Regulators'].append(link)

            try:
                new_relations[rel_id].Nodes['Targets'].sort(key=self.link_id)
            except KeyError: pass

            try:
                new_relations[rel_id].Nodes['Regulators'].sort(key=self.link_id)
            except KeyError: pass
            
            for rel in new_relations.values():
                regulator_target = rel.get_regulators_targets()
                for pair in regulator_target:
                    try:
                        ref_count = rel['RelationNumberOfReferences'][0]
                    except KeyError: ref_count = 0
                    new_graph.add_edge(pair[0], pair[1], relation=rel, weight=float(ref_count))
                    # print (newGraph.get_edge_data(pair[0], pair[1]))

            self.IDtoRelation.update(new_relations)  # must be kept since Resnet relation may not be binary

        if add2self:
            self.Graph = nx.compose(self.Graph,new_graph)

        return new_graph

    def load_graph_from_oql(self, oql_query: str, relation_props: list=None, entity_props: list=None, get_links=True, add2self=True):
        
        if isinstance(entity_props,(list, set)):
            entity_props = set(['Name']+list(entity_props)) if entity_props else set()
            #specify empty list explicitly to avoid retreival of any properties
            #can be used to retrieve only IDs for speed
        else: entity_props = {'Name'}
        
        if get_links:
            if isinstance(relation_props,(set,list)):
                relation_props = set(relation_props+['Name','RelationNumberOfReferences']) if relation_props else set()
                #specify empty list explicitly to avoid retreival of any properties
                #can be used to retrieve only IDs for speed
            else: relation_props = {'Name','RelationNumberOfReferences'}
            
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

    def _get_obj_ids_by_props(self, propValues: list, search_by_properties=None, get_childs=True,
                             only_obj_types=[]):
        if search_by_properties is None: search_by_properties = ['Name','Alias']
        
        query_node = OQL.get_entities_by_props(propValues, search_by_properties, only_obj_types)
        target_ids = self._obj_id_by_oql(query_node)
        
        if get_childs:
            child_ids = list()
            for i in target_ids:
                try:
                    child_ids = child_ids + self.ID2Children[i]
                except KeyError:
                    query_ontology = OQL.get_childs([i],['id'])
                    children = list(self._obj_id_by_oql(query_ontology))
                    self.ID2Children[i] = children
                    child_ids = child_ids + children

            target_ids.update(child_ids)
            
        return target_ids


    def get_children(self, parent_ids:list):
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
                    ppi_keeper = nx.compose(ppi_keeper, new_ps_relations)

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
            execution_time = "{}".format(str(timedelta(seconds=time.time() - iter_start)))
            print('Iteration %d out of %d was completed in %s' % (s,number_of_splits, execution_time))
            s += 1
        return accumulate_network

    def _iterate_oql(self, oql_query:str, id_set:set, REL_PROPS:list=None, ENTITY_PROPS:list=None):
        # oql_query MUST contain string placeholder called {ids} 
        entire_graph = ResnetGraph()
        id_list = list(id_set)
        step = 500
        for i in range(0,len(id_list), step):
            ids = id_list[i:i+step]
            oql_query_with_ids = oql_query.format(ids=','.join(map(str,ids)))
            iter_graph = self.load_graph_from_oql(oql_query_with_ids,REL_PROPS,ENTITY_PROPS)
            entire_graph = nx.compose(entire_graph,iter_graph)
        return entire_graph

    def _iterate_oql2(self, oql_query:str, id_set1:set, id_set2:set, REL_PROPS:list=None, ENTITY_PROPS:list=None):
        # oql_query MUST contain 2 string placeholders called {ids1} and {ids2}
        entire_graph = ResnetGraph()
        id_list1 = list(id_set1)
        id_list2 = list(id_set2)
        step = 500
        for i1 in range(0,len(id_list1), step):
            ids1 = id_list1[i1:i1+step]
            for i2 in range(0,len(id_list2), step):
                ids2 = id_list2[i2:i2+step]
                oql_query_with_ids = oql_query.format(ids1=','.join(map(str,ids1)),ids2=','.join(map(str,ids2)))
                iter_graph = self.load_graph_from_oql(oql_query_with_ids,REL_PROPS,ENTITY_PROPS)
                entire_graph = nx.compose(entire_graph,iter_graph)
        return entire_graph
    
    def get_pathway_member_ids(self, PathwayIds: list, search_pathways_by=None, only_entities=None,
                               with_properties=None):
        if with_properties is None:
            with_properties = ['Name', 'Alias']
        if only_entities is None:
            only_entities = []
        if search_pathways_by is None:
            search_pathways_by = ['id']

        if search_pathways_by[0] in ['id', 'Id', 'ID']:
            oql_query = 'SELECT Entity WHERE MemberOf (SELECT Network WHERE id = (' + ','.join(
                [str(i) for i in PathwayIds]) + '))'
        else:
            property_names, values = OQL.get_search_strings(search_pathways_by, PathwayIds)
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
        pathway_query = 'SELECT Network WHERE Name = (' + OQL.join_with_quotes(',', PathwayNames) + ')'
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


    def get_pathway_components(self, prop_vals: list, search_by_property: str, retrieve_rel_properties:set=set(), retrieve_ent_properties:set=set()):

        ent_props = retrieve_ent_properties
        ent_props.add('Name')

        rel_props = retrieve_rel_properties
        rel_props.update(['Name','RelationNumberOfReferences'])
         
        rel_query = 'SELECT Relation WHERE MemberOf (SELECT Network WHERE {propName} = {pathway})'
        accumulate_pathways = ResnetGraph()
        for pathway_prop in prop_vals:

            rel_q = rel_query.format(propName=search_by_property,pathway=pathway_prop)
            pathway_relations_ids_only = self.load_graph_from_oql(rel_q, relation_props=[],get_links=True,add2self=False)
            new_relations = pathway_relations_ids_only.subtract(self.Graph)
            exist_relations = self.Graph.intersect(pathway_relations_ids_only)
            # must subtract from self.Graph to keep references
    
            if new_relations:
                new_relation_ids = new_relations._relations_ids()
                oql_query = OQL.get_relations_by_props(list(new_relation_ids),['id'])
                new_relations = self.load_graph_from_oql(oql_query,relation_props=list(rel_props))

            pathway_relations = nx.compose(exist_relations,new_relations)
            
            exist_node_ids = set(exist_relations)
            if exist_node_ids:
                existing_nodes_str = ','.join(map(str,exist_node_ids))
                ent_q = 'SELECT Entity WHERE MemberOf (SELECT Network WHERE {propName} = {pathway}) AND NOT (id = ({exist_nodes}))'
                ent_q = ent_q.format(propName=search_by_property,pathway=pathway_prop,exist_nodes=existing_nodes_str)
            else:
                ent_q = 'SELECT Entity WHERE MemberOf (SELECT Network WHERE {propName} = {pathway})'
                ent_q = ent_q.format(propName=search_by_property,pathway=pathway_prop)

            pathway_nodes = self.load_graph_from_oql(ent_q,entity_props=list(ent_props), get_links=False)
            pathway_graph = nx.compose(pathway_nodes, pathway_relations)
            accumulate_pathways = nx.compose(accumulate_pathways,pathway_graph)

        return accumulate_pathways


    def to_rnef(self, in_graph=None,add_rel_props:dict=None,add_pathway_props:dict=None):
        # add_rel_props structure {PropName:[PropValues]}
        if not isinstance(in_graph,ResnetGraph): in_graph=self.Graph
        all_rnef_props = set(self.RNEFnameToPropType.keys()) #set(list(self.RNEFnameToPropType.keys())+REF_PROPS)
        return in_graph.to_rnef(list(all_rnef_props),add_rel_props,add_pathway_props)

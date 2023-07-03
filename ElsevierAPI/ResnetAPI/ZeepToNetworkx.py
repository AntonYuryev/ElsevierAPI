import networkx as nx
from  .PathwayStudioGOQL import OQL
from  .PathwayStudioZeepAPI import DataModel
from  .NetworkxObjects import PSObject,PSRelation,len,REGULATORS,TARGETS,EFFECT
from  .ResnetGraph import ResnetGraph,REFCOUNT
import math
import time
from datetime import timedelta


class PSNetworx(DataModel):
    '''
    PSNetworx is not aware of retreived properties.\n
    Property names for retrieval must be passed to all functions explicitly
    '''
    def __init__(self, *args,**kwargs):
        '''
        Input
        -----
        APIconfig = args[0]
        no_mess - default False, if True your script becomes more verbose\n
        connect2server - default True, set to False to run script using data in __pscache__ files instead of database
        '''
        my_kwargs = dict(kwargs)
        super().__init__(*args,**my_kwargs)
        self.dbid2relation = dict()  # {relID:{node_id1,node_id2,PSRelation}} needs to be - Resnet relations may not be binary
        self.Graph = ResnetGraph()


    @staticmethod
    def _zeep2psobj(zeep_objects):
        '''
        Returns
        -------
        id2entity = {db_id:PSObject}
        '''
        dbid2entity = dict()
        if type(zeep_objects) != type(None):
            for o in zeep_objects.Objects.ObjectRef:
                ps_obj = PSObject.from_zeep(o)
                dbid2entity[ps_obj.dbid()] = ps_obj

            for prop in zeep_objects.Properties.ObjectProperty:
                dbid = prop.ObjId
                #prop_id = prop.PropId
                prop_name = prop.PropName
                if type(prop.PropValues) != type(None):
                    values = prop.PropValues.string
                    #id2entity[db_id][str(prop_id] = values
                    dbid2entity[dbid][prop_name] = values
                    try:
                        prop_display_name = prop.PropDisplayName
                        dbid2entity[dbid][prop_display_name] = values
                    except AttributeError:
                        continue

        return dbid2entity

    
    @staticmethod
    def execution_time(execution_start,remaining_iterations=None,number_of_iterations=None):
        '''
        Input
        -----
        if "number_of_iterations" is supplied assumes that "execution_start" is global start
        otherwise assumes "execution_start" is the start of the current iteration if "remaining_iterations" is supplied
        '''
        delta = time.time() - execution_start
        if type(number_of_iterations) != type(None):
            passed_iterations = number_of_iterations-remaining_iterations
            if passed_iterations:
                remaining_time = delta*float(remaining_iterations/passed_iterations)
            else:
                remaining_time = 0
            return "{}".format(str(timedelta(seconds=delta))), "{}".format(str(timedelta(seconds=remaining_time)))
        elif type(remaining_iterations) != type(None):
                remaining_time = delta*remaining_iterations
                return "{}".format(str(timedelta(seconds=delta))), "{}".format(str(timedelta(seconds=remaining_time)))
        else:
            return "{}".format(str(timedelta(seconds=delta)))
    

    def __psrel2dict(self,rels:dict):
        '''
        Input
        -----
        rels = {dbid:PSRelation}
        '''
        for dbid, rel in rels.items():
            try:
                self.dbid2relation[dbid].merge_rel(rel)
            except:
                self.dbid2relation[dbid] = rel
            
        
    def _load_graph(self, zeep_relations, zeep_objects, add2self=True, merge_data=False):
        new_graph = ResnetGraph()
        # loading entities and their properties
        id2psobj = self._zeep2psobj(zeep_objects)
        new_graph.add_nodes_from([(n.uid(),n.items()) for n in id2psobj.values()])

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
                direction = l['Dir'] # 0:no arrow, 1:arrow to entity from relation, -1:arrow from entity to relation
                graph_node = id2psobj[l['EntityId']]
                # link = (l['EntityId'], direction, l[EFFECT])
                link = (graph_node.uid(), direction, l[EFFECT])

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
            
            [new_graph.add_rel(rel,merge=False) for rel in new_relations.values()]
            
        if add2self:
            if merge_data:
                self.Graph.add_graph(new_graph,merge=True)
                self.__psrel2dict(new_relations)
            else:
                self.Graph = self.Graph.compose(new_graph)
                self.dbid2relation.update(new_relations)

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
                return self._load_graph(zeep_relations,zeep_objects,add2self)
            else: return ResnetGraph()
        else:
            zeep_objects = self.get_data(oql_query, list(entity_props), getLinks=False)
            return self._load_graph(None, zeep_objects,add2self)


    def _db_id_by_oql(self, oql_query: str):
        zeep_entities = self.get_data(oql_query, retrieve_props=['Name'], getLinks=False)
        if type(zeep_entities) != type(None):
            obj_ids = set([x['Id'] for x in zeep_entities.Objects.ObjectRef])
            return obj_ids
        else:
            return set()


    def find_drugs(self, for_targets_with_ids: list, REL_PROPS: list, ENTITY_PROPS: list):
        oql_query = OQL.drugs4(for_targets_with_ids)
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
        rel_props = {'Name', REFCOUNT} 
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

    '''
    def get_ppi(self, interactors_dbids:set, REL_PROPS:list, ENTITY_PROPS:list,add2self=True):
        """
        PPI relation types: Binding, DirectRegulation, ProtModification
        """
        splitter = list() #holds lists of ids splits 
        splitter.append(list(interactors_dbids))
        number_of_splits = int(math.log2(len(interactors_dbids)))
        ppi_keeper = ResnetGraph()
        for s in range(1, number_of_splits):
            new_splitter = list()
            half = int(len(splitter[0]) / 2)
            #futures = list()
            for split in splitter:
                uq_list1 = split[0:half]
                uq_list2 = split[half:]
                
                oql_query = OQL.get_ppi(set(uq_list1), set(uq_list2))
                new_ps_relations = self.load_graph_from_oql(oql_query,REL_PROPS,ENTITY_PROPS,add2self=add2self)
                ppi_keeper = ppi_keeper.compose(new_ps_relations)
                new_splitter.append(uq_list1)
                new_splitter.append(uq_list2)

            splitter = new_splitter
            s += 1

        return ppi_keeper
        '''

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
            oql_query = oql_query.format(ids1=','.join(list(map(str,ids1))), ids2=','.join(list(map(str,ids2))))
            if isinstance(connect_by_rel_types,list):
                oql_query = oql_query +  'AND objectType = ({rel_types})'.format(rel_types=','.join(connect_by_rel_types))

            network_iter = self.load_graph_from_oql(oql_query, REL_PROPS,ENTITY_PROPS)
            accumulate_network = nx.compose(accumulate_network, network_iter)

            splitter = new_splitter
            executiontime = self.execution_time(iter_start)
            print('Iteration %d out of %d was completed in %s' % (s,number_of_splits, executiontime))
            s += 1
        return accumulate_network
    
    
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

        rel_props = ['Name', REFCOUNT, 'DOI', 'PMID', 'Source']
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
            relation_properties = ['Name', REFCOUNT]

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
        ResnetGraph containing graphs merged from all pathways found with "by_pathway_props" "in_prop_type"
        """
         
        rel_query = 'SELECT Relation WHERE MemberOf (SELECT Network WHERE {propName} = \'{pathway}\')'
        subgraph_relation_dbids = set()
        loaded_node_ids = set(self.Graph.dbids4nodes())
        loaded_relation_ids = set(self.Graph.dbids4nodes())
        for pathway_prop in by_pathway_props:
            rel_q = rel_query.format(propName=in_prop_type,pathway=pathway_prop)
            pathway_id_only_graph = self.load_graph_from_oql(rel_q,[],[],get_links=True,add2self=False)
            pathway_relation_dbids = set(pathway_id_only_graph.relation_dbids())
            subgraph_relation_dbids.update(pathway_relation_dbids)
            new_relation_ids = pathway_relation_dbids.difference(loaded_relation_ids)
    
            if new_relation_ids:
                oql_query = OQL.get_relations_by_props(list(new_relation_ids),['id'])
                new_relation_graph = self.load_graph_from_oql(oql_query,relprops2load)
                # new_relation_graph is now fully loaded with desired rel_props
                # still need desired props for new nodes that did not exist in self.Graph
                new_nodes_ids = set(new_relation_graph.dbids4nodes()).difference(loaded_node_ids)
                if new_nodes_ids:
                    entity_query = OQL.get_objects(list(new_nodes_ids))
                    self.load_graph_from_oql(entity_query,[],entprops2load,get_links=False)
        
        my_dbid2relation = dict(self.dbid2relation) 
        # need to make a copy of self.dbid2relation because it is mutating during multithreaded retreival
        rels4subgraph = [rel for i,rel in my_dbid2relation.items() if i in subgraph_relation_dbids]
        return_subgraph = self.Graph.subgraph_by_rels(rels4subgraph)
        return return_subgraph
    

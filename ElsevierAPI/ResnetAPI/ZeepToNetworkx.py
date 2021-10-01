import networkx as nx
import ElsevierAPI.ResnetAPI.PathwayStudioGOQL as OQL
from ElsevierAPI.ResnetAPI.PathwayStudioZeepAPI import DataModel
from ElsevierAPI.ResnetAPI.NetworkxObjects import PSObject, PSRelation, REF_PROPS, REF_ID_TYPES
from ElsevierAPI.ResnetAPI.ResnetGraph import ResnetGraph
from xml.dom import minidom
import xml.etree.ElementTree as et

REL_PROPS = ['Effect','Mechanism']

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

    def _load_graph(self, zeep_relations, zeep_objects):
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

            for rel in new_relations.values():
                regulator_target = rel.get_regulators_targets()
                for pair in regulator_target:
                    ref_count = rel['RelationNumberOfReferences'][0]
                    new_graph.add_edge(pair[0], pair[1], relation=rel, weight=float(ref_count))
                    # print (newGraph.get_edge_data(pair[0], pair[1]))

            self.IDtoRelation.update(new_relations)  # must be kept since Resnet relation may not be binary

        self.Graph = nx.compose(new_graph, self.Graph)
        return new_graph

    def load_graph_from_oql(self, oql_query: str, relation_props: list=None, entity_props: list=None, get_links=True):
        entity_props = set(['Name']+list(entity_props)) if isinstance(entity_props,(list, set)) else {'Name'}
        if get_links:
            if isinstance(relation_props,list):
                relation_props = set(relation_props+['Name','RelationNumberOfReferences'])
            else: relation_props = {'Name','RelationNumberOfReferences'}
            
            zeep_relations = self.get_data(oql_query, list(relation_props), getLinks=True)
            if type(zeep_relations) != type(None):
                obj_id_list = list(set([x['EntityId'] for x in zeep_relations.Links.Link]))
                zeep_objects = self.get_object_properties(obj_id_list, list(entity_props))
                return self._load_graph(zeep_relations, zeep_objects)
            else: return ResnetGraph()
        else:
            zeep_objects = self.get_data(oql_query, list(entity_props), getLinks=False)
            return self._load_graph(None, zeep_objects)

    def _obj_id_by_oql(self, oql_query: str):
        zeep_entities = self.get_data(oql_query, retrieve_props=['Name'], getLinks=False)
        if type(zeep_entities) != type(None):
            obj_ids = set([x['Id'] for x in zeep_entities.Objects.ObjectRef])
            return obj_ids
        else:
            return set()

    def _get_obj_ids_by_props(self, propValues: list, search_by_properties=None, get_childs=True,
                             only_obj_types=None):

        if only_obj_types is None: only_object_types = [] 
        if search_by_properties is None: search_by_properties = ['Name','Alias']
        
        query_node = OQL.get_entities_by_props(propValues, search_by_properties, only_object_types)
        target_ids = self._obj_id_by_oql(query_node)
        
        if get_childs:
            child_ids = list()
            for i in target_ids:
                try:
                    child_ids = child_ids + self.ID2Children[i]
                except KeyError:
                    query_ontology = OQL.get_childs([i],['id'])
                    self.ID2Children[i] = list(self._obj_id_by_oql(query_ontology))
                    child_ids = child_ids + self.ID2Children[i]

            target_ids.update(child_ids)
            
        return target_ids


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
        import math
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


    def get_objects_from_folders(self, FolderIds: list, property_names=None, with_layout=False):
        if property_names is None: property_names = ['Name']
        if not hasattr(self,'id2folders'): self.id2folders = self.load_folder_tree()
        if not hasattr(self,'id2pathways'): self.id2pathways = dict()
        if not hasattr(self,'id2groups'): self.id2groups = dict()
        id2objects = dict()
        for f in FolderIds:
            folder_name = self.id2folders[f][0]['Name']
            zeep_objects = self.get_folder_objects_props(f, property_names)
            id2objs = self._zeep2psobj(zeep_objects)
            for Id, psObj in id2objs.items():
                if psObj['ObjTypeName'][0] == 'Pathway':
                    try:
                        self.id2pathways[Id].add_unique_property('Folders', folder_name)
                    except KeyError:
                            psObj['Folders'] = [folder_name]
                            self.id2pathways[Id] = psObj
                    if with_layout:
                        psObj['layout'] = self.get_layout(Id)
                if psObj['ObjTypeName'][0] == 'Group':
                    try:
                        self.id2groups[Id].add_unique_property('Folders', folder_name)
                    except KeyError:
                            psObj['Folders'] = [folder_name]
                            self.id2groups[Id] = psObj

            id2objects.update(id2objs)
        return id2objects


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


    def get_pathway_components(self, prop_vals: list, search_by_property: str, retrieve_rel_properties:list=None, retrieve_ent_properties:list=None):
        if not isinstance(retrieve_ent_properties,list): retrieve_ent_properties = ['Name']
        else: retrieve_ent_properties = list(set(retrieve_ent_properties.append('Name')))

        if not isinstance(retrieve_rel_properties,list): retrieve_rel_properties = ['Name','RelationNumberOfReferences']
        else: retrieve_rel_properties = list(set(retrieve_rel_properties + ['Name','RelationNumberOfReferences']))
         
        ent_query = 'SELECT Entity WHERE MemberOf (SELECT Network WHERE {propName} = {pathway})'
        rel_query = 'SELECT Relation WHERE MemberOf (SELECT Network WHERE {propName} = {pathway})'
        
        accumulate_pathways = ResnetGraph()
        for pathway_prop in prop_vals:
            ent_q = ent_query.format(propName=search_by_property,pathway=pathway_prop)
            pathway_nodes = self.load_graph_from_oql(ent_q,entity_props=retrieve_ent_properties, get_links=False)

            rel_q = rel_query.format(propName=search_by_property,pathway=pathway_prop)
            pathway_relations = self.load_graph_from_oql(rel_q,relation_props=retrieve_rel_properties, get_links=True)

            pathway_graph = nx.compose(pathway_nodes, pathway_relations)
            accumulate_pathways = nx.compose(pathway_graph, accumulate_pathways)

        return accumulate_pathways

    def to_rnef(self, in_graph=None,add_rel_props:dict=None,add_pathway_props:dict=None):
        # add_rel_props structure {PropName:[PropValues]}
        if not isinstance(in_graph,ResnetGraph): in_graph=self.Graph
        all_rnef_props = set(list(self.RNEFnameToPropType.keys())+REF_PROPS)
        return in_graph.to_rnef(list(all_rnef_props),add_rel_props,add_pathway_props)

    def get_all_pathways(self, property_names=None):
        if property_names is None: property_names = ['Name']
        print('retrieving identifiers of all pathways from database')

        if (len(self.id2folders)) == 0: self.load_folder_tree()

        self.id2pathways = dict()
        urn2pathway = dict()
        for folderList in self.id2folders.values():
            for folder in folderList:
                zeep_objects = self.get_folder_objects_props(folder['Id'], property_names)
                ps_objects = self._zeep2psobj(zeep_objects)
                for Id, psObj in ps_objects.items():
                    if psObj['ObjTypeName'][0] == 'Pathway':
                        try:
                            self.id2pathways[Id].add_unique_property('Folders', folder['Name'])
                        except KeyError:
                            psObj['Folders'] = [folder['Name']]
                            self.id2pathways[Id] = psObj
                            urn2pathway[psObj['URN'][0]] = psObj

        print('Found %d pathways in the database' % (len(self.id2pathways)))
        return urn2pathway


    def get_pathway(self, pathwayId,path_urn=None,path_name=None,rel_props:list=None, ent_props:list=None,
                    xml_format='RNEF',put2folder:str=None, add_rel_props:dict=None, add_pathway_props:dict=None, as_batch=True):
    # add_rel_props, add_pathway_props structure - {PropName:[PropValues]}
        if not isinstance(rel_props,list): rel_props = REF_PROPS+REF_ID_TYPES

        if hasattr(self,'id2pathways'):
            if not isinstance(path_urn,str):
                try:
                    path_urn = self.id2pathways[pathwayId]['URN'][0]
                    path_name = self.id2pathways[pathwayId]['Name'][0]
                except KeyError:
                    print('Pathway collection does not have %s pathway with URN %s' % (path_name,path_urn))

        if not isinstance(path_urn,str):
            print('Pathway has no URN specifed!!!! ')
            path_urn = 'no_urn'
        
        if not isinstance(path_name,str):
            print('Pathway has no Name specifed!!!! ')
            path_name = 'no_name'


        pathway_graph = self.get_pathway_components([pathwayId],'id',retrieve_rel_properties=rel_props,
                                                    retrieve_ent_properties=ent_props)
        pathway_graph.count_references()

        graph_xml = self.to_rnef(pathway_graph,add_rel_props,add_pathway_props)
        import xml.etree.ElementTree as et
        rnef_xml = et.fromstring(graph_xml)
        rnef_xml.set('name', path_name)
        rnef_xml.set('urn', path_urn)
        rnef_xml.set('type', 'Pathway')

        lay_out = et.Element('attachments')
        lay_out.append(et.fromstring(self.get_layout(pathwayId)))
        rnef_xml.append(lay_out)
        
        batch_xml = et.Element('batch')
        batch_xml.insert(0,rnef_xml)
                   
        if xml_format == ['SBGN']:
            pathway_xml = et.tostring(batch_xml,encoding='utf-8',xml_declaration=True).decode("utf-8")
            pathway_xml = minidom.parseString(pathway_xml).toprettyxml(indent='   ')
            from ElsevierAPI.ResnetAPI.rnef2sbgn import rnef2sbgn_str
            pathway_xml = rnef2sbgn_str(pathway_xml, classmapfile='ElsevierAPI/ResnetAPI/rnef2sbgn_map.xml')
        else:
            if isinstance(put2folder,str):
                resnet = et.Element('resnet')
                xml_nodes = et.SubElement(resnet, 'nodes')
                folder_local_id = 'F0'
                xml_node_folder = et.SubElement(xml_nodes, 'node', {'local_id':folder_local_id, 'urn': 'urn:agi-folder:xxxxx_yyyyy_zzzzz'})
                et.SubElement(xml_node_folder, 'attr', {'name': 'NodeType', 'value': 'Folder'})
                et.SubElement(xml_node_folder, 'attr', {'name': 'Name', 'value': put2folder})
                pathway_local_id = 'P0'
                xml_node_pathway = et.SubElement(xml_nodes, 'node', {'local_id':pathway_local_id, 'urn': path_urn})
                et.SubElement(xml_node_pathway, 'attr', {'name': 'NodeType', 'value': 'Pathway'})
                xml_controls = et.SubElement(resnet, 'controls')
                xml_control = et.SubElement(xml_controls, 'control', {'local_id':'CFE1'})
                et.SubElement(xml_control, 'attr', {'name':'ControlType', 'value':'MemberOf'})
                et.SubElement(xml_control, 'link', {'type':'in', 'ref':pathway_local_id})
                et.SubElement(xml_control, 'link', {'type':'out', 'ref':folder_local_id})
                batch_xml.append(resnet)
            
            if as_batch:
                pathway_xml = et.tostring(batch_xml,encoding='utf-8',xml_declaration=True).decode("utf-8")
                pathway_xml = minidom.parseString(pathway_xml).toprettyxml(indent='   ')
            else:
                pathway_xml = et.tostring(rnef_xml,encoding='utf-8',xml_declaration=True).decode("utf-8")
                pathway_xml = str(minidom.parseString(pathway_xml).toprettyxml(indent='   '))
                pathway_xml = pathway_xml[pathway_xml.find('\n')+1:]
                #minidom does not work without xml_declaration

        print('\"%s\" pathway downloaded: %d nodes, %d edges supported by %d references' % 
            (path_name, pathway_graph.number_of_nodes(),pathway_graph.number_of_edges(),pathway_graph.size(weight="weight")))

        return pathway_graph, str(pathway_xml)


    def get_group(self, group_id,group_urn=None,group_name=None, ent_props:list=None,put2folder:str=None,as_batch=True):
        if hasattr(self,'id2groups'):
            if not isinstance(group_urn,str):
                try:
                    group_urn = self.id2groups[group_id]['URN'][0]
                    group_name = self.id2groups[group_id]['Name'][0]
                except KeyError:
                    print('Pathway collection does not have %s pathway with URN %s' % (group_name,group_urn))

        if not isinstance(group_urn,str):
            print('Pathway has no URN specifed!!!!')
            group_urn = 'no_urn'
        
        if not isinstance(group_name,str):
            print('Pathway has no Name specifed!!!!')
            group_name = 'no_name'

        group_graph = self.load_graph_from_oql('SELECT Entity WHERE MemberOf (SELECT Group WHERE Name = \'{name}\')'.format(name = group_name),entity_props=ent_props,get_links=False)
       
        rnef_xml = et.fromstring(self.to_rnef(group_graph))
        rnef_xml.set('name', group_name)
        rnef_xml.set('urn', group_urn)
        rnef_xml.set('type', 'Group')
        
        batch_xml = et.Element('batch')
        batch_xml.insert(0,rnef_xml)
                   
        if isinstance(put2folder,str):
            folder_resnet = et.Element('resnet')
            xml_nodes = et.SubElement(folder_resnet, 'nodes')
            folder_local_id = 'F0'
            xml_node_folder = et.SubElement(xml_nodes, 'node', {'local_id':folder_local_id, 'urn': 'urn:agi-folder:xxxxx_yyyyy_zzzzz'})
            et.SubElement(xml_node_folder, 'attr', {'name': 'NodeType', 'value': 'Folder'})
            et.SubElement(xml_node_folder, 'attr', {'name': 'Name', 'value': put2folder})
            pathway_local_id = 'P0'
            xml_node_pathway = et.SubElement(xml_nodes, 'node', {'local_id':pathway_local_id, 'urn':group_urn})
            et.SubElement(xml_node_pathway, 'attr', {'name': 'NodeType', 'value': 'Group'})
            xml_controls = et.SubElement(folder_resnet, 'controls')
            xml_control = et.SubElement(xml_controls, 'control', {'local_id':'CFE1'})
            et.SubElement(xml_control, 'attr', {'name':'ControlType', 'value':'MemberOf'})
            et.SubElement(xml_control, 'link', {'type':'in', 'ref':pathway_local_id})
            et.SubElement(xml_control, 'link', {'type':'out', 'ref':folder_local_id})
            batch_xml.append(folder_resnet)
        
        if as_batch:
            group_xml = et.tostring(batch_xml,encoding='utf-8',xml_declaration=True).decode("utf-8")
            group_xml = minidom.parseString(group_xml).toprettyxml(indent='   ')
        else:
            group_xml = et.tostring(rnef_xml,encoding='utf-8',xml_declaration=True).decode("utf-8")
            group_xml = str(minidom.parseString(group_xml).toprettyxml(indent='   '))
            group_xml = group_xml[group_xml.find('\n')+1:]
            #minidom does not work without xml_declaration

        print('\"%s\" group downloaded: %d nodes' % (group_name, group_graph.number_of_nodes()))
        return group_graph, str(group_xml)

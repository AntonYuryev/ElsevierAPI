import networkx as nx
import pandas as pd
import xml.etree.ElementTree as et
from ElsevierAPI.ResnetAPI.NetworkxObjects import PSObject, PSRelation

NO_RNEF_NODE_PROPS = ['Id','URN','ObjClassId','ObjTypeId','ObjTypeName','OwnerId','DateCreated','DateModified']
NO_RNEF_REL_PROPS = NO_RNEF_NODE_PROPS + ['RelationNumberOfReferences', '# of Total References', 'Name']

class ResnetGraph (nx.MultiDiGraph):
    def __init__(self):
        super().__init__()

    def get_entity_ids(self, SearchValues: list, in_graph=None, search_by_properties: list=None):
        search_by_properties = ['ObjTypeName'] if search_by_properties is None else search_by_properties
        graph = in_graph if isinstance(in_graph, ResnetGraph) else self
        all_ids = set()
        for i, node in graph.nodes(data=True):
            for propName in search_by_properties:
                if PSObject(node).has_str_property(propName, SearchValues):
                    all_ids.add(i)
                    break
        return list(all_ids)

    def _get_node(self, nodeId, in_graph=None):
        graph = in_graph if isinstance(in_graph, ResnetGraph) else self
        #dic = {k: v for k, v in graph.nodes[nodeId].items()}
        #ps_obj = PSObject(dic)
        return PSObject({k: v for k, v in graph.nodes[nodeId].items()})

    def _get_nodes(self, nodeIds: list, in_graph=None):
        graph = in_graph if isinstance(in_graph, ResnetGraph) else self
        node_list = list()
        for n in nodeIds:
            node_list.append(PSObject({k: v for k, v in graph.nodes[n].items()}))
        return node_list

    def __get_relations(self, in_graph=None):
        graph = in_graph if isinstance(in_graph, ResnetGraph) else self
        graph_relations = set()
        for regulatorID, targetID, rel in graph.edges.data('relation'):
            graph_relations.add(rel)

        return graph_relations

    def get_properties(self, IDList: set, PropertyName):
        id2props = {x: y[PropertyName] for x, y in self.nodes(data=True) if x in IDList}
        return id2props

    def get_neighbors(self, EntityIDs: set, in_graph=None, only_neighbors_with_ids=None):
        if only_neighbors_with_ids is None:
            only_neighbors_with_ids = []
        graph = in_graph if isinstance(in_graph, ResnetGraph) else self
        id2neighbors = dict()
        for Id in EntityIDs:
            if Id in graph:
                neighbors_ids = set([x for x in nx.all_neighbors(graph, Id)])
                id2neighbors[Id] = list(neighbors_ids)
        if len(only_neighbors_with_ids) > 0:
            filtered_id2neighbors = dict()
            for k, v in id2neighbors.items():
                filtered_neighbors_ids = [i for i in v if i in only_neighbors_with_ids]
                filtered_id2neighbors[k] = filtered_neighbors_ids
            return filtered_id2neighbors
        else:
            return id2neighbors

    def expand_nodes(self, EntityIDs: set, in_graph=None, only_neighbors_with_ids=None):
        only_neighbors_with_ids = [] if only_neighbors_with_ids is None else only_neighbors_with_ids
        graph = in_graph if isinstance(in_graph, ResnetGraph) else self
        neighbors = set()
        for Id in EntityIDs:
            if Id in graph:
                neighbors.update(set([x for x in nx.all_neighbors(graph, Id)]))
        if len(only_neighbors_with_ids) > 0:
            return [i for i in list(neighbors) if i in only_neighbors_with_ids]
        else:
            return list(neighbors)

    def get_regulome(self, StartEntityIDs: set, in_graph=None):
        graph = in_graph if isinstance(in_graph, ResnetGraph) else self
        all_trees = nx.DiGraph()
        for Id in StartEntityIDs:
            if Id in graph:
                t = nx.bfs_tree(graph, Id)
                all_trees = nx.compose(all_trees, t)
        return all_trees

    def get_subgraph(self, between_node_ids: list, and_node_ids: list, in_graph=None):
        graph = in_graph if isinstance(in_graph, ResnetGraph) else self
        subgraph = ResnetGraph()
        for n1 in between_node_ids:
            for n2 in and_node_ids:
                if graph.has_edge(n1, n2):
                    for i in range(0, len(graph[n1][n2])):
                        subgraph.add_edge(n1, n2, relation=graph[n1][n2][i]['relation'],
                                          weight=graph[n1][n2][i]['weight'])
                if graph.has_edge(n2, n1):
                    for i in range(0, len(graph[n2][n1])):
                        subgraph.add_edge(n2, n1, relation=graph[n2][n1][i]['relation'],
                                          weight=graph[n2][n1][i]['weight'])

        return subgraph

    def count_references(self, in_graph=None):
        if not isinstance(in_graph, ResnetGraph): in_graph = self
        references = set()
        for regulatorID, targetID, rel in in_graph.edges.data('relation'):
            rel.load_references()
            references.update(rel.References.values())
        return references


    def set_edge_property(self, nodeId1, nodeId2, PropertyName, PropertyValues: list, bothDirs=True):
        if self.has_edge(nodeId1, nodeId2):
            for i in range(0, len(self[nodeId1][nodeId2])):
                self[nodeId1][nodeId2][i]['relation'][PropertyName] = PropertyValues
        if bothDirs:
            if self.has_edge(nodeId2, nodeId1):
                for i in range(0, len(self[nodeId2][nodeId1])):
                    self[nodeId2][nodeId1][i]['relation'][PropertyName] = PropertyValues


    def print_triples(self, fileOut, relPropNames, in_graph=None, access_mode='w', printHeader=True):
        graph = in_graph if isinstance(in_graph, ResnetGraph) else self
        with open(fileOut, access_mode, encoding='utf-8') as f:
            if printHeader:
                header = '\t'.join(relPropNames) + '\t' + "Regulators Id" + '\t' + "Targets Id"
                f.write(header + '\n')

        for regulatorID, targetID, rel in graph.edges.data('relation'):
            f.write(rel.triple2str(relPropNames))




    def print_references(self, fileOut, relPropNames, entity_prop_names=None, in_graph=None,
                         access_mode='w', printHeader=True, RefNumPrintLimit=0, col_sep: str='\t'):

        entity_prop_names = [] if entity_prop_names is None else entity_prop_names
        rel_props = relPropNames
        in_graph = in_graph if isinstance(in_graph, ResnetGraph) else self

        with open(fileOut, access_mode, encoding='utf-8') as f:
            if printHeader:
                header = col_sep.join(rel_props) + col_sep + "Regulators Id" + col_sep + "Targets Id"
                target_header = [''] * len(entity_prop_names)
                reg_header = [''] * len(entity_prop_names)

                for i in range(0, len(entity_prop_names)):
                    reg_header[i] = 'Regulator:' + entity_prop_names[i]
                    target_header[i] = 'Target:' + entity_prop_names[i]

                if len(reg_header) > 0:
                    header = col_sep.join(reg_header) + col_sep + header + col_sep + col_sep.join(target_header)
                f.write(header + '\n')

            if in_graph.number_of_edges() > 0:
                if len(entity_prop_names) == 0:
                    for regulatorID, targetID, rel in in_graph.edges.data('relation'):
                        reference_view_triple = str(rel.triple2str(rel_props))
                        f.write(reference_view_triple)
                else:
                    for regulatorID, targetID, rel in in_graph.edges.data('relation'):
                        reg = self._get_node(regulatorID, in_graph)
                        target = self._get_node(targetID, in_graph)
                        reg_props_str = reg.data2str(entity_prop_names, col_sep=col_sep)
                        target_props_str = target.data2str(entity_prop_names, col_sep=col_sep)
                        rel_props_str_list = dict(
                            rel.triple2str(rel_props, return_dict=True, RefNumPrintLimit=RefNumPrintLimit,
                                           col_sep=col_sep))

                        reference_table_view = str()
                        for row in rel_props_str_list.values():
                            reference_table_view += reg_props_str[0:len(reg_props_str) - 1]
                            reference_table_view += col_sep + col_sep.join(row) + col_sep + target_props_str
                        f.write(reference_table_view)
            else:
                for node_id in in_graph.nodes:
                    n = self._get_node(node_id, in_graph)
                    node_prop_str = n.data2str(entity_prop_names, col_sep=col_sep)
                    f.write(node_prop_str)


    def ref2pandas (self, relPropNames: list, entity_prop_names=None, in_graph=None, RefNumPrintLimit=0):
        entity_prop_names = [] if entity_prop_names is None else entity_prop_names
        rel_props = relPropNames
        graph = in_graph if isinstance(in_graph, ResnetGraph) else self
        
        target_header = [''] * len(entity_prop_names)
        reg_header = [''] * len(entity_prop_names)
        for i in range(0, len(entity_prop_names)):
            reg_header[i] = 'Regulator:' + entity_prop_names[i]
            target_header[i] = 'Target:' + entity_prop_names[i]

        cols = reg_header + rel_props + target_header + ["Regulators Id","Targets Id"]
        reference_pandas = pd.DataFrame(columns = cols)

        if graph.number_of_edges() > 0:
            if len(entity_prop_names) == 0:
                for regulatorID, targetID, rel in in_graph.edges.data('relation'):
                    relation_pandas = rel.to_pandas(rel_props,RefNumPrintLimit)
                    reference_pandas = pd.concat([reference_pandas, relation_pandas], axis=0)
            else:
                for regulatorID, targetID, rel in graph.edges.data('relation'):
                    reg = self._get_node(regulatorID, in_graph)
                    target = self._get_node(targetID, in_graph)
                    reg_pandas = reg.to_pandas(entity_prop_names)
                    reg_pandas = reg_pandas.add_prefix('Regulator:')
                    target_pandas = target.to_pandas(entity_prop_names)
                    target_pandas = target_pandas.add_prefix('Target:')
                    relation_pandas = rel.to_pandas(rel_props,RefNumPrintLimit)

                    for col in reg_pandas.columns:
                        relation_pandas[col] = reg_pandas.at[0,col]
                    
                    for col in target_pandas.columns:
                        relation_pandas[col] = target_pandas.at[0,col]
                    
                    reference_pandas = pd.concat([reference_pandas, relation_pandas], axis=0)
        else:
            for node_id in graph.nodes:
                n = self._get_node(node_id, in_graph)
                node_pandas = n.to_pandas(entity_prop_names)
                reference_pandas = pd.concat([reference_pandas, node_pandas], axis=0)

        return reference_pandas


    def dump_entities(self, fileOut, PropNames, entity_ids=None):
        entity_ids = [] if entity_ids is None else entity_ids
        header = '\t'.join(PropNames)
        node_dict = dict([(i, v) for i, v in self.nodes(data=True)])
        with open(fileOut, 'w', encoding='utf-8') as f:
            f.write(header + '\n')
            if len(entity_ids) == 0:
                for ent in node_dict.values():
                    f.write(ent.data2str(PropNames))
            else:
                for ent in node_dict.values():
                    if ent['Id'] in entity_ids:
                        f.write(ent.data2str(PropNames))


    def rename_rel_property(self, oldPropertyName='MedlineTA', newPropertyName='Journal'):
        for regulatorID, targetID, rel in self.edges.data('relation'):
            try:
                rel[newPropertyName] = rel.pop(oldPropertyName)
            except KeyError:
                for prop in rel.PropSetToProps.values():
                    try:
                        prop[newPropertyName] = prop.pop(oldPropertyName)
                    except KeyError:
                        continue
                continue


    def to_rnef(self, RNEFnameToPropType: list, in_graph=None):
        if not isinstance(in_graph, ResnetGraph): in_graph = self
        resnet = et.Element('resnet')
        xml_nodes = et.SubElement(resnet, 'nodes')
        local_id_counter = 0
        for nodeId, n in in_graph.nodes(data=True):
            local_id = n['URN'][0]
            xml_node = et.SubElement(xml_nodes, 'node', {'local_id': local_id, 'urn': n['URN'][0]})
            et.SubElement(xml_node, 'attr', {'name': str('NodeType'), 'value': str(n['ObjTypeName'][0])})
            for prop_name, prop_values in n.items():
                if prop_name in RNEFnameToPropType:
                    if prop_name not in NO_RNEF_NODE_PROPS:
                        for prop_value in prop_values:
                            et.SubElement(xml_node, 'attr', {'name': str(prop_name), 'value': str(prop_value)})

            local_id_counter += 1

        xml_controls = et.SubElement(resnet, 'controls')

        graph_relations = self.__get_relations(in_graph)
        for rel in graph_relations:
            control_id = rel['URN'][0]
            xml_control = et.SubElement(xml_controls, 'control', {'local_id': control_id})
            et.SubElement(xml_control, 'attr', {'name': str('ControlType'), 'value': str(rel['ObjTypeName'][0])})
            # adding links
            regulators = rel.Nodes['Regulators']
            try:
                targets = rel.Nodes['Targets']
                for r in regulators:
                    regulator_local_id = in_graph.nodes[r[0]]['URN'][0]
                    et.SubElement(xml_control, 'link', {'type': 'in', 'ref': regulator_local_id})

                for t in targets:
                    target_local_id = in_graph.nodes[t[0]]['URN'][0]
                    et.SubElement(xml_control, 'link', {'type': 'out', 'ref': target_local_id})
            except KeyError:
                for r in regulators:
                    regulator_local_id = in_graph.nodes[r[0]]['URN'][0]
                    et.SubElement(xml_control, 'link', {'type': 'in-out', 'ref': regulator_local_id})

            # adding non-reference properties
            for prop_name, prop_values in rel.items():
                if prop_name in RNEFnameToPropType:
                    if prop_name not in NO_RNEF_REL_PROPS:
                        for prop_value in prop_values:
                            et.SubElement(xml_control, 'attr', {'name': str(prop_name), 'value': str(prop_value)})

            # adding references
            references = list(set(rel.References.values()))
            for i in range(0, len(references)):
                for prop_name, prop_values in references[i].items():
                    for prop_value in prop_values:
                        if prop_name in RNEFnameToPropType:
                            et.SubElement(
                                xml_control, 'attr',{'name': str(prop_name), 'value': str(prop_value), 'index': str(i)})
                for ref_id_type,ref_id in references[i].Identifiers.items():
                    et.SubElement(
                        xml_control, 'attr',{'name': str(ref_id_type), 'value': str(ref_id), 'index': str(i)})

        xml_str = et.tostring(resnet, encoding='utf-8').decode("utf-8")
        return xml_str

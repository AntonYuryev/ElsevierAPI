import networkx as nx
import pandas as pd
import xml.etree.ElementTree as et
from ElsevierAPI.ResnetAPI.NetworkxObjects import PSObject

NO_RNEF_NODE_PROPS = ['Id','URN','ObjClassId','ObjTypeId','ObjTypeName','OwnerId','DateCreated','DateModified']
NO_RNEF_REL_PROPS = NO_RNEF_NODE_PROPS + ['RelationNumberOfReferences', '# of Total References', 'Name']

class ResnetGraph (nx.MultiDiGraph):
    def __init__(self):
        super().__init__()

    def add_graph(self, other: "ResnetGraph"):
        self.update(other)
      
    def get_entity_ids(self, SearchValues: list, search_by_properties: list=None):
        if search_by_properties is None: search_by_properties = ['ObjTypeName']
        all_ids = set()
        for i, node in self.nodes(data=True):
            for propName in search_by_properties:
                if PSObject(node).has_str_property(propName, SearchValues):
                    all_ids.add(i)
                    break
        return list(all_ids)

    def _get_node(self, nodeId):
        return PSObject({k: v for k, v in self.nodes[nodeId].items()})

    def _get_nodes(self, nodeIds: list):
        node_list = list()
        for n in nodeIds:
            node_list.append(PSObject({k: v for k, v in self.nodes[n].items()}))
        return node_list

    def __relations(self):
        graph_relations = set()
        for regulatorID, targetID, rel in self.edges.data('relation'):
            graph_relations.add(rel)

        return graph_relations

    def get_properties(self, IDList: set, PropertyName):
        id2props = {x: y[PropertyName] for x, y in self.nodes(data=True) if x in IDList}
        return id2props

    def get_neighbors(self, node_ids: set, only_neighbors_with_ids=None):
        if not isinstance(only_neighbors_with_ids,list): only_neighbors_with_ids = list()
        neighbors = set()
        for Id in node_ids:
                neighbors.update(set([x for x in nx.all_neighbors(self, Id)]))
        if len(only_neighbors_with_ids) > 0:
            return [i for i in list(neighbors) if i in only_neighbors_with_ids]
        else:
            return list(neighbors)

    def get_regulome(self, StartEntityIDs: set):
        all_trees = nx.DiGraph()
        for Id in StartEntityIDs:
                t = nx.bfs_tree(self, Id)
                all_trees = nx.compose(all_trees, t)
        return all_trees

    def get_subgraph(self,between_node_ids:list,and_node_ids:list,by_relation_type=None,with_effect=None,in_direction=None)->"ResnetGraph":
        subgraph = ResnetGraph()
        for n1 in between_node_ids:
            for n2 in and_node_ids:
                if not isinstance(in_direction,str) or in_direction == '>':
                    if self.has_edge(n1, n2):
                        for i in range(0, len(self[n1][n2])):
                            rel = self[n1][n2][i]['relation']
                            if isinstance(by_relation_type,list) and rel['ObjTypeName'] not in by_relation_type: continue
                            if isinstance(with_effect,list):
                                try: ef = rel['Effect'] 
                                except KeyError:
                                    rel['Effect'] = 'unknown'
                                    ef = rel['Effect']
                                if ef not in with_effect: continue
                                
                            subgraph.add_edge(n1, n2, relation=rel,weight=self[n1][n2][i]['weight'])
                            subgraph.add_nodes_from([(n1, self.nodes[n1]),(n2, self.nodes[n2])])
                if not isinstance(in_direction,str) or in_direction == '<':
                    if self.has_edge(n2, n1):
                        for i in range(0, len(self[n2][n1])):
                            rel = self[n2][n1][i]['relation']
                            if isinstance(by_relation_type,list) and rel['ObjTypeName'] not in by_relation_type: continue
                            if isinstance(with_effect,list):
                                try: ef = rel['Effect'] 
                                except KeyError:
                                    rel['Effect'] = 'unknown'
                                    ef = rel['Effect']
                                if ef not in with_effect: continue
                                    
                            subgraph.add_edge(n2, n1, relation=rel,weight=self[n2][n1][i]['weight'])
                            subgraph.add_nodes_from([(n1, self.nodes[n1]),(n2, self.nodes[n2])])
            
        return subgraph

    def get_neighbors_graph(self, node_ids:set, only_neighbors_with_ids=None):
        neighbors_ids = self.get_neighbors(node_ids,only_neighbors_with_ids)
        return self.get_subgraph(node_ids,neighbors_ids)

    def count_references(self):
        references = set()
        for regulatorID, targetID, rel in self.edges.data('relation'):
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


    def print_triples(self, fileOut, relPropNames, access_mode='w', printHeader=True):
        
        with open(fileOut, access_mode, encoding='utf-8') as f:
            if printHeader:
                header = '\t'.join(relPropNames) + '\t' + "Regulators Id" + '\t' + "Targets Id"
                f.write(header + '\n')

        for regulatorID, targetID, rel in self.edges.data('relation'):
            f.write(rel.triple2str(relPropNames))




    def print_references(self, fileOut, relPropNames, entity_prop_names=None,
                         access_mode='w', printHeader=True, RefNumPrintLimit=0, col_sep: str='\t', debug=False):

        if entity_prop_names is None: entity_prop_names = []
        
        with open(fileOut, access_mode, encoding='utf-8') as f:
            if printHeader:
                header = col_sep.join(relPropNames)
                if debug:
                    header = header + col_sep + "Regulators Id" + col_sep + "Targets Id"

                target_header = [''] * len(entity_prop_names)
                reg_header = [''] * len(entity_prop_names)
                for i in range(0, len(entity_prop_names)):
                    reg_header[i] = 'Regulator:' + entity_prop_names[i]
                    target_header[i] = 'Target:' + entity_prop_names[i]

                if len(reg_header) > 0:
                    header = col_sep.join(reg_header) + col_sep + header + col_sep + col_sep.join(target_header)
                
                f.write(header + '\n')

            if self.number_of_edges() > 0:
                if len(entity_prop_names) == 0:
                    for regulatorID, targetID, rel in self.edges.data('relation'):
                        reference_view_triple = str(rel.triple2str(relPropNames))
                        f.write(reference_view_triple)
                else:
                    # rel = PSRelation()
                    for regulatorID, targetID, rel in self.edges.data('relation'):
                        reg = self._get_node(regulatorID)
                        target = self._get_node(targetID)
                        reg_props_str = reg.data2str(entity_prop_names, col_sep=col_sep)
                        target_props_str = target.data2str(entity_prop_names, col_sep=col_sep)
                        rel_props_str_list = dict(
                            rel.triple2str(relPropNames, return_dict=True, RefNumPrintLimit=RefNumPrintLimit,
                                           col_sep=col_sep, add_entities=debug))

                        reference_table_view = str()
                        for row in rel_props_str_list.values():
                            reference_table_view += reg_props_str[0:len(reg_props_str) - 1]
                            reference_table_view += col_sep + col_sep.join(row) + col_sep + target_props_str
                        f.write(reference_table_view)
            else:
                for node_id in self.nodes:
                    n = self._get_node(node_id)
                    node_prop_str = n.data2str(entity_prop_names, col_sep=col_sep)
                    f.write(node_prop_str)


    def ref2pandas (self, relPropNames: list, entity_prop_names=None, RefNumPrintLimit=0):
        entity_prop_names = [] if entity_prop_names is None else entity_prop_names
        rel_props = relPropNames
        
        
        target_header = [''] * len(entity_prop_names)
        reg_header = [''] * len(entity_prop_names)
        for i in range(0, len(entity_prop_names)):
            reg_header[i] = 'Regulator:' + entity_prop_names[i]
            target_header[i] = 'Target:' + entity_prop_names[i]

        cols = reg_header + rel_props + target_header + ["Regulators Id","Targets Id"]
        reference_pandas = pd.DataFrame(columns = cols)

        if self.number_of_edges() > 0:
            if len(entity_prop_names) == 0:
                for regulatorID, targetID, rel in self.edges.data('relation'):
                    relation_pandas = rel.to_pandas(rel_props,RefNumPrintLimit)
                    reference_pandas = pd.concat([reference_pandas, relation_pandas], axis=0)
            else:
                for regulatorID, targetID, rel in self.edges.data('relation'):
                    reg = self._get_node(regulatorID)
                    target = self._get_node(targetID)
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
            for node_id in self.nodes:
                n = self._get_node(node_id)
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


    def to_rnef(self, RNEFprops: list):
        resnet = et.Element('resnet')
        xml_nodes = et.SubElement(resnet, 'nodes')
        local_id_counter = 0
        for nodeId, n in self.nodes(data=True):
            local_id = n['URN'][0]
            xml_node = et.SubElement(xml_nodes, 'node', {'local_id': local_id, 'urn': n['URN'][0]})
            et.SubElement(xml_node, 'attr', {'name': str('NodeType'), 'value': str(n['ObjTypeName'][0])})
            for prop_name, prop_values in n.items():
                if prop_name in RNEFprops:
                    if prop_name not in NO_RNEF_NODE_PROPS:
                        for prop_value in prop_values:
                            et.SubElement(xml_node, 'attr', {'name': str(prop_name), 'value': str(prop_value)})

            local_id_counter += 1

        xml_controls = et.SubElement(resnet, 'controls')

        graph_relations = self.__relations()
        for rel in graph_relations:
            control_id = rel['URN'][0]
            xml_control = et.SubElement(xml_controls, 'control', {'local_id': control_id})
            et.SubElement(xml_control, 'attr', {'name': str('ControlType'), 'value': str(rel['ObjTypeName'][0])})
            # adding links
            regulators = rel.Nodes['Regulators']
            try:
                targets = rel.Nodes['Targets']
                for r in regulators:
                    regulator_local_id = self.nodes[r[0]]['URN'][0]
                    et.SubElement(xml_control, 'link', {'type': 'in', 'ref': regulator_local_id})

                for t in targets:
                    target_local_id = self.nodes[t[0]]['URN'][0]
                    et.SubElement(xml_control, 'link', {'type': 'out', 'ref': target_local_id})
            except KeyError:
                for r in regulators:
                    regulator_local_id = self.nodes[r[0]]['URN'][0]
                    et.SubElement(xml_control, 'link', {'type': 'in-out', 'ref': regulator_local_id})

            # adding non-reference properties
            for prop_name, prop_values in rel.items():
                if prop_name in RNEFprops:
                    if prop_name not in NO_RNEF_REL_PROPS:
                        for prop_value in prop_values:
                            et.SubElement(xml_control, 'attr', {'name': str(prop_name), 'value': str(prop_value)})

            # adding references
            references = list(set(rel.References.values()))
            for i in range(0, len(references)):
                for prop_name, prop_values in references[i].items():
                    for prop_value in prop_values:
                        if prop_name in RNEFprops:
                            et.SubElement(
                                xml_control, 'attr',{'name': str(prop_name), 'value': str(prop_value), 'index': str(i)})
                for ref_id_type,ref_id in references[i].Identifiers.items():
                    et.SubElement(
                        xml_control, 'attr',{'name': str(ref_id_type), 'value': str(ref_id), 'index': str(i)})

        xml_str = et.tostring(resnet, encoding='utf-8').decode("utf-8")
        return xml_str

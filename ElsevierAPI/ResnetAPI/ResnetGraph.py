import networkx as nx
import pandas as pd
import os
import xml.etree.ElementTree as et
from ElsevierAPI.ResnetAPI.NetworkxObjects import PSObject, PSRelation

NO_RNEF_NODE_PROPS = ['Id','URN','ObjClassId','ObjTypeId','ObjTypeName','OwnerId','DateCreated','DateModified']
NO_RNEF_REL_PROPS = NO_RNEF_NODE_PROPS + ['RelationNumberOfReferences', '# of Total References', 'Name']

class ResnetGraph (nx.MultiDiGraph):
    def __init__(self):
        super().__init__()

    def add_graph(self, other: "ResnetGraph"):
        self.update(other)
      
    def get_entity_ids(self, SearchValues:list, search_by_properties:list=None):
        if search_by_properties is None: search_by_properties = ['ObjTypeName']
        all_ids = set()
        for i, node in self.nodes(data=True):
            for propName in search_by_properties:
                if PSObject(node).has_property(propName, SearchValues):
                    all_ids.add(i)
                    break
        return list(all_ids)

    def weight(self):
        return self.size(weight="weight")

    def subtract(self, other: "ResnetGraph"):
        #only self graph is analyzed 
        # works faster if other graph is bigger than self
        unique2self = ResnetGraph()
        edges_from_other = other._relations_ids()
        for n1,n2,e in self.edges(data=True):
            if e['relation']['Id'][0] not in edges_from_other:
                unique2self.add_edge(n1, n2, relation=e['relation'],weight=e['weight'])
                unique2self.add_nodes_from([(n1, self.nodes[n1]),(n2, self.nodes[n2])])
        return unique2self

    def intersect (self, other: "ResnetGraph"):
        #only self graph is analyzed 
        # works faster if other graph is bigger than self
        intersection = ResnetGraph()
        edges_from_other = other._relations_ids()
        for n1,n2,e in self.edges(data=True):
            if e['relation']['Id'][0] in edges_from_other:
                intersection.add_edge(n1, n2, relation=e['relation'],weight=e['weight'])
                intersection.add_nodes_from([(n1, self.nodes[n1]),(n2, self.nodes[n2])])
        return intersection

    def get_objects(self, SearchValues: list, search_by_properties: list=None):
        if search_by_properties is None: search_by_properties = ['ObjTypeName']
        all_objects = set()
        for i, node in self.nodes(data=True):
            for propName in search_by_properties:
                ps_obj = PSObject(node)
                if ps_obj.has_property(propName, SearchValues):
                    all_objects.add(ps_obj)
                    break
        return list(all_objects)

    def get_relations(self, search_values:list, search_by_properties: list=None):
        if search_by_properties is None: search_by_properties = ['ObjTypeName']
        relations2return = set()
        for regulatorID, targetID, rel in self.edges.data('relation'):
            for prop_name in search_by_properties:
                if rel.has_property(prop_name, search_values):
                    relations2return.add(rel)
                    break
        return relations2return

    def get_prop2obj_dic(self, search_by_property:str):
        to_return = dict()
        for i, n in self.nodes(data=True):
            try:
                for v in n[search_by_property]:
                    try:
                        to_return[v].append(PSObject(n))
                    except KeyError:
                        to_return[v] = [PSObject(n)]
            except KeyError:
                continue
        return to_return

    def _get_node(self, nodeId):
        return PSObject({k: v for k, v in self.nodes[nodeId].items()})

    def _get_nodes(self, node_ids:list=None):
        if not isinstance(node_ids,list):
            return [PSObject(ddict) for id,ddict in self.nodes(data=True)]
        else:
            return [PSObject(ddict) for id,ddict in self.nodes(data=True) if id in node_ids]

    def find_relations(self, reg_id, targ_id, rel_type, effect:str=None, mechanism:str=None):
        rel_set = [rel for regulatorID, targetID, rel in self.edges.data('relation') 
                if regulatorID == reg_id and targetID == targ_id and rel['ObjTypeName'][0]==rel_type
        ]
        if isinstance(effect,str):
            rel_set = [x for x in rel_set if x['Effect'] == effect]
        
        if isinstance(mechanism,str):
            rel_set = [x for x in rel_set if x['Mechanism'] == mechanism]

        return rel_set

    def _relations_ids(self):
        return {rel['Id'][0] for regulatorID, targetID, rel in self.edges.data('relation')}

    def find_regulators(self, for_relation:PSRelation, filter_by:list=None, in_properties:list=None):
        regulators_ids = for_relation.get_regulator_ids()
        regulators = self.get_objects(regulators_ids,['Id'])

        if isinstance(filter_by,list):
            search_in_properties = in_properties if isinstance(in_properties,list) else ['ObjTypeName']
            must_values = set(filter_by)
            filtered_regulators = set()
            for prop_name in search_in_properties:
                filtered_regulators.update([x for x in regulators if not set(x[prop_name]).isdisjoint(must_values)]) 
            return filtered_regulators
        else:
            return regulators

    def get_properties(self, IDList: set, PropertyName):
        id2props = {x: y[PropertyName] for x, y in self.nodes(data=True) if x in IDList}
        return id2props

    def get_neighbors(self, node_ids: set, only_neighbors_with_ids=None):
        if not isinstance(only_neighbors_with_ids,list): only_neighbors_with_ids = list()
        neighbors = set()
        for i in [n for n in node_ids if self.has_node(n)]:
            neighbors.update(set([x for x in nx.all_neighbors(self, i)]))
                
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

    def subgraph_by_relprops(self, search_values:list, in_properties:list=None):
        if in_properties is None: in_properties = ['ObjTypeName']
        search_value_set = set(search_values)
        subgraph = ResnetGraph()
        for prop_type in in_properties:
            for n1, n2, rel in self.edges.data('relation'):
                if set(rel[prop_type]).intersection(search_value_set):
                    subgraph.add_edge(n2, n1, relation=rel,weight=rel['RelationNumberOfReferences'])
                    subgraph.add_nodes_from([(n1, self.nodes[n1]),(n2, self.nodes[n2])])
        return subgraph


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
            
    def count_references(self, weight_by_prop_name:str=None, proval2weight:dict=None):
        references = set()
        if isinstance(weight_by_prop_name,str):
            for regulatorID, targetID, rel in self.edges.data('relation'):
                rel.load_references()
                rel._weight2ref(weight_by_prop_name,proval2weight)
                references.update(rel.References.values())
        else:
            for regulatorID, targetID, rel in self.edges.data('relation'):
                rel.load_references()
                references.update(rel.References.values())
        return references

    def count_references_between(self, between_node_ids: list, and_node_ids: list, weight_by_prop_name:str=None, proval2weight:dict=None):
        sub_graph = self.get_subgraph(between_node_ids, and_node_ids)
        return sub_graph.count_references(weight_by_prop_name,proval2weight)

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


    def print_references(self, fileOut:str, relPropNames:list, entity_prop_names=[],access_mode='w',
                          printHeader=True, RefNumPrintLimit=0, col_sep:str='\t', debug=False, single_rel_row=False):
        
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
                if not entity_prop_names:
                    for regulatorID, targetID, rel in self.edges.data('relation'):
                        reference_view_triple = str(rel.triple2str(relPropNames,as1row=single_rel_row))
                        f.write(reference_view_triple)
                else:
                    # rel = PSRelation()
                    for regulatorID, targetID, rel in self.edges.data('relation'):
                        reg = self._get_node(regulatorID)
                        target = self._get_node(targetID)
                        reg_props_str = reg.data2str(entity_prop_names, col_sep=col_sep)
                        target_props_str = target.data2str(entity_prop_names, col_sep=col_sep)

                        reference_table_view = str()
                        if single_rel_row:
                            reference_table_view = rel.to1row(relPropNames, RefNumPrintLimit=RefNumPrintLimit,
                                                                col_sep=col_sep, add_entities=debug)

                            reference_table_view = reference_table_view[0:len(reference_table_view) - 1]#remove end of line character
                            reference_table_view = reg_props_str[0:len(reg_props_str) - 1]+col_sep+reference_table_view
                            reference_table_view += col_sep + target_props_str
                        else:
                            rel_props_str_list = dict(rel.to_table_dict(
                                                        relPropNames, RefNumPrintLimit=RefNumPrintLimit,add_entities=debug))
                            
                            for row in rel_props_str_list.values():
                                reference_table_view += reg_props_str[0:len(reg_props_str) - 1]
                                reference_table_view += col_sep + col_sep.join(row) + col_sep + target_props_str

                        f.write(reference_table_view)
            else:
                for node_id in self.nodes:
                    n = self._get_node(node_id)
                    node_prop_str = n.data2str(entity_prop_names, col_sep=col_sep)
                    f.write(node_prop_str)


    def ref2pandas (self, relPropNames:list, entity_prop_names=[], RefNumPrintLimit=0) -> 'pd.DataFrame':
        temp_fname = '__temp__.tsv'
        self.print_references(temp_fname,relPropNames,entity_prop_names,RefNumPrintLimit=RefNumPrintLimit)
        to_return = pd.read_csv(temp_fname,sep='\t',header=0,index_col=False, dtype='unicode')
        os.remove(temp_fname)
        return to_return

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


    def to_rnef(self, RNEFprops: list,add_rel_props:dict=None,add_pathway_props:dict=None):
        # add_rel_props,add_pathway_props structure - {PropName:[PropValues]}
        resnet = et.Element('resnet')
        if isinstance(add_pathway_props,dict):
            pathway_props = et.SubElement(resnet, 'properties')
            for prop_name,prop_val in add_pathway_props.items():
                    for val in prop_val:
                        et.SubElement(pathway_props, 'attr', {'name':str(prop_name), 'value':str(val)})

        xml_nodes = et.SubElement(resnet, 'nodes')
        local_id_counter = 0
        for nodeId, n in self.nodes(data=True):
            local_id = n['URN'][0]
            xml_node = et.SubElement(xml_nodes, 'node', {'local_id': local_id, 'urn': n['URN'][0]})
            et.SubElement(xml_node, 'attr', {'name': 'NodeType', 'value': str(n['ObjTypeName'][0])})
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
            xml_control = et.SubElement(xml_controls, 'control', {'local_id':control_id})
            et.SubElement(xml_control, 'attr', {'name':'ControlType', 'value':str(rel['ObjTypeName'][0])})
            
            # adding links
            regulators = rel.Nodes['Regulators']
            try:
                targets = rel.Nodes['Targets']
                for r in regulators:
                    regulator_local_id = self.nodes[r[0]]['URN'][0]
                    et.SubElement(xml_control, 'link', {'type':'in', 'ref':regulator_local_id})

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
                            et.SubElement(xml_control, 'attr', {'name':str(prop_name), 'value':str(prop_value)})

            if isinstance(add_rel_props,dict):
                for prop_name,prop_val in add_rel_props.items():
                    for val in prop_val:
                        et.SubElement(xml_control, 'attr', {'name':str(prop_name), 'value':str(val)})

            # adding references
            references = list(set(rel.References.values()))
            index = 0
            for ref in references:
                ref_props = {k:v for k,v in ref.items() if k in RNEFprops}
                for textref, sentence_props in ref.snippets.items():
                    et.SubElement(xml_control, 'attr',{'name': str('TextRef'), 'value': textref, 'index': str(index)})
                    for sentprop_name, sentprop_values in sentence_props.items():
                        for v in sentprop_values:
                            et.SubElement(xml_control, 'attr',{'name': str(sentprop_name), 'value': str(v), 'index': str(index)})

                        for prop_name, prop_values in ref_props.items():
                            for prop_value in prop_values:
                                    et.SubElement( xml_control, 'attr',{'name': str(prop_name), 'value': str(prop_value), 'index': str(index)})
                        for ref_id_type,ref_id in ref.Identifiers.items():
                            et.SubElement(xml_control, 'attr',{'name': str(ref_id_type), 'value': str(ref_id), 'index': str(index)})

                        index += 1

        xml_str = et.tostring(resnet, encoding='utf-8').decode("utf-8")
        return xml_str

    @staticmethod
    def get_att_set(prop_name:str,ps_objects:list):
        return set([i for sublist in [x[prop_name] for x in ps_objects] for i in sublist])

    def read_rnef(self, rnef_file:str, new_node_id=100000, new_control_id= 1000000):
        tree = et.parse(rnef_file)
        resnets = tree.findall('./batch/resnet')
        nodel_local_ids = dict()
        for resnet in resnets:
            for node in resnet.findall('./nodes/node'):
                node_urn = node.get('urn')
                local_id = node.get('local_id')
                node_objs = list(self.get_objects([node_urn],['URN']))
                if len(node_objs) > 0:
                    node_obj = node_objs[0]
                else:  
                    new_node_id += 1
                    node_id = new_node_id
                    node_obj = PSObject({'Id':[node_id]})

                for attr in node.findall('attr'):
                    prop_id = attr.get('name')
                    prop_value = attr.get('value')
                    if prop_id == 'NodeType': prop_id = 'ObjTypeName'
                    node_obj.append_property(prop_id, prop_value)

                self.add_node(node_id,node_obj.items())
                nodel_local_ids[local_id] = node_obj
                    

            for rel in resnet.findall('./controls/control'):
                regulators = list()
                targets = list()
                for link in rel.findall('link'):
                    link_type = link.get('type')
                    link_ref = link.get('ref')
                    node_obj = nodel_local_ids[link_ref]
                    if link_type == 'out': targets.append(node_obj)
                    else: regulators.append(node_obj)

                ps_rel = PSRelation(dict())
                effect = rel.find('./Effect')
                effect_val = 'unknown' if type(effect) == type(None) else effect.text
                mechanism = rel.find('./Mechanism')
                if type(mechanism) != type(None): mechanism = mechanism.text

                for reg in regulators:
                    try: 
                        ps_rel.Nodes['Regulators'].append(tuple(reg['Id'][0], '0', effect_val))
                    except KeyError:
                        ps_rel.Nodes['Regulators'] = [tuple(reg['Id'][0], '0', effect_val)]
                    
                for targ in targets:
                    try: 
                        ps_rel.Nodes['Targets'].append(tuple(targ['Id'][0], '0', effect_val))
                    except KeyError:
                        ps_rel.Nodes['Targets'] = [tuple(targ['Id'][0], '0', effect_val)]
                            
                for attr in rel.findall('attr'):
                    prop_id = attr.get('name')
                    prop_value = attr.get('value')
                    if prop_id == 'ControlType': prop_id = 'ObjTypeName'

                    index = attr.get('index')
                    if type(index) == type(None):
                        ps_rel.append_property(prop_id, prop_value)
                    else:
                        try:
                            props = ps_rel.PropSetToProps[index]
                            try:
                                props[prop_id].append(prop_value)
                            except KeyError:
                                 props[prop_id] = [prop_value]
                        except KeyError:
                            ps_rel.PropSetToProps[index] = {prop_id:prop_value}

                ref_count = len(ps_rel.PropSetToProps)
                if targets:
                    for r in regulators:
                        for t in targets:
                            existing_rel = self.find_relations(r['Id'][0],t['Id'][0],rel['ObjTypeName'][0],effect,mechanism)
                            if existing_rel:
                                for e in existing_rel:
                                    e.copy(rel)
                                    self.add_edge(r['Id'][0],t['Id'][0],relation=e, weight=float(len(e.PropSetToProps)))
                            else:
                                new_control_id += 1
                                rel['Id'] = [new_control_id]
                                self.add_edge(r['Id'][0],t['Id'][0],relation=rel, weight=float(ref_count))

                else:
                    for i in range(0, len(regulators)):
                        r = regulators[i]
                        for j in range(i, len(regulators)):
                            t = regulators[j]
                            existing_rel = self.find_relations(r['Id'][0],t['Id'][0],rel['ObjTypeName'][0],effect,mechanism)
                            if existing_rel:
                                for e in existing_rel:
                                    e.copy(rel)
                                    self.add_edge(r['Id'][0],t['Id'][0],relation=e, weight=float(len(e.PropSetToProps)))
                            else:
                                new_control_id += 1
                                rel['Id'] = [new_control_id]
                                self.add_edge(r['Id'][0],t['Id'][0],relation=rel, weight=float(ref_count))


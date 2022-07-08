import networkx as nx
import pandas as pd
import time
from datetime import timedelta
import os
from xml.dom import minidom
import xml.etree.ElementTree as et
from .NetworkxObjects import PSObject,PSRelation,REGULATORS,TARGETS,CHILDS,REFCOUNT
from ..ETM_API.references import PUBYEAR,EFFECT
from itertools import combinations
from ..pandas.panda_tricks import df

#NO_RNEF_NODE_PROPS = {'Id','URN','ObjClassId','ObjTypeId','ObjTypeName','OwnerId','DateCreated','DateModified'}
#NO_RNEF_REL_PROPS = NO_RNEF_NODE_PROPS | {'RelationNumberOfReferences', '# of Total References', 'Name'}

class ResnetGraph (nx.MultiDiGraph):
    def __init__(self):
        super().__init__()
        self.urn2obj = dict() #lookup for combining nodes in database graph and rnef graph
        self.urn2rel = dict() #lookup for combining relations in database graph and rnef graph

    def copy(self):
        cp = ResnetGraph()
        cp.update(self)
        cp.urn2obj = dict(self.urn2obj)
        cp.urn2rel = dict(self.urn2rel)
        return cp

    def load_urn_dicts(self):
        self.urn2obj = {o['URN'][0]:o for i,o in self.nodes(data=True)}
        self.urn2rel = {self.rel_urn(rel):rel for r,t,rel in self.edges(data='relation')}

    @staticmethod
    def execution_time(execution_start):
        return "{}".format(str(timedelta(seconds=time.time() - execution_start)))

###############    ADD ADD ADD    #############################
    def add1node(self,node:PSObject):
        self.add_nodes_from([(int(node['Id'][0]), node.items())])
        self.urn2obj[node['URN'][0]] = node

    def add_nodes(self,nodes:list): # nodes = [PSObject]
        self.add_nodes_from([(int(n['Id'][0]), n.items()) for n in nodes])
        [self.urn2obj.update({n['URN'][0]:n}) for n in nodes]


    def property2node(self, node_id,prop_name:str,prop_values:list):
        nx.set_node_attributes(self,{node_id:{prop_name:prop_values}})

    def __add_rel(self, rel:PSRelation):
        """
        nodes connected by rel must exist in the graph
        """
        ref_count = rel.get_reference_count()
        rel_urn = self.rel_urn(rel)
        if rel.is_directional():
            for regulator_id in rel.Nodes[REGULATORS]:
                for target_id in rel.Nodes[TARGETS]:
                    self.add_edge(regulator_id[0], target_id[0], relation=rel,weight=ref_count,key=rel_urn)
                    self.urn2rel[rel_urn] = rel
        else:
            reg_pairs = combinations(rel.Nodes[REGULATORS],2)
            [self.add_edge(pair[0][0], pair[1][0], relation=rel,weight=ref_count,key=rel_urn) for pair in reg_pairs]
            self.urn2rel[rel_urn] = rel


    def add_rel(self,rel:PSRelation):
        """
        nodes connected by relation must exist in the graph. written for read_rnef
        use copy_rel to move relation from one graph to another
        """
        rel_urn = self.rel_urn(rel)
        try:
            self.urn2rel[rel_urn].copy(rel)
            #existing_rel.copy(rel)
            #self.remove_relation(rel)
            self.__add_rel(self.urn2rel[rel_urn])
        except KeyError:
            self.urn2rel[rel_urn] = rel
            self.__add_rel(rel)
            

    def copy_rel(self,rel:PSRelation,from_graph:"ResnetGraph"):
        self.__add_rel(rel)
        add_nodes_with_ids = rel.get_regulator_ids() + rel.get_target_ids()
        nodes2add = list()
        for i in add_nodes_with_ids:
            #node_urn = from_graph.nodes[i]['URN'][0]
            node2add = from_graph.nodes(data=True)[i]
            try:
                self.urn2obj[node2add['URN'][0]]
            except KeyError:
                nodes2add.append((i,node2add))
                self.urn2obj[node2add['URN'][0]] = node2add

        self.add_nodes_from(nodes2add)

    
    def add_triple(self, regulator:PSObject,target:PSObject,rel_props: dict or PSRelation,refs=[],is_directional=True):
        """
        adds nodes and their relation to graph
        """
        self.add_nodes([regulator,target])

        if isinstance(rel_props,PSRelation):
            self.add_rel(rel_props)
            return
        else:
            rel = PSRelation.make_rel(regulator,target,rel_props,refs,is_directional)
            rel['Id'] = [self.number_of_edges()]
            self.rel_name(rel)
            rel_urn = self.rel_urn(rel)
            self.add_edge(regulator['Id'][0], target['Id'][0], relation=rel,weight=rel.get_reference_count(),key=rel_urn)
    

    def add_graph(self, other: "ResnetGraph"):
        self.update(other)

################## SET SET SET ##########################################
    def set_annotation(self, urn2value:dict, new_prop_name:str):
        """
            adds new property to exsisting nodes. Existing values of 'with_new_prop' will be replaced
            urn2value = {map2prop_value:[with_prop_values]}
        """
        id2value = dict() #remapping urn to ids
        for urn, value in urn2value.items():
            try:
                node2annotate = self.urn2node(urn)
                id2value[node2annotate['Id'][0]] = {new_prop_name:value}
            except KeyError: continue

        nx.set_node_attributes(self, id2value)
        #print('%d nodes were annotated with attributes "%s"' % (len(id2value),new_prop_name))


    def set_edge_property(self, nodeId1, nodeId2, PropertyName, PropertyValues: list, bothDirs=True):
        if self.has_edge(nodeId1, nodeId2):
            for i in range(0, len(self[nodeId1][nodeId2])):
                self[nodeId1][nodeId2][i]['relation'][PropertyName] = PropertyValues
        if bothDirs:
            if self.has_edge(nodeId2, nodeId1):
                for i in range(0, len(self[nodeId2][nodeId1])):#gives error if edge does not exist
                    self[nodeId2][nodeId1][i]['relation'][PropertyName] = PropertyValues


    def add_node_annotation(self, with_new_prop:str, map2prop:str, using_map:dict):
        """
            using_map = {map2prop_value:[annotations]}
        """
        annotation_counter = 0
        for i, node in self.nodes(data=True):
            try:
                # collecting all values to add 
                map_by_values = node[map2prop]
                annotate_with_values = set()
                for v in map_by_values:
                    try:
                        values2add = using_map[v]
                        annotate_with_values.update(values2add)
                    except KeyError: continue

                if annotate_with_values:
                    annotation_counter +=1
                    try:
                        # case when node has with_new_prop 
                        merged_annotation = set(node[with_new_prop]) | annotate_with_values
                        # set(node[with_new_prop]) - current existing annotation in the node
                        # annotate_with_values - new annotation
                        nx.set_node_attributes(self, {node['Id'][0]:{with_new_prop:list(merged_annotation)}})
                    except KeyError:
                        # case when node has no with_new_prop 
                        nx.set_node_attributes(self, {node['Id'][0]:{with_new_prop:list(annotate_with_values)}})
            except KeyError: continue
        print('%d nodes were annotated "%s" values out of %d "%s" values used for mapping' %
                (annotation_counter,with_new_prop,len(using_map), map2prop) )


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

################## SET-GET SET-GET SET-GET ###############################
    def rel_name(self, rel:PSRelation):
        try:
            return rel['Name'][0]
        except KeyError:
            regulators, targets = self.find_nodes(rel)
            reg_names = ','.join([r['Name'][0] for r in regulators])
            targ_names = ','.join([t['Name'][0] for t in targets])
            arrow = '--->'
            try:
                effect = rel[EFFECT]
                if effect == 'positive': arrow[2]='+'
                elif effect == 'negative': arrow[3]='|'
            except KeyError: pass
            name = reg_names+arrow+targ_names
            try:
                name += ':'+rel['Mechanism']
            except KeyError: pass
            rel['Name'] = [name]
            return name


    def rel_urn(self, rel:PSRelation):
        try:
            urn = rel['URN'][0]
        except KeyError:
            regulators, targets = self.find_nodes(rel)
            reg_urns = 'in:'.join([r['URN'][0] for r in regulators])
            urn = 'urn:agi-'+rel['ObjTypeName'][0]+':'+reg_urns
            if targets:
                urn += ':'+'out:'.join([r['URN'][0] for r in targets])
            try:
                urn += ':'+rel[EFFECT][0]
            except KeyError: pass
            try:
                urn += ':'+rel['Mechanism'][0]
            except KeyError: pass
            rel.set_property('URN', urn)
        return urn
    

    def load_references(self, weight_by_prop_name:str=None, proval2weight:dict=None):
        """
        uses rel.References
        does not work with REFCOUNT or 'RelationNumberOfReferences'
        """
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


    def load_references_between(self, between_node_ids: list, and_node_ids: list, weight_by_prop_name:str=None, proval2weight:dict=None):
        sub_graph = self.get_subgraph(between_node_ids, and_node_ids)
        return sub_graph.load_references(weight_by_prop_name,proval2weight)

##################### DEL DEL DEL ######################################
    def remove_relation(self, rel:PSRelation):
        rel_urn = self.rel_urn(rel)
        if rel.is_directional():
            for regulator_id in rel.Nodes[REGULATORS]:
                for target_id in rel.Nodes[TARGETS]:
                    if self.has_edge(regulator_id[0], target_id[0], key=rel_urn):
                        self.remove_edge(regulator_id[0], target_id[0], key=rel_urn)

        else:
            reg_pairs = combinations(rel.Nodes[REGULATORS],2)
            for pair in reg_pairs:
                if self.has_edge(regulator_id[0], target_id[0], key=rel_urn):
                    self.remove_edge(pair[0][0], pair[1][0], key=rel_urn)

        self.urn2rel.pop(rel_urn)

    def clear_graph(self):
        super().clear()

    def clear(self):
        super().clear()
        self.urn2obj.clear()
        self.urn2rel.clear()

######################    GET GET GET ######################################
    def get_node_ids(self, SearchValues:list, search_by_properties:list=None):
        if search_by_properties is None: search_by_properties = ['ObjTypeName']
        all_ids = set()
        for i, node in self.nodes(data=True):
            for propName in search_by_properties:
                if PSObject(node).is_annotated(propName, SearchValues):
                    all_ids.add(i)
                    break
        return list(all_ids)


    def weight(self):
        """
        Returns the number of edges or total of all edge weights
        """
        return self.size(weight="weight")

    def _get_node(self, nodeId:int):
        return PSObject({k:v for k,v in self.nodes[nodeId].items()})

    def urn2node(self, urn:str):
        try:
            return self.urn2obj[urn]
        except KeyError:
            raise KeyError

    def _get_nodes(self, node_ids=[]):
        if not node_ids:
            return [PSObject(ddict) for i,ddict in self.nodes(data=True)]
        else:
            nodes = list()
            for i in node_ids:
                nodes.append(PSObject(self.nodes(data=True)[i]))
            return nodes
            #return [PSObject(ddict) for id,ddict in self.nodes(data=True) if id in in_ids]

    def get_obj_by_prop(self,prop_value:str, prop_name='Name'): # use this function to find object by Name
        obj_ids = [i for i,o in self.nodes(data=True) if prop_value in o[prop_name]]
        return self._get_nodes(obj_ids)


    def get_objects(self, SearchValues: list, search_by_properties: list=None):
        if search_by_properties is None: search_by_properties = ['ObjTypeName']
        all_objects = set()
        for i, node in self.nodes(data=True):
            for propName in search_by_properties:
                ps_obj = PSObject(node)
                if ps_obj.is_annotated(propName, SearchValues):
                    all_objects.add(ps_obj)
                    break
        return list(all_objects)


    def __relations(self):
        """
        Returns list of all relation in the graph [PSRelation]
        """
        return list(r for n1,n2,r in self.edges.data('relation'))


    def _relation4(self, regulator_id, target_id):
        return [e['relation'] for i,e in dict(self[regulator_id][target_id]).items()]


    def get_relations(self, with_values:list, in_properties=[]):
        if not in_properties: in_properties = ['ObjTypeName']
        relations2return = set()
        for regulatorID, targetID, rel in self.edges.data('relation'):
            for prop_name in in_properties:
                if rel.has_property(prop_name, with_values):
                    relations2return.add(rel)
                    break
        return relations2return

    
    def find_relations(self, reg_id, targ_id, rel_type, effect='', mechanism=''):
        all_relations = self._relation4(reg_id, targ_id)
        my_rel = [rel for rel in all_relations if rel['ObjTypeName'][0]==rel_type]

        if effect:
            my_rel = [x for x in my_rel if x['Effect'] == effect]
        
        if mechanism:
            my_rel = [x for x in my_rel if x['Mechanism'] == mechanism]

        return my_rel


    def get_prop2obj_dic(self, search_by_property:str, filter_by_values=[], case_insensitive=False):
        search_value2obj = dict()
        objid2search_values = dict()
        if filter_by_values:
            allowed_values = list(map(lambda x:str(x).lower(),filter_by_values)) if case_insensitive else filter_by_values
            for nodeid, n in self.nodes(data=True):
                try:
                    node_prop_values = list(map(lambda x:str(x).lower(),n[search_by_property])) if case_insensitive else n[search_by_property] 
                    matched_values = [x for x in allowed_values if x in node_prop_values]
                    if matched_values:
                        objid2search_values[nodeid] = matched_values
                        for v in matched_values:
                            try:
                                search_value2obj[v].append(PSObject(n))
                            except KeyError:
                                search_value2obj[v] = [PSObject(n)]
                except KeyError: continue
        else:
            for nodeid, n in self.nodes(data=True):
                try:
                    all_values = list(map(lambda x: str(x).lower(),n[search_by_property])) if case_insensitive else n[search_by_property] 
                    objid2search_values[nodeid] = all_values
                    for v in all_values:
                        try:
                            search_value2obj[v].append(PSObject(n))
                        except KeyError:
                            search_value2obj[v] = [PSObject(n)]
                except KeyError: continue
        return search_value2obj, objid2search_values


    def get_props2obj_dic(self, propValues:list, prop_names:list, case_insensitive=False):
        propval2objs = dict()
        objid2propval = dict()
        for prop_name in prop_names:
            p2o, i2p = self.get_prop2obj_dic(prop_name, propValues,case_insensitive)
            for p,objs in p2o.items():
                try:
                    mapped_objs = set(propval2objs[p])
                    mapped_objs.update(objs)
                    propval2objs[p] = list(mapped_objs)
                except KeyError:
                    propval2objs[p] = objs

            for id,prop_vals in i2p.items():
                try:
                    mapped_values = set(objid2propval[id])
                    mapped_values.update(prop_vals)
                    objid2propval[id] = list(mapped_values)
                except KeyError:
                    objid2propval[id] = prop_vals

        return propval2objs, objid2propval


    def rel4pair(self, entity1prop, entity2prop, search_property='Name', filter_rel_with:dict={}, with_children=False):
        """
            entity1prop - property value to find first entity
            entity2prop - property value to find second entity
            search_property - property name to find pair entities by entity1prop,entity2prop
            filter_rel_with = {propName:[values]}
            if with_children - entities must have property CHILD to include ontology children

            Output - [PSRelation]
        """
        # both dicts must be {propName:[values]}
        entities1 = self.get_objects([entity1prop],search_by_properties=[search_property])
        entities2 = self.get_objects([entity2prop],search_by_properties=[search_property])
        entity1ids = [x['Id'][0] for x in entities1]
        entity2ids = [x['Id'][0] for x in entities2]
        if with_children:
            for e1 in entities1:
                try:
                    entity1ids = entity1ids+e1[CHILDS]
                except KeyError:
                    pass

            for e2 in entities2:
                try:
                    entity2ids = entity2ids+e2[CHILDS]
                except KeyError:
                    pass

        ids = set(entity1ids+entity2ids)
        if filter_rel_with:
            return [rel for regulatorID, targetID, rel in self.edges.data('relation') 
                if regulatorID in ids and targetID in ids and rel.has_value_in(filter_rel_with)]
        else:
            return [rel for regulatorID, targetID, rel in self.edges.data('relation') 
                if regulatorID in ids and targetID in ids]


    def recent_refs(self, entity1, entity2, search_property='Name', ref_limit=5,with_children=False):
        relations = self.rel4pair(entity1, entity2, search_property,with_children=with_children)
        references = set()
        for r in relations:
            references.update(list(r.References.values()))

        references = list(references)

        def sortkey(ref:dict):
            try:
                year = ref[PUBYEAR][0]
            except KeyError:
                try:
                    year = ref['Start'][0]
                    year = year[-4:]
                except KeyError: year = '1812'
            return tuple([year] + list(ref.Identifiers.values()))

        references.sort(key=lambda x: sortkey(x), reverse=True)
        total_refs = len(references)
        recent_refs = references[:ref_limit] if ref_limit else references
        ref_ids = dict()
        for ref in recent_refs:
            id_type,identifier = ref._identifier()
            try:
                ref_ids[id_type].append(identifier)
            except KeyError:
                ref_ids[id_type] = [identifier]
        return total_refs, ref_ids, recent_refs
    

    def _relations_ids(self):
        return {rel['Id'][0] for regulatorID, targetID, rel in self.edges.data('relation')}


    def find_nodes(self, for_relation:PSRelation, filter_by:list=None, in_properties:list=None):
        regulators_ids = for_relation.get_regulator_ids()
        regulators = [PSObject(o) for i,o in self.nodes(data=True) if i in regulators_ids]#self.get_objects(regulators_ids,['Id'])
        target_ids = for_relation.get_target_ids()
        targets = [PSObject(o) for i,o in self.nodes(data=True) if i in target_ids]
        #self.get_objects(target_ids,['Id'])

        if isinstance(filter_by,list):
            search_in_properties = in_properties if isinstance(in_properties,list) else ['ObjTypeName']
            must_values = set(filter_by)
            filtered_regulators = set()
            filtered_targets = set()
            for prop_name in search_in_properties:
                filtered_regulators.update([x for x in regulators if not set(x[prop_name]).isdisjoint(must_values)]) 
                filtered_targets.update([x for x in targets if not set(x[prop_name]).isdisjoint(must_values)]) 
            return list(filtered_regulators), list(filtered_targets)
        else:
            return regulators, targets

    
    def get_properties(self, ids:set, property_name:str):
        """
        # returns id2props = {id:[prop_values]}
        """
        id2props = {x: y[property_name] for x, y in self.nodes(data=True) if x in ids}
        return id2props


 


    def get_neighbors(self, node_ids:set, only_neighbors_with_ids=[]):
        """
        returns list of neighbor IDs
        """
        neighbor_ids = set()
        for i in [n for n in node_ids if self.has_node(n)]:
            neighbor_ids.update(set([x for x in nx.all_neighbors(self, i)]))
                
        if only_neighbors_with_ids:
            return [i for i in list(neighbor_ids) if i in only_neighbors_with_ids]
        else:
            return list(neighbor_ids)


    def get_neighbors_graph(self, for_node_ids:set, only_neighbors_with_ids=None):
        neighbors_ids = self.get_neighbors(for_node_ids,only_neighbors_with_ids)
        return self.get_subgraph(for_node_ids,neighbors_ids)


    def get_neighbors_rels(self, for_node_ids:set, only_neighbors_with_ids=None):
        neighbor_graph = self.get_neighbors_graph(for_node_ids,only_neighbors_with_ids)
        return neighbor_graph.__relations()


    def get_neighbors_refs(self, for_node_ids:set, only_neighbors_with_ids=None):
        neighbor_graph = self.get_neighbors_graph(for_node_ids,only_neighbors_with_ids)
        return neighbor_graph.load_references()

    @staticmethod
    def get_att_set(prop_name:str,ps_objects:list):
        return set([i for sublist in [x[prop_name] for x in ps_objects] for i in sublist])

    def get_children_ids(self, for_node_ids:list, at_depth=1):
        """
        needs node['Child Ids'] annotation
        use after APISessions::load_children
        """
        children_ids = set()
        level_parent_ids = for_node_ids
        for level in range(0,at_depth):
            level_children_ids = set()
            [level_children_ids.update(p[CHILDS]) for p in self._get_nodes(level_parent_ids)]
            level_parent_ids = level_children_ids
            children_ids.update(level_children_ids)

        return children_ids


    def __get_rels_between(self,node1_id, node2_id, from_relation_types=[]):
        try:
            edges = dict(self[node1_id][node2_id])
            if from_relation_types:
                return [v['relation'] for v in edges.values() if v['relation']['ObjTypeName'][0] in from_relation_types]
            else:
                return [r['relation'] for r in edges.values()]
        except KeyError:
            return list()


    def get_regulators(self, only_objtype=[], min_targets=1):
        if only_objtype:
            return {x for x, y in self.nodes(data=True) if ((self.out_degree(x) >= min_targets) & (y['ObjTypeName'][0] in only_objtype))}
        else:
            return {x for x in self.nodes() if self.out_degree(x) >= min_targets}


    def get_targets(self, only_objtype=[], min_regulators=1):
        if only_objtype:
            return {x for x, y in self.nodes(data=True) if ((self.in_degree(x) >= min_regulators) & (y['ObjTypeName'][0] in only_objtype))}
        else:
            return {x for x in self.nodes() if self.in_degree(x) >= min_regulators}


#################################   WRITE-DUMP, WRITE-DUMP, WRITE-DUMP ##############################
    def print_triples(self, fileOut, relPropNames, access_mode='w', printHeader=True, add_entities=False, as1row=False):
        with open(fileOut, access_mode, encoding='utf-8') as f:
            if printHeader:
                header = '\t'.join(relPropNames) + '\t' + "Regulators Id" + '\t' + "Targets Id"
                f.write(header + '\n')

            for regulatorID, targetID, rel in self.edges.data('relation'):
                f.write(rel.triple2str(relPropNames,add_entities=add_entities, as1row=as1row))

    def __get_ref_list(self, relPropNames:list, entity_prop_names=[], RefNumPrintLimit=0,col_sep:str='\t',single_rel_row=False):
        references = list()
        if not entity_prop_names:
            for regulatorID, targetID, rel in self.edges.data('relation'):
                reference_view_triple = str(rel.triple2str(relPropNames,as1row=single_rel_row))
                references.append(reference_view_triple)
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
                                                        col_sep=col_sep)

                    reference_table_view = reference_table_view[0:len(reference_table_view) - 1]#remove end of line character
                    reference_table_view = reg_props_str[0:len(reg_props_str) - 1]+col_sep+reference_table_view
                    reference_table_view += col_sep + target_props_str
                else:
                    rel_props_str_list = dict(rel.to_table_dict(
                                                relPropNames, RefNumPrintLimit=RefNumPrintLimit))
                    
                    for row in rel_props_str_list.values():
                        reference_table_view += reg_props_str[:-1]
                        reference_table_view += col_sep + col_sep.join(row) + col_sep + target_props_str[:-1]
                        references.append(reference_table_view)

            return references


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


    def to_rnef(self, ent_props:list, rel_props:list, add_rel_props:dict=None,add_pathway_props:dict=None):
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
            if ent_props:
                for prop_name, prop_values in n.items():
                    if prop_name in ent_props:
                        for prop_value in prop_values:
                            et.SubElement(xml_node, 'attr', {'name': str(prop_name), 'value': str(prop_value)})
 

            local_id_counter += 1

        xml_controls = et.SubElement(resnet, 'controls')

        graph_relations = self.__relations()
        for rel in graph_relations:
            control_id = self.rel_urn(rel)
            xml_control = et.SubElement(xml_controls, 'control', {'local_id':control_id})
            et.SubElement(xml_control, 'attr', {'name':'ControlType', 'value':str(rel['ObjTypeName'][0])})
            
            # adding links
            regulators = rel.Nodes[REGULATORS]
            try:
                targets = rel.Nodes[TARGETS]
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

            # non-reference properties
            for prop_name, prop_values in rel.items():
                if prop_name in rel_props:
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
                ref_props = {k:v for k,v in ref.items() if k in rel_props}
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


    def to_rnef_file(self, fname:str, ent_prop2print:set, rel_prop2print:set, add_rel_props:dict=None,add_pathway_props:dict=None):
        rnef_str = self.to_rnef(ent_prop2print,rel_prop2print, add_rel_props,add_pathway_props)
        rnef_str = '<batch>\n'+str(rnef_str).strip()+'</batch>\n'
        pretty_xml = str(minidom.parseString(rnef_str).toprettyxml(indent='  '))
        with open(fname,'w',encoding='utf-8') as f: f.write(pretty_xml)


    def add_row2(self,to_df:pd.DataFrame,from_relation_types:list, between_node_id, and_node_id, from_properties:list,
                cell_sep=';'):
        """
        Appends row to_df generated from valid triples (between_node_id, and_node_id, rel) according to specifications in from_properties

        from_properties = [ 'N1:PropName', 'N2:PropName,'R:PropName']
        N1 prefix is for between_node_id properties, N2 prefix is for and_node_id properties, R prefix is for relation properties
        all properties will be added in the order of from_properties
        multiple values from property are joined by cell_sep=';'
        """
        node1 = self._get_node(between_node_id)
        node2 = self._get_node(and_node_id)
        relations = self.__get_rels_between(between_node_id,and_node_id,from_relation_types)
        relations += self.__get_rels_between(and_node_id,between_node_id,from_relation_types)

        row_template = ['']*len(from_properties)
        for i in range(0, len(from_properties)):
            if from_properties[i][0] == 'N':
                node4cell = node1 if from_properties[i][1] =='1' else node2
                cell = node4cell.prop_values(from_properties[i][3:],sep=cell_sep)
                row_template[i] = cell
    
        row_counter = 0
        for rel in relations:
            row = row_template
            for i in range(0, len(from_properties)):
                if from_properties[i][:2] == 'R:':
                    cell = rel._props2str(from_properties[i][2:])
                    row[i] = cell
                    to_df.loc[len(to_df.index)] = row
            row_counter += 1

        if not row_counter:
            to_df.loc[len(to_df.index)] = row_template


################################# READ READ READ ##########################################
    def _parse_nodes_controls(self, resnet:et.Element, new_node_count=1):
        nodel_local_ids = dict()
        new_node_id = new_node_count
        for node in resnet.findall('./nodes/node'):
            node_urn = node.get('urn')
            local_id = node.get('local_id')
            try:
                node_obj = self.urn2node(node_urn)
            except KeyError:
                node_obj = PSObject({'Id':[new_node_id],'URN':[node_urn]})
                new_node_id += 1

            [node_obj.update_with_value(attr.get('name'), attr.get('value')) for attr in node.findall('attr')]
            node_obj['ObjTypeName'] = node_obj.pop('NodeType')
            self.add1node(node_obj)
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
            effect_val = '' if type(effect) == type(None) else effect.text

            for reg in regulators:
                try: 
                    ps_rel.Nodes[REGULATORS].append((reg['Id'][0], '0', effect_val))
                except KeyError:
                    ps_rel.Nodes[REGULATORS] = [(reg['Id'][0], '0', effect_val)]
                
            for targ in targets:
                try: 
                    ps_rel.Nodes[TARGETS].append((targ['Id'][0], '0', effect_val))
                except KeyError:
                    ps_rel.Nodes[TARGETS] = [(targ['Id'][0], '0', effect_val)]
                        
            for attr in rel.findall('attr'):
                prop_id = attr.get('name')
                prop_value = attr.get('value')
            
                index = attr.get('index')
                if type(index) == type(None):
                    ps_rel.update_with_value(prop_id, prop_value)
                else:
                    try:
                        props = ps_rel.PropSetToProps[index]
                        try:
                            props[prop_id].append(prop_value)
                        except KeyError:
                                props[prop_id] = [prop_value]
                    except KeyError:
                        ps_rel.PropSetToProps[index] = {prop_id:[prop_value]}

            ps_rel['ObjTypeName'] = ps_rel.pop('ControlType')
            self.add_rel(ps_rel)
        return new_node_id


    def read_rnef(self, rnef_file:str, last_new_node_id=1):
        try:
            root = et.parse(rnef_file).getroot()
            print ('Loading graph from file %s' % rnef_file)
        except FileNotFoundError:
            raise FileNotFoundError

        for resnet in root.findall('resnet'):
            last_new_node_id = self._parse_nodes_controls(resnet,last_new_node_id)


    @classmethod
    def fromRNEF(cls,rnef_file:str):
        try:
            g = ResnetGraph()
            g.read_rnef(rnef_file)
            return g
        except FileNotFoundError:
            raise FileNotFoundError


########################  SUBGRAPH SUBGRAPH SUBGRAPH #####################################
    def subtract(self, other: "ResnetGraph"):
        #only self graph is analyzed 
        # works faster if other graph is bigger than self
        unique2self = ResnetGraph()
        edges_from_other = other._relations_ids()
        for n1,n2,e in self.edges(data=True):
            if e['relation']['Id'][0] not in edges_from_other:
                unique2self.add_edge(n1, n2, relation=e['relation'],weight=e['weight'], key=self.rel_urn(e['relation']))
                unique2self.add_nodes_from([(n1, self.nodes[n1]),(n2, self.nodes[n2])])
        return unique2self


    def intersect (self, other: "ResnetGraph"):
        #only self graph is analyzed 
        # works faster if other graph is bigger than self
        intersection = ResnetGraph()
        edges_from_other = other._relations_ids()
        for n1,n2,e in self.edges(data=True):
            if e['relation']['Id'][0] in edges_from_other:
                intersection.add_edge(n1, n2, relation=e['relation'],weight=e['weight'],key=self.rel_urn(e['relation']))
                intersection.add_nodes_from([(n1, self.nodes[n1]),(n2, self.nodes[n2])])
        return intersection


    def filter_references(self, keep_prop2values:dict, rel_types=[]):
        """
        # prop_names2values = {prop_name:[values]}
        """
        filtered_graph = self.copy()
        if rel_types:
            for regulatorID, targetID, rel in filtered_graph.edges.data('relation'):        
                if rel['ObjTypeName'][0] in rel_types:
                    rel.filter_references(keep_prop2values)
        else:
            for regulatorID, targetID, rel in filtered_graph.edges.data('relation'):
                rel.filter_references(keep_prop2values)

        return filtered_graph
  
                    
    def regulatory_network(self,for_targets_with_urns:list,network_name:str):
        """
        returns subgraph with targets from 'for_targets_with_urns' and their regulators 
        funtion is used by SNEA
        """
        start = time.time()
        targets_urns = set(for_targets_with_urns)
        input_target_ids = [node['Id'][0] for node in self.urn2obj.values() if node['URN'][0] in targets_urns]

        reg_subgraph = ResnetGraph()
        rel2copy = [e for r,t,e in self.edges(data='relation') if t in set(input_target_ids)]
        [reg_subgraph.copy_rel(rel,self) for rel in rel2copy]
            
        reg_subgraph.name = network_name
        print('%s subnetwork was selected in %s' % (network_name, self.execution_time(start)))
        return reg_subgraph


    def make_simple(self):
        """
        returns graph with only one relation between pair of nodes.
        keeps relation with biggest reference count.
        all other relations are merged into the most referenced one
        """
        simple_g = ResnetGraph()
        added_edges = set()
        for regulator_id, target_id, edges in self.edges(data=True):
            if (regulator_id, target_id) in added_edges: 
                continue
            reg2target_edges = dict(self[regulator_id][target_id])
            reg2target_rels = [e['relation'] for i,e in reg2target_edges.items()]
            if len(reg2target_rels) == 1:
                simple_g.copy_rel(reg2target_rels[0],self)
                added_edges.add((regulator_id, target_id))
            else:
                reg2target_rels.sort(key=lambda x: x[REFCOUNT][0], reverse=True)

                best_rel = PSRelation()
                best_rel_idx = 0
                for i in range(0, len(reg2target_rels)):
                    rel = reg2target_rels[i]
                    try:
                        eff = rel['Effect'][0]
                        if eff in ['positive', 'negative']:
                            best_rel = rel
                            best_rel_idx = i
                            break
                        else:
                            continue
                    except KeyError:
                        continue
                
                reg2target_rels.pop(best_rel_idx)            
                [best_rel.copy(rel) for rel in reg2target_rels]
                    
                simple_g.copy_rel(best_rel,self)
                added_edges.add((regulator_id, target_id))

        return simple_g
                
                    
    def regulome_dict(self, min_size=2):
        """
        returns {regulator_id: [targets]};
        len([targets]) >= min_size
        """
        regulator_ids = list(self.get_regulators(min_targets=min_size))
        subnetworks = dict()
        for regulator_id in regulator_ids:
            target_ids = list(self.neighbors(regulator_id))
            targets = self._get_nodes(target_ids)
            subnetworks[regulator_id] = targets
            nx.set_node_attributes(self, {regulator_id:{'# targets':[len(target_ids)]}})

        print('Generated %d regulome subnetworks with more than %d targets from %s' 
                        % (len(subnetworks), min_size, self.name))
        return subnetworks


    def downstream_relations(self,node_id:int):
        return [rel for r,t,rel in self.edges.data('relation') if r == node_id]


    def get_regulome(self, start_node_ids: set):
        """
        returns composition of bfs_trees for all ids in start_node_ids
        """
        all_trees = nx.DiGraph()
        for Id in start_node_ids:
            if not self.has_node(Id): continue
            t = nx.bfs_tree(self, Id)
            all_trees = nx.compose(all_trees, t)
        return all_trees

    def subgraph_by_rel_urns(self, rel_urns:set):
        """
        #not tested!!
        """
        subgraph = ResnetGraph()
        rels2add = [r for u,r in self.urn2rel.items() if u in rel_urns]
        [subgraph.copy_rel(rel,self) for rel in rels2add]
        return subgraph

    def subgraph_by_rel_ids(self, rels:list):
        """
        rels = [PSRelation]
        """
        subgraph = ResnetGraph()
        [subgraph.copy_rel(rel,self) for rel in rels]
        return subgraph


    def subgraph_by_relprops(self, search_values:list, in_properties:list=None):
        if in_properties is None: in_properties = ['ObjTypeName']
        search_value_set = set(search_values)
        subgraph = ResnetGraph()
        for prop_type in in_properties:
            for n1, n2, rel in self.edges.data('relation'):
                if not set(rel[prop_type]).isdisjoint(search_value_set):
                    subgraph.add_edge(n2, n1, relation=rel,weight=rel.get_reference_count(),key=self.rel_urn(rel))
                    subgraph.add_nodes_from([(n1, self.nodes[n1]),(n2, self.nodes[n2])])
        return subgraph


    def get_subgraph(self,between_node_ids:list,and_node_ids:list,by_relation_type=None,with_effect=None,in_direction=None)->"ResnetGraph":
        subgraph = ResnetGraph()
        for n1 in between_node_ids:
            for n2 in and_node_ids:
                if not isinstance(in_direction,str) or in_direction == '>':
                    if self.has_edge(n1, n2):
                        for i in range(0, len(self[n1][n2])): #gives error if edge does not exist
                            rel = self[n1][n2][i]['relation']
                            if isinstance(by_relation_type,list) and rel['ObjTypeName'] not in by_relation_type: continue
                            if isinstance(with_effect,list):
                                try: ef = rel['Effect'] 
                                except KeyError:
                                    rel['Effect'] = 'unknown'
                                    ef = rel['Effect']
                                if ef not in with_effect: continue
                                
                            subgraph.add_edge(n1, n2, relation=rel,weight=self[n1][n2][i]['weight'],key=self.rel_urn(rel))
                            subgraph.add_nodes_from([(n1, self.nodes[n1]),(n2, self.nodes[n2])])
                if not isinstance(in_direction,str) or in_direction == '<':
                    if self.has_edge(n2, n1):
                        for i in range(0, len(self[n2][n1])):#gives error if edge does not exist
                            rel = self[n2][n1][i]['relation']
                            if isinstance(by_relation_type,list) and rel['ObjTypeName'] not in by_relation_type: continue
                            if isinstance(with_effect,list):
                                try: ef = rel['Effect'] 
                                except KeyError:
                                    rel['Effect'] = 'unknown'
                                    ef = rel['Effect']
                                if ef not in with_effect: continue
                                    
                            subgraph.add_edge(n2, n1, relation=rel,weight=self[n2][n1][i]['weight'],key=self.rel_urn(rel))
                            subgraph.add_nodes_from([(n1, self.nodes[n1]),(n2, self.nodes[n2])])
            
        return subgraph
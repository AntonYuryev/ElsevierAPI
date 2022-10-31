import networkx as nx
from pip import List
from ..pandas.panda_tricks import df
from datetime import timedelta
import os, math, time
from xml.dom import minidom
import xml.etree.ElementTree as et

from .NetworkxObjects import PSObject,PSRelation,len,REGULATORS,TARGETS,CHILDS,REFCOUNT
from ..ETM_API.references import PUBYEAR,EFFECT,TITLE,pubmed_hyperlink, make_hyperlink
from ..ETM_API.etm import ETMstat
from itertools import combinations

RESNET = 'resnet'
PHYSICAL_INTERACTIONS = ['Binding','DirectRegulation','ProtModification','PromoterBinding','ChemicalReaction']
PROTEIN_TYPES = ['Protein','FunctionalClass','Complex']

#NO_RNEF_NODE_PROPS = {'Id','URN','ObjClassId','ObjTypeId','ObjTypeName','OwnerId','DateCreated','DateModified'}
#NO_RNEF_REL_PROPS = NO_RNEF_NODE_PROPS | {'RelationNumberOfReferences', '# of Total References', 'Name'}

class ResnetGraph (nx.MultiDiGraph):
    pass
    def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.urn2obj = dict() #lookup for combining nodes in database graph and rnef graph
            self.urn2rel = dict() #lookup for combining relations in database graph and rnef graph


    def copy(self)->'ResnetGraph':
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


    def add_nodes(self,nodes:list or dict):
        """
        Input
        -----
        nodes: [PSObject] or {id:PSObject}
        """
        if isinstance(nodes,dict):
            self.add_nodes_from([(k,v.items()) for k,v in nodes.items()])
            new_urn_dict = {n['URN'][0]:PSObject(n) for n in nodes.values()}
            self.urn2obj.update(new_urn_dict)
        else:
            self.add_nodes_from([(int(n['Id'][0]), n.items()) for n in nodes])
            new_urn_dict = {n['URN'][0]:PSObject(n) for n in nodes}
            self.urn2obj.update(new_urn_dict)


    def property2node(self, node_id,prop_name:str,prop_values:list):
        nx.set_node_attributes(self,{node_id:{prop_name:prop_values}})


    def copy_node_annotation(self, from_:'ResnetGraph', prop_name:str):
        annotations = {i:r for i,r in from_.nodes(data=prop_name)}
        nx.set_node_attributes(self,annotations,prop_name)

    def __add_rel(self, rel:PSRelation):
        """
        # nodes connected by rel must exist in the graph
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
            self.urn2rel[rel_urn].merge_rel(rel)
            self.__add_rel(self.urn2rel[rel_urn])
        except KeyError:
            self.urn2rel[rel_urn] = rel
            self.__add_rel(rel)
            

    def copy_rel(self,rel:PSRelation,from_graph:"ResnetGraph"):
        ids4nodes2add = rel.entities_ids()
        nodes2add = list()
        for node_id in ids4nodes2add:
            node2add = from_graph._get_node(node_id)
            try:
                self.urn2obj[node2add.urn()]
            except KeyError:
                nodes2add.append((node_id,node2add.items()))
                self.urn2obj[node2add.urn()] = node2add

        self.add_nodes_from(nodes2add)
        self.__add_rel(rel)

    
    def add_triple(self,regulator:PSObject,target:PSObject,rel:dict or PSRelation,refs=[],is_directional=True):
        """
        adds nodes and their relations to graph
        """
        self.add_nodes([regulator,target])

        if isinstance(rel,PSRelation):
            self.add_rel(rel)
            return
        else:
            rel = PSRelation.make_rel(regulator,target,rel,refs,is_directional)
            rel['Id'] = [self.number_of_edges()]
            self.rel_name(rel)
            rel_urn = self.rel_urn(rel)
            self.add_edge(regulator['Id'][0], target['Id'][0], relation=rel,weight=rel.get_reference_count(),key=rel_urn)


    def add_graph(self, other:"ResnetGraph"):
        self.update(other)
        self.urn2obj.update(other.urn2obj)
        self.urn2rel.update(other.urn2rel)
        self.references = dict()

################## SET SET SET ##########################################
    def set_node_annotation(self, urn2value:dict, new_prop_name:str):
        """
            adds new property to exsisting nodes. Existing values of 'new_prop_name' will be replaced
            urn2value = {urn:[with_prop_values]}
        """
        id2value = dict() #remapping urn to ids
        for urn, value in urn2value.items():
            try:
                node2annotate = self.urn2node(urn)
                id2value[node2annotate['Id'][0]] = value
            except KeyError: continue

        nx.set_node_attributes(self,id2value,new_prop_name)
        #print('%d nodes were annotated with attributes "%s"' % (len(id2value),new_prop_name))


    def set_edge_property(self,nodeId1,nodeId2,rel_urn:str,prop_name,prop_values:list):
        annotate_values = set(prop_values)
        if self.has_edge(nodeId1, nodeId2,rel_urn):
            try:
                my_prop_values = self[nodeId1][nodeId2][rel_urn]['relation'][prop_name]
                annotate_values.update(my_prop_values)
            except KeyError:
                pass
            self[nodeId1][nodeId2][rel_urn]['relation'][prop_name] = list(annotate_values)


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


    def copy_node_annotation(self, from_property:str, in_other:"ResnetGraph", as_propname=''):
        urn2value = {n.urn():n[from_property] for i,n in in_other.nodes(data=True)}
        prop_name = as_propname if as_propname else from_property
        self.set_node_annotation(urn2value, prop_name)


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


    def closeness(self):
        """
        Annotates
        ---------
        nodes with property 'Closeness' calculated by nx.harmonic_centrality -\n 
        average length of the shortest path between the node and all other nodes in the graph
        Returns
        -------
        {node_id:Closeness}
        """
        centrality_dic = nx.harmonic_centrality(self)
        norm_const = self.number_of_nodes()-1
        centrality_dic_norm = {i:v/norm_const for i,v in centrality_dic.items()}
        nx.set_node_attributes(self, centrality_dic_norm, 'Closeness')
        return centrality_dic_norm


    def rank_regulator(self, regulator_id:int, target_weights:dict, max_distance=5):
        '''
        Input
        -----
        target_weights = {node_id:weight}
        Returns
        -------
        regulator rank for regulator_id
        '''
        regulator_rank = 0
        regulator_tree = nx.bfs_tree(self, regulator_id)

        for level in range (1,max_distance+1):
            targets_on_level = nx.descendants_at_distance(regulator_tree, regulator_id, level)
            if not targets_on_level: break
            for target_id in targets_on_level:
                try:
                    target_weight = target_weights[target_id]
                    regulator_rank = regulator_rank + target_weight/math.pow(level,2)
                except KeyError:
                    continue
        return regulator_rank


    def rank_regulators(self, node_weights:dict, add2prop:str, max_distance=5):
        '''
        Input
        -----
        node_weights = {node_id:weight}
        Returns
        -------
        {regulator_id:rank}
        '''
        regulator_ranks = dict()
        regulator_trees = nx.DiGraph()
        for Id in node_weights.keys():
            if self.has_node(Id):
                source_weight = node_weights[Id]
                tree = nx.bfs_tree(self, Id, reverse=True)
                regulator_trees = nx.compose(regulator_trees,tree)
                scored_regulators = set()
                for level in range (1,max_distance+1):
                    regulators = nx.descendants_at_distance(tree, Id, level)
                    regulators_on_level = regulators.difference(scored_regulators)
                    for regulator_id in regulators_on_level:
                        try:
                            rank = float(regulator_ranks[regulator_id])
                            regulator_ranks[regulator_id] = rank + source_weight/math.pow(level,2)
                        except KeyError:
                            regulator_ranks[regulator_id] = source_weight/math.pow(level,2)
                    scored_regulators.update(regulators)
                    
        nx.set_node_attributes(self,regulator_ranks,add2prop)

        #network adjustment to boost regulators regulating other regulators on the same level of regulator_trees
        for regulator_id in self.nodes():
            my_target_ids = self.neighbors(regulator_id)
            neighborhood_weight = 0.0
            for target_id in my_target_ids:
                if regulator_trees.has_edge(target_id,regulator_id):
                    continue #if edge was visited
                    # regulator_trees is reversed 
                target = self.nodes[target_id]
                try:
                    target_weight = target[add2prop]
                    neighborhood_weight += target_weight
                except KeyError:
                    continue
            try:
                node_rank = regulator_ranks[regulator_id]
                regulator_ranks[regulator_id] = node_rank + 0.5*neighborhood_weight
            except KeyError:
                # this should not happen
                regulator_ranks[regulator_id] = 0.5*neighborhood_weight

        nx.set_node_attributes(self,regulator_ranks,add2prop)
        return regulator_ranks


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
        Input
        -----
        references are annotated with 'weight' property according to specification in proval2weight
        """
        graph_references = set()
        if isinstance(weight_by_prop_name,str):
            for regulatorID, targetID, rel in self.edges.data('relation'):
                rel.load_references()
                rel._weight2ref(weight_by_prop_name,proval2weight)
                graph_references.update(rel.References.values())
        else:
            for regulatorID, targetID, rel in self.edges.data('relation'):
                rel.load_references()
                graph_references.update(rel.References.values())

        return graph_references


    def citation_index(self)->dict:
        """
        Returns
        -------
        dictionary {identifier:Reference} with all references in graph.\n 
        Reference objects in returned dictionary are annotated by Reference['Citation index'] property
        """
        graph_references = dict()
        for regulatorID, targetID, rel in self.edges.data('relation'):
            rel_refs = rel._get_refs()
            for ref in rel_refs:
                id_type, ref_id = ref.get_doc_id()
                if ref_id:
                    try:
                        ref_counter = graph_references[id_type+':'+ref_id]
                        count = ref_counter['Citation index'][0]
                        ref_counter['Citation index'] = [count + 1]
                    except KeyError:
                        ref['Citation index'] = [1]
                        graph_references[id_type+':'+ref_id] = ref

        return graph_references


    def refs2df(self,df_name:str,for_rel_types=list(),with_effect=list()):
        sub_graph = self
        if for_rel_types:
            sub_graph = sub_graph.subgraph_by_relprops(for_rel_types)
        if with_effect:
            sub_graph = sub_graph.subgraph_by_relprops(with_effect,['Effect'])
        
        references = sub_graph.citation_index() # annotates 
        ref_df = ETMstat.external_counter2pd(set(references.values()),stat_prop='Citation index')
        ref_df._name_ = df_name
        return ref_df


    def add_recent_refs(self,to_df:df,between_col:str,and_col:str,map2prop='Name'):
        annoate_df = df.copy_df(to_df)
        for row in to_df.index:
            node1_prop = annoate_df.loc[row][between_col]
            node2_prop = annoate_df.loc[row][and_col]
            total_refs, ref_ids, ps_references = self.recent_refs(node1_prop,node2_prop,map2prop,with_children=True)
            dois = str()
            try:
                pmids = ref_ids['PMID']
                refcount = pubmed_hyperlink([pmids],total_refs)
            except KeyError:
                try:
                    dois = ref_ids['DOI']
                    refcount = make_hyperlink(dois[0],url='http://dx.doi.org/', display_str=str(total_refs))
                except KeyError:
                    refcount = str(total_refs)
            annoate_df.loc[row]['Recent pubs'] = refcount
        return annoate_df


    def load_references_between(self, between_node_ids: list, and_node_ids: list, weight_by_prop_name:str=None, proval2weight:dict=None):
        '''
        Input
        -----
        proval2weight = {value:weight} for values in weight_by_prop_name
        Returns
        -------
        {References}
        '''
        sub_graph = self.get_subgraph(between_node_ids, and_node_ids)
        return sub_graph.load_references(weight_by_prop_name,proval2weight)


    def snippets2df(self, df_name='Snippets'):
        annotated_refs = self.citation_index()
        header = ['Concept','Entity','Citation index','PMID','DOI',PUBYEAR,TITLE,'Snippets']
        row_tuples = set()

        for regulatorID, targetID, rel in self.edges.data('relation'):
            regulator_name = self.nodes[regulatorID]['Name'][0]
            target_name = self.nodes[targetID]['Name'][0]
            rel_refs = rel._get_refs()
            for ref in rel_refs:
                id_type, identifier = ref.get_doc_id()
                try:
                    annotated_ref = annotated_refs[id_type+':'+identifier]
                    ref_list = annotated_ref.to_list(['PMID','DOI'],True,[PUBYEAR,TITLE],['Citation index'],True)
                    row_tuples.add(tuple([regulator_name,target_name] + ref_list))
                except KeyError:
                    continue
                    
        snippet_df = df.from_rows(row_tuples,header)
        snippet_df.sort_values(['Concept','Entity','Citation index'], inplace=True, ascending=False)
        snippet_df._name_ = df_name
        return snippet_df

##################### DEL DEL DEL ######################################
    def remove_nodes_by_prop(self, property_values:list, prop_names:list=['ObjTypeName']):
        node_ids = self.get_node_ids(property_values,prop_names)
        self.remove_nodes_from(node_ids)
        print("%d nodes with %s were removed" % (len(node_ids), ','.join(property_values)))


    def remove_nodes_by_degree(self,min_degree=0,max_degree=1000000,only_with_prop:list=['ObjTypeName'],having_values:list=[]):
        if having_values:
            ids2remove = set(self.get_node_ids(having_values,only_with_prop))
            ids2remove = {x for x in self.nodes() if self.degree(x) > max_degree and x in ids2remove}
        else:
            ids2remove = {x for x in self.nodes() if self.degree(x) > max_degree}

        if min_degree:
            ids2remove = {x for x in self.nodes() if self.degree(x) < min_degree and x in ids2remove}
        self.remove_nodes_from(ids2remove)
        if min_degree:
            print('%d nodes with %d < degree < %d were removed' % (len(ids2remove),min_degree,max_degree))
        else:
            print('%d nodes with degree > %d were removed' % (len(ids2remove),max_degree))


    def remove_nodes_by_outdegree(self,min_degree=0,max_degree=1000000,only_with_prop:list=['ObjTypeName'],having_values:list=[]):
        if having_values:
            ids2remove = set(self.get_node_ids(having_values,only_with_prop))
            ids2remove = {x for x in self.nodes() if self.out_degree(x) > max_degree and x in ids2remove}
        else:
            ids2remove = {x for x in self.nodes() if self.out_degree(x) > max_degree}

        if min_degree:
            ids2remove = {x for x in self.nodes() if self.out_degree(x) < min_degree and x in ids2remove}
        self.remove_nodes_from(ids2remove)
        if min_degree:
            print('%d nodes with %d < outdegree < %d were removed' % (len(ids2remove),min_degree,max_degree))
        else:
            print('%d nodes with outdegree > %d were removed' % (len(ids2remove),max_degree))

        
    def remove_nodes_by_indegree(self,min_degree=0,max_degree=1000000,only_with_prop:list=['ObjTypeName'],having_values:list=[]):
        if having_values:
            ids2remove = set(self.get_node_ids(having_values,only_with_prop))
            ids2remove = {x for x in self.nodes() if self.in_degree(x) > max_degree and x in ids2remove}
        else:
            ids2remove = {x for x in self.nodes() if self.in_degree(x) > max_degree}

        if min_degree:
            ids2remove = {x for x in self.nodes() if self.in_degree(x) < min_degree and x in ids2remove}
        self.remove_nodes_from(ids2remove)
        if min_degree:
            print('%d nodes with %d < indegree < %d were removed' % (len(ids2remove),min_degree,max_degree))
        else:
            print('%d nodes with indegree > %d were removed' % (len(ids2remove),max_degree))


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


    def clear_resnetgraph(self):
        super().clear()
        self.urn2obj.clear()
        self.urn2rel.clear()


######################   GET GET GET   ######################################
    def get_node_ids(self, SearchValues:list, search_by_properties:list=['ObjTypeName']):
        all_ids = set()
        for i, node in self.nodes(data=True):
            for propName in search_by_properties:
                if PSObject(node).is_annotated(propName, SearchValues):
                    all_ids.add(i)
                    break
        return list(all_ids)


    def node_id2urn(self, SearchValues=list(), search_by_properties:list=['ObjTypeName']):
        all_ids = dict()
        if SearchValues:
            for i, node in self.nodes(data=True):
                for propName in search_by_properties:
                    ps_obj = PSObject(node)
                    if ps_obj.is_annotated(propName, SearchValues):
                        all_ids[i] = ps_obj.urn()
                        break
            return all_ids
        else:
            return {i:u[0] for i,u in self.nodes(data='URN')}

    def weight(self):
        """
        Returns
        -------
        the number of edges or total of all edge weights
        """
        return self.size(weight="weight")


    def _get_node(self, nodeId:int):
        return PSObject({k:v for k,v in self.nodes[nodeId].items()})


    def urn2node(self, urn:str):
        try:
            return PSObject(self.urn2obj[urn])
        except KeyError:
            raise KeyError


    def _get_nodes(self, node_ids=[]):
        if not node_ids:
            return [PSObject(ddict) for i,ddict in self.nodes(data=True)]
        else:
            id_set = set(node_ids)
            return [PSObject(ddict) for i,ddict in self.nodes(data=True) if i in id_set]

            #may be faster?
            '''
            nodes = list()
            for i in node_ids:
                try:
                    node = self.get_objects()
                    node = PSObject(self.nodes(data=True)[i])
                    nodes.append(PSObject(node))
                except KeyError:
                    continue
            return nodes
            '''

    def get_obj_by_prop(self,prop_value:str, prop_name='Name'): # use this function to find object by Name
        obj_ids = [i for i,o in self.nodes(data=True) if prop_value in o[prop_name]]
        return self._get_nodes(obj_ids)


    def get_objects(self, SearchValues:list, search_by_properties=['ObjTypeName']):
        all_objects = set()
        for i, node in self.nodes(data=True):
            for propName in search_by_properties:
                ps_obj = PSObject(node)
                if ps_obj.is_annotated(propName, SearchValues):
                    all_objects.add(ps_obj)
                    break
        return list(all_objects)


    def _relations(self):
        """
        Returns
        -------
        list of all relations in graph, [PSRelation]
        """
        return list(r for n1,n2,r in self.edges.data('relation'))


    def _relation4(self, regulator_id, target_id):
        try:
            return [e['relation'] for i,e in dict(self[regulator_id][target_id]).items()]
        except KeyError:
            return []


    def get_relations(self, with_values:list, in_properties:list=['ObjTypeName']):
        relations2return = set()
        for regulatorID, targetID, rel in self.edges.data('relation'):
            for prop_name in in_properties:
                if rel.is_annotated(prop_name, with_values):
                    relations2return.add(rel)
                    break
        return relations2return

    
    def __find_relations(self, reg_id, targ_id, rel_types=list(), effect=list(), mechanism=list(), any_direction=False):
        my_rels = self._relation4(reg_id, targ_id)
        if any_direction:
            my_rels = my_rels + self._relation4(targ_id, reg_id)

        if rel_types:
            my_rels = [rel for rel in my_rels if rel['ObjTypeName'][0] in rel_types]

        if effect:
            try:
                my_rels = [x for x in my_rels if x['Effect'][0] in effect]
            except KeyError:
                return list()
        
        if mechanism:
            try:
                my_rels = [x for x in my_rels if x['Mechanism'][0] in mechanism]
            except KeyError:
                return list()

        return my_rels


    def relation_exist(self,between_ids:list,and_ids:list,with_reltypes=list(),with_effects=list(),mechanism=list(),any_direction=False):
        for node_id1 in between_ids:
            for node_id2 in and_ids:
                my_rels = self.__find_relations(node_id1,node_id2,with_reltypes,with_effects,mechanism,any_direction)
                if my_rels:
                    return True
                else:
                    continue

        return False


    def find_relations(self, between_ids,and_ids,with_reltypes=list(),with_effects=list(),mechanism=list(),any_direction=False):
        my_rels = set()
        for node_id1 in between_ids:
            for node_id2 in and_ids:
                my_rels.update(self.__find_relations(node_id1,node_id2,with_reltypes,with_effects,mechanism,any_direction))
          
        return my_rels


    def _effect_counts__(self):
        positive_rel = self.get_relations(['positive'],[EFFECT])
        negative_rel = self.get_relations(['negative'],[EFFECT])

        positive_refs = set()
        [positive_refs.update(rel._get_refs()) for rel in positive_rel]
        negative_refs = set()
        [negative_refs.update(rel._get_refs()) for rel in negative_rel]

        return positive_refs,negative_refs


    def effect_stats(self, between_node_ids:list, and_node_ids:list):
        between_graph = self.get_subgraph(between_node_ids, and_node_ids)
        positive_refs, negative_refs = between_graph._effect_counts__()
        return between_graph,positive_refs,negative_refs


    def get_prop2obj_dic(self, search_by_property:str, filter_by_values=[], case_insensitive=False):
        '''
        Returns
        -------
        propval2objs = {search_by_property_value:[PSObject]}
        objid2propval = {id:[search_by_property_value]}
        '''
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


    def rel4pair(self, entity1prop:str, entity2prop:str, search_property='Name', filter_rel_with:dict={}, with_children=False):
        """
        Input
        -----
        entity1prop - property value to find first entity
        entity2prop - property value to find second entity
        search_property - property name to find pair entities by entity1prop,entity2prop
        filter_rel_with = {propName:[values]}
        if with_children - entities must have property CHILD to include ontology children

        Returns
        -------
        list of PSRelation objects
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


    def recent_refs(self, entity1prop:str, entity2prop:str, search_property='Name', ref_limit=5,with_children=False):
        """
        Returns
        -------
        total_refs, ref_ids={id_type:[identifiers]}, recent_refs={Reference}
        """
        relations = self.rel4pair(entity1prop, entity2prop, search_property,with_children=with_children)
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


    def find_nodes(self,for_relation:PSRelation,filter_by:list=None,in_properties:list=None):
        """
        Returns
        -------
        regulators = [PSObject], targets = - [PSObject]
        """
        regulators_ids = for_relation.regulator_ids()
        regulators = [PSObject(o) for i,o in self.nodes(data=True) if i in regulators_ids]#self.get_objects(regulators_ids,['Id'])
        target_ids = for_relation.target_ids()
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

    
    def get_properties(self, property_name:str, for_ids=set()):
        """
        Returns
        -------
        id2props = {id:[prop_values]}
        """
        if for_ids:
            return {x: y[property_name] for x, y in self.nodes(data=True) if x in for_ids}
        else:
            return {x: y[property_name] for x, y in self.nodes(data=True)}


    def get_neighbors(self, of_node_with_ids:set, only_with_ids=[]):
        """
        Returns
        -------
        list of neighbor IDs
        """
        neighbor_ids = set()
        for i in [n for n in of_node_with_ids if self.has_node(n)]:
            neighbor_ids.update(set([x for x in nx.all_neighbors(self, i)]))
                
        if only_with_ids:
            return [i for i in list(neighbor_ids) if i in only_with_ids]
        else:
            return list(neighbor_ids)


    def get_neighbors_graph(self, for_node_ids:set, only_neighbors_with_ids=None):
        neighbors_ids = self.get_neighbors(for_node_ids,only_neighbors_with_ids)
        return self.get_subgraph(for_node_ids,neighbors_ids)


    def get_neighbors_rels(self, for_node_ids:set, only_neighbors_with_ids=None):
        neighbor_graph = self.get_neighbors_graph(for_node_ids,only_neighbors_with_ids)
        return neighbor_graph._relations()

    def get_neighbors_refs(self, for_node_ids:set, only_neighbors_with_ids=None):
        neighbor_graph = self.get_neighbors_graph(for_node_ids,only_neighbors_with_ids)
        return neighbor_graph.load_references()

    @staticmethod
    def get_att_set(prop_name:str,ps_objects:list):
        return set([i for sublist in [x[prop_name] for x in ps_objects] for i in sublist])


    def get_children_ids(self, for_node_ids:list, at_depth=1):
        """
        Input
        -----
        nodes in self must have ['Child Ids'] annotation
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


    def regulators(self, only_objtype=[], min_targets=1):
        if only_objtype:
            return {x for x, y in self.nodes(data=True) if ((self.out_degree(x) >= min_targets) & (y['ObjTypeName'][0] in only_objtype))}
        else:
            return {x for x in self.nodes() if self.out_degree(x) >= min_targets}


    def targets(self, only_objtype=[], min_regulators=1):
        if only_objtype:
            return {x for x, y in self.nodes(data=True) if ((self.in_degree(x) >= min_regulators) & (y['ObjTypeName'][0] in only_objtype))}
        else:
            return {x for x in self.nodes() if self.in_degree(x) >= min_regulators}

    def unconnected_node_ids(self):
        return {x for x in self.nodes() if self.degree(x) == 0}

    
    def root_nodes(self):
        return [i for i in self.nodes if self.in_degree(i)==0]


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


    def ref2pandas (self, relPropNames:list, entity_prop_names=[], RefNumPrintLimit=0) -> 'df':
        temp_fname = '__temp__.tsv'
        self.print_references(temp_fname,relPropNames,entity_prop_names,RefNumPrintLimit=RefNumPrintLimit)
        to_return = df.read(temp_fname,header=0,index_col=False, dtype='unicode')
        os.remove(temp_fname)
        return to_return

    '''
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
    '''

    def _2rnef_s(self,ent_props:list,rel_props:list,add_rel_props=dict(),add_pathway_props=dict()):
        # add_rel_props,add_pathway_props structure - {PropName:[PropValues]}
        resnet = et.Element('resnet')
        if add_pathway_props:
            pathway_props = et.SubElement(resnet, 'properties')
            for prop_name,prop_val in add_pathway_props.items():
                    for val in prop_val:
                        et.SubElement(pathway_props, 'attr', {'name':str(prop_name), 'value':str(val)})

        xml_nodes = et.SubElement(resnet, 'nodes')
        for nodeId, n in self.nodes(data=True):
            try:
                local_id = n['URN'][0]
                xml_node = et.SubElement(xml_nodes, 'node', {'local_id': local_id, 'urn': n['URN'][0]})
                et.SubElement(xml_node, 'attr', {'name': 'NodeType', 'value': str(n['ObjTypeName'][0])})
                if ent_props:
                    for prop_name, prop_values in n.items():
                        if prop_name in ent_props:
                            for prop_value in prop_values:
                                et.SubElement(xml_node, 'attr', {'name': str(prop_name), 'value': str(prop_value)})
            except KeyError:
                continue

        xml_controls = et.SubElement(resnet, 'controls')
        graph_relations = self._relations()
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

        xml_str = str(et.tostring(resnet, encoding='utf-8').decode("utf-8"))
        return xml_str


    def _2rnef_file(self,fname:str,ent_prop2print:list,rel_prop2print:list,add_rel_props=dict(),add_pathway_props=dict()):
        rnef_str = self._2rnef_s(ent_prop2print,rel_prop2print, add_rel_props,add_pathway_props)
        rnef_str = '<batch>\n'+str(rnef_str).strip()+'</batch>\n'
        pretty_xml = str(minidom.parseString(rnef_str).toprettyxml(indent='  '))
        with open(fname,'w',encoding='utf-8') as f: f.write(pretty_xml)


    def __2rnef(self,ent_prop2print:List,rel_prop2print:list,add_rel_props=dict(),with_section_size=1000):
        resnet_sections = ResnetGraph()
        graph_rnef_str = str()
        if self.number_of_edges():
            for regulatorID, targetID, e in self.edges(data='relation'):
                resnet_sections.copy_rel(e,self)
                if resnet_sections.number_of_edges() == with_section_size:
                    rnef_str = resnet_sections._2rnef_s(ent_prop2print,rel_prop2print,add_rel_props)
                    rnef_str = str(minidom.parseString(rnef_str).toprettyxml(indent='  '))
                    rnef_str = rnef_str[rnef_str.find('\n')+1:]
                    graph_rnef_str += rnef_str
                # file_out.write(rnef_str)
                    resnet_sections.clear_resnetgraph()

            rnef_str = resnet_sections._2rnef_s(ent_prop2print,rel_prop2print,add_rel_props)
            rnef_str = str(minidom.parseString(rnef_str).toprettyxml(indent='  '))
            rnef_str = rnef_str[rnef_str.find('\n')+1:]
            return graph_rnef_str + rnef_str
        else:
            all_nodes = self._get_nodes()
            for sec in range(0, len(all_nodes), with_section_size):
                section_nodes = all_nodes[sec:sec+with_section_size]
                section_graph = ResnetGraph()
                section_graph.add_nodes(section_nodes)
                rnef_str = section_graph._2rnef_s(ent_prop2print,rel_prop2print,add_rel_props)
                rnef_str = str(minidom.parseString(rnef_str).toprettyxml(indent='  '))
                rnef_str = rnef_str[rnef_str.find('\n')+1:]
                graph_rnef_str += rnef_str
            return graph_rnef_str


    def _2rnef(self,file_out,ent_prop2print:list,rel_prop2print:list,add_rel_props=dict(),with_section_size=1000):
        my_rnef_str = self.__2rnef(ent_prop2print,rel_prop2print,add_rel_props,with_section_size)
        file_out.write(my_rnef_str)


    def add_row2(self,to_df:df,from_relation_types:list, between_node_id, and_node_id, from_properties:list,
                cell_sep=';'):
        """
        Input
        -----
        from_properties format: [ 'N1:PropName', 'N2:PropName,'R:PropName']
        N1 prefix is for 'between_node_id' properties, N2 prefix is for 'and_node_id' properties, R prefix is for relation properties\n
        All properties are added in the order of 'from_properties'\n
        Multiple values from property are joined by cell_sep=';'\n
        Appends row 'to_df' generated from triples (between_node_id, and_node_id, rel) according to specifications in from_properties
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


    def tree4(self,root_node_id:int):
        tree_rn = ResnetGraph()
        tree = nx.bfs_tree(self,root_node_id)
        for r,t in tree.edges():
            rels = self._relation4(r,t)
            [tree_rn.copy_rel(rel,self) for rel in rels]
        return tree_rn


    def largest_tree(self):
        largest_tree = nx.DiGraph()
        root_nodes = self.root_nodes()
        for root in root_nodes:
            tree = nx.bfs_tree(self,root)
            if len(tree) > len(largest_tree):
                largest_tree = tree

        largest_tree_rn = ResnetGraph()
        for r,t in largest_tree.edges():
            rels = self._relation4(r,t)
            [largest_tree_rn.copy_rel(rel,self) for rel in rels]
        
        return largest_tree_rn


    def get_simple_paths2(self,end_tree_ids:list, with_only_objtypes=[]):
        """
        Returns
        -------
        {end_tree_id:{len(path):[all_simple_paths]}}
        """
        def _is_valid(path_obj:list):
            for obj in path_obj:
                if obj.objtype() not in with_only_objtypes:
                    return False
            return True

        id2paths = dict()
        root_ids = self.root_nodes()
        for end_id in end_tree_ids:
            end_paths = dict() # {len(path):[paths]}
            for root_id in root_ids:   
                paths = [p for p in nx.all_simple_paths(self,root_id,end_id)]
                for path in paths:
                    path_obj = self.path2obj(path)
                    if _is_valid(path_obj):
                        try:
                            end_paths[len(path)].append(path_obj)
                        except KeyError:
                            end_paths[len(path)] = [path_obj]
    
            end_paths_sorted = dict(sorted(end_paths.items()))
            id2paths[end_id] = end_paths_sorted
 
        return id2paths
            

    def path2obj(self,path:list, merge_rel2parent=True):
        """
        Input
        -----
        path - list of node ids specifying path
        """
        path_objs = list()
        if merge_rel2parent:
            for i in range(1,len(path)):
                parent_id = path[i-1]
                child_id = path[i]

                parent = self._get_node(parent_id)
                #child = self._get_node(child_id)
                path_rels = self.__find_relations(parent_id,child_id)
                #print('')
                [parent.merge_with_rel(r) for r in path_rels]
                path_objs.append(parent)

            path_objs.append(self._get_node(path[-1]))
        else:
            path_objs.append(self._get_node(path[0]))
            for i in range(1,len(path)):
                parent_id = path[i-1]
                child_id = path[i]

                parent = self._get_node(parent_id)
                child = self._get_node(child_id)
                path_rels = self.__find_relations(parent_id,child_id)
                [child.merge_with_rel(r) for r in path_rels]
                path_objs.append(child)

        return path_objs


    def nested_dic(self,top_level_class_prop:str,rel_key_prop:str,node_key_prop:str,rel_props:list,node_props:list,end_tree_ids:list):
    
        """
        Returns
        -------
        {level1_prop:{node_props:prop_values,{regulator_props:prop_values},{target_props:prop_values}},\n
        where "regulator_props" values are added from relation if node is regulator and\n "target_props" values are added from relation if node is target
        
        """
        def deref_multi(data:dict, keys:list):
            return deref_multi(data[keys[0]], keys[1:]) if keys else data

        def add_element(to_dict:dict, at_path:list, add_dict:dict, using_key:str):
            path_element = deref_multi(to_dict,at_path)
            try:
                exist_dict = path_element[using_key]
            except KeyError:
                path_element[using_key] = add_dict
                
            return at_path + [using_key]

        propval2objs, objid2propval = self.get_prop2obj_dic(top_level_class_prop)
        dict2return = dict()
        for value,nodes in propval2objs.items():
            dict2return[value] = dict()
            for node in nodes:
                node_prop_dict = node.props2dict(node_props)
                node_key = node._prop2str(node_key_prop)
                dict2return[value][node_key] = node_prop_dict

                #node_tree = self.tree4(node.id())
                #tree_leafs = {i for i in node_tree if node_tree.out_degree(i)==0}
                for leaf in end_tree_ids:
                    #curent_path = [value,node_key]
                    #shortest_path = nx.shortest_path(node_tree,node.id(),leaf)
                    #shortest_paths = nx.single_source_shortest_path(node_tree,node.id(),leaf)
                    shortest_paths = [p for p in nx.all_simple_paths(self,node.id(),leaf)]

                    for shortest_path in shortest_paths:
                        curent_path = [value,node_key]
                        for i in range(1,len(shortest_path)):
                            parent_id = shortest_path[i-1]
                            child_id = shortest_path[i]
                            
                            path_rels = self.__find_relations(parent_id,child_id)         
                            for rel in path_rels:
                                #curent_path = [value,node_key]
                                path_rel_key = rel_key_prop+':'+rel._prop2str(rel_key_prop)
                                if not path_rel_key: break

                                path_rel_props = rel.props2dict(rel_props)
                                curent_path = add_element(dict2return,curent_path,path_rel_props,path_rel_key)

                                child_node = self._get_node(child_id)
                                child_prop_dict = child_node.props2dict(node_props)
                                child_key = child_node._prop2str(node_key_prop)
                                try:
                                    child_top_level_class = child_node[top_level_class_prop][0]
                                    curent_path = add_element(dict2return,curent_path,{},child_top_level_class)
                                except KeyError:
                                    pass

                                if not child_key: break
                                curent_path = add_element(dict2return,curent_path,child_prop_dict,child_key)
                        
        return dict2return

########################  SUBGRAPH SUBGRAPH SUBGRAPH #####################################
    def subtract(self, other: "ResnetGraph"):
        #only self graph is analyzed 
        # works faster if other graph is bigger than self
        unique2self = ResnetGraph()
        edges_from_other = other._relations_ids()
        for n1,n2,e in self.edges(data=True):
            if e['relation']['Id'][0] not in edges_from_other:
                unique2self.add_edge(n1, n2, relation=e['relation'],weight=e['weight'], key=self.rel_urn(e['relation']))
                unique2self.add_nodes([self.nodes[n1],self.nodes[n2]])
        return unique2self


    def intersect (self, other: "ResnetGraph"):
        #only self graph is analyzed 
        # works faster if other graph is bigger than self
        intersection = ResnetGraph()
        edges_from_other = other._relations_ids()
        for n1,n2,e in self.edges(data=True):
            if e['relation']['Id'][0] in edges_from_other:
                intersection.add_edge(n1, n2, relation=e['relation'],weight=e['weight'],key=self.rel_urn(e['relation']))
                intersection.add_nodes([self.nodes[n1],self.nodes[n2]])
        return intersection


    def filter_references(self, keep_prop2values:dict, rel_types=[]):
        """
        Input
        -----
        keep_prop2values = {prop_name:[values]}
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
        Returns 
        -------
        subgraph with targets from 'for_targets_with_urns' and their regulators 
        function is used for SNEA
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
        Returns 
        -------
        graph with only one relation between pair of nodes.
        keeps relation with biggest reference count.
        all other relations are merged into the most referenced one
        """
        def find_best_rel(rels:list):
            rels.sort(key=lambda x: x[REFCOUNT][0], reverse=True)
            for index,rel in enumerate(rels):
                try:
                    eff = rel['Effect'][0]
                    if eff in ['positive', 'negative']:
                        return index
                    else:
                        continue
                except KeyError:
                    continue
            
            # to merge Binding to DirectRegulation
            for index,rel in enumerate(rels):
                if rel.is_directional():
                    return index

            return 0
        
        simple_g = ResnetGraph()
        for regulator_id, target_id, edges in self.edges(data=True):
            if not simple_g.has_edge(regulator_id, target_id): 
                reg2target_edges = dict(self[regulator_id][target_id])
                reg2target_rels = [e['relation'] for i,e in reg2target_edges.items()]

                best_rel_index = find_best_rel(reg2target_rels)
                best_rel = reg2target_rels[best_rel_index].copy()
                reg2target_rels.pop(best_rel_index)
                [best_rel.merge_rel(rel) for rel in reg2target_rels]
                    
                simple_g.copy_rel(best_rel,self)

        return simple_g
    

    def node_targets(self,regulator_id:int):
        target_ids = list(self.neighbors(regulator_id))
        return self._get_nodes(target_ids)

                    
    def regulome_dict(self, min_size=2):
        """
        Returns 
        -------
        {regulator_id: [targets]}, where len([targets]) >= min_size
        """
        regulator_ids = list(self.regulators(min_targets=min_size))
        subnetworks = dict()
        len_targets = dict()
        for regulator_id in regulator_ids:
            target_ids = list(self.neighbors(regulator_id))
            targets = self._get_nodes(target_ids)
            subnetworks[regulator_id] = targets
            len_targets[regulator_id] = len(target_ids)
        
        nx.set_node_attributes(self, len_targets, '# targets')

        print('Generated %d regulome subnetworks with more than %d targets from %s' 
                        % (len(subnetworks), min_size, self.name))
        return subnetworks


    def downstream_relations(self,node_id:int,with_types=list()):
        my_rels = [rel for r,t,rel in self.edges.data('relation') if r == node_id]
        if with_types:
            my_rels = [r for r in my_rels if r['ObjTypeName'] in with_types]
        return my_rels


    def upstream_relations(self,node_id:int,with_types=list()):
        my_rels = [rel for r,t,rel in self.edges.data('relation') if t == node_id]
        if with_types:
            my_rels = [r for r in my_rels if r['ObjTypeName'] in with_types]
        return my_rels


    def downstream_targets(self,of_node_id:int,linkedby_reltypes=list()):
        target_ids = [t for r,t,rel in self.edges.data('relation') if r == of_node_id and rel['ObjTypeName'][0] in linkedby_reltypes]
        return target_ids


    def get_regulome(self, start_node_ids: set):
        """
        Returns 
        -------
        composition of bfs_trees for all ids in start_node_ids
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


    def subgraph_by_rels(self, rels:list)->'ResnetGraph':
        """
        Input
        -----
        rels = [PSRelation]
        """
        subgraph = ResnetGraph()
        [subgraph.copy_rel(rel,self) for rel in rels]
        return subgraph


    def subgraph_by_relprops(self, search_values:list, in_properties:list=['ObjTypeName']):
        search_value_set = set(search_values)
        subgraph = ResnetGraph()
        for prop_type in in_properties:
            for n1, n2, rel in self.edges.data('relation'):
                if not set(rel[prop_type]).isdisjoint(search_value_set):
                    subgraph.copy_rel(rel,self)
                  #  subgraph.add_edge(n1, n2, relation=rel,weight=rel.get_reference_count(),key=self.rel_urn(rel))
                   # subgraph.add_nodes([self.nodes[n1],self.nodes[n2]])
        return subgraph


    def get_subgraph(self,between_node_ids:list,and_node_ids:list,by_relation_type=None,with_effect=None,in_direction=None)->"ResnetGraph":
        subgraph = ResnetGraph()
        for n1 in between_node_ids:
            for n2 in and_node_ids:
                if not isinstance(in_direction,str) or in_direction == '>':
                    if self.has_edge(n1, n2):
                        for r in self[n1][n2].values(): #gives error if edge does not exist
                            rel = r['relation']
                            if isinstance(by_relation_type,list) and rel['ObjTypeName'] not in by_relation_type: continue
                            if isinstance(with_effect,list):
                                try: ef = rel['Effect'] 
                                except KeyError:
                                    rel['Effect'] = 'unknown'
                                    ef = rel['Effect']
                                if ef not in with_effect: continue
                            subgraph.copy_rel(rel,self)
                         #   subgraph.add_edge(n1, n2, relation=rel,weight=r['weight'],key=self.rel_urn(rel))
                         #   subgraph.add_nodes([self.nodes[n1],self.nodes[n2]])
                if not isinstance(in_direction,str) or in_direction == '<':
                    if self.has_edge(n2, n1):
                        for r in self[n2][n1].values():#gives error if edge does not exist
                            rel = r['relation']
                            if isinstance(by_relation_type,list) and rel['ObjTypeName'] not in by_relation_type: continue
                            if isinstance(with_effect,list):
                                try: ef = rel['Effect'] 
                                except KeyError:
                                    rel['Effect'] = 'unknown'
                                    ef = rel['Effect']
                                if ef not in with_effect: continue
                            subgraph.copy_rel(rel,self)
                           # subgraph.add_edge(n2, n1, relation=rel,weight=r['weight'],key=self.rel_urn(rel))
                           # subgraph.add_nodes([self.nodes[n1],self.nodes[n2]])
            
        return subgraph


    def centrality(self):
        '''
        Returns
        -------
        nodes : dictionary
        Dictionary of nodes with degree centrality as the value.
        '''
        return nx.degree_centrality(self)


class PSPathway(PSObject):
    pass
    def __init__(self,dic=dict(),g=ResnetGraph()):
        super().__init__(dic)
        self.graph = ResnetGraph()
        
        self.update(dic)
        self.graph.update(g)
        self.graph.urn2obj.update(g.urn2obj)
        self.graph.urn2rel.update(g.urn2rel)


    @classmethod
    def from_resnet(cls, resnet:et.Element, add_annotation=dict()):
        pathway = PSPathway()
        pathway['Name'] = [resnet.get('name')]
        pathway['URN'] = [resnet.get('urn')]
        pathway.update(add_annotation)
        [pathway.update_with_value(p.get('name'),p.get('value')) for p in resnet.findall('properties/attr')]

        if pathway.graph._parse_nodes_controls(resnet):
            pathway[RESNET] = et.tostring(resnet, encoding='unicode', method='xml')
            return pathway
        else:
            print('Invalid <resnet> section')
            return None

    @classmethod
    def from_pathway(cls,other:"PSPathway"):
        return PSPathway(other,other.graph)

    def number_of_nodes(self, obj_type=''):
        if obj_type:
            try:
                return self['#'+obj_type]
            except KeyError:
                obj_count = len([o for o in self.graph.urn2obj.values() if obj_type in o['ObjTypeName']])
                self['#'+obj_type] = obj_count
                return obj_count
        else:
            return len(self.graph.urn2obj)


    def merge_pathway(self, other:"PSPathway"):
        self.merge_obj(other)
        self.graph.add_graph(other.graph)
        self.graph.urn2obj.update(other.graph.urn2obj)
        self.graph.urn2rel.update(other.graph.urn2rel)


import networkx as nx
from ..pandas.panda_tricks import df
from datetime import timedelta
import os, math, time
from .rnef2sbgn import minidom
import xml.etree.ElementTree as et

from .NetworkxObjects import PSObject,PSRelation,len,REGULATORS,TARGETS,CHILDS,REFCOUNT
from ..ETM_API.references import PUBYEAR,EFFECT,TITLE,pubmed_hyperlink, make_hyperlink
from ..ETM_API.etm import ETMstat
from itertools import combinations

RESNET = 'resnet'
PHYSICAL_INTERACTIONS = ['Binding','DirectRegulation','ProtModification','PromoterBinding','ChemicalReaction']
PROTEIN_TYPES = ['Protein','FunctionalClass','Complex']
PS_REF_COULUMN = 'Number of reference. Link opens recent publications in PubMed'
NUMBER_OF_TARGETS = '# targets'


def execution_time(execution_start):
    return "{}".format(str(timedelta(seconds=time.time() - execution_start)))

class ResnetGraph (nx.MultiDiGraph):
    pass
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.urn2obj = dict() #lookup for combining nodes in database graph and rnef graph
        self.urn2rel = dict() #lookup for combining relations in database graph and rnef graph
        self.rnef_node_count = 0

    def copy(self)->'ResnetGraph':
        cp = super().copy(self)
        cp.urn2obj = dict(self.urn2obj)
        cp.urn2rel = dict(self.urn2rel)
        cp.rnef_node_count = self.rnef_node_count
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
        non-directional relations are duplicated in all possible directions
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
        nodes connected by relation must exist in the graph.
        non-directional relations are duplicated in all possible directions
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
        '''
        copies rel and its nodes specified in rel.Nodes from_graph to self\n
        re-assigns new Id[0] and rel.Nodes if node with the same URN is present in self
        '''
        rel2add = rel.copy()
        nodes2add = list()
        for s, node_tup_list in rel.Nodes.items():
        # rel.Nodes = {"Regulators':[(entityID, 0, effect)], "Targets':[(entityID, 1, effect)]}
            for i, node_tup in enumerate(node_tup_list):
                from_graph_node_id = node_tup[0]
                my_node = from_graph._get_node(from_graph_node_id)
                try:
                    exist_node = self.urn2node(my_node.urn())
                    exist_node_id = exist_node.id()
                    rel2add.Nodes[s][i] = (exist_node_id,node_tup[1],node_tup[2])
                except KeyError:
                    if rel.is_from_rnef():
                        my_node['Id'][0] = self.rnef_node_count                        
                        rel2add.Nodes[s][i] = (self.rnef_node_count,node_tup[1],node_tup[2])
                        nodes2add.append((self.rnef_node_count,my_node.items()))
                        self.rnef_node_count += 1
                    else:
                        nodes2add.append((from_graph_node_id,my_node.items()))
                    
                    self.urn2obj[my_node.urn()] = my_node
                    
        self.add_nodes_from(nodes2add)
        self.add_rel(rel2add)
        return

    
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
        self.add_nodes(dict(other.nodes(data=True)))
        other_rels = other._relations()
        [self.copy_rel(rel,other) for rel in other_rels]


    def compose(self,other:"ResnetGraph"):
        '''
        Returns
        -------
        new graph created by nx.compose\n
        other graph attributes take precedent
        '''
        composed_graph = ResnetGraph(nx.compose(self,other))
        composed_graph.urn2obj.update(other.urn2obj)
        composed_graph.urn2rel.update(other.urn2rel)
        return composed_graph

################## SET SET SET ##########################################
    def set_node_annotation(self, urn2value:dict, new_prop_name:str):
        """
        Input
        -----
        urn2value = {urn:[with_prop_values]}

        Adds
        ----
        new property to exsisting nodes. Existing values of 'new_prop_name' will be replaced
        """
        id2value = dict() #remapping urn to ids
        for urn, value in urn2value.items():
            try:
                node2annotate = self.urn2node(urn)
                id2value[node2annotate['Id'][0]] = value
            except KeyError: continue

        nx.set_node_attributes(self,id2value,new_prop_name)
        #print('%d nodes were annotated with attributes "%s"' % (len(id2value),new_prop_name))


    def set_edge_property(self,between_node_id,and_node_id,with_rel_urn:str,prop_name,prop_values:list):
        annotate_values = set(prop_values)
        if self.has_edge(between_node_id, and_node_id,with_rel_urn):
            try:
                my_prop_values = self[between_node_id][and_node_id][with_rel_urn]['relation'][prop_name]
                annotate_values.update(my_prop_values)
            except KeyError:
                pass
            self[between_node_id][and_node_id][with_rel_urn]['relation'][prop_name] = list(annotate_values)


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


    def centrality(self):
        '''
        Returns
        -------
        nodes : dictionary
        Dictionary of nodes with degree centrality as the value.
        '''
        return nx.degree_centrality(self)


    def rank_regulator(self, regulator_id:int, target_weights:dict, max_distance=5):
        '''
        Input
        -----
        target_weights = {node_id:weight}

        Returns
        -------
        regulator rank for regulator_id
        '''
        regulator_rank = float(0.0)
        regulator_tree = nx.bfs_tree(self, regulator_id)

        for level in range (1,max_distance+1):
            targets_on_level = nx.descendants_at_distance(regulator_tree, regulator_id, level)
            if not targets_on_level: break
            for target_id in targets_on_level:
                try:
                    target_weight = float(target_weights[target_id])
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
                refcount = pubmed_hyperlink(pmids,total_refs)
            except KeyError:
                try:
                    dois = ref_ids['DOI']
                    refcount = make_hyperlink(dois[0],url='http://dx.doi.org/', display_str=str(total_refs))
                except KeyError:
                    refcount = str(total_refs)
            annoate_df.loc[row][PS_REF_COULUMN] = refcount
        return annoate_df

    
    def recent_refs2df(self, ref_limit=5,pathway=''):
        self.load_references()
        return_df = df(name='RefCount',columns=['Regulator','Target',PS_REF_COULUMN,'Number of References','Pathway'])
        row_counter = 0
        for regulator_id, target_id,rel in self.edges(data='relation'):
            regulator = self._get_node(regulator_id)
            target = self._get_node(target_id)
            regulator_name = regulator.name()
            target_name = target.name()
            total_refs, ref_ids, ps_references = self.__recent_refs([rel],ref_limit)
            dois = str()
            try:
                pmids = ref_ids['PMID']
                refcount = pubmed_hyperlink(pmids,total_refs)
            except KeyError:
                try:
                    dois = ref_ids['DOI']
                    refcount = make_hyperlink(dois[0],url='http://dx.doi.org/', display_str=str(total_refs))
                except KeyError:
                    refcount = str(total_refs)
            return_df.loc[row_counter] = [regulator_name,target_name,refcount,total_refs,pathway]
            row_counter += 1
        return_df.sort_values(by='Number of References',ascending=False, inplace=True)
        snippets_df = self.snippets2df()
        return return_df,snippets_df


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
        snippet_df.set_hyperlink_color(['PMID','DOI'])
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
        try:
            return PSObject({k:v for k,v in self.nodes[nodeId].items()})
        except KeyError:
            raise KeyError


    def urn2node(self, urn:str):
        try:
            return PSObject(self.urn2obj[urn])
        except KeyError:
            raise KeyError


    def _get_nodes(self, node_ids:list=[]):
        if not node_ids:
            return [PSObject(ddict) for i,ddict in self.nodes(data=True)]
        else:
            id_set = set(node_ids)
            return [PSObject(ddict) for i,ddict in self.nodes(data=True) if i in id_set]    
            '''
            #may be faster?
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
        return [o for i,o in self.nodes(data=True) if prop_value in o[prop_name]]
        return self._get_nodes(obj_ids)


    def get_objects(self, with_properties:list=['ObjTypeName'], only_with_values=[]):
        return_node_ids = set()
        if only_with_values:
            my_values = set(only_with_values)
            for prop_name in with_properties:
                id2prop = nx.get_node_attributes(self,prop_name)
                my_ids = [i for i,v in id2prop.items() if not my_values.isdisjoint(v)]
                return_node_ids.update(my_ids)
        else:
            for prop_name in with_properties:
                id2prop = nx.get_node_attributes(self,prop_name)
                return_node_ids.update(id2prop.keys())
            
        return self._get_nodes(list(return_node_ids))


    def _relations(self):
        """
        Returns
        -------
        [PSRelation] - list of all relations in graph
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
        [positive_refs.update(rel._get_refs(sort_by='')) for rel in positive_rel]
        negative_refs = set()
        [negative_refs.update(rel._get_refs(sort_by='')) for rel in negative_rel]

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


    def __rels4pair(self, entity1prop:str, entity2prop:str, search_property='Name', filter_rel_with:dict={}, with_children=False):
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
        entities1 = self.get_objects(with_properties=[search_property],only_with_values=[entity1prop])
        entities2 = self.get_objects(with_properties=[search_property],only_with_values=[entity2prop])
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

    @staticmethod
    def __recent_refs(relations:list,ref_limit=5):
        '''
        Input
        -----
        [PSRelations] - list of PSRelation objects
        '''
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


    def recent_refs(self, entity1prop:str, entity2prop:str, search_property='Name', ref_limit=5,with_children=False):
        """
        Returns
        -------
        total_refs, ref_ids={id_type:[identifiers]}, recent_refs={Reference}
        """
        relations = self.__rels4pair(entity1prop, entity2prop, search_property,with_children=with_children)
        return self.__recent_refs(relations,ref_limit)
    

    def _relations_ids(self):
        return {rel['Id'][0] for regulatorID, targetID, rel in self.edges.data('relation')}


    def find_nodes(self,for_relation:PSRelation,filter_by:list=None,in_properties:list=None):
        """
        Returns
        -------
        regulators = [PSObject], targets = [PSObject]
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


    def get_neighbors_graph(self, for_node_ids:set, only_neighbors_with_ids:list=[],by_relation_types:list=[],with_effects:list=[],in_direction=None):
        neighbors_ids = self.get_neighbors(for_node_ids,only_neighbors_with_ids)
        return self.get_subgraph(for_node_ids,neighbors_ids,by_relation_types,with_effects,in_direction)


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

    
    def prop2outdegree(self,prop_name='Name'):
        return {p[0]:self.out_degree(i) for i,p in self.nodes(data=prop_name)}


################################# WRITE-DUMP, WRITE-DUMP, WRITE-DUMP ##############################
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


    def ref2pandas (self, relPropNames:list, entity_prop_names=[], RefNumPrintLimit=0,single_rel_row=False) -> 'df':
        temp_fname = '__temp__.tsv'
        self.print_references(temp_fname,relPropNames,entity_prop_names,RefNumPrintLimit=RefNumPrintLimit,single_rel_row=single_rel_row)
        to_return = df.read(temp_fname,header=0,index_col=False, dtype='unicode')
        os.remove(temp_fname)
        return to_return


    def rnef(self,ent_props:list,rel_props:list,add_rel_props:dict={},add_pathway_props:dict={}):
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
                    try:
                        regulator_local_id = self.nodes[r[0]]['URN'][0]
                        et.SubElement(xml_control, 'link', {'type':'in', 'ref':regulator_local_id})
                    except KeyError:
                        continue

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


    def rnef2sections(self,ent_prop2print:list,rel_prop2print:list,add_rel_props=dict(),with_section_size=1000):
        """
        resolves printing graphs with and without edges\n
        splits RNEF into <resnet> sections with_section_size <control> elements
        
        Returns
        -------
        sectioned RNEF string
        """
        resnet_sections = ResnetGraph()
        graph_rnef_str = str()
        if self.number_of_edges():
            for regulatorID, targetID, e in self.edges(data='relation'):
                resnet_sections.copy_rel(e,self)
                if resnet_sections.number_of_edges() == with_section_size:
                    rnef_str = resnet_sections.rnef(ent_prop2print,rel_prop2print,add_rel_props)
                    rnef_str = str(minidom.parseString(rnef_str).toprettyxml(indent='  '))
                    rnef_str = rnef_str[rnef_str.find('\n')+1:]
                    graph_rnef_str += rnef_str
                    resnet_sections.clear_resnetgraph()

            rnef_str = resnet_sections.rnef(ent_prop2print,rel_prop2print,add_rel_props)
            rnef_str = str(minidom.parseString(rnef_str).toprettyxml(indent='  '))
            rnef_str = rnef_str[rnef_str.find('\n')+1:]
            return graph_rnef_str + rnef_str
        else:
            all_nodes = self._get_nodes()
            for sec in range(0, len(all_nodes), with_section_size):
                section_nodes = all_nodes[sec:sec+with_section_size]
                section_graph = ResnetGraph()
                section_graph.add_nodes(section_nodes)
                rnef_str = section_graph.rnef(ent_prop2print,rel_prop2print,add_rel_props)
                rnef_str = str(minidom.parseString(rnef_str).toprettyxml(indent='  '))
                rnef_str = rnef_str[rnef_str.find('\n')+1:]
                graph_rnef_str += rnef_str
            return graph_rnef_str


    def rnefsections2stream(self,file_out,ent_prop2print:list,rel_prop2print:list,add_rel_props=dict(),with_section_size=1000):
        my_rnef_str = self.rnef2sections(ent_prop2print,rel_prop2print,add_rel_props,with_section_size)
        file_out.write(my_rnef_str)


    def rnef2file(self,fname:str,ent_prop2print:list,rel_prop2print:list,add_rel_props:dict={},add_pathway_props:dict={}):
        rnef_str = self.rnef(ent_prop2print,rel_prop2print, add_rel_props,add_pathway_props)
        rnef_str = '<batch>\n'+str(rnef_str).strip()+'</batch>\n'
        pretty_xml = str(minidom.parseString(rnef_str).toprettyxml(indent='  '))
        with open(fname,'w',encoding='utf-8') as f: 
            f.write(pretty_xml)


    def add_row2(self,to_df:df,from_relation_types:list,between_node_id,and_node_id,from_properties:list,cell_sep=';'):
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
                cell = node4cell.prop_values2str(from_properties[i][3:],sep=cell_sep)
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
    def _parse_nodes_controls(self, resnet:et.Element):
        nodel_local_ids = dict()
        for node in resnet.findall('./nodes/node'):
            node_urn = node.get('urn')
            local_id = node.get('local_id')
            try:
                node_obj = self.urn2node(node_urn)
            except KeyError:
                node_obj = PSObject({'Id':[self.rnef_node_count],'URN':[node_urn]})
                self.rnef_node_count += 1

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
            self.add_rel(ps_rel) #ps_rel does not have 'Id' property.This is used by PSRelation.is_from_rnef()
        return self.rnef_node_count


    def __read_rnef(self, rnef_file:str, last_new_node_id=1):
        try:
            print ('\nLoading graph from file %s' % rnef_file,flush=True)
            root = et.parse(rnef_file).getroot()  
        except FileNotFoundError:
            raise FileNotFoundError

        for resnet in root.findall('resnet'):
            self._parse_nodes_controls(resnet)


    @classmethod
    def fromRNEF(cls,rnef_file:str):
        try:
            g = ResnetGraph()
            start = time.time()
            g.__read_rnef(rnef_file)
            print('File %s with %d edges and %d nodes was loaded in %s' 
        % (rnef_file,g.number_of_edges(),g.number_of_nodes(),execution_time(start)))
            return g
        except FileNotFoundError:
            raise FileNotFoundError


    def tree4(self,root_node_id:int,reverse=False):
        tree_rn = ResnetGraph()
        tree = nx.bfs_tree(self,root_node_id,reverse)
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
        def _is_valid(path_objs:list):
            for obj in path_objs:
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
                    path_objs = self.path2obj(path)
                    if _is_valid(path_objs):
                        try:
                            end_paths[len(path)].append(path_objs)
                        except KeyError:
                            end_paths[len(path)] = [path_objs]
    
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
        {top_level_class_prop:{node_props:prop_values,{regulator_props:prop_values},{target_props:prop_values}},\n
        where "regulator_props" values are added from relation if node is regulator and\n 
        "target_props" values are added from relation if node is a target
        top_level_class_prop defines first  in every path
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
        subgraph with targets from 'for_targets_with_urns' and their regulators\n 
        function is used for SNEA
        """
        start = time.time()
        targets_urns = set(for_targets_with_urns)
        input_target_ids = [n['Id'][0] for u,n in self.urn2obj.items() if u in targets_urns]

        reg_subgraph = ResnetGraph()
        rel2copy = [e for r,t,e in self.edges(data='relation') if t in set(input_target_ids)]
        [reg_subgraph.copy_rel(rel,self) for rel in rel2copy]
            
        reg_subgraph.name = network_name
        print('%s subnetwork with %d edges was selected in %s from network with %d edges' % (network_name, reg_subgraph.number_of_edges(), self.execution_time(start),self.number_of_edges()))
        return reg_subgraph


    def make_simple(self, rel_type_rank:list=[]):
        """
        Returns 
        -------
        graph with only one edge between nodes.\n
        keeps relation with the biggest reference count.\n
        all other relations are merged into the most referenced one
        """
        def find_best_effect(rels:list):
            rels.sort(key=lambda x: x[REFCOUNT][0], reverse=True)
            for index,rel in enumerate(rels):
                try:
                    eff = rel.effect()
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

       # REL_TYPE_RANK = ['DirectRegulation','Binding','ProtModification','PromoterBinding',
        #'ChemicalReaction','Expression','Biomarker','QuantitativeChange','StateChange','MolSynthesis','MolTransport','CellExpression','Regulation','FunctionalAssosiation']
        
        def find_best_type(rels:list):
            for rel_type in rel_type_rank:
                for rel in rels:
                    if rel.objtype() == rel_type:
                        if rel_type == 'Binding':
                            return 'DirectRegulation'
                        else:
                            return rel_type
            return ''
        
        simple_g = ResnetGraph()
        for regulator_id, target_id in self.edges():
            if not simple_g.has_edge(regulator_id, target_id): 
                reg2target_edges = dict(self[regulator_id][target_id])
                reg2target_rels = [e['relation'] for e in reg2target_edges.values()]

                if len(reg2target_rels) > 1:
                    best_rel_index = find_best_effect(reg2target_rels)
                    best_rel = reg2target_rels[best_rel_index].copy()
                    
                    best_rel_type = find_best_type(reg2target_rels)
                    if best_rel_type:
                        best_rel['ObjTypeName'][0] = best_rel_type

                    reg2target_rels.pop(best_rel_index)
                    [best_rel.merge_rel(rel) for rel in reg2target_rels]
                else:
                    best_rel = reg2target_rels[0]

                simple_g.copy_rel(best_rel,self)

        print('%d edges in graph were simplified' % (self.number_of_edges()-simple_g.number_of_edges()))
        simple_g.rnef_node_count = self.rnef_node_count

        return simple_g
    

    def node_targets(self,regulator_id:int):
        target_ids = list(self.neighbors(regulator_id))
        return self._get_nodes(target_ids)

                    
    def regulome_dict(self, only_objtype:list, min_size=2):
        """
        Returns 
        -------
        {regulator_id : [targets]}, where len([targets]) >= min_size
        """
        regulator_ids = list(self.regulators(only_objtype,min_targets=min_size))
        subnetworks = dict()
        len_targets = dict()
        for regulator_id in regulator_ids:
            target_ids = list(self.neighbors(regulator_id))
            targets = self._get_nodes(target_ids)
            subnetworks[regulator_id] = targets
            len_targets[regulator_id] = len(target_ids)
        
        nx.set_node_attributes(self, len_targets, NUMBER_OF_TARGETS)

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


    def upstream_regulators(self,of_node_id:int,linkedby_reltypes=list()):
        regulator_ids = [r for r,t,rel in self.edges.data('relation') if t == of_node_id and rel['ObjTypeName'][0] in linkedby_reltypes]
        return regulator_ids


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
        if search_values:
            for prop_type in in_properties:
                for n1, n2, rel in self.edges.data('relation'):
                    try:
                        if not set(rel[prop_type]).isdisjoint(search_value_set):
                            subgraph.copy_rel(rel,self)
                    except KeyError:
                        continue
        else:
            for prop_name in in_properties:
                for n1, n2, rel in self.edges.data('relation'):
                    try:
                        my_rel = rel[prop_type]
                        subgraph.copy_rel(my_rel,self)
                    except KeyError:
                        continue
        return subgraph


    def subgraph_by_nodeprops(self,has_properties:list,with_values=set()):
        '''
        Returns
        -------
        neighborhood of nodes that has_properties if with_values is empty\n
        otherwise neighborhood of nodes that has_properties with_values
        '''
        my_nodes = self.get_objects(has_properties,with_values)
        my_nodes_ids = [n.id() for n in my_nodes]
        return self.get_neighbors_graph(my_nodes_ids)


    def get_subgraph(self,between_node_ids:list,and_node_ids:list,by_relation_types:list=[],with_effect:list=[],in_direction=None)->"ResnetGraph":
        subgraph = ResnetGraph()
        for n1 in between_node_ids:
            for n2 in and_node_ids:
                if not isinstance(in_direction,str) or in_direction == '>':
                    if self.has_edge(n1, n2):
                        for r in self[n1][n2].values(): #gives error if edge does not exist
                            rel = r['relation']
                            if by_relation_types and rel.objtype() not in by_relation_types: continue
                            if with_effect:
                                try: ef = rel.effect()
                                except KeyError:
                                    rel['Effect'] = ['unknown']
                                    ef = rel.effect()
                                if ef not in with_effect: continue
                            subgraph.copy_rel(rel,self)
                         #   subgraph.add_edge(n1, n2, relation=rel,weight=r['weight'],key=self.rel_urn(rel))
                         #   subgraph.add_nodes([self.nodes[n1],self.nodes[n2]])
                if not isinstance(in_direction,str) or in_direction == '<':
                    if self.has_edge(n2, n1):
                        for r in self[n2][n1].values():#gives error if edge does not exist
                            rel = r['relation']
                            if by_relation_types and rel.objtype() not in by_relation_types: continue
                            if with_effect:
                                try: ef = rel.effect() 
                                except KeyError:
                                    rel['Effect'] = ['unknown']
                                    ef = rel.effect()
                                if ef not in with_effect: continue
                            subgraph.copy_rel(rel,self)
                           # subgraph.add_edge(n2, n1, relation=rel,weight=r['weight'],key=self.rel_urn(rel))
                           # subgraph.add_nodes([self.nodes[n1],self.nodes[n2]])
            
        return subgraph


    def ontology_resnet(self,members:list,add2parent:PSObject):
        """
        Input
        -----
        members - [PSObjects]
        """
        resnet = ResnetGraph()
        parent_id = add2parent.id()
        resnet.add1node(add2parent)
        resnet.add_nodes(members)
        for m in members:
            child_id = m['Id'][0]
            rel = PSRelation({'ObjTypeName':['MemberOf'],'Relationship':['is-a'],'Ontology':['Pathway Studio Ontology']})
            rel.Nodes['Regulators'] = [(child_id,0,0)]
            rel.Nodes['Targets'] = [(parent_id,1,0)]
            rel.append_property('Id',int(str(parent_id)+str(child_id)))
            rel.append_property('URN', self.rel_urn(rel)) #str(child_id)+'0') #fake URN
            resnet.add_edge(child_id,parent_id,relation=rel)

        return resnet


    def ontology_graph(self):
        '''
        Input
        -----
        nodes in self are annotated with [CHILDS] property
        '''
        ontology_graph = ResnetGraph()
        for parent in self._get_nodes():
            try:
                child_ids = parent[CHILDS]
                children = self._get_nodes(child_ids)
                resnet = self.ontology_resnet(children,parent)
                ontology_graph.add_graph(resnet)
            except KeyError:
                continue

        return ontology_graph


    def all_paths_from(self,child_id:int):
        '''
        self must be ontology_graph
        returns [[PSObject]] - list of all pathes to node_id
        '''
        parent_tree = self.tree4(child_id)
        top_parent_ids = {x for x in parent_tree.nodes() if parent_tree.out_degree(x) == 0}
        all_ontology_paths = list()
        for parent_id in top_parent_ids:
            ontology_id_paths = list(nx.all_simple_paths(parent_tree,child_id,parent_id))
            for path in ontology_id_paths:
                ontology_obj_path = self.path2obj(list(path))
                all_ontology_paths.append(ontology_obj_path)

        return all_ontology_paths


    def direct_indirect_targets(self,of_node_id,
    direct_reltypes:list=['DirectRegulation','Binding','ChemicalReaction','ProtModification','PromoterBinding'],
    indirect_reltypes=['Regulation','MolTransport','Expression','MolSynthesis']):
        direct_neighbors_ids = set()
        indirect_neighbors_ids = set()
        #neighbor_ids = list(self[of_node_id])
        for neighbor_id in self[of_node_id]:
            for urn in self[of_node_id][neighbor_id]:
                #rel = self[of_node_id][neighbor_id][urn]['relation']
                objtype  = self[of_node_id][neighbor_id][urn]['relation']['ObjTypeName'][0]
                if objtype in direct_reltypes:
                    direct_neighbors_ids.add(neighbor_id)
                elif objtype in indirect_reltypes:
                    indirect_neighbors_ids.add(neighbor_id)

        return direct_neighbors_ids, indirect_neighbors_ids


    def split2puddles(self):
        '''
        Input: Connected Graph
        ----------------------

        Returns
        -------
        communities_subgraphs - [nx.Graph]:\n
        list of subgraphs with non-overlapping set of nodes generated by \n
        due to small-world structure of ResnetGraph the returned list usually has several small subgraphs and 1 big subgraph\n
        big subgraph can be 80-90% of the input ResnetGraph
        '''
        my_graph = nx.Graph(self)
        small_components_rels = set()
        communities = list()
        components = nx.connected_components(my_graph)

        large_comminity_count = 0
        for i, component in enumerate(components):
            component_graph = my_graph.subgraph(component)
            large_comminity_count += 1
            puddles = list(nx.community.asyn_fluidc(component_graph,100))
            print('Created %d communities from component with %d nodes and %d edges' % 
            (len(puddles),component_graph.number_of_nodes(),component_graph.number_of_edges()))
            communities += puddles
        else:
            component_rels = [e for r,t,e in component_graph.edges(data='relation')]
            small_components_rels.update(component_rels)

        print('Found %d components' % i)

        communities.sort(key=len,reverse=True)
        communities_subgraphs = list()
        for community in communities:
            community_subgraph = self.subgraph(community)
            communities_subgraphs.append(community_subgraph)

        small_components_graph = self.subgraph_by_rels(small_components_rels)
        communities_subgraphs.append(small_components_graph)

        return communities_subgraphs


    @staticmethod
    def trisect(graph:'ResnetGraph'):
        '''
        Input: Connected Graph
        ----------------------

        Returns
        -------
        tuple containing nx.Graph objects: kl_partion1_graph, kl_partion2_graph, overlap_graph, where\n
        kl_partion1_graph, kl_partion2_graph - subgraphs with non-overlaping set of nodes made by kernighan_lin_bisection algorithm\n
        overlap_graph - nx.Graph containing edges common between two non-overlapping sectors from above
        '''
        nodes4partion1,nodes4partion2 = nx.community.kernighan_lin_bisection(nx.Graph(graph))
        kl_partion1 = graph.subgraph(nodes4partion1)
        kl_partion2 = graph.subgraph(nodes4partion2)

        partion1_rels = {e for r,t,e in kl_partion1.edges(data='relation')}
        partion2_rels = {e for r,t,e in kl_partion2.edges(data='relation')}

        kl_partion1_graph = graph.subgraph_by_rels(partion1_rels)
        kl_partion2_graph = graph.subgraph_by_rels(partion2_rels)

        common_rels = set(graph._relations()).difference(partion1_rels|partion2_rels)
        overlap_graph = graph.subgraph_by_rels(common_rels)

        return kl_partion1_graph,kl_partion2_graph, overlap_graph


    def components(self):
        '''
        Returns
        -------
        [nx.Graph] - List of unconnected subgraphs that do not have common edges
        '''
        components = nx.connected_components(nx.Graph(self))
        components_graphs = list()
        for c, component in enumerate(components):
            component_graph = self.subgraph(component)
            components_graphs.append(component_graph)

        return components_graphs


    def split2threads(self, max_edge_count=10000):
        '''
        Returns
        -------
        tuple of 3 lists - thread1, thread2, for_thread3:\n
        multithread1 - has non-overlapping subgraphs with number of edges less than max_edge_count
        multithread2 - has non-overlapping subgraphs with number of edges less than max_edge_count, which overlap with thread1
        if for_thread3 is not empty - subgraphs with number of edges more than max_edge_count. Threads in for_thread3 must be split further
        '''
        partition1, partition2, overlap = ResnetGraph.trisect(self)
        for_multithread2 = [overlap]
        multithread1 = [partition1, partition2]

        need_split = True
        while need_split:
            partitions2split = list(multithread1)
            for partition in partitions2split:
                if partition.number_of_edges() > max_edge_count:
                    p1,p2, o = ResnetGraph.trisect(partition)
                    multithread1 += [p1,p2]
                    for_multithread2.append(o)

            if len(multithread1) == len(partitions2split):
                # no split happened
                break
            else:
                continue

        for_multithread3 = list()
        multithread2 = list()
        for g in for_multithread2:
            if g.number_of_edges() > max_edge_count:
                for_multithread3.append(g)
            else:
                multithread2.append(g)

        return multithread1, multithread2, for_multithread3


    def split(self, max_edge_count=10000):
        '''
        Returns
        -------
        multithreads - [list], where each "list" contains non-overlapping subgraphs with number of edges less than max_edge_count\n
        subgraphs from different lists can overlap\n
        each list can be imported into graph database using multithreading without deadlocks
        '''
        t1,t2,leftover = ResnetGraph.split2threads(self,max_edge_count)
        multithreads = [t1,t2] # list of lists of graphs

        while leftover:
            for g in leftover:
                t2add1,t2add2,leftover = ResnetGraph.split2threads(g,max_edge_count)
                multithreads.append(t2add1)
                multithreads.append(t2add2)

        return multithreads


    def splitOLD(self, max_edge_count=10000):
        #nodes_per_subgraph = int(len(self)/number_of_subgraphs)
        isolated_subgraphs = list()
        rels_in_subgraphs = set()
        common_edges_leftovers = list()

        my_graph = nx.Graph(self)
        components = nx.connected_components(my_graph)
        small_components_graph = nx.Graph()
        for c, component in enumerate(components):
            component_graph = self.subgraph(component)
            if component_graph.number_of_edges() > max_edge_count:
               isolated_partitions, leftovers = self.split_multi()
               isolated_subgraphs += isolated_partitions
               common_edges_leftovers += leftovers
            else:
                small_components_graph.add_graph(component_graph)
                component_rels = [e for r,t,e in component_graph.edges(data='relation')]
                rels_in_subgraphs.update(component_rels)

        print('Graph has %d components' % c)
        isolated_subgraphs.append(small_components_graph)
        all_relations = set(self._relations())
        leftover_rels = all_relations.difference(rels_in_subgraphs)  #{r for r in all_relations if r not in rels_in_subgraphs}
        leftover_graph = self.subgraph_by_rels(leftover_rels)
        isolated_subgraphs.sort(key=len,reverse=True)

        isolated_subgraphs2 = list()
        common_edges_leftovers2 = list()
        for g in common_edges_leftovers:
            isolated_partitions, leftovers = self.split_multi()
            isolated_subgraphs2 += isolated_partitions
            common_edges_leftovers2 += leftovers

        print('Graph was split into %d isolated subgraphs and one leftover graph with %d nodes and %d edges' %
                        (len(isolated_subgraphs),len(leftover_graph),leftover_graph.number_of_edges()))
        #subgraphs.append(leftover_graph)

        return isolated_subgraphs
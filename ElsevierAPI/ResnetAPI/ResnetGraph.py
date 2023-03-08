import networkx as nx
from networkx.exception import NetworkXError
from ..pandas.panda_tricks import df
from datetime import timedelta
import os, math, time
from .rnef2sbgn import minidom
#import xml.etree.ElementTree as et
from lxml import etree as et
import glob
from .NetworkxObjects import PSObject,PSRelation,len,REGULATORS,TARGETS,CHILDS,REFCOUNT
from ..ETM_API.references import pubmed_hyperlink, make_hyperlink
from ..ETM_API.references import PUBYEAR,EFFECT,TITLE,REFERENCE_PROPS,CITATION_INDEX
from ..ETM_API.etm import ETMstat,IDENTIFIER_COLUMN
from itertools import combinations,product

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
        self.urn2rel = dict() #lookup for combining relations in database graph and to_rnef graph

    def copy(self)->'ResnetGraph':
        cp = super().copy(self)
        cp.urn2rel = dict(self.urn2rel)
        return cp

    def load_urn_dicts(self):
        self.urn2rel = {self.rel_urn(rel):rel for r,t,rel in self.edges(data='relation')}

    @staticmethod
    def execution_time(execution_start):
        return "{}".format(str(timedelta(seconds=time.time() - execution_start)))

###############    ADD ADD ADD    #############################
    def add_psobj(self,node:PSObject):
        self.add_nodes_from([(int(node.uid()), node.items())])


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
            for regulator_uid in rel.Nodes[REGULATORS]:
                for target_uid in rel.Nodes[TARGETS]:
                    self.add_edge(regulator_uid[0], target_uid[0], relation=rel,weight=ref_count,key=rel_urn)
                    self.urn2rel[rel_urn] = rel
        else:
            try:
                reg_pairs = combinations(rel.Nodes[REGULATORS],2)
            except KeyError:
                reg_pairs = combinations(rel.Nodes[TARGETS],2)
            [self.add_edge(pair[0][0], pair[1][0], relation=rel,weight=ref_count,key=rel_urn) for pair in reg_pairs]
            [self.add_edge(pair[1][0], pair[0][0], relation=rel,weight=ref_count,key=rel_urn) for pair in reg_pairs]
            self.urn2rel[rel_urn] = rel


    def add_psobj(self,node:PSObject):
        node_uid = node.uid()
        try:
            exist_node = PSObject(self.nodes[node_uid])
            node.merge_obj(exist_node)
            self.add_node(node_uid,**node)
        except KeyError:
            self.add_node(node_uid,**node)
            

    def add_psobjs(self,nodes:set,merge=True):
        '''
        Input
        -----
        nodes - {PSObject} set
        '''
        if merge:
            existing_nodes = list()
            new_nodes = list()

            my_uids = set(self.nodes())
            [(new_nodes,existing_nodes)[n.uid() in my_uids].append(n) for n in nodes]
            [n.merge_obj(self._psobj(n.uid())) for n in existing_nodes]
            nodes2add = existing_nodes+new_nodes
        else:
            nodes2add = nodes

        self.add_nodes_from([(n.uid(),n.items()) for n in nodes2add])


    def add_rel(self,rel:PSRelation,merge=True):
        """
        nodes connected by relation must exist in the graph.\n
        non-directional relations are duplicated in all possible directions\n
        use "add_rel" to move relation from one graph to another or to add nodes together with "rel"
        """
        rel_urn = self.rel_urn(rel) # will raise error if nodes for relation do not exists in graph
        if merge:
            try:
                self.urn2rel[rel_urn].merge_rel(rel)
                self.__add_rel(self.urn2rel[rel_urn])
            except KeyError:
                self.__add_rel(rel)
        else:
            self.__add_rel(rel)


    def add_psrels(self,rels:set(),merge=True):
        """
        Input
        -----
        rels - {PSRelation}
        nodes connected by relation must exist in the graph.\n
        non-directional relations are duplicated in all possible directions\n
        """
        if merge:
            existing_rels = set()
            new_rels = set()
            [(new_rels,existing_rels)[self.rel_urn(rel) in self.urn2rel].add(rel) for rel in rels]
            [rel.merge_rel(self.urn2rel[rel.urn()]) for rel in existing_rels]
            rels2add = existing_rels|new_rels
        else:
            rels2add = rels

        # __add_rel will raise error if nodes for relation do not exists in graph
        [self.__add_rel(rel) for rel in rels2add]
        
        
    def __copy_rel(self,rel:PSRelation,from_graph:"ResnetGraph"):
        '''
        copies rel and its nodes specified in rel based on their URNs rather than database identifiers from ['Id'][0]\n
        re-assigns node identifiers in rel.Nodes to make them consistent with relations in "self" if necessary
        '''
        rel2add = rel.copy()
        nodes2add = list()
        for s, node_tup_list in rel.Nodes.items():
            for i, node_tup in enumerate(node_tup_list):
                from_graph_node_id = node_tup[0]
                my_node = from_graph._psobj(from_graph_node_id)
                try:
                    exist_node = self.urn2node(my_node.urn())
                    exist_node_id = exist_node.uid()
                    rel2add.Nodes[s][i] = (exist_node_id,node_tup[1],node_tup[2])
                except KeyError:
                    nodes2add.append((my_node.uid(),my_node.items()))
                   # self.urn2obj[my_node.urn()] = my_node
                    
        self.add_nodes_from(nodes2add)
        self.add_rel(rel2add)
        return

    
    def add_triple(self,regulator:PSObject,target:PSObject,rel:dict or PSRelation,refs=[],is_directional=True):
        """
        adds nodes and their relations to graph
        """
        self.add_psobjs([regulator,target])

        if isinstance(rel,PSRelation):
            self.add_rel(rel)
            return
        else:
            rel = PSRelation.make_rel(regulator,target,rel,refs,is_directional)
            #rel['Id'] = [self.number_of_edges()]
            self.rel_name(rel)
            rel_urn = self.rel_urn(rel)
            self.add_edge(regulator.uid(), target.uid(), relation=rel,weight=rel.get_reference_count(),key=rel_urn)


    def add_graph(self, other:"ResnetGraph",merge=True):
        '''
        slow. use it only to merge ResnetGraphs from different RNEF files or RNEF ResnetGraph with database ResnetGraph
        node ids and relation ids in self take precedent.\n
        To merge graphs from database use ResnetGraph.compose()
        '''
        self.add_psobjs( other._get_nodes(),merge)
        self.add_psrels(other.__psrels(),merge)

        
    def compose(self,other:"ResnetGraph")->"ResnetGraph":
        '''
        Returns
        -------
        new graph created by nx.compose\n
        other graph attributes take precedent
        '''
        composed_graph = ResnetGraph(nx.compose(self,other))
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
        uid2value = dict() #remapping urn to ids
        for urn, value in urn2value.items():
            try:
                node2annotate = self.urn2node(urn)
                uid2value[node2annotate.uid()] = value
            except KeyError: continue

        nx.set_node_attributes(self,uid2value,new_prop_name)
        #print('%d nodes were annotated with attributes "%s"' % (len(id2value),new_prop_name))


    def add_edge_annotation(self,between_node_uid,and_node_uid,for_rel_with_urn:str,prop_name,prop_values:list):
        annotate_values = set(prop_values)
        if self.has_edge(between_node_uid, and_node_uid,for_rel_with_urn):
            try:
                my_prop_values = self[between_node_uid][and_node_uid][for_rel_with_urn]['relation'][prop_name]
                annotate_values.update(my_prop_values)
            except KeyError:
                pass
            self[between_node_uid][and_node_uid][for_rel_with_urn]['relation'][prop_name] = list(annotate_values)


    def set_edge_annotation(self,between_node_uid,and_node_uid,for_rel_with_urn:str,prop_name,prop_values:list):
        if self.has_edge(between_node_uid, and_node_uid,for_rel_with_urn):
            self[between_node_uid][and_node_uid][for_rel_with_urn]['relation'][prop_name] = prop_values


    def set_rel_annotation(self,for_rel:PSRelation,prop_name,prop_values:list):
        for ruid in for_rel.regulator_uids():
            for tuid in for_rel.target_uids():
                    self[ruid][tuid][for_rel.urn()]['relation'][prop_name] = prop_values


    def add_node_annotation(self, with_new_prop:str, map2prop:str, using_map:dict):
        """
        using_map = {map2prop_value:[annotations]}
        """
        annotation_counter = 0
        for i, node in self._get_nodes():
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
                        nx.set_node_attributes(self, {node.uid():{with_new_prop:list(merged_annotation)}})
                    except KeyError:
                        # case when node has no with_new_prop 
                        nx.set_node_attributes(self, {node.uid():{with_new_prop:list(annotate_with_values)}})
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
        centrality_dic_norm = {uid:v/norm_const for uid,v in centrality_dic.items()}
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


    def rank_regulator(self, regulator:PSObject, target_weights:dict, max_distance=5):
        '''
        Input
        -----
        target_weights = {node_id:weight}

        Returns
        -------
        regulator rank for regulator_id
        '''
        regulator_rank = float(0.0)
        regulator_tree = nx.bfs_tree(self, regulator.uid())

        for level in range (1,max_distance+1):
            targets_on_level = nx.descendants_at_distance(regulator_tree, regulator.uid(), level)
            if not targets_on_level: break
            for target_uid in targets_on_level:
                try:
                    target_weight = float(target_weights[target_uid])
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
            if self.has_node_with_id(Id):
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
            return rel.urn()
        except KeyError:
            regulators, targets = self.find_nodes(for_relation=rel)
            reg_urns = [r.urn() for r in regulators]
            tar_urns = [r.urn() for r in targets]
            return rel.make_urn(reg_urns,tar_urns)
            

    def load_references(self,weight_by_property='',using_value2weight:dict=None):
        """
        Input
        -----
        references are annotated with 'weight_by_property' 'using_value2weight' specifications
        """
        graph_references = set()
        if weight_by_property:
            for regulatorID, targetID, rel in self.edges.data('relation'):
                rel.load_references()
                rel._weight2ref(weight_by_property,using_value2weight)
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
                        count = int(ref_counter[CITATION_INDEX][0])
                        ref_counter[CITATION_INDEX] = [count + 1]
                    except KeyError:
                        ref[CITATION_INDEX] = [1]
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

        clinvar_pmids = [['10447503'],['10592272'],['10612825'],['11125122'],['26619011'],['25741868'],['26582918'],['28492532']]
        clinvar_hyperlinks = list(map(pubmed_hyperlink,clinvar_pmids))
        ref_df = ref_df.remove_rows_by(clinvar_hyperlinks,IDENTIFIER_COLUMN)
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
            regulator = self._psobj(regulator_id)
            target = self._psobj(target_id)
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


    def load_references_between(self,nodes:list,and_nodes:list,weight_by_property='',using_value2weight:dict=None):
        '''
        Input
        -----
        nodes,and_nodes - [PSObject]
        proval2weight = {value:weight} for values in weight_by_prop_name

        Returns
        -------
        {References}
        '''
        sub_graph = self.get_subgraph(nodes, and_nodes)
        return sub_graph.load_references(weight_by_property,using_value2weight)


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
        node_ids = self.uids4nodes(property_values,prop_names)
        self.remove_nodes_from(node_ids)
        print("%d nodes with %s were removed" % (len(node_ids), ','.join(property_values)))


    def remove_nodes_by_degree(self,min_degree=0,max_degree=1000000,only_with_prop:list=['ObjTypeName'],having_values:list=[]):
        if having_values:
            ids2remove = set(self.uids4nodes(having_values,only_with_prop))
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
            ids2remove = set(self.uids4nodes(having_values,only_with_prop))
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
            ids2remove = set(self.uids4nodes(having_values,only_with_prop))
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
                if self.has_edge(pair[0][0], pair[1][0], key=rel_urn):
                    self.remove_edge(pair[0][0], pair[1][0], key=rel_urn)

        self.urn2rel.pop(rel_urn)


    def remove_edges4psobjs(self,with_value:str, in_property='Name'):
        psobjs = self._psobjs_with(with_value,in_property)
        psobjs_uids = self.uids(psobjs)
        edges2remove = [(r,t) for r,t in self.edges() if (r in psobjs_uids) or (t in psobjs_uids)]
        [self.remove_edge(r,t) for r,t in edges2remove]
        print(f'Removed {len(edges2remove)} edges connected to node with {in_property} = {with_value}')


    def clear_resnetgraph(self):
        super().clear()
        self.urn2rel.clear()


######################   GET GET GET   ######################################
    def dbids4nodes(self, with_values=None, in_properties:list=['ObjTypeName']):
        '''
        Return
        ------
        if with_values is not supplied returns dbids for nodes from entire graph 
        '''
        if isinstance(with_values,(set,list)):
            all_dbids = set()
            for node in self._get_nodes():
                for propName in in_properties:
                    if PSObject(node).is_annotated(propName, with_values):
                        all_dbids.add(node.dbid())
                        break
            return list(all_dbids)
        else:
            return [id[0] for i,id in self.nodes(data='Id')]
    
    
    def uids4nodes(self, with_values:list, in_properties=['ObjTypeName']):
        all_ids = set()
        for node in self._get_nodes():
            for propName in in_properties:
                if PSObject(node).is_annotated(propName, with_values):
                    all_ids.add(node.uid())
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



    def urn2node(self, urn:str):
        try:
            node_uid = PSObject.urn2uid(urn)
            return self._psobj(node_uid)
        except KeyError:
            raise KeyError
    

    def urn2id(self, urn:str):
        try:
            my_node = self.urn2node(urn) 
            return my_node.uid()
        except KeyError:
            raise KeyError


    def _get_node(self, with_uid:int):
        return PSObject(self.nodes[with_uid])


    def _get_nodes(self, with_uids=set()):
        '''
        Input
        -----
        with_uids - {uids}

        Return
        ------
        [PSObject] for graph nodes with with_uids or nodes of entire graph if with_uids is empty
        '''
        if with_uids:
            return [PSObject(ddict) for i,ddict in self.nodes(data=True) if i in with_uids] 
        else:
            return [PSObject(ddict) for i,ddict in self.nodes(data=True)]
               
            '''
            #may be faster?
            nodes = list()
            for i in node_ids:
                try:
                    node = self.psobjs_with()
                    node = PSObject(self.nodes(data=True)[i])
                    nodes.append(PSObject(node))
                except KeyError:
                    continue
            return nodes
            '''

    def _psobjs_with(self,value,in_property='Name'): # use this function to find object by Name
        '''
        Input
        -----
        value - str,int,float
        '''
        node_collector = list()
        for node in self._get_nodes():
            try:
                node_property_values = node[in_property]
                if value in node_property_values:
                    node_collector.append(node)
            except KeyError:
                continue
                
        return node_collector


    def psobj_with_dbids(self,dbids:set):
        return [PSObject(o) for i,o in self.nodes(data=True) if o['Id'][0] in dbids]


    def psobjs_with(self, with_properties:list=['ObjTypeName'], only_with_values=[]):
        '''
        Returns
        -------
        [PSObject], if "only_with_values" is empty returns all nodes annotated "with_properties"
        '''
        return_node_uids = set()
        if only_with_values:
            my_values = set(only_with_values)
            for prop_name in with_properties:
                id2prop = nx.get_node_attributes(self,prop_name)
                my_uids = [i for i,v in id2prop.items() if not my_values.isdisjoint(v)]
                return_node_uids.update(my_uids)
        else:
            for prop_name in with_properties:
                id2prop = nx.get_node_attributes(self,prop_name)
                return_node_uids.update(id2prop.keys())
            
        return self._get_nodes(list(return_node_uids)) if return_node_uids else []


    def find_parents4(self,psobjects:list):
        '''
        Returns
        -------
        parents2return - {PSObject}
        childuid2parents - {child_uid:{PSObject}}
        '''
        childuid2parents = dict()
        parents2return = set()
        all_parents = self.psobjs_with([CHILDS])
        for parent in all_parents:
            has_children = set(parent[CHILDS]).intersection(psobjects)
            for child in has_children:
                try:
                    childuid2parents[child.uid()].update(parent)
                except KeyError:
                    childuid2parents[child.uid()] = {parent}
                parents2return.add(child)

        return parents2return, childuid2parents




    def __psrels(self):
        """
        Returns
        -------
        [PSRelation] - list of all relations in graph
        """
        return [r for n1,n2,r in self.edges.data('relation')]
    

    def relation_dbids(self):
        return [r.dbid() for n1,n2,r in self.edges.data('relation')]


    def _psrels4(self, regulator_uid, target_uid):
        try:
            edges = dict(self[regulator_uid][target_uid]).values()
            return [e['relation'] for e in edges]
        except KeyError:
            return []


    def psrels_with(self, with_values:list=[], in_properties:list=['ObjTypeName']):
        '''
        Input
        -----
        if "with_values" is empty will return True if self has any value in "in_properties"
        '''
        relations2return = set()
        if with_values:
            for regulatorID, targetID, rel in self.edges.data('relation'):
                for prop_name in in_properties:
                    if rel.is_annotated(prop_name, with_values):
                        relations2return.add(rel)
                        break

        return relations2return

    
    def __find_relations(self, reg_uid, targ_uid, rel_types=list(), with_effects=list(), mechanism=list(), any_direction=False):
        my_rels = self._psrels4(reg_uid, targ_uid)
        if any_direction:
            my_rels = my_rels + self._psrels4(targ_uid, reg_uid)

        if rel_types:
            my_rels = [rel for rel in my_rels if rel.objtype() in rel_types]

        if with_effects:
            try:
                my_rels = [x for x in my_rels if x.effect() in with_effects]
            except KeyError:
                return list()
        
        if mechanism:
            try:
                my_rels = [x for x in my_rels if x.mechnaism() in mechanism]
            except KeyError:
                return list()

        return my_rels


    def __relation_exist(self,between_ids:list,and_ids:list,with_reltypes=list(),with_effects=list(),mechanism=list(),any_direction=False):
        for node_id1 in between_ids:
            for node_id2 in and_ids:
                my_rels = self.__find_relations(node_id1,node_id2,with_reltypes,with_effects,mechanism,any_direction)
                if my_rels:
                    return True
                else:
                    continue

        return False


    def relation_exist(self,between_psobjects1:list,and_psobjects2:list,with_reltypes=list(),with_effects=list(),mechanism=list(),any_direction=False):
        uids1 = self.uids(between_psobjects1)
        uids2 = self.uids(and_psobjects2)
        for node_id1 in uids1:
            for node_id2 in uids2:
                my_rels = self.__find_relations(node_id1,node_id2,with_reltypes,with_effects,mechanism,any_direction)
                if my_rels:
                    return True
                else:
                    continue

        return False


    def find_relations(self, between_uids,and_uids,with_reltypes=list(),with_effects=list(),mechanism=list(),any_direction=False):
        '''
        Return
        ------
        Empty set if no relation exists between_uids and_uids  
        '''
        my_rels = set()
        for node_uid1 in between_uids:
            for node_uid2 in and_uids:
                my_rels.update(self.__find_relations(node_uid1,node_uid2,with_reltypes,with_effects,mechanism,any_direction))
          
        return my_rels


    def _effect_counts__(self):
        positive_rels = self.psrels_with(['positive'],[EFFECT])
        negative_rels = self.psrels_with(['negative'],[EFFECT])

        positive_refs = set()
        [positive_refs.update(rel._get_refs(sort_by='')) for rel in positive_rels]
        negative_refs = set()
        [negative_refs.update(rel._get_refs(sort_by='')) for rel in negative_rels]

        # hack to count effect when references are not loaded:
        if not positive_refs:
            if not negative_refs:
                return positive_rels,negative_rels

        return positive_refs,negative_refs


    def effect_stats(self, between_nodes:list, and_nodes:list):
        '''
        Input
        -----
        between_nodes, and_nodes - PSObject

        Returns
        -------
        tuple: between_graph - graph between_node_ids and_node_ids\n
        positive_refs - set(Reference) of all relations with Effect = positive\n
        negative_refs - set(Reference) of all relations with Effect = negative
        '''
        between_graph = self.get_subgraph(between_nodes, and_nodes)
        positive_refs, negative_refs = between_graph._effect_counts__()
        return between_graph,positive_refs,negative_refs


    def effect_vote(self,regulator:PSObject,target:PSObject, any_direction=False):
        '''
        Return
        ------
        Consensus Effect sign for the relations between "regulator" and "target":\n
        1 - positive
        -1 - negative
        0 - unknown
        '''
        positive_counter = set()
        negative_counter= set()
        pair_rels = self.__find_relations(
            regulator.uid(),target.uid(),
            with_effects=['positive','negative'],
            any_direction=any_direction
            )

        for rel in pair_rels:
            refs = rel._get_refs(sort_by='')
            if rel.effect() == 'positive':
                positive_counter.update(refs) if refs else positive_counter.add(rel)
            elif rel.effect() == 'negative':
                negative_counter.update(refs) if refs else negative_counter.add(rel)

        if len(positive_counter) > len(negative_counter):
            return 1
        elif len(positive_counter) < len(negative_counter):
            return -1
        return 0


    def net_regulator_effect(self,regulator_id:int,target_ids:list,vote_effect_by_ref=False):
        activated_target_ids = list()
        inhibited_target_ids = list()
        if vote_effect_by_ref:
            for target_id in target_ids:
                between_graph = self.get_subgraph([regulator_id], [target_id],in_direction='>')
                positive_refs, negative_refs = between_graph._effect_counts__()
                if len(positive_refs) > len(negative_refs):
                    activated_target_ids.append(target_id) 
                elif len(positive_refs) < len(negative_refs):
                    inhibited_target_ids.append(target_id)
        else:
            for target_id in target_ids:
                rels = self._psrels4(regulator_id,target_id)
                net_effect = 0
                for rel in rels:
                    try:
                        rel_effect = rel.effect()
                        if rel_effect == 'positive':
                            net_effect += 1
                        elif rel_effect == 'negative':
                            net_effect += -1
                    except KeyError:
                        continue

                if net_effect > 0:
                    activated_target_ids.append(target_id)
                elif  net_effect < 0:
                    inhibited_target_ids.append(target_id)

        return activated_target_ids,inhibited_target_ids


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
        entities1 = self.psobjs_with(with_properties=[search_property],only_with_values=[entity1prop])
        entities2 = self.psobjs_with(with_properties=[search_property],only_with_values=[entity2prop])
        all_entities = set(entities1+entities2)
        if with_children:
            [all_entities.update(e.childs()) for e in all_entities]
        
        all_uids = self.uids(all_entities)
        all_rels = [rel for regulatorID, targetID, rel in self.edges.data('relation') 
                if regulatorID in all_uids and targetID in all_uids]
        
        if filter_rel_with:
            return [rel for rel in all_rels if rel.has_value_in(filter_rel_with)]
        else:
            return all_rels


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
    

    def _relations_dbids(self):
        return {rel.dbid() for regulatorID, targetID, rel in self.edges.data('relation')}


    def find_nodes(self,for_relation:PSRelation,filter_by:list=None,in_properties:list=None):
        """
        Returns
        -------
        regulators = [PSObject], targets = [PSObject]
        """
        regulators_ids = for_relation.regulator_uids()
        regulators = [PSObject(o) for i,o in self.nodes(data=True) if i in regulators_ids]
        target_ids = for_relation.target_uids()
        targets = [PSObject(o) for i,o in self.nodes(data=True) if i in target_ids]

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


    def has_node(self, node:PSObject):
        return super().has_node(node.uid())
    

    def has_node_with_id(self, node_id:int):
        return super().has_node(node_id)

    def _psobj(self, node_uid:int):
        try:
            return PSObject({k:v for k,v in self.nodes[node_uid].items()})
        except KeyError:
            raise KeyError

    def __psobjs(self,with_uids=None):
        '''
        Input
        ----
        returns [PSObject] for all nodes in self "with_uids"
        if "with_uids" is None returns all nodes in self
        '''
        if isinstance(with_uids,(set,list)):
            return {PSObject(n) for u,n in self.nodes(data=True) if u in with_uids}
        else:
            return {PSObject(n) for u,n in self.nodes(data=True)}

    @staticmethod
    def uids(nodes:list):
        '''
        Input
        -----
        nodes = [PSObject]
        '''
        return [n.uid() for n in nodes]
    
    @staticmethod
    def urns(nodes:list):
        '''
        Input
        -----
        nodes = [PSObject]
        '''
        return [n.urn() for n in nodes]
    
    @staticmethod
    def dbids(nodes:list):
        '''
        Input
        -----
        nodes = [PSObject]
        '''
        return [n.dbid() for n in nodes]


    def names(nodes:list):
        '''
        Input
        -----
        nodes = [PSObject]
        '''
        return [n.name() for n in nodes]
    

    def dbid2uid(self,dbids:list=[]):
        '''
        Input
        -----
        nodes = [PSObject]

        Returns
        -------
        {node.dbid():node.uid()}
        '''
        if dbids:
            nodes_with_dbids = self.psobj_with_dbids(dbids)
            return {n.dbid():n.uid() for n in nodes_with_dbids}
        else:
            all_nodes = self._get_nodes()
            return {n.dbid():n.uid() for n in all_nodes}
    

    def uid2dbid(self,for_uids:list=[]):
        '''
        Input
        -----
        nodes = [PSObject]

        Returns
        -------
        {node.uid():node.dbid()}
        '''
        if for_uids:
            nodes_with_uids = self.__psobjs(for_uids)
            return {n.uid():n.dbid() for n in nodes_with_uids}
        else:
            all_nodes = self.__psobjs()
            return {n.uid():n.dbid() for n in all_nodes}
        
    
    def node_urns(self):
        return self.urns(self.__psobjs())
    
    def rel_urns(self):
        return self.urns(self.__psrels())


    def get_neighbors(self,of_nodes:set,allowed_neigbors:list=[]):
        """
        Input
        -----
        of_nodes = {PSObject}
        allowed_neigbors = [PSObject], optional

        Returns
        -------
        list of both upstream and downstream PSObjects
        """
        neighbor_uids = set()
        for n in of_nodes:
            if self.has_node(n):
                n_neighbors_uids = list(nx.all_neighbors(self, n.uid()))
                neighbor_uids.update(n_neighbors_uids)

        if allowed_neigbors:
            allowed_neigbors_uids = self.uids(allowed_neigbors)
            neighbor_uids = neighbor_uids.intersection(allowed_neigbors_uids)
        
        return self.__psobjs(neighbor_uids)
 


    def get_neighbors_rels(self, for_node_ids:set, only_neighbors_with_ids=None):
        neighbor_graph = self.neighborhood(for_node_ids,only_neighbors_with_ids)
        return neighbor_graph.__psrels()

    def get_neighbors_refs4(self, psobjects:set, only_neighbors_with_ids=None):
        neighbor_graph = self.neighborhood(psobjects,only_neighbors_with_ids)
        return neighbor_graph.load_references()

    @staticmethod
    def get_att_set(prop_name:str,ps_objects:list):
        return set([i for sublist in [x[prop_name] for x in ps_objects] for i in sublist])


    @staticmethod
    def children4(parents:list,at_depth=1):
        """
        Input
        -----
        parents in self must have [CHILD] annotation
        use after APISessions::load_children4(parents)

        Return
        ------
        {PSObject}
        """
        children = set()
        level_parents = parents
        for level in range(0,at_depth):
            level_childrens = set()
            [level_childrens.update(p.childs()) for p in level_parents]
            level_parents = level_childrens
            children.update(level_childrens)

        return children


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
            return {x for x,y in self.nodes(data=True) if ((self.out_degree(x) >= min_targets) & (y['ObjTypeName'][0] in only_objtype))}
        else:
            return {x for x in self.nodes() if self.out_degree(x) >= min_targets}


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

    '''
    def __get_ref_list(self, relPropNames:list, entity_prop_names=[], RefNumPrintLimit=0,col_sep:str='\t',single_rel_row=False):
        references = list()
        if not entity_prop_names:
            for regulatorID, targetID, rel in self.edges.data('relation'):
                reference_view_triple = str(rel.triple2str(relPropNames,as1row=single_rel_row))
                references.append(reference_view_triple)
        else:
            # rel = PSRelation()
            for regulatorID, targetID, rel in self.edges.data('relation'):
                reg = self._psobj(regulatorID)
                target = self._psobj(targetID)
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
        '''

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
                        reg = self._psobj(regulatorID)
                        target = self._psobj(targetID)
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
                    n = self._psobj(node_id)
                    node_prop_str = n.data2str(entity_prop_names, col_sep=col_sep)
                    f.write(node_prop_str)


    def ref2pandas (self, relPropNames:list, entity_prop_names=[], RefNumPrintLimit=0,single_rel_row=False) -> 'df':
        temp_fname = '__temp__.tsv'
        self.print_references(temp_fname,relPropNames,entity_prop_names,RefNumPrintLimit=RefNumPrintLimit,single_rel_row=single_rel_row)
        to_return = df.read(temp_fname,header=0,index_col=False, dtype='unicode')
        os.remove(temp_fname)
        return to_return


    def to_rnefstr(self,ent_props:list,rel_props:list,add_rel_props:dict={},add_pathway_props:dict={}):
        # add_rel_props,add_pathway_props structure - {PropName:[PropValues]}
        print_snippets = set(REFERENCE_PROPS).intersection(rel_props)
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
        graph_relations = self.__psrels()
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
            if print_snippets:
                references = list(set(rel.References.values()))
                ref_index = 0
                for ref in references:
                    # each snippet has its own index in RNEF
                    for textref, sentence_props in ref.snippets.items():
                        et.SubElement(xml_control, 'attr',{'name': str('TextRef'), 'value': textref, 'index': str(ref_index)})
                        for sentprop_name, sentprop_values in sentence_props.items():
                            if sentprop_name in print_snippets:
                                for v in sentprop_values:
                                    et.SubElement(xml_control, 'attr',{'name':str(sentprop_name), 'value':str(v), 'index':str(ref_index)})

                        for prop_name, prop_values in ref.items():
                            if prop_name in print_snippets:
                                for prop_value in prop_values:
                                    et.SubElement( xml_control, 'attr',{'name': str(prop_name), 'value': str(prop_value), 'index': str(ref_index)})
                            
                        for ref_id_type,ref_id in ref.Identifiers.items():
                            if ref_id_type in print_snippets:
                                et.SubElement(xml_control, 'attr',{'name':str(ref_id_type), 'value':str(ref_id), 'index':str(ref_index)})
                        ref_index+=1

        xml_str = str(et.tostring(resnet, encoding='utf-8').decode("utf-8"))
        return xml_str
    

    def __2resnet(self,resnet:et.ElementTree,ent_props:list,rel_props:list,add_rel_props:dict={},add_pathway_props:dict={}):
        '''
        Input
        -----
        add_rel_props,add_pathway_props structure - {PropName:[PropValues]}
        '''
        print_snippets = set(REFERENCE_PROPS).intersection(rel_props)
        #resnet = et.Element('resnet')
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
        graph_relations = self.__psrels()
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
            if print_snippets:
                references = list(set(rel.References.values()))
                ref_index = 0
                for ref in references:
                    # each snippet has its own index in RNEF
                    for textref, sentence_props in ref.snippets.items():
                        et.SubElement(xml_control, 'attr',{'name': str('TextRef'), 'value': textref, 'index': str(ref_index)})
                        for sentprop_name, sentprop_values in sentence_props.items():
                            if sentprop_name in print_snippets:
                                for v in sentprop_values:
                                    et.SubElement(xml_control, 'attr',{'name':str(sentprop_name), 'value':str(v), 'index':str(ref_index)})

                        for prop_name, prop_values in ref.items():
                            if prop_name in print_snippets:
                                for prop_value in prop_values:
                                    et.SubElement( xml_control, 'attr',{'name': str(prop_name), 'value': str(prop_value), 'index': str(ref_index)})
                            
                        for ref_id_type,ref_id in ref.Identifiers.items():
                            if ref_id_type in print_snippets:
                                et.SubElement(xml_control, 'attr',{'name':str(ref_id_type), 'value':str(ref_id), 'index':str(ref_index)})
                        ref_index+=1


    def to_rnefstr(self,ent_props:list,rel_props:list,add_rel_props:dict={},add_pathway_props:dict={}):
        resnet = et.Element('resnet')
        self.__2resnet(resnet,ent_props,rel_props,add_rel_props,add_pathway_props)
        xml_str = str(et.tostring(resnet, encoding='utf-8').decode("utf-8"))
        return xml_str


    def write_rnef(self,to_file:str,ent_props:list,rel_props:list,add_rel_props:dict={},add_pathway_props:dict={}):
        with et.xmlfile(to_file) as xf:
            with xf.element('batch'):
                resnet = et.Element('resnet')
                self.__2resnet(resnet,ent_props,rel_props,add_rel_props,add_pathway_props)
                xf.write(resnet,pretty_print=True)
    

    def rnef2sections(self,file_out,ent_prop2print:list,rel_prop2print:list,add_rel_props=dict(),with_section_size=1000):
        """
        resolves printing graphs with and without edges\n
        splits RNEF into <resnet> sections with_section_size <control> elements
        used to write large graphs to file by redcuing the length of xml string
        
        Returns
        -------
        sectioned RNEF string
        """
        if self.number_of_edges():
            resnet_sections_rels = set()
            for regulatorID, targetID, rel in self.edges.data('relation'):
                resnet_sections_rels.add(rel)
                if len(resnet_sections_rels) == with_section_size:
                    section_graph = self.subgraph_by_rels(resnet_sections_rels)
                    rnef_str = section_graph.to_rnefstr(ent_prop2print,rel_prop2print,add_rel_props)
                    rnef_str = str(minidom.parseString(rnef_str).toprettyxml(indent='  '))
                    rnef_str = rnef_str[rnef_str.find('\n')+1:]
                    file_out.write(rnef_str)
                    resnet_sections_rels.clear()

            # printing leftovers
            section_graph = self.subgraph_by_rels(resnet_sections_rels)
            rnef_str = section_graph.to_rnefstr(ent_prop2print,rel_prop2print,add_rel_props)
            rnef_str = str(minidom.parseString(rnef_str).toprettyxml(indent='  '))
            rnef_str = rnef_str[rnef_str.find('\n')+1:]
            file_out.write(rnef_str)
            return 
        else:
            all_nodes = self._get_nodes()
            for sec in range(0, len(all_nodes), with_section_size):
                section_nodes = all_nodes[sec:sec+with_section_size]
                section_graph = ResnetGraph()
                section_graph.add_psobjs(section_nodes)
                rnef_str = section_graph.to_rnefstr(ent_prop2print,rel_prop2print,add_rel_props)
                rnef_str = str(minidom.parseString(rnef_str).toprettyxml(indent='  '))
                rnef_str = rnef_str[rnef_str.find('\n')+1:]
                file_out.write(rnef_str)
            return


    def rnef2file(self,fname:str,ent_prop2print:list,rel_prop2print:list,add_rel_props:dict={},with_section_size=0):
        graph_copy = self.copy() # copying graph to enable using the function in multithreaded file writing
        with open(fname,'w',encoding='utf-8') as f:
            print(f'Writing "{graph_copy.name}" graph to file {fname}')
            if with_section_size:
                print(f'Writing graph "{self.name}" to {fname} file in resnet section of size {with_section_size}')
                f.write('<batch>\n')
                graph_copy.rnef2sections(f,ent_prop2print,rel_prop2print,add_rel_props,with_section_size)
                f.write('</batch>')
            else:
                print(f'Writing graph "{self.name}" to {fname} file in one resnet section')
                self.write_rnef(fname,ent_prop2print,rel_prop2print,add_rel_props)


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
        node1 = self._psobj(between_node_id)
        node2 = self._psobj(and_node_id)
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
                node_psobj = self.urn2node(node_urn)
            except KeyError:
                node_psobj = PSObject({'URN':[node_urn]})

            [node_psobj.update_with_value(attr.get('name'), attr.get('value')) for attr in node.findall('attr')]
            node_psobj['ObjTypeName'] = node_psobj.pop('NodeType')
            self.add_psobj(node_psobj)
            nodel_local_ids[local_id] = node_psobj
                
        for rel in resnet.findall('./controls/control'):
            regulators = list()
            targets = list()
            for link in rel.findall('link'):
                link_type = link.get('type')
                link_ref = link.get('ref')
                link_psobj = nodel_local_ids[link_ref]
                if link_type == 'out': 
                    targets.append(link_psobj)
                else: 
                    regulators.append(link_psobj)

            ps_rel = PSRelation(dict())
            effect = rel.find('./Effect')
            effect_val = 'unknown' if type(effect) == type(None) else effect.text

            for reg in regulators:
                try: 
                    ps_rel.Nodes[REGULATORS].append((reg.uid(), '0', effect_val))
                except KeyError:
                    ps_rel.Nodes[REGULATORS] = [(reg.uid(), '0', effect_val)]
                
            for targ in targets:
                try: 
                    ps_rel.Nodes[TARGETS].append((targ.uid(), '1', effect_val))
                except KeyError:
                    ps_rel.Nodes[TARGETS] = [(targ.uid(), '1', effect_val)]
                        
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
            #ps_rel does not have 'Id' property.This is used by PSRelation.is_from_rnef()
        return


    def __read_rnef(self, rnef_file:str):
        try:
            print ('\nLoading graph from file %s' % rnef_file,flush=True)
            root = et.parse(rnef_file).getroot()  
        except FileNotFoundError:
            raise FileNotFoundError
        except OSError:
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
            #print('Cannot find %s file' % rnef_file)
            raise FileNotFoundError


    @classmethod
    def fromRNEFdir(cls,path2dir:str,merge=True):
        start = time.time()
        listing = glob.glob(path2dir+'*.rnef')
        if listing:
            combo_g = ResnetGraph.fromRNEF(listing[0])
            for i in range(1,len(listing)):
                g = ResnetGraph.fromRNEF(listing[i])
                combo_g.add_graph(g,merge)
            
            print('Graph %d edges and %d nodes was loaded from "%s" with %d files in %s' 
            % (combo_g.number_of_edges(),combo_g.number_of_nodes(),path2dir,len(listing),execution_time(start)))
        else:
            combo_g = ResnetGraph()
            print('Cannot find "%s" directory' % path2dir)
        
        return combo_g
        

    def tree4(self,root:PSObject,reverse=False):
        tree_rn = ResnetGraph()
        try:
            tree = nx.bfs_tree(self,root.uid(),reverse)
            for r,t in tree.edges():
                rels = self._psrels4(r,t)
                [tree_rn.add_rel(rel,self) for rel in rels]
            return tree_rn
        except NetworkXError:
            print(f'{root.name()} nodes was not found in graph {self.name}')
            return ResnetGraph()
        


    def largest_tree(self):
        largest_tree = nx.DiGraph()
        root_nodes = self.root_nodes()
        for root in root_nodes:
            tree = nx.bfs_tree(self,root)
            if len(tree) > len(largest_tree):
                largest_tree = tree

        largest_tree_rn = ResnetGraph()
        for r,t in largest_tree.edges():
            rels = self._psrels4(r,t)
            [largest_tree_rn.add_rel(rel,self) for rel in rels]
        
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
                    path_objs = self.idpath2objs(path)
                    if _is_valid(path_objs):
                        try:
                            end_paths[len(path)].append(path_objs)
                        except KeyError:
                            end_paths[len(path)] = [path_objs]
    
            end_paths_sorted = dict(sorted(end_paths.items()))
            id2paths[end_id] = end_paths_sorted
 
        return id2paths
            

    def idpath2objs(self,uid_path:list, merge_rel2parent=True):
        """
        Input
        -----
        path - list of node ids specifying path
        """
        path_objs = list()
        if merge_rel2parent:
            for i in range(1,len(uid_path)):
                parent_uid = uid_path[i-1]
                child_uid = uid_path[i]
                parent = self._psobj(parent_uid)
                path_rels = self.__find_relations(parent_uid,child_uid)
                [parent.merge_with_rel(r) for r in path_rels]
                path_objs.append(parent)

            path_objs.append(self._psobj(uid_path[-1]))
        else:
            path_objs.append(self._psobj(uid_path[0]))
            for i in range(1,len(uid_path)):
                parent_uid = uid_path[i-1]
                child_uid = uid_path[i]

                parent = self._psobj(parent_uid)
                child = self._psobj(child_uid)
                path_rels = self.__find_relations(parent_uid,child_uid)
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
            '''
            Returns
            -------
            end element in "data" using path specified by "keys"
            '''
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

                for leaf in end_tree_ids:
                    shortest_paths = [p for p in nx.all_simple_paths(self,node.uid(),leaf)]

                    for shortest_path in shortest_paths:
                        curent_path = [value,node_key]
                        for i in range(1,len(shortest_path)):
                            parent_id = shortest_path[i-1]
                            child_id = shortest_path[i]
                            
                            path_rels = self.__find_relations(parent_id,child_id)         
                            for rel in path_rels:
                                path_rel_key = rel_key_prop+':'+rel._prop2str(rel_key_prop)
                                if not path_rel_key: break

                                path_rel_props = rel.props2dict(rel_props)
                                curent_path = add_element(dict2return,curent_path,path_rel_props,path_rel_key)

                                child_node = self._psobj(child_id)
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
    def subgraph(self,node_uids:list):
        '''
        Returns subgraph made by nx.subgraph containg all edges between node_ids
        '''
        sub_g = ResnetGraph(nx.subgraph(self,node_uids))
        sub_g.load_urn_dicts()
        return sub_g


    def neighborhood(self,psobjects:set,only_neighbors:list=[],by_relation_types:list=[],with_effects:list=[],in_direction=''):
        '''
        Input
        -----
        of_nodes, only_neighbors - [PSObject]
        '''
        neighbors = self.get_neighbors(psobjects,only_neighbors)
        return self.get_subgraph(psobjects,neighbors,by_relation_types,with_effects,in_direction)

   
    def subtract(self, other: "ResnetGraph"):
        #only self graph is analyzed 
        # works faster if other graph is bigger than self
        unique2self = ResnetGraph()
        edges_from_other = set(other.__psrels())
        for n1,n2,e in self.edges(data='relation'):
            if e not in edges_from_other:
                unique2self.add_psobjs([self.nodes[n1],self.nodes[n2]],merge=False)
                unique2self.add_edge(n1, n2, relation=e['relation'],weight=e['weight'], key=self.rel_urn(e['relation']))
        return unique2self


    def intersect (self, other: "ResnetGraph"):
        #only self graph is analyzed 
        # works faster if other graph is bigger than self
        intersection = ResnetGraph()
        edges_from_other = other.__psrels()
        for n1,n2,e in self.edges(data='relation'):
            if e in edges_from_other:
                intersection.add_psobjs([self.nodes[n1],self.nodes[n2]],merge=False)
                intersection.add_edge(n1, n2, relation=e['relation'],weight=e['weight'],key=self.rel_urn(e['relation']))
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
  
                    
    def regulatory_network_urn(self,for_targets_with_urns:list,network_name:str):
        """
        Input
        -----
        for_targets - [PSObject]

        Returns 
        -------
        subgraph with targets from 'for_targets_with_urns' and their regulators\n 
        function is used for SNEA
        """
        start = time.time()
        input_target_uids = {PSObject.urn2uid(urn) for urn in for_targets_with_urns}

        rel2copy = [e for r,t,e in self.edges(data='relation') if t in input_target_uids]
        reg_subgraph = self.subgraph_by_rels(rel2copy)      
        reg_subgraph.name = network_name
        print('%s subnetwork with %d edges was selected in %s from network with %d edges' % 
              (network_name, reg_subgraph.number_of_edges(), self.execution_time(start),self.number_of_edges()))
        return reg_subgraph
    

    def regulatory_network(self,for_targets:list,network_name:str):
        """
        Input
        -----
        for_targets - [PSObject]

        Returns 
        -------
        subgraph with targets from 'for_targets_with_urns' and their regulators\n 
        function is used for SNEA
        """
        start = time.time()
        input_target_uids = set(self.uids(for_targets))

        rel2copy = [e for r,t,e in self.edges(data='relation') if t in input_target_uids]
        reg_subgraph = self.subgraph_by_rels(rel2copy)      
        reg_subgraph.name = network_name
        print('%s subnetwork with %d edges was selected in %s from network with %d edges' % 
              (network_name, reg_subgraph.number_of_edges(), self.execution_time(start),self.number_of_edges()))
        return reg_subgraph


    def clean_version_number(self,in_property:str):
        uid2value = dict()
        for uid, psobj in self.nodes(data=True):
            try:
                identifiers = psobj[in_property]
                no_version_set = set()
                for i in identifiers:
                    version_pos = str(i).rfind('.')
                    id_end = version_pos if version_pos > 0 else len(i)
                    no_version_set.add(i[:id_end])
                uid2value[uid] = list(no_version_set)
            except KeyError:
                continue

        graph2return = self.copy()
        nx.set_node_attributes(graph2return,uid2value,in_property)
        return graph2return


    def make_simple(self, rel_type_rank:list=[]):
        """
        Input
        -----
        specify rel_type_rank if graph contains relations with different types connecting the same pair of nodes\n
        Examples:\n
        ['DirectRegulation','Binding','ProtModification','Regulation']\n
        ['PromoterBinding','Expression','Regulation']\n
        ['Biomarker','StateChange']\n
        ['Biomarker','QuantitativeChange']\n
        [MolSynthesis','Regulation']\n
        [MolTransport','Regulation']\n
        [MolTransport','CellExpression']\n
        ['Regulation','FunctionalAssosiation']\n
        
        Returns
        -------
        graph with only one edge between nodes.\n
        keeps relation with the biggest reference count.\n
        all other relations are merged into the most referenced one
        if rel_type_rank is specified the new relation type is assigned accordingly 
        """
        def best_effect_index(rels:list):
            rels.sort(key=lambda x: x[REFCOUNT][0], reverse=True)
            for index,rel in enumerate(rels):
                if rel.effect() != 'unknown':
                    return index
                else:
                    continue
            
            # to merge Binding to DirectRegulation
            for index,rel in enumerate(rels):
                if rel.is_directional():
                    return index
            return 0


        def find_best_type(rels:list):
            if rel_type_rank:
                my_rel_types = {r.objtype() for r in rels}
                for rel_type in rel_type_rank:
                    if rel_type in my_rel_types:
                        return 'DirectRegulation' if rel_type == 'Binding' else rel_type
            return ''
        print(f'Simplifying {self.name}')
        simple_g = self.copy()
        for regulator_uid, target_uid in self.edges():
            reg2target_rels = simple_g._psrels4(regulator_uid,target_uid)
            if len(reg2target_rels) > 1:
                best_rel_index = best_effect_index(reg2target_rels)
                best_rel_type = find_best_type(reg2target_rels)

                best_rel = reg2target_rels[best_rel_index].copy()
                if best_rel_type:
                    simple_g.set_rel_annotation(best_rel,'ObjTypeName',[best_rel_type])
                    
                reg2target_rels.pop(best_rel_index)
                [best_rel._add_refs(rel.References.values()) for rel in reg2target_rels]
                [simple_g.remove_relation(rel) for rel in reg2target_rels]

        print('%d redundant edges in graph "%s" were removed by simplification' % 
              (self.number_of_edges()-simple_g.number_of_edges(),self.name))
        return simple_g
    
                    
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
            if target_ids:
                targets = self._get_nodes(target_ids)
                subnetworks[regulator_id] = targets
                len_targets[regulator_id] = len(target_ids)
        
        nx.set_node_attributes(self, len_targets, NUMBER_OF_TARGETS)

        print('Generated %d regulome subnetworks with more than %d targets from %s' 
                        % (len(subnetworks), min_size, self.name))
        return subnetworks


    def upstream_relations(self,node_id:int,with_types=list()):
        my_rels = [rel for r,t,rel in self.edges.data('relation') if t == node_id]
        if with_types:
            my_rels = [r for r in my_rels if r['ObjTypeName'] in with_types]
        return my_rels


    def find_targets(self,of_regulators:list=[],targets_objtype:list=[],linkedby_reltypes:list=[], min_regulators=1):
        '''
        Input
        -----
        of_regulators - [PSObject]

        Return
        ------
        tuple [PSObject] [PSRelation] [PSObject], where\n
        0 - regulators
        1 - relations between regulators and targets
        2 - targets
        '''
        regulator_uids = set()
        target_uids = set()
        if of_regulators:
            regulators_uids = self.uids(of_regulators)
            targets_uids = [t for r,t in self.edges() if r in regulators_uids]
        else:
            for r,t in self.edges():
                regulator_uids.add(r)
                target_uids.add(t)

        if targets_objtype:
            targets = self.__psobjs(target_uids)
            targets = {t for t in targets if t.objtype() in targets_objtype}
            target_uids = self.uids(targets)

        relation2targets = set()
        if linkedby_reltypes:
            filtered_regulator_uids = set()
            filtered_target_uids = set()
            for r,t,rel in self.edges.data('relation'):
                if rel.objtype() in linkedby_reltypes:
                    if r in regulator_uids:
                        if t in targets_uids:
                            filtered_regulator_uids.add(r)
                            filtered_target_uids.add(t)
                            relation2targets.add(rel)
            regulator_uids = filtered_regulator_uids
            target_uids = filtered_target_uids
        else:
            [relation2targets.add(rel) for r,t,rel in self.edges.data('relation') if r in regulator_uids and t in target_uids]

        if min_regulators > 1:
            reg2target_subgraph = self.subgraph_by_rels(relation2targets)
            target_uids = [uid for uid in reg2target_subgraph.nodes() if reg2target_subgraph.in_degree(uid) > min_regulators]
            filtered_regulator_uids = set()
            filtered_rels = set()
            for r,t,rel in self.edges.data('relation'):
                if t in targets_uids:
                    filtered_regulator_uids.add(r)
                    filtered_rels.add(rel)
            regulator_uids = filtered_regulator_uids
            relation2targets = filtered_rels
        
        return self.__psobjs(regulator_uids), relation2targets, self.__psobjs(target_uids)
    

    def downstream_relations(self,node_id:int,with_types=list()):
        my_rels = [rel for r,t,rel in self.edges.data('relation') if r == node_id]
        if with_types:
            my_rels = [r for r in my_rels if r['ObjTypeName'] in with_types]
        return my_rels


    def upstream_regulators(self,of_node_id:int,linkedby_reltypes=list()):
        regulator_ids = [r for r,t,rel in self.edges.data('relation') if t == of_node_id and rel['ObjTypeName'][0] in linkedby_reltypes]
        return regulator_ids


    def get_regulome(self, start_nodes:set):
        """
        Input
        -----
        start_nodes - {PSobject}

        Returns 
        -------
        composition of bfs_trees for all ids in start_node_ids
        """
        all_trees = nx.DiGraph()
        for node in start_nodes:
            if not self.has_node(node): 
                continue
            t = nx.bfs_tree(self, node.uid())
            all_trees = nx.compose(all_trees, t)
        
        regulome_edges = all_trees.edges()
        regulome_rels = [rel for r,t,rel in self.edges.data('relation') if (r,t) in regulome_edges]
        return self.subgraph_by_rels(regulome_rels)


    def subgraph_by_rel_urns(self, rel_urns:set):
        """
        #not tested!!
        """
        subgraph = ResnetGraph()
        rels2add = [r for u,r in self.urn2rel.items() if u in rel_urns]
        [subgraph.add_rel(rel,self) for rel in rels2add]
        return subgraph


    def subgraph_by_rels(self, rels:list)->'ResnetGraph':
        """
        Input
        -----
        rels = [PSRelation]
        """
        subgraph = ResnetGraph()
        if rels:
            rels_uids = set()
            [rels_uids.update(rel.entities_uids()) for rel in rels]
            rel_nodes = self._get_nodes(rels_uids)
            subgraph.add_psobjs(set(rel_nodes))
            subgraph.add_psrels(rels)
        return subgraph


    def subgraph_by_relprops(self, search_values:list, in_properties:list=['ObjTypeName']):
        '''
        Return
        ------
        if "search_values" is empty returns subgraph with all relations annotated by "in_properties"
        '''
        search_value_set = set(search_values)
        subgraph = ResnetGraph()
        my_rels = set()
        if search_values:
            for prop_type in in_properties:
                for n1, n2, rel in self.edges.data('relation'):
                    try:
                        if not set(rel[prop_type]).isdisjoint(search_value_set):
                            my_rels.add(rel)
                            #subgraph.add_rel(rel,self)
                    except KeyError:
                        continue
        else:
            for prop_name in in_properties:
                for n1, n2, rel in self.edges.data('relation'):
                    try:
                        my_rel = rel[prop_type]
                        my_rels.add(my_rel)
                        #subgraph.add_rel(my_rel,self)
                    except KeyError:
                        continue

        print(f'Select subgraph with {len(my_rels)} relations with {in_properties} from graph with {self.number_of_edges()} edges')
        return self.subgraph_by_rels(my_rels)


    def subgraph_by_nodeprops(self,has_properties:list,with_values=set()):
        '''
        Returns
        -------
        neighborhood of nodes that has_properties if with_values is empty\n
        otherwise neighborhood of nodes that has_properties with_values
        '''
        my_nodes = self.psobjs_with(has_properties,with_values)
        my_nodes_ids = [n.uid() for n in my_nodes]
        return self.neighborhood(my_nodes_ids)


    def get_subgraph(self,between_nodes:list,and_nodes:list,by_relation_types:list=[],with_effect:list=[],in_direction='',urn_ids=False)->"ResnetGraph":
        '''
        Input
        -----
        between_nodes,and_nodes - [PSObject]
        in_direction in ['>','<',None], defaults to None
        '''
        reg_uids = self.uids(between_nodes)
        tar_uids = self.uids(and_nodes)

        edges = set(self.edges())
        my_rels = set()
        if not in_direction or in_direction == '>':
            edges4subgraph = set(product(reg_uids,tar_uids)).intersection(edges)
            rels4subgraph = {self[reg][targ][urn]['relation'] for (reg,targ) in edges4subgraph for urn in self[reg][targ]}
            my_rels.update(rels4subgraph)
        
        if not in_direction or in_direction == '<':
            edges4subgraph = set(product(tar_uids,reg_uids)).intersection(edges)
            rels4subgraph = {self[reg][targ][urn]['relation'] for (reg,targ) in edges4subgraph for urn in self[reg][targ]}
            my_rels.update(rels4subgraph)

        if by_relation_types:
            my_rels = {r for r in my_rels if r.objtype() in by_relation_types}

        if with_effect:
            my_rels = {r for r in my_rels if r.effect() in with_effect}

        subgraph = self.subgraph_by_rels(my_rels)
        return subgraph

 
    def __ontology_resnet(self,children:list,add2parent:PSObject):
        """
        Input
        -----
        members - [PSObjects]
        """
        resnet = ResnetGraph()
        resnet.add_psobj(add2parent)
        resnet.add_psobjs(children)
        for child in children:
            #child_id = m.uid()
            rel = PSRelation({'ObjTypeName':['MemberOf'],'Relationship':['is-a'],'Ontology':['Pathway Studio Ontology']})
            rel.Nodes[REGULATORS] = [(child.uid(),0,0)]
            rel.Nodes[TARGETS] = [(add2parent.uid(),1,0)]
            resnet.add_rel(rel,self)

        return resnet


    def ontology_graph(self):
        '''
        Input
        -----
        nodes in self are annotated with [CHILD_UIDS] property
        '''
        ontology_graph = ResnetGraph()
        parents = self.psobjs_with([CHILDS])
        [ontology_graph.add_graph(self.__ontology_resnet(p[CHILDS],p)) for p in parents]
        return ontology_graph


    def all_paths_from(self,child:PSObject):
        '''
        self must be ontology_graph
        returns [[PSObject]] - list of all pathes to node_id
        '''
        parent_tree = self.tree4(child)
        top_parent_uids = {x for x in parent_tree.nodes() if not parent_tree.out_degree(x)}
        all_ontology_paths = list()
        for parent_uid in top_parent_uids:
            ontology_uid_paths = list(nx.all_simple_paths(parent_tree,child.uid(),parent_uid))
            for path in ontology_uid_paths:
                ontology_obj_path = self.idpath2objs(list(path))
                all_ontology_paths.append(ontology_obj_path)

        return all_ontology_paths


    def direct_indirect_targets(self,of_node_id,
    direct_reltypes:list=['DirectRegulation','Binding','ChemicalReaction','ProtModification','PromoterBinding'],
    indirect_reltypes=['Regulation','MolTransport','Expression','MolSynthesis']):
        direct_neighbors_ids = set()
        indirect_neighbors_ids = set()
        for neighbor_id in self[of_node_id]:
            for urn in self[of_node_id][neighbor_id]:
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

        common_rels = set(graph.__psrels()).difference(partion1_rels|partion2_rels)
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
        all_relations = set(self.__psrels())
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
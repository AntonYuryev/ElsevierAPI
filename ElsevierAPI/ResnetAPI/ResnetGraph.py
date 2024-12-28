import networkx as nx
from networkx.exception import NetworkXError
from ..pandas.panda_tricks import df,np
from datetime import timedelta
import os, math, time, torch,glob,csv
from xml.dom import minidom
from lxml import etree as et
from .NetworkxObjects import PSObject,PSRelation,len, DIRECT, INDIRECT, DBID
from .NetworkxObjects import REGULATORS,TARGETS,CHILDS,REFCOUNT,STATE,DIRECT_RELTYPES,OBJECT_TYPE
from ..ETM_API.references import Reference, pubmed_hyperlink, make_hyperlink
from ..ETM_API.references import PUBYEAR,EFFECT,TITLE,REFERENCE_PROPS,INT_PROPS,PS_CITATION_INDEX,SENTENCE_PROPS,SENTENCE
from ..ETM_API.RefStats import RefStats,IDENTIFIER_COLUMN
from itertools import product
from typing import Optional
from torch_geometric.data import HeteroData
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor,as_completed
from ..utils import execution_time, execution_time2,str2str,unpack,normalize

RESNET = 'resnet'
PHYSICAL_INTERACTIONS = ['Binding','DirectRegulation','ProtModification','PromoterBinding','ChemicalReaction']
NONDIRECTIONAL_RELATIONS = ['Binding','FunctionalAssociation','Paralog','Metabolization','CellExpression']
PROTEIN_TYPES = ['Protein','FunctionalClass','Complex']
ANATOMICAL_CONCEPTS = ['Cell','Organ','Tissue']
PS_REF_COULUMN = 'Number of reference. Link opens recent publications in PubMed'
NUMBER_OF_TARGETS = '# targets'
CLOSENESS = 'Closeness'
CONSISTENCY = 'Consistency coefficient'
NUMERICAL_PROPS = [CLOSENESS]+list(INT_PROPS)
SENTENCE_PROPS_SET = set(SENTENCE_PROPS+['TextRef'])

RNEF_DISCLAIMER = str('Disclaimer: please refer to our Terms and Conditions on authorized use of Elsevier data. https://www.elsevier.com/legal/elsevier-website-terms-and-conditions?dgcid=RN_AGCM_Sourced_300005028')
MAX_RNEF_THREADS = 4


class ResnetGraph (nx.MultiDiGraph):
  pass

  def __init__(self, *args, **kwargs):
      super().__init__(*args, **kwargs)
      self.urn2rel = dict() #lookup for combining relations in database graph and to_rnef graph

  def copy(self)->'ResnetGraph':
      cp = ResnetGraph(super().copy())
      cp.urn2rel = dict(self.urn2rel)
      return cp

  def load_urn_dicts(self):
      self.urn2rel = {rel['relation'].urn():rel['relation'] for r,t,rel in self.edges(data=True)}

  @staticmethod
  def execution_time(execution_start):
      return "{}".format(str(timedelta(seconds=time.time() - execution_start)))

###############    ADD ADD ADD    #############################
  def set_node_attributes(self,values,name):
    nx.function.set_node_attributes(self,values,name)
    for ruid,tuid,rel in self.edges.data('relation'):
      if REGULATORS in rel.Nodes:
        for i,reg in enumerate(rel.Nodes[REGULATORS]):
          rel.Nodes[REGULATORS][i] = self._get_node(reg.uid())
      if TARGETS in rel.Nodes:
        for i,tar in enumerate(rel.Nodes[TARGETS]):
          rel.Nodes[REGULATORS][i] = self._get_node(tar.uid())


  def property2node(self, node_uid,prop_name:str,prop_values:list):
      self.set_node_attributes(self,{node_uid:{prop_name:prop_values}})


  def copy_node_annotation(self, from_property:str, in_other:"ResnetGraph", as_propname=''):
      urn2value = {n.urn():n[from_property] for i,n in in_other.nodes(data=True)}
      prop_name = as_propname if as_propname else from_property
      self.set_node_annotation(urn2value, prop_name)


  def __add_rel(self, rel:PSRelation,refresh_urn=False):
      """
      # nodes connected by rel must exist in the graph
      non-directional relations are duplicated in all possible directions
      """
      rel_urn = rel.urn(refresh_urn)
    #   if rel_urn == 'urn:agi-Binding:in-out:urn:agi-cas:165689-82-7:in-out:urn:agi-llid:834489':
    #       print('')
      uid_pairs = rel.get_regulators_targets()
      if uid_pairs:
          [self.add_edge(pair[0], pair[1], relation=rel,
                        weight=rel.count_refs(),key=rel_urn) for pair in uid_pairs]   
          self.urn2rel[rel_urn] = rel
      else:
          print(f'No nodes exist in the graph for relation {rel_urn} or relation is self-loop')


  def add_psobj(self,node:PSObject):
      node_uid = node.uid()
      try:
          exist_node = PSObject(self.nodes[node_uid])
          merged_node = node.merge_obj(exist_node)
          self.add_node(node_uid,**merged_node)
      except KeyError:
          self.add_node(node_uid,**node)
          

  def add_psobjs(self,nodes:set[PSObject],merge=True):
      '''
      Input
      -----
      nodes - {PSObject} set
      '''
      if not nodes: return
      nodes2add = nodes
      if merge:
          existing_nodes = list()
          new_nodes = list()
          my_uids = set(self.nodes())
          if my_uids:
              [(new_nodes,existing_nodes)[n.uid() in my_uids].append(n) for n in nodes]
              existing_nodes = [n.merge_obj(self._psobj(n.uid())) for n in existing_nodes]
              nodes2add = existing_nodes+new_nodes

      self.add_nodes_from([(n.uid(),n.items()) for n in nodes2add])


  def relations(self):
      for r,t,urn in self.edges(keys=True):
          source, target = PSObject(r), PSObject(t)
          rel = self[source][target][urn]['relation']
          yield source, target, rel


  def add_rel(self,rel:PSRelation,merge=True):
      """
      nodes connected by relation must exist in the graph.\n
      non-directional relations are duplicated in all possible directions\n
      use "__copy_rel" to move relation from one graph to another or to add nodes together with "rel"
      """
      rel_urn = rel.urn() # will raise error if nodes for relation do not exists in graph
      if merge:
          try:
              self.urn2rel[rel_urn] = self.urn2rel[rel_urn].merge_rel(rel)
              self.__add_rel(self.urn2rel[rel_urn])
          except KeyError:
              self.__add_rel(rel)
      else:
          self.__add_rel(rel)


  def add_psrels(self,rels:set[PSRelation],merge=True):
      """
      Input
      -----
      rels - {PSRelation}
      nodes connected by relation must exist in the graph.\n
      non-directional relations are duplicated in all possible directions\n
      """
      rel_nodes = [x.regulators()+x.targets() for x in rels]
      rel_nodes = unpack(rel_nodes)
      self.add_psobjs(rel_nodes)
      if merge:
          existing_rels = set()
          new_rels = set()
          [(new_rels,existing_rels)[rel.urn() in self.urn2rel].add(rel) for rel in rels]
          existing_rels = [rel.merge_rel(self.urn2rel[rel.urn()]) for rel in existing_rels]
          rels2add = set(existing_rels)|new_rels
      else:
          rels2add = rels
      # __add_rel will raise error if nodes for relation do not exists in graph
      [self.__add_rel(rel) for rel in rels2add]
  

  @classmethod
  def from_rels(cls,rels:set[PSRelation]|list[PSRelation]):
      newG = ResnetGraph()
      my_rels = list(rels) if isinstance(rels,set) else list(set(rels))
      newG.add_psrels(list(my_rels))
      return newG


  def __copy_rel(self,rel:PSRelation):
      '''
      copies rel and its nodes specified in rel if they do not exist in self. 
      Otherwise copies node from self to rel  
      '''
      rel2add = rel.copy()
      nodes2add = list()
      for s, node_list in rel.Nodes.items():
          for i, node_from_graph in enumerate(node_list):
              node_from_graph_uid = node_from_graph.uid()
              if node_from_graph_uid in self.nodes():
                  rel2add.Nodes[s][i] = self._psobj(node_from_graph_uid)
              else:
                  nodes2add.append(node_from_graph)
                  
      self.add_psobjs(nodes2add)
      self.add_rel(rel2add)
      return

  
  def add_triple(self,regulator:PSObject,target:PSObject,props_or_rel:(dict|PSRelation),refs=[],is_directional=True):
      """
      adds nodes and their relations to graph
      """
      self.add_psobjs({regulator,target})

      if isinstance(props_or_rel,PSRelation):
          self.add_rel(props_or_rel)
          return
      else:
          assert(isinstance(props_or_rel,dict))
          rel = PSRelation.make_rel(regulator,target,props_or_rel,refs,is_directional)
          #self.rel_name(rel)
          self.__add_rel(rel)
        

  def add_graph(self, other:"ResnetGraph",merge=True):
      '''
      slow. use only to merge ResnetGraphs from different RNEF files 
      or RNEF ResnetGraph with database ResnetGraph
      or update ResnetGraph with database properties

      node ids and relation ids in self take precedent.\n
      To merge graphs from database use ResnetGraph.compose()
      '''
      self.add_psobjs(set(other._get_nodes()),merge)
      self.add_psrels(other._psrels(),merge)

      
  def compose(self,other:"ResnetGraph")->"ResnetGraph":
      '''
      Returns
      -------
      new graph created by nx.compose\n
      other graph attributes take precedent
      '''
      composed_graph = ResnetGraph(nx.operators.binary.compose(self,other))
      composed_graph.urn2rel = dict(self.urn2rel)
      composed_graph.urn2rel.update(other.urn2rel)
      return composed_graph
  

  def clone_node(self,n:PSObject,replace_with:PSObject,flip_effect=0,set_reltype=''):
    '''
    Input
    -----
    flip_effect - ACTIVATED, REPRESSED, UNKNOWN_STATE (1,-1,0)
    if flip_effect == 0 assigns "unknown" Effect to cloned relation
    '''
    new_rels = list()
    for rel in self.get_neighbors_rels({n}):
      new_rel = rel.copy()
      if flip_effect < 0:
        new_rel.flip_effect()
      elif not flip_effect:
        new_rel[EFFECT] = ['unknown']

      if set_reltype:
        new_rel[OBJECT_TYPE] = [set_reltype]

      for nodes in new_rel.Nodes.values():
        if n in nodes:
          index = nodes.index(n) 
          nodes[index] = replace_with

      new_rel.urn(refresh=True)
      new_rels.append(new_rel)
    self.add_psrels(new_rels,merge=False)
    return new_rels


################## SET SET SET ##########################################
  def set_node_annotation(self, urn2values:dict[str,list], new_prop_name:str):
      """
      Input
      -----
      urn2values = {urn:[with_prop_values]}

      Adds
      ----
      new property to exsisting nodes. Existing values of 'new_prop_name' will be replaced
      """
      uid2values = {PSObject.urn2uid(k):v for k,v in urn2values.items()}
      nx.function.set_node_attributes(self,uid2values,new_prop_name)
      # set_node_attributes() does not update nodes that do not exist in Graph
      return


  def add_edge_annotation(self,between_node_uid,and_node_uid,for_rel_with_urn:str,prop_name,prop_values:list):
      annotate_values = set(prop_values)
      if self.has_edge(between_node_uid, and_node_uid,for_rel_with_urn):
          try:
              my_prop_values = self[between_node_uid][and_node_uid][for_rel_with_urn]['relation'][prop_name]
              annotate_values.update(my_prop_values)
          except KeyError:
              pass
          
          self[between_node_uid][and_node_uid][for_rel_with_urn]['relation'][prop_name] = list(annotate_values)


  def set_edge_annotation(self,between_node_uid:int,and_node_uid:int,for_rel_with_urn:str,prop_name,prop_values:list):
      if self.has_edge(between_node_uid, and_node_uid,for_rel_with_urn):
          self[between_node_uid][and_node_uid][for_rel_with_urn]['relation'][prop_name] = prop_values


  def refprop2rel(self,ref_prop:str, relprop:str,min_max=0):
      for r,t,urn in self.edges(keys=True):
          self[r][t][urn]['relation']._refprop2rel(ref_prop, relprop,min_max)


  def set_rel_annotation(self,for_rel:PSRelation,prop_name,prop_values:list):
      for ruid,tuid in for_rel.get_regulators_targets():
          self[ruid][tuid][for_rel.urn()]['relation'][prop_name] = prop_values


  def add_rel_annotation(self,for_rel:PSRelation,prop_name,prop_values:list):
      for ruid,tuid in for_rel.get_regulators_targets():
          self.add_edge_annotation(ruid,tuid,for_rel.urn(),prop_name,prop_values)


  def add_refs2(self,rel:PSRelation,refs:list[Reference]):
      for ruid,tuid in rel.get_regulators_targets():
          self[ruid][tuid][rel.urn()]['relation']._add_refs(refs)


  def merge_rel(self,rel:PSRelation,with_rels:list[PSRelation]):
      for ruid,tuid in rel.get_regulators_targets():
          for r in with_rels:
              merged_rel = self[ruid][tuid][rel.urn()]['relation'].merge_rel(r)
              self[ruid][tuid][rel.urn()]['relation'] = merged_rel


  def add_node_annotation(self, with_new_prop:str, map2prop:str, using_map:dict):
      """
      using_map = {map2prop_value:[annotations]}
      """
      annotation_counter = 0
      for i, node in self._get_nodes():
          if map2prop in node:
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
                  if with_new_prop in node:
                      # case when node has with_new_prop 
                      merged_annotation = set(node[with_new_prop]) | annotate_with_values
                      nx.function.set_node_attributes(self, {node.uid():{with_new_prop:list(merged_annotation)}})
                  else:
                      # case when node has no with_new_prop 
                      nx.function.set_node_attributes(self, {node.uid():{with_new_prop:list(annotate_with_values)}})

      print('%d nodes were annotated "%s" values out of %d "%s" values used for mapping' %
              (annotation_counter,with_new_prop,len(using_map), map2prop) )


  def rename_rel_property(self, oldPropertyName='MedlineTA', newPropertyName='Journal'):
      for ruid,tuid,rel in self.edges.data():
          try:
              rel[newPropertyName] = rel['relation'].pop(oldPropertyName)
          except KeyError:
              for prop in rel.PropSetToProps.values():
                  try:
                      prop[newPropertyName] = prop.pop(oldPropertyName)
                  except KeyError:
                      continue
              continue


############################   LABELING LABELS LABELING   ######################################

  def closeness(self)->dict[int,float]:
      """
      Annotates
      ---------
      nodes with property 'Closeness' calculated by nx.closeness_centrality -\n 
      average length of the shortest path between the node and all other nodes in the graph
      
      Returns
      -------
      {node_id:Closeness}
      """
      #closness = nx.algorithms.centrality.closeness_centrality(self)
      return dict(nx.algorithms.centrality.closeness_centrality(self)) 


  def centrality(self):
      '''
      Returns
      -------
      nodes : dictionary
      Dictionary of nodes with degree centrality as the value.
      '''
      return nx.algorithms.centrality.degree_centrality(self)


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
      regulator_tree = nx.algorithms.traversal.breadth_first_search.bfs_tree(self, regulator.uid())

      for level in range (1,max_distance+1):
          targets_on_level = nx.algorithms.traversal.breadth_first_search.descendants_at_distance(regulator_tree, regulator.uid(), level)
          if not targets_on_level: break
          for target_uid in targets_on_level:
              try:
                  target_weight = float(target_weights[target_uid])
                  regulator_rank = regulator_rank + target_weight/math.pow(level,2)
              except KeyError:
                  continue
      return regulator_rank


  def rank_regulators(self,node_weights:dict[int,float],add2prop:str,max_distance=5):
      '''
      Input
      -----
      node_weights = {node_id:weight}\n
      add2prop - node property name to write regulator rank to
      
      Returns
      -------
      {regulator_id:rank}
      '''
      regulator_ranks = dict()
      regulator_trees = nx.DiGraph()
      for Id in node_weights.keys():
          if self.has_node_with_id(Id):
              source_weight = node_weights[Id]
              tree = nx.algorithms.traversal.breadth_first_search.bfs_tree(self, Id, reverse=True)
              regulator_trees = nx.algorithms.operators.binary.compose(regulator_trees,tree)
              scored_regulators = set()
              for level in range (1,max_distance+1):
                  regulators = nx.algorithms.traversal.breadth_first_search.descendants_at_distance(tree, Id, level)
                  regulators_on_level = regulators.difference(scored_regulators)
                  for regulator_id in regulators_on_level:
                      try:
                          rank = float(regulator_ranks[regulator_id])
                          regulator_ranks[regulator_id] = rank + source_weight/math.pow(level,2)
                      except KeyError:
                          regulator_ranks[regulator_id] = source_weight/math.pow(level,2)
                  scored_regulators.update(regulators)
                  
      nx.function.set_node_attributes(self,regulator_ranks,add2prop)

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

      nx.function.set_node_attributes(self,regulator_ranks,add2prop)
      return regulator_ranks


  def propagate_state(self, seed:PSObject):
      '''
      Input
      -----
      seed must have attribute STATE equal to ACTIVATED or REPRESSED
      '''
      seed_uid = seed.uid()
      #downstream_of_edges = self.out_edges(seed_uid,data='relation')
      nx.set_node_attributes(self,{seed_uid:seed.state()},STATE)
      seed_tree = self.tree4(seed)
      regulator_uids = [seed_uid]
      max_depth = 6
      for level in range (1,max_depth):
          targets_uids = nx.algorithms.traversal.breadth_first_search.descendants_at_distance(seed_tree,seed_uid,level)
          for ruid in regulator_uids:
              regulator = self._get_node(ruid)
              regulator_state = regulator.state()
              if regulator_state:
                  for tuid in targets_uids:
                      rels = self._psrels4(ruid, tuid)
                      target = self._get_node(tuid)
                      for rel in rels:
                          rel_sign = rel.effect_sign()
                          if rel_sign:
                              add2target_state = rel_sign*regulator_state
                              target_state = target.state()
                              target_state += add2target_state
                              nx.set_node_attributes(self,{tuid:target_state},STATE)
          regulator_uids = targets_uids

      annoated_with_state = self.psobjs_with([STATE])
      uid2state = {o.uid():o.state() for o in annoated_with_state}
      return uid2state
          


################## SET-GET SET-GET SET-GET ###############################
  def load_references(self,relpval2weight:dict[str,dict[str,float]]=dict(),weight_name='relweight')->set[Reference]:
    """
    input:
      relpval2weight = {property_name:{prop_value:weight}}, 
      where property_name - PSRelation property used to assign weight to references (usually property_name = OBJECT_TYPE)
    output:
        graph references annotated with 'weight_name' property with values specified in 'relpval2weight'
    """
    my_rels = self._psrels()
    if relpval2weight:
      assert(len(relpval2weight) == 1)
      weight_by_property,using_value2weight = next(iter(relpval2weight.items()))
      [rel.set_weight2ref(weight_by_property,using_value2weight,weight_name) for rel in my_rels]

    graph_references = set()
    [graph_references.update(rel.refs()) for rel in my_rels]
    return graph_references


  def add_node_weight2ref(self,regurn2weight:dict,tarurn2weight:dict,weight_name='nodeweight'):
      '''
      input:
          nodes must have "regulator weight" and "target weight" properties.  
          If both regulator and target node has no "weight" property references with no previously annotated weight will receive zero weight property
      '''
      for r,t,urn,rel in self.edges.data('relation',keys=True):
          regulator_weight = regurn2weight.get(self._get_node(r).urn(),0.0)
          target_weight = tarurn2weight.get(self._get_node(t).urn(),0.0)
          self[r][t][urn]['relation']._set_weight2ref(regulator_weight+target_weight,weight_name)


  def add_weights2neighbors(self,of_nodes:list[PSObject],with_name:str,under_new_name:str,in_direction=''):
      '''
      input:
          in_direction = ['<','>','']
      '''
      assert (in_direction in ['<','>',''])
      urn2weight = {o.urn():o.get_prop(with_name,if_missing_return=0.0) for o in of_nodes}
      urn2weight4annotation = dict()
      for r,t,urn,rel in self.edges.data('relation',keys=True):
          target_weight = regulator_weight = 0.0
          reg_urn = self._get_node(r).urn()
          tar_urn = self._get_node(t).urn()
          if in_direction in ['>','']:
              regulator_weight = urn2weight.get(reg_urn,0)
              if regulator_weight:
                  urn2weight4annotation[tar_urn] = [regulator_weight]
          if in_direction in ['<','']:
              target_weight = urn2weight.get(tar_urn,0)
              if target_weight:
                  urn2weight4annotation[reg_urn] = [target_weight]

      self.set_node_annotation(urn2weight4annotation,under_new_name)


  def citation_index(self)->dict[str,Reference]:
      """
      Returns
      -------
      dictionary {id_type:ref_id:Reference} with all references in graph.\n 
      Reference objects in return dictionary are annotated by PS_CITATION_INDEX property\n
      undirected duplicates are not counted 
      """
      my_graph = self.remove_undirected_duplicates()
      graph_references = dict()
      for r,t,rel in my_graph.edges.data('relation'):
          assert(isinstance(rel,PSRelation))
          for ref in rel.refs():
              id_type, ref_id = ref.get_doc_id()
              if ref_id:
                  try:
                      ref = graph_references[id_type+':'+ref_id]
                      count = int(ref[PS_CITATION_INDEX][0])
                      ref[PS_CITATION_INDEX] = [count + 1]
                      graph_references[id_type+':'+ref_id] = ref
                  except KeyError:
                      ref[PS_CITATION_INDEX] = [1]
                      graph_references[id_type+':'+ref_id] = ref

      sorted_dict = dict(sorted(graph_references.items(),reverse=True,key=lambda item: item[1][PS_CITATION_INDEX][0]))
      return sorted_dict
  

  def author_index(self)->dict[str,int]:
      citation_index = self.citation_index()
      my_refs = set(citation_index.values())
      author_idx = dict()
      for ref in my_refs:
          citation_count = int(ref[PS_CITATION_INDEX][0])
          authors = ref.author_list()
          for au in authors:
              au_clean = au.strip(',.')
              if au:
                  try:
                      author_idx[au_clean] += citation_count
                  except KeyError:
                      author_idx[au_clean] = citation_count
      author_idx = dict(sorted(author_idx.items(), key=lambda item: item[1],reverse=True))
      return author_idx
  

  def number_of_snippets(self):
      return sum([rel.number_of_snippets() for r,t,rel in self.iterate()])


  def textrefs(self):
      all_textrefs = set()
      [all_textrefs.update(rel.textrefs()) for r,t,rel in self.iterate()]
      return all_textrefs
  

  def bibliography(self,df_name:str,for_rel_types=list(),with_effect=list())->df:
      '''
      Return
      ------
      df['Resnet Citation index','Citation','Identifier type','Document identifier: PMID or DOI']\n
      where values of 'Document identifier: PMID or DOI' column are hyperlinked to respective sources
      '''
      sub_graph = self
      if for_rel_types:
          sub_graph = sub_graph.subgraph_by_relprops(for_rel_types)
      if with_effect:
          sub_graph = sub_graph.subgraph_by_relprops(with_effect,['Effect'])
      
      references = sub_graph.citation_index() # annotates
      ref_df = RefStats.external_counter2pd(set(references.values()),stat_prop=PS_CITATION_INDEX)
      ref_df._name_ = df_name

      clinvar_pmids = [['10447503'],['10592272'],['10612825'],['11125122'],['26619011'],
                      ['25741868'],['26582918'],['28492532'],['24033266'],['18414213'],
                      ['26467025'],['24728327']]
      clinvar_hyperlinks = list(map(pubmed_hyperlink,clinvar_pmids))
      ref_df = ref_df.remove_rows_by(clinvar_hyperlinks,IDENTIFIER_COLUMN)
      return ref_df


  def add_recent_refs(self,_2df:df,between_col:str,and_col:str,map2prop='Name'):
      annotated_df = _2df.dfcopy()
      for row in _2df.index:
          node1_prop = annotated_df.loc[row][between_col]
          node2_prop = annotated_df.loc[row][and_col]
          total_refs,ref_ids,_ = self.recent_refs(node1_prop,node2_prop,map2prop,with_children=True)
          dois = str()
          try:
              pmids = ref_ids['PMID']
              refcount = pubmed_hyperlink(pmids,str(total_refs))
          except KeyError:
              try:
                  dois = ref_ids['DOI']
                  refcount = make_hyperlink(dois[0],url='http://dx.doi.org/', display_str=str(total_refs))
              except KeyError:
                  refcount = str(total_refs)
          annotated_df.loc[row][PS_REF_COULUMN] = refcount
      return annotated_df

  
  def recent_refs2df(self, ref_limit=5,pathway=''):
      self.load_references()
      kwargs = {
          'name' : 'RefCount',
          'columns' : ['Regulator','Target',PS_REF_COULUMN,'Number of References','Pathway']
      }

      return_df = df(**kwargs)
      row_counter = 0
      my_graph = self.remove_undirected_duplicates()
      for regulator_id, target_id,rel in my_graph.edges.data():
          regulator = self._psobj(regulator_id)
          target = self._psobj(target_id)
          regulator_name = regulator.name()
          target_name = target.name()
          total_refs, ref_ids,_ = self.__recent_refs([rel['relation']],ref_limit)
          dois = str()
          try:
              pmids = ref_ids['PMID']
              refcount = pubmed_hyperlink(pmids,str(total_refs))
          except KeyError:
              try:
                  dois = ref_ids['DOI']
                  refcount = make_hyperlink(dois[0],url='http://dx.doi.org/', display_str=str(total_refs))
              except KeyError:
                  refcount = str(total_refs)
          return_df.iat[row_counter] = [regulator_name,target_name,refcount,total_refs,pathway]
          row_counter += 1
      return_df = return_df.sortrows(by='Number of References')
      return return_df


  def load_references_between(self,nodes:list,and_nodes:list,weight_by_property='',using_value2weight:dict[str,float]=dict()):
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
      return sub_graph.load_references({weight_by_property:using_value2weight})


  def snippets2df(self, df_name='Snippets',add_nodetype=False, 
                  add_rel_props:list=['Name','RelType',EFFECT],
                  ref_identifiers:list = ['PMID','DOI'],
                  ref_biblio_props:list = [PUBYEAR,TITLE],
                  ref_sentence_props:list=[]):
      '''
      input:
          add_rel_props - can be either empty or ['Name','RelType',EFFECT]
          ref_sentence_props - additional reference properties 
      Returns
      -------
      df with columns: "Concept","Entity",{}'Concept','Concept Type','Entity','Entity Type'},add_rel_props,add_ref_props,PS_CITATION_INDEX,"PMID","DOI",PUBYEAR,TITLE,"Snippets",\n
      where "Concept" contains graph regulators, "Entity" contains graph targets
      '''
      my_graph = self.remove_undirected_duplicates()
      my_graph.load_references()
      if add_nodetype:
          header = ['Concept','Concept Type','Entity','Entity Type']
      else:
          header = ['Concept','Entity']    
  
      header += add_rel_props
      header += [PS_CITATION_INDEX]+ref_sentence_props+ref_identifiers+ref_biblio_props #+['Snippets']
      
      annotated_refs = my_graph.citation_index()
      rows = list()
      for regulatorID, targetID, rel in my_graph.edges.data('relation'):
          regulator_name = my_graph.nodes[regulatorID]['Name'][0]
          target_name = my_graph.nodes[targetID]['Name'][0]
          assert(isinstance(rel,PSRelation))
          reltype = rel.objtype()
          relname = rel.name()
          releffect = rel.effect()
          for ref in rel.refs():
              id_type, identifier = ref.get_doc_id()
              try:
                  annotated_ref = annotated_refs[id_type+':'+identifier]
                  ref_citation_idex = str(annotated_ref[PS_CITATION_INDEX][0])
                  ref_list = ref.to_list(ref_identifiers,False,ref_biblio_props,ref_sentence_props,True)
                  row = [regulator_name,my_graph.nodes[regulatorID][OBJECT_TYPE][0],target_name,my_graph.nodes[targetID][OBJECT_TYPE][0]] if add_nodetype else [regulator_name,target_name]
                  if add_rel_props:
                      row += [relname,reltype,releffect]
                  row += [ref_citation_idex]+ref_list
                  rows.append(row)
              except KeyError:
                  continue

      snippet_df = df.from_rows(rows,header)
      snippet_df[PS_CITATION_INDEX] = snippet_df[PS_CITATION_INDEX].astype(int)
      if PUBYEAR in snippet_df.columns:
          snippet_df[PUBYEAR] = snippet_df[PUBYEAR].astype(int)
      snippet_df = snippet_df.sortrows([PS_CITATION_INDEX,'Concept','Entity'],ascending=[False,True,True])

      clinvar_pmids = [['10447503'],['10592272'],['10612825'],['11125122'],['26619011'],
                      ['25741868'],['26582918'],['28492532'],['24033266'],['18414213'],
                      ['26467025'],['24728327']]
      clinvar_hyperlinks = list(map(pubmed_hyperlink,clinvar_pmids))
      snippet_df = snippet_df.remove_rows_by(clinvar_hyperlinks,'PMID')

      snippet_df._name_ = df_name
      snippet_df.set_hyperlink_color(ref_identifiers)
      return snippet_df

##################### DEL DEL DEL ######################################
  def remove_nodes_by_prop(self, property_values:list[str], prop_names:list[str]=[OBJECT_TYPE]):
      node_ids = self.uids4nodes(property_values,prop_names)
      self.remove_nodes_from(node_ids)
      print("%d nodes with %s were removed" % (len(node_ids), ','.join(property_values)))


  def remove_nodes_by_degree(self,min_degree=0,max_degree=1000000,only_with_prop:list=[OBJECT_TYPE],having_values:list=[]):
      if having_values:
          only_ids = set(self.uids4nodes(having_values,only_with_prop))
          ids2remove = {x for x in self.nodes() if self.degree(x) > max_degree and x in only_ids}
          ids2remove.update({x for x in self.nodes() if self.degree(x) < min_degree and x in only_ids})
      else:
          ids2remove = {x for x in self.nodes() if self.degree(x) > max_degree}
          ids2remove.update({x for x in self.nodes() if self.degree(x) < min_degree})

      self.remove_nodes_from(ids2remove)
      if min_degree:
          print('%d nodes with degree less than %d and more than %d were removed' 
                % (len(ids2remove),min_degree,max_degree))
      else:
          print(f'{len(ids2remove)} nodes with degree > {max_degree} were removed')


  def remove_nodes_by_outdegree(self,min_degree=0,max_degree=1000000,only_with_prop:list=[OBJECT_TYPE],having_values:list=[]):
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

      
  def remove_nodes_by_indegree(self,min_degree=0,max_degree=1000000,only_with_prop:list=[OBJECT_TYPE],having_values:list=[]):
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


  def remove_nodes_by_targets(self,min_target=0,max_target=1000000,only_with_prop:list=[OBJECT_TYPE],having_values:list=[]):
      uids2remove = set()
      for uid in self.nodes():
          targets = list(self.neighbors(uid))
          target_count = len(targets)
          if target_count > max_target:
              if min_target:
                  if target_count < min_target: 
                      uids2remove.add(uid)
              else:
                  uids2remove.add(uid)
      
      if having_values:
          valid_uids = set(self.uids4nodes(having_values,only_with_prop))
          uids2remove = uids2remove.intersection(valid_uids)

      self.remove_nodes_from(uids2remove)
      if min_target:
          print('%d nodes with %d < outdegree < %d were removed' % (len(uids2remove),min_target,max_target))
      else:
          print('%d nodes with outdegree > %d were removed' % (len(uids2remove),max_target))


  def remove_nodes_by_regulators(self,min_regulator=0,max_regulator=1000000,only_with_prop:list=[OBJECT_TYPE],having_values:list=[]):
      uids2remove = set()
      for uid in self.nodes():
          regulators = list(self.predecessors(uid))
          regulator_count = len(regulators)
          if regulator_count > max_regulator:
              if min_regulator:
                  if regulator_count < min_regulator:
                      uids2remove.add(uid)
              else:
                  uids2remove.add(uid)
      
      if having_values:
          valid_uids = set(self.uids4nodes(having_values,only_with_prop))
          uids2remove = uids2remove.intersection(valid_uids)

      self.remove_nodes_from(uids2remove)
      if min_regulator:
          print('%d nodes with %d < outdegree < %d were removed' % (len(uids2remove),min_regulator,max_regulator))
      else:
          print('%d nodes with outdegree > %d were removed' % (len(uids2remove),max_regulator))


  def remove_relation(self, rel:PSRelation):
      rel_urn = rel.urn()
      for pair in rel.get_regulators_targets():
          if self.has_edge(pair[0], pair[1], key=rel_urn):
              self.remove_edge(pair[0], pair[1], key=rel_urn)
      self.urn2rel.pop(rel_urn,'')


  def remove_edges4psobjs(self,with_value:str, in_property='Name'):
      psobjs = self._psobjs_with(with_value,in_property)
      psobjs_uids = self.uids(psobjs)
      edges2remove = [(r,t) for r,t in self.edges() if (r in psobjs_uids) or (t in psobjs_uids)]
      [self.remove_edge(r,t) for r,t in edges2remove]
      print(f'Removed {len(edges2remove)} edges connected to node with {in_property} = {with_value}')


  def remove_edges(self,with_property:str,having_values:list):
      prop2values = {with_property:having_values}
      edges2remove = [(r,t,urn) for r,t,urn,rel in self.edges.data('relation',keys=True) if rel.has_value_in(prop2values)]
      [self.remove_edge(r,t,urn) for r,t,urn in edges2remove]
      print(f'Removed {len(edges2remove)} edges with {having_values} in {with_property} property')


  def remove_references(self,with_prop2values:dict):
      '''
      Input
      -----
      prop2values = {prop_name:[values]}
      ''' 
      norefs_rels = list()
      noduids4norefs_rels = set()
      affected_relcount = 0
      for r,t,rel in self.edges.data('relation'):
          if isinstance(rel,PSRelation):
              noref_rel = rel.remove_references(with_prop2values)
              if noref_rel:
                  norefs_rels.append(noref_rel)
                  noduids4norefs_rels.add(r)
                  noduids4norefs_rels.add(t)
                  if len(noref_rel.refs()) < len(rel.refs()):
                      affected_relcount += 1

      psobjs4norefs_rels = set(self._get_nodes(noduids4norefs_rels))
      norefs_graph = ResnetGraph()
      #norefs_graph.add_psobjs(psobjs4norefs_rels)
      norefs_graph.add_psrels(set(norefs_rels),merge=False)

      deleted_refcount = len(self.load_references()) - len(norefs_graph.load_references())
      print(f'{deleted_refcount} references were deleted in {affected_relcount} relations')

      deleted_relcount = self.number_of_edges()-norefs_graph.number_of_edges()
      print(f'{deleted_relcount} relations without references were removed')
      return norefs_graph


  def clear_resnetgraph(self):
      super().clear()
      self.urn2rel.clear()


######################   GET GET GET   ######################################
  def iterate(self):
    for r,t,rel in self.edges.data('relation'):
      assert(isinstance(rel,PSRelation))
      yield self._get_node(r), self._get_node(t), rel


  def targets_of(self,n:PSObject):
      for r,t,rel in self.edges(n.uid(),data='relation'):
          assert(isinstance(rel,PSRelation))
          yield self._get_node(r), self._get_node(t), rel


  def get_node_attributes(self,prop_name):
      '''
      Returns
      -------
      Dictionary of attributes keyed by node.
      '''
      return nx.function.get_node_attributes(self,prop_name)
  

  def dbids4nodes(self, with_values=None, in_properties:list=[OBJECT_TYPE]):
      '''
      Return
      ------
      if with_values is not supplied returns dbids for nodes from entire graph 
      '''
      if isinstance(with_values,(set,list)):
          all_dbids = set()
          for node in self._get_nodes():
              for propName in in_properties:
                  if PSObject(node).is_annotated(propName, list(with_values)):
                      all_dbids.add(node.dbid())
                      break
          return list(all_dbids)
      else:
          return [n['Id'][0] for i,n in self.nodes(data=True) if 'Id' in n.keys()]
  
  
  def uids4nodes(self, with_values:list, in_properties=[OBJECT_TYPE]):
      all_ids = set()
      for node in self._get_nodes():
          for propName in in_properties:
              if PSObject(node).is_annotated(propName, with_values):
                  all_ids.add(node.uid())
                  break
      return list(all_ids)


  def node_id2urn(self, SearchValues=list(), search_by_properties:list=[OBJECT_TYPE]):
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
      '''
      raises KeyError if PSObject does not exist in self
      '''
      try:
          node_uid = PSObject.urn2uid(urn)
          return self._psobj(node_uid)
      except KeyError:
          raise KeyError
  

  def urn2uid(self, urn:str):
      try:
          my_node = self.urn2node(urn) 
          return my_node.uid()
      except KeyError:
          raise KeyError


  def _get_node(self, with_uid:int):
      return PSObject(self.nodes[with_uid])


  def _get_nodes(self, with_uids:list[int]=[])->list[PSObject]:
      '''
      output:
          [PSObject] for graph nodes with_uids
          if with_uids is empty returns all nodes of self
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

  def _psobjs_with(self,value:str|int,in_property='Name')->list[PSObject]:
      '''
      # finds objects by Name or other property
      Input
      -----
      value - str,int,float
      '''
      my_uids = [uid for uid,values in self.nodes(data=in_property) if value in values]
      return self._get_nodes(list(my_uids)) if my_uids else list()


  def psobj_with_dbids(self,dbids:set)->list[PSObject]:
      my_uids = [uid for uid,ndbids in self.nodes(data=DBID) if ndbids[0] in dbids]
      return self._get_nodes(list(my_uids)) if my_uids else list()


  def psobjs_with(self, with_properties:list=[OBJECT_TYPE], only_with_values:list=[])->list[PSObject]:
      '''
      finds objects by Name,Alias,ObjTypeName or other list of properties\n
      ouput:
          [PSObject], if "only_with_values" is empty returns all nodes annotated "with_properties"
      '''
      return_node_uids = set()
      if only_with_values:
          my_values = set(only_with_values)
          for prop_name in with_properties:
              id2prop = nx.function.get_node_attributes(self,prop_name)
              my_uids = [i for i,v in id2prop.items() if not my_values.isdisjoint(v)]
              return_node_uids.update(my_uids)
      else:
          for prop_name in with_properties:
              id2prop = nx.function.get_node_attributes(self,prop_name)
              return_node_uids.update(id2prop.keys())
      
      #_get_nodes return all nodes if return_node_uids but we need to return empty list here
      return self._get_nodes(list(return_node_uids)) if return_node_uids else list()


  def __node_stats(self):
      my_nodes = self._get_nodes()
      node_stats = dict()
      for n in my_nodes:
          try:
              node_stats[n.objtype()] += 1
          except KeyError:
              node_stats[n.objtype()] = 1

      return dict(sorted(node_stats.items(), key=lambda item: item[1]))


  def __rel_stats(self,_4prop=OBJECT_TYPE):
      '''
      Return
      ------
      dict - rel_stats[rel.get_prop(_4prop)] = (total_count,abs_count,ref1count,ref2count,ref3count)
      '''
      rel_stats = dict()
      my_rels = self._psrels()
      # need to collect unique by URN relations here to not double count non-directional duplications
      for rel in my_rels:
          #reltype = rel.objtype()
          if _4prop == PUBYEAR:
              prop_value = rel._1st_ref().pubyear()
          else:    
              prop_value = rel.get_prop(_4prop)

          if prop_value:
              refcount = rel.count_refs()
              try:
                  total_count,abs_count,ref1,ref2,ref3 = rel_stats[prop_value]
                  total_count += 1          
                  if refcount>1:ref1 += 1
                  if refcount>2:ref2 += 1
                  if refcount>3:ref3 += 1
                  if rel.is_from_abstract():
                      abs_count += 1
                  rel_stats[prop_value] = (total_count,abs_count,ref1,ref2,ref3) 
              except KeyError:
                  total_count = 1
                  ref1 = 1 if refcount>1 else 0
                  ref2 = 1 if refcount>2 else 0
                  ref3 = 1 if refcount>3 else 0
                  abs_count = 1 if rel.is_from_abstract() else 0
                  rel_stats[prop_value] = (total_count,abs_count,ref1,ref2,ref3)

      return dict(sorted(rel_stats.items(), key=lambda item: item[1],reverse=True))


  def timeline(self):
      timeline = self.__rel_stats(PUBYEAR)
      rows = list()
      for year, stats in timeline.items():
          rows.append([year]+list(stats))

      year_col_name = 'Year'
      header = [year_col_name,'Total Counts','Abstract Counts','>1 ref counts','>2 ref counts','>3 ref counts']
      stats_df = df.from_rows(rows,header)
      stats_df[year_col_name] = stats_df[year_col_name].astype(int)
      stats_df = stats_df.sortrows(by=[year_col_name])
      return stats_df


  def get_stats(self):
      '''
      Return
      ------
      {NodeType:Count}, {RelationType:(TotalCount,AbstractCount)}
      '''
      node_stats = self.__node_stats()
      rel_stats = self.__rel_stats()
      
      rows = list()
      for reltype, stats in rel_stats.items():
          rows.append(['Relation',reltype]+list(stats))

      for nodetype, count in node_stats.items():
          rows.append(['Node',nodetype,count,'','','',''])

      header = ['Object','Type','Total Counts','Abstract Counts','>1 ref counts','>2 ref counts','>3 ref counts']
      stats_df = df.from_rows(rows,header)
    #   stats_df.loc[len(stats_df)] = stats_df.sum(axis=0,numeric_only=True)
      return stats_df


  def find_parents4(self,psobjects:list):
      '''
      Input
      -----
      [PSObject] - PSObject list annotated with CHILD property
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


  def _psrels(self):
      """
      Returns
      -------
      {PSRelation} - list of all unique relations in graph
      \n bi-directional duplicates from non-directional relations are removed
      """
      return {PSRelation.copy(r) for n1,n2,r in self.edges.data('relation')}
  

  def relation_dbids(self):
      return [r.dbid() for n1,n2,r in self.edges.data('relation')]


  def _psrels4(self, regulator_uid, target_uid)->list[PSRelation]:
      try:
          edges = dict(self[regulator_uid][target_uid]).values()
          return [e['relation'] for e in edges]
      except KeyError:
          return []


  def psrels_with(self, with_values:list=[], in_properties:list=[OBJECT_TYPE])->set[PSRelation]:
      '''
      Return
      -----
      if "with_values" is empty will return unique set of relations for entire graph with no non-directional duplicates 
      '''
      if with_values:
          relations2return = set()
          for r,t,rel in self.iterate():
              for prop_name in in_properties:
                  if rel.is_annotated(prop_name, with_values):
                      relations2return.add(rel.copy())
                      break
      else:
          relations2return = self._psrels()

      return relations2return

  
  def __find_relations(self, reg_uid, targ_uid, rel_types:list=[], with_effects:list=[], mechanism:list=[], any_direction=False):
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
              my_rels = [x for x in my_rels if x.mechanism() in mechanism]
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
      [positive_refs.update(rel.refs()) for rel in positive_rels]
      negative_refs = set()
      [negative_refs.update(rel.refs()) for rel in negative_rels]

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
          refs = rel.refs()
          if rel.effect() == 'positive':
              positive_counter.update(refs) if refs else positive_counter.add(rel)
          elif rel.effect() == 'negative':
              negative_counter.update(refs) if refs else negative_counter.add(rel)

      if len(positive_counter) > len(negative_counter):
          return 1
      elif len(positive_counter) < len(negative_counter):
          return -1
      return 0


  def net_regulator_effect(self,regulator:PSObject,targets:list,vote_effect_by_ref=False):
      '''
      Input: targets - [PSObject]
      -----
      '''
      activated_targets = list()
      inhibited_targets = list()
      r_uid = regulator.uid()
      if vote_effect_by_ref:
          for target in targets:
              rels = self._psrels4(r_uid,target.uid())
              positive_rels = [r for r in rels if r.effect()=='positive']
              negative_rels = [r for r in rels if r.effect()=='negative']

              positive_refs = set()
              negative_refs = set()
              [positive_refs.update(rel.refs()) for rel in positive_rels]
              [negative_refs.update(rel.refs()) for rel in negative_rels]

              if len(positive_refs) > len(negative_refs):
                  activated_targets.append(target)
              elif len(positive_refs) < len(negative_refs):
                  inhibited_targets.append(target)
      else:
          for target in targets:
              rels = self._psrels4(r_uid,target.uid())
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
                  activated_targets.append(target)
              elif  net_effect < 0:
                  inhibited_targets.append(target)

      return activated_targets,inhibited_targets


  def __targets4(self,regulator:PSObject,linkedby_reltypes:list=[]):
      r_uid = regulator.uid()
      if linkedby_reltypes:
          taregt_uids = [t for r,t,rel in self.edges(r_uid,data='relation') if rel.objtype() in linkedby_reltypes]
      else:
          taregt_uids = [t for r,t in self.edges(r_uid)]
      return self._get_nodes(taregt_uids)


  def __mean_effect(self,of_regulator:PSObject|int,with_reltypes:set=DIRECT_RELTYPES,cutoff = 0.8):
      '''
      Return
      ------
      mean effect of the regulator on its targets
      used to predict if drug is inhibitor or agonist.  
      cutoff - the minimal proportion of targets that must have the same effect sign
      '''
      if isinstance(of_regulator,int):
          of_regulator = self._get_node(of_regulator)
      targets = self.__targets4(of_regulator,list(with_reltypes))

      activated_targets,inhibited_targets = self.net_regulator_effect(of_regulator,targets,True)
      targets_count = len(activated_targets)+len(inhibited_targets)
      if len(activated_targets) > cutoff*targets_count:
          return 'positive' 
      elif len(inhibited_targets) > cutoff*targets_count:
          return 'negative'
      else:
          return 'unknown'


  def predict_effect4(self,_4enttypes:list,_4reltypes:list):
      '''
      Return
      ------
      predicts effect for all targets of a regulator based on the ResnetGraph.__mean_effect(regulator)
      graph copy with predicted Effect _4enttypes with _4reltypes\n
      [PSRelation] that were modified
      '''
      my_copy = self.copy()
      counter = 0
      lazy_dict = dict()
      predicted_rels = list()
      for r,t,urn,rel in my_copy.edges(data='relation',keys=True):
          if rel.objtype() in _4reltypes and rel.effect() == 'unknown':
              regulator = my_copy._get_node(r)
              if regulator.objtype() in _4enttypes:
                  try:
                      predicted_effect = lazy_dict[r]
                  except KeyError:
                      predicted_effect = my_copy.__mean_effect(regulator,cutoff=0.8)
                      lazy_dict[r] = predicted_effect
                      
                  if predicted_effect != 'unknown':
                      my_copy[r][t][urn]['relation'][EFFECT] = [predicted_effect]
                      counter += 1
                      predicted_rels.append(my_copy[r][t][urn]['relation'])
      print(f'Assigned Effect to {counter} relations in {my_copy.name} with {my_copy.number_of_edges()} edges')
      return my_copy,predicted_rels


  def predict_effect(self,mean_effect_cutoff=0.8):
      '''
      Return
      ------
      predicts effect for all targets of a regulator based on the ResnetGraph.__mean_effect(regulator)
      graph copy with predicted Effect _4enttypes with _4reltypes\n
      [PSRelation] that were modified
      '''
      my_copy = self.copy()
      counter = 0
      lazy_dict = dict()
      predicted_rels = list()
      for r,t,urn,rel in my_copy.edges(data='relation',keys=True):
          if rel.effect() == 'unknown':
              regulator = my_copy._get_node(r)
              try:
                  predicted_effect = lazy_dict[r]
              except KeyError:
                  predicted_effect = my_copy.__mean_effect(regulator,cutoff=mean_effect_cutoff)
                  lazy_dict[r] = predicted_effect
                  
              if predicted_effect != 'unknown':
                  my_copy[r][t][urn]['relation'][EFFECT] = [predicted_effect]
                  counter += 1
                  predicted_rels.append(my_copy[r][t][urn]['relation'])
      print(f'Assigned Effect to {counter} relations in {my_copy.name} with {my_copy.number_of_edges()} edges')
      return my_copy,predicted_rels


  def get_prop2obj_dic(self, search_by_property:str, filter_by_values=[], case_insensitive=False)->tuple[dict[str,list[PSObject]],dict[int,list[str]]]:
      '''
      Returns
      -------
      propval2objs = {search_by_property_value:[PSObject]}
      objuid2propval = {id:[search_by_property_value]}\n
      if filter_by_values empty returns dictionary keyed by all values in search_by_property
      '''
      search_value2obj = defaultdict(list)
      objid2search_values = defaultdict(list)
      if filter_by_values:
          allowed_values = set(map(lambda x:str(x).lower(),filter_by_values)) if case_insensitive else filter_by_values
          for n_uid, n in self.nodes(data=True):
              try:
                  node_prop_values = list(map(lambda x:str(x).lower(),n[search_by_property])) if case_insensitive else n[search_by_property]
                  matched_values = list(set(node_prop_values).intersection(allowed_values))
                  if matched_values:
                      objid2search_values[n_uid] = matched_values
                      [search_value2obj[v].append(PSObject(n)) for v in matched_values]
              except KeyError: continue
      else:
          for n_uid, n in self.nodes(data=True):
              try:
                  all_values = list(map(lambda x: str(x).lower(),n[search_by_property])) if case_insensitive else n[search_by_property] 
                  objid2search_values[n_uid] = all_values
                  [search_value2obj[v].append(PSObject(n)) for v in all_values]
              except KeyError: continue
      return search_value2obj, objid2search_values


  def props2obj_dict(self, propValues:list, prop_names:list, case_insensitive=False) -> tuple[dict[str, list[PSObject]], dict[int, list[str]]]:
      '''
      Returns
      -------
      propval2objs = {search_by_property_value:[PSObject]}
      objuid2propval = {id:[search_by_property_value]}\n
      if propValues empty returns dictionary keyed by all values in prop_names
      '''
      propval2objs = dict()
      uid2propval = dict()
      for prop_name in prop_names:
          p2o, i2p = self.get_prop2obj_dic(prop_name, propValues,case_insensitive)
          for p,objs in p2o.items():
              try:
                  mapped_objs = set(propval2objs[p])
                  mapped_objs.update(objs)
                  propval2objs[p] = list(mapped_objs)
              except KeyError:
                  propval2objs[p] = objs

          for uid,prop_vals in i2p.items():
              try:
                  mapped_values = set(uid2propval[uid])
                  mapped_values.update(prop_vals)
                  uid2propval[uid] = list(mapped_values)
              except KeyError:
                  uid2propval[uid] = prop_vals

      return propval2objs, uid2propval


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
      
      all_uids = self.uids(list(all_entities))
      all_rels = [rel for regulatorID, targetID, rel in self.edges.data('relation') 
              if regulatorID in all_uids and targetID in all_uids]
      
      if filter_rel_with:
          return [rel for rel in all_rels if rel.has_value_in(filter_rel_with)]
      else:
          return all_rels


  @staticmethod
  def __recent_refs(relations:list[PSRelation],ref_limit=5)->tuple[int,dict[str,list],list[Reference]]:
      '''
      Input
      -----
      [PSRelations] - list of PSRelation objects
      '''
      references = set()
      [references.update(r.refs()) for r in relations]
      references = list(references)

      def sortkey(ref:Reference):
          try:
              year = int(ref[PUBYEAR][0])
          except KeyError:
              try:
                  year = ref['Start'][0]
                  year = int(year[-4:])
              except KeyError: year = 1812
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
      return {rel.dbid() for _,_,rel in self.edges.data('relation')}


  def find_nodes(self,for_relation:PSRelation,filter_by:list=[],in_properties:list=[])->tuple[list[PSObject],list[PSObject]]:
      """
      Returns
      -------
      regulators = [PSObject], targets = [PSObject]
      """
      regulators = for_relation.regulators()
      targets = for_relation.targets()

      if filter_by:
          search_in_properties = in_properties if in_properties else [OBJECT_TYPE]
          must_values = set(filter_by)
          filtered_regulators = set()
          filtered_targets = set()
          for prop_name in search_in_properties:
              filtered_regulators.update([x for x in regulators if not set(x[prop_name]).isdisjoint(must_values)]) 
              filtered_targets.update([x for x in targets if not set(x[prop_name]).isdisjoint(must_values)]) 
          return list(filtered_regulators), list(filtered_targets)
      else:
          return regulators, targets

  
  def get_properties(self,property_name:str, for_ids=set()):
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
  

  def __has_nodetype(self, nodetype:str):
      for uid in self.nodes():
          psobj = self._get_node(uid)
          if psobj.objtype() == nodetype:
              return True
      return False
  

  def has_nodetypes(self,nodetypes:list):
      for type in nodetypes:
          if not self.__has_nodetype(type):
              return False
      return True


  def __has_reltype(self, reltype:str):
      for r,t,psrel in self.edges.data('relation'):
          if psrel.objtype() == reltype:
              return True
      return False
  

  def has_reltypes(self,reltypes:list):
      for type in reltypes:
          if not self.__has_reltype(type):
              return False
      return True


  def has_node_with_id(self, node_id:int):
      return super().has_node(node_id)


  def _psobj(self, node_uid:int):
      '''
      raises KeyError if node_uid does not exist in self
      '''
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
  def uids(nodes:list[PSObject])->list[int]:
      '''
      Input
      -----
      nodes = [PSObject]
      '''
      return [n.uid() for n in nodes]
  

  @staticmethod
  def urns(nodes:list[PSObject]):
      '''
      Input
      -----
      nodes = [PSObject]
      '''
      return [n.urn() for n in nodes]
  

  @staticmethod
  def dbids(nodes:list[PSObject])->list[int]:
      '''
      Input
      -----
      nodes = [PSObject]
      '''
      return [n.dbid() for n in nodes]


  @staticmethod
  def names(nodes:list[PSObject]):
      '''
      Input
      -----
      nodes = [PSObject]
      '''
      return [n.name() for n in nodes]
  

  @staticmethod
  def classes(nodes:list[PSObject]):
      '''
      Input
      -----
      nodes = [PSObject]
      '''
      return {n.get_prop('Class',0,'') for n in nodes}
  

  @staticmethod
  def find_object(my_set:set[PSObject], urn:str):
      """
      Finds the first object in the set that matches the given condition.
      Returns a reference to the object, or None if no match is found.
      """
      return next((o for o in my_set if o.urn() == urn), PSObject())


  @staticmethod
  def childs(nodes:list[PSObject]):
      childs = set()  
      [childs.update(n.childs()) for n in nodes]
      return list(map(PSObject,childs))
  

  def dbid2uid(self,dbids:list[int]=[])->dict[int,int]:
      '''
      Input
      -----
      nodes = [PSObject]

      Returns
      -------
      {node.dbid():node.uid()}
      '''
      if dbids:
          nodes_with_dbids = self.psobj_with_dbids(set(dbids))
          return {n.dbid():n.uid() for n in nodes_with_dbids}
      else:
          all_nodes = self._get_nodes()
          return {n.dbid():n.uid() for n in all_nodes}
  

  def uid2dbid(self,for_uids:list[int]=[])->dict[int,int]:
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
      return self.urns(list(self.__psobjs()))
  

  def rel_urns(self):
      return self.urns(list(self._psrels()))
  

  def node_props(self,_4props:list[str], in_psobjs:list[PSObject]=[])->set[int|str]:
      psobjs = in_psobjs if in_psobjs else self.psobjs_with(_4props)
      all_props = set()
      for p in _4props:
          for o in psobjs:
              all_props.update(o[p])
      return all_props


  def get_neighbors(self,of_nodes:set[PSObject],allowed_neigbors:list[PSObject]=[]):
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
              n_neighbors_uids = list(nx.function.all_neighbors(self, n.uid()))
              neighbor_uids.update(n_neighbors_uids)

      if allowed_neigbors:
          allowed_neigbors_uids = self.uids(allowed_neigbors)
          neighbor_uids = neighbor_uids.intersection(allowed_neigbors_uids)
      
      return self.__psobjs(neighbor_uids)


  def get_neighbors_rels(self, _4objs:set[PSObject], only_neighbors:list[PSObject]=[]):
      '''
      Input
      -----
      _4psobs - {PSObject}
      only_neighbors - [PSObject]
      '''
      neighbor_graph = self.neighborhood(_4objs,only_neighbors)
      return neighbor_graph._psrels()


  def get_neighbors_refs4(self, _4psobs:set[PSObject], only_neighbors:list[PSObject]=[]):
      '''
      Input
      -----
      _4psobs - {PSObject}
      only_neighbors - [PSObject]
      '''
      neighbor_graph = self.neighborhood(_4psobs,only_neighbors)
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
              return [v['relation'] for v in edges.values() if v['relation'][OBJECT_TYPE][0] in from_relation_types]
          else:
              return [r['relation'] for r in edges.values()]
      except KeyError:
          return list()


  def regulators(self, only_objtype=[], min_targets=1):
      if only_objtype:
          return {PSObject(y) for x,y in self.nodes(data=True) if ((self.out_degree(x) >= min_targets) & (y[OBJECT_TYPE][0] in only_objtype))}
      else:
          return {PSObject(y) for x,y in self.nodes(data=True) if self.out_degree(x) >= min_targets}


  def unconnected_node_ids(self):
      return {PSObject(x) for x in self.nodes() if self.degree(x) == 0}

  
  def root_nodes(self)->list[int]:
      return [i for i in self.nodes if self.in_degree(i)==0]

  
  def prop2outdegree(self,prop_name='Name'):
      return {p[0]:self.out_degree(i) for i,p in self.nodes(data=prop_name)}


################################# WRITE-DUMP, WRITE-DUMP, WRITE-DUMP ##############################
  def print_triples(self, fileOut, relPropNames, access_mode='w', printHeader=True, add_entities=False, as1row=False):
      with open(fileOut, access_mode, encoding='utf-8') as f:
          if printHeader:
              header = '\t'.join(relPropNames) + '\t' + "Regulators Id" + '\t' + "Targets Id"
              f.write(header + '\n')

          for _,_, rel in self.edges.data('relation'):
              assert(isinstance(rel,PSRelation))
              f.write(rel.triple2str(relPropNames,add_entities=add_entities, as1row=as1row))


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
                      assert(isinstance(rel,PSRelation))
                      reference_view_triple = str(rel.triple2str(relPropNames,as1row=single_rel_row))
                      f.write(reference_view_triple)
              else:
                  for regulatorID, targetID, rel in self.edges.data('relation'):
                      assert(isinstance(rel,PSRelation))
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


  def ref2pandas (self, relPropNames:list, entity_prop_names=[], RefNumPrintLimit=0,single_rel_row=False) -> df:
      temp_fname = '__temp__.tsv'
      self.print_references(temp_fname,relPropNames,entity_prop_names,RefNumPrintLimit=RefNumPrintLimit,single_rel_row=single_rel_row)
      to_return = df.read(temp_fname,header=0,index_col=False, dtype='unicode')
      os.remove(temp_fname)
      return to_return


  def __2resnet(self,resnet:et._Element,ent_props:list,rel_props:list,
                add_rel_props:dict[str,list[str]]={},add_pathway_props:dict[str,list[str]]={},
                delete_nodes=False):
      '''
      Input
      -----
      add_rel_props,add_pathway_props structure - {PropName:[PropValues]}
      '''
      def _2b_printed(prop_name:str,prop_list:list):
          return prop_name in prop_list if prop_list else True

      if add_pathway_props:
          pathway_props = et.SubElement(resnet, 'properties',attrib=None, nsmap=None)
          for prop_name,prop_val in add_pathway_props.items():
              for val in prop_val:
                  et.SubElement(pathway_props, 'attr', {'name':str(prop_name), 'value':str(val)},nsmap=None)
      
      node_attr = {'delete':'true'} if delete_nodes else None
      xml_nodes = et.SubElement(resnet,'nodes',attrib=node_attr,nsmap=None)  
      for nodeId, n in self.nodes(data=True):
          try: 
              local_id = n['URN'][0]
              xml_node = et.SubElement(xml_nodes, 'node', {'local_id': local_id, 'urn': n['URN'][0]},nsmap=None)
              et.SubElement(xml_node, 'attr', {'name': 'NodeType', 'value': str(n[OBJECT_TYPE][0])},nsmap=None)
              for prop_name, prop_values in n.items():
                  if _2b_printed(prop_name, ent_props):
                      for prop_value in prop_values:
                          et.SubElement(xml_node, 'attr', {'name': str(prop_name), 'value': str(prop_value)},nsmap=None)
          except KeyError:
              continue

      xml_controls = et.SubElement(resnet, 'controls',attrib=None,nsmap=None)
      graph_relations = self._psrels()
      for rel in graph_relations:
          control_id = rel.urn()
          xml_control = et.SubElement(xml_controls, 'control', {'local_id':control_id},nsmap=None)
          et.SubElement(xml_control, 'attr', {'name':'ControlType', 'value':str(rel[OBJECT_TYPE][0])},nsmap=None)
          
          # adding links
          if TARGETS in rel.Nodes:
              linktype4reg = 'in'
              for t in rel.Nodes[TARGETS]:
                  et.SubElement(xml_control, 'link', {'type': 'out', 'ref': t.urn()},nsmap=None)
          else:
              linktype4reg = 'in-out'

          for r in rel.Nodes[REGULATORS]:
              try:
                  et.SubElement(xml_control, 'link', {'type':linktype4reg, 'ref':r.urn()},nsmap=None)
              except IndexError or KeyError:
                  continue
          # non-reference properties
          for prop_name, prop_values in rel.items():
              if _2b_printed(prop_name,rel_props):
                  for prop_value in prop_values:
                      et.SubElement(xml_control, 'attr', {'name':str(prop_name), 'value':str(prop_value)},nsmap=None)

          for prop_name,prop_val in add_rel_props.items():
              for val in prop_val:
                  et.SubElement(xml_control, 'attr', {'name':str(prop_name), 'value':str(val)},nsmap=None)

          # adding references
          snippet_props = set(REFERENCE_PROPS).intersection(rel_props)
          print_snippets = True if not rel_props else True if snippet_props else False
          if print_snippets:
            references = list(set(rel.refs()))
            ref_index = 0
            for ref in references:
              # each snippet has its own index in RNEF
              for textref, snippet in ref._snippets():
                et.SubElement(xml_control, 'attr',{'name': str('TextRef'), 'value': textref, 'index': str(ref_index)},nsmap=None)
                for sentprop_name, sentprop_values in snippet.items():
                  if _2b_printed(sentprop_name,list(snippet_props)):
                    v_str = ','.join(sentprop_values)
                    et.SubElement(xml_control, 'attr',{'name':str(sentprop_name), 'value':str(v_str), 'index':str(ref_index)},nsmap=None)

                for prop_name, prop_values in ref.items():
                    if _2b_printed(prop_name,list(snippet_props)):
                        for prop_value in prop_values:
                            et.SubElement( xml_control, 'attr',{'name': str(prop_name), 'value': str(prop_value), 'index': str(ref_index)},nsmap=None)
                    
                for ref_id_type,ref_id in ref.Identifiers.items():
                    if _2b_printed(ref_id_type,list(snippet_props)):
                        et.SubElement(xml_control, 'attr',{'name':str(ref_id_type), 'value':str(ref_id), 'index':str(ref_index)},nsmap=None)
                ref_index+=1

              if not ref.snippets:
                  textref = ref._make_textref()
                  et.SubElement(xml_control, 'attr',{'name': str('TextRef'), 'value': textref, 'index': str(ref_index)},nsmap=None)
                  for prop_name, prop_values in ref.items():
                      if prop_name in snippet_props:
                          for prop_value in prop_values:
                              et.SubElement( xml_control, 'attr',{'name': str(prop_name), 'value': str(prop_value), 'index': str(ref_index)},nsmap=None)
                      
                  for ref_id_type,ref_id in ref.Identifiers.items():
                      if ref_id_type in snippet_props:
                          et.SubElement(xml_control, 'attr',{'name':str(ref_id_type), 'value':str(ref_id), 'index':str(ref_index)},nsmap=None)
                  
                    
  def to_rnefstr(self,ent_props:list,rel_props:list,add_rel_props:dict={},add_pathway_props:dict={},delete_nodes=False):
      resnet_attr = {'refonly':'true'} if delete_nodes else None
      resnet = et.Element('resnet',resnet_attr,nsmap=None)
      self.__2resnet(resnet,ent_props,rel_props,add_rel_props,add_pathway_props,delete_nodes)
      xml_str = et.tostring(resnet)
      return xml_str


  def __2rnef(self,to_file:str,ent_props:list,rel_props:list,add_rel_props:dict={},add_pathway_props:dict={},delete_nodes=False):
      with et.xmlfile(to_file,encoding='utf-8',buffered=False) as xf:
          xf.write(et.Comment(RNEF_DISCLAIMER),pretty_print=True)
          resnet_attr = {'refonly':'true'} if delete_nodes else None
          with xf.element('batch'):
              resnet = et.Element('resnet',attrib=resnet_attr,nsmap=None)
              self.__2resnet(resnet,ent_props,rel_props,add_rel_props,add_pathway_props,delete_nodes)
              xf.write(resnet,pretty_print=True)
  

  def __2rnef_secs(self,xmlfile:et.xmlfile,ent_prop2print:list,rel_prop2print:list,add_rel_props=dict(),with_section_size=1000,delete_nodes=False):
      """
      splits RNEF into <resnet> sections with_section_size <control> elements
      used to write large graphs to file by redcuing the length of xml string
      \nwith_section_size = 1000 recommended to avoid memory problems
      \nresolves printing graphs with and without edges\n
      """
      resnet_attr = {'refonly':'true'} if delete_nodes else None
      if self.number_of_edges():
          resnet_sections_rels = set()
          for regulatorID, targetID, rel in self.edges.data('relation'):
              resnet_sections_rels.add(rel)
              if len(resnet_sections_rels) == with_section_size:
                  section_graph = self.subgraph_by_rels(list(resnet_sections_rels))
                  resnet = et.Element('resnet',attrib=resnet_attr,nsmap=None)
                  section_graph.__2resnet(resnet,ent_prop2print,rel_prop2print,add_rel_props,delete_nodes=delete_nodes)
                  xmlfile.write(resnet,pretty_print=True)
                  resnet_sections_rels.clear()
                  
          # printing leftovers
          section_graph = self.subgraph_by_rels(list(resnet_sections_rels))
          resnet = et.Element('resnet',attrib=resnet_attr,nsmap=None)
          section_graph.__2resnet(resnet,ent_prop2print,rel_prop2print,add_rel_props,delete_nodes=delete_nodes)
          xmlfile.write(resnet,pretty_print=True)
          return 
      else:
          all_nodes = self._get_nodes()
          for sec in range(0, len(all_nodes), with_section_size):
              section_nodes = all_nodes[sec:sec+with_section_size]
              section_graph = ResnetGraph()
              section_graph.add_psobjs(set(section_nodes))
              rnef_str = section_graph.to_rnefstr(ent_prop2print,rel_prop2print,add_rel_props,delete_nodes=delete_nodes)
              rnef_str = str(minidom.parseString(rnef_str).toprettyxml(indent='  '))
              rnef_str = rnef_str[rnef_str.find('\n')+1:]
              xmlfile.write(rnef_str)
              return


  def remove_undirected_duplicates(self):
      relset = set()
      to_return = self.copy()
      for r,t,rel in self.edges.data('relation'):
          assert(isinstance(rel,PSRelation))
          if rel in relset:
              to_return.remove_edge(r,t,rel.urn())
          else:
              relset.add(rel)

      return to_return
      

  def dump2rnef(self,fname:str,ent_prop2print:list=['Name'],rel_prop2print:list=[],add_rel_props:dict={},with_section_size=0,delete_nodes=False):
      '''
      Input
      -----
      if fname is empty will create file with graph self.name

      Dumps
      -----
      graph into RNEF file with <resnet> sections of size "with_section_size".
      if "with_section_size" is zero dumps graph into RNEF file with single <resnet>. Large graphs must be dumped with resnet sections
      single <resnet> section reduces size of RNEF file, but may slow down import of large RNEF files into Pathway Studio
      '''
      rnef_fname = fname if fname else self.name
      if fname[-5:] != '.rnef':
          rnef_fname += '.rnef'

      message = f'Writing graph "{self.name}" with {self.number_of_nodes()} nodes, {self.number_of_edges()} edges to {rnef_fname} file'
      graph_copy = self.remove_undirected_duplicates() 
      # copying graph to enable using the function in multithreaded file writing
      
      if with_section_size:
          with et.xmlfile(rnef_fname,encoding='utf-8',buffered=False) as xf:  
              print(message + f' in resnet section of size {with_section_size}')
              xf.write(et.Comment(RNEF_DISCLAIMER),pretty_print=True) 
              with xf.element('batch'):
                  graph_copy.__2rnef_secs(xf,ent_prop2print,rel_prop2print,add_rel_props,with_section_size,delete_nodes)
      else:
          print(message + f' in one resnet section')
          graph_copy.__2rnef(rnef_fname,ent_prop2print,rel_prop2print,add_rel_props,delete_nodes=delete_nodes)


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
  
      rows2add = list()
      for rel in relations:
          row = row_template
          for i in range(0, len(from_properties)):
              if from_properties[i][:2] == 'R:':
                  cell = rel._props2str(from_properties[i][2:])
                  row[i] = cell
                  rows2add.append(row)

      if not rows2add: 
          rows2add = [row_template]
      
      new_data = df.from_rows(rows2add, header=to_df.columns.to_list())
      to_df = to_df.append_df(new_data)
      return to_df

################################# READ READ READ ##########################################
  @staticmethod
  def _parse_nodes_controls(resnet:et._Element, prop2values:dict=dict(),
      only_relprops:set=set(),only4objs:set[PSObject]=set(),on_both_ends=True)->tuple[set[PSObject],set[PSRelation]]:
      '''
      prop2values = {prop_name:[values]} - filter to load relations only for nodes with desired properties,\n
      merge - should be false if graph is loaded from RNEF file with single <resnet> section
      '''
      def validate(my_regulators:list[PSObject], my_targets:list[PSObject]):
          is_valid = False
          if only4objs:
              is_valid = set(my_regulators).issubset(only4objs)
              if not is_valid: return False
              if on_both_ends:
                  is_valid = set(my_targets).issubset(only4objs)
              if not is_valid: return False

          if prop2values:
              if any(n.has_value_in(prop2values) for n in my_regulators + my_targets):
                  return True
          
          return True
  
      nodel_local_ids = dict()
      for node in resnet.findall('./nodes/node',''):
          node_urn = node.get('urn')
          local_id = node.get('local_id')
          node_psobj = PSObject({'URN':[node_urn]})
          [node_psobj.update_with_value(attr.get('name'), attr.get('value')) for attr in node.findall('attr')]
          node_psobj[OBJECT_TYPE] = node_psobj.pop('NodeType')

          exist_obj = ResnetGraph.find_object(only4objs,node_urn)
          if exist_obj:
              node_psobj = exist_obj = node_psobj.merge_obj(exist_obj)

          nodel_local_ids[local_id] = node_psobj
  
      valid_nodes = set()
      valid_rels = list()
      if only_relprops: only_relprops.add('ControlType')
      
      for rel_tag in resnet.findall('./controls/control',''):
          regulators = list()
          targets = list()
          psobjs4rel = list()
          for link in rel_tag.findall('link'):
              link_ref = link.get('ref')
              link_psobj = nodel_local_ids[link_ref]
              assert(isinstance(link_psobj,PSObject))
              psobjs4rel.append(link_psobj)
              link_type = link.get('type')
              if link_type == 'out': 
                  targets.append(link_psobj)
              else: 
                  regulators.append(link_psobj)

          is_valid_rel = validate(regulators,targets)
          if is_valid_rel:
              valid_nodes.update(psobjs4rel)
              ps_rel = PSRelation(dict())
              [ps_rel.Nodes[REGULATORS].append(reg) for reg in regulators]
              [ps_rel.Nodes[TARGETS].append(targ) for targ in targets]

              for attr in rel_tag.findall('attr'):
                  prop_id = attr.get('name')
                  if only_relprops and prop_id not in only_relprops: continue
                  prop_value = attr.get('value')    
                  index = attr.get('index')
                  if type(index) == type(None):
                      # case if RNEF is generated by MedScan
                      if prop_id in SENTENCE_PROPS_SET: # MedScan does not generate "index" if relation has only one reference
                          propid = SENTENCE if prop_id == 'msrc' else prop_id
                          if '1' in ps_rel.PropSetToProps:
                              ps_rel.PropSetToProps['1'][propid] = [prop_value]
                          else:
                              ps_rel.PropSetToProps['1'] = {propid:[prop_value]}
                      else:
                          ps_rel.update_with_value(prop_id, prop_value)
                  else:
                      # property has index
                      propid = SENTENCE if prop_id == 'msrc' else prop_id
                      if index in ps_rel.PropSetToProps:
                          try:
                              ps_rel.PropSetToProps[index][propid].append(prop_value)
                          except KeyError:
                              ps_rel.PropSetToProps[index][propid] = [prop_value]
                      else:
                          ps_rel.PropSetToProps[index] = {propid:[prop_value]}

              ps_rel[OBJECT_TYPE] = ps_rel.pop('ControlType')
              ps_rel.refs()
              ps_rel.urn()
              valid_rels.append(ps_rel)
      #ps_rel does not have 'Id' property.This is used by PSRelation.is_from_rnef()
      return valid_nodes, set(valid_rels)


  @staticmethod
  def __read_rnef(rnef_file:str,prop2values:dict=dict(),only_relprops:set=set(),
                  no_mess=False,only4objs:set[PSObject]=set(),on_both_ends=True)->tuple[set[PSObject],set[PSRelation]]:
      '''
      Input
      -----
      prop2values={prop_name:[values]} - filter to load relations only for nodes with desired properties
      '''
      nodes = set()
      rels = set()
      if not no_mess:
          print ('\nLoading graph from file %s' % rnef_file,flush=True)
      with open(rnef_file, "rb") as f:
          context = et.iterparse(f, tag="resnet")
          for action, elem in context:
              resnet_nodes,resnet_rels = ResnetGraph._parse_nodes_controls(elem,prop2values,only_relprops,only4objs,on_both_ends)
              nodes.update(resnet_nodes)
              rels.update(resnet_rels)
              elem.clear()
          del context
      return nodes,rels


  @classmethod
  def fromRNEF(cls,rnef_file:str,
                prop2values:dict=dict(),only_relprops:set=set(),merge=False,no_mess=False,
                only4objs:set[PSObject]=set(),on_both_ends=True):
      '''
      Input
      -----
      set merge=True if graph loaded from multiple RNEF files with multiple <resnet> sections
      prop2values={prop_name:[values]} - filter to load relations only for nodes with desired properties

      Raises FileNotFoundError if "rnef_file" is not found
      '''
      try:
          start = time.time()
          g = ResnetGraph()
          nodes,rels = g.__read_rnef(rnef_file,prop2values,only_relprops,no_mess,only4objs,on_both_ends)
          g.name = f'from {rnef_file}'
          #g.add_psobjs(nodes,merge)
          g.add_psrels(rels,merge)

          if not no_mess:
              print('File %s with %d edges and %d nodes was loaded in %s' 
              % (rnef_file,g.number_of_edges(),g.number_of_nodes(),execution_time(start)))
          return g
      except FileNotFoundError:
          raise FileNotFoundError


  @classmethod
  def fromRNEFflist(cls,flist:list[str],prop2values:dict=dict(),
                  only_relprops:set=set(),merge=True):
      '''
      Input
      -----
      set merge=True if graph loaded form multiple RNEF files with multiple <resnet> sections
      prop2values={prop_name:[values]} - filter to load relations only for nodes with desired properties
      '''
      max_workers = min(32,len(flist),(os.cpu_count() or 1) + 4)
      combo_g = ResnetGraph()
      with ThreadPoolExecutor(max_workers=max_workers,thread_name_prefix='readRNEFdir') as e:
          futures = list()
          [futures.append(e.submit(combo_g.__read_rnef,f,prop2values,only_relprops,True)) for f in flist]
          for f in as_completed(futures):
              nodes,rels = f.result()
              #combo_g.add_psobjs(nodes,merge)
              combo_g.add_psrels(rels,merge)
          e.shutdown()    
      return combo_g


  @classmethod
  def readRNEFflist_in_batches(cls,flist:list[str],batch_size=0):
    batch_size = batch_size if batch_size else min(32,len(flist),(os.cpu_count() or 1) + 4)
    for i in range(0, len(flist), batch_size):
      batch_end = min(i + batch_size, len(flist))
      batch = flist[i:batch_end]
      # Process the batch of files and yield data
      yield cls.fromRNEFflist(batch)


  @classmethod
  def fromRNEFdir(cls,path2dir:str,prop2values:dict=dict(),
                  only_relprops:set=set(),merge=True,include_subdirs=False):
      '''
      Input
      -----
      set merge=True if graph loaded form multiple RNEF files with multiple <resnet> sections
      prop2values={prop_name:[values]} - filter to load relations only for nodes with desired properties

      Raises FileNotFoundError if "rnef_file" is not found
      '''
      start = time.time()
      real_path = os.path.join(path2dir, '')
      listing = glob.glob(f'{real_path}/*.rnef')+glob.glob(os.path.join(real_path, '**/*.rnef'),recursive=include_subdirs)
      combo_g = ResnetGraph()
      if listing:
          combo_g = ResnetGraph.fromRNEFflist(listing,prop2values,only_relprops,merge)
          print('Graph (%d edges, %d nodes) was loaded from "%s" with %d files in %s' 
          % (combo_g.number_of_edges(),combo_g.number_of_nodes(),real_path,len(listing),execution_time(start)))
          combo_g.name = f'{os.path.basename(os.path.normpath(real_path))}'
      else:
          print('Cannot find "%s" directory' % real_path)
      
      return combo_g
      

  @classmethod
  def read_dir_in_batches(cls,path2dir:str,include_subdirs=True,subdirs_only=False):
      def soft_pop(l:list,index:int):
          try:
              return l.pop(index)
          except IndexError:
              return ('','')

      def redistribute_files_by_size(dir_flist:list[str],batch_size:int):
          flist = [(f,os.path.getsize(f)) for f in dir_flist]
          flist = sorted(flist, key=lambda x: x[1])
          number_of_batches = math.ceil(len(flist)/batch_size)
          batches = list()
          [batches.append(list()) for _ in range(number_of_batches)] # [[]]*number_of_batches does not work here
          while flist:
              for i in range(0,number_of_batches):
                  batches[i] += [soft_pop(flist,i),soft_pop(flist,-(i+1))]
              for i in range(0,number_of_batches):
                  batches[i] += [soft_pop(flist,i),soft_pop(flist,-(i+1))]

          flist = [i for sub_list in batches for i in sub_list if i != ('','') ]
          flist = [x[0] for x in flist]
          return flist

      def get_flist(path2dir,include_subdirs,subdirs_only):
          real_path = os.path.join(path2dir, '')
          rnef_files = []
          for root, dirs, files in os.walk(real_path):
              if not include_subdirs and root != real_path:
                  continue
              if subdirs_only and root == real_path:
                  continue
              
              for filename in files:
                  if filename.lower().endswith(".rnef"):
                      full_path = os.path.join(root, filename)
                      rnef_files.append(full_path)

          return rnef_files

      flist = get_flist(path2dir,include_subdirs,subdirs_only)
      batch_size = min(32,len(flist),(os.cpu_count() or 1) + 4)
      flist = redistribute_files_by_size(flist,batch_size)
      process_start = time.time()
      number_of_iteration = math.ceil((len(flist))/batch_size)
      for i in range(0, len(flist), batch_size):
          batch_end = min(i + batch_size, len(flist))
          batch = flist[i:batch_end]
          # Process the batch of files and yield data
          batch_start = time.time()
          yield cls.fromRNEFflist(batch)
          remaining_iter = math.ceil((len(flist) - batch_end)/batch_size)
          time_passed, time_remained = execution_time2(process_start,remaining_iter,number_of_iteration)
          print(f'Processed {batch_size} files out of {len(flist)} in {execution_time(batch_start)}')
          print(f'Overall processing time: {time_passed}\nEstimated remaining time for {len(flist)-batch_end} out of {len(flist)} files: {time_remained}\n')


  def tree4(self,root:PSObject,reverse=False):
      tree_rn = ResnetGraph()
      try:
          tree = nx.algorithms.traversal.breadth_first_search.bfs_tree(self,root.uid(),reverse)
          assert(isinstance(tree,nx.DiGraph))
          tree_nodes = self._get_nodes(list(tree.nodes()))
          tree_rn.add_psobjs(set(tree_nodes))
          for r,t in tree.edges():
              rels = self._psrels4(r,t)
              [tree_rn.__copy_rel(rel) for rel in rels]
          return tree_rn
      except NetworkXError:
          print(f'{root.name()} nodes was not found in graph {self.name}')
          return ResnetGraph()
      


  def largest_tree(self):
      largest_tree = nx.DiGraph()
      root_nodes = self.root_nodes()
      for root in root_nodes:
          tree = nx.algorithms.traversal.breadth_first_search.bfs_tree(self,root)
          if len(tree) > len(largest_tree):
              largest_tree = tree

      largest_tree_rn = ResnetGraph()
      largest_tree_rn.add_psobjs(set(self._get_nodes(largest_tree.nodes())))
      for r,t in largest_tree.edges():
          rels = self._psrels4(r,t)
          [largest_tree_rn.__copy_rel(rel) for rel in rels]
      
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
              paths = [p for p in nx.simple_paths.all_simple_paths(self,root_id,end_id)]
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
          

  def idpath2objs(self,uid_path:list,rel2parent=True):
      """
      Input
      -----
      path - list of node uids specifying the path in the graph\n
      if "rel2parent" is True relation between two nodes is merged to predecessor node\n
      otherwise relation is merged to successor

      Return
      ------
      [PSObject] - list of PSObject in th path with relations merged either to predecessor or successor nodes 
      """
      path_objs = list()
      if rel2parent:
          for i in range(1,len(uid_path)):
              parent_uid = uid_path[i-1]
              child_uid = uid_path[i]
              parent = self._psobj(parent_uid)
              path_rels = self.__find_relations(parent_uid,child_uid)
              [parent.merge_obj(rel.rel2psobj()) for rel in path_rels]
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
              [child.merge_obj(rel.rel2psobj()) for rel in path_rels]
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
          if using_key not in path_element:
              path_element[using_key] = add_dict # ???????
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
                  shortest_paths = [p for p in nx.simple_paths.all_simple_paths(self,node.uid(),leaf)]

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
                              child_top_level_class = child_node.get_prop(top_level_class_prop)
                              if child_top_level_class:
                                  curent_path = add_element(dict2return,curent_path,{},child_top_level_class)

                              if not child_key: break
                              curent_path = add_element(dict2return,curent_path,child_prop_dict,child_key)
                      
      return dict2return

########################  SUBGRAPH SUBGRAPH SUBGRAPH #####################################
  def subgraph(self,node_uids:list):
      '''
      Returns subgraph made by nx.subgraph containg all edges between node_ids
      '''
      sub_g = ResnetGraph(super().subgraph(node_uids))
      subg_rels = sub_g._psrels()
      all_uids = set()
      [all_uids.update(r.entities_uids()) for r in subg_rels] 
      all_nodes = set(self._get_nodes(all_uids))
      # need to collect nodes from ChemicalReaction
      sub_g.add_psobjs(all_nodes)
      sub_g.load_urn_dicts()
      return sub_g


  def neighborhood(self,_4psobs:set[PSObject],only_neighbors:list[PSObject]=[],
                    only_reltypes:list[str]=[],with_effects:list[str]=[],in_direction=''):
      '''
      Input
      -----
      _4psobs, only_neighbors - [PSObject]\n
      in_direction = ['>','<',None], defaults to None
      '''
      neighbors = self.get_neighbors(_4psobs,only_neighbors)
      return self.get_subgraph(list(_4psobs),list(neighbors),only_reltypes,with_effects,in_direction)

  
  def subtract(self, other: "ResnetGraph"):
      '''
      only self graph is analyzed, works faster if other graph is bigger than self
      '''
      edges_from_other = other._psrels()
      nodes2add  = set()
      rels2add = set()
      for n1,n2,rel in self.iterate():
          if rel not in edges_from_other:
              nodes2add.update([n1,n2])
              rels2add.add(rel)

      unique2self = ResnetGraph()
      #unique2self.add_psobjs(nodes2add)
      unique2self.add_psrels(rels2add)
      return unique2self


  def intersect(self, other: "ResnetGraph"):
      #only self graph is analyzed 
      # works faster if other graph is bigger than self
      intersection = ResnetGraph()
      edges_from_other = other._psrels()
      for n1,n2,e in self.edges(data='relation'):
          if e in edges_from_other:
              intersection.add_psobjs({self.nodes[n1],self.nodes[n2]},merge=False)
              intersection.__add_rel(e['relation'])

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
              assert(isinstance(rel,PSRelation))       
              if rel.objtype() in rel_types:
                  rel.filter_references(keep_prop2values)
      else:
          for regulatorID, targetID, rel in filtered_graph.edges.data('relation'):
              rel.filter_references(keep_prop2values)

      return filtered_graph

                  
  def regulatory_network_urn(self,for_targets_with_urns:list[str],network_name:str):
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
          if in_property in psobj:
              identifiers = psobj[in_property]
              no_version_set = set()
              for i in identifiers:
                  version_pos = str(i).rfind('.')
                  id_end = version_pos if version_pos > 0 else len(i)
                  no_version_set.add(i[:id_end])
              uid2value[uid] = list(no_version_set)

      graph2return = self.copy()
      nx.function.set_node_attributes(graph2return,uid2value,in_property)
      return graph2return


  @staticmethod
  def __bestrel(rels:list,ranks:list[list[str]]=[]) -> PSRelation:
      '''
      Input
      -----
      rels - [PSRelation]
      ranks - [[str,..]], where str - PSRelation obtypes sorted by rank. Order of ranking lists may be important.\n
      Examples of ranks:\n
      [['DirectRegulation','Binding','ProtModification','Regulation'],
      ['PromoterBinding','Expression','Regulation'],
      ['Biomarker','StateChange'],
      ['Biomarker','QuantitativeChange','FunctionalAssosiation'],
      ['MolSynthesis','Regulation'],
      ['MolTransport','Regulation'],
      ['MolTransport','CellExpression'],
      ['Regulation','FunctionalAssosiation']]
      '''
      for rank_list in ranks:
          my_rels = [r for r in rels if r.objtype() in rank_list]
          if len(my_rels) > 1:
              my_rels.sort(key=lambda x: int(x[REFCOUNT][0]), reverse=True)
              my_effect = my_rels[0].effect()
              if my_effect == 'unknown':
                  for rel in my_rels[1:]:
                      if rel.effect() != 'unknown':
                          my_effect = rel.effect()
                          break

              for reltype in rank_list:
                  for rel in my_rels:
                      if rel.objtype() == reltype:
                          rel[EFFECT] = [my_effect]
                          return rel
                      
          elif len(my_rels) == 1:
              return my_rels[0]

      # we are here because rels cannot be ranked by any ranks
      # case when relations are annotated with pX from Reaxys
      rels.sort(key=lambda x: x.pX(), reverse=True)
      if rels[0].pX() > 0.0:
          best_rel = rels[0]
      else:
          rels.sort(key=lambda x: x.count_refs(), reverse=True)
          best_rel = rels[0] 
          for rel in rels:
              if rel.is_directional():
                  best_rel = rel
                  break
              
      if best_rel.effect() == 'unknown':
          rels.sort(key=lambda x: x.count_refs(), reverse=True)
          for rel in rels:
              rel_effect = rel.effect()
              if rel_effect != 'unknown':
                  best_rel[EFFECT] = [rel_effect]
                  return best_rel
      return best_rel
              

  def __set_bestrel(self, from_rels:list[PSRelation],ranks:list[list[str]]=[]):
      '''
      Input
      -----
      from_rels = [PSRelation]\n
      ranks - [[str,..]], where str - PSRelation obtypes sorted by rank.
      '''
      best_rel = ResnetGraph.__bestrel(from_rels,ranks)
      refs2add = [ref for rel in from_rels for ref in rel.refs()]
      # refs2add must be list because some refs may come from the same article and have the same __hash__
      best_rel._add_refs(refs2add)

      [self.remove_relation(rel) for rel in from_rels]
      self.__add_rel(best_rel,refresh_urn=True)


  def make_simple(self, rel_type_rank:list[str]=[]):
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
      print(f'Simplifying {self.name} graph')
      ranks = [rel_type_rank] if rel_type_rank else []
      simple_g = self.copy()
      for ruid, tuid in self.edges():
          reg2target_rels = simple_g._psrels4(ruid,tuid)
          if len(reg2target_rels) > 1: # need simplification
              simple_g.__set_bestrel(reg2target_rels,ranks)
          
      print('%d redundant edges in graph "%s" were removed by simplification' % 
            (self.number_of_edges()-simple_g.number_of_edges(),self.name))

      # now removing non-directional relations created by direction duplication in __add_rel
      duplicate_rels = [r for u,v,urn,r in simple_g.edges(keys=True,data='relation') if urn not in simple_g.urn2rel.keys()]
      size_before = simple_g.number_of_edges()
      [simple_g.remove_relation(r) for r in duplicate_rels]

      print('%d additional non-directional edges in graph "%s" were removed as duplicates of removed relation' % 
            (size_before-simple_g.number_of_edges(),self.name))
      return simple_g
  

  def curate(self, ranks:list[list[str]]=[]):
      """
      Input
      -----
      default_curation_rules = [
      ['DirectRegulation','Binding','ProtModification','Regulation'],
      ['PromoterBinding','Expression','Regulation'],
      ['Biomarker','StateChange'],
      ['Biomarker','QuantitativeChange','FunctionalAssosiation'],
      ['MolSynthesis','Regulation'],
      ['MolTransport','Regulation'],
      ['MolTransport','CellExpression'],
      ['Regulation','FunctionalAssosiation']
      ]
      
      Returns
      -------
      graph with only one edge between nodes.\n
      keeps relation with the biggest reference count.\n
      all other relations are merged into the most referenced one
      if rel_type_rank is specified the new relation type is assigned accordingly 
      """
      curation_rules = [
      ['DirectRegulation','Binding','ProtModification','Regulation'],
      ['PromoterBinding','Expression','Regulation'],
      ['Biomarker','StateChange'],
      ['Biomarker','QuantitativeChange','FunctionalAssosiation'],
      ['MolSynthesis','Regulation'],
      ['MolTransport','Regulation'],
      ['MolTransport','CellExpression'],
      ['Regulation','FunctionalAssosiation']
      ] if not ranks else ranks

      print(f'Curating {self.name}')
      curated_g = self.copy()
      for ruid, tuid in self.edges():
          # if regulator_uid ==PSObject.urn2uid('urn:agi-cas:106266-06-2'):
          #     if  target_uid == PSObject.urn2uid('urn:agi-llid:3351'):
          #         print('')
          reg2target_rels = curated_g._psrels4(ruid,tuid)
          if len(reg2target_rels) > 1: # need simplification
              curated_g.__set_bestrel(reg2target_rels,curation_rules)
              
      print('%d redundant edges in graph "%s" were removed by simplification' % 
          (self.number_of_edges()-curated_g.number_of_edges(),self.name))

      # now removing non-directional relations created by direction duplication in __add_rel
      duplicate_rels = [r for u,v,urn,r in curated_g.edges(keys=True,data='relation') if urn not in curated_g.urn2rel.keys()]
      size_before = curated_g.number_of_edges()
      [curated_g.remove_relation(r) for r in duplicate_rels]

      print('%d additional non-directional edges in graph "%s" were removed as duplicates of removed relation' % 
          (size_before-curated_g.number_of_edges(),self.name))
      return curated_g
  

  def clean(self):
      """
      Merges relations of same type between regulator and target to remove duplicates by Effect,Mechanism,ChangeType
      """
      print(f'Cleaning {self.name}')
      curated_g = self.copy()
      for ruid, tuid in self.edges():
          reg2target_rels = curated_g._psrels4(ruid,tuid)
          rtype2rels = dict()
          for rel in reg2target_rels:
              reltype = rel.objtype()
              try:
                  rtype2rels[reltype].append(rel)
              except KeyError:
                  rtype2rels[reltype] = [rel]

          for rtype, sametype_rels in rtype2rels.items():
              if len(sametype_rels) > 1: # need simplification
                  curated_g.__set_bestrel(sametype_rels,[[rtype]])
              
      print('%d redundant edges in graph "%s" were removed by merging relations of the same type' % 
          (self.number_of_edges()-curated_g.number_of_edges(),self.name))
  
      return curated_g
  

  def regulome_dict(self, only_objtype:list[str], min_size=2):
      """
      Returns
      -------
      {regulator_id : [targets]}, where len([targets]) >= min_size
      """
      regulators = list(self.regulators(only_objtype,min_targets=min_size))
      subnetworks = dict()
      len_targets = dict()
      for regulator in regulators:
          regulator_id = regulator.uid()
          target_ids = list(self.neighbors(regulator_id))
          if target_ids:
              targets = self._get_nodes(target_ids)
              subnetworks[regulator_id] = targets
              len_targets[regulator_id] = len(target_ids)
      
      nx.function.set_node_attributes(self, len_targets, NUMBER_OF_TARGETS)
      print('Generated %d regulome subnetworks with more than %d targets from %s' 
                      % (len(subnetworks), min_size, self.name))
      return subnetworks


  def upstream_relations(self,node_uid:int,with_types=list()):
      my_rels = [rel for r,t,rel in self.edges.data('relation') if t == node_uid]
      if with_types:
          my_rels = [r for r in my_rels if r[OBJECT_TYPE] in with_types]
      return my_rels


  def find_targets(self,of_regulators:list=[],targets_objtype:list=[],linkedby_reltypes:list=[], min_regulators=1):
      '''
      Input
      -----
      of_regulators - [PSObject]

      Return
      ------
      tuple [PSObject] [PSRelation] [PSObject],\nwhere
      0 - regulators
      1 - relations between regulators and targets
      2 - targets
      '''
      regulators_uids = set()
      targets_uids = set()
      if of_regulators:
          regulators_uids = self.uids(of_regulators)
          targets_uids = [t for r,t in self.edges() if r in regulators_uids]
      else:
          for r,t in self.edges():
              regulators_uids.add(r)
              targets_uids.add(t)

      if targets_objtype:
          targets = self.__psobjs(targets_uids)
          targets = {t for t in targets if t.objtype() in targets_objtype}
          targets_uids = self.uids(list(targets))

      relation2targets = set()
      if linkedby_reltypes:
          filetered_edges = [(r,t,rel) for r,t,rel in self.edges.data('relation') if r in regulators_uids and t in targets_uids and rel.objtype() in linkedby_reltypes]
          filtered_regulator_uids = set()
          filtered_target_uids = set()

          for e in filetered_edges:
              filtered_regulator_uids.add(e[0])
              filtered_target_uids.add(e[1])
              relation2targets.add(e[2])
          
          regulators_uids = filtered_regulator_uids
          targets_uids = filtered_target_uids
      else:
          [relation2targets.add(rel) for r,t,rel in self.edges.data('relation') if r in regulators_uids and t in targets_uids]

      if min_regulators > 1:
          reg2target_subgraph = self.subgraph_by_rels(list(relation2targets))
          targets_uids = [uid for uid in reg2target_subgraph.nodes() if reg2target_subgraph.in_degree(uid) > min_regulators]
          filtered_regulator_uids = set()
          filtered_rels = set()
          for r,t,rel in self.edges.data('relation'):
              if t in targets_uids:
                  filtered_regulator_uids.add(r)
                  filtered_rels.add(rel)
          regulators_uids = filtered_regulator_uids
          relation2targets = filtered_rels
      
      return self.__psobjs(regulators_uids), relation2targets, self.__psobjs(targets_uids)
  

  def downstream_relations(self,node_id:int,with_types=list()):
      my_rels = [rel for r,t,rel in self.edges.data('relation') if r == node_id]
      if with_types:
          my_rels = [r for r in my_rels if r[OBJECT_TYPE] in with_types]
      return my_rels


  def upstream_regulators(self,of_node_id:int,linkedby_reltypes=list()):
      regulator_ids = [r for r,t,rel in self.edges.data('relation') if t == of_node_id and rel[OBJECT_TYPE][0] in linkedby_reltypes]
      return regulator_ids


  def get_regulome(self, start_nodes:set[PSObject]):
      """
      Input
      -----
      start_nodes - {PSobject}

      Returns 
      -------
      ResnetGraph() - composition of bfs_trees orginated from every node in "start_nodes"
      """
      all_trees = nx.DiGraph()
      for node in start_nodes:
          if  self.has_node(node): 
              t = nx.algorithms.traversal.breadth_first_search.bfs_tree(self, node.uid())
              all_trees = nx.algorithms.operators.binary.compose(all_trees, t)
      
      regulome_edges = list(all_trees.edges())
      regulome_rels = [rel for r,t,rel in self.edges.data('relation') if (r,t) in regulome_edges]
      return self.subgraph_by_rels(regulome_rels)


  def subgraph_by_relurns(self, rel_urns:set):
      """
      #not tested!!
      """
      subgraph = ResnetGraph()
      rels2add = [r for u,r in self.urn2rel.items() if u in rel_urns]
      uids2add = set()
      [uids2add.update(r.entities_uids()) for r in rels2add]
      subgraph.add_psobjs(set(self._get_nodes(self._get_nodes(uids2add))))
      [subgraph.add_rel(rel) for rel in rels2add]
      return subgraph


  def subgraph_by_rels(self, rels:list)->'ResnetGraph':
      """
      Input
      -----
      rels = [PSRelation]
      """
      node_uids = set()
      [node_uids.update(rel.entities_uids()) for rel in rels]
      subgraph = ResnetGraph()
      if node_uids:
          rel_nodes = self._get_nodes(node_uids)
          subgraph.add_psobjs(set(rel_nodes),merge=False) # has to add original objects from self 
          subgraph.add_psrels(set(rels),merge=False)

      return subgraph


  def subgraph_by_relprops(self, search_values:list, in_properties:list=[OBJECT_TYPE]):
      '''
      Return
      ------
      if "search_values" is empty returns subgraph with all relations annotated by "in_properties"
      '''
      search_value_dict = {p:search_values for p in in_properties}
      my_rels = set()
      if search_values:
          for n1, n2, rel in self.iterate():
            if rel.has_value_in(search_value_dict):
                my_rels.add(rel)
      else:
        for prop_name in in_properties:
          for n1, n2, rel in self.iterate():
            if rel.has_properties({prop_name}):
                my_rels.add(rel)

      print(f'Select subgraph with {len(my_rels)} relations with {in_properties} from graph with {self.number_of_edges()} edges')
      return self.subgraph_by_rels(list(my_rels))


  def subgraph_by_refcount(self,min_refcount:int, max_refcount:int):
      my_rels = self._psrels()
      need_rels = [rel for rel in my_rels if (min_refcount <= rel.count_refs() <= max_refcount)]
      return self.subgraph_by_rels(need_rels)
  

  def subgraph_by_nodeprops(self,has_properties:list,with_values=set()):
      '''
      Returns
      -------
      neighborhood of nodes that has_properties if with_values is empty\n
      otherwise neighborhood of nodes that has_properties with_values
      '''
      my_nodes = self.psobjs_with(has_properties,with_values)
      return self.neighborhood(set(my_nodes))


  def get_subgraph(self,between_nodes:list[PSObject],and_nodes:list[PSObject],
          by_relation_types:list[str]=[],with_effect:list[str]=[],in_direction='')->"ResnetGraph":
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

      subgraph = self.subgraph_by_rels(list(my_rels))
      return subgraph


  def __ontology_resnet(self,children:list,add2parent:PSObject):
      """
      Input
      -----
      children - [PSObjects]
      """
      resnet = ResnetGraph()
      resnet.add_psobj(add2parent)
      resnet.add_psobjs(set(children))
      for child in children:
          #child_id = m.uid()
          rel = PSRelation({OBJECT_TYPE:['MemberOf'],'Relationship':['is-a'],'Ontology':['Pathway Studio Ontology']})
          rel.Nodes[REGULATORS] = [child]
          rel.Nodes[TARGETS] = [add2parent]
          resnet.add_rel(rel)

      return resnet


  def ontology_graph(self):
      '''
      Input
      -----
      nodes in self must be annotated with [CHILDS] property

      Return
      ------
      ResnetGraph with edges (child,parent,relation=[MemberOf, is_a, Pathway Studio Ontology]
      '''
      ontology_graph = ResnetGraph()
      parents = self.psobjs_with([CHILDS])
      for p in parents:
          ontology_graph.add_psobj(p)
          ontology_graph.add_psobjs(p[CHILDS])
          for child in p[CHILDS]:
              rel = PSRelation({OBJECT_TYPE:['MemberOf'],'Relationship':['is-a'],'Ontology':['Pathway Studio Ontology']})
              rel.Nodes[REGULATORS] = [child]
              rel.Nodes[TARGETS] = [p]
              ontology_graph.add_rel(rel,merge=False)
      return ontology_graph


  def all_paths_from(self,child:PSObject)->list[list[PSObject]]:
      '''
      Input
      -----
      self must be ontology_graph

      Returns
      -------
      [[PSObject]] - list of all pathes (lists) from child\n
      used to retrieve all ontology branches containing input child
      '''
      parent_tree = self.tree4(child)
      top_parent_uids = {x for x in parent_tree.nodes() if not parent_tree.out_degree(x)}
      all_ontology_paths = list()
      for parent_uid in top_parent_uids:
          ontology_uid_paths = list(nx.simple_paths.all_simple_paths(parent_tree,child.uid(),parent_uid))
          for path in ontology_uid_paths:
              ontology_obj_path = self.idpath2objs(list(path))
              all_ontology_paths.append(ontology_obj_path)

      return all_ontology_paths


  def direct_indirect_targets(self,of_node_id:int
                      )->tuple[set[tuple[int,float,float]],set[tuple[int,float,float]]]:
      '''
      output:
          tuple: {direct_neighbors}, {indirect_neighbors_ids}\n
          where boths sets = {(target_uid,pX,consistency)}, pX = -1.0 if absent
      '''
      direct_targets = set()
      indirect_targets = set()
      for neighbor_id in self[of_node_id]:
          for urn in self[of_node_id][neighbor_id]:
              rel = self[of_node_id][neighbor_id][urn]['relation']
              assert(isinstance(rel,PSRelation))
              rel_is_direct = rel.isprimarytarget()
              pX = float(rel.get_prop('Affinity',0,-1.0))
              consistency = round(float(rel.get_prop(CONSISTENCY,0,0.0)),3)
              if rel_is_direct == DIRECT:
                  direct_targets.add((neighbor_id,pX,consistency))
              elif rel_is_direct == INDIRECT:
                  indirect_targets.add((neighbor_id,pX,consistency))
              else:
                  continue

      direct_targets = list(direct_targets)
      direct_targets.sort(key=lambda x:(float(x[1]),self._get_node(x[0]).name()),reverse=True)
      indirect_targets = list(indirect_targets)
      indirect_targets.sort(key=lambda x:(float(x[1]),self._get_node(x[0]).name()),reverse=True)
      return direct_targets, indirect_targets


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
      components = nx.algorithms.components.connected_components(my_graph)

      large_comminity_count = 0
      #component_graph = ResnetGraph()
      for i, component in enumerate(components):
          component_graph = my_graph.subgraph(component)
          large_comminity_count += 1
          puddles = list(nx.algorithms.community.asyn_fluidc(component_graph,100))
          print('Created %d communities from component with %d nodes and %d edges' % 
          (len(puddles),component_graph.number_of_nodes(),component_graph.number_of_edges()))
          communities += puddles

          component_rels = [e for r,t,e in component_graph.edges(data='relation')]
          small_components_rels.update(component_rels)

          print('Found %d components' % i)

      communities.sort(key=len,reverse=True)
      communities_subgraphs = list()
      for community in communities:
          community_subgraph = self.subgraph(community)
          communities_subgraphs.append(community_subgraph)

      small_components_graph = self.subgraph_by_rels(list(small_components_rels))
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
      nodes4partion1,nodes4partion2 = nx.algorithms.community.kernighan_lin_bisection(nx.Graph(graph))
      kl_partion1 = graph.subgraph(nodes4partion1)
      kl_partion2 = graph.subgraph(nodes4partion2)

      partion1_rels = {e for r,t,e in kl_partion1.edges(data='relation')}
      partion2_rels = {e for r,t,e in kl_partion2.edges(data='relation')}

      kl_partion1_graph = graph.subgraph_by_rels(list(partion1_rels))
      kl_partion2_graph = graph.subgraph_by_rels(list(partion2_rels))

      common_rels = graph._psrels().difference(partion1_rels|partion2_rels)
      overlap_graph = graph.subgraph_by_rels(list(common_rels))

      return kl_partion1_graph,kl_partion2_graph, overlap_graph


  def components(self,only_with_nodetypes=list(),only_with_reltypes=list()):
      '''
      Returns
      -------
      [ResnetGraph] - List of unconnected subgraphs that have no common edges sorted by size
      '''
      components = nx.algorithms.components.connected_components(nx.Graph(self))
      components_graphs = list()
      for component in components:
          component_graph = self.subgraph(list(component))
          if component_graph.has_nodetypes(only_with_nodetypes):
              if component_graph.has_reltypes(only_with_reltypes):
                  components_graphs.append(component_graph)

      components_graphs = sorted(components_graphs, key=lambda x:len(x), reverse=True)
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


  def __remap_rels(self,uid_remap:dict[int,list[int]])->list[PSRelation]:
      '''
      Input
      -----
      uid_remap = {current_node_uid:[new_node_uids]}
      '''
      def map_node_tuples(node_tuples:list):
          '''
          Return
          ------
          original node_tuples if no tuples were remapped
          '''
          mapped_node_tuples = list() # list of lists with remapped node tuples [[((new_uid(),'0' or '1' for direction,effect))]]
          for t in node_tuples:
              try:
                  mapped_uids = uid_remap[t[0]]
                  remapped_tuples = [(m_uid,t[1],t[2]) for m_uid in mapped_uids]
                  mapped_node_tuples.append(remapped_tuples)
              except KeyError:
                  mapped_node_tuples.append([t])
          # mapped_node_tuples - [[mappings from tuple1],[mappings from tuple2],] 
          return list(product(*mapped_node_tuples)) # all possible combinations from all tuples mappings
          # each combination produces new relation
          
      remaped_rels = list()
      for n1,n2,rel in self.edges.data('relation'):
          assert(isinstance(rel,PSRelation))
          regulator_combinations = map_node_tuples(rel.Nodes[REGULATORS])
          new_rels = list()
          if TARGETS in rel.Nodes:
              target_combinations = map_node_tuples(rel.Nodes[TARGETS])
              for reg_tuples in regulator_combinations:
                  for targ_tuples in target_combinations:
                      new_rel = rel.copy()
                      new_rel.Nodes[REGULATORS] = list(reg_tuples)
                      new_rel.Nodes[TARGETS] = list(targ_tuples)
                      new_rels.append(new_rel)
          else:
              for reg_tuples in regulator_combinations:
                  new_rel = rel.copy()
                  new_rel.Nodes[REGULATORS] = list(reg_tuples)
                  new_rels.append(new_rel)
          
          assert(new_rels)
          remaped_rels += new_rels

      return remaped_rels
  

  def remap_graph(self,props2newobjs:dict[str,list[PSObject]],map_by_props:list[str],
                  case_insensitive=True, stricrlen=4)->tuple['ResnetGraph',set[PSObject]]:
      '''
      Input
      -----
      props2objs - {propvalue:[PSObject]}, propvalue are tokenized if its lentgh > stricrlen\n
      if case_insensitive props2obj keys must be in lowercase\n

      Tokenized match
      ------------
      chars in '-/\\():.[]' are replaced by whitespace for matching\n

      Return
      ------
      ResnetGraph with nodes from "props2objs" and edges from self, set of unmapped_nodes
      '''
      def tokenize(s:str):
          stop_tokens = '-/\\():.[]'  
          return s.translate(str.maketrans(stop_tokens, ' ' * len(stop_tokens))) if len(s) > stricrlen else s

      def make_lower(s:str):
          return s.lower() if len(s) > stricrlen else s

      my_props2newobjs = {tokenize(k).replace(' ',''):v for k,v in props2newobjs.items()}
      if case_insensitive:
          my_props2newobjs = {make_lower(k):v for k,v in my_props2newobjs.items()}
      
      uid_remap = dict() # {old_uid:[new_uids]}
      unmapped_nodes = set()
      remapped_nodes = set()
      for graph_node in self._get_nodes():
          node_props = graph_node.get_props(map_by_props)
          node_props = list(map(lambda x:tokenize(x),node_props))
          if case_insensitive:
              node_props = list(map(lambda x:make_lower(x),node_props))

          mapped_new_nodes = set()
          for normalized_node_prop in node_props:
              try:
                  mapped_new_nodes.update(my_props2newobjs[normalized_node_prop])
              except KeyError:
                  continue

          if mapped_new_nodes:
              if len(mapped_new_nodes) > 1: #make sure that mapping is done between objects with same object type
                  mapped_new_nodes = [n for n in mapped_new_nodes if n.objtype() == graph_node.objtype()]

              new_nodes = [graph_node.copy().merge_obj(n,replace_identity=True) for n in mapped_new_nodes]
              uid_remap[graph_node.uid()] = ResnetGraph.uids(new_nodes)
              remapped_nodes.update(new_nodes)
          else:
              unmapped_nodes.add(graph_node)

      mapped_graph = ResnetGraph()
      #mapped_graph.add_psobjs(remapped_nodes|unmapped_nodes,merge=False)
      mapped_graph.add_psrels(set(self.__remap_rels(uid_remap)))
      return mapped_graph, unmapped_nodes
  

  def replace_nodes(self,psobj_pairs:list[tuple[PSObject,PSObject]]):
      '''
      Input
      -----
      psobj_pairs - [(node_for_replacement,replace_by_node)]

      Return
      ------
      replaced_graph, replaced_graph.neighborhood(of replaced nodes)
      '''
      new_graph_nodes = set(self._get_nodes())
      old_nodes = {n[0] for n in psobj_pairs}
      new_graph_nodes = new_graph_nodes.difference(old_nodes)
      new_nodes = {n[1] for n in psobj_pairs}
      new_graph_nodes.update(new_nodes)

      uid_remap = defaultdict(list)
      [uid_remap[p[0].uid()].append(p[1].uid()) for p in psobj_pairs]

      replaced_graph = ResnetGraph()
      #replaced_graph.add_psobjs(new_graph_nodes,merge=False)
      replaced_graph.add_psrels(set(self.__remap_rels(dict(uid_remap))))
      print(f'Replaced {len(psobj_pairs)} node pairs')
      return replaced_graph, replaced_graph.neighborhood(new_nodes)


############################# EMBEDDING EMBEDDING EMBEDDING #######################################
  def rn2tensor(self,node_stats: Optional[dict[str,int]]=None,node_features: Optional[list[list[int]]]=None,
                edge_stats: Optional[dict[str,int]]=None,idx2objs: Optional[dict[int,PSObject]]=None):
      '''
      Input
      -----
      node_stats - {objtype:counts}. If empty will be created from self\n
      edge_stats - {reltype:counts}. If empty will be created from self\n
      idx2objs - {node_index:PSobject}. Dictionary of indexed graph nodes. If empty will be created from self\n

      Return
      ------
      {index:PSobject} - dictionary of indexed graph nodes 
      node_features - [[0,1,0,0]...[1,0,0,0]] - 1-hot encoded node types. feature index = index of objtype in node_stats
      edges - [r,t,w,[0,1,0,0]] - regulator index, target index, edge weight, 1-hot encoded edge types\n
      node_stats - {objtype:counts}, sorted by count in descending order
      edge_stats - {reltype:counts}, sorted by count in descending order, edge feature index = index of reltype in edge_stats
      '''
      if node_stats is None:
          node_stats = self.__node_stats() 
      
      if edge_stats is None:
          edge_stats = self.__rel_stats()

      urn2idx = dict()
      if isinstance(idx2objs,dict):
          urn2idx = {o.urn():i for i,o in idx2objs.items()}
      else:
          idx2objs = dict()

      def objidx(obj:PSObject):
          try:
              return urn2idx[obj.urn()]
          except KeyError:
              newidx = len(urn2idx)
              urn2idx[obj.urn()] = newidx
              idx2objs[newidx] = obj
              return newidx
      
      if node_features is None:
          node_features_dict = list(node_stats.keys())
          nodes = self._get_nodes()
          node_features = [[]]*2*len(nodes)
          for n in nodes:
              nobjtype = n.objtype()
              n_features = [0]*len(node_features_dict)
              objtype_idx = node_features_dict.index(nobjtype)
              n_features[objtype_idx] = 1

              active_newidx = objidx(n.make_active())
              node_features[active_newidx] = n_features
              repressed_newidx = objidx(n.make_repressed())
              node_features[repressed_newidx] = n_features
      
      edge_features_dict = list(edge_stats.keys())
      edge_list = list()
      for r,t,e in self.edges.data(True):
          regulator = self._get_node(r)
          target = self._get_node(t)
          r_active_idx = objidx(regulator.make_active())
          r_repressed_idx = objidx(regulator.make_repressed())
          t_active_idx = objidx(target.make_active())
          t_repressed_idx = objidx(target.make_repressed())
          
          psrel = e['relation']
          if isinstance(psrel,PSRelation):
              reltype = psrel.objtype()
              if reltype == 'ClinicalTrial':
                  effect_sign = -1 # hacking ClinicalTrial that do not have Effect sign in Resnet
              else:
                  effect_sign = psrel.effect_sign()

              relweight = psrel.count_refs()

              reltype = psrel.objtype()
              e_features = [0]*len(edge_features_dict)
              reltype_idx = edge_features_dict.index(reltype)
              e_features[reltype_idx] = 1

              if effect_sign > 0:
                  edge_list.append([r_active_idx,t_active_idx,relweight,e_features])
                  edge_list.append([r_repressed_idx,t_repressed_idx,relweight,e_features])
              elif effect_sign < 0:
                  edge_list.append([r_active_idx,t_repressed_idx,relweight,e_features])
                  edge_list.append([r_repressed_idx,t_active_idx,relweight,e_features])
              else:
                  edge_list.append([r_active_idx,t_active_idx,relweight,e_features])
                  edge_list.append([r_repressed_idx,t_repressed_idx,relweight,e_features])
                  edge_list.append([r_active_idx,t_repressed_idx,relweight,e_features])
                  edge_list.append([r_repressed_idx,t_active_idx,relweight,e_features])
              
      return edge_list,node_features,idx2objs,node_stats,edge_stats
  
  
  def rn2hd4dp(self)->tuple[HeteroData,dict[int,PSObject]]:
      '''
      Return
      ------
      HeteroData object with edge_weight and edge_label for drug repurposing usecase, where\n
      edge_weight is L2 normalized\n
      edge_label = 1 for positive edges, -1 for negative edges, 0 for test edges\n
      idx2obj - {index:PSobject} - dictionary of indexed graph nodes\n

      [SmallMol,ClinicalTrial,Disease] edge types are converted into [SmallMol,Regulation,Disease] with edge_label = 1\n
      [SmallMol,Regulation,Disease] edges with effect positive are labeled -1\n
      all other [SmallMol,Regulation,Disease] edges have label = 0
      '''
      node_stats = self.__node_stats()
      # assigning indexes to objects
      urn2idx = dict()
      idx2objs = dict()
      def objidx(obj:PSObject):
          try:
              return urn2idx[obj.urn()]
          except KeyError:
              newidx = len(urn2idx)
              urn2idx[obj.urn()] = newidx
              idx2objs[newidx] = obj
              return newidx
          
      def efsign(rel:PSRelation):# hacking ClinicalTrial that do not have Effect sign in Resnet
          return -1 if rel.objtype() == 'ClinicalTrial' else rel.effect_sign()

      # initializing urn2idx to count total number of nodes with different states
      for r,t,rel in self.edges.data(data='relation'):
          regulator = self._get_node(r)
          target = self._get_node(t)
          assert(isinstance(rel,PSRelation))
          effect_sign = efsign(rel)
          if effect_sign > 0:
              objidx(regulator.make_active())
              objidx(target.make_active())
          elif effect_sign < 0:
              objidx(regulator.make_active())
              objidx(target.make_repressed())
          else:
              objidx(regulator.make_active())
              objidx(target.make_active())
              objidx(regulator.make_repressed())
              objidx(target.make_repressed())
          
      # collecting edge data into {triple_type:list} dicts
      edges_reg = defaultdict(list)
      edges_tar = defaultdict(list)
      edge_weights = defaultdict(list)
      edge_labels = defaultdict(list)

      for r,t,rel in self.edges.data(data='relation'):
          regulator = self._get_node(r)
          target = self._get_node(t)
          rtype = regulator.objtype()
          ttype = target.objtype()
          r_active_idx = objidx(regulator.make_active())
          r_repressed_idx = objidx(regulator.make_repressed())
          t_active_idx = objidx(target.make_active())
          t_repressed_idx = objidx(target.make_repressed())

          assert(isinstance(rel,PSRelation))
          effect_sign = efsign(rel)
          if effect_sign > 0:
              new_edges = [[r_active_idx,r_repressed_idx], [t_active_idx,t_repressed_idx]]
          elif effect_sign < 0:
                  new_edges = [[r_active_idx,r_repressed_idx], [t_repressed_idx,t_active_idx]]
          else:
              new_edges = [[r_active_idx,r_active_idx,r_repressed_idx,r_repressed_idx], [t_active_idx,t_repressed_idx,t_active_idx,t_repressed_idx]]

          reltype = rel.objtype()
          new_reltype = 'Regulation' if reltype == 'ClinicalTrial' else reltype
          triple_type = (rtype,new_reltype,ttype)
          edges_reg[triple_type] += new_edges[0]
          edges_tar[triple_type] += new_edges[1]

          relweight = rel.count_refs()
          new_edges_count = len(new_edges[0])
          edge_weights[triple_type] += [relweight]*new_edges_count
          
          if rtype == 'SmallMol' and ttype == 'Disease':
              if reltype == 'ClinicalTrial':
                  edge_labels[triple_type] += [1]*new_edges_count
              else:
                  if reltype == 'Regulation':  
                      if effect_sign > 0:# negative edges are drug-disease toxicities
                          edge_labels[triple_type] += [-1]*new_edges_count
                      else:# all other drug-disease pairs are for testing
                          edge_labels[triple_type] += [0]*new_edges_count

      # creating HeteroData to return
      data = HeteroData()
      number_of_nodes = len(urn2idx)
      node_types = list(node_stats.keys())
      for nt in node_types:
          data[nt].x = [0]*number_of_nodes

      for idx,obj in idx2objs.items():
          data[obj.objtype()].x[idx] = 1

      triple_types = list(edges_reg.keys())
      for triple_type in triple_types:
          #regulator, relation, target = triple_type
          edges = torch.tensor([edges_reg[triple_type],edges_tar[triple_type]])
          edge_weight = torch.from_numpy(np.asarray(edge_weights[triple_type]))

          data[triple_type].edge_index = edges
          data[triple_type].edge_weight = edge_weight
          if triple_type in edge_labels.keys():
              edge_label = torch.from_numpy(np.asarray(edge_labels[triple_type]))
              data[triple_type].edge_label = edge_label
      
      # normalizing edge_weight:
      all_edge_weights = torch.cat([data[tt].edge_weight for tt in triple_types], dim=0)
      all_edge_weights = all_edge_weights.to(torch.float32)
      normalized_weights = torch.nn.functional.normalize(all_edge_weights, p=2, dim=0)
      split_index = [0]+[data[tt].edge_index.shape[1] for tt in triple_types]
      split_intervals = torch.cumsum(torch.tensor(split_index), dim=0)

      for i in range(0,len(triple_types)):
          triple_type = triple_types[i]
          start, end = split_intervals[i], split_intervals[i+1]
          data[triple_type].edge_weight = normalized_weights[start:end]

      return data, idx2objs


  def merge_refs(self,_2refs:dict[str,Reference]):
    '''
    input:
      _2refs = {'idtype:id':Reference}
    '''
    my_refs = self.load_references()
    ref_dic_size = len(set(_2refs.values()))
    print(f'Adding {len(my_refs)} references to {ref_dic_size} references dictionary')
    merged_counter = 0
    for ref in my_refs:
      was_merged = False
      ref_identifiers = list()
      for id_type,id in ref.Identifiers.items():
        identifier = id_type+':'+id
        try:
          _2refs[identifier]._merge(ref)
          was_merged = True
          break
        except KeyError:
          ref_identifiers.append(identifier)
          continue

      if not was_merged:
        _2refs.update({i:ref for i in ref_identifiers})
      else:
        merged_counter += 1

    print(f'{len(set(_2refs.values()))-ref_dic_size} reference added to dictionary')
    print(f'{merged_counter} references were merged')
    return
        

  def neo4j_df(self)->tuple[set[PSObject],df,df,set[PSObject]]:                
    my_nodes = set(self._get_nodes())
    refsets = list()
   # node_rows = [str2str(o) for o in my_nodes]

    '''
    all_refs = dict()
    self.merge_refs(all_refs)
    all_refs = set(all_refs.values())

    if print_refs:
      for ref in all_refs:
        assert(isinstance(ref,Reference))
        doi_id, refdic = ref.todict(as_str2str=True)
        refdic['URN'] = doi_id
        refdic[OBJECT_TYPE] = 'Reference'
        refdic['Name'] = refdic.pop(TITLE)
        refdic['Number of snippets:int'] = ref.number_of_snippets()
        node_rows.append(refdic)
    '''

    rel_rows = list()
    for r,t,rel in self.iterate():
      relname = rel.name()
      refset_id = "'"+relname+"'"
      reldic,_ = rel.todict()
      reldic.pop(REFCOUNT)
      reldic = str2str(reldic)
      r_urn = r.urn()
      t_urn = t.urn()
      reldic[':START_ID'] = r_urn
      reldic[':END_ID'] = t_urn
      reldic['Name'] = relname
      rel_rows.append(reldic)

      refset = {'URN':refset_id,
                'Name':relname,
                OBJECT_TYPE:'ReferenceSet',
                "Number of references:int":rel.count_refs(), 
                "Number of snippets:int": rel.number_of_snippets()}
      
      refsets.append(refset)
      for ref in rel.refs():
        ref_id = ref.doi_id()
        refrel_urn = f'urn:els-BelongsTo:in-out:{ref_id}:in:{refset_id}'
        refrel_name = f'{ref_id}---BelongTo--->{relname}'
        refrel = {':START_ID':ref_id,':END_ID':refset_id, OBJECT_TYPE:'BelongsTo','URN':refrel_urn,'Name':refrel_name}
        rel_rows.append(refrel)

      reltype = rel.objtype()
      urn1 = f'urn:els-{reltype}:in-out:{r_urn}:in:{relname}'
      name1 = f'{r_urn}----{reltype}--->{relname}'
      refset_rel1 = {':START_ID':r.urn(), ':END_ID':refset_id, OBJECT_TYPE:reltype,'URN':urn1,'Name':name1}
      
      urn2 = f'urn:els-{reltype}:in-out:{relname}:in:{t_urn}'
      name2 = f'{relname}----{reltype}--->{t_urn}'
      refset_rel2 = {':START_ID':refset_id, ':END_ID':t.urn(), OBJECT_TYPE:reltype,'URN':urn2,'Name':name2}
      rel_rows += [refset_rel1,refset_rel2]

    refset_df = df(refsets)
    refset_df = refset_df.dfcopy(rename2={OBJECT_TYPE:':LABEL','URN':':ID'})
    refset_df = refset_df.move_cols({'URN:ID':0})

    relations_df = df(rel_rows)
    relations_df = relations_df.dfcopy(rename2={OBJECT_TYPE:':TYPE'})
    relations_df = relations_df.move_cols({':START_ID':0,':END_ID':1})
    return my_nodes,refset_df,relations_df
  

  def rn2neo4jDump(self,save2dir:str,sep='|'):
    nodopG = self.remove_undirected_duplicates()
    nodes_df,relations_df = nodopG.neo4j_df()
    extension = '.txt' if sep== '|' else '.tsv' if sep== '\t' else '.csv'
    #nodes_df.to_csv(save2dir+f'{self.name}Nodes.neo4j'+extension,sep=sep,index=False,quoting=csv.QUOTE_MINIMAL)
    relations_df.to_csv(save2dir+f'{self.name}Relations.neo4j'+extension,sep=sep,index=False,quoting=csv.QUOTE_MINIMAL)


  @staticmethod
  def _make_map(using_props:list, from_nodes:list=[],norm=True)->dict[str,dict[str,dict[str,list[PSObject]]]]:
    '''
    output:
      {objectype:{propname:{propval:PSObject}}}, where propval is in lowercase()
    '''
    mapdic = defaultdict(dict)
    for n in from_nodes:
      if 'protein:prophash' in n.urn(): continue
      for prop in using_props:
        try:
          vals = n.get_props([prop])
          vals = [x.lower() for x in vals]
          val_dic = {v:n for v in vals}
          try:
            mapdic[n.objtype()][prop].update(val_dic)
          except KeyError:
            mapdic[n.objtype()].update({prop:val_dic})
        except KeyError:
            continue

    if norm:
      dcopy = mapdic.copy()
      for t,pvo in dcopy.items():
        for p,vo in pvo.items():
          mapdic[t][p].update({normalize(v):o for v,o in vo.items()})

    return dict(mapdic)

  def make_map(self,using_props:list, from_nodes:list=[])->dict[str,dict[str,dict[str,list[PSObject]]]]:
    '''
    output:
      {objectype:{propname:{propval:PSObject}}}, where propval is in lowercase()
    '''
    nodes = from_nodes if from_nodes else self._get_nodes()
    return self._make_map(using_props,nodes)
  
    '''
    output:
      {objectype:{propname:{propval:PSObject}}}, where propval is in lowercase()
    '''
    nodes = from_nodes if from_nodes else self._get_nodes()
    mapdic = dict()
    for n in nodes:
      for prop in using_props:
        try:
          vals = n.get_props([prop])
          vals = [x.lower() for x in vals]
          n_type = n.objtype()
          try:
            [mapdic[n_type][prop][v].append(n) for v in vals]
          except KeyError:
            try:
              for v in vals: 
                mapdic[n_type][prop][v] =[n]
            except KeyError:
                mapdic[n_type] = {prop:{vals[0]:[n]}}
              #  for i in range(1,len(vals)):
              #    mapdic[n_type][prop].update({vals[i]:[n.urn()]})
        except KeyError:
            continue
        
    return dict(mapdic)

'''
class Regulome(ResnetGraph):
    def __init__(self,regulator:PSObject,regulome:ResnetGraph):
        super().__init__(regulome)
        self.urn2rel = dict(regulome.urn2rel)
        self.regulator = PSObject(regulator)


        def activation_score(self,score_name:str):
            """
            Input
            -----
            self must be anno
            """
            regulator_activation_score = 0.0
            effect_target_counter = 0
            for target in targets_with_anno:
                target_exp_value, target_exp_pvalue = target[sample_annotation_name][0]
                #target_exp_pvalue = target[annotation][0][1]
                if target_exp_pvalue < 0.05 or (str(target_exp_pvalue) == str(np.nan) and abs(target_exp_value) >= 1.0):
                    reg2target_rels = my_graph._psrels4(reg_uid,target.uid())
                    rel = reg2target_rels[0]
                    sign = rel.effect_sign()
                    if sign != 0:
                        regulator_activation_score = regulator_activation_score + sign*target_exp_value
                        effect_target_counter += 1

            if len(targets_with_anno) >= self.min_subnet_size and effect_target_counter > 0:
                regulator_activation_score = regulator_activation_score/math.sqrt(effect_target_counter)
            else:
                regulator_activation_score = np.nan
            return regulator_activation_score, effect_target_counter
'''
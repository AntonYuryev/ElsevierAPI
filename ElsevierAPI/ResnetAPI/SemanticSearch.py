import time,math,os
from .PathwayStudioGOQL import OQL
from ..ETM_API.references import Reference, NODEWEIGHT, RELWEIGHT
from ..pandas.panda_tricks import np, pd, df,MAX_TAB_LENGTH
from .ResnetGraph import ResnetGraph,PSObject,EFFECT,defaultdict
from .ResnetAPISession import APISession,CURRENT_SPECS
from .ResnetAPISession import DO_NOT_CLONE,BELONGS2GROUPS,NO_REL_PROPERTIES,REFERENCE_IDENTIFIERS,DBID
from ..ETM_API.RefStats import SBSstats
from ..utils import execution_time,sortdict,ThreadPoolExecutor,unpack


COUNTS = 'counts'
PS_BIBLIOGRAPHY = 'EBKGrefs'
TM_BIBLIOGRAPHY = 'TMrefs'
CHILDREN_COUNT = 'Number of ontology children'
RANK = 'Literature rank' # PRTS = Probability of Technical and Regulatory Success
ONTOLOGY_ANALYSIS = 'ontology'
WEIGHTED = 'weighted '
MAX_CHILDS = 11 # (10 children + 1 top-level concept)
TOTAL_REFCOUNT = 'Total Semantic refcount'
RELEVANT_CONCEPTS = 'Relevant concepts'
REFCOUNT_COLUMN = '#sentences with co-occurrence. Link opens most relevant articles in PubMed'
INFO_WORKSHEET = 'info'
INPUT_WOKSHEET = 'Input'
PHENOTYPE_WORKSHEET = 'Phenotype'
TARGET_WEIGHT = 'target weight' # default name for the target node weight for relation reference count 
REGULATOR_WEIGHT = 'regulator weight' # default name for the regulator node weight for relation reference count 


class SemanticSearch (APISession):
  __refCache__='reference_cache.tsv'
  __cntCache__='counts_cache.tsv'
  __concept_name__ = 'Concept name'
  __mapped_by__ = 'mapped_by'
  __resnet_name__ = 'Resnet name'
  _col_name_prefix = WEIGHTED
  _col_name_prefix2 = "RefCount to "
  __temp_id_col__ = 'entity_IDs'
  __rel_dir__ = '' # relation directions: allowed values: '>', '<',''
  __iter_size__ = 1000 #controls the iteration size in sematic reference count
  __print_refs__ = True
  data_dir = ''
  iteration_step = 1000
  min_etm_relevance = 0.0
  __colname4GV__  = ''
  TMsearch = 'SBSsearch'

  
  def __init__(self,*args,**kwargs):
      '''
      input:
          APIconfig - args[0]
          kwargs:
          what2retrieve - defaults NO_REL_PROPERTIES
          what2retrieve other options:
          [DATABASE_REFCOUNT_ONLY,REFERENCE_IDENTIFIERS,BIBLIO_PROPERTIES,SNIPPET_PROPERTIES,ONLY_REL_PROPERTIES,ALL_PROPERTIES]
          connect2server - default True, set to False to run script using data in __pscache__ files instead of database
          max_ontology_parent - default: 10
      '''
      my_kwargs = {
              'what2retrieve':REFERENCE_IDENTIFIERS,
              'max_ontology_parent': 10,
              'init_refstat':True
              }
      my_kwargs.update(kwargs)
      session_kwargs, parameters = self.split_kwargs(my_kwargs)

      super().__init__(*args,**session_kwargs)
      self.params = dict(parameters)
      self.PageSize = 1000
      self.RefCountPandas = df() # stores result of semantic retreival
      self.RefCountPandas._name_ = COUNTS
      self.__Ontology__ = set()

      self.__only_map2types__ = list()
      self.__connect_by_rels__= list()
      self.__rel_effect__ = list()
      self.columns2drop = [self.__temp_id_col__] # columns to drop before printing pandas
      self.__boost_with__ = list()
      self.relprop2weight = dict() # {propname:{propval:weight}}
      self.ConceptsHaveWeights = False
      self.max_ontology_parent = self.params.get('max_ontology_parent',10)

      search_kwargs = {'limit':10,'min_relevance':self.min_etm_relevance}
      self.RefStats = SBSstats(self.APIconfig,**search_kwargs)

      self.report_pandas=dict() # stores pandas used to generate report file
      self.raw_data = dict() # stores raw data pandas used for generation pandas in self.report_pandas
      ontology_file = my_kwargs.get('ontology_file','')
      if ontology_file:
        self.ontology_loader = ThreadPoolExecutor(thread_name_prefix=f'Loading ontology')
        print('Submiting future for')
        self.__ontology_future__ = self.ontology_loader.submit(self._load_ontology, ontology_file)
      else:
        self.__ontology_future__ = None


  def _clone(self, **kwargs):
      '''
      output:
          SemanticSearch object copy of self with copy of self.Graph
      '''
      my_kwargs = dict(kwargs)
      my_kwargs['init_refstat'] = False # no need for RefStat in cloned session
      new_session = SemanticSearch(self.APIconfig,**my_kwargs)
      rel_props = my_kwargs.pop('rel_props',[]) # saving input props to add after __retrieve__
      new_session.__retrieve__(my_kwargs.get('what2retrieve',CURRENT_SPECS))
      new_session.add_rel_props(rel_props)
      new_session.add_ent_props(my_kwargs.pop('ent_props',[]))
      new_session.data_dir = self.data_dir
      new_session.__connect_by_rels__ = self.__connect_by_rels__
      new_session.__rel_effect__ = self.__rel_effect__
      new_session.__rel_dir__ = self.__rel_dir__
      new_session.__boost_with__ = self.__boost_with__
      new_session.relprop2weight = self.relprop2weight
      new_session.iteration_step = self.iteration_step
      new_session.TMsearch = self.TMsearch
      new_session.Graph = self.Graph  # need to copy self.Graph to enable using cache in link2concept
      new_session.add2self = False # controls whether to add fetched relation to self.Graph 
      if self.__rel_effect__:
        new_session.add_rel_props([EFFECT])
      
      new_session.neo4j = self.neo4j
      # cloning is done to avoid adding irrelevant references to self.Graph to avoid in PS_Bibliography worksheet
      # therefore new_session.add2self is set to false
      return new_session
  

  def reset(self):
      score_columns = [col for col in self.RefCountPandas.columns if self._col_name_prefix in col]
      for col in score_columns:
        self.RefCountPandas[col] = 0.0
      print('DataFrame was reset')


  def debug(self):
     return self.params['debug']


  def report_df(self,dfname:str)->df:
    return self.report_pandas[dfname]


  def raw_df(self,dfname:str)->df:
    return self.raw_data[dfname]

  def idtype(self):
    return 'URN' if self.useNeo4j() else DBID 

  def __concept(self,colname:str):
    return colname[len(self._col_name_prefix+self._col_name_prefix2):]
  
  def _refcount_colname(self,concept_name:str):
      '''
      output: 
        self._col_name_prefix2 + concept_name
      '''
      return self._col_name_prefix2 + concept_name
  
  def _weighted_refcount_colname(self,concept_name:str):
      return self._col_name_prefix+self._col_name_prefix2 + concept_name
  
  def _linkedconcepts_colname(self,concept_name:str):
      return 'Linked '+ concept_name + ' count'
  
  def _concept_size_colname(self,concept_name:str):
      return concept_name + ' count'
  
  def _refcount_columns(self,concept_name:str):
      '''
      output:
        [self._refcount_colname(concept_name),
         self._weighted_refcount_colname(concept_name),
         self._linkedconcepts_colname(concept_name),
         self._concept_size_colname(concept_name)]
      '''
      return [self._refcount_colname(concept_name),
              self._weighted_refcount_colname(concept_name),
              self._linkedconcepts_colname(concept_name),
              self._concept_size_colname(concept_name)]


  def drop_refcount_columns(self,counts_df=df.from_pd(pd.DataFrame())):
      """
      Removes
      -----
      columns with refcount prefix before redoing semantic search with weighted reference count
      """
      d = self.RefCountPandas if counts_df.empty else counts_df
      ref_count_columns = self.__refcount_columns(d)
      d.drop(columns=ref_count_columns, inplace=True) #defaults to deleting rows, column must be specified
      print('%d refcount columns were dropped' % len(ref_count_columns))


  def report_name(self):
      return 'semantic refcount'


  def report_path(self,extension=''):
    if extension:
      ext = extension if extension.find('.')>-1 else '.'+extension
    else:
      ext = ''
    return self.data_dir+self.report_name()+ext


  def add_entities(self, entities:list[PSObject]):
    self.Entities4df.update(entities)
    self.explodedEntities4df.update(entities)
    children = unpack(e.childs() for e in entities)
    self.explodedEntities4df.update(children)
    self.Graph.add_psobjs(entities+children)
    

  '''
  def df_entities(self,df2score:df=None,include_children = True):
    """
    output:
      list of all ids from the tuples in column self.__temp_id_col__].\n
      if df2score is empty returns all ids from self.RefCountPandas
    """
    if df2score is None:
      if include_children:
        return self.explodedEntities4df
      else:
        return self.Entities4df
    else:
      if self.__temp_id_col__ in df2score.columns:
        list_of_tuples = list(df2score[self.__temp_id_col__])          
        ids = unpack(list_of_tuples)
        idtype = self.idtype()
        if include_children:
           return set([e for e in self.explodedEntities4df if e.id(idtype) in ids])
        else:
          return set([e for e in self.Entities4df if e.id(idtype) in ids])
      else:
        print(f'DataFrame has no {self.__temp_id_col__} column)')
        return []
      '''

  def entities(self,from_df=df(),id_column='Name',map_by_property='Name'):
      '''
      Return
      ------
      propval2objs = {map_by_property:[PSObject]}
      '''
      my_df = self.RefCountPandas if from_df.empty else from_df
      entity_names = my_df[id_column].to_list()
      propval2objs,_ = self.Graph.get_prop2obj_dic(map_by_property,entity_names)
      return propval2objs


  def expand_entities(self, in_df:df,with_objtypes:list,linked_by=list()):
      my_session = self._clone_session(what2retrieve=NO_REL_PROPERTIES)
      my_session.entProps = []

      oql = 'SELECT Relation WHERE NeighborOf (SELECT Entity WHERE id = ({ids}))'
      oql += f' AND NeighborOf (SELECT Entity WHERE objectType = ({with_objtypes}))'
      if linked_by:
          oql += f' AND objectType = ({linked_by})'

      my_dbids = ResnetGraph.dbids(in_df.entities())
      reqname = f'Expanding entities with {with_objtypes} linked by {linked_by}'
      expanded_entity_graph = my_session.iterate_oql(oql,my_dbids,request_name=reqname)

      def __expand_row(cell):
          row_dbids = set(list(cell)) # row_dbids is a tuple
          enities2expand = expanded_entity_graph.psobj_with_ids(row_dbids)
          enity_neighbors = expanded_entity_graph.get_neighbors(set(enities2expand))
          row_dbids.update(ResnetGraph.dbids(list(enity_neighbors)))
          return tuple(list(row_dbids)), ','.join(ResnetGraph.names(list(enity_neighbors))),len(enity_neighbors)

      my_df_copy = in_df.dfcopy()
      my_df_copy[[self.__temp_id_col__,'ExpandedBy','Expand Size']] = my_df_copy[self.__temp_id_col__].apply(lambda x: pd.Series(__expand_row(x)))
      return my_df_copy


  def load_pandas(self,from_entity_df:df,prop_names_in_header=True,use_cache=False,
                  map2type:list=[],max_childs=MAX_CHILDS,expand2=[],with_rels=[]):
      '''
      Input
      -----
      if prop_names_in_header=True uses column headers as mapping property\n
      iterates through all columns in for mapping
      '''
      self.__only_map2types__ = map2type
      if use_cache:
          try: 
              self.read_cache()
              refcount_df = self.RefCountPandas
          except FileNotFoundError:
              print ('Cannot find cache file %s.  You will have to map entities on input identifiers again!' % self.__cntCache__)
              refcount_df = self.__map_entities(from_entity_df,prop_names_in_header)
      else:
          refcount_df = self.__map_entities(from_entity_df,prop_names_in_header)
      
      if max_childs:
          refcount_df = self.remove_high_level_entities(refcount_df,max_childs)

      if expand2:
          refcount_df = self.expand_entities(refcount_df,expand2,with_rels)
      
      #self.Entities4df.update(self.df_entities(refcount_df))
      return refcount_df
  

  def add_temp_id(self,to_df:df,map2columns:list[str]=['URN','Name'],max_childs=MAX_CHILDS,max_threads=50):
    '''
    input:
      if max_childs=0 (=ALL_CHILDS) loads __temp_id_col__ column for all rows in "to_df"
    output:
      df where all rows have values __temp_id_col__, rows that cannot be mapped on DB identifiers are removed
    '''
    new_pd = to_df.dfcopy()
    tempid_col = self.__temp_id_col__
    mapped_by = self.__mapped_by__
    resnet_name = self.__resnet_name__

    if self.__temp_id_col__ not in new_pd.columns:
      new_pd[tempid_col] = [np.nan] * len(new_pd)

    for map2column in map2columns:
      mapping_values = new_pd.loc[new_pd[tempid_col].isna(), map2column].tolist()
      if not mapping_values: break
      mapped_entities = self._props2psobj(mapping_values,[map2column],get_childs=False)
      [o.update_with_value(mapped_by,map2column) for o in mapped_entities]
    
      children,_ = self.load_children4(mapped_entities,max_childs=max_childs,max_threads=max_threads)
      for c in children:
        self.Graph.nodes[c.uid()][mapped_by] = [map2column] 
        #child can me tagged with parent mapped_by but it will not have values in tempid_col 
        self.Graph.nodes[c.uid()][resnet_name] = [c.name()]

      graph_psobjects = self.Graph._get_nodes(ResnetGraph.uids(mapped_entities))
      if max_childs:
        parents4df = [o for o in graph_psobjects if o.number_of_children() <= max_childs] 
        print(f'{len(graph_psobjects)-len(parents4df)} entities with > ontology children {max_childs}  were removed from further calculations')
      else:
        parents4df = graph_psobjects
      prop2tempids = {o.get_prop(map2column):tuple(o.child_dbids()+[o.dbid()]) for o in parents4df}
      assert(map2column in new_pd.columns)
      new_pd = new_pd.set_index(map2column)
      new_pd[tempid_col] = new_pd[tempid_col].fillna(new_pd.index.to_series().map(prop2tempids))
      new_pd = new_pd.reset_index()

    # sometimes dbids are not available :(
    clean_df = df.from_pd(new_pd.dropna(subset=[tempid_col],inplace=False),to_df._name_)
    clean_df.copy_format(to_df)
    removed_rows = len(new_pd)-len(clean_df)
    if removed_rows:
      names_with_nan = new_pd.loc[new_pd[tempid_col].isna(), 'Name'].tolist()
      print(f'{removed_rows} for {names_with_nan} entities were removed from worksheet because database identifiers cannot be found')
    clean_df.add_entities(parents4df)
    return clean_df
  

  def load_df_neo4j(self,from_entities:list[PSObject]|set[PSObject],max_childs=MAX_CHILDS):
    assert self.useNeo4j(), 'Neo4j database connection was not specified'
    entities_with_children = self.neo4j._load_children_(from_entities,max_childs=max_childs)

    children = set()
    parents4df = []
    for parent in entities_with_children:
      if parent.number_of_children() <= max_childs:
        childs = parent.childs()
        children.update(childs)
        parents4df.append(parent)

    for c in children:
      c[self.__mapped_by__] = 'Name'
      c[self.__resnet_name__] = c.name()

    self.Graph.add_psobjs(set(parents4df)|children)
    remove_count = len(from_entities) - len(parents4df)
    print(f'{remove_count} out of {len(from_entities)} entities with > {max_childs} ontology children were excluded from worksheet')

    urns4tempid = [tuple(o.child_urns()+[o.urn()]) for o in parents4df]
    names = [o.name() for o in parents4df] # "graph_psobjects" and "from_entities" may not have same 'Name' !!!!
    urns = [o.urn() for o in parents4df]
    objtypes = [o.objtype() for o in parents4df]
    refcount_df = df.from_dict({'Name':names,'ObjType':objtypes,'URN':urns,self.__temp_id_col__:urns4tempid})
    refcount_df.add_entities(parents4df)
    return refcount_df


  def load_df(self,from_entities:list[PSObject]|set[PSObject],max_childs=MAX_CHILDS,max_threads=50): 
      # do not use self.max_ontology_parent instead of max_childs to avoid problems for session cloning  
      '''
      input:
       if max_childs == 0 no children are loaded
      output:
          refcount_df with columns 'Name',
          if max_childs > 0 refcount_df will have "self.__temp_id_col__" column
      '''
      [o.update_with_value(self.__mapped_by__,o.name()) for o in from_entities] 
      if max_childs:
        if from_entities[0].is_from_rnef():# database ids were not loaded
          self.load_dbids4(from_entities)

        # load_children4() adds empty PSObject in CHILDS property if parents has children > max_childs 
        children,_ = self.load_children4(from_entities, max_childs=max_childs, max_threads=max_threads)
        for c in children:
          self.Graph.nodes[c.uid()][self.__mapped_by__] = 'Name'
          self.Graph.nodes[c.uid()][self.__resnet_name__] = c.name()

        # load_children4() adds empty PSObject in CHILDS property if parents has children > max_childs 
        graph_psobjects = self.Graph._get_nodes(ResnetGraph.uids(from_entities))
          # child_dbids() will return empty list if o[CHILDS] has empty PSObject in CHILDS property
        parents4df = [o for o in graph_psobjects if o.number_of_children() <= max_childs]
        print(f'{len(graph_psobjects)-len(parents4df)} entities with > {max_childs} ontology children were removed from further calculations')

        child_dbids = [tuple(o.child_dbids()+[o.dbid()]) for o in parents4df]
        names = [o.name() for o in parents4df] # "graph_psobjects" and "from_entities" may not have same 'Name' !!!!
        assert(len(names) == len(child_dbids))
        urns = [o.urn() for o in parents4df]
        objtypes = [o.objtype() for o in parents4df]
        refcount_df = df.from_dict({'Name':names,'ObjType':objtypes,'URN':urns,self.__temp_id_col__:child_dbids})
      else:
        names = [o.name() for o in from_entities]
        urns = [o.urn() for o in from_entities]
        objtypes = [o.objtype() for o in from_entities]
        refcount_df = df.from_dict({'Name':names,'ObjType':objtypes,'URN':urns})
          
      refcount_df.add_entities(parents4df)
      return refcount_df


  def __map_entities(self,EntityPandas:df,prop_names_in_header=False):
      start_mapping_time  = time.time()
      if not prop_names_in_header:
          print('Entity infile does not have header with property names:\nwill use 1st column for mapping as Name and then as Alias')
          EntityPandas.rename(columns={0:'Name'}, inplace=True)
          EntityPandas['Alias'] = EntityPandas['Name']
          self.add_ent_props(['Alias'])
      
      map2types = self.__only_map2types__
      PropName2Prop2psobj = dict()
      RemainToMap = EntityPandas.dfcopy()
      mapped_count = 0
      all_mapped_objs = set() # {PSObject}
      for propName in EntityPandas.columns:
        identifiers = list(map(str,list(RemainToMap[propName])))
        prop2dbid2psobj = self.map_prop2entities(identifiers, propName, map2types, get_childs=True)
        PropName2Prop2psobj[propName] = prop2dbid2psobj
        mapped_count = mapped_count + len(prop2dbid2psobj)
        RemainToMap = RemainToMap[~RemainToMap[propName].isin(list(prop2dbid2psobj.keys()))]
        for dbid2objs in prop2dbid2psobj.values():
          all_mapped_objs.update(dbid2objs)
      print('%d rows out of %d were mapped to %d entities that have at least one connection in Resnet using %s' % 
            (mapped_count,len(EntityPandas),self.Graph.number_of_nodes(),','.join(EntityPandas.columns)))

      def get_entIDs(col_name,cell_value):
          prop2dbid2psobj = PropName2Prop2psobj[col_name]
          try:
              mapped_objs = list(map(PSObject,prop2dbid2psobj[str(cell_value)].values()))
              resnet_names = ','.join(ResnetGraph.names(mapped_objs))
              all_obj_dbids = set(ResnetGraph.dbids(mapped_objs + ResnetGraph.childs(mapped_objs)))
              return resnet_names,str(cell_value),tuple(all_obj_dbids)
          except KeyError:
              return None,None,None
      
      df_name = EntityPandas._name_ if EntityPandas._name_ else COUNTS
      df2return = df()
      df2return._name_ = df_name
      for propName in EntityPandas.columns:
          MappedEntitiesByProp = EntityPandas.dfcopy()
          MappedEntitiesByProp = df.apply_and_concat(MappedEntitiesByProp,propName,get_entIDs,
                                      [self.__resnet_name__,self.__mapped_by__,self.__temp_id_col__])
          MappedEntitiesByProp = df.from_pd(MappedEntitiesByProp[MappedEntitiesByProp[self.__temp_id_col__].notnull()])
          df2return = df2return.append_df(MappedEntitiesByProp)
      ex_time = execution_time(start_mapping_time)
      print ('Mapped %d out of %d identifiers to entities in database in %s\n' % 
            (len(df2return), len(EntityPandas.index), ex_time))
      
      df2return.add_entities(all_mapped_objs)
      return df2return


  def map_prop2entities(self,propValues:list,propName:str,map2types=[],
                        get_childs=False,MinConnectivity=1,max_childs=MAX_CHILDS,
                        max_threads=50)->dict[str,dict[int,PSObject]]:
      """
      Returns
      -------
      prop2psobj = {prop_val:{dbid:PSObject}}
      where PSObject are annotated with self.__mapped_by__,CHILD_DBIDS,CHILD_UIDS properties
      """
      ent_props = self.entProps
      if propName not in ent_props: 
          ent_props.append(propName)

      step = 950 # must be slightly less than 1000 to accomodate names with commas
      iteration_counter = math.ceil(len(propValues) / step)

      print('Will use %d %s identifiers to find entities in %d iterations' % 
            (len(propValues), propName, iteration_counter))

      dbid2entity = dict()
      for i in range(0, len(propValues), step):
          start_time = time.time()
          propval_chunk = propValues[i:i + step]
          query_node = OQL.get_entities_by_props(propval_chunk, [propName], map2types, MinConnectivity)
          zeep_entities = self.get_data(query_node, ent_props, getLinks=False)
          if type(zeep_entities) != type(None):
              dbid2entity_chunk = self._zeep2psobj(zeep_entities)
              dbid2entity.update(dbid2entity_chunk)
              ex_time = execution_time(start_time)
              print("%d in %d iteration found %d entities for %d %s identifiers in %s" % 
                  (i / step + 1, iteration_counter, len(dbid2entity_chunk), len(propval_chunk), propName, ex_time))

      # lazy_child_dict = dict()
      propValues_set = list(map(lambda x: str(x).lower(),propValues))
      def my_prop_vals(psobj:PSObject):
          lower_case_values = set(map(lambda x: str(x).lower(),psobj[propName]))
          intersection = list()
          for propval in lower_case_values:
              try:
                  pval_idx = propValues_set.index(propval)
                  intersection.append(propValues[pval_idx])
              except ValueError:
                  continue
          return intersection

      prop2psobj = dict()
      for psobj in dbid2entity.values():
          prop_values = my_prop_vals(psobj)
          mapped_by_propvalue = propName + ':' + ','.join(prop_values)
          psobj.update_with_value(self.__mapped_by__, mapped_by_propvalue)
          for prop_val in prop_values: 
              try:
                  prop2psobj[prop_val][psobj.dbid()] = psobj
              except KeyError:
                  prop2psobj[prop_val] = {psobj.dbid(): psobj}

      if get_childs:
          children,_ = self.load_children4(list(dbid2entity.values()),
                                                                max_childs=max_childs,
                                                                max_threads=max_threads)
          dbid2entity.update({child.dbid():child for child in children})
          for c in children:
              self.Graph.nodes[c.uid()][self.__mapped_by__] = 'Name'

      self.Graph.add_nodes_from([(v.uid(), v.items()) for k, v in dbid2entity.items()])
      print('%d out of %d %s identifiers were mapped on entities in the database' % 
            (len(prop2psobj), len(propValues), propName))
      return prop2psobj
  

  def set_how2connect(self,**kwargs):
    '''
    kwargs:
        'connect_by_rels' - desired relation types for connection. Defaults to []
        'with_effects' - desired relation effect for connection [positive,negative]. Deafults to [] (any)
        'in_dir' allowed values: 
            '>' - from entities in row to concept in column
            '<' - from concept in column to entities in row
            '' - any direction (default)
        'boost_with_reltypes' - if not empty adds to refcount references from other relation types listed in [boost_with_reltypes] regardless Effect sign and direction.  Equivalent to merging relations for refcount if at least one relation specified by other parameters exist between entities and concept
        'step' = step for iterate_oql() function. Defaults to 500
        'nodeweight_prop' - if true concepts must be annotated with NODE_WEIGHT. Defaults to False
    '''
    self.__connect_by_rels__ = kwargs.get('connect_by_rels',[])
    self.__rel_effect__ = kwargs.get('with_effects',[])
    self.__rel_dir__ = kwargs.get('in_dir','')
    self.__boost_with__ = kwargs.get('boost_by_reltypes',[])
    self.iteration_step = kwargs.get('step',500)
    self.ConceptsHaveWeights = kwargs.get('nodeweight_prop',False)

    connect_by_rels = self.__connect_by_rels__+self.__boost_with__

    def connect_psapi(nodes1:list[PSObject], nodes2:list[PSObject]):
      return self.connect_entities(set(nodes1),set(nodes2),connect_by_rels,self.__rel_effect__,self.__rel_dir__,step=self.iteration_step)
    
    def connect_neo4j(nodes1:list[PSObject], nodes2:list[PSObject]):
      return self.connect_entities_neo4j(nodes1,nodes2,connect_by_rels,self.__rel_effect__,self.__rel_dir__)
    
    return connect_neo4j if self.useNeo4j() else connect_psapi
  

  def connect(self,my_df:df,concepts:list[PSObject],how2connect) -> ResnetGraph:
    '''
    output:
      connects all entities in my_df with all concepts in one request.
      combines relations from database with relations from self.Graph non-existing in database
      unique to self.Graph relations can come from in RNEF cache or cretaed by metabolite2inhibit
    '''
    start_time = time.time()
    df_entities = my_df.entities()
    DBconnectionsG = how2connect(concepts,df_entities)
    print('%d nodes were linked by %d relations from database supported by %d references in %s' %
        (len(DBconnectionsG),DBconnectionsG.number_of_edges(),DBconnectionsG.weight(),execution_time(start_time)))
    CacheconnectionsG = self.Graph.get_subgraph(concepts,df_entities,self.__connect_by_rels__,self.__rel_effect__,self.__rel_dir__)
    if CacheconnectionsG:
      print('Additional %d nodes were linked by %d relations from cache supported by %d references in %s' %
          (len(CacheconnectionsG),CacheconnectionsG.number_of_edges(),CacheconnectionsG.weight(),execution_time(start_time)))
      
    connectionsG = CacheconnectionsG.compose(DBconnectionsG)
    if self.useNeo4j():
      self.load_references(connectionsG)
    return connectionsG
  

  def __annotate_rels(self, from_graph:ResnetGraph, with_concept_name:str):
      for regulatorID, targetID, rel in from_graph.edges.data('relation'):
          try:
              entity_search_attrs = self.Graph.nodes[regulatorID][self.__mapped_by__]
          except KeyError:
              try:
                  entity_search_attrs = self.Graph.nodes[targetID][self.__mapped_by__]
              except KeyError:
                  continue

          self.Graph.set_edge_annotation(regulatorID, targetID,rel.urn(),self.__mapped_by__, entity_search_attrs)
          self.Graph.set_edge_annotation(regulatorID, targetID,rel.urn(),self.__concept_name__, [with_concept_name])
      return


  def init_concept(self,my_df:df, concept_name:str):
      in_df = my_df.dfcopy()
      refcount_column = self._refcount_colname(concept_name)
      weighted_refcount_column = self._weighted_refcount_colname(concept_name) 
      linked_count_column = self._linkedconcepts_colname(concept_name)
      concept_size_column = self._concept_size_colname(concept_name)

      in_df.insert(len(in_df.columns),weighted_refcount_column,[float(0)]*len(in_df))
      in_df.insert(len(in_df.columns),refcount_column,[0]*len(in_df))
      in_df.insert(len(in_df.columns),linked_count_column,[0]*len(in_df))
      in_df.insert(len(in_df.columns),concept_size_column,[1]*len(in_df))
      return in_df


  def __link2concept(self,ConceptName:str,concepts:list[PSObject],to_entities:df|pd.DataFrame,
                     how2connect)->tuple[ResnetGraph,df]:
    """
    input:
      to_entities.columns must have self.__temp_id_col__\n
      concepts - [PSObject]\n
      how2connect - function with instructions how to connect "concepts","to_entities"
    output:
      linked_row_count -int\n
      linked_entitiess - {PSObject}\n
      copy of "to_entities" df with 4 columns added:
        "weighted Refcount to {ConceptName}", 
        "Refcount to {ConceptName}", 
        "Linked {ConceptName} count", 
        "{ConceptName} count"
    """
    my_df = to_entities.dfcopy() if isinstance(to_entities,df) else df.from_pd(to_entities)

    concept_size = len(concepts)-1 if concepts else 0
    if (len(concepts) > 501 and len(my_df) > 500):
        print(f'"{ConceptName}" concept has {concept_size} ontology children! Linking may take a while, be patient' )
    else:
        print(f'\nLinking row entities to "{ConceptName}" column which has {concept_size} ontology children')

    weighted_refcount_column = self._weighted_refcount_colname(ConceptName) 
    refcount_column = self._refcount_colname(ConceptName)
    linked_count_column = self._linkedconcepts_colname(ConceptName)
    concept_size_column = self._concept_size_colname(ConceptName)

    concepts_uids = set(ResnetGraph.uids(concepts))
    concepts_count = len(concepts_uids)
    try:
        my_df.insert(len(my_df.columns),weighted_refcount_column,[float(0)]*len(my_df))
        my_df.insert(len(my_df.columns),refcount_column,[0]*len(my_df))
        my_df.insert(len(my_df.columns),linked_count_column,[0]*len(my_df))
        my_df.insert(len(my_df.columns),concept_size_column,[concepts_count]*len(my_df)) # assuming all concepts have at least one member 
        # occurence of linked concepts. Format str(#linked concepts)/str(#total concepts)
    except ValueError:
        print('%s column already exists in dataframe!!!' % refcount_column)
        pass

    linked_row_count = 0
    #linked_entities = set()
    start_time  = time.time()
    id_type = 'URN' if self.useNeo4j() else DBID
    connection_graph = self.connect(my_df,concepts, how2connect)
    if connection_graph.number_of_edges()>0:
      self.__annotate_rels(connection_graph, ConceptName)
      ref_counter = set()
      for idx in my_df.index:
        entities_ids = set(my_df.at[idx,self.__temp_id_col__])
        row_has_connection = False
        row_entities = connection_graph.psobj_with_ids(entities_ids,id_type)
        if row_entities:
          if self.__rel_dir__ =='>':
            row_has_connection = connection_graph.relation_exist(row_entities,concepts,
                                          self.__connect_by_rels__,self.__rel_effect__,[],False)
          elif self.__rel_dir__ =='<':
            row_has_connection = connection_graph.relation_exist(concepts,row_entities,
                                          self.__connect_by_rels__,self.__rel_effect__,[],False)
          else:
            row_has_connection = connection_graph.relation_exist(concepts,row_entities,
                                          self.__connect_by_rels__,self.__rel_effect__,[],True)

        if not row_has_connection:
          continue

        row_subgraph = connection_graph.get_subgraph(row_entities, concepts)
        connected_concepts_uids = [c for c in row_subgraph.nodes() if row_subgraph.degree(c) and c in concepts_uids]
        my_df.at[idx,linked_count_column] = len(connected_concepts_uids) # used to calculate concept incidence at normalization step
        #it measures the occurence of concepts linked to row entities among all input concepts
        references = list(row_subgraph.load_references(self.relprop2weight,postgres=self.postgres()))
        # placeholder for possible future use of Scopus citation index:
        # references = [self.RefStats.citation_index(r) for r in references]
        ref_weights = [1.0]*len(references)
        scopus_score = max([x.get('Relation score',[0.0])[0] for x in references])

        if self.ConceptsHaveWeights: # add ref_nodeweights to ref_weights
        # assumes concepts are annotated by REGULATOR_WEIGHT and/or TARGET_WEIGHT
        # currently does not differentiate between regulators and targets because concepts can be upstream and downstream from entities
          regulatorurn2weight = {o.urn():o.get_prop(REGULATOR_WEIGHT) for o in concepts}
          targeturn2weight = {o.urn():o.get_prop(TARGET_WEIGHT) for o in concepts}
          row_subgraph.add_node_weight2ref(regulatorurn2weight,targeturn2weight)
          ref_nodeweights = [r.get_weight(NODEWEIGHT) for r in references]       
          ref_weights = [x + y for x, y in zip(ref_weights, ref_nodeweights)]
        
        if self.relprop2weight:# add ref_relweights to ref_weights
          ref_relweights = [r.get_weight(RELWEIGHT) for r in references]    
          ref_weights =  [x + y for x, y in zip(ref_weights, ref_relweights)]

        row_score = float(sum(ref_weights)) * (1 + scopus_score/100)

        number_of_children = len(list(my_df.at[idx,self.__temp_id_col__]))
        entities_uids = ResnetGraph.uids(row_entities)
        connected_entities_count = len([c for c in row_subgraph.nodes() if row_subgraph.degree(c) and c in entities_uids])
        corrected_row_score = row_score * (1+connected_entities_count/number_of_children)  # boost by multiple component connectivity
        corrected_row_score /= math.sqrt(number_of_children) # normalize by number of entity components
        # correction by the number of connected concepts is done by self.normalize function
        my_df.at[idx,weighted_refcount_column] = corrected_row_score
        my_df.at[idx,refcount_column] = len(references)
        
        ref_counter.update(references)
        linked_row_count += 1

      effecStr = ','.join(self.__rel_effect__) if len(self.__rel_effect__)>0 else 'all'
      relTypeStr = ','.join(self.__connect_by_rels__) if len(self.__connect_by_rels__)>0 else 'all'
      exec_time = execution_time(start_time)
      if linked_row_count:
        print("Concept \"%s\" is linked to %d entities by %s relations of type \"%s\" supported by %d references with effect \"%s\" in %s" %
              (ConceptName,linked_row_count, connection_graph.number_of_edges(),relTypeStr,len(ref_counter),effecStr,exec_time))
      elif connection_graph:
        print("Concept \"%s\" has no links of type \"%s\" with effect \"%s\" to entities in boosted graph" %
            (ConceptName,relTypeStr,effecStr))
    else:  # connection_graph is empty
        print("Concept \"%s\" has no links to entities" % (ConceptName))

    return connection_graph, my_df


  def link2concept(self,to_concept_named:str,concepts:list[PSObject],to_entities:df|pd.DataFrame,
                   how2connect,clone2retrieve=DO_NOT_CLONE)->tuple[ResnetGraph,df]:
    """
    input:
      to_entities.columns must have self.__temp_id_col__\n
      concepts - [PSObject]\n
      how2connect - function with instructions how to connect "concepts","to_entities"
    output:
      linked_row_count - int\n
      linked_entities - list[PSObject]\n
      df copy of "to_entities" with 4 columns added: 
        "weighted Refcount to {to_concept_named}", 
        "Refcount to {to_concept_named}", 
        "Linked {to_concept_named} count", 
        "{to_concept_named} count"
    """
    if clone2retrieve:
      my_session = self._clone(to_retrieve=clone2retrieve,init_refstat=False) # careful to_retrieve MPSVI ClinicalTrial
      connection_graph,return_df = my_session.__link2concept(to_concept_named,concepts,to_entities,how2connect)
      my_session.close_connection()
      return connection_graph,return_df
    else:
      return self.__link2concept(to_concept_named,concepts,to_entities,how2connect)


  def link2RefCountPandas(self,to_concept_named:str,concepts:list,how2connect=set_how2connect,clone2retrieve=DO_NOT_CLONE):
      '''
      Input
      -----
      concepts - [PSObject]
      wrapper for backward compatibility
      '''
      connectionG,self.RefCountPandas = self.link2concept(to_concept_named,concepts,
                                            self.RefCountPandas,how2connect,clone2retrieve)
      linked_row_count = len(connectionG._get_nodes(ResnetGraph.uids(self.RefCountPandas)))
      return linked_row_count


  def flush_dump(self):
      self.flush_dump_files()
      open(self.__cntCache__, 'w').close()
      open(self.__refCache__, 'w').close()


  def add2report(self,table:df):
      '''
      input:
        uses table._name_ attribute to call worksheet in the report
      '''
      if not table.empty:
        table_name = table._name_
        if not table_name:
            table_name = 'Sheet' + str(len(self.report_pandas))
        self.report_pandas[str(table_name)] = table


  def add2raw(self,table:df):
      if not table.empty:
        assert(table._name_)
        self.raw_data[table._name_] = table


  def find_ref(self,ref:Reference):
      id_type, identifier = ref.get_doc_id()
      try:
          return self.references[id_type+':'+identifier]
      except KeyError:
          return dict()


  def read_cache(self):# returns last concept linked in cache
      try:
          self.RefCountPandas = df.from_pd(pd.read_csv(self.__cntCache__,sep='\t',header=0))
          last_col_name = list(self.RefCountPandas.columns)[-1]
          return last_col_name[len(self._col_name_prefix):]
      except FileNotFoundError:
          print('Cache was not found! Will start processing all input files from scratch')
          return FileNotFoundError

###########################  SCORE  ############################   SCORE  #######################
  def __refcount_columns(self,counts_df=df(),column_prefix=''):
      d = self.RefCountPandas if counts_df.empty else counts_df
      refcount_column_prefix = column_prefix if column_prefix else self._col_name_prefix # self._col_name_prefix=WEIGHTED
      to_return = [col for col in d.columns if refcount_column_prefix in col]
      return to_return
  

  def set_rank4(self, concept:str, in_df:df, rank = None):
    '''
      assigns in_df.max_colrank()+1 to in_df.col2rank if rank is None
    '''
    column_name = self._weighted_refcount_colname(concept)
    in_df.set_rank(column_name,rank)


  def rank2weight(self,col2rank:dict[str,int]):
      '''
      input:
        rank is in ascending order: 0 or 1 most important
      '''
      unique_ranks =  sorted(list(set(col2rank.values())))
      number_of_weights = len(unique_ranks)
      weight_step = 1.0/number_of_weights
      r2w = {rank:round((1.0-i*weight_step),2) for i,rank in enumerate(unique_ranks)}
      column2weight = {c:r2w[r] for c,r in col2rank.items()}
      return sortdict(column2weight,by_key=False,reverse=True)
      

  def make_count_df(self,from_rawdf=df(),with_name=COUNTS,sort_by_1stcol=True):
      '''
      Returns
      -------
      df with_name from_df with formatted CHILDREN_COUNT column, soreted by first refcount_column
      '''
      my_df = self.RefCountPandas if from_rawdf.empty else from_rawdf
      count_df = my_df.dfcopy()
      if self.__temp_id_col__ in count_df.columns:
        def sortlen(x):
          return 0 if x is np.nan else max(len(list(x)),1)
        count_df[CHILDREN_COUNT] = count_df[self.__temp_id_col__].apply(lambda x: sortlen(x))
        count_df.add_column_format(CHILDREN_COUNT,'align','center')

      if sort_by_1stcol:
        def __col2sort__(my_df:df):
          refcount_columns = self.__refcount_columns(my_df)
          if refcount_columns: return refcount_columns[0]
          else:
            for c in my_df.columns.tolist():
              if my_df.is_numeric(c): return c

            return str(my_df.columns[1])

        sort_by = __col2sort__(my_df)
        count_df = count_df.sortrows(by=sort_by)

      if self.__colname4GV__ in count_df.columns: # moves column with GVs to the end 
          count_df.move_cols({self.__colname4GV__:len(count_df.columns)})

      count_df._name_ = with_name
      print (f'Created "{count_df._name_}" raw_data count worksheet with {len(count_df)} rows' )
      return count_df
      

  def __adjust4concept_incidence(self,raw_df:df):
      """
      input:
          assumes "raw_df" has following column format for each concept: 
          [weighted Refcount to {concept},Refcount to {concept},Linked {concept} count, {concept} count]

          values in "columns2norm" must start with WEIGHTED
          if "columns2norm" is empty uses list of columns with prefix WEIGHTED
      output:
          copy of raw_df with follwoing columns:
          [weighted Refcount to {concept},Refcount to {concept}, {concept} incidence
      """

      incidence_df = df(name = raw_df._name_)
      rawdf_colnames = raw_df.columns.to_list()
      empty_refcount_cols = list()
      i = 0
      while i < len(rawdf_colnames): # in range() controls i and cannot be used here
          colname = rawdf_colnames[i]
          if colname.startswith(self._col_name_prefix): # startswith(WEIGHTED)
            linked_concepts_count_col = rawdf_colnames[i+2]
            max_concept_count = max(raw_df[linked_concepts_count_col])
            if max_concept_count >= 1:# ensures that refcount column with all zeros is not copied
              incidence_df[colname] = raw_df[colname]
              refcount_colname = rawdf_colnames[i+1]
              incidence_df[refcount_colname] = raw_df[refcount_colname]
          
              concepts_count_col = rawdf_colnames[i+3]
              if max_concept_count > 1: # boosts multicomponent concepts
                incidence_df[colname] = raw_df[colname]*(1+raw_df[linked_concepts_count_col]/raw_df[concepts_count_col])
              
              concept_name = self.__concept(colname)
              incidence_colname = concept_name + ' incidence'
              temp_df = df()
              temp_df['percentage'] = 100*raw_df[linked_concepts_count_col]/raw_df[concepts_count_col]
              temp_df['percentage'] = temp_df['percentage'].round(decimals=2)
              temp_df['percentage_str'] = temp_df['percentage'].apply(lambda x: f"{x:.0f}")+'%'
              incidence_df[incidence_colname] = raw_df[linked_concepts_count_col].astype(str)
              incidence_df[incidence_colname] = incidence_df[incidence_colname] + '/' +raw_df[concepts_count_col].astype(str)
              incidence_df[incidence_colname] = incidence_df[incidence_colname] +'('+temp_df['percentage_str']+')'
            else:
              empty_refcount_cols.append(colname)
            i += 4
          else:
            incidence_df[colname] = raw_df[colname]
            i += 1

      incidence_df.copy_format(raw_df)
      return incidence_df,empty_refcount_cols


  def normalize(self,raw_df_named:str,to_df_named:str,entity_column='Name',
                drop_empty_columns=True,nozero_rank=True,add_pvalue=True)->tuple[df,df]:
    """
    input:
      df with name "raw_df_named" must be in self.raw_data and must have columns:
      [weighted Refcount to {concept},Refcount to {concept},Linked {concept} count, {concept} count]

    output:
      rankedf4report with name "to_df_named" in self.report_pandas[to_df_named]
      normdf4raw with name "norm.to_df_named" in self.raw_data[norm.to_df_named]
    """
    try:
      if self.raw_data[raw_df_named].empty:
        print(f'{raw_df_named} worksheet is empty')
        return df(),df()
    except KeyError:
      print(f'No worksheet named "{raw_df_named}" is available in raw_data')
      return df(),df()
    
    incidence_df,empty_cols = self.__adjust4concept_incidence(self.raw_data[raw_df_named])
    col2rank = {k:v for k,v in incidence_df.col2rank.items() if k not in empty_cols}
    rankedf4report,normdf4raw = self._normalize(incidence_df,col2rank,entity_column,
                                                drop_empty_columns,nozero_rank,add_pvalue)
    rankedf4report._name_ = to_df_named
    self.add2report(rankedf4report)
    return rankedf4report,normdf4raw


  def _normalize(self,incidence_df:df,col2rank:dict=dict(),entity_column='Name',
                 drop_empty_columns=True,nozero_rank=True,add_pvalue=True):

    if not col2rank: col2rank = incidence_df.col2rank
    weight_row = {entity_column:'WEIGHTS:'}
    col2weight = self.rank2weight(col2rank)
    weight_row.update(col2weight)
    weights_df = df([weight_row])
    refcount_header = weights_df.columns.to_list()
    refcount_cols = list(col2weight.keys())
    
    #calculating weighted cumulative score
    normdf4raw = incidence_df.dfcopy(refcount_header)
    normdf4raw = normdf4raw.l2norm(refcount_cols)
    weights = weights_df.loc[0, refcount_cols].values.tolist()
    for i in normdf4raw.index:
      row_scores = normdf4raw.loc[i,refcount_cols].values.tolist()
      assert(len(weights) == len(row_scores))
      weighted_sum = np.nan_to_num(sum(s*w for s,w in zip(row_scores, weights)))
      normdf4raw.loc[i,RANK] = weighted_sum

    normdf4raw = normdf4raw.l2norm([RANK])
    # moving all other columns including CHILDREN_COUNT to normdf4raw 
    # except count columns that were used to create incidence column
    added_columns = set(normdf4raw.columns)
    for col in list(incidence_df.columns):
      if col not in added_columns and ' count' not in col:
        normdf4raw[col] = incidence_df[col]

    if nozero_rank:
      normdf4raw = df.from_pd(normdf4raw.loc[normdf4raw[RANK] >= 0.001]) # removes rows with all zeros
      print(f'Removed {len(incidence_df)-len(normdf4raw)} rows out of {len(incidence_df)} from normalized worksheet with score=0')

    pvals_columns = []
    if add_pvalue:
      rank_pval = RANK+' empirpvalue'
      normdf4raw[rank_pval] = df.calculate_pvalues(normdf4raw[RANK])
      rank_expopval = RANK + ' expopvalue'
      normdf4raw[rank_expopval] = df.calculate_expo_pvalues(normdf4raw[RANK])
      normdf4raw = normdf4raw.sortrows(by=[RANK,rank_pval,entity_column], ascending=[False, True,True])
      pvals_columns = [rank_pval,rank_expopval]


    rename_cols = dict()
    header4rankedf = [entity_column]
    for c in refcount_cols:
      refcount_colname = c
      if c.startswith(self._col_name_prefix):
        # since ranking is done will copy  to rankedf "Refcount to" columns instead of WEIGHTED
        refcount_colname = c[len(self._col_name_prefix):] 
      rename_cols[c] = refcount_colname
      header4rankedf.append(refcount_colname) # adding "Disease model regulatory score"
    
    #prettyfying weighter_df header:
    weights_df = df.from_pd(weights_df.map(lambda x: f'{x:,.4f}' if isinstance(x,float) else x))
    # renaming weight_df columns from "weighted " to "RefCount to ":
    weights_df = weights_df.dfcopy(rename2=rename_cols)

    # re-ordering normdf4raw colums for clarity:
    normdf4raw_cols = normdf4raw.columns.to_list()
    forbidden_cols = set(refcount_header+header4rankedf)
    forbidden_cols.update([self.__temp_id_col__,'URN',RANK,'ObjType',self.__colname4GV__,CHILDREN_COUNT]+pvals_columns)
    # refcount_header has "WEIGHTED ..." columns
    # header4rankedf has "Refcount to ..." columns
    other_cols4rankedf = [x for x in normdf4raw_cols if x not in forbidden_cols]
    header4rankedf += other_cols4rankedf
    header4rankedf +=  [RANK]+pvals_columns+[CHILDREN_COUNT,'URN','ObjType']
   # if self.__colname4GV__ in normdf4raw_cols:
    #    header4rankedf.append(self.__colname4GV__)
    # at this point: header4rankedf = [entity_column]+refcount_colnames+other_cols4rankedf+[RANK,'URN','ObjType',self.__colname4GV__]

    rankedf4report = normdf4raw.dfcopy(only_columns=header4rankedf)
    # now finish normdf4raw by adding weights to df and df to raw data:
    normdf4raw = weights_df.append_df(normdf4raw)
    normdf4raw._name_ = 'norm.'+incidence_df._name_
    self.add2raw(normdf4raw)

    assert(weights_df.columns.to_list() == rankedf4report.columns.to_list()[:len(weights_df.columns)])
    rankedf4report = rankedf4report.sortrows(by=RANK)
    rankedf4report = weights_df.append_df(rankedf4report)# append will add columns missing in weights_str 

    rankedf4report.copy_format(incidence_df)
    rankedf4report.add_column_format(RANK,'align','center')
    if drop_empty_columns:
      rankedf4report = rankedf4report.drop_empty_columns()
    return rankedf4report,normdf4raw


  def clear(self):
    super().clear()
    self.RefCountPandas = df()


  def score_concepts(self,*args,**kwargs)->tuple[ResnetGraph,df]:
    '''
    input:
      args[0] - [PSObject]
      args[1] - df::df2score will score self.RefCountPandas if len(args) == 1\n
      kwargs:
          'concept_name' - string is used to create column name 'Refcount to column_name'
          'connect_by_rels' - desired relation types for connection. Defaults to []
          'with_effects' - desired relation effect for connection [positive,negative]. Deafults to [] (any)
          'in_dir' allowed values: 
              '>' - from entities in row to concept in column
              '<' - from concept in column to entities in row
              '' - any direction (default)
          'boost_with_reltypes' - if not empty adds to refcount references from other relation types listed in [boost_with_reltypes] regardless Effect sign and direction.  Equivalent to merging relations for refcount if at least one relation specified by other parameters exist between entities and concept
          'step' = step for iterate_oql() function. Defaults to 500
          'ConceptsHaveWeights' - indicates the name of the property holding node weight. Deafaults to '' (no weight property
          'clone2retrieve' - defaults to DO_NOT_CLONE
          'column_rank' - column rank to calculate combined score. Set column_rank to -1 to skip ranking
    '''
    concepts = args[0]
    df2score = args[1] if len(args) > 1 else self.RefCountPandas
    if concepts:
      concept_name = kwargs.get('concept_name','Concepts')
      clone2retrieve = kwargs.get('clone2retrieve',DO_NOT_CLONE)
      self.add2self = kwargs.get('add_relevance_concept_column',False)
      my_step = kwargs.get('step',500) # saving old step size
      if len(concepts) > 500:
        #print('reducing step to 250 for large number of concepts (%d)' % len(concepts))
        # to avoid table lock in Oracle:
        kwargs['step'] = 250
      how2connect = self.set_how2connect(**kwargs)
      connectionG,scored_df = self.link2concept(concept_name,concepts,df2score,how2connect,clone2retrieve)
      linked_row_count = len(connectionG._get_nodes(scored_df.uids()))
      print(f'{linked_row_count} rows linked to column {concept_name}')
      kwargs['step'] = my_step
      if linked_row_count:
        # pop ensures next column to have new name assigned explicityly in the upstream function
        rank = kwargs.get('column_rank',None)
        self.set_rank4(concept_name, scored_df, rank) 
        # rank must be set here because different classes use score_concepts or score_concept
      return connectionG,scored_df
    else:
      return ResnetGraph(),df2score
    

  def __explode(self, concepts:list[PSObject], using_concept_params:dict):
    if isinstance(using_concept_params,dict): # concepts have specific parameters (weight, no_childs)
      # parameter json must be in written as "inflammation":[1.0,"no_childs"]

      self.ConceptsHaveWeights = True # signal for __link2concept
      # normalizing concept names for matching:
      concept_params = {k.lower():v for k,v in using_concept_params.items()} 
      concepts_need_childs = []
      for concept in concepts:
        this_concept_params = concept_params[concept.name().lower()]
        if len(this_concept_params) == 1: # concept has only weight parameter
          concepts_need_childs.append(concept)
        else:
          assert(this_concept_params[1] == 'no_childs'), 'Second parameter for concept must be "no_childs"'

      _, concepts_with_childs = self._load_childs_(concepts_need_childs)

      # adding all concepts annoated with children into one array:
      all_concepts = set(concepts_with_childs)
      [all_concepts.add(c) for c in concepts if c not in concepts_with_childs]
      # now concepts in all_concepts are annoated with property CHILDS 

      # now assigning nodeweights from input paramters to all concepts and its children 
      # and collect them for scoring:
      exploded_concepts = set()
      for concept in all_concepts:
        node_weight = concept_params[concept.name().lower()][0]
        concept.set_property(TARGET_WEIGHT,node_weight)
        concept.set_property(REGULATOR_WEIGHT,node_weight)
        exploded_concepts.add(concept)
        for child in concept.childs():
          child.set_property(TARGET_WEIGHT,node_weight)
          child.set_property(REGULATOR_WEIGHT,node_weight)
          exploded_concepts.add(child)
    else: # by default all concepts undergo term expansion:
        children,exploded_concepts = self._load_childs_(concepts)
        exploded_concepts.update(children)
    print(f'{len(concepts)} concepts were exploded to {len(exploded_concepts)}')
    return exploded_concepts
  
  
  def score_concept(self,*args,**kwargs)->tuple[ResnetGraph,df,list[PSObject]]:
    '''
    output:
      connection graph, df with new column, annotated concepts
    input:
      args[0] - key name in self.params that holds the list of concept names\n
      args[1] - df to score.  Wil use self.RefCountPandas if len(args) < 2\n
      kwargs:
        'concept_name' - string is used to create column name 'Refcount to column_name'
        'connect_by_rels' - desired relation types for connection. Defaults to []
        'with_effects' - desired relation effect for connection [positive,negative]. Deafults to [] (any)
        'in_dir' allowed values: 
            '>' - from entities in row to concept in column
            '<' - from concept in column to entities in row
            '' - any direction (default)
        'boost_with_reltypes' - if not empty adds to refcount references from other relation types listed in [boost_with_reltypes] regardless Effect sign and direction.  Equivalent to merging relations for refcount if at least one relation specified by other parameters exist between entities and concept
        'step' = step for iterate_oql() function. Defaults to 500
        'ConceptsHaveWeights' - indicates the name of the property holding node weight. Deafaults to '' (no weight property
        'clone2retrieve' - defaults to DO_NOT_CLONE
        'column_rank' - column rank to calculate combined score
    '''
    start = time.time()
    df2score = args[1] if len(args) > 1 else self.RefCountPandas
    concept_params = self.params.get(args[0],[])
    if not concept_params:
      return ResnetGraph(),df2score,[]
      
    concept_names = list(concept_params.keys())  if isinstance(concept_params,dict) else concept_params
    if concept_names:
      concept_name = args[0]
      print(f'\n\nLinking concepts from "{concept_name}" parameter to "{df2score._name_}" worksheet with {len(df2score)} rows',flush=True)
      print(f'Results will be added to column "{kwargs['concept_name']}"')
     
      if self.useNeo4j():
        input_concepts = self.neo4j.get_nodes('','Name',concept_names)
      else:
        request_name = f'Loading "{concept_name}" from script parameters'
        oql = OQL.get_entities_by_props(concept_names,['Name'])
        input_concepts = self.process_oql(oql,request_name)._get_nodes()

      if input_concepts:
        print(f'Found {len(input_concepts)} {concept_name} in ontology for {len(concept_params)} {concept_name} in parameters')
        exploded_concepts = self.__explode(input_concepts,concept_params)
        connectionG,scored_df = self.score_concepts(exploded_concepts,df2score,**kwargs)
        self.ConceptsHaveWeights = False

        if connectionG:
          if kwargs.get('add_relevance_concept_column',False):
            scored_df = self.add_relevant_concepts(scored_df,{concept_name:list(exploded_concepts)},connectionG)
          
          annotated_inconcepts = [c for c in exploded_concepts if c in input_concepts]
          for c in annotated_inconcepts:
            connectivity = connectionG.connectivity(c,with_children=True)
            c.set_property('Local connectivity',connectivity)
            c.set_property('rank',scored_df.max_colrank())

          print(f'Linked {len(exploded_concepts)} expanded from "{concept_name} to  worksheet "{df2score._name_}" by {connectionG.size()} references')
          print(f'Linking was done in {execution_time(start)}',flush=True)
          return connectionG,scored_df,annotated_inconcepts
        else:
          print(f'No entities were linked to {concept_names}')
          return ResnetGraph(),df2score,[]
      else:
          print(f'No concepts found for {concept_names}. Check your spelling !!!')
          self.ConceptsHaveWeights = False
          return ResnetGraph(),df2score,[]
    else:
      print('No concept name was provided!!!!')
      self.ConceptsHaveWeights = False
      return ResnetGraph(),df2score,[]


##################  ANNOTATE  ############################## ANNOTATE ############################
  def tm_refcount_colname(self,between_column:str,and_concepts:str|list):
    concept_str = and_concepts if isinstance(and_concepts,str) else ','.join(and_concepts)
    return REFCOUNT_COLUMN + ' between '+between_column+' and '+concept_str

  def tm_doi_colname(self,between_names_in_col,and_concepts):
      return self.RefStats.doi_column(between_names_in_col,and_concepts)


  def refs2report(self,to_df_named:str,input_names:list,entity_name_col:str='Name',
                  add2query=[],add2report=True):
    """
    input:
      self.report_pandas[to_df_named] must exist and have column 'Name'
    output:
        copy of self.report_pandas[to_df_named] with added columns:
          RefStats.refcount_column(between_names_in_col,and_concepts)
    """
    multithread = False if self.params.get('debug',False) else True
    names2hyperlinks = self.RefStats.reflinks(self.report_df(to_df_named),entity_name_col,input_names,add2query,multithread)
    refcountcol = self.tm_refcount_colname(entity_name_col,input_names)
    my_df = self.RefStats.add_reflinks(names2hyperlinks,refcountcol,self.report_df(to_df_named),entity_name_col)
    #my_df = self.RefStats.add_refs(self.report_df(to_df_named),entity_name_col,input_names,add2query,multithread)
    if add2report:
      self.add2report(my_df)
    return my_df


  def add_tm_bibliography_df(self,suffix=''):
    """
    Adds
    ----
    df with ETM_BIBLIOGRAPHY name to self.report_pandas from self.RefStats.counter2df()
    """
    if self.params.get('add_bibliography',True):
      biblio_df = self.RefStats.counter2df()
      biblio_df._name_ = TM_BIBLIOGRAPHY
      if suffix: biblio_df._name_ += '-'+suffix
      biblio_df._name_ = biblio_df._name_[:31]
      self.add2report(biblio_df)
      return biblio_df._name_


  def add_graph_bibliography(self,suffix='',from_graph=ResnetGraph()):
    """
    adds:
      df with PS_BIBLIOGRAPHY-suffix name to self.report_pandas
    """
    my_graph = from_graph if from_graph else self.Graph
    ref_df_name = PS_BIBLIOGRAPHY
    if suffix: ref_df_name += '4'+suffix
    ref_df_name = ref_df_name[:31]
    ref_df = my_graph.bibliography(ref_df_name)
    self.add2report(ref_df)
    return


  def id2paths(self,_4df:df, for_entities_in_column='Name', 
                map_by_graph_property='Name',ontology_depth=3)->dict[str,str]:
      '''
      output:
          dict{name: 'parent1->parent2->entity'} for each entity in _4df
      '''

      name2child_objs = self.entities(_4df,for_entities_in_column,map_by_graph_property)
      return self.ontopaths2(name2child_objs,ontology_depth)


  def add_group_annotation(self,group_names:list, _2graph: ResnetGraph = ResnetGraph(),_4groups=True):
    if group_names:
      super().add_group_annotation(group_names,_2graph,_4groups)
    else:
      if self.__ontology_future__:
        if self.__ontology_future__.result():
          urns2values = defaultdict(list)
          [[urns2values[child.urn()].append(group.name()) for child in group.childs()] for group in self.__Ontology__]
          self.Graph.set_node_annotation(dict(urns2values),BELONGS2GROUPS)
      else:
         print('No ontology file was specified !!!')


  def add_groups(self,_2df:df,group_names:list,
      for_entities_in_column='Name',map_by_graph_property='Name',_4groups=True):
      '''
      output:
          copy of _2df with added column 'Groups' containing groups from group_names to which entity belongs
      '''
      self.add_group_annotation(group_names,_4groups=_4groups)
      prop2objs = self.entities(_2df,for_entities_in_column,map_by_graph_property)
      prop2groups = dict()
      for prop, psobjs in prop2objs.items():
        prop_groups = set()
        for psobj in psobjs:
          obj_groups = psobj.propvalues(BELONGS2GROUPS)
          if obj_groups:
            prop_groups.update(obj_groups)
        prop2groups[prop] = ',\n'.join(prop_groups)

      new_df = _2df.dfcopy()
      new_df['Groups'] = new_df[for_entities_in_column].map(prop2groups)
   #   new_df['Groups'] = new_df['Groups'].fillna('')
      self.add2report(new_df)
      return new_df
  

  def add_groups_from_file(self,_2df:df,group_file:str,
    for_entities_in_column='Name',map_by_graph_property='Name',_4groups=False):
    with open(group_file, 'r') as f:
      group_names = [line.strip() for line in f.readlines()]
    return self.add_groups(_2df, group_names, for_entities_in_column, map_by_graph_property,_4groups)
    

  def add_entity_annotation(self,_2column:str,in_df:df,from_node_property='', for_entities_in_column='Name',map_by_graph_property='Name'):
      node_property = from_node_property if from_node_property else _2column
      prop2objs = self.entities(in_df,for_entities_in_column,map_by_graph_property)
      prop2values = dict()
      for prop, psobjs in prop2objs.items():
          prop_values = set()
          for o in psobjs:
              try:
                  o_props = o[node_property]
                  prop_values.update(o_props)
              except KeyError:
                  continue
          
          if prop_values:
              prop2values[prop] = ','.join(prop_values)

      new_df = in_df.dfcopy()
      new_df = new_df.merge_dict(prop2values,_2column,for_entities_in_column)
      self.add2report(new_df)
      return new_df


  def _load_ontology(self, ontology_file:dict={}):
    '''
    input:
      ontology_file = {NodeType:path2file_with_groups}
    '''
    if not self.__Ontology__ and ontology_file:
      assert isinstance(ontology_file,dict), 'Format for ontology file input: {{NodeType:path2file}}'
      for node_type, fpath in ontology_file.items():
        group_names = {x.strip() for x in open(fpath, 'r').readlines()}
        print(f'Loading {len(group_names)} {node_type} groups listed in {fpath}')
        if self.useNeo4j():
          parents = self.neo4j.get_nodes(node_type,'Name',list(group_names))
          self.__Ontology__.update(self.neo4j._load_children_(parents))
        else:
          oql = 'SELECT Entity WHERE Name = ({props})'
          new_session = self._clone_session() #to avoid mutating self.Graph in parallel calculations
          request_name = f'Loading ontology {len(group_names)} groups from file {ontology_file}'
          ontology_graph = new_session.__iterate__(oql,group_names,request_name,step=500)
          if isinstance(ontology_graph,ResnetGraph):
            ontology_groups = ontology_graph._get_nodes()
            _, parents_with_children = new_session.load_children4(ontology_groups, max_threads=50)
            self.__Ontology__.update(parents_with_children)
          new_session.close_connection()
           
    return self.__Ontology__


  def add_ontology_df(self,_4df:df,add2report=True):
      '''
      input:
          df with name = for_df_name must be in self.report_pandas
      output:
          Adds worksheet to "report_pandas" containing statistics for ontology groups listed in 'ElsevierAPI/ResnetAPI/ontology/disease4ontology_analysis.txt'
      '''
      self._load_ontology()

      try:
          scored_disease_names = _4df[self.__resnet_name__].to_list()
      except KeyError:
          scored_disease_names = _4df['Name'].to_list()  # df made by load_df does not have self.__resnet_name__

      all_diseases = self.Graph._psobjs_with('Disease','ObjTypeName')
      scored_diseases = [x for x in all_diseases if x.name() in scored_disease_names]
      grant_total = len(scored_diseases)
      rows = [['All diseases',grant_total,100]]
      all_scored_children_counter = set()
      for parent_disease in self.__Ontology__:
          parent_ontology = [parent_disease] + parent_disease.childs()
          scored_children_in_parent_ontology = set(parent_ontology).intersection(scored_diseases)
          group_stat = len(scored_children_in_parent_ontology)
          rows.append([parent_disease.name(),group_stat, 100*group_stat/grant_total])
          all_scored_children_counter.update(scored_children_in_parent_ontology)

      stats_colname = 'Number of indications'
      other_count = len(scored_diseases)-len(all_scored_children_counter)
      rows.append(['Other diseases',other_count, 100*other_count/grant_total])
      ontology_df = df.from_rows(rows,['Ontology category',stats_colname,'Percentage'])
      ontology_df._name_ = ONTOLOGY_ANALYSIS+'4'+ _4df._name_ if _4df._name_ else ONTOLOGY_ANALYSIS

      ontology_df = ontology_df.sortrows(by=stats_colname)
      print('Created "Ontology analysis" table')

      if add2report:
          self.add2report(ontology_df)
      return ontology_df


  def remove_high_level_entities(self,from_df:df, max_childs=MAX_CHILDS):
    return_df = from_df.dfcopy()
    if max_childs:
        return_df = df.from_pd(return_df[return_df.apply(lambda x: len(x[self.__temp_id_col__]) <= max_childs, axis=1)])
        print(f'{len(from_df)-len(return_df)} entities with > {max_childs-1} ontology children were removed from further calculations')
        return_df.copy_format(from_df)
        return_df._name_ = from_df._name_
    return return_df


  def add_relevant_concepts(self,to_df:df,name2concepts:dict[str,list[PSObject]],from_graph:ResnetGraph=None):
    '''
    input:
      name2concepts = {concept_name:[Concepts]}
    output:
      df with new column "Relavant concept_names" containing [Concepts] linked to row_entities in the self.Graph
    '''
    def _relevant_concept_colname(concept_name:str):
      return 'Relevant ' + concept_name
    
    myG = from_graph if from_graph else self.Graph
    my_df = to_df.dfcopy()
    idtype = self.idtype()

    for name, concepts in name2concepts.items():
      column = _relevant_concept_colname(name)
      print(f'Adding {column} column to {to_df._name_}')
      for i in my_df.index:
        ids = list(to_df.at[i,self.__temp_id_col__])
        row_uids = list(myG.dbid2uid(ids,idtype).values())
        strs4cell = []
        for concept in concepts:
          concep_rels = list(myG.get_rels_between(row_uids, [concept.uid()]))
          if concep_rels:
            concep_rels.sort(key=lambda x: 
                             (
                              x.get_prop(REGULATOR_WEIGHT,if_missing_return=1),
                              x.count_refs()
                              ),
                              reverse=True)
            e2i_refs = set()
            [e2i_refs.update(rel.refs()) for rel in concep_rels]
            strs4cell.append(f'{concept.name()} ({str(len(e2i_refs))})')
        my_df.loc[i,column] = ','.join(strs4cell)
    return my_df


##################  PRINT ############################## PRINT ##################################
  def add_params2df(self):
    internal_params  = {"skip","debug", "consistency_correction4target_rank",
      "add_bibliography","strict_mode", "target_types","max_ontology_parent",
      "init_refstat","max_childs","add_closeness",'propagate_target_state_in_model',
      'DTfromDB','BBBP','use_in_children','ontology_file'
                        }
    rows = []
    for type, name2weights in self.params.items():
      if name2weights:
        if type not in internal_params:
          if isinstance(name2weights, dict):
            [rows.append([type,name,value]) for name,value in name2weights.items()]
          else:
            rows.append([type,type,name2weights])
        
    param_df = df.from_rows(rows,['Parameter','Concept','Value'])
    param_df = param_df.sortrows(by=['Value','Parameter','Concept'])
    param_df._name_ = INPUT_WOKSHEET
    param_df.tab_format['tab_color'] = 'yellow'
    self.add2report(param_df)


  def add_infodf(self):
     self.add2report(df.info_df())


  def cleand_df(self,report_df:df):
    df2print = report_df.dfcopy()
    report_df_columns = df2print.columns.to_list()
    my_drop = [c for c in self.columns2drop if c in report_df_columns]
    clean_pd = df2print.drop(columns=my_drop)
    clean_df = df.from_pd(clean_pd.loc[(clean_pd!=0).any(axis=1)].fillna('')) # removes rows with all zeros
    clean_df_columns = clean_df.columns.to_list()

    stop_index = clean_df_columns.index('Name') + 1 if 'Name' in clean_df_columns else 1
    new_column_order = clean_df_columns[:stop_index]

    hyperlinked_cols = [c for c in clean_df_columns if c.startswith(REFCOUNT_COLUMN)]
    if hasattr(self,'RefStats'):
      [hyperlinked_cols.append(c) for c in self.RefStats.refcols 
       if c in clean_df_columns and c not in new_column_order]
    clean_df.set_hyperlink_color(hyperlinked_cols)

    centered_columns = hyperlinked_cols
    if CHILDREN_COUNT in clean_df_columns: 
      centered_columns.append(CHILDREN_COUNT)
    [clean_df.add_column_format(c,'align','center') for c in centered_columns]
    new_column_order += centered_columns

    ranked_cols = [c for c in clean_df_columns if 'rank' in c.lower() and c not in new_column_order]
    new_column_order += ranked_cols
    
    other_columns = [c for c in clean_df_columns if c not in new_column_order]
    new_column_order += other_columns

    clean_df = clean_df.dfcopy(new_column_order)
    clean_df.copy_format(report_df)
    
    if any(s in report_df._name_ for s in {'ibliography','nippets','refs','snpt',INPUT_WOKSHEET,PHENOTYPE_WORKSHEET}):
      clean_df.make_header_horizontal()
    else:
      clean_df.make_header_vertical()
    return clean_df


  def add2writer(self,writer:pd.ExcelWriter,ws_prefix='',df_names:list[str]=[]):
    '''
    input:
      self.columns2drop
    '''
    if df_names:
      my_pandas = {k:self.report_pandas[k] for k in df_names if k in self.report_pandas}
    else:
      my_pandas = self.report_pandas

    def truncate_ws(ws_name:str):
      return ws_prefix+'_'+ws_name[:MAX_TAB_LENGTH-len(ws_prefix)-1] if ws_prefix else ws_name[:MAX_TAB_LENGTH]
  
    default_ws_order = [INPUT_WOKSHEET,PHENOTYPE_WORKSHEET,'Drugs']
    my_worksheets = list(my_pandas.keys())
    my_ws_order = [w for w in default_ws_order if w in my_worksheets]
    [my_ws_order.append(w) for w in my_worksheets if w not in my_ws_order]

    # gathering worksheet statistics into Info worksheet:
    info_rows = [['Worksheets in this file:','=SHEETS()']]
    [info_rows.append([f'=HYPERLINK("#\'{truncate_ws(w)}\'!A1", "{truncate_ws(w)}")',1]) for w in my_ws_order]
    infodf = df.from_rows(info_rows,['Info','Counts'],dfname=INFO_WORKSHEET)
    infodf.set_hyperlink_color(['Info'])
    self.add2report(infodf)
    my_ws_order = [INFO_WORKSHEET] + my_ws_order

    # writing all worksheets:
    for ws_name in my_ws_order:
      if ws_name in self.report_pandas:
        clean_df = self.cleand_df(self.report_pandas[ws_name])
        clean_df.df2excel(writer,sheet_name=truncate_ws(ws_name),float_format='%.3f')


  def addraw2writer(self, writer:pd.ExcelWriter, ws_prefix=''):
    if hasattr(self,'raw_data'):
      for ws_name, rawdf in self.raw_data.items():
        sh_name = ws_prefix+'_'+ws_name if ws_prefix else ws_name
        assert(isinstance(rawdf,df))
        raw_df = rawdf.dfcopy()
        raw_df.make_header_vertical()
        raw_df.df2excel(writer, sheet_name=sh_name[:30])


  def print_report(self, path2report:str, ws_prefix=''):
      writer = pd.ExcelWriter(path2report, engine='xlsxwriter')
      self.add2writer(writer,ws_prefix)
      writer.close()
      print('Report is in "%s"' % path2report)


  def print_rawdata(self, path2report:str, ws_prefix=''):
      writer = pd.ExcelWriter(path2report, engine='xlsxwriter')
      self.addraw2writer(writer,ws_prefix)
      writer.close()
      print('Raw data is in "%s"' % path2report)


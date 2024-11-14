import time, math, os, glob, json, threading
import networkx as nx
from zeep import exceptions
from xml.dom import minidom
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import timedelta
from collections import defaultdict

from .ZeepToNetworkx import PSNetworx, len
from .ResnetGraph import ResnetGraph,df,REFCOUNT,CHILDS,DBID,PSObject,PSRelation
from .PathwayStudioGOQL import OQL
from .Zeep2Experiment import Experiment
from ..ETM_API.references import PS_BIBLIO_PROPS,PS_SENTENCE_PROPS,PS_REFIID_TYPES,RELATION_PROPS,ALL_PSREL_PROPS
from ..ScopusAPI.scopus import loadCI, SCOPUS_CI
from ..utils import unpack,execution_time,execution_time2,load_api_config

TO_RETRIEVE = 'to_retrieve'
BELONGS2GROUPS = 'belongs2groups'
ALL_CHILDS = 0
# options for TO_RETRIEVE argument:
NO_REL_PROPERTIES = -2
CURRENT_SPECS = -1 # keep current retrieval properties
DO_NOT_CLONE = 0
DATABASE_REFCOUNT_ONLY = 1 #retrieves only RelationNumberOfReferences
REFERENCE_IDENTIFIERS = 2 # retrieves only PS_ID_TYPES to count reference in runtime
BIBLIO_PROPERTIES = 3 # retrieves only PS_BIBLIO_PROPS to generate biblio_str for reference output
SNIPPET_PROPERTIES = 4 # retrieves properties required to dump supporting references using self.Graph.snippets2df()
ONLY_REL_PROPERTIES = 5
ALL_PROPERTIES = 10

NO_RNEF_REL_PROPS={'RelationNumberOfReferences','Name','URN'}

APISESSION_KWARGS = {'what2retrieve','connect2server','no_mess','data_dir',TO_RETRIEVE,
                  'use_cache','load_model','ent_props','rel_props'}
MAX_SESSIONS = 50 # by default sessions_max=200 in Oracle 
MAX_PAGE_THREADS = 80
MAX_OQLSTR_LEN = 65000 # string properties can form long oql queries exceeding 65000 chars limit
ONTOLOGY_CACHE = os.path.join(os.getcwd(),'ENTELLECT_API/ElsevierAPI/ResnetAPI','__pscache__','ontology_cache.json')

class APISession(PSNetworx):
    '''
    manages graph retrieval from database and loading cache files
    '''
    pass
    entProps = ['Name'] # Name is highly recommended to make sense of the data
    relProps = ['URN'] # URNs are used as MultiDiGraph keys in self.Graph
    ResultPos = 0
    PageSize = 1000
    reference_cache_size = 1000000 # max number of reference allowed in self.Graph. Clears self.Graph if exceeded
    resnet_size = 1000 # number of <node><control> sections in RNEF dump
    max_rnef_size = 100000000 # max size of RNEF XML dump file. If dump file exceeds max_file_size new file is opened with index++
    max_sessions = MAX_SESSIONS
    data_dir = ''
    sep = '\t'
    dump_oql_queries = False
    max_threads4ontology = 50

######################################  CONFIGURATION  ######################################
    def __init__(self,*args,**kwargs):
        '''
        Input
        -----
        APIconfig = args[0], defaults to "ElsevierAPI/APIconfig.json"\n
        what2retrieve - defaults NO_REL_PROPERTIES, other options -\n
        [DATABASE_REFCOUNT_ONLY,REFERENCE_IDENTIFIERS,BIBLIO_PROPERTIES,SNIPPET_PROPERTIES,ONLY_REL_PROPERTIES,ALL_PROPERTIES]
        no_mess - default True, if False your script becomes more verbose
        connect2server - default True, set to False to run script using data in __pscache__ files instead of database
        '''
        self.GOQLquery = str()
        self.DumpFiles = []

        my_kwargs = {'what2retrieve':NO_REL_PROPERTIES,
                            'ent_props' : ['Name'],
                            'rel_props' : ['URN'],
                            'data_dir' : '',
                            'use_cache' : False,
                            'oql_queries' : [],
                            'add2self':True
                            }

        ent_props = kwargs.pop('ent_props',[])
        rel_props = kwargs.pop('rel_props',[])
        my_kwargs.update(kwargs)
        super().__init__(*args, **my_kwargs)
        
        self.set_dir(my_kwargs.get('data_dir',''))
        self.use_cache = my_kwargs.get('use_cache',False)# if True signals to use graph data from cache files instead of retrieving data from database using GOQL queries 

        self.__retrieve(my_kwargs['what2retrieve']) #__retrieve overides self.relProps, self.entProps
        # properties have to updated after self.__retrieve
        self.add_ent_props(ent_props)
        self.add_rel_props(rel_props)

        self.add2self = my_kwargs.get('add2self',True) # if False will not add new graph to self.Graph
        self.print_rel21row = False
        self.getLinks = True
        self.__IsOn1st_page = True # to print header in dump file
        self.dump_folder = str()
        self.my_oql_queries = list(my_kwargs.get('oql_queries',[])) # [(oql_query,request_name),...] list of GOQL queries tuples for fetching data from database
        self.ResultPos = int()
        self.ResultSize = int()
        self.ontology_cache = dict() # {urn:[urns]}

    @staticmethod
    def _what2retrieve(what2retrieve:int):
        if what2retrieve == NO_REL_PROPERTIES:
            return ['URN']
        elif what2retrieve == DATABASE_REFCOUNT_ONLY:
            return ['URN',REFCOUNT]
        elif what2retrieve == REFERENCE_IDENTIFIERS:
            return ['URN']+PS_REFIID_TYPES
        elif what2retrieve == BIBLIO_PROPERTIES:
            return ['URN']+PS_REFIID_TYPES+list(PS_BIBLIO_PROPS)
        elif what2retrieve == SNIPPET_PROPERTIES:
            return ['URN']+PS_REFIID_TYPES+list(PS_BIBLIO_PROPS)+PS_SENTENCE_PROPS
        elif what2retrieve == ONLY_REL_PROPERTIES:
            return ['URN']+list(RELATION_PROPS)
        elif what2retrieve == ALL_PROPERTIES:
            return ['URN']+list(ALL_PSREL_PROPS)
        return []


    def __retrieve(self,what2retrieve=CURRENT_SPECS):
        if isinstance(what2retrieve,int):
            retrieve_props = self._what2retrieve(what2retrieve)
            if retrieve_props:
                self.relProps = retrieve_props
            else:
                return # do nothing, keep CURRENT_SPECS
        else:
            print(f'Invalid what2retrieve: {what2retrieve}')
             

    @staticmethod
    def split_kwargs(kwargs:dict):
        '''
        Returns
        -------
        session_kwargs,other_parameters
        '''
        session_kwargs = dict()
        other_parameters = dict()
        [(other_parameters,session_kwargs)[k in APISESSION_KWARGS].update({k:v}) for k,v in kwargs.items()]
        return session_kwargs,other_parameters


    def _clone_session(self,**kwargs):
        """
        kwargs:
            what2retrieve - desired relation properties in clone session. Defaults to CURRENT_SPECS
            load_model - defaults to False to avoid database model reloading
            copy_graph - copy self.Graph. Defaults to False
        """
        my_kwargs = dict(kwargs)
        my_kwargs['what2retrieve'] = my_kwargs.pop(TO_RETRIEVE,CURRENT_SPECS)
        my_kwargs['load_model'] = False
        if my_kwargs['what2retrieve'] == CURRENT_SPECS:  
            new_session = APISession(self.APIconfig,**my_kwargs)
            new_session.entProps = list(self.entProps)
            new_session.relProps = list(self.relProps)
        else:
            new_session = APISession(self.APIconfig,**my_kwargs)
            new_session.entProps = list(self.entProps)
        
        new_session.data_dir = self.data_dir
        new_session.dump_folder = self.dump_folder
        new_session.no_mess = my_kwargs.get('no_mess',self.no_mess)
        new_session.max_threads4ontology = self.max_threads4ontology
        new_session.max_sessions = self.max_sessions
        new_session.IdtoObjectType = dict(self.IdtoObjectType)
        new_session.IdToPropType = dict(self.IdToPropType)
        new_session.RNEFnameToPropType = dict(self.RNEFnameToPropType)

        copy_graph = my_kwargs.pop('copy_graph',False)
        if copy_graph:
            new_session.Graph = self.Graph
            new_session.dbid2relation = dict(self.dbid2relation)

        return new_session


    def __set_get_links(self,oql_query=''):
        my_query = oql_query if oql_query else self.GOQLquery
        return my_query[7:15] == 'Relation'


    def __replace_goql(self, oql_query:str):
        self.GOQLquery = oql_query
        self.getLinks = self.__set_get_links()
        self.ResultRef = None


    def __add2self(self,other:'APISession'):
        '''
        needs relations to have dbid
        works only with another graph from database
        '''
        if self.add2self:
            self.Graph = self.Graph.compose(other.Graph)
            self.dbid2relation.update(other.dbid2relation)


    def set_dir(self,indir:str):
        self.data_dir = os.path.join(indir, '')


    def start_download_from(self,result_position:int):
        self.ResultPos = result_position


    def add_rel_props(self, add_props:list):
        self.relProps = self.relProps+[i for i in add_props if i not in self.relProps]
        my_sent_props = set(PS_SENTENCE_PROPS).intersection(add_props)
        if my_sent_props:
            ms = 5
            if self.max_sessions > ms:
                print(f'{type(self).__name__} decreases multithreading to {ms} \
to retreive {my_sent_props} properties')
                self.max_sessions = ms
            if 'TextRef' not in self.relProps:
                self.relProps.append('TextRef')
            

    def add_ent_props(self, add_props: list):
        self.entProps = self.entProps+[i for i in add_props if i not in self.entProps]

    
    def add_dump_file(self, dump_fname, replace_main_dump=False):
        if replace_main_dump:
            self.DumpFiles = []
        self.DumpFiles.append(dump_fname)

#########################  MULTITHREADED PROCESSING, RETRIEVAL   ##################################
    def __init_session(self, first_iteration=0,max_result=0, request_name=''):
        '''
        Return
        ------
        empty ResnetGraph if number of results > max_result
        '''
        self.__set_get_links()
        obj_props = self.relProps if self.getLinks else self.entProps
        page_size = max_result if max_result else self.PageSize
        zeep_data, (self.ResultRef, self.ResultSize, self.ResultPos) = self.init_session(self.GOQLquery,
                                                                                        page_size,obj_props,
                                                                                 getLinks = self.getLinks)
        
        if max_result and self.ResultSize > max_result:
            return ResnetGraph()
    
        if first_iteration > 0:
            self.ResultPos = first_iteration
            print (f'Resuming retrieval from {first_iteration} position' ,flush=True)
            return self.__nextpage__(self.ResultPos)
        else:
            if request_name and not max_result:  # to suppress messages from load_children4
                print(f'query "{request_name}" found {self.ResultSize} results',flush=True)

        if not isinstance(zeep_data, type(None)):
            if self.getLinks:
                obj_dbids = list(set([x['EntityId'] for x in zeep_data.Links.Link]))
                zeep_objects = self.get_object_properties(obj_dbids, self.entProps)
                return self._load_graph(zeep_data, zeep_objects,self.add2self)
            else:
                return self._load_graph(None, zeep_data,self.add2self)
        else:
            return ResnetGraph()

        
    def __nextpage__(self,result_pos:int):
        '''
        Input
        -----
        result_pos - resut position to begin download. Provide explicitly to initiate multithreading
        '''
        with threading.Lock():
          # current_pos does not change during multithreading initiation!!!! 
          if int(self.ResultPos) < int(self.ResultSize):
              obj_props = self.relProps if self.getLinks else self.entProps
              zeep_data, self.ResultSize, current_pos = self.get_session_page(self.ResultRef, result_pos, self.PageSize,
                                                                        self.ResultSize,obj_props,getLinks=self.getLinks)                                                          
              if not isinstance(zeep_data, type(None)):
                  if self.getLinks and len(zeep_data.Links.Link) > 0:
                      obj_dbids = list(set([x['EntityId'] for x in zeep_data.Links.Link]))
                      zeep_objects = self.get_object_properties(obj_dbids, self.entProps)
                      return self._load_graph(zeep_data, zeep_objects,self.add2self)
                  else:
                      return self._load_graph(None, zeep_data,self.add2self)
              else:
                  return ResnetGraph()
          else: return ResnetGraph()


    def __thread__(self,pages:int,process_name='oql_results'):
        max_workers = pages if pages < MAX_PAGE_THREADS else MAX_PAGE_THREADS
        if not self.no_mess:
            print(f'Retrieval starts from {self.ResultPos} result')
            remaining_retrieval = self.ResultSize - self.ResultPos
            next_download_size = min(pages*self.PageSize, remaining_retrieval)
            print(f'Begin retrieving next {next_download_size} results in {max_workers} threads',flush=True)
        
        entire_graph = ResnetGraph()
        result_pos = self.ResultPos
        with ThreadPoolExecutor(max_workers, thread_name_prefix=process_name) as e:
            futures = list()
            for i in range(0,pages):
                futures.append(e.submit(self.__nextpage__,result_pos))
                result_pos += self.PageSize
            
            for future in as_completed(futures):
                try:
                    entire_graph = entire_graph.compose(future.result())
                except exceptions.TransportError:
                    raise exceptions.TransportError
            e.shutdown()
        
        self.ResultPos = result_pos
        return entire_graph


    def process_oql(self,oql_query:str,request_name='',max_result=0, debug=False) -> ResnetGraph|int:
        '''
        Return
        ------
        if max_result is not 0 and number of results exceeds max_result returns int = self.ResultSize\n
        otherwise returns ResnetGraph with query results
        '''
        with threading.Lock():
          self.__replace_goql(oql_query)
          start_time = time.time()
          return_type = 'relations' if self.getLinks else 'entities'
          entire_graph = self.__init_session(max_result=max_result,request_name=request_name)
          if debug: return entire_graph

          if max_result and self.ResultSize > max_result:
              return self.ResultSize

          if entire_graph:
              pages = int(self.ResultSize / self.PageSize)
              if pages:
                  my_request_name = request_name if request_name else self.GOQLquery[:100]+'...'
                  if not self.no_mess:
                      print('\n\"%s\"\nrequest found %d %s.\n%d is retrieved. Remaining %d results will be retrieved in %d iterations' % 
                  (my_request_name,self.ResultSize,return_type,self.ResultPos,(self.ResultSize-self.PageSize),pages))
                  try:
                      iterations_graph = self.__thread__(pages,process_name=request_name)
                      entire_graph = entire_graph.compose(iterations_graph)
                  except exceptions.TransportError:
                      raise exceptions.TransportError('Table lock detected!!! Aborting operation!!!')

                  if not self.no_mess:
                      print('"%s"\nretrieved %d nodes and %d edges in %s by %d parallel iterations' % 
                      (my_request_name, entire_graph.number_of_nodes(), entire_graph.number_of_edges(),
                              execution_time(start_time), pages+1),flush=True)

          self.ResultRef = ''
          self.ResultPos = 0
          self.ResultSize = 0
          self.__IsOn1st_page = True
          return entire_graph


    def __annotate_dbid_graph__(self,id_only_graph:ResnetGraph,use_cache=True,request_name='',step=1000):
        '''
        loads
        -----
        new relations to self.Graph by ids in id_only_graph
        '''
        print(f'Retrieving annotations for graph with {len(id_only_graph)} nodes and {id_only_graph.number_of_edges()} relations')
        req_name = f'Annotations 4 "{request_name}"'
        need_reldbids = id_only_graph._relations_dbids()
        need_nodedbids = id_only_graph.dbids4nodes()
        
        return_subgraph = ResnetGraph()
        relations_dbids2retreive = set(need_reldbids)
        nodes_dbids2retreive = set(need_nodedbids)
        
        if self.Graph.number_of_nodes() > 0: # check because downstream functions do not like empty self.Graph
            rels4subgraph = [rel for dbid,rel in self.dbid2relation.items() if dbid in need_reldbids]
            return_subgraph = self.Graph.subgraph_by_rels(rels4subgraph) 
            nodes_in_cache = self.Graph.psobj_with_dbids(set(need_nodedbids))
            if nodes_in_cache:
                return_subgraph.add_psobjs(set(nodes_in_cache))
            
            if use_cache:
                relations_dbids2retreive = need_reldbids.difference(self.dbid2relation.keys())
                nodes_dbids2retreive = set(need_nodedbids).difference(set(self.Graph.dbids4nodes()))

        if id_only_graph.number_of_nodes() > 0:
            exist_nodes = len(need_nodedbids)-len(nodes_dbids2retreive)
            exist_rels = len(need_reldbids)-len(relations_dbids2retreive)
            print(f'{exist_nodes} nodes and {exist_rels} relations were downloaded from database previously')
            print('%d nodes and %d relations will be loaded from database' % 
                        (len(nodes_dbids2retreive),len(relations_dbids2retreive)))
        
        add2return = ResnetGraph()
        if relations_dbids2retreive:
            rel_query = 'SELECT Relation WHERE id = ({ids})'
            rels_add2return = self.__iterate__(rel_query,relations_dbids2retreive,req_name,step=step)
            add2return = rels_add2return
            nodes_dbids2retreive = nodes_dbids2retreive.difference(add2return.dbids4nodes())
            
        if nodes_dbids2retreive:
            entity_query = 'SELECT Entity WHERE id = ({ids})'
            nodes_add2return = self.__iterate__(entity_query,nodes_dbids2retreive,req_name,step=step)
            add2return = add2return.compose(nodes_add2return)
        
        return_subgraph = return_subgraph.compose(add2return)
        assert(return_subgraph.number_of_nodes() == len(need_nodedbids))
       # assert(return_subgraph.number_of_edges() >= len(relation_dbids2return))
       # privately owned relations are not returned by API. 
       # therefore, return_subgraph.number_of_edges() may have less relations than "relation_dbids2return"
       
       # return_subgraph may also contain more relations than "relation_dbids2return" 
       # due to duplication of non-directional relations into both directions
        return return_subgraph


    def choose_process(self,oql:str,request_name:str,download=False,lock=None):
        '''
        if download does not close the last batch file\n
        use self.close_rnef_dump() to close batch
        '''
        if download:
            self.download_oql(oql,request_name,lock=lock,close_batch=False)
        else:
            return self.process_oql(oql,request_name)


    def __iterate__(self,oql_query:str,dbids:set,request_name='',step=1000,download=False):
        '''
        Input
        -----
        oql_query MUST contain string placeholder called {ids} to iterate dbids\n
        oql_query MUST contain string placeholder called {props} to iterate other database properties
        "step" - controls duration of request to avoid timeout during multithreading

        uses self.entProp and self.relProp to retrieve properties\n
        use self.add2self and self.merge2self to control caching
 
        if download does not close the last batch file\n
        use self.close_rnef_dump() to close batch
        '''
        if not dbids: return ResnetGraph()
        if oql_query.find('{ids',20) > 0:
            my_oql = oql_query
            iteration_size =  min(step,1000)
            def join_list (l:list):
                return ','.join(list(map(str,l)))
        else:
            my_oql = oql_query.replace('{props','{ids')
            # string properties can form long oql queries exceeding 65500 chars limit
            max_id_len = max([len(str(s)) for s in dbids])
            iteration_size = min(step, 1000, int(MAX_OQLSTR_LEN/max_id_len))
            def join_list (l:list):
                return OQL.join_with_quotes(l)

        number_of_iterations = math.ceil(len(dbids)/iteration_size)  
        # new sessions must be created to avoid ResultRef pointer overwriting
        entire_graph = ResnetGraph()
        thread_name = f'Iterate_{len(dbids)}_ids_by_{iteration_size}'
        print(f'Retrieval of {len(dbids)} objects will be done in {self.max_sessions} parallel sessions' )
        dbids_list = list(dbids)
        lock = threading.Lock() if download else None
        with ThreadPoolExecutor(max_workers=self.max_sessions, thread_name_prefix=thread_name) as e:
            future_sessions = list()
            futures = list()
            for i in range(0,len(dbids_list), iteration_size):
                iter_dbids = dbids_list[i:i+iteration_size]
                oql_query_with_dbids = my_oql.format(ids=join_list(iter_dbids))
                iter_name = f'Iteration #{int(i/iteration_size)+1} in {number_of_iterations} for "{request_name}" retrieving {len(iter_dbids)} ids'
                if self.dump_oql_queries:
                    self.my_oql_queries.append((oql_query_with_dbids,iter_name))   
                future_session = self._clone_session()
                futures.append(e.submit(future_session.choose_process,oql_query_with_dbids,iter_name,download,lock))
                future_sessions.append(future_session)
            
            if not download:
                for future in as_completed(futures):
                    entire_graph = entire_graph.compose(future.result())
            
            for session in future_sessions:
                self.__add2self(session)
                session.close_connection()
            e.shutdown()
        
        if self.dump_oql_queries:
            with open(self.data_dir+"iterations_oql.json","w") as o:
                json.dump(self.my_oql_queries,o)
        return entire_graph


    def __iterate_oql__(self,oql_query:str,ids_or_props:set,request_name='',step=1000,download=False):
        '''
        Input
        -----
        oql_query MUST contain string placeholder called {ids} to iterate dbids\n
        oql_query MUST contain string placeholder called {props} to iterate other database properties
        "step" - controls duration of request to avoid timeout during multithreading

        Returns
        -------
        ResnetGraph containing only dbids and no other properties
        function does not add retrieved graph to self.Graph\n
        and must be used togeter with APISession::__annotate_dbid_graph__ \n
        '''
        req_name = 'Retrieve DB identifiers 4 "'+request_name+'"'
        my_ent_props = list(self.entProps)
        my_rel_props = list(self.relProps)

        self.relProps.clear()
        self.entProps.clear()
        old_add2self = self.add2self
        self.add2self = False

        dbidonly_graph = self.__iterate__(oql_query,ids_or_props,req_name,step,download)

        self.entProps = my_ent_props
        self.relProps = my_rel_props
        self.add2self = old_add2self 
        return dbidonly_graph
    

    def __iterate2__(self,oql_query:str,ids1:list,ids2:list,
                     request_name='iteration1',step=500,download=False):
        '''
        Input
        -----
        oql_query MUST contain string placeholder called {ids} to iterate dbids\n
        oql_query MUST contain string placeholder called {props} to iterate other database properties
        '''
        
        if oql_query.find('{ids',20) > 0:
            id_oql = oql_query
            hint = r'{ids}'
            iteration1step = min(step, 1000) # cannot exceed 1000
            iteration2step = min(2*iteration1step, 1000) # cannot exceed 1000
            def join_list (l:list):
                return ','.join(list(map(str,l)))
        else:
            id_oql = oql_query.replace('{props','{ids')
            hint = r'{props}'
            # string properties can form long oql queries exceeding 65500 chars limit
            max_id_len = max(len(str(s)) for s in ids1)
            iteration1step = min(step, 1000, int(MAX_OQLSTR_LEN/max_id_len))
            iteration2step = min(2*iteration1step, 1000)
            def join_list (l:list):
                return OQL.join_with_quotes(l)
                             
        number_of_iterations = math.ceil(len(ids1)/iteration1step) * math.ceil(len(ids2)/iteration2step)
        if not number_of_iterations: return ResnetGraph()
        print(f'Iterating {len(ids1)} entities with {len(ids2)} entities for {request_name}')
        if number_of_iterations > 2:
            print(f'Query will be executed in {number_of_iterations} iterations by {self.max_sessions} sessions')

        entire_graph = ResnetGraph()
        if len(ids1) <= len(ids2):
            for i1 in range(0,len(ids1), iteration1step):
                iter_ids1 = ids1[i1:i1+iteration1step]
                iter_oql_query = id_oql.format(ids1=join_list(iter_ids1),ids2=hint)
                iter_graph = self.__iterate__(iter_oql_query,set(ids2),'iteration2',iteration2step,download)
                entire_graph = entire_graph.compose(iter_graph)
                print(f'Processed {min([i1+iteration1step,len(ids1)])} in {len(ids1)} identifiers against {len(ids2)} identifiers')
        else:
            for i2 in range(0,len(ids2), iteration1step):
                iter_ids2 = ids2[i2:i2+iteration1step]
                iter_oql_query = id_oql.format(ids1=hint,ids2=join_list(iter_ids2))
                iter_graph = self.__iterate__(iter_oql_query,set(ids1),'iteration2',iteration2step,download)
                entire_graph = entire_graph.compose(iter_graph)
                print(f'Processed {min([i2+iteration1step,len(ids2)])} in {len(ids2)} identifiers against {len(ids1)} identifiers')
            
        return entire_graph
    

    def __iterate_oql2__(self,oql_query:str, ids_or_props1:set, ids_or_props2:set, request_name='', step=500):
        '''
        Input
        -----
        oql_query MUST contain 2 string placeholders called {ids1} and {ids2} to iterate dbids1/2 as dbids
        oql_query MUST contain 2 string placeholders called {props1} and {props2} to iterate other database properties

        Return
        ------
        ResnetGraph containing only dbids and no other properties
        function does not add retrieved graph to self.Graph\n
        and must be used togeter with APISession::__annotate_dbid_graph__ \n
        '''
        my_ent_props = list(self.entProps)
        my_rel_props = list(self.relProps)

        self.relProps.clear()
        self.entProps.clear()
        old_add2self = self.add2self
        self.add2self = False

        dbidonly_graph = self.__iterate2__(oql_query,list(ids_or_props1),list(ids_or_props2),request_name,step)

        self.entProps = my_ent_props
        self.relProps = my_rel_props
        self.add2self = old_add2self
        return dbidonly_graph

 
    def iterate_oql(self,oql_query:str,search_values:set[str]|set[int],use_cache=True,request_name='',step=1000):
        """
        # oql_query MUST contain string placeholder called {ids} if iterable id_set contains dbids\n
        # oql_query MUST contain string placeholder called {props} if iterable id_set contain property values other than database id
        """
        print('Processing "%s" request\n' % request_name)
        dbid_only_graph = self.__iterate_oql__(oql_query,search_values,request_name,step)
        self.close_connection()
        if dbid_only_graph:
            annotated_graph = self.__annotate_dbid_graph__(dbid_only_graph,use_cache,request_name,step)
            return annotated_graph
        else:
            return ResnetGraph()
   

    def iterate_download(self,oql_query:str,search_values:set[str]|set[int],request_name='',step=1000):
        """
        # oql_query MUST contain string placeholder called {ids} if iterable id_set contains dbids\n
        # oql_query MUST contain string placeholder called {props} if iterable id_set contain property values other than database id
        """
        return self.__iterate_oql__(oql_query,search_values,request_name,step,download=True)


    def iterate_oql2(self,oql_query:str,search_values1:set[str]|set[int],search_values2:set[str]|set[int]
                     ,use_cache=True,request_name='',step=500)->ResnetGraph:
        """
        # oql_query MUST contain 2 string placeholders called {ids1} and {ids2}
        # oql_query MUST contain 2 string placeholder called {props1} and {props2}\n
        if iterable id_sets contain property values other than database id
        """
        print(f'Processing "{request_name}" request')
        dbid_only_graph = self.__iterate_oql2__(oql_query,search_values1,search_values2,request_name,step)
        self.close_connection()
        if dbid_only_graph:
            #sleep(10)
            print('Will retrieve annotated graph with %d entities and %d relations' 
                        %(dbid_only_graph.number_of_nodes(),dbid_only_graph.number_of_edges()))
            old_max = self.max_sessions
            self.max_sessions = 10 # reducing self.max_sessions to allow sessions from __iterate_oql2__ to close
            annotated_graph = self.__annotate_dbid_graph__(dbid_only_graph,use_cache,request_name)
            self.max_sessions = old_max
            return annotated_graph
        else:
            return ResnetGraph()
     

    def __iterate_oqls__(self,oqls:list[str]):
        max_workers = len(oqls) if len(oqls) < MAX_SESSIONS else MAX_SESSIONS
        accumulate_graph = ResnetGraph()
        with ThreadPoolExecutor(max_workers=max_workers,thread_name_prefix=f'Iterate {len(oqls)} OQLs') as e:
            futures = list()
            for i, oql in enumerate(oqls):
                job_name = f'OQL#{i}'
                new_session = self._clone_session() # need to clone since self.dbid2relation mutates
                new_session.add2self = False
                futures.append(e.submit(new_session.process_oql,oql,request_name=job_name))

            for f in as_completed(futures):
                accumulate_graph = accumulate_graph.compose(f.result())
            e.shutdown()
        return accumulate_graph
            

    def iterate_oql3(self,reg2targ:list,oql_template:str):
        '''
        Input
        -----
        reg2targ - [(PSObject.PSObject)..], where 1st element is regulator, 2nd is target
        oql_template - must have two placeholders {reg_urn},{targ_urn}
        '''
        oql_queries = list()
        [oql_queries.append(oql_template.format(reg_urn=rt[0].urn(),targ_urn=rt[1].urn())) for rt in reg2targ]
        return self.__iterate_oqls__(oql_queries)


    def clear(self):
        self.Graph.clear_resnetgraph()
        self.dbid2relation.clear()


    def flush_dump_files(self):
        for f in self.DumpFiles:
            open(f, 'w').close()
            if not self.no_mess:
                print('File "%s" was cleared before processing' % f)


    @staticmethod
    def reopen(fname):
        open(fname, "w", encoding='utf-8').close()
        return open(fname, "a", encoding='utf-8')

    def get_result_size(self,oql_query:str):
        zeep_relations, (self.ResultRef, self.ResultSize, self.ResultPos) = self.init_session(oql_query, PageSize=1,
                                                                                              property_names=[])
        return self.ResultSize



###########################GET GET RETRIEVE GET RETRIEVE GET GET ################################################
    def __get_saved_results(self,psobjs:list[PSObject]):
        result_names = [o.name() for o in psobjs]
        oql_query = f'SELECT Result WHERE Name = ({OQL.join_with_quotes(result_names)})'
        results = self.load_graph_from_oql(oql_query,self.relProps,
                                                self.entProps,get_links=False)._get_nodes()
        
        def __has_relations(result):
            # undestands if search results contain entities or relations 
            # load  self.ResultRef
            self.ResultRef  = 'ID='+str(result.dbid())
            zeepdata,self.ResultSize, pos = self.get_session_page(self.ResultRef,0, 1,1,[],getLinks=True)
            if  not isinstance(zeepdata, type(None)):
                    first_obj = zeepdata.Objects.ObjectRef[0]
            else:
                return False
            
            return True if first_obj['ObjClassId'] == 3 else False

        self.PageSize  = 10000
        result_graph = ResnetGraph()
        for result in results:
            if __has_relations(result):
                self.getLinks = True
            else:
                self.getLinks = False                 
           
            result_graph = self.__nextpage__(0)
            number_of_iterations = int(self.ResultSize/self.PageSize)

            if number_of_iterations:
                process_name = f'Retrieve {result.name()}'
                with ThreadPoolExecutor(number_of_iterations, thread_name_prefix=process_name) as e:
                    futures = list()
                    result_pos = self.PageSize
                    for i in range(0,number_of_iterations):          
                        futures.append(e.submit(self.__nextpage__,result_pos))
                        result_pos += self.PageSize
                
                    for future in as_completed(futures):
                        try:
                            result_graph = result_graph.compose(future.result())
                        except exceptions.TransportError:
                            raise exceptions.TransportError
                    e.shutdown()
                
            print (f'Retreived {len(result_graph)} objects and {result_graph.number_of_edges()} from {result.name()}')
        return result_graph
    

    def __get_group_graph(self,psobjs:list):
        dbids = [str(o.dbid()) for o in psobjs]
        dbids_str = ",".join(dbids)
        oql_query = f'SELECT Entity WHERE MemberOf (SELECT Group WHERE Id = ({dbids_str}))'
        return self.load_graph_from_oql(oql_query,entity_props=self.entProps,get_links=False)
    

    def graph4obj(self,folder_obj:PSObject):
        '''
        Input
        -----
        obj must have type 'Pathway', 'Group', 'attributesSearch'
        '''
        if folder_obj.objtype() == 'Pathway':
            return self.pathway_components([folder_obj.dbid()],'id',self.relProps,self.entProps)
        elif folder_obj.objtype() == 'Group':
            return self.__get_group_graph([folder_obj])
        elif folder_obj.objtype() == 'attributesSearch':
            return self.__get_saved_results([folder_obj])
        else:
            print(f'Cannot return graph for {folder_obj.objtype()}')
            return ResnetGraph()


    def get_network(self,_4dbids:set,connect_by_reltypes:list=[],download=False):
        if download:
            oql_query = 'SELECT Relation WHERE NeighborOf (SELECT Entity WHERE id = ({ids1})) AND NeighborOf (SELECT Entity WHERE id = ({ids2}))'
            if connect_by_reltypes:
                oql_query += f' AND objectType = ({str(",".join(connect_by_reltypes))})'

            splitter = list() #holds lists of ids for splitting
            splitter.append(list(_4dbids)) # inititalizing with orginal list
            number_of_splits = int(math.log2(len(_4dbids)))
            print(f'Will load network of {len(_4dbids)} nodes in {number_of_splits} iterations using {2**number_of_splits} splits')
            retreival_start = time.time()

            step = 1000
            number_step_iterations = ((len(_4dbids)/step)**2)/self.max_sessions # estimating number of iterations in __iterate2__()
            print(f'Estimated number of iterations with step {step}: {number_step_iterations}')

            for s in range(1, number_of_splits):
                splitter_splits = list()
                half = int(len(splitter[0]) / 2)
                #iter_start = time.time()
                ids1 = set()
                ids2 = set()
                for i, split in enumerate(splitter):
                    uq_list1 = split[:half]
                    uq_list2 = split[half:]
                    if len(splitter) <= number_step_iterations: # it is faster to download large splits individually
                        self.__iterate2__(oql_query,list(uq_list1),list(uq_list2),
                        f'Download network: split #{i+1} in {len(splitter)} in iteration #{s} out of {number_of_splits}',
                        download=True,step=step)
                    else:
                        ids1.update(uq_list1)
                        ids2.update(uq_list2)
                    splitter_splits.append(uq_list1)
                    splitter_splits.append(uq_list2)

                if len(splitter) > number_step_iterations: # it is faster to accumulate small splits and process them all at once
                    self.__iterate2__(oql_query,list(ids1),list(ids2),f'Download network: iteration #{s} out of {number_of_splits}',
                                      download=True,step=step)
                splitter = splitter_splits
                remaining_iterations = number_of_splits-s
                executiontime, remaintime = execution_time2(retreival_start,remaining_iterations,number_of_splits)
                print(f'Network retrieval iteration {s} out of {number_of_splits} was completed in {executiontime}')
                print(f'Estimated remaining time for network retrieval: {remaintime}\n\n')
                s += 1
            self.close_rnef_dump(self.dump_folder)
        else:
            return super().get_network(_4dbids,connect_by_reltypes,self.relProps,self.entProps)


    def connect_nodes(self,node_dbids1:set,node_dbids2:set,
                        by_relation_type=[],with_effect=[],in_direction='',
                        use_relation_cache=True, step=500,download=False):
        """
        Input
        -----
        in_direction must be '>' or '<'

        Returns
        -----
        ResnetGraph containing relations between node_dbids1 and node_dbids2
        """
        oql_query = r'SELECT Relation WHERE '
        if in_direction:
            if in_direction == '>':
                oql_query = oql_query + 'NeighborOf downstream (SELECT Entity WHERE id = ({ids1})) AND NeighborOf upstream (SELECT Entity WHERE id = ({ids2}))'
            elif in_direction == '<':
                oql_query = oql_query + 'NeighborOf upstream (SELECT Entity WHERE id = ({ids1})) AND NeighborOf downstream (SELECT Entity WHERE id = ({ids2}))'
            else: print('Direction sign is not recorgnized!')
        else:
            oql_query = oql_query + 'NeighborOf (SELECT Entity WHERE id = ({ids1})) AND NeighborOf (SELECT Entity WHERE id = ({ids2}))'

        if with_effect:
            effect_str = ','.join(with_effect)
            oql_query = oql_query + ' AND Effect = ('+effect_str+')'
        else:
            effect_str = 'any'

        if by_relation_type:
            reltype_str = ','.join(by_relation_type)
            oql_query = oql_query + ' AND objectType = ('+reltype_str+')'
        else:
            reltype_str = 'all'

        my_dbids2 = node_dbids2.copy()
        common_ids = node_dbids1.intersection(node_dbids2)
        if common_ids and not in_direction:
            print(f'Removed {len(common_ids)} common identifiers between two input lists from second input list to avoid false expand')
            my_dbids2 = node_dbids2.difference(node_dbids1)

        req_name = f'Connecting {len(node_dbids1)} and {len(my_dbids2)} entities by {reltype_str} relations and {effect_str} effects'
         # step is set 250 to mitigate possible timeout during multithreading
        if download:
           self.__iterate2__(f'{oql_query}',list(node_dbids1),list(my_dbids2),req_name,step,download=True)
           return ResnetGraph()
        else:
            return self.iterate_oql2(f'{oql_query}',node_dbids1,my_dbids2,use_relation_cache,request_name=req_name,step=step)
   

    def connect_entities(self,objs:list[PSObject],and_obj:list[PSObject],
                         by_relation_type=[],with_effect=[],in_direction=''):
        dbids1 = set(ResnetGraph.dbids(objs))
        dbids2 = set(ResnetGraph.dbids(and_obj))
        assert(dbids1)
        assert(dbids2)
        return self.connect_nodes(dbids1,dbids2,by_relation_type,with_effect,in_direction)
    

    def get_group_members(self, group_names:list):
        groups = OQL.join_with_quotes(group_names)
        oql_query = f'SELECT Entity WHERE MemberOf (SELECT Group WHERE Name = ({groups}))'
        req_name = 'Find members of groups: ' + ','.join(group_names)
        graph2return = self.process_oql(oql_query, request_name=req_name)
        if isinstance(graph2return,ResnetGraph):
            group_size = graph2return.number_of_nodes()
            if group_size:
                print(f'loaded {group_size} members from {group_names}')
            else:
                print(f'{group_names} groups are empty or do not exist in database')
        return graph2return
    

    def add_group_annotation(self,group_names:list,_2graph=ResnetGraph()):
        urns2values = PSObject()
        for group_name in group_names:
            group_graph = self.get_group_members([group_name])
            if isinstance(group_graph,ResnetGraph):
                group_members = group_graph._get_nodes()
                [urns2values[o.urn()].append(group_name) for o in group_members]

        if _2graph: # my_graph = _2graph if _2graph else self.Graph does not work here
            _2graph.set_node_annotation(urns2values,BELONGS2GROUPS)
        else:
            self.Graph.set_node_annotation(urns2values,BELONGS2GROUPS)
        return


    def map_props2objs(self,using_values:list,in_properties:list[str],
                       case_insensitive=False,only_objtypes:list=[])->tuple[dict[str,list[PSObject]],dict[int,list[str]]]:
        '''
        output:\n
        propval2objs = {prop_value:[PSObject]}, where prop_value is from 'using_values'\n
        objid2propval = {node_uid:[prop_values]}
        '''
        propval2objs,uid2propval = self.Graph.props2obj_dict(using_values, in_properties,case_insensitive)
        need_db_mapping = set(using_values).difference(propval2objs.keys())

        current_ent_props = self.entProps
        self.add_ent_props(in_properties)

        step = 1000
        iteration_counter = math.ceil(len(need_db_mapping) / step)
        print('Will use %d %s identifiers to find entities in %d iterations' % 
             (len(need_db_mapping), ','.join(in_properties), iteration_counter))

        prop_name_str,_ = OQL.get_search_strings(in_properties,using_values)
        oql = f'SELECT Entity WHERE ({prop_name_str})' + ' = ({props})'
        if only_objtypes:
            oql += f' AND objectType = ({only_objtypes})'
        request_name = f'Mapping {len(using_values)} {prop_name_str}'
        mapping_graph = self.iterate_oql(oql,set(using_values),request_name=request_name)

        p2obj,id2p = mapping_graph.props2obj_dict(list(need_db_mapping),in_properties,case_insensitive)
        propval2objs.update(p2obj)
        for id,prop_vals in id2p.items():
            try:
                mapped_values = set(uid2propval[id])
                mapped_values.update([str(x).lower() for x in prop_vals])
                uid2propval[id] = list(mapped_values)
            except KeyError:
                uid2propval[id] = prop_vals

        print("%d entities were mapped to %d attributes in %s properties" % 
                (len(uid2propval), len(using_values),','.join(in_properties)))

        self.entProps = current_ent_props
        return propval2objs, uid2propval

################################# ONTOLOGY  ONTOLOGY ############################################
    def __load_children(self,parent:PSObject,min_connectivity=0,depth=0,max_childs=ALL_CHILDS)->tuple[PSObject,list[PSObject]]:
        '''
        Input
        -----
        if max_childs=0 finds all children for parent in database

        Return
        ------
        Tuple: parent:PSObject, children:list[PSObject]
        '''
        parent_dbid = parent.dbid()
        if parent_dbid:
            query_ontology = OQL.get_childs([parent_dbid],['id'],depth=depth)
            if min_connectivity:
                query_ontology += f' AND Connectivity >= {min_connectivity}'

            my_session = self._clone_session(what2retrieve=NO_REL_PROPERTIES)
            #request_name = f'Retrieve ontology children for {parent.name()}'
            request_name = ''
            children_graph = my_session.process_oql(query_ontology,request_name,max_result=max_childs)
        #      if parent.name() == 'physical illness':
        #          print('')
            my_session.close_connection()
            if isinstance(children_graph,int):
                # fake list of empty children for consistency to enable downstream remove_high_level_entities()
                return parent, [PSObject()]*children_graph 
            else:
                return parent, children_graph._get_nodes()
        else:
            print(f'Parent {parent.name()} has no database identifier. Its children cannot be loaded')
            return parent,[]
        

    def _load_children4(self,parents:list[PSObject],min_connectivity=0,depth=0,
                       max_childs=ALL_CHILDS,max_threads=50)->tuple[set[PSObject],list[PSObject]]:
        '''
        Input:
            if "parents" is empty will use all nodes from self.Graph
            depth - ontology depth to get children from. If depth=0, all children from all depths are returned
            max_childs - maximum number of children allowed in parent - cutoff for high level ontology concepts
            min_connectivity - avoid from using.  min_connectivity > 0 slows down children retrival significantly

        Output:
            {PSObject},[PSObject] - set of all found children, list of all parents with children\n
            if max_childs > 0, parents with number of children exeeding max_childs are excluded
        
        Updates:
            parents in self.Graph with CHILDS properties as [PSObject] children
        '''
        process_start = time.time()
        print(f'Loading ontology children for {len(parents)} entities')
        # by default sessions_max=200 in Oracle
        thread_name_prefix = f'Childs 4 {len(parents)} prnts in {max_threads} thrds-'
        max_threaded_time = 0
        need_children = list()
        child_counter = set()
        iteration_counter = 0
        for p, parent in enumerate(parents):
            if CHILDS in self.Graph.nodes[parent.uid()]:
                child_counter.update(self.Graph.nodes[parent.uid()][CHILDS])
            else:
                need_children.append(parent)
                if len(need_children) >= max_threads or p == len(parents)-1:
                    max_workers = min(len(need_children),max_threads)
                    thread_start = time.time()
                    iteration_counter += 1
                    thread_name = thread_name_prefix+str(iteration_counter) + 'iter'
                    with ThreadPoolExecutor(max_workers=max_workers, thread_name_prefix=thread_name) as e:
                        futures = list()
                        [futures.append(e.submit(self.__load_children,nc,min_connectivity,depth,max_childs)) for nc in need_children]

                        for future in as_completed(futures):
                            f_parent, children = future.result()
                            # children has empty PSObjects if len(children) <= max_childs
                            nx.function.set_node_attributes(self.Graph, {f_parent.uid():{CHILDS:children}})
                            if max_childs == ALL_CHILDS or len(children) <= max_childs:
                                # case when children has real PSObjects
                                self.Graph.add_psobjs(children) # set_node_attributes() does not update nodes that do not exist in Graph
                                child_counter.update(children)
                    e.shutdown()

                    threaded_time = time.time() - thread_start
                    if threaded_time > max_threaded_time:
                        max_threaded_time = threaded_time
                        if max_threaded_time > 240:
                            print(f'Longest {max_threads}-threaded time {"{}".format(str(timedelta(seconds=max_threaded_time)))} \
is close to Apache default 5 minutes transaction timeout !!!')
                    if not self.no_mess:
                        print (f'{p+1} {len(parents)} parents were processed in {execution_time(process_start)}')
                    need_children.clear()
        
        print(f'{len(child_counter)} children for {len(parents)} parent entities were found in database \
in {execution_time(process_start)}')
        print(f'Longest {max_threads}-threaded time: {"{}".format(str(timedelta(seconds=max_threaded_time)))}')

        parent_with_children = self.Graph.psobjs_with([CHILDS])
        parent_with_children = [p for p in parent_with_children if p.childs() and p in parents]
        child_counter = {x for x in child_counter if x} #sometimes they are empty?
        return child_counter, parent_with_children


    def load_children4(self,parents:list[PSObject],min_connectivity=0,depth=0,
                       max_childs=ALL_CHILDS,max_threads=50)->tuple[set[PSObject],list[PSObject]]:
        '''
            Input:
                if "parents" is empty will use all nodes from self.Graph
                depth - ontology depth to get children from. If depth=0, all children from all depths are returned
                max_childs - maximum number of children allowed in parent - cutoff for high level ontology concepts
                min_connectivity - avoid from using.  min_connectivity > 0 slows down children retrival significantly

            Output:
                {PSObject},[PSObject] - set of all found children, list of all parents with children\n
                if max_childs > 0, parents with number of children exeeding max_childs are excluded
            
            Updates:
                parents in self.Graph with CHILDS properties as [PSObject] children
        '''
        if parents:
            my_parents = list(set(parents))
            self.Graph.add_psobjs(set(parents))
        else:
            my_parents = self.Graph._get_nodes()
    
        return self._load_children4(my_parents,min_connectivity,depth,max_childs,max_threads)
    

    def child_graph(self, propValues:list, search_by_properties=[],include_parents=True):
        if not search_by_properties: search_by_properties = ['Name','Alias']
        oql_query = OQL.get_childs(propValues,search_by_properties,include_parents=include_parents)
        prop_val_str = ','.join(propValues)
        request_name = f'Find ontology children for {prop_val_str}'
        ontology_graph = self.process_oql(oql_query,request_name)
        if isinstance(ontology_graph, ResnetGraph):
            print('Found %d ontology children' % len(ontology_graph))
        return ontology_graph


    def ontology_graph(self,members:list,add2parent:PSObject):
        """
        Input
        -----
        members - [PSObjects]
        """
        ontology_graph = ResnetGraph()
        parent_id = add2parent.uid()
        ontology_graph.add_psobj(add2parent)
        ontology_graph.add_psobjs(set(members))
        for m in members:
            child_id = m.uid()
            rel = PSRelation({'ObjTypeName':['MemberOf'],'Relationship':['is-a'],'Ontology':['Pathway Studio Ontology']})
            rel.Nodes['Regulators'] = [(child_id,0,0)]
            rel.Nodes['Targets'] = [(parent_id,1,0)]
            rel[DBID].append(child_id)
            rel['URN'].append(child_id.urn()) #fake URN
            ontology_graph.add_edge(child_id,parent_id,relation=rel)

        self.add_rel_props(['Relationship','Ontology'])
        return ontology_graph


    def __ontology_graph(self,for_childs=[],depth=1):
        """
        Adds
        ----
        ontology parents to_graph\n
        annotates each parent with CHILDS property

        Return
        ------
        ResnetGraph with edges (child,parent,relation=[MemberOf, is_a, 'Pathway Studio Ontology']
        """
        my_childs = for_childs if for_childs else self.Graph._get_nodes()
        child_dbids = ResnetGraph.dbids(my_childs)

        my_session = self._clone_session(what2retrieve=NO_REL_PROPERTIES)
        my_session.max_threads4ontology = 50
        get_parent_query = 'SELECT Entity WHERE InOntology (SELECT Annotation WHERE Ontology=\'Pathway Studio Ontology\' AND Relationship=\'is-a\') inRange {steps} over (SELECT OntologicalNode WHERE id = ({ids}))'
        oql_query = get_parent_query.format(steps=str(depth),ids=','.join(list(map(str,child_dbids))))
        request_name = f'Find parents of {len(my_childs)} nodes with depth {depth}'
        parent_graph = my_session.process_oql(oql_query,request_name)
        if isinstance(parent_graph, ResnetGraph):
            my_session.load_children4(parent_graph._get_nodes(),depth=1)
            # loading only immediate children for all parents to not create shortcuts path in ontology graph
            my_session.close_connection()
        return my_session.Graph.ontology_graph()


    def load_ontology(self, parent_ontology_groups:list)->defaultdict[str,list[str]]:
        """
        Input
        -----
        [List of Ontology Parent Names]
        
        Return
        ------
        self.child2parent = {child_name:[ontology_parent_names]}
        """
        self.child2parent = defaultdict(list)
        self.add_ent_props(['Name'])
        for group_name in parent_ontology_groups:
            self.child_graph([group_name],['Name'],include_parents=False)
            for id, child in self.Graph.nodes(data=True):
                child_name = child['Name'][0]
                self.child2parent[child_name].append(group_name)

        return self.child2parent
 


    def ontopaths2(self,parent2childs:dict[str,list[PSObject]], ontology_depth:int)->dict[str,str]:
        '''
        Input
        -----
        name2childs = {name:[PSObject]}, where name = parent concept with ontology children listed in [PSObject]\n
        name2childs is generated by SemanticSearch:entities() 

        Return
        ------
        {name:paths}, where paths has format: path1<EOF>path2<EOF>...
        '''
        list_of_objlists = list(parent2childs.values())
        children = unpack(list_of_objlists)
        ontology_graph = self.__ontology_graph(children,ontology_depth)

        name2paths = dict()
        path_sep = '->'
        for concept_name, childs in parent2childs.items():
            for child in childs:
                all_parent_paths = ontology_graph.all_paths_from(child)
                for path in all_parent_paths:
                    parent_path_names = [x.name() for x in path]
                    try:
                        name2paths[concept_name] += '\n'+path_sep+path_sep.join(parent_path_names[1:])
                        # first elemnt in component path == concept_name
                    except KeyError:
                        name2paths[concept_name] = path_sep+path_sep.join(parent_path_names[1:])

        return name2paths


###################################  WRITE DUMP  CACHE, WRITE, DUMP  CACHE, WRITE DUMP  CACHE ##############################################
    def to_csv(self, file_out, in_graph=ResnetGraph(), access_mode='w'):
        debug = not self.no_mess
        my_graph = in_graph if in_graph else self.Graph
        my_graph.print_references(file_out, self.relProps, self.entProps, access_mode, printHeader=self.__IsOn1st_page,
                                col_sep=self.sep,debug=debug,single_rel_row=self.print_rel21row)


    def to_pandas (self, in_graph=None, RefNumPrintLimit=0)-> 'df':
        if not isinstance(in_graph, ResnetGraph): in_graph = self.Graph
        return df(in_graph.ref2pandas(self.relProps,self.entProps,RefNumPrintLimit,single_rel_row=self.print_rel21row))


    @staticmethod
    def pretty_xml(xml_string:str, remove_declaration = False):
        '''
        xml_string must have xml declration
        '''
        pretty_xml = str(minidom.parseString(xml_string).toprettyxml(indent='   '))
        if remove_declaration:
            pretty_xml = pretty_xml[pretty_xml.find('\n')+1:]
        return pretty_xml


    @staticmethod
    def filename4(name: str) -> str:
        """
        Normalizes a filename by replacing illegal characters.

        Args:
            name: The filename to normalize.

        Returns:
            The normalized filename with illegal characters replaced.
        """
        replacements = {'>': '-', '<': '-', '|': '-', '/': '-', ':': '_'}
        return "".join(replacements.get(char, char) for char in name)


    def __find_folder(self,folder_name:str,in_parent_folder='',root_dump_folder=''):
        '''
        Input:
            if self.data_dir has several dump directories with the same subdirectory names 
            use root_folder to direct search in specific subdirectory in self.data_dir
        Return:
            self.data_dir string if path does not exist
        '''
        folder_legal_name = self.filename4(folder_name)
        parent_folder_legal_name = self.filename4(in_parent_folder)
        root_folder_legal_name = self.filename4(root_dump_folder)

        search_dir = os.path.join(self.data_dir,root_folder_legal_name)    
        longest_path = search_dir
        for dirpath, dirnames,_ in os.walk(os.path.abspath(search_dir)):
            if folder_legal_name in dirnames and (not in_parent_folder or parent_folder_legal_name in dirpath):
                path = os.path.join(dirpath, folder_legal_name)
                if len(path) > len(longest_path):
                    longest_path = path
        return longest_path
   

    def dump_path(self,folder_name:str, parent_folder='',root_dump_folder='') -> str:
        """Creates a directory, optionally within a parent folder.

        Args:
            folder_name: The database name of the folder to create.
            parent_folder: (Optional) The name of the parent folder.

        Returns:
            The full path of the created directory.
        """
        if parent_folder:
            my_root = '' if parent_folder == root_dump_folder else root_dump_folder
            parent_path = self.__find_folder(parent_folder,'',my_root)
        else:
            root_dump_folder_legal_name = self.filename4(root_dump_folder)
            parent_path = os.path.join(self.data_dir,root_dump_folder_legal_name)

        if folder_name != root_dump_folder:
            folder_legal_name = self.filename4(folder_name)
            full_path = os.path.join(parent_path, folder_legal_name)
        else:
           full_path = parent_path

        try:
            os.makedirs(full_path, exist_ok=True)  # Create directory, including parents if needed
            # If the directory already exists, the function will do nothing and will not raise an error. 
            # It will only create the directory and any missing parent directories if they don't already exist.
        except OSError as e: # Handle potential errors (e.g., permissions issues)
            pass  # Or raise a more specific exception

        return full_path


    def __dump_base_name(self,folder_name:str):
        return 'content of '+ APISession.filename4(folder_name)+'_'


    def __dumpfiles(self, in_folder:str,in2parent_folder='', root_folder='', with_extension='rnef'):
        folder_path = self.dump_path(in_folder,in2parent_folder,root_folder)
        base_name = self.__dump_base_name(in_folder)
        listing = glob.glob(os.path.join(folder_path, base_name+'*.' + with_extension))
        return folder_path, base_name, len(listing)


    @staticmethod
    def __dumpfile_has_closing_batch_tag(dump_file:str):
        try:
            with open(dump_file, 'rb') as f:  # Open in binary mode to avoid potential line ending issues
                f.seek(-2, 2)  # Seek to the second-to-last byte (to handle newline variations)
                while f.read(1) != b'\n':  # Move backward until a newline is found
                    if not f.tell():
                        break
                    f.seek(-2, 1)
                last_line = f.readline().decode('utf-8').strip()  # Read and decode the last line
                return "</batch>" == last_line  # Check if "</batch>" is present
        except FileNotFoundError:
            return False # directory has no RNEF 


    def was_downloaded(self,folder_name:str,in_parent_folder='',root_folder=''):
        parent_folder = '' if folder_name == in_parent_folder else in_parent_folder
        last_dump_file = self.__make_dumpfile_path(folder_name,parent_folder,root_folder)
        return APISession.__dumpfile_has_closing_batch_tag(last_dump_file)


    def __make_dumpfile_path(self,of_folder:str,in2parent_folder='',root_folder='',new=False,start_index=0):
        '''
        Input
        -----
        new - forces to return name for new dumpfile
        filecount - if not zero forces to create dump file with index = filecount+1
        '''
        folder_path,of_folder_base_name,file_index = self.__dumpfiles(of_folder,in2parent_folder,root_folder,'rnef')
        if start_index:
            file_index = start_index
        if new:
            file_index += 1
        else:
            last_dumpfile = os.path.join(folder_path,of_folder_base_name+str(file_index)+'.rnef')
            if APISession.__dumpfile_has_closing_batch_tag(last_dumpfile):
                file_index += 1

        return os.path.join(folder_path,of_folder_base_name+str(file_index)+'.rnef')


    def close_rnef_dump(self,for_folder='',in_parent_folder='',root_folder='',check_last_tag=False):
        parent_folder = '' if for_folder == in_parent_folder else in_parent_folder
        last_dump_file = self.__make_dumpfile_path(for_folder,parent_folder,root_folder)
        if os.path.exists(last_dump_file):
            if check_last_tag:
                if APISession.__dumpfile_has_closing_batch_tag(last_dump_file):
                    return
            f = open(last_dump_file,'a',encoding='utf-8')
            f.write('</batch>')
            f.close()


    def _2rnefs(self,graph=ResnetGraph(),add_rel_props:dict={},add_pathway_props:dict={}):
        '''
        used for dumping pathway, group and results objects by FolderContent
        Returns
        -------
        graph RNEF XML with single <resnet> section and session properties for nodes and edges
        '''
        my_graph = graph if graph.number_of_nodes()>0 else self.Graph
        rel_props = [p for p in self.relProps if p not in NO_RNEF_REL_PROPS]
        return my_graph.to_rnefstr(self.entProps,rel_props,add_rel_props,add_pathway_props)
    

    def rnefs2dump(self,rnef_xml:str,to_folder='',in2parent_folder='',root_folder='',
        can_close=True,lock=None):
        '''
        # set can_close=False to continue dumping into last dump file
        Dumps
        -----
        "rnef_xml" string into 'to_folder' inside 'in2parent_folder' located in "self.data_dir"
        if size of dump file exceeds "max_rnef_size", "rnef_xml" is splitted into several RNEF files\n
        dump RNEF files are named as: 'content of to_folder#', where # - dump file number
        '''
        if lock is None: lock = threading.Lock()

        with lock:
            write2 = self.__make_dumpfile_path(to_folder,in2parent_folder,root_folder)
            if Path(write2).exists():
                file_size = os.path.getsize(write2)
                if file_size < self.max_rnef_size:
                    with open(write2,'a',encoding='utf-8') as f:
                        f.write(rnef_xml)
                        f.flush()
                        if os.path.getsize(write2) > self.max_rnef_size and can_close:
                            # need to create new dump file ASAP for next thread to write into it
                            new_write2 = self.__make_dumpfile_path(to_folder,in2parent_folder,root_folder,new=True)
                            with open(new_write2,'w',encoding='utf-8') as new_f:
                                new_f.write('<batch>\n')
                                new_f.flush()
                            f.write('</batch>')
                        f.flush()
                else:
                    # file_size >= self.max_rnef_size
                    if can_close:
                        # need to create new dump file ASAP for next thread to write into it
                        new_write2 = self.__make_dumpfile_path(to_folder,in2parent_folder,root_folder,new=True)
                        with open(new_write2,'w',encoding='utf-8') as new_f:
                            new_f.write('<batch>\n')
                            new_f.write(rnef_xml)
                            new_f.flush()

                        with open(write2,'a',encoding='utf-8') as f:
                            f.write('</batch>')
                            f.flush()
                    else:
                        with open(write2,'a',encoding='utf-8') as f:
                            f.write(rnef_xml)
                            f.flush()                        
            else:
                write2 = self.__make_dumpfile_path(to_folder,in2parent_folder,root_folder,new=True)
                with open(write2,'w',encoding='utf-8') as new_f:
                    new_f.write('<batch>\n')
                    new_f.write(rnef_xml)
                    new_f.flush() 


    def _dump2rnef(self,graph=ResnetGraph(),to_folder='',in_parent_folder='',root_folder='',can_close=True,lock=None):
        '''
        Input
        -----
        to_folder - name of the folder with downloaded RNEF data\n
        in_parent_folder - name of the folder in "self.data_dir" containinig "to_folder". Use "in_parent_folder" to create structured dump mirroring the folder structure in database\n
        By default "in_parent_folder" does not exist and "to_folder" is in "self.data_dir"\n
        can_close - <batch> element in dump RNEF file can be closed if file size exceeds self.max_rnef_size to complete current RNEF XML and start  new file if necessary\n
        set "can_close" to True to allow adding additional elements to RNEF XML, such as folder structure at the end of the dump

        Dumps
        -----
        large graph objects into several RNEF XML files
        
        Return
        -------
        execution time
        '''
        dump_start = time.time()
        my_graph = graph.copy() if graph else self.Graph.copy()
        # making input graph copy to release it for multithreading
        if my_graph:
            my_graph = my_graph.remove_undirected_duplicates()
            section_rels = set()
            for r, t, e in my_graph.edges(data='relation'):
                section_rels.add(e)
                if len(section_rels) == self.resnet_size:
                    resnet_section = my_graph.subgraph_by_rels(list(section_rels))
                    rnef_str = resnet_section.to_rnefstr(ent_props=self.entProps,rel_props=self.relProps)
                    rnef_str = self.pretty_xml(rnef_str,remove_declaration=True)
                    # dumps section
                    self.rnefs2dump(rnef_str,to_folder,in_parent_folder,root_folder,can_close,lock)
                    resnet_section.clear_resnetgraph()
                    section_rels.clear()
        
            # dumps leftover resnet_section with size < self.resnet_size
            resnet_section = my_graph.subgraph_by_rels(list(section_rels))
            rnef_str = resnet_section.to_rnefstr(ent_props=self.entProps,rel_props=self.relProps)
            rnef_str = self.pretty_xml(rnef_str,remove_declaration=True)
            self.rnefs2dump(rnef_str,to_folder,in_parent_folder,root_folder,can_close,lock)

            if my_graph.number_of_edges() == 0:
                rnef_str = my_graph.to_rnefstr(ent_props=self.entProps,rel_props=self.relProps)
                rnef_str = self.pretty_xml(rnef_str,remove_declaration=True)
                self.rnefs2dump(rnef_str,to_folder,in_parent_folder,root_folder,can_close,lock)

            if not self.no_mess:
                print('RNEF dump of "%s" graph into %s folder was done in %s' % 
                    (my_graph.name,to_folder,execution_time(dump_start)),flush=True)
                
        return time.time()-dump_start


    def download_oql(self,oql_query,request_name:str,resume_page=0,threads=25,close_batch=True,lock=None):
        # Use for oql_query producing large results
        # 50 threads crashed connection
        '''
        Input
        -----
        threads - number of threads to use in one download iteration\n
        Inreasing "threads" accelerates download but slows down writing cache file to disk\n
        if cache file writing time exceeds 5 min (Apache default timeout) download will crash
        
        Dumps
        -----
        results of oql_query to self.dump_folder in self.data_dir.\n
        Splits dump into small files smaller than self.max_rnef_size\n
        use 
        '''
        self.__replace_goql(oql_query)
        reference_counter = 0
        start_time = time.time()
        return_type = 'relations' if self.getLinks else 'entities'
        resume_pos = resume_page*self.PageSize
        iterations_graph = self.__init_session(resume_pos)
        self.clear() # to save RAM self.Graph is cleared
        if iterations_graph:
            number_of_iterations = int(self.ResultSize / self.PageSize)
            if number_of_iterations:
                result_size = self.ResultSize
                print(f'\n\nrequest "{request_name}" found {result_size} {return_type}. Begin download in {number_of_iterations} iterations')
                with ThreadPoolExecutor(1,f'{request_name} Dump{number_of_iterations}iters') as e:
                    # max_workers = 1 otherwise _dump2rnef gets locked
                    for i in range(0,number_of_iterations,threads):
                        iterations_graph = iterations_graph.compose(self.__thread__(threads))
                        e.submit(self._dump2rnef, iterations_graph.copy(), self.dump_folder,'',True,lock)
                        page_ref_count = iterations_graph.weight()
                        reference_counter += page_ref_count
                        remaining_iterations = number_of_iterations-i-threads
                        exec_time, remaining_time = execution_time2(start_time,remaining_iterations,number_of_iterations) 
                        print("With %d in %d iterations, %d %s in %d results with %d references retrieved in %s using %d threads" % 
                        (i+threads,number_of_iterations,self.ResultPos,return_type,result_size,reference_counter,exec_time,threads))
                        if i < number_of_iterations:
                            print(f'Estimated remaining retrieval time: {remaining_time}')
                        self.clear()
                        iterations_graph.clear()
                    e.shutdown()
            else:
                self._dump2rnef(iterations_graph, to_folder=self.dump_folder,lock=lock)
            
            if close_batch:
                self.close_rnef_dump(for_folder=self.dump_folder)
            self.clear()
            self.ResultRef = str()
            self.ResultPos = int()
            self.ResultSize = int()
            self.__IsOn1st_page = True
            return ResnetGraph() # fake return for consistency with process_oql
        else:
           return ResnetGraph() # fake return for consistency with process_oql 
    

    def common_neighbors(self,with_entity_types:list,of_dbids1:list,reltypes12:list,dir12:str,
                                and_dbids3:list,reltypes23:list,dir23:str):
        """
            Input
            -----
            dir12,dir23 must be '<','>' or empty string

            Returns
            -------
            ResnetGraph with 1-->2-->3 relations, where nodes in position 1 are from "of_dbids1" and nodes in position 3 are from "and_dbids3"\n
            if nodes in positon 1 and positon 3 are neighbors returns graph to complete 3-vertex cliques
        """
        sel_rels = 'SELECT Relation WHERE objectType = ({})'
        sel_rels12 = sel_rels.format(','.join(reltypes12))
        sel_rels23 = sel_rels.format(','.join(reltypes23))
        

        find_neighbors_oql = r'SELECT Entity objectType = ('+','.join(with_entity_types)+')'

        find_neighbors_oql = find_neighbors_oql+' AND Connected by ('+sel_rels12+') to (SELECT Entity WHERE id = ({ids1}))'
        find_neighbors_oql += ' AND Connected by ('+sel_rels23+') to (SELECT Entity WHERE id = ({ids2}))'

        r_n = f'Find entities common between {len(of_dbids1)} and {len(and_dbids3)} entities' 
        neighbors_entity_graph = self.iterate_oql2(f'{find_neighbors_oql}',set(of_dbids1),set(and_dbids3),request_name=r_n)
        neighbors_dbids = list(neighbors_entity_graph.dbids4nodes())
        print('Found %d common neighbors' % len(neighbors_dbids))

        if not neighbors_dbids: return ResnetGraph()
        if dir12 == '>':
            sel_12_graph_oql = sel_rels12 + r' AND downstream NeighborOf (SELECT Entity WHERE id = ({ids1})) AND upstream NeighborOf (SELECT Entity WHERE id = ({ids2}))'
        elif dir12 == '<':
            sel_12_graph_oql = sel_rels12 + r' AND upstream NeighborOf (SELECT Entity WHERE id = ({ids1})) AND downstream NeighborOf (SELECT Entity WHERE id = ({ids2}))'
        else:
            sel_12_graph_oql = sel_rels12 + r' AND NeighborOf (SELECT Entity WHERE id = ({ids1})) AND NeighborOf (SELECT Entity WHERE id = ({ids2}))'
        rn12 = f'Find 1-2 relations between {len(of_dbids1)} entities in position 1 and {len(neighbors_dbids)} common neighbors'

        if dir23 == '>':
            sel_23_graph_oql = sel_rels23 + r' AND upstream NeighborOf (SELECT Entity WHERE id = ({ids1})) AND downstream NeighborOf (SELECT Entity WHERE id = ({ids2}))'
        elif dir12 == '<':
            sel_23_graph_oql = sel_rels23 + r' AND downstream NeighborOf (SELECT Entity WHERE id = ({ids1})) AND upstream NeighborOf (SELECT Entity WHERE id = ({ids2}))'
        else:
            sel_23_graph_oql = sel_rels23 + r' AND NeighborOf (SELECT Entity WHERE id = ({ids1})) AND NeighborOf (SELECT Entity WHERE id = ({ids2}))'
        rn23 = f'Find 2-3 relations between {len(and_dbids3)} entities in position 3 and {len(neighbors_dbids)} common neighbors'
       # graph23 = self.iterate_oql2(f'{sel_23_graph_oql}',and_dbids3,neighbors_dbids,request_name=rn23)

        with ThreadPoolExecutor(max_workers=1, thread_name_prefix='FindCommonNeighbors') as executor:
            future12 = executor.submit(self.iterate_oql2,f'{sel_12_graph_oql}',set(of_dbids1),set(neighbors_dbids), True,rn12) 
            future23 = executor.submit(self.iterate_oql2,f'{sel_23_graph_oql}',set(and_dbids3),set(neighbors_dbids), True,rn23)
            
            graph12 = future12.result()
            graph23 = future23.result()
            executor.shutdown()

        graph12 = graph12.compose(graph23)
        return graph12


    def common_neighbors_with_effect(self,with_entity_types:list,of_dbids1:list,reltypes12:list,dir12:str,
                                and_dbids3:list,reltypes23:list,dir23:str,id_type='id')->"ResnetGraph":
        # written to circumvent bug in Patwhay Studio
        """
            Input
            -----
            dir12,dir23 must be '<','>' or empty string
            specify "id_type" == 'id' to process databse identifiers in "of_dbids1" and "and_dbids3", otherwise they will be treated as list of properties

            Returns
            -------
            ResnetGraph with 1-->2-->3 relations, where nodes in position 1 are from "of_dbids1" and nodes in position 3 are from "and_dbids3"\n
            if nodes in positon 1 and positon 3 are neighbors returns graph to complete 3-vertex cliques
        """
        hint = 'ids' if id_type == 'id' else 'props'
        sel_rels = 'SELECT Relation WHERE objectType = ({}) AND NOT Effect = unknown'
        sel_rels12 = sel_rels.format(','.join(reltypes12))
        sel_rels23 = sel_rels.format(','.join(reltypes23))

        find_neighbors_oql = 'SELECT Entity objectType = ('+','.join(with_entity_types)+')'
        find_neighbors_oql = find_neighbors_oql+' AND Connected by ('+sel_rels12+') to (SELECT Entity WHERE '+id_type+' = ({'+ hint +'1}))'
        find_neighbors_oql += ' AND Connected by ('+sel_rels23+') to (SELECT Entity WHERE '+id_type+' = ({'+hint+'2}))'

        r_n = f'Find entities common between {len(of_dbids1)} and {len(and_dbids3)} entities' 
        neighbors_entity_graph = self.iterate_oql2(f'{find_neighbors_oql}',set(of_dbids1),set(and_dbids3),request_name=r_n)
        
        neighbors = neighbors_entity_graph._get_nodes()
        if not neighbors: return ResnetGraph()
        if id_type == 'id':
            neighbors_ids = ResnetGraph.dbids(neighbors)
        else:
            neighbors_ids = [n[id_type][0] for n in neighbors]
        print('Found %d common neighbors' % len(neighbors_ids))

        if dir12 == '>':
            sel_12_graph_oql = sel_rels12 + ' AND downstream NeighborOf (SELECT Entity WHERE '+id_type+' = ({'+hint+'+1})) AND upstream NeighborOf (SELECT Entity WHERE '+id_type+' = ({'+hint+'2}))'
        elif dir12 == '<':
            sel_12_graph_oql = sel_rels12 + ' AND upstream NeighborOf (SELECT Entity WHERE '+id_type+r' = ({'+hint+'1})) AND downstream NeighborOf (SELECT Entity WHERE '+id_type+' = ({'+hint+'2}))'
        else:
            sel_12_graph_oql = sel_rels12 + ' AND NeighborOf (SELECT Entity WHERE '+id_type+' = ({'+hint+'1})) AND NeighborOf (SELECT Entity WHERE '+id_type+' = ({'+hint+'2}))'
        rn12 = f'Find 1-2 relations between {len(of_dbids1)} entities in position 1 and {len(neighbors_ids)} common neighbors'
       

        if dir23 == '>':
            sel_23_graph_oql = sel_rels23 + ' AND upstream NeighborOf (SELECT Entity WHERE '+id_type+' = ({'+hint+'1})) AND downstream NeighborOf (SELECT Entity WHERE '+id_type+' = ({'+hint+'2}))'
        elif dir12 == '<':
            sel_23_graph_oql = sel_rels23 + ' AND downstream NeighborOf (SELECT Entity WHERE '+id_type+' = ({'+hint+'1})) AND upstream NeighborOf (SELECT Entity WHERE '+id_type+' = ({'+hint+'}2))'
        else:
            sel_23_graph_oql = sel_rels23 + ' AND NeighborOf (SELECT Entity WHERE '+id_type+' = ({'+hint+'1})) AND NeighborOf (SELECT Entity WHERE '+id_type+' = ({'+hint+'2}))'
        rn23 = f'Find 2-3 relations between {len(and_dbids3)} entities in position 3 and {len(neighbors_ids)} common neighbors'

        graph12 = self.iterate_oql2(sel_12_graph_oql,set(of_dbids1),set(neighbors_ids),True,rn12)
        graph23 = self.iterate_oql2(sel_23_graph_oql,set(and_dbids3),set(neighbors_ids),True,rn23)
        graph123 = graph12.compose(graph23)
        # to remove relations with no Effect annotation:
        rels_with_effect = list(graph123.psrels_with(['positive','negative'],['Effect']))
        graph123_effect = graph123.subgraph_by_rels(rels_with_effect)
        if self.add2self:
            self.Graph = self.Graph.compose(graph123_effect)
        return graph123_effect


    def gv2gene(self,gv_dbids:list):
        """
        Input
        -----
        list of GV ids

        Returns
        -------
        {gv_id:[gene_names]}
        """
        prot2gvs_graph = ResnetGraph()
        print ('Finding genes for %d genetic variants' % len(gv_dbids))
        number_of_iterations = int(len(gv_dbids)/1000)+1
        for i in range(0, len(gv_dbids),1000):
            chunk = gv_dbids[i: i+1000]
            oql_query = OQL.expand_entity(PropertyValues=chunk, SearchByProperties=['id'], 
                                expand_by_rel_types=['GeneticChange'],expand2neighbors=['Protein'])
            request_name = f'{str(int(i/1000)+1)} iteration in {str(number_of_iterations)} to find genes linked to GVs'
            gvs = self.process_oql(oql_query,request_name)
            if isinstance(gvs,ResnetGraph):
                prot2gvs_graph.add_graph(gvs)

        # making gvid2genes for subsequent annotation
        gvid2genes = dict()
        for gv_id, protein_id, rel in prot2gvs_graph.edges.data('relation'):
            protein_node = prot2gvs_graph._psobj(protein_id)
            gene_name = protein_node['Name'][0]
            gvid2genes[gv_id].append(gene_name)

        return gvid2genes
    
 
    def _props2psobj(self, propValues: list, search_by_properties=[], get_childs=True,
                             only_obj_types=[], add2self=True):
        '''
        Annotates
        ---------
        nodes in self.Graph with CHILDS property if get_childs == True  

        Returns
        -------
        {PSObject} containing all parents and children combined
        '''
        prop_names = search_by_properties if search_by_properties else ['Name','Alias']
        prop_names_str = OQL.join_prop_names(prop_names)
        oql_query =  f'SELECT Entity WHERE ({prop_names_str})'+ ' = ({props})'
        if only_obj_types:
            obj_types_str = OQL.join_with_quotes(only_obj_types)
            oql_query = oql_query + ' AND objectType = (' + obj_types_str + ')'

        rn = f'Retrieving objects from database using {len(propValues)} values for {prop_names} properties'
        # string properties can form long oql queries exceeding 65500 chars limit
        max_id_len = max(len(s) for s in propValues)
        iteration_size = min(1000, int(MAX_OQLSTR_LEN/max_id_len))
        my_psobjs = self.iterate_oql(oql_query,set(propValues),request_name=rn,step=iteration_size)._get_nodes()

        if get_childs:
            children, parents_with_children = self.load_children4(my_psobjs)
            my_psobjs = self.Graph._get_nodes(ResnetGraph.uids(my_psobjs))+list(children)
        elif add2self:
            self.Graph.add_psobjs(set(my_psobjs))
        
        return my_psobjs


    def __dbid4prop(self,prop:str,psobjs:list[PSObject])->tuple[set[PSObject],set[PSObject]]:
        prop_values = unpack([list(o.get_props(prop)) for o in psobjs])
        prop2objs,_ = self.map_props2objs(prop_values,[prop])
        prop2dbid = dict()
        for propval,objs in prop2objs.items():
            for obj in objs:
                prop2dbid[propval] = obj.dbid()

        mapped_objs = set()
        notmapped_objs = set(psobjs)
        for psobj in psobjs:
            for propval in psobj.get_props(prop):
                try:
                    my_dbid = prop2dbid[propval]
                    psobj.update_with_value(DBID,my_dbid)
                    mapped_objs.add(psobj)
                    notmapped_objs.discard(psobj)
                except KeyError:
                    continue
        return mapped_objs, notmapped_objs


    def load_dbids4(self,psobjs:list[PSObject]):
        '''
        Return
        ------
        mapped_objs,no_dbid_objs - [PSObject],[PSObject]
        mapping is done first by Name then by URN
        '''
        
        print(f'Reterieving database identifiers for {len(psobjs)} entities using Name identifier')
        kwargs = {TO_RETRIEVE:NO_REL_PROPERTIES}
        my_session = self._clone_session(**kwargs)
        mapped_objs, notmapped_objs = my_session.__dbid4prop('Name',psobjs)
        
        if notmapped_objs:
            mo, notmapped_objs = my_session.__dbid4prop('URN',notmapped_objs)
            mapped_objs.update(mo)
        
        if notmapped_objs:
            [x.set_property('Alias',x.name()) for x in notmapped_objs]
            mo,_ = my_session.__dbid4prop('Alias',notmapped_objs)
            mapped_objs.update(mo)

        urn2dbids = {n.urn():n[DBID] for n in mapped_objs}
        self.Graph.set_node_annotation(urn2dbids,DBID)
        print(f'Loaded {len(mapped_objs)} database identitiers for {len(psobjs)} entities')
        my_session.close_connection()
        return mapped_objs,notmapped_objs
    

    def update_graph(self,g:ResnetGraph,with_node_props:list[str],with_rel_props:list[str],inplace=True):
        my_graph = g if inplace else g.copy()
        new_session = self._clone_session()
        new_session.add2self = False
        if with_node_props:
            new_session.add_ent_props(with_node_props)
            oql = 'SELECT Entity WHERE URN = ({props})'
            objs = g._get_nodes()
            urns = {r.urn() for r in objs}
            start = time.time()
            chunk_len = 10000
            number_of_iterations = math.ceil(len(urns)/chunk_len)
            for i in range(0,len(urns),chunk_len):
                req_name = f'Updating {chunk_len} entities with {len(with_node_props)} properties'
                new_g = new_session.iterate_oql(oql,urns,use_cache=False,request_name=req_name)
                new_g.load_references()
                my_graph.add_graph(new_g)
                time_passed, remaining_time = execution_time2(start,number_of_iterations-i,number_of_iterations)
                print(f'Retrieved {i+chunk_len} entities out of {len(objs)} in {time_passed}')
                print(f'Estimated remaining update time: {remaining_time}')

        if with_rel_props:
            new_session.add_rel_props(with_rel_props)
            oql = 'SELECT Relation WHERE URN = ({props})'
            objs = g.psrels_with()
            urns = list({r.urn() for r in objs})
            start = time.time()
            chunk_len = 10000
            number_of_iterations = math.ceil(len(urns)/chunk_len)
            iteration_counter = 0
            for i in range(0,len(urns),chunk_len):
                update_len = min(len(urns)-i,chunk_len)
                req_name = f'Updating {update_len} in {len(objs)} Relations with {len(with_rel_props)} properties'
                new_g = new_session.iterate_oql(oql,set(urns[i:i+chunk_len]),use_cache=False,request_name=req_name)
                new_g.load_references()
                #new_g_rels = new_g.psrels_with()
                my_graph.add_graph(new_g)
                iteration_counter +=1
                remaining_iteration = number_of_iterations-iteration_counter
                time_passed, remaining_time = execution_time2(start,remaining_iteration,number_of_iterations)
                retreived_len = min(len(urns),i+chunk_len)
                print(f'Retrieved {retreived_len} relations out of {len(objs)} in {time_passed}')
                print(f'Estimated remaining update time: {remaining_time}')

        return None if inplace else my_graph
    

    def child_update(self,g:ResnetGraph,make_new_rels=False):
        '''
            make_new_rels - if True will create for all child entities the same relations as their parent have
        '''
        myG = g.copy() if g else self.Graph.copy()
        graph_entities = myG._get_nodes()
        children, parents_with_children = self.load_children4(graph_entities)
        myG.add_psobjs(parents_with_children+list(children))

        visited_parents = set()
        add2g = ResnetGraph()
        for r,t,rel in myG.iterate():
            for i,parent in enumerate([r,t]):
                if parent not in visited_parents and parent in parents_with_children:
                    visited_parents.add(parent)
                    n_neighborhood = myG.neighborhood([parent])
                    n_neighbors = n_neighborhood._get_nodes()
                    n_neighbors.remove(parent)
                    parent_childs = parent.childs()
                    connections = self.connect_entities(n_neighbors,parent_childs)
                    add2g = add2g.compose(connections)
                    if make_new_rels:
                        for n1,n2,parent_rel in n_neighborhood.iterate():
                            # must use dict(parent_rel) to make new rel
                            if n1 == parent:
                                for c in parent_childs:
                                    add2g.add_triple(c,n2,dict(parent_rel),parent_rel.refs(),parent_rel.is_directional())
                            else:
                                for c in parent_childs:
                                    add2g.add_triple(n1,c,dict(parent_rel),parent_rel.refs(),parent_rel.is_directional())
                
        return myG.compose(add2g)


##################### EXPERIMENT EXPERIMENT EXPERIMENT EXPERIMENT ###########################
    def map_experiment(self, exp:Experiment):
        identifier_name = exp.identifier_name()
        identifier_names = identifier_name.split(',')
        identifiers = exp.list_identifiers()
        map2objtypes = exp['ObjTypeName']
        identifier2objs, objid2prop = self.map_props2objs(identifiers,identifier_names,map2objtypes)
        return exp.map(identifier2objs)
    

    def scopusCI(self):
        references = self.Graph.load_references()
        refs_with_ci, no_ci_refs = loadCI(self.APIconfig,references)

        graph_with_ci = ResnetGraph()
        graph_with_ci.add_psobjs(set(self.Graph._get_nodes()))
        for n1,n2,rel in self.Graph.edges(data='relation'):
            e_refs = set(rel.refs())
            e_refs_with_ci = {r for r in refs_with_ci if r in e_refs}
            e_refs.update(e_refs_with_ci)
            graph_with_ci.add_rel(rel)

        return graph_with_ci
    

    def map_graph(self,graph:ResnetGraph, map_by=['Name'])->tuple[ResnetGraph,set[PSObject]]:
        '''
        Input
        -----
        graph - ResnetGraph where nodes may have arbitrary URNs
        map_by - list of database properties for mapping. Nodes in "graph" must have same properties  

        Return
        ------
        ResnetGraph with nodes with database URNs
        '''
        print(f'Mapping input graph with {len(graph)} entities using {map_by} identifiers')
        kwargs = {TO_RETRIEVE:NO_REL_PROPERTIES}
        my_session = self._clone_session(**kwargs)
        my_session.add_ent_props(map_by)

        graph_psobjs = graph.psobjs_with(map_by)
        obj_props = graph.node_props(map_by,graph_psobjs)
        my_session._props2psobj(list(obj_props),map_by,get_childs=False)
        props2objs,uid2propval  = my_session.Graph.props2obj_dict([],map_by,case_insensitive=True)

        return graph.remap_graph(props2objs,map_by)
    

    def get_ppi(self, interactors:set[PSObject], fromRNEFs:list=[]):
        '''
        Return
        ------
        ResnetGraph with relation types: Binding, DirectRegulation, ProtModification 
        '''
        ppi_keeper = ResnetGraph()
        if fromRNEFs:
            for fname in fromRNEFs:
                ppi_keeper = ppi_keeper.compose(ResnetGraph.fromRNEF(fname,only4objs=interactors))
        else:
            interactors_dbids = ResnetGraph.dbids(list(interactors))
            splitter = list() #holds lists of ids splits 
            splitter.append(list(interactors_dbids))
            number_of_splits = int(math.log2(len(interactors_dbids)))
            
            for s in range(1, number_of_splits):
                new_splitter = list()
                half = int(len(splitter[0]) / 2)
                with ThreadPoolExecutor(max_workers=number_of_splits,thread_name_prefix='Retrieve PPI') as e:
                    futures = list()
                    job_name = f'split {s} in {number_of_splits}'
                    for split in splitter:
                        uq_list1 = split[0:half]
                        uq_list2 = split[half:]
                        oql_query = "SELECT Relation WHERE objectType = (Binding, DirectRegulation, ProtModification) "
                        oql_query += "AND NeighborOf (SELECT Entity WHERE id = ({ids1})) "
                        oql_query += "AND NeighborOf (SELECT Entity WHERE id = ({ids2}))"
                        new_session = self._clone_session() # need to clone since self.dbid2relation mutates
                        new_session.add2self = False
                        futures.append(e.submit(new_session.iterate_oql2,oql_query,set(uq_list1),set(uq_list2),request_name=job_name))
                        new_splitter.append(uq_list1)
                        new_splitter.append(uq_list2)
                    splitter = new_splitter
                    s += 1

                    for f in as_completed(futures):
                        ppi_keeper = ppi_keeper.compose(f.result())
                    e.shutdown()
        print(f'Retrieved protein-protein interaction network with {ppi_keeper.number_of_nodes()} nodes and {ppi_keeper.number_of_edges()} edges')
        return ppi_keeper


def open_api_session(api_config_file='',what2retrieve=1) -> APISession:
  APIconfig = load_api_config(api_config_file)
  return APISession(APIconfig,what2retrieve=what2retrieve)
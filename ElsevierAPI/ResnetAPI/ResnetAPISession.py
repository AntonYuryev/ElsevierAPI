import time, math, os, glob, json, threading
import networkx as nx
from zeep import exceptions
from xml.dom import minidom
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import timedelta

from .ZeepToNetworkx import PSNetworx, len
from .ResnetGraph import ResnetGraph,df,REFCOUNT,CHILDS
from .NetworkxObjects import DBID,PSObject,PSRelation
from .PathwayStudioGOQL import OQL
from .Zeep2Experiment import Experiment
from ..ETM_API.references import SENTENCE_PROPS,PS_BIBLIO_PROPS,PS_SENTENCE_PROPS,PS_REFIID_TYPES,RELATION_PROPS,ALL_PS_PROPS
from ..ScopusAPI.scopus import loadCI, SCOPUS_CI

TO_RETRIEVE = 'to_retrieve'
BELONGS2GROUPS = 'belongs2groups'
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
MAX_PAGE_THREADS = 100


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
    max_threads4ontology = 4
    # 4 threads perform a bit faster than 8 threads 
    # 288 disease parents out of 288 were processed in 0:19:51.732684 by 8 threads
    # 288 disease parents out of 288 were processed in 0:18:52.715937 by 4 threads
    

######################################  CONFIGURATION  ######################################
    def __init__(self,*args,**kwargs):
        '''
        Input
        -----
        APIconfig = args[0]\nwhat2retrieve - defaults NO_REL_PROPERTIES, other options -\n
        [DATABASE_REFCOUNT_ONLY,REFERENCE_IDENTIFIERS,BIBLIO_PROPERTIES,SNIPPET_PROPERTIES,ONLY_REL_PROPERTIES,ALL_PROPERTIES]
        no_mess - default True, if False your script becomes more verbose
        connect2server - default True, set to False to run script using data in __pscache__ files instead of database
        '''
        self.GOQLquery = str()
        self.DumpFiles = []

        supported_kwargs = {'what2retrieve':NO_REL_PROPERTIES,
                            'ent_props' : ['Name'],
                            'rel_props' : ['URN'],
                            'data_dir' : '',
                            'use_cache' : False,
                            'oql_queries' : []
                            }

        ent_props = kwargs.pop('ent_props',[])
        rel_props = kwargs.pop('rel_props',[])
        supported_kwargs.update(kwargs)
        super().__init__(*args, **supported_kwargs)
        
        self.set_dir(supported_kwargs.get('data_dir',''))
        self.use_cache = supported_kwargs.get('use_cache',False)# if True signals to use graph data from cache files instead of retrieving data from database using GOQL queries 

        self.__retrieve(supported_kwargs['what2retrieve']) #__retrieve overides self.relProps, self.entProps
        # properties have to updated after self.__retrieve
        self.add_ent_props(ent_props)
        self.add_rel_props(rel_props)

        self.add2self = True # if False will not add new graph to self.Graph
        self.print_rel21row = False
        self.getLinks = True
        self.__IsOn1st_page = True # to print header in dump file
        self.dump_folder = str()
        self.my_oql_queries = list(supported_kwargs.get('oql_queries',[])) # [(oql_query,request_name),...] list of GOQL queries tuples for fetching data from database
        self.ResultPos = int()
        self.ResultSize = int()

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
            return ['URN']+ALL_PS_PROPS
        return []


    def __retrieve(self,what2retrieve=CURRENT_SPECS):
        retrieve_props = self._what2retrieve(what2retrieve)
        if retrieve_props:
            self.relProps = retrieve_props
        else:
            return # do nothing, keep CURRENT_SPECS
             

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
        used by self.process_oql to clone session
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


    def set_dir(self,dir:str):
        self.data_dir = dir
        if dir and dir[-1] != '/':
            self.data_dir += '/'


    def start_download_from(self,result_position:int):
        self.ResultPos = result_position


    def add_rel_props(self, add_props:list):
        self.relProps = self.relProps+[i for i in add_props if i not in self.relProps]
        if set(SENTENCE_PROPS).intersection(add_props) and 'TextRef' not in self.relProps:
            self.relProps.append('TextRef')


    def add_ent_props(self, add_props: list):
        self.entProps = self.entProps+[i for i in add_props if i not in self.entProps]

    
    def add_dump_file(self, dump_fname, replace_main_dump=False):
        if replace_main_dump:
            self.DumpFiles = []
        self.DumpFiles.append(dump_fname)

#########################  MULTITHREADED PROCESSING, RETRIEVAL   ##################################
    def __init_session(self, first_iteration=0,max_result=0):
        '''
        Return
        ------
        empty ResnetGraph if number of results > max_result
        '''
        self.__set_get_links()
        obj_props = self.relProps if self.getLinks else self.entProps
        zeep_data, (self.ResultRef, self.ResultSize, self.ResultPos) = self.init_session(self.GOQLquery,
                                                                                        self.PageSize,
                                                                                        obj_props,
                                                                                        getLinks = self.getLinks)
        if first_iteration > 0:
            self.ResultPos = first_iteration
            print ('Resuming retrieval from %d position' % first_iteration)
            return self.__nextpage__(self.ResultPos)

        if max_result and self.ResultSize > max_result:
            return ResnetGraph()

        if type(zeep_data) != type(None):
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
        # current_pos does not change during multithreading initiation!!!!  
        if int(self.ResultPos) < int(self.ResultSize):
            obj_props = self.relProps if self.getLinks else self.entProps
            zeep_data, self.ResultSize, current_pos = self.get_session_page(self.ResultRef, result_pos, self.PageSize,
                                                                       self.ResultSize,obj_props,getLinks=self.getLinks)                                                          
            if type(zeep_data) != type(None):
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


    def process_oql(self,oql_query,request_name='',max_result=0) -> ResnetGraph|int:
        '''
        Return
        ------
        if max_result is not 0 and number of results exceeds max_result returns int = self.ResultSize\n
        otherwise returns ResnetGraph with query results
        '''
        self.__replace_goql(oql_query)
        start_time = time.time()
        return_type = 'relations' if self.getLinks else 'entities'
        entire_graph = self.__init_session(max_result=max_result)

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
                            self.execution_time(start_time)[0], pages+1),flush=True)

        self.ResultRef = ''
        self.ResultPos = 0
        self.ResultSize = 0
        self.__IsOn1st_page = True
        self.clone = False
        return entire_graph


    def __annotate_dbid_graph__(self,id_only_graph:ResnetGraph,use_cache=True,request_name=''):
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
            rels_add2return = self.__iterate__(rel_query,relations_dbids2retreive,req_name,step=1000)
            add2return = rels_add2return
            nodes_dbids2retreive = nodes_dbids2retreive.difference(add2return.dbids4nodes())
            
        if nodes_dbids2retreive:
            entity_query = 'SELECT Entity WHERE id = ({ids})'
            nodes_add2return = self.__iterate__(entity_query,nodes_dbids2retreive,req_name,step=1000)
            add2return = add2return.compose(nodes_add2return)
        
        return_subgraph = return_subgraph.compose(add2return)
        assert(return_subgraph.number_of_nodes() == len(need_nodedbids))
       # assert(return_subgraph.number_of_edges() >= len(relation_dbids2return))
       # privately owned relations are not returned by API. 
       # therefore, return_subgraph.number_of_edges() may have less relations than "relation_dbids2return"
       # 
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
            def join_list (l:list):
                return ','.join(list(map(str,l)))
        else:
            my_oql = oql_query.replace('{props','{ids')
            def join_list (l:list):
                return OQL.join_with_quotes(l)
            

        number_of_iterations = math.ceil(len(dbids)/step)  
        # new sessions must be created to avoid ResultRef pointer overwriting
        entire_graph = ResnetGraph()
        thread_name = f'Iterate_{len(dbids)}_ids_by_{step}'
        print(f'Retrieval of {len(dbids)} objects will be done in {self.max_sessions} parallel sessions' )
        dbids_list = list(dbids)
        lock = threading.Lock() if download else None
        with ThreadPoolExecutor(max_workers=self.max_sessions, thread_name_prefix=thread_name) as e:
            future_sessions = list()
            futures = list()
            for i in range(0,len(dbids_list), step):
                iter_dbids = dbids_list[i:i+step]
                oql_query_with_dbids = my_oql.format(ids=join_list(iter_dbids))
                iter_name = f'Iteration #{int(i/step)+1} out of {number_of_iterations} for "{request_name}" retrieving {len(iter_dbids)} ids'
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
        req_name = 'Retrieving database identifiers for '+request_name
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
            iteration1step = min(step,950) # string properties can form long oql queries exceeding zeep limits
            iteration2step = min(2*iteration1step, 950)
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
                print(f'Processed {min([i1+iteration1step,len(ids1)])} out of {len(ids1)} identifiers against {len(ids2)} identifiers')
        else:
            for i2 in range(0,len(ids2), iteration1step):
                iter_ids2 = ids2[i2:i2+iteration1step]
                iter_oql_query = id_oql.format(ids1=hint,ids2=join_list(iter_ids2))
                iter_graph = self.__iterate__(iter_oql_query,set(ids1),'iteration2',iteration2step,download)
                entire_graph = entire_graph.compose(iter_graph)
                print(f'Processed {min([i2+iteration1step,len(ids2)])} out of {len(ids2)} identifiers against {len(ids1)} identifiers')
            
        return entire_graph
    

    def __iterate_oql2__(self, oql_query:str, ids_or_props1:set, ids_or_props2:set, request_name='', step=500):
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

 
    def iterate_oql(self,oql_query:str,dbid_set:set,use_cache=True,request_name='',step=1000):
        """
        # oql_query MUST contain string placeholder called {ids} if iterable id_set contains dbids\n
        # oql_query MUST contain string placeholder called {props} if iterable id_set contain property values other than database id
        """
        print('Processing "%s" request\n' % request_name)
        dbid_only_graph = self.__iterate_oql__(oql_query,dbid_set,request_name,step)
        if dbid_only_graph:
            annoated_graph = self.__annotate_dbid_graph__(dbid_only_graph,use_cache,request_name)
            return annoated_graph
        else:
            return ResnetGraph()
   

    def iterate_download(self,oql_query:str,dbid_set:set,request_name='',step=1000):
        """
        # oql_query MUST contain string placeholder called {ids} if iterable id_set contains dbids\n
        # oql_query MUST contain string placeholder called {props} if iterable id_set contain property values other than database id
        """
        return self.__iterate_oql__(oql_query,dbid_set,request_name,step,download=True)


    def iterate_oql2(self,oql_query:str,dbid_set1:set,dbid_set2:set,use_cache=True,request_name='',step=500):
        """
        # oql_query MUST contain 2 string placeholders called {ids1} and {ids2}
        # oql_query MUST contain 2 string placeholder called {props1} and {props2}\n
        if iterable id_sets contain property values other than database id
        """
        print(f'Processing "{request_name}" request')
        dbid_only_graph = self.__iterate_oql2__(oql_query,dbid_set1,dbid_set2,request_name,step)
        
        if dbid_only_graph:
            print('Will retrieve annotated graph with %d entities and %d relations' 
                        %(dbid_only_graph.number_of_nodes(),dbid_only_graph.number_of_edges()))
            annotated_graph = self.__annotate_dbid_graph__(dbid_only_graph,use_cache,request_name)
            return annotated_graph
        else:
            return ResnetGraph()
     

    def clear(self):
        self.Graph.clear_resnetgraph()
        self.dbid2relation.clear()
        #self.parent_dbid2children.clear()


    def flush_dump_files(self):
        for f in self.DumpFiles:
            open(f, 'w').close()
            if not self.no_mess:
                print('File "%s" was cleared before processing' % f)


    @staticmethod
    def reopen(fname):
        open(fname, "w", encoding='utf-8').close()
        return open(fname, "a", encoding='utf-8')

    def get_result_size(self,oql_query):
        zeep_relations, (self.ResultRef, self.ResultSize, self.ResultPos) = self.init_session(oql_query, PageSize=1,
                                                                                              property_names=[])
        return self.ResultSize



###########################GET GET RETRIEVE GET RETRIEVE GET GET ################################################
    def __get_saved_results(self,psobjs:list):
        result_names = [o.name() for o in psobjs]
        oql_query = f'SELECT Result WHERE Name = ({OQL.join_with_quotes(result_names)})'
        results = self.load_graph_from_oql(oql_query,self.relProps,
                                                self.entProps,get_links=False)._get_nodes()
        
        def __has_relations(result):
            # undestands if search results contain entities or relations 
            # load  self.ResultRef
            self.ResultRef  = 'ID='+str(result.dbid())
            zeepdata,self.ResultSize, pos = self.get_session_page(self.ResultRef,0, 1,1,[],getLinks=True)
            first_obj = zeepdata.Objects.ObjectRef[0]
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


    def get_network(self,_4dbids:set,connect_by_reltypes:list=[],download=False,threads=25):
        if download:
            oql_query = 'SELECT Relation WHERE NeighborOf (SELECT Entity WHERE id = ({ids1})) AND NeighborOf (SELECT Entity WHERE id = ({ids2}))'
            if connect_by_reltypes:
                oql_query += f' AND objectType = ({str(",".join(connect_by_reltypes))})'

            splitter = list() #holds lists of ids splits 
            splitter.append(list(_4dbids))
            number_of_splits = int(math.log2(len(_4dbids)))
            print('Will load network of %d nodes in %d iterations' % (len(_4dbids),number_of_splits))
            for s in range(1, number_of_splits):
                new_splitter = list()
                half = int(len(splitter[0]) / 2)
                iter_start = time.time()
                ids1 = set()
                ids2 = set()
                for split in splitter:
                    uq_list1 = split[:half]
                    uq_list2 = split[half:]
                    ids1.update(uq_list1)
                    ids2.update(uq_list2)
                    new_splitter.append(uq_list1)
                    new_splitter.append(uq_list2)

                self.__iterate2__(oql_query,list(ids1),list(ids2),f'Download network: iteration #{s}',download=True)
                splitter = new_splitter
                executiontime = self.execution_time(iter_start)[0]
                print(f'Iteration {s} out of {number_of_splits} was completed in {executiontime}')
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

        req_name = f'Connecting {len(node_dbids1)} and {len(node_dbids2)} entities by {reltype_str} relations and {effect_str} effects'
         # step is set 250 to mitigate possible timeout during multithreading
        if download:
           self.__iterate2__(f'{oql_query}',list(node_dbids1),list(node_dbids2),req_name,step,download=True)
           return ResnetGraph()
        else:
            return self.iterate_oql2(f'{oql_query}',node_dbids1,node_dbids2,use_relation_cache,request_name=req_name,step=step)
   

    def get_group_members(self, group_names:list):
        groups = OQL.join_with_quotes(group_names)
        oql_query = f'SELECT Entity WHERE MemberOf (SELECT Group WHERE Name = ({groups}))'
        req_name = 'Find members of groups: ' + ','.join(group_names)
        graph2return = self.process_oql(oql_query, request_name=req_name)
        if isinstance(graph2return,ResnetGraph):
            if graph2return.number_of_nodes() > 0:
                print('%s groups are empty or do not exist in databse' % str(group_names))
            else:
                print('loaded %d members from %s' % (graph2return.number_of_nodes(),str(group_names)))
        return graph2return
    

    def add_group_annotation(self,group_names:list,_2graph=ResnetGraph()):
        urns2value = PSObject()
        for group_name in group_names:
            group_graph = self.get_group_members([group_name])
            group_members = group_graph._get_nodes()
            [urns2value.append_property(o.urn(),group_name) for o in group_members]

        if _2graph:
            _2graph.set_node_annotation(urns2value,BELONGS2GROUPS)
        else:
            self.Graph.set_node_annotation(urns2value,BELONGS2GROUPS)
        # my_graph = _2graph if _2graph else self.Graph does not work here
        return


    def map_props2objs(self,using_values:list,in_properties:list,map2types=[],case_insensitive=False):
        """
        Returns
        -------
        propval2objs = {prop_value:[PSObject]}\n
        objid2propval = {node_id:[prop_values]},\n  where prop_value is from 'using_values'
        """
        propval2objs,objid2propval = self.Graph.props2obj_dict(using_values, in_properties,case_insensitive)
        need_db_mapping = set(using_values).difference(propval2objs.keys())

        current_ent_props = self.entProps
        self.entProps = in_properties

        step = 1000
        iteration_counter = math.ceil(len(need_db_mapping) / step)
        print('Will use %d %s identifiers to find entities in %d iterations' % 
             (len(need_db_mapping), ','.join(in_properties), iteration_counter))

        for i in range(0, len(need_db_mapping), step):
            last = i + step
            propval_chunk = using_values[i:last]
            query_node = OQL.get_entities_by_props(propval_chunk,in_properties,map2types)
            request_name = 'Mapping {start}-{end} values out of {total}'
            end = min(last,len(using_values))
            request_name = request_name.format(start=str(i), end=str(end), total = len(using_values))
            self.process_oql(query_node,request_name)

        p2obj,id2p = self.Graph.props2obj_dict(list(need_db_mapping),in_properties,case_insensitive)
        propval2objs.update(p2obj)
        for id,prop_vals in id2p.items():
            try:
                mapped_values = set(objid2propval[id])
                mapped_values.update([str(x).lower() for x in prop_vals])
                objid2propval[id] = list(mapped_values)
            except KeyError:
                objid2propval[id] = prop_vals

        print("%d entities were mapped to %d attributes in %s properties" % 
                (len(objid2propval), len(using_values),','.join(in_properties)))

        self.entProps = current_ent_props
        return propval2objs, objid2propval

################################# ONTOLOGY  ONTOLOGY ############################################
    def __load_children(self,parent:PSObject,min_connectivity=0,depth=0,max_childs=0):
            '''
            Input
            -----
            if max_childs=0 finds all children for parent in database
            '''
            parent_dbid = parent.dbid()
            if parent_dbid:
                query_ontology = OQL.get_childs([parent_dbid],['id'],depth=depth)
                if min_connectivity:
                    query_ontology += f' AND Connectivity >= {min_connectivity}'

                my_session = self._clone_session()
                request_name = f'Retrieve ontology children for {parent.name()}'
                children_graph = my_session.process_oql(query_ontology,request_name,max_result=max_childs)
          #      if parent.name() == 'physical illness':
          #          print('')
                my_session.close_connection()
                if isinstance(children_graph,int):
                    # fake list of empty children to enable downstream remove_high_level_entities()
                    return parent, [PSObject()]*children_graph 
                else:
                    return parent, children_graph._get_nodes()
            else:
                print(f'Parent {parent.name()} has no database identifier. Its children cannot be loaded')
                return parent,[]
        

    def load_children4(self,parents:list,add2self=True,min_connectivity=0,depth=0,max_childs=0,max_threads=10):
        '''
        Input
        -----
        parents - [PSObject]\n
        depth - ontology depth to get children from. If depth=0, all children from all depths are returned\n
        max_childs - mximum number of children allowed in each parent.  Use it as cutoff for high level ontology concepts

        Return
        ------
        {PSObject},[PSObject] - set of all found children, list of all parents with children\n
        if max_childs > 0, parents with number of children exeeding max_childs are excluded
    
        Updates
        -------
        parents in self.Graph with CHILDS properties as [PSObject] children
        '''
        process_start = time.time()
        if parents:
            annotate_parents = list(set(parents))
            if add2self:
                self.Graph.add_psobjs(set(parents)) 
        else:
            annotate_parents = self.Graph._get_nodes()

        print(f'Loading ontology children for {len(annotate_parents)} entities')
        #max_threads = self.max_threads4ontology
        # by default sessions_max=200 in Oracle 
        # for Disease optimal max_threads = 5
        # for Protein+FunctionalCLass+Complex optimal max_threads = 50 
        thread_name_prefix = f'Childs 4 {len(annotate_parents)} prnts in {max_threads} thrds-'
        # 4 threads perform a bit faster than 8 threads 
        # 288 disease parents out of 288 were processed in 0:19:51.732684 by 8 threads
        # 288 disease parents out of 288 were processed in 0:18:52.715937 by 4 threads
        max_threaded_time = 0
        need_children = list()
        child_counter = set()
        iteration_counter = 0
        for p, parent in enumerate(annotate_parents):
            try:
                children = self.Graph.nodes[parent.uid()][CHILDS]
                child_counter.update(children)
            except KeyError:
                need_children.append(parent)
                if len(need_children) >= max_threads or p == len(annotate_parents)-1:       
                    thread_start = time.time()
                    iteration_counter += 1
                    thread_name = thread_name_prefix+str(iteration_counter) + 'iter'
                    with ThreadPoolExecutor(max_workers=len(need_children), thread_name_prefix=thread_name) as e:
                        futures = list()
                        [futures.append(e.submit(self.__load_children,nc,min_connectivity,depth,max_childs)) for nc in need_children]

                        for future in as_completed(futures):
                            f_parent, children = future.result()
                            if max_childs:
                                if len(children) > max_childs: # children has empty PSObjects
                                    nx.function.set_node_attributes(self.Graph, {f_parent.uid():{CHILDS:children}})
                                    continue
                            self.Graph.add_psobjs(children)
                            nx.function.set_node_attributes(self.Graph, {f_parent.uid():{CHILDS:children}})
                            child_counter.update(children)
                                # set_node_attributes() does not update nodes that do not exist in Graph
                        e.shutdown()
       
                    threaded_time = time.time() - thread_start
                    if threaded_time > max_threaded_time:
                        max_threaded_time = threaded_time
                        if max_threaded_time > 240:
                            print(f'Longest {max_threads}-threaded time {"{}".format(str(timedelta(seconds=max_threaded_time)))} \
                                  is close to Apache default 5 minutes transaction timeout !!!')
                    if not self.no_mess:
                        print (f'{p+1} parents out of {len(annotate_parents)} were processed in {self.execution_time(process_start)[0]}')
                    need_children.clear()
        
        print(f'{len(child_counter)} children for {len(annotate_parents)} parent entities were found in database \
              in {self.execution_time(process_start)[0]}')
        print(f'Longest {max_threads}-threaded time: {"{}".format(str(timedelta(seconds=max_threaded_time)))}')
        return child_counter, self.Graph.psobjs_with([CHILDS])


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
            rel.append_property(DBID,child_id)
            rel.append_property('URN', child_id.urn()) #fake URN
            ontology_graph.add_edge(child_id,parent_id,relation=rel)

        self.add_rel_props(['Relationship','Ontology'])
        return ontology_graph


    def __ontology_graph(self,for_childs=[],depth=1):
        """
        Adds
        ----
        ontology parents to_graph\n
        annotates each parent with CHILDS property

        Returns
        -------
        [PSObject] - list of parents with CHILDS property
        """
        my_childs = for_childs if for_childs else self.Graph._get_nodes()
        child_dbids = ResnetGraph.dbids(my_childs)

        my_session = self._clone_session(what2retrieve=NO_REL_PROPERTIES)
        my_session.max_threads4ontology = 25
        get_parent_query = 'SELECT Entity WHERE InOntology (SELECT Annotation WHERE Ontology=\'Pathway Studio Ontology\' AND Relationship=\'is-a\') inRange {steps} over (SELECT OntologicalNode WHERE id = ({ids}))'
        oql_query = get_parent_query.format(steps=str(depth),ids=','.join(list(map(str,child_dbids))))
        request_name = f'Find parents of {len(my_childs)} nodes with depth {depth}'
        parent_graph = my_session.process_oql(oql_query,request_name)
        if isinstance(parent_graph, ResnetGraph):
            my_session.load_children4(parent_graph._get_nodes(),depth=1)
            # loading only immediate children for all parents to not create shortcuts path in ontology graph
            my_session.close_connection()
        return my_session.Graph.ontology_graph()


    def load_ontology(self, parent_ontology_groups:list):
        """
        Input
        -----
        [List of Ontology Parent Names]
        
        Return
        ------
        self.child2parent = {child_name:[ontology_parent_names]}
        """
        self.child2parent = dict()
        self.add_ent_props(['Name'])
        for group_name in parent_ontology_groups:
            childs = self.child_graph([group_name],['Name'],include_parents=False)
            for id, child in self.Graph.nodes(data=True):
                child_name = child['Name'][0]
                try:
                    self.child2parent[child_name].append(group_name)
                except KeyError:
                    self.child2parent[child_name] = [group_name]


    def ontopaths2(self,name2childs:dict, ontology_depth:int):
        '''
        Input
        -----
        name2childs = {name:[PSObject]}, where name is the concept name with ontology children in [PSObject]
        name2childs is generated by SemanticSearch class

        Return
        ------
        name2paths = {name:[path]}, where path has format 
        '''
        list_of_objlists = list(name2childs.values())
        children = PSObject.unpack(list_of_objlists)
        ontology_graph = self.__ontology_graph(children,ontology_depth)

        name2paths = dict()
        path_sep = '->'
        for concept_name, childs in name2childs.items():
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
    def filename4(folder_or_object_named:str):
        new_str = list(folder_or_object_named)
        for i in range(0, len(new_str)):
            if new_str[i] in {'>', '<', '|','/'}:
                new_str[i] = '-'
            if new_str[i] in {':'}:
                new_str[i] = '_'
        return "".join(new_str)


    def __find_folder(self, folder_name:str, in_parent_folder=''):
        folder_name_ = self.filename4(folder_name)
        parent_name = self.filename4(in_parent_folder)
        longest_path = str()
        for dirpath, dirnames, filenames in os.walk(os.path.abspath(self.data_dir)):
            for dirname in dirnames:
                if dirname == folder_name_:
                    if len(dirpath+dirname)+1 > len(longest_path):
                        if in_parent_folder:
                            if parent_name in dirpath:
                                longest_path = dirpath+'/'+dirname +'/'
                        else:
                            longest_path = dirpath+'/'+dirname +'/'

        return longest_path

    
    def dump_path(self, of_folder:str, in2parent_folder=''):
        my_dir = self.__find_folder(of_folder,in2parent_folder)
        if not my_dir:
            if in2parent_folder:
                parent_path = self.__find_folder(in2parent_folder)
                if parent_path:
                    my_dir = parent_path+self.filename4(of_folder)+'/'
                else:
                    my_dir = self.data_dir+self.filename4(in2parent_folder)+'/'+self.filename4(of_folder)+'/'
            else:
                my_dir = self.data_dir+self.filename4(of_folder)+'/'
            try: 
               os.mkdir(my_dir)
            except FileExistsError: 
                pass
        return my_dir


    def __dump_base_name(self,folder_name:str):
        return 'content of '+ APISession.filename4(folder_name)+'_'

    def __dumpfiles(self, in_folder:str,in2parent_folder='', with_extension='rnef'):
        folder_path = self.dump_path(in_folder,in2parent_folder)
        listing = glob.glob(folder_path+'*.' + with_extension)
        return folder_path, self.__dump_base_name(in_folder), len(listing)


    def __make_dumpfile_name(self,of_folder:str,in2parent_folder='',new=False,filecount=0):
        '''
        Input
        -----
        new - forces to return name for new dumpfile
        filecount - if not zero forces to create dump file with index = filecount+1
        '''
        folder_path,of_folder_base_name,file_index = self.__dumpfiles(of_folder,in2parent_folder)
        if filecount:
            file_index = filecount
        if new:
            file_index += 1
    
        return folder_path + of_folder_base_name+str(file_index)+'.rnef'


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
    

    def rnefs2dump(self,rnef_xml:str,to_folder='',in2parent_folder='',can_close=True, lock=None):
        '''
        Dumps
        -----
        "rnef_xml" string into 'to_folder' inside 'in2parent_folder' located in "self.data_dir"
        if size of dump file exceeds "max_rnef_size", "rnef_xml" is splitted into several RNEF files\n
        dump RNEF files are named as: 'content of to_folder#',
        where # - dump file number
        keep can_close = False to continue dumping into one file
        '''
        if lock is None:
            lock = threading.Lock()

        with lock:
            write2 = self.__make_dumpfile_name(to_folder,in2parent_folder)
            if Path(write2).exists():
                file_size = os.path.getsize(write2)
                if file_size < self.max_rnef_size:
                    with open(write2,'a',encoding='utf-8') as f:
                        f.write(rnef_xml)
                        f.flush()
                        if os.path.getsize(write2) > self.max_rnef_size and can_close:
                            # need to create new dump file ASAP for next thread to write into it
                            new_write2 = self.__make_dumpfile_name(to_folder,in2parent_folder,new=True)
                            with open(new_write2,'w',encoding='utf-8') as new_f:
                                new_f.write('<batch>\n')
                                new_f.flush()

                            f.write('</batch>')
                        f.flush()
                else:
                    # file_size >= self.max_rnef_size
                    if can_close:
                        # need to create new dump file ASAP for next thread to write into it
                        new_write2 = self.__make_dumpfile_name(to_folder,in2parent_folder,new=True)
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
                write2 = self.__make_dumpfile_name(to_folder,in2parent_folder,new=True)
                with open(write2,'w',encoding='utf-8') as new_f:
                    new_f.write('<batch>\n')
                    new_f.write(rnef_xml)
                    new_f.flush()


    def _dump2rnef(self,graph=ResnetGraph(),to_folder='',in_parent_folder='',can_close=True,lock=None):
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
        my_graph = graph.copy() if graph.number_of_nodes()>0 else self.Graph.copy() 
        # making input graph copy to release it for multithreading
        section_rels = set()
        for regulatorID, targetID, e in my_graph.edges(data='relation'):
            section_rels.add(e)
            if len(section_rels) == self.resnet_size:
                resnet_section = my_graph.subgraph_by_rels(list(section_rels))
                rnef_str = resnet_section.to_rnefstr(ent_props=self.entProps,rel_props=self.relProps)
                rnef_str = self.pretty_xml(rnef_str,remove_declaration=True)
                # dumps section
                self.rnefs2dump(rnef_str,to_folder,in_parent_folder,can_close,lock)
                resnet_section.clear_resnetgraph()
                section_rels.clear()
    
        # dumps leftover resnet_section with size < self.resnet_size
        resnet_section = my_graph.subgraph_by_rels(list(section_rels))
        rnef_str = resnet_section.to_rnefstr(ent_props=self.entProps,rel_props=self.relProps)
        rnef_str = self.pretty_xml(rnef_str,remove_declaration=True)
        self.rnefs2dump(rnef_str,to_folder,in_parent_folder,can_close,lock)

        if my_graph.number_of_edges() == 0:
            rnef_str = my_graph.to_rnefstr(ent_props=self.entProps,rel_props=self.relProps)
            rnef_str = self.pretty_xml(rnef_str,remove_declaration=True)
            self.rnefs2dump(rnef_str,to_folder,in_parent_folder,can_close,lock)

        if not self.no_mess:
            print('RNEF dump of "%s" graph into %s folder was done in %s' % 
                  (my_graph.name,to_folder,self.execution_time(dump_start)[0]),flush=True)
            
        return time.time()-dump_start


    def close_rnef_dump(self,for_folder='',in_parent_folder=''):
        parent_folder = '' if for_folder == in_parent_folder else in_parent_folder
        last_dump_file = self.__make_dumpfile_name(for_folder,parent_folder)
        f = open(last_dump_file,'a',encoding='utf-8')
        f.write('</batch>')
        f.close()


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
        resume_pos = resume_page*self.PageSize #self.add2self=False
        iterations_graph = self.__init_session(resume_pos)
        self.clear() # to save RAM self.Graph is cleared
        number_of_iterations = int(self.ResultSize / self.PageSize)
        if number_of_iterations:
            print('\n\nrequest "%s" found %d %s. Begin download in %d %d-thread iterations' % 
                (request_name,self.ResultSize,return_type, number_of_iterations/threads,threads))
            with ThreadPoolExecutor(1,f'{request_name} Dump{number_of_iterations}iters') as e:
                # max_workers = 1 otherwise _dump2rnef gets locked
                # dump_futures = list()
                for i in range(0,number_of_iterations,threads):
                    iterations_graph = iterations_graph.compose(self.__thread__(threads))
                    e.submit(self._dump2rnef, iterations_graph.copy(), self.dump_folder,'',True,lock)
                    page_ref_count = iterations_graph.weight()
                    reference_counter += page_ref_count
                    remaining_iterations = number_of_iterations-i-threads
                    exec_time, remaining_time = self.execution_time(start_time,remaining_iterations,number_of_iterations) 
                    print("With %d in %d iterations, %d %s out of %d results with %d references retrieved in %s using %d threads" % 
                    (i+threads,number_of_iterations,self.ResultPos,return_type,self.ResultSize,reference_counter,exec_time,threads))
                    print('Estimated remaining retrieval time %s: '% remaining_time)
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
        self.clone = False
        return ResnetGraph() # fake return for consistency with process_oql
    

    def common_neighbors(self,with_entity_types:list,of_dbids1:list,reltypes12:list,effect12:list,dir12:str,
                                and_dbids3:list,reltypes23:list,effect23:list,dir23:str):
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
        sel_rels12 = sel_rels+' AND Effect = ({})'.format(','.join(effect12)) if effect12 else sel_rels
        sel_rels23 = sel_rels+' AND Effect = ({})'.format(','.join(effect23)) if effect23 else sel_rels

        sel_rels12 = sel_rels12.format(','.join(reltypes12))
        sel_rels23 = sel_rels23.format(','.join(reltypes23))

        find_neighbors_oql = r'SELECT Entity objectType = ('+','.join(with_entity_types)+')'

        find_neighbors_oql = find_neighbors_oql+' AND Connected by ('+sel_rels12+') to (SELECT Entity WHERE id = ({ids1}))'
        find_neighbors_oql += ' AND Connected by ('+sel_rels23+') to (SELECT Entity WHERE id = ({ids2}))'

        r_n = f'Find entities common between {len(of_dbids1)} and {len(and_dbids3)} entities' 
        neighbors_entity_graph = self.iterate_oql2(f'{find_neighbors_oql}',of_dbids1,and_dbids3,request_name=r_n)
        neighbors_dbids = list(neighbors_entity_graph.dbids4nodes())
        print('Found %d common neighbors' % len(neighbors_dbids))

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
            request_name = f'{str(int(i/1000)+1)} iteration out of {str(number_of_iterations)} to find genes linked to GVs'
            gvs = self.process_oql(oql_query,request_name)
            if isinstance(gvs,ResnetGraph):
                prot2gvs_graph.add_graph(gvs)

        # making gvid2genes for subsequent annotation
        gvid2genes = dict()
        for gv_id, protein_id, rel in prot2gvs_graph.edges.data('relation'):
            protein_node = prot2gvs_graph._psobj(protein_id)
            gene_name = protein_node['Name'][0]
            try:
                gvid2genes[gv_id].append(gene_name)
            except KeyError:
                gvid2genes[gv_id] = [gene_name]

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
        my_psobjs = self.iterate_oql(oql_query,propValues,request_name=rn,step=950)._get_nodes()
        # step size is reduced to 950 to mitigate potential excessive length of oql_query string due to zeep package limitations

        if get_childs:
            children, parents_with_children = self.load_children4(my_psobjs,add2self)
            my_psobjs = self.Graph._get_nodes(ResnetGraph.uids(my_psobjs))+list(children)
        elif add2self:
            self.Graph.add_psobjs(set(my_psobjs))
        
        return my_psobjs

##################### EXPERIMENT EXPERIMENT EXPERIMENT EXPERIMENT ###########################
    def load_dbids4(self,psobjs:list):
        '''
        Return
        ------
        mapped_objs,no_dbid_objs - [PSObject],[PSObject]
        mapping is done first by Name then by URN
        '''
        print(f'Reterieving database identifiers for {len(psobjs)} entities using Name identifier')
        kwargs = {TO_RETRIEVE:NO_REL_PROPERTIES}
        my_session = self._clone_session(**kwargs)
        names = ResnetGraph.names(psobjs)
        db_nodes = self._props2psobj(names,['Name'],get_childs=False)
        name2dbid = {n.name():n.dbid() for n in db_nodes}
        mapped_objs = list()
        no_dbid_objs = list()
        for psobj in psobjs:
            try:
                my_dbid = name2dbid[psobj.name()]
                psobj.update_with_value(DBID,my_dbid)
                mapped_objs.append(psobj)
            except KeyError:
                no_dbid_objs.append(psobj)
                continue
        
        if no_dbid_objs:
            urns_need_mapping = ResnetGraph.urns(no_dbid_objs)
            db_nodes = self._props2psobj(urns_need_mapping,['URN'],get_childs=False)
            urn2dbid = {n.urn():n.dbid() for n in db_nodes}
            mapped_by_urn = list()
            for psobj in no_dbid_objs:
                try:
                    my_dbid = urn2dbid[psobj.urn()]
                    psobj.update_with_value(DBID,my_dbid)
                    mapped_by_urn.append(psobj)
                except KeyError:
                    continue

            mapped_objs += mapped_by_urn
            no_dbid_objs = [obj for obj in no_dbid_objs if obj not in mapped_by_urn]

        urn2dbids = {n.urn():n[DBID] for n in mapped_objs}
        self.Graph.set_node_annotation(urn2dbids,DBID)
        print(f'Loaded {len(mapped_objs)} database identitiers for {len(psobjs)} entities')
        my_session.close_connection()
        return mapped_objs,no_dbid_objs


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
    

    def map_graph(self,graph:ResnetGraph, map_by=['Name']):
        '''
        Input
        -----
        graph - ResnetGraph where nodes may have arbitrary URNs
        map_by - list of database properties for mapping. Nodes in "graph" must have same properties  

        Return
        ------
        "graph" copy with nodes with database URNs
        '''
        print(f'Mapping input graph with {len(graph)} entities using {map_by} identifiers')
        kwargs = {TO_RETRIEVE:NO_REL_PROPERTIES}
        my_session = self._clone_session(**kwargs)
        my_session.add_ent_props(map_by)

        graph_psobjs = graph.psobjs_with(map_by)
        obj_props = graph.props(map_by,graph_psobjs)
        my_session._props2psobj(list(obj_props),map_by,get_childs=False)
        props2objs,uid2propval  = my_session.Graph.props2obj_dict([],map_by,case_insensitive=True)

        return graph.remap_graph(props2objs,map_by)
    

    def get_ppi(self, interactors_dbids:set):
        '''
        Return
        ------
        ResnetGraph with relation types: Binding, DirectRegulation, ProtModification 
        '''
        splitter = list() #holds lists of ids splits 
        splitter.append(list(interactors_dbids))
        number_of_splits = int(math.log2(len(interactors_dbids)))
        ppi_keeper = ResnetGraph()
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




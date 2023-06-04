import time, math, os, glob
import networkx as nx
from zeep import exceptions
from xml.dom import minidom
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import timedelta

from .ZeepToNetworkx import PSNetworx, len
from .ResnetGraph import ResnetGraph,PSObject,PSRelation,df,REFCOUNT,CHILDS
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

SESSION_KWARGS = {'what2retrieve','connect2server','no_mess','data_dir',TO_RETRIEVE,'use_cache'}


class APISession(PSNetworx):
    '''
    manages graph retrieval from database and loading cache files
    '''
    pass
    APIconfig = dict()
    ResultRef = None
    ResultPos = 0
    ResultSize = None
    PageSize = 1000
    __IsOn1st_page = True # to print header in dump file
    GOQLquery = str()
    __getLinks = True

    relProps = ['URN'] # URNs are used as MultiDiGraph keys in self.Graph
    entProps = ['Name'] # Name is highly recommended to make sense of the data

    add2self = True # if False will not add new graph to self.Graph
    merge2self = True #if False will replace nodes and relations properties in self.Graph by properties from new graph
    
    use_cache = False # if True signals functions to use graph data from cache files instead of retrieving data from database using GOQL queries 
    # all cached networks are written into cache_dir as files named cache_name.rnef:
    cache_dir = 'ElsevierAPI/ResnetAPI/__pscache__/'
    reference_cache_size = 1000000 # max number of reference allowed in self.Graph. Clears self.Graph if exceeded
    resnet_size = 1000 # number of <node><control> sections in RNEF dump
    max_rnef_size = 100000000 # max size of RNEF XML dump file. If dump file exceeds max_file_size new file is opened with index++
    
    data_dir = ''
    DumpFiles = []
    sep = '\t'
    print_rel21row = False
    

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
        self.APIconfig= dict(args[0])
        my_kwargs = {'what2retrieve':NO_REL_PROPERTIES}
        my_kwargs.update(kwargs)
        super().__init__(self.APIconfig, **my_kwargs)
        self.__retrieve(my_kwargs['what2retrieve'])
        self.DumpFiles = []
        self.set_dir(kwargs.get('data_dir',''))
        self.use_cache = my_kwargs.get('use_cache',False)
    

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
        session_kwargs = {k:v for k,v in kwargs.items() if k in SESSION_KWARGS}
        other_parameters = {k:v for k,v in kwargs.items() if k not in SESSION_KWARGS}
        return session_kwargs,other_parameters


    def _clone_session(self,**kwargs):
        """
        used by self.process_oql to clone session 
        """
        my_kwargs = dict(kwargs)
        my_kwargs['what2retrieve'] = my_kwargs.pop(TO_RETRIEVE,CURRENT_SPECS)
        if my_kwargs['what2retrieve'] == CURRENT_SPECS:  
            new_session = APISession(self.APIconfig,**my_kwargs)
            new_session.entProps = list(self.entProps)
            new_session.relProps = list(self.relProps)
        else:
            new_session = APISession(self.APIconfig,**my_kwargs)
            new_session.entProps = list(self.entProps)
        
        new_session.data_dir = self.data_dir
        new_session.no_mess = my_kwargs.get('no_mess',self.no_mess)
        return new_session


    def __set_get_links(self,oql_query=''):
        my_query = oql_query if oql_query else self.GOQLquery
        return my_query[7:15] == 'Relation'


    def __replace_goql(self, oql_query:str):
        self.GOQLquery = oql_query
        self.__getLinks = self.__set_get_links()
        self.ResultRef = None


    def __add2self(self,other:'APISession'):
        '''
        needs relations to have dbid
        works only with another graph from database
        '''
        if self.add2self:
            if self.merge2self:
                self.Graph.add_graph(other.Graph)
            else:
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
    def __init_session(self, first_iteration=0):
        self.__set_get_links()
        obj_props = self.relProps if self.__getLinks else self.entProps
        zeep_data, (self.ResultRef, self.ResultSize, self.ResultPos) = self.init_session(self.GOQLquery,
                                                                                        self.PageSize,
                                                                                        obj_props,
                                                                                        getLinks = self.__getLinks)
        if first_iteration > 0:
            self.ResultPos = first_iteration
            print ('Resuming retrieval from %d position' % first_iteration)
            return self.__get_next_page(self.ResultPos)                                                  

        if type(zeep_data) != type(None):
            if self.__getLinks:
                obj_dbids = list(set([x['EntityId'] for x in zeep_data.Links.Link]))
                zeep_objects = self.get_object_properties(obj_dbids, self.entProps)
                return self._load_graph(zeep_data, zeep_objects,self.add2self,self.merge2self)
            else:
                return self._load_graph(None, zeep_data,self.add2self,self.merge2self)
        else:
            return ResnetGraph()

        
    def __get_next_page(self,result_pos:int):
        '''
        Input
        -----
        result_pos - resut position to begin download.  Has to be provided explicitly to initiate multithreading
        '''
        # current_pos does not change during multithreading initiation!!!!  
        if self.ResultPos < self.ResultSize:
            obj_props = self.relProps if self.__getLinks else self.entProps
            zeep_data, self.ResultSize, current_pos = self.get_session_page(self.ResultRef, result_pos, self.PageSize,
                                                                       self.ResultSize,obj_props,getLinks=self.__getLinks)                                                          
            if type(zeep_data) != type(None):
                if self.__getLinks and len(zeep_data.Links.Link) > 0:
                    obj_dbids = list(set([x['EntityId'] for x in zeep_data.Links.Link]))
                    zeep_objects = self.get_object_properties(obj_dbids, self.entProps)
                    return self._load_graph(zeep_data, zeep_objects,self.add2self,self.merge2self)
                else:
                    return self._load_graph(None, zeep_data,self.add2self,self.merge2self)
            else:
                return ResnetGraph()
        else: return ResnetGraph()


    def __thread__(self,number_of_iterations:int,max_threads=0,process_name='process_oql'):
        max_workers = max_threads if max_threads else number_of_iterations
        if not self.no_mess:
            print(f'Retrieval starts from {self.ResultPos} result')
            remaining_retrieval = self.ResultSize - self.ResultPos
            next_download_size = min(number_of_iterations*self.PageSize, remaining_retrieval)
            print(f'Begin retrieving next {next_download_size} results in {max_workers} threads',flush=True)
        
        entire_graph = ResnetGraph()
        result_pos = self.ResultPos
        with ThreadPoolExecutor(max_workers, thread_name_prefix=process_name) as e:
            futures = list()
            for i in range(0,number_of_iterations):          
                futures.append(e.submit(self.__get_next_page,result_pos))
                result_pos += self.PageSize
            
            for future in as_completed(futures):
                try:
                    entire_graph = entire_graph.compose(future.result())
                except exceptions.TransportError:
                    raise exceptions.TransportError
            
        self.ResultPos = result_pos
        return entire_graph


    def process_oql(self, oql_query, request_name='') -> ResnetGraph:
        self.__replace_goql(oql_query)
        start_time = time.time()
        return_type = 'relations' if self.__getLinks else 'entities'
        entire_graph = self.__init_session()
        number_of_iterations = int(self.ResultSize / self.PageSize)
        my_request_name = request_name if request_name else self.GOQLquery[:100]+'...'
    
        if number_of_iterations:
            if not self.no_mess:
                print('\n\"%s\"\nrequest found %d %s.\n%d is retrieved. Remaining %d results will be retrieved in %d parallel iterations' % 
            (my_request_name,self.ResultSize,return_type,self.ResultPos,(self.ResultSize-self.PageSize),number_of_iterations))
            try:
                iterations_graph = self.__thread__(number_of_iterations,process_name=request_name)
                entire_graph = entire_graph.compose(iterations_graph)
            except exceptions.TransportError:
                raise exceptions.TransportError('Table lock detected!!! Aborting operation!!!')

        if not self.no_mess:
            print('"%s"\nretrieved %d nodes and %d edges in %s by %d parallel iterations' % 
            (my_request_name, entire_graph.number_of_nodes(), entire_graph.number_of_edges(),
                    self.execution_time(start_time), number_of_iterations+1),flush=True)

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
        req_name = f'Retrieving annotations for "{request_name}"'
        relation_dbids2return = id_only_graph._relations_dbids()
        node_dbids2return = id_only_graph.dbids4nodes()
        
        rels4subgraph = [rel for dbid,rel in self.dbid2relation.items() if dbid in relation_dbids2return]
        return_subgraph = self.Graph.subgraph_by_rels(rels4subgraph)
        return_subgraph.add_psobjs(self.Graph.psobj_with_dbids(node_dbids2return))
        
        if use_cache:
            relation_dbids2retreive = relation_dbids2return.difference(self.dbid2relation.keys())
            nodes_dbids2retreive = set(node_dbids2return).difference(set(self.Graph.dbids4nodes()))
        else:
            relation_dbids2retreive = set(relation_dbids2return)
            nodes_dbids2retreive = set(node_dbids2return)

        if relation_dbids2retreive or nodes_dbids2retreive:
            exist_nodes = len(node_dbids2return)-len(nodes_dbids2retreive)
            exist_rels = len(relation_dbids2return)-len(relation_dbids2retreive)
            print(f'{exist_nodes} nodes and {exist_rels} relations were downloaded from database previously')
            print('%d nodes and %d relations will be loaded from database' % 
                        (len(nodes_dbids2retreive),len(relation_dbids2retreive)))
        
        add2return = ResnetGraph()
        if relation_dbids2retreive:
            rel_query = 'SELECT Relation WHERE id = ({ids})'
            rels_add2return = self.__iterate__(rel_query,relation_dbids2retreive,req_name,step=1000)
            add2return = rels_add2return
            nodes_dbids2retreive = nodes_dbids2retreive.difference(add2return.dbids4nodes())
            
        if nodes_dbids2retreive:
            entity_query = 'SELECT Entity WHERE id = ({ids})'
            nodes_add2return = self.__iterate__(entity_query,nodes_dbids2retreive,req_name,step=1000)
            add2return = add2return.compose(nodes_add2return)
        
        return_subgraph = return_subgraph.compose(add2return)
        assert(return_subgraph.number_of_nodes() == len(node_dbids2return))
       # assert(return_subgraph.number_of_edges() >= len(relation_dbids2return))
       # privately owned relations are not returned by API. 
       # therefore, return_subgraph.number_of_edges() may have less relations than "relation_dbids2return"
       # 
       # return_subgraph may also contain more relations than "relation_dbids2return" 
       # due to duplication of non-directional relations to both directions
        return return_subgraph
 

    def __iterate__(self,oql_query:str,dbids:set,request_name='',step=1000):
        '''
        Input
        -----
        oql_query MUST contain string placeholder called {ids} to iterate dbids\n
        oql_query MUST contain string placeholder called {props} to iterate other database properties
        "step" - controls duration of request to avoid timeout during multithreading

        uses self.entProp and self.relProp to retrieve properties\n
        use self.add2self and self.merge2self to control caching
        '''
        if not dbids: return ResnetGraph()
        if oql_query.find('{ids',20) > 0:
            my_oql = oql_query
            def join_list (l:list):
                return ','.join(map(str,l))
        else:
            my_oql = oql_query.replace('{props','{ids')
            def join_list (l:list):
                return OQL.join_with_quotes(l) 
            
        number_of_iterations = math.ceil(len(dbids)/step)  
        # new sessions must be created to avoid ResultRef pointer overwriting
        future_sessions = list()
        entire_graph = ResnetGraph()
        thread_name = f'Iterate_{len(dbids)}_ids'
        print(f'Retrieval of {len(dbids)} objects will be done in {number_of_iterations} parallel iterations' )
        dbids_list = list(dbids)
        with ThreadPoolExecutor(max_workers=number_of_iterations, thread_name_prefix=thread_name) as e:   
            futures = list()
            for i in range(0,len(dbids_list), step):
                iter_dbids = dbids_list[i:i+step]
                oql_query_with_dbids = my_oql.format(ids=join_list(iter_dbids))
                iter_name = f'Iteration #{int(i/step)+1} out of {number_of_iterations} for "{request_name}" retrieving {len(iter_dbids)} ids'
                future_session = self._clone_session()
                futures.append(e.submit(future_session.process_oql,oql_query_with_dbids,iter_name))
                future_sessions.append(future_session)
            
            for future in as_completed(futures):
                entire_graph = entire_graph.compose(future.result())

        for session in future_sessions:
            self.__add2self(session)
            session.close_connection()
            
        return entire_graph 


    def __iterate_oql__(self,oql_query:str,ids_or_props:set,request_name='',step=1000):
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

        dbidonly_graph = self.__iterate__(oql_query,ids_or_props,req_name,step)

        self.entProps = my_ent_props
        self.relProps = my_rel_props
        self.add2self = old_add2self 
        return dbidonly_graph

    '''
    def __iterate2_singlethread_(self, oql_query:str, ids1:list, ids2:list,request_name='',threads=1):
        """
        Input
        -----
        oql_query MUST contain string placeholder called {ids} to iterate dbids\n
        oql_query MUST contain string placeholder called {props} to iterate other database properties
        """
        if oql_query.find('{ids',20) > 0:
            step = 1000
            my_oql = oql_query
            def join_list (l:list):
                return ','.join(map(str,l))
        else:
            step = 950 # string properties can for too long oql queries
            my_oql = oql_query.replace('{props','{ids')
            def join_list (l:list):
                return OQL.join_with_quotes(l)
                             
        number_of_iterations = math.ceil(len(ids1)/step) * math.ceil(len(ids2)/step)
        print('Connecting %d with %d entities' % (len(ids1), len(ids2)))
        if number_of_iterations > 2:
            print('Query will be executed in %d iterations' % number_of_iterations)
        entire_id_graph = ResnetGraph()
        iteration_counter = 1
        start  = time.time()
        for i1 in range(0,len(ids1), step):
            iter_ids1 = ids1[i1:i1+step]
            for i2 in range(0,len(ids2), step):
                iter_ids2 = ids2[i2:i2+step]
                oql_query_with_dbids = my_oql.format(ids1=join_list(iter_ids1),ids2=join_list(iter_ids2))
                iter_graph = self.process_oql(oql_query_with_dbids,request_name)
                entire_id_graph.add_graph(iter_graph)
                if number_of_iterations > 2:
                    print('Iteration %d out of %d performed in %s' % 
                        (iteration_counter,number_of_iterations,self.execution_time(start)))
                iteration_counter +=1

        return entire_id_graph
    '''

    def __iterate2__(self, oql_query:str, ids1:list, ids2:list,request_name='iteration1',step=500):
        '''
        Input
        -----
        oql_query MUST contain string placeholder called {ids} to iterate dbids\n
        oql_query MUST contain string placeholder called {props} to iterate other database properties
        '''
        
        if oql_query.find('{ids',20) > 0:
            id_oql = oql_query
            hint = r'{ids}'
            iteration1step = min(step, 1000)
            iteration2step = min(2*iteration1step, 1000)
            def join_list (l:list):
                return ','.join(map(str,l))
        else:
            id_oql = oql_query.replace('{props','{ids')
            hint = r'{props}'
            iteration1step = min(step,950) # string properties can form long oql queries exceeding zeep limits
            iteration2step = min(2*iteration1step, 950)
            def join_list (l:list):
                return OQL.join_with_quotes(l)
                             
        number_of_iterations = math.ceil(len(ids1)/step) * math.ceil(len(ids2)/iteration1step)
        if not number_of_iterations: return ResnetGraph()
        print(f'Iterating {len(ids1)} entities with {len(ids2)} entities for {request_name}')
        if number_of_iterations > 2:
            print('Query will be executed in %d iterations in parallel threads' % number_of_iterations)

        entire_graph = ResnetGraph()
        
        if len(ids1) <= len(ids2):
            number_of_iterations = math.ceil(len(ids1)/iteration1step)
            for i1 in range(0,len(ids1), iteration1step):
                iter_ids1 = ids1[i1:i1+iteration1step]
                iter_oql_query = id_oql.format(ids1=join_list(iter_ids1),ids2=hint)
                iter_graph = self.__iterate__(iter_oql_query,ids2,'iteration2',iteration2step)
                entire_graph = entire_graph.compose(iter_graph)
                print(f'Processed {min([i1+iteration1step,len(ids1)])} out of {len(ids1)} identifiers against {len(ids2)} identifiers')
        else:
            number_of_iterations = math.ceil(len(ids2)/iteration1step)
            for i2 in range(0,len(ids2), iteration1step):
                iter_ids2 = ids2[i2:i2+iteration1step]
                iter_oql_query = id_oql.format(ids1=hint,ids2=join_list(iter_ids2))
                iter_graph = self.__iterate__(iter_oql_query,ids1,'iteration2',iteration2step)
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

        entire_graph = self.__iterate2__(oql_query,list(ids_or_props1),list(ids_or_props2),request_name,step)

        self.entProps = my_ent_props
        self.relProps = my_rel_props
        self.add2self = old_add2self
        return entire_graph

 
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
   

    def iterate_oql2(self,oql_query:str,dbid_set1:set,dbid_set2:set,use_cache=True,request_name='',step=500):
        """
        # oql_query MUST contain 2 string placeholders called {ids1} and {ids2}
        # oql_query MUST contain 2 string placeholder called {props1} and {props2}\n
        if iterable id_sets contain property values other than database id
        """
        print('Processing "%s" request' % request_name)
        dbid_only_graph = self.__iterate_oql2__(oql_query,dbid_set1,dbid_set2,request_name,step)
        
        if dbid_only_graph:
            print('Will retrieve annotated graph with %d entities and %d relations' 
                        %(dbid_only_graph.number_of_nodes(),dbid_only_graph.number_of_edges()))
            annoated_graph = self.__annotate_dbid_graph__(dbid_only_graph,use_cache,request_name)
            return annoated_graph
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
    def _get_saved_results(self, results_names:list, get_links=True):

        #self = self.__clone() if self.clone else self

        query_names = OQL.join_with_quotes(results_names)
        oql_query = 'SELECT Result WHERE Name = ({names})'
        oql_query = oql_query.format(names=query_names)
        result_graph = self.load_graph_from_oql(oql_query,self.relProps,
                                                self.entProps,get_links=False)
        self.ResultPos = 0
        self.ResultSize  = 1
        self.PageSize  = 10000
        entire_graph = ResnetGraph()
        start = time.time()

        self.__getLinks = get_links
        for id, name in result_graph.nodes(data='Name'):
            self.ResultRef  = 'ID='+str(id)
            while self.ResultPos < self.ResultSize:
                page_graph = self.__get_next_page(self.ResultPos)
                entire_graph = nx.compose(page_graph, entire_graph)
                iteration = math.ceil(self.ResultPos / self.PageSize)
                total_iterations = math.ceil(self.ResultSize/self.PageSize)
                print ('%d out of %d iterations fetching %d objects from %s search results performed in %s' % 
                (iteration, total_iterations, self.ResultSize, name[0], self.execution_time(start)[0]))

        return entire_graph


    def _get_ppi_graph(self,fout):
        # Build PPI network between proteins
        graph_proteins = self.Graph.dbids4nodes(['Protein'])
        print('Retrieving PPI network for %d proteins' % (len(graph_proteins)))
        start_time = time.time()
        ppi_graph = self.get_ppi(graph_proteins, self.relProps, self.entProps)
        self.to_csv(self.data_dir+fout, in_graph=ppi_graph)
        print("PPI network with %d relations was retrieved in %s ---" % (
            ppi_graph.size(), self.execution_time(start_time)[0]))

        return ppi_graph


    def _get_network_graph(self, for_interactor_dbids:set, connect_by_rel_types:list=None):
        #self = self.__clone() if self.clone else self
        network2return = self.get_network(for_interactor_dbids,connect_by_rel_types,
                                            self.relProps,self.entProps)
        return network2return


    def connect_nodes(self,node_dbids1:set,node_dbids2:set,
                        by_relation_type=[],with_effect=[],in_direction='',
                        use_relation_cache=True, step=500):
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
        return self.iterate_oql2(f'{oql_query}',node_dbids1,node_dbids2,use_relation_cache,request_name=req_name,step=step)
    # step is set 250 to mitigate possible timeout during multithreading


    def get_group_members(self, group_names:list):
        groups = OQL.join_with_quotes(group_names)
        oql_query = f'SELECT Entity WHERE MemberOf (SELECT Group WHERE Name = ({groups}))'
        req_name = 'Find members of groups: ' + ','.join(group_names)
        graph2return = self.process_oql(oql_query, request_name=req_name)
        if len(graph2return) == 0:
            print('%s groups are empty or do not exist in databse' % str(group_names))
        else:
            print('loaded %d members from %s' % (graph2return.number_of_nodes(),str(group_names)))
        return graph2return
    

    def add_group_annotation(self,group_names:list,graph=ResnetGraph()):
        my_graph = graph if graph else self.Graph
        urns2value = PSObject()
        for group_name in group_names:
            group_graph = self.get_group_members([group_name])
            group_members = group_graph._get_nodes()
            [urns2value.append_property(o.urn(),group_name) for o in group_members]

        my_graph.set_node_annotation(urns2value,BELONGS2GROUPS)


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

        p2obj,id2p = self.Graph.props2obj_dict(need_db_mapping,in_properties,case_insensitive)
        propval2objs.update(p2obj)
        for id,prop_vals in id2p.items():
            try:
                mapped_values = set(objid2propval[id])
                mapped_values.update(map(lambda x:str(x).lower(),prop_vals))
                objid2propval[id] = list(mapped_values)
            except KeyError:
                objid2propval[id] = prop_vals

        print("%d entities were mapped to %d attributes in %s properties" % 
                (len(objid2propval), len(using_values),','.join(in_properties)))

        self.entProps = current_ent_props
        return propval2objs, objid2propval

################################# ONTOLOGY  ONTOLOGY ############################################
    def __load_children(self,parent:PSObject,min_connectivity=0):
            parent_dbid = parent.dbid()
            if parent_dbid:
                query_ontology = OQL.get_childs([parent_dbid],['id'])
                if min_connectivity:
                    query_ontology += f' AND Connectivity >= {min_connectivity}'

                my_session = self._clone_session()
                children_graph = my_session.process_oql(query_ontology,f'Retrieve ontology children for {parent.name()}')
                my_session.close_connection()
                
                return parent, children_graph._get_nodes()
            else:
                print(f'Parent {parent.name()} has no database identifier. Its children cannot be loaded')
                return parent,[]
        

    def load_children4(self,parents:list,add2self=True,min_connectivity=0):
        '''
        Input
        -----
        parents - [PSObject]

        Return
        ------
        {PSObject} - set of all found children\n
    
        Updates
        -------
        parents in self.Graph with CHILDS properties as [PSObject] children
        '''
        process_start = time.time()
        if parents:
            annotate_parents = list(set(parents))
            if add2self:
                self.Graph.add_psobjs(parents)    
        else:
            annotate_parents = self.Graph._get_nodes()

        print(f'Loading ontology children for {len(annotate_parents)} entities')

        max_threads = 50 
        # 4 threads perform a bit faster than 8 threads 
        # 288 parents out of 288 were processed in 0:19:51.732684 by 8 threads
        # 288 parents out of 288 were processed in 0:18:52.715937 by 4 threads
        max_threaded_time = 0
        need_children = list()
        child_counter = set()
        for p, parent in enumerate(annotate_parents):
            try:
                children = self.Graph.nodes[parent.uid()][CHILDS]
                child_counter.update(children)
            except KeyError:
                need_children.append(parent)
                if len(need_children) >= max_threads or p == len(annotate_parents)-1:       
                    thread_start = time.time()
                    with ThreadPoolExecutor(max_workers=len(need_children), thread_name_prefix='load_children') as e:
                        futures = list()
                        [futures.append(e.submit(self.__load_children, nc)) for nc in need_children]

                        for future in as_completed(futures):
                            f_parent, children = future.result()
                            self.Graph.add_psobjs(children)
                            nx.function.set_node_attributes(self.Graph, {f_parent.uid():{CHILDS:children}})
                            child_counter.update(children)
                            # set_node_attributes() does not update nodes that do not exist in Graph
       
                    threaded_time = time.time() - thread_start
                    if threaded_time > max_threaded_time:
                        max_threaded_time = threaded_time
                        if max_threaded_time > 240:
                            print(f'Longest {max_threads}-threaded time {"{}".format(str(timedelta(seconds=max_threaded_time)))} \
                                  is close to Apache default 5 minutes transaction timeout !!!')
                    if not self.no_mess:
                        print (f'{p+1} parents out of {len(annotate_parents)} were processed in {self.execution_time(process_start)}')
                    need_children.clear()
        
        print(f'{len(child_counter)} children for {len(annotate_parents)} parent entities were found in database \
              in {self.execution_time(process_start)}')
        print(f'Longest {max_threads}-threaded time: {"{}".format(str(timedelta(seconds=max_threaded_time)))}')
        return child_counter


    def child_graph(self, propValues:list, search_by_properties=[],include_parents=True):
        if not search_by_properties: search_by_properties = ['Name','Alias']
        oql_query = OQL.get_childs(propValues,search_by_properties,include_parents=include_parents)
        prop_val_str = ','.join(propValues)
        request_name = f'Find ontology children for {prop_val_str}'
        ontology_graph = self.process_oql(oql_query,request_name)
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
        ontology_graph.add_psobjs(members)
        for m in members:
            child_id = m.uid()
            rel = PSRelation({'ObjTypeName':['MemberOf'],'Relationship':['is-a'],'Ontology':['Pathway Studio Ontology']})
            rel.Nodes['Regulators'] = [(child_id,0,0)]
            rel.Nodes['Targets'] = [(parent_id,1,0)]
            rel.append_property('Id',child_id)
            rel.append_property('URN', child_id.urn()) #fake URN
            ontology_graph.add_edge(child_id,parent_id,relation=rel)

        self.add_rel_props(['Relationship','Ontology'])
        return ontology_graph


    def add_parents(self,for_childs=[],depth=1):
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

        get_parent_query = 'SELECT Entity WHERE InOntology (SELECT Annotation WHERE Ontology=\'Pathway Studio Ontology\' AND Relationship=\'is-a\') inRange {steps} over (SELECT OntologicalNode WHERE id = ({ids}))'
        oql_query = get_parent_query.format(steps=str(depth),ids=','.join(map(str,child_dbids)))
        request_name = f'Find parents of {len(my_childs)} nodes with depth {depth}'
        parent_graph = self.process_oql(oql_query,request_name)
        
        self.load_children4(parent_graph._get_nodes())
        return self.Graph._get_nodes(list(parent_graph))


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


###################################  WRITE DUMP  CACHE, WRITE, DUMP  CACHE, WRITE DUMP  CACHE ##############################################
    def to_csv(self, file_out, in_graph: ResnetGraph=None, access_mode='w'):
        debug = not self.no_mess
        if not isinstance(in_graph,ResnetGraph): in_graph = self.Graph
        in_graph.print_references(file_out, self.relProps, self.entProps, access_mode, printHeader=self.__IsOn1st_page,
                                col_sep=self.sep,debug=debug,single_rel_row=self.print_rel21row)


    def to_pandas (self, in_graph=None, RefNumPrintLimit=0)-> 'df':
        if not isinstance(in_graph, ResnetGraph): in_graph = self.Graph
        return df(in_graph.ref2pandas(self.relProps,self.entProps,RefNumPrintLimit,single_rel_row=self.print_rel21row))


    @staticmethod
    def pretty_xml(xml_string:str, no_declaration = False):
        pretty_xml = str(minidom.parseString(xml_string).toprettyxml(indent='   '))
        if no_declaration:
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
            try: os.mkdir(my_dir) 
            except FileExistsError: pass
        return my_dir


    def __count_dumpfiles(self, in_folder:str,in2parent_folder='', with_extension='rnef'):
        folder_path = self.dump_path(in_folder,in2parent_folder)
        listing = glob.glob(folder_path+'*.' + with_extension)
        return folder_path,len(listing)


    def __dump_base_name(self,of_folder:str):
        return 'content of '+ APISession.filename4(of_folder)+'_'


    def __make_dumpfile_name(self,of_folder:str,in2parent_folder='',new=False,filecount=0):
        folder_path,file_count = self.__count_dumpfiles(of_folder,in2parent_folder)
        file_count += int(new)
        if filecount: file_count = filecount
        return folder_path + self.__dump_base_name(of_folder)+str(file_count)+'.rnef'


    def _2rnefs(self,graph=ResnetGraph(),add_rel_props:dict={},add_pathway_props:dict={}):
        '''
        used for dumping pathway, group and results objects by FolderContent
        Returns
        -------
        graph RNEF XML with single <resnet> section and session properties for nodes and edges
        '''
        my_graph = graph if graph else self.Graph
        rel_props = [p for p in self.relProps if p not in NO_RNEF_REL_PROPS]
        return my_graph.to_rnefstr(self.entProps,rel_props,add_rel_props,add_pathway_props)
    

    def rnefs2dump(self,rnef_xml:str,to_folder='',in2parent_folder='',can_close=True):
        '''
        Dumps
        -----
        "rnef_xml" string into 'to_folder' inside 'in2parent_folder' located in "self.data_dir"
        if size of dump file exceeds "max_rnef_size", "rnef_xml" is splitted into several RNEF files\n
        dump RNEF files are named as: 'content of to_folder#',
        where # - dump file number
        keep can_close = False to continue dumping
        '''
        write2 = self.__make_dumpfile_name(to_folder,in2parent_folder)
        if Path(write2).exists():
            f = open(write2,'a',encoding='utf-8')
            f.write(rnef_xml)
            file_size = os.path.getsize(write2)
            if file_size > self.max_rnef_size and can_close:
                f.write('</batch>')
                f.close()

                new_write2 = self.__make_dumpfile_name(to_folder,in2parent_folder,new=True)
                f = open(new_write2,'w',encoding='utf-8')
                f.write('<batch>\n')
            f.close()
        else:
            write2 = self.__make_dumpfile_name(to_folder,in2parent_folder,new=True)
            with open(write2,'w',encoding='utf-8') as f:
                f.write('<batch>\n'+rnef_xml)


    def __dump2rnef(self,graph=ResnetGraph(),to_folder='',in_parent_folder='',can_close=True):
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
        '''
        dump_start = time.time()
        my_graph = graph.copy() if graph else self.Graph.copy() 
        # making input graph copy to release it for multithreading
        section_rels = set()
        for regulatorID, targetID, e in my_graph.edges(data='relation'):
            section_rels.add(e)
            if len(section_rels) == self.resnet_size:
                resnet_section = my_graph.subgraph_by_rels(section_rels)
                rnef_str = resnet_section.to_rnefstr(ent_props=self.entProps,rel_props=self.relProps)
                rnef_str = self.pretty_xml(rnef_str,no_declaration=True)
                # dumps section
                self.rnefs2dump(rnef_str,to_folder,in_parent_folder,can_close)
                resnet_section.clear_resnetgraph()
                section_rels.clear()
    
        # dumps leftover resnet_section with size < self.resnet_size
        resnet_section = my_graph.subgraph_by_rels(section_rels)
        rnef_str = resnet_section.to_rnefstr(ent_props=self.entProps,rel_props=self.relProps)
        rnef_str = self.pretty_xml(rnef_str,no_declaration=True)
        self.rnefs2dump(rnef_str,to_folder,in_parent_folder,can_close)

        if my_graph.number_of_edges() == 0:
            rnef_str = my_graph.to_rnefstr(ent_props=self.entProps,rel_props=self.relProps)
            rnef_str = self.pretty_xml(rnef_str,no_declaration=True)
            self.rnefs2dump(rnef_str,to_folder,in_parent_folder,can_close)

        if not self.no_mess:
            print('RNEF dump of "%s" graph into %s folder was done in %s' % 
                  (my_graph.name,to_folder,self.execution_time(dump_start)),flush=True)
            
        return time.time()-dump_start


    def close_rnef_dump(self,for_folder='',in_parent_folder=''):
        parent_folder = '' if for_folder == in_parent_folder else in_parent_folder
        last_dump_file = self.__make_dumpfile_name(for_folder,parent_folder)
        f = open(last_dump_file,'a',encoding='utf-8')
        f.write('</batch>')
        f.close()


    def download_oql(self,oql_query,request_name:str,resume_page=0,threads=25):
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
        results of oql_query to folder in self.data_dir named "request_name.rnef".\n
        Splits dump into small files smaller than self.max_rnef_size 
        '''
        self.__replace_goql(oql_query)
        reference_counter = 0
        start_time = time.time()
        return_type = 'relations' if self.__getLinks else 'entities'
        resume_pos = resume_page*self.PageSize
        iterations_graph = self.__init_session(resume_pos)
        self.clear()
        number_of_iterations = int(self.ResultSize / self.PageSize)
        if number_of_iterations:
            print('\n\nrequest "%s" found %d %s. Begin download in %d %d-thread iterations' % 
                (request_name,self.ResultSize,return_type, number_of_iterations/threads,threads))
            for i in range(0,number_of_iterations,threads):
                iterations_graph = iterations_graph.compose(self.__thread__(threads))
                self.__dump2rnef(iterations_graph, request_name)
                if not self.no_mess:
                    page_ref_count = iterations_graph.weight()
                    reference_counter += page_ref_count
                    remaining_iterations = number_of_iterations-i-threads
                    exec_time, remaining_time = self.execution_time(start_time,remaining_iterations,number_of_iterations) 
                    print("With %d in %d iterations, %d %s out of %d results with %d references retrieved in %s using %d threads" % 
                    (i+threads,number_of_iterations,self.ResultPos,return_type,self.ResultSize,reference_counter,exec_time,threads))
                    print('Estimated remaining retrieval time %s: '% remaining_time)
                self.clear()
                iterations_graph.clear()
            
        self.close_rnef_dump(request_name)
        self.clear()
        self.ResultRef = ''
        self.ResultPos = 0
        self.ResultSize = 0
        self.__IsOn1st_page = True
        self.clone = False
        return

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
            future12 = executor.submit(self.iterate_oql2,f'{sel_12_graph_oql}',of_dbids1,neighbors_dbids, True,rn12) 
            future23 = executor.submit(self.iterate_oql2,f'{sel_23_graph_oql}',and_dbids3,neighbors_dbids, True,rn23)
            
            graph12 = future12.result()
            graph23 = future23.result()

        graph12.add_graph(graph23)
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
        neighbors_entity_graph = self.iterate_oql2(f'{find_neighbors_oql}',of_dbids1,and_dbids3,request_name=r_n)
        
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

        graph12 = self.iterate_oql2(sel_12_graph_oql,of_dbids1,neighbors_ids,True,rn12)
        graph23 = self.iterate_oql2(sel_23_graph_oql,and_dbids3,neighbors_ids,True,rn23)
        graph123 = graph12.compose(graph23)
        # to remove relations with no Effect annotation:
        rels_with_effect = graph123.psrels_with(['positive','negative'],['Effect']) 
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
            prot2gvs_graph.add_graph(self.process_oql(oql_query,request_name))

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
            children = self.load_children4(my_psobjs,add2self)
            my_psobjs = self.Graph._get_nodes(ResnetGraph.uids(my_psobjs))+list(children)
        elif add2self:
            self.Graph.add_psobjs(my_psobjs) 
        
        return my_psobjs

##################### EXPERIMENT EXPERIMENT EXPERIMENT EXPERIMENT ###########################
    def load_dbids4(self,psobjs:list):
        '''
        Return
        ------
        mapped_objs,no_dbid_objs - [PSObject],
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
                psobj.update_with_value('Id',my_dbid)
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
                    psobj.update_with_value('Id',my_dbid)
                    mapped_by_urn.append(psobj)
                except KeyError:
                    continue

            mapped_objs += mapped_by_urn
            no_dbid_objs = [obj for obj in no_dbid_objs if obj not in mapped_by_urn]

        urn2dbids = {n.urn():n['Id'] for n in mapped_objs}
        self.Graph.set_node_annotation(urn2dbids,'Id')
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
        graph_with_ci.add_psobjs(self.Graph._get_nodes())
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




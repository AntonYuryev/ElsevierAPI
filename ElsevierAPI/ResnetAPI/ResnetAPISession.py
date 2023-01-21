import time
import math
import os
import glob
import networkx as nx
from .ZeepToNetworkx import PSNetworx, len
from .ResnetGraph import ResnetGraph,PSObject,PSRelation,CHILDS,df,REFCOUNT
from .PSPathway import PSPathway
from .PathwayStudioGOQL import OQL
from .Zeep2Experiment import Experiment,Sample
from ..ETM_API.references import SENTENCE_PROPS,PS_BIBLIO_PROPS,PS_SENTENCE_PROPS,PS_REFIID_TYPES,RELATION_PROPS,ALL_PS_PROPS
from xml.dom import minidom
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

# options for "to_retrieve" parameter
NO_REL_PROPERTIES = -2
CURRENT_SPECS = -1 # keeps current retreival properties
DO_NOT_CLONE = 0
DATABASE_REFCOUNT_ONLY = 1 #retrieves only RelationNumberOfReferences
REFERENCE_IDENTIFIERS = 2 # retrieves only PS_ID_TYPES to count reference in runtime
BIBLIO_PROPERTIES = 3 # retrieves only PS_BIBLIO_PROPS to generate biblio_str for reference output
SNIPPET_PROPERTIES = 4 # retrieves only PS_SENTENCE_PROPS to generate biblio_str for reference output
ONLY_REL_PROPERTIES = 5
ALL_PROPERTIES = 10

NO_RNEF_REL_PROPS={'RelationNumberOfReferences','Name','URN'}


class APISession(PSNetworx):
    '''
    manages retreival from database
    '''
    pass
    ResultRef = None
    ResultPos = 0
    ResultSize = None
    PageSize = 1000
    GOQLquery = str()
    add2self = True # if False will not add new graph to self.Graph
    __IsOn1st_page = True # to print header in dump file
    relProps = ['URN'] # URN is required for proper adding relations to ResnetGraph
   # __relPropCache__ = relProps # for temporary caching relProps if props2load() was used
    entProps = ['Name'] # Name is highly recommended to make sense of the data

    reference_cache_size = 1000000 # max number of reference allowed in self.Graph. Clears self.Graph if exceeded
    resnet_size = 1000 # number of <node><control> sections in RNEF dump
    max_rnef_size = 100000000 # max size of RNEF XML dump file. If dump file exceeds max_file_size new file is opened with index++
    cache_dir = 'ElsevierAPI/ResnetAPI/__pscache__/'
    cache_dict = PSObject() # {cache_name:[GOQLquery]}, [GOQLquery] - list of GOQL queries used to make the network in cache
    # all cached networks are written into cache_dir as file named cache_name.rnef

    DumpFiles = []
    csv_delimeter = '\t'
    print_rel21row = False
    clear_graph_cache = False
    __getLinks = True
    APIconfig = dict()
    data_dir = ''
    debug = True

    def __init__(self,APIconfig:dict,what2retrieve=NO_REL_PROPERTIES):
        super().__init__(APIconfig['ResnetURL'],APIconfig['PSuserName'],APIconfig['PSpassword'])
        self.APIconfig = APIconfig
        self.DumpFiles = []
        self.__retrieve(what2retrieve)

######################################  CONFIGURATION  ######################################
    def __retrieve(self,what2retrieve=CURRENT_SPECS):
        if what2retrieve == NO_REL_PROPERTIES:
            self.relProps = ['URN']
        elif what2retrieve == DATABASE_REFCOUNT_ONLY:
            self.relProps = ['URN',REFCOUNT]
        elif what2retrieve == REFERENCE_IDENTIFIERS:
            self.relProps = ['URN']+PS_REFIID_TYPES
        elif what2retrieve == BIBLIO_PROPERTIES:
            self.relProps = ['URN']+PS_REFIID_TYPES+list(PS_BIBLIO_PROPS)
        elif what2retrieve == SNIPPET_PROPERTIES:
            self.relProps = ['URN']+PS_REFIID_TYPES+list(PS_BIBLIO_PROPS)+PS_SENTENCE_PROPS
        elif what2retrieve == ONLY_REL_PROPERTIES:
            self.relProps = ['URN']+list(RELATION_PROPS)
        elif what2retrieve == ALL_PROPERTIES:
            self.relProps = ['URN']+ALL_PS_PROPS
        else:
            return #keeps current self.relProps
            #use 


    def _clone(self, to_retreive=CURRENT_SPECS):
        """
        used by self.process_oql to clone session 
        """
        if to_retreive == CURRENT_SPECS:
            new_session = APISession(self.APIconfig)
            new_session.entProps = list(self.entProps)
            new_session.relProps = list(self.relProps)
        else:
            new_session = APISession(self.APIconfig, to_retreive)
            new_session.__retrieve(to_retreive)
            new_session.entProps = list(self.entProps)
        return new_session


    def __set_get_links(self,oql_query=''):
        my_query = oql_query if oql_query else self.GOQLquery
        return my_query[7:15] == 'Relation'


    def __replace_goql(self, oql_query:str):
        self.GOQLquery = oql_query
        self.__getLinks = self.__set_get_links()
        self.ResultRef = None


    def set_dir(self,dir:str):
        self.data_dir = dir
        if dir and dir[-1] != '/':
            self.data_dir += '/'

    def _dir(self):
        return self.data_dir

  #  def set_temp_relProp(self,new_props:list):
  #      self.__relPropCache__ = list(self.relProps)
  #      self.relProps = list(new_props)
        # use restore_relProps() at the end of the function

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

#########################  PAGE BY PAGE PROCESSING, RETRIEVAL   ##################################
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
            return self.__get_next_page()                                                  

        if type(zeep_data) != type(None):
            if self.__getLinks:
                obj_ids = list(set([x['EntityId'] for x in zeep_data.Links.Link]))
                zeep_objects = self.get_object_properties(obj_ids, self.entProps)
                return self._load_graph(zeep_data, zeep_objects,self.add2self)
            else:
                return self._load_graph(None, zeep_data,self.add2self)
        else:
            return ResnetGraph()

        
    def __get_next_page(self, no_mess=True):
        if self.ResultPos < self.ResultSize:
            if not no_mess:
                print('fetched %d results out of %d' % (self.ResultPos, self.ResultSize))

            obj_props = self.relProps if self.__getLinks else self.entProps
            zeep_data, self.ResultSize, self.ResultPos = self.get_session_page(self.ResultRef, self.ResultPos, self.PageSize,
                                                                       self.ResultSize, obj_props ,getLinks=self.__getLinks)                                                                      
            if type(zeep_data) != type(None):
                if self.__getLinks and len(zeep_data.Links.Link) > 0:
                    obj_ids = list(set([x['EntityId'] for x in zeep_data.Links.Link]))
                    zeep_objects = self.get_object_properties(obj_ids, self.entProps)
                    return self._load_graph(zeep_data, zeep_objects,self.add2self)
                else:
                    return self._load_graph(None, zeep_data,self.add2self)
            else:
                return ResnetGraph()
        else: return ResnetGraph()


    def process_oql(self, oql_query, request_name='', flush_dump=False, debug=False, no_mess=True, 
                    iteration_limit=1) -> ResnetGraph:

        if flush_dump and self.ResultPos == 0:
            self.flush_dump_files(no_mess)
       
        self.__replace_goql(oql_query)
        global_start = time.time()
        entire_graph = ResnetGraph()
        reference_counter = 0
        
        start_time = time.time()
        return_type = 'relations' if self.__getLinks else 'entities'
        page_graph = self.__init_session()
        number_of_iterations = math.ceil(self.ResultSize / self.PageSize)

        my_request_name = request_name if request_name else self.GOQLquery[:100]
        if number_of_iterations > 1:
            print('\n\"%s\"\nrequest found %d %s. Begin retrieval in %d iterations' % 
                (my_request_name,self.ResultSize,return_type, number_of_iterations))
            if debug:
                print("Processing GOQL query:\n\"%s\"\n" % (self.GOQLquery))
            if no_mess:
                print('Progress report is suppressed. Retrieval may take long time - be patient!')

        while page_graph:
            page_ref_count = page_graph.size(weight='weight')
            reference_counter += page_ref_count
            exec_time = self.execution_time(start_time)
            iteration = math.ceil(self.ResultPos / self.PageSize) 
            if not no_mess:           
                edge_count = page_graph.number_of_edges()
                node_count = page_graph.number_of_nodes()
                print("Iteration %d in %d retrieved %d relations for %d nodes supported by %d references in %s seconds" %
                     (iteration, number_of_iterations,edge_count,node_count,page_ref_count,exec_time))
                start_time = time.time()
            
            self.Graph.add_graph(page_graph)
            
            if len(self.DumpFiles) > 0:
                self.to_csv(self.DumpFiles[0], in_graph=page_graph, access_mode='a',debug=True)
                self.__IsOn1st_page = False
                if number_of_iterations > 1:
                    print("With %d in %d iterations, %d %s in %d with %d references saved into \"%s\" file. Retrieval time: %s" % 
                        (iteration,number_of_iterations,self.ResultPos,return_type,self.ResultSize,reference_counter,self.DumpFiles[0],exec_time))
            else:
                if number_of_iterations > 1:
                    # progress report:
                    print("With %d in %d iterations, %d %s in %d with %d references retrieved in: %s" % 
                        (iteration,number_of_iterations,self.ResultPos,return_type,self.ResultSize,reference_counter,exec_time))

            entire_graph = page_graph.compose(entire_graph)

            if debug and number_of_iterations >= iteration_limit: break
            if self.clear_graph_cache: self.clear()
            page_graph = self.__get_next_page(no_mess)

        if debug: print("GOQL query:\n \"%s\"\n was executed in %s in %d iterations" % 
                 (self.GOQLquery, self.execution_time(global_start), number_of_iterations))
            
        if number_of_iterations == 1: 
            print('"%s"\nwas retrieved in %s by %d iteration' % 
                (my_request_name, self.execution_time(global_start), number_of_iterations))
        else:
            print('"%s"\nwas retrieved in %s by %d iterations' % 
                (my_request_name, self.execution_time(global_start), number_of_iterations))

        self.ResultRef = ''
        self.ResultPos = 0
        self.ResultSize = 0
        self.__IsOn1st_page = True
        self.clone = False
        #self.what2retrieve = CURRENT_SPECS
        return entire_graph


    def __annotated_graph_by_ids__(self,id_only_graph:ResnetGraph,use_cache=True):
        '''
        loads
        -----
        new relations to self.Graph by ids in id_only_graph
        '''
        relation_ids2return = id_only_graph._relations_ids()
        if relation_ids2return:
            if use_cache:
                relation_ids2retreive = relation_ids2return.difference(self.id2relation.keys())
                nodes_ids2retreive = set(id_only_graph.nodes()).difference(set(self.Graph.nodes()))
            else:
                relation_ids2retreive = list(relation_ids2return)
                nodes_ids2retreive = list(id_only_graph.nodes())

            print('Retreiving additional graph by id with %d nodes and %d relations' % (len(nodes_ids2retreive),len(relation_ids2retreive)))
            if nodes_ids2retreive:
                entity_query = 'SELECT Entity WHERE id = ({ids})'
                self.__iterate_id__(entity_query,list(nodes_ids2retreive))

            if relation_ids2retreive:
                rel_query = 'SELECT Relation WHERE id = ({ids})'
                self.__iterate_id__(rel_query,list(relation_ids2retreive))
    
            # now all relation_ids2return should be in self.id2relation
            rels4subgraph = [rel for i,rel in self.id2relation.items() if i in relation_ids2return]
            assert(len(rels4subgraph) == len(relation_ids2return))
            return_subgraph = self.Graph.subgraph_by_rels(rels4subgraph)
            return return_subgraph
        else:
            # case when id_only_graph has only nodes
            nodes_ids2return = list(id_only_graph.nodes())
            if use_cache:
                nodes_ids2retreive = set(nodes_ids2return).difference(set(self.Graph.nodes()))
            else:
                nodes_ids2retreive = set(nodes_ids2return)

            if nodes_ids2retreive:
                entity_query = 'SELECT Entity WHERE id = ({ids})'
                self.__iterate_id__(entity_query,list(nodes_ids2retreive))

            nodes_graph = ResnetGraph()
            nodes_graph.add_nodes({i:node for i,node in self.Graph.nodes(data=True) if i in nodes_ids2return})
            assert(len(nodes_graph) == len(nodes_ids2return))
            return nodes_graph


    def __iterate_id__(self,oql_query:str,id_list:list):
        entire_graph = ResnetGraph()
        step = 1000
        for i in range(0,len(id_list), step):
            ids = id_list[i:i+step]
            oql_query_with_ids = oql_query.format(ids=','.join(map(str,ids)))
            iter_graph = self.process_oql(oql_query_with_ids)
            entire_graph.add_graph(iter_graph)
        return entire_graph 


    def __iterate_oql4id__(self,oql_query:str,id_set:set):
        # oql_query MUST contain string placeholder called {ids} 
        my_ent_props = list(self.entProps)
        my_rel_props = list(self.relProps)

        self.relProps.clear()
        self.entProps.clear()
        self.add2self = False
        entire_id_graph = ResnetGraph()
        entire_id_graph = self.__iterate_id__(oql_query, list(id_set))

        self.entProps = my_ent_props
        self.relProps = my_rel_props
        self.add2self = True
        return entire_id_graph


    def __iterate_oql4id_s__(self,oql_query:str,prop_set:set):
        '''
            # oql_query MUST contain string placeholder called {props} 
        '''
        my_ent_props = list(self.entProps)
        my_rel_props = list(self.relProps)

        self.relProps.clear()
        self.entProps.clear()
        self.add2self = False

        entire_id_graph = ResnetGraph()
        prop_list = list(prop_set)
        step = 950 # to mitigate oql_query length limitation 
        for i in range(0,len(prop_list), step):
            props = prop_list[i:i+step]
            oql_query_with_props = oql_query.format(props=OQL.join_with_quotes(props))
            iter_graph = self.process_oql(oql_query_with_props)
            entire_id_graph.add_graph(iter_graph)

        self.entProps = my_ent_props
        self.relProps = my_rel_props
        self.add2self = True
        return entire_id_graph


    def __iterate_oql4id_2__(self, oql_query:str, id_set1:set, id_set2:set):
        '''
        # oql_query MUST contain 2 string placeholders called {ids1} and {ids2}
        '''
        my_ent_props = list(self.entProps)
        my_rel_props = list(self.relProps)

        self.relProps.clear()
        self.entProps.clear()
        self.add2self = False

        entire_id_graph = ResnetGraph()
        id_list1 = list(id_set1)
        id_list2 = list(id_set2)
        step = 1000
        number_of_iterations = math.ceil(len(id_set1)/step) * math.ceil(len(id_set2)/step)
        print('\nConnecting %d with %d entities' % (len(id_set1), len(id_set2)))
        if number_of_iterations > 2:
            print('Query will be executed in %d iterations' % number_of_iterations)

        iteration_counter = 1
        start  = time.time()
        for i1 in range(0,len(id_list1), step):
            ids1 = id_list1[i1:i1+step]
            for i2 in range(0,len(id_list2), step):
                ids2 = id_list2[i2:i2+step]
                oql_query_with_ids = oql_query.format(ids1=','.join(map(str,ids1)),ids2=','.join(map(str,ids2)))
                iter_graph = self.process_oql(oql_query_with_ids)
                entire_id_graph.add_graph(iter_graph)
                if number_of_iterations > 2:
                    print('Iteration %d out of %d performed in %s' % 
                        (iteration_counter,number_of_iterations,self.execution_time(start)))
                iteration_counter +=1

        self.entProps = my_ent_props
        self.relProps = my_rel_props
        self.add2self = True
        return entire_id_graph


    def __iterate_oql4id_2s__(self, oql_query:str, id_set1:set, id_set2:set):
        '''
        # oql_query MUST contain 2 string placeholders called {props1} and {props2}
        '''
        my_ent_props = list(self.entProps)
        my_rel_props = list(self.relProps)

        self.relProps.clear()
        self.entProps.clear()
        self.add2self = False

        entire_id_graph = ResnetGraph()
        id_list1 = list(id_set1)
        id_list2 = list(id_set2)
        step = 950
        number_of_iterations = math.ceil(len(id_set1)/step) * math.ceil(len(id_set2)/step)
        print('\nConnecting %d with %d entities' % (len(id_set1), len(id_set2)))
        if number_of_iterations > 2:
            print('Query will be executed in %d iterations' % number_of_iterations)

        iteration_counter = 1
        start  = time.time()
        for i1 in range(0,len(id_list1), step):
            ids1 = id_list1[i1:i1+step]
            for i2 in range(0,len(id_list2), step):
                ids2 = id_list2[i2:i2+step]
                oql_query_with_ids = oql_query.format(props1=','.join(map(str,ids1)),props2=','.join(map(str,ids2)))
                iter_graph = self.process_oql(oql_query_with_ids)
                entire_id_graph.add_graph(iter_graph)
                if number_of_iterations > 2:
                    print('Iteration %d out of %d performed in %s' % 
                        (iteration_counter,number_of_iterations,self.execution_time(start)))
                iteration_counter +=1

        self.entProps = my_ent_props
        self.relProps = my_rel_props
        self.add2self = True
        return entire_id_graph


    def iterate_oql(self, oql_query:str,id_set:set,use_cache=True,request_name='',iterate_ids=True):
        """
        # oql_query MUST contain string placeholder called {ids} if iterate_ids==True\n
        # oql_query MUST contain string placeholder called {props} if iterate_ids==False
        """
        print('Processing "%s" request\n' % request_name)
        if iterate_ids:
            id_only_graph = self.__iterate_oql4id__(oql_query,id_set)
        else:
            id_only_graph = self.__iterate_oql4id_s__(oql_query,id_set)

        return self.__annotated_graph_by_ids__(id_only_graph,use_cache)


    def iterate_oql2(self,oql_query:str,id_set1:set,id_set2:set,use_cache=True,request_name='',iterate_ids=True):
        """
        # oql_query MUST contain 2 string placeholders called {ids1} and {ids2}
        """
        print('Processing "%s" request\n' % request_name)
        if iterate_ids:
            id_only_graph = self.__iterate_oql4id_2__(oql_query,id_set1,id_set2)
        else:
            id_only_graph = self.__iterate_oql4id_2s__(oql_query,id_set1,id_set2)
            
        print('Will retrieve annotated graph with %d entities and %d relations' %(id_only_graph.number_of_nodes(),id_only_graph.number_of_edges()))
        return self.__annotated_graph_by_ids__(id_only_graph,use_cache)


    def clear(self):
        self.Graph.clear_resnetgraph()
        self.id2relation.clear()
        self.ID2Children.clear()


    def flush_dump_files(self, no_mess=True):
        for f in self.DumpFiles:
            open(f, 'w').close()
            if not no_mess:
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
                page_graph = self.__get_next_page()
                entire_graph = nx.compose(page_graph, entire_graph)
                iteration = math.ceil(self.ResultPos / self.PageSize)
                total_iterations = math.ceil(self.ResultSize/self.PageSize)
                print ('%d out of %d iterations fetching %d objects from %s search results performed in %s' % 
                (iteration, total_iterations, self.ResultSize, name[0], self.execution_time(start)))

        return entire_graph


    def _get_ppi_graph(self,fout):
        # Build PPI network between proteins
        graph_proteins = self.Graph.get_node_ids(['Protein'])
        print('Retrieving PPI network for %d proteins' % (len(graph_proteins)))
        start_time = time.time()
        ppi_graph = self.get_ppi(graph_proteins, self.relProps, self.entProps)
        self.to_csv(self.data_dir+fout, in_graph=ppi_graph)
        print("PPI network with %d relations was retrieved in %s ---" % (
            ppi_graph.size(), self.execution_time(start_time)))

        return ppi_graph


    def _get_network_graph(self, for_interactor_ids:set, connect_by_rel_types:list=None):
        #self = self.__clone() if self.clone else self
        network2return = self.get_network(for_interactor_ids,connect_by_rel_types,
                                            self.relProps,self.entProps)
        return network2return


    def connect_nodes(self,node_ids1:set,node_ids2:set,
                        by_relation_type=[],with_effect=[],in_direction='',
                        use_relation_cache=True):
        """
        Input
        -----
        in_direction must be '>' or '<'

        Returns
        -----
        ResnetGraph containing relations between node_ids1 and node_ids2
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

        req_name = f'Connecting {str(len(node_ids1))} and {str(len(node_ids2))} entities by {reltype_str} relations and {effect_str} effects'
        return self.iterate_oql2(f'{oql_query}',node_ids1,node_ids2,use_relation_cache,request_name=req_name)


    def get_group_members(self, group_names:list):
        oql_query = 'SELECT Entity WHERE MemberOf (SELECT Group WHERE Name = ({group}))'
        groups = OQL.join_with_quotes(group_names)
        req_name = 'Find members of groups: ' + ','.join(group_names)
        graph2return = self.process_oql(oql_query.format(group=groups), request_name=req_name)
        if len(graph2return) == 0:
            print('%s groups are empty or do not exist in databse' % str(group_names))
        else:
            print('loaded %d members from %s' % (graph2return.number_of_nodes(),str(group_names)))
        return graph2return


    def map_props2objs(self,using_values:list,in_properties:list,map2types=[],case_insensitive=False):
        """
        Returns
        -------
        propval2objs = {prop_value:[PSObject]}\n
        objid2propval = {node_id:[prop_values]},\n  where prop_value is from 'using_values'
        """
        propval2objs,objid2propval = self.Graph.get_props2obj_dic(using_values, in_properties,case_insensitive)
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

        p2obj,id2p = self.Graph.get_props2obj_dic(need_db_mapping,in_properties,case_insensitive)
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
    def child_graph(self, propValues:list, search_by_properties=[],include_parents=True):
        if not search_by_properties: search_by_properties = ['Name','Alias']
        oql_query = OQL.get_childs(propValues,search_by_properties,include_parents=include_parents)
        request_name = 'Find ontology children for {}'.format(','.join(propValues))
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
        parent_id = add2parent.id()
        ontology_graph.add1node(add2parent)
        ontology_graph.add_nodes(members)
        for m in members:
            child_id = m['Id'][0]
            rel = PSRelation({'ObjTypeName':['MemberOf'],'Relationship':['is-a'],'Ontology':['Pathway Studio Ontology']})
            rel.Nodes['Regulators'] = [(child_id,0,0)]
            rel.Nodes['Targets'] = [(parent_id,1,0)]
            rel.append_property('Id',child_id)
            rel.append_property('URN', str(child_id)+'0') #fake URN
            ontology_graph.add_edge(child_id,parent_id,relation=rel)

        self.add_rel_props(['Relationship','Ontology'])
        return ontology_graph


    def load_children(self,for_parent_ids=[],add_new_nodes=False):
        """
        slower than self.child_graph() but annotates each parent with CHILDS property
        """
        annotate_parent_ids = set(for_parent_ids) if for_parent_ids else set(self.Graph)
        self.get_children(annotate_parent_ids)

        if add_new_nodes:
            [self.Graph.property2node(parent_id,CHILDS,self.ID2Children[parent_id]) for parent_id in annotate_parent_ids]
        else:
            for parent_id in annotate_parent_ids:
                all_child_ids = set(self.ID2Children[parent_id])
                graph_child_ids = all_child_ids.intersection(set(self.Graph))
                self.Graph.property2node(parent_id,CHILDS,list(graph_child_ids))


    def add_parents(self,for_child_ids=[],depth=1):
        """
        Adds
        ----
        ontology parents to_graph\n
        annotates each parent with CHILDS property

        Returns
        -------
        [PSObject] - list of parents with CHILDS property
        """
        if not for_child_ids:
            for_child_ids = list(self.Graph)

        get_parent_query = 'SELECT Entity WHERE InOntology (SELECT Annotation WHERE Ontology=\'Pathway Studio Ontology\' AND Relationship=\'is-a\') inRange {steps} over (SELECT OntologicalNode WHERE id = ({ids}))'
        oql_query = get_parent_query.format(steps=str(depth),ids=','.join(map(str,for_child_ids)))
        request_name = 'Find parents of {count} nodes with depth {d}'.format(count = str(len(for_child_ids)), d=str(depth))
        parent_graph = self.process_oql(oql_query,request_name)

      #  if to_graph:
       #     graph_registry[index].add_graph(parent_graph)
        
        self.load_children(for_parent_ids=list(parent_graph))
        return self.Graph._get_nodes(list(parent_graph))



    def load_ontology(self, parent_ontology_groups:list):
        """
        Input
        -----
        [Parent Name]
        loads self.child2parent = {child_name:[ontology_parent_names]}
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

###################################   WRITE, DUMP   ##############################################
    def to_csv(self, file_out, in_graph: ResnetGraph=None, access_mode='w', debug=False):
        if not isinstance(in_graph,ResnetGraph): in_graph = self.Graph
        in_graph.print_references(file_out, self.relProps, self.entProps, access_mode, printHeader=self.__IsOn1st_page,
                                col_sep=self.csv_delimeter,debug=debug,single_rel_row=self.print_rel21row)


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
    def __make_file_name(folder_or_object_name: str):
        new_str = list(folder_or_object_name)
        for i in range(0, len(new_str)):
            if new_str[i] in {'>', '<', '|','/'}:
                new_str[i] = '-'
            if new_str[i] in {':'}:
                new_str[i] = '_'
        return "".join(new_str)


    def __find_folder(self, folder_name:str):
        folder_name_ = self.__make_file_name(folder_name)
        longest_path = str()
        for dirpath, dirnames, filenames in os.walk(self.data_dir):
            for dirname in dirnames:
                if dirname == folder_name_:
                    if len(dirpath+dirname)+1 > len(longest_path):
                        longest_path = dirpath+'/'+dirname +'/'

        return longest_path


    def dump_path(self, of_folder:str, in2parent_folder=''):
        # finding directory
        my_dir = self.__find_folder(of_folder)
        if not my_dir:
            if in2parent_folder:
                parent_dir = self.__find_folder(in2parent_folder)
                if parent_dir:
                    my_dir = parent_dir+self.__make_file_name(of_folder)+'/'
                    os.mkdir(my_dir) 
                else:
                    my_dir = self.data_dir+self.__make_file_name(in2parent_folder)+'/'+self.__make_file_name(of_folder)+'/'
                    os.mkdir(my_dir)
            else:
                my_dir = self.data_dir+self.__make_file_name(of_folder)+'/'
                os.mkdir(my_dir)

        return my_dir


    def __count_dumpfiles(self, in_folder:str,in2parent_folder='', with_extension='rnef'):
        folder_path = self.dump_path(in_folder,in2parent_folder)
        listing = glob.glob(folder_path+'*.' + with_extension)
        if listing:
            return folder_path,len(listing)
        else:
            return folder_path,1


    def __dump_base_name(self,of_folder:str):
        return 'content of '+ APISession.__make_file_name(of_folder)+'_'


    def __make_dumpfile_name(self,of_folder:str,in2parent_folder='',new=False):
        folder_path,file_count = self.__count_dumpfiles(of_folder,in2parent_folder)
        if new:
            file_count += 1
        return folder_path + self.__dump_base_name(of_folder)+str(file_count)+'.rnef'


    def str2rnefdump(self,rnef_xml:str,to_folder='',in2parent_folder='',can_close=True):
        '''
        Dumps
        -----
        "rnef_xml" into 'to_folder' inside 'in2parent_folder' located in "self.data_dir"
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
            with open(write2,'w',encoding='utf-8') as f:
                f.write('<batch>\n'+rnef_xml)


    def to_rnef(self,graph=ResnetGraph(),add_rel_props:dict={},add_pathway_props:dict={}):
        '''
        used for dumping pathway, group and results objects by FolderContent
        Returns
        -------
        graph RNEF XML with single <resnet> section and session properties for nodes and edges
        '''
        my_graph = graph if graph else self.Graph
        rel_props = [p for p in self.relProps if p not in NO_RNEF_REL_PROPS]
        return my_graph.rnef(self.entProps,rel_props,add_rel_props,add_pathway_props)


    def graph2rnefdump(self,graph=ResnetGraph(),to_folder='',in_parent_folder='',can_close=True):
        '''
        Dumps
        -----
        large graph objects into several RNEF XML files
        '''
        my_graph = graph if graph else self.Graph
        rel_props = [p for p in self.relProps if p not in NO_RNEF_REL_PROPS]

        resnet_sections = ResnetGraph()
        for regulatorID, targetID, e in my_graph.edges(data='relation'):
            resnet_sections.copy_rel(e,my_graph)
            if resnet_sections.number_of_edges() == self.resnet_size:
                rnef_str = resnet_sections.rnef(ent_props=self.entProps,rel_props=rel_props)
                rnef_str = self.pretty_xml(rnef_str,no_declaration=True)
                self.str2rnefdump(rnef_str,to_folder,in_parent_folder,can_close)
                resnet_sections.clear_resnetgraph()

        rnef_str = resnet_sections.rnef(ent_props=self.entProps,rel_props=self.relProps)
        rnef_str = self.pretty_xml(rnef_str,no_declaration=True)
        self.str2rnefdump(rnef_str,to_folder,in_parent_folder,can_close)


    def close_rnef_dump(self,to_folder='',in2parent_folder=''):
        last_dump_file = self.__make_dumpfile_name(to_folder,in2parent_folder)
        f = open(last_dump_file,'a',encoding='utf-8')
        f.write('</batch>')
        f.close()


    def download_oql(self,oql_query,request_name:str,resume_page=0):
        '''
        Use for oql_query producing large results

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
        page_graph = self.__init_session(resume_pos)
        number_of_iterations = math.ceil(self.ResultSize / self.PageSize)
        iterations_str = 'iteration' 
        if number_of_iterations > 1:
            iterations_str += 's'
            print('\n\nrequest found %d %s. Begin download in %d iterations' % 
                (self.ResultSize,return_type, number_of_iterations))
            
        while page_graph:
            page_ref_count = page_graph.weight()
            reference_counter += page_ref_count
            exec_time = self.execution_time(start_time)
            iteration = math.ceil(self.ResultPos / self.PageSize)
            print("With %d in %d %s, %d %s in %d with %d references retrieved in: %s" % 
            (iteration,number_of_iterations,iterations_str,self.ResultPos,return_type,self.ResultSize,reference_counter,exec_time))
                         
            self.graph2rnefdump(in_graph=page_graph, to_folder=request_name)
            self.clear()
            page_graph = self.__get_next_page()

        self.close_rnef_dump(to_folder=request_name)
        self.clear()
        self.ResultRef = ''
        self.ResultPos = 0
        self.ResultSize = 0
        self.__IsOn1st_page = True
        self.clone = False
        return


    def load_cache(self,cache_name:str,oql_queries:list=[],ent_props=['Name'],rel_props=['URN'])->'ResnetGraph':
        cache_file = self.cache_dir+cache_name+'.rnef'
        try:
            return ResnetGraph.fromRNEF(cache_file)
        except FileNotFoundError:
            database_graph = ResnetGraph.fromRNEFdir(self.cache_dir+cache_name)
            if not database_graph:
                print('Cache %s was not found\nBegin downloading network from database' % cache_name)
                my_session = self._clone(REFERENCE_IDENTIFIERS) #need identifiers to make graph simple
                my_session.add_ent_props(ent_props)
                my_session.add_rel_props(rel_props)
                my_session.set_dir(self.cache_dir)
                for oql_query in oql_queries:
                    self.download_oql(oql_query,cache_name)
                my_session.clear()
                database_graph = ResnetGraph.fromRNEFdir(self.cache_dir+cache_name)
                #  network4oql = my_session.process_oql(oql_query, f'Fetching {cache_name}')
                #  downloaded_network = downloaded_network.compose(network4oql)
            
            downloaded_network = database_graph.make_simple()
            network_psobj = PSPathway({'Name':[cache_name],'OQLs':oql_queries},downloaded_network)
            network_xml_str = network_psobj.to_xml(ent_props,rel_props)
            with open(cache_file,'w',encoding='utf-8') as f: 
                f.write(network_xml_str)
            
            print('%s with %d edges and %d nodes was downloaded from database in %s' 
                % (cache_name,downloaded_network.number_of_edges(),downloaded_network.number_of_nodes(),self.execution_time(start)))
            
            return downloaded_network

    
    #####################   EXPERIMENT EXPERIMENT EXPERIMENT EXPERIMENT ###########################
    def map_experiment(self, exp:Experiment):
        identifier_name = exp.identifier_name()
        identifier_names = identifier_name.split(',')
        identifiers = exp.list_identifiers()
        map2objtypes = exp['ObjTypeName']
        identifier2objs, objid2prop = self.map_props2objs(identifiers,identifier_names,map2objtypes)

        unmapped_identifiers = set(identifiers).difference(set(identifier2objs.keys()))
        print('Following experiment identifiers were not mapped:\n%s' % '\n'.join(unmapped_identifiers))

        identifiers_columns = exp.identifiers.columns.to_list()
        mapped_entities_pd = df(columns = [identifier_name,'URN'])
        row = 0
        for i, psobj_list in identifier2objs.items():
            for psobj in psobj_list:
                mapped_entities_pd.at[row,identifier_name] = i
                mapped_entities_pd.at[row,'URN'] = psobj['URN'][0]
                row += 1
        
        samples_pd = exp.samples2pd()
        new_data = mapped_entities_pd.merge(samples_pd,'outer', on=identifier_name)
        new_data = new_data.drop('URN_y', axis=1)
        new_data.rename(columns={'URN_x':'URN'},inplace=True)
        mapped_exp = Experiment(exp)
        mapped_exp.identifiers = new_data[identifiers_columns]

        column_counter = len(identifiers_columns)
        for sample in exp.get_samples():
            mapped_sample = Sample(sample)
            mapped_sample.data['value'] = samples_pd.iloc[:,column_counter]
            if sample.has_pvalue():
                column_counter += 1
                mapped_sample.data['pvalue'] = samples_pd.iloc[:,column_counter]
            mapped_exp._add_sample(mapped_sample)

        print('%d out of %d entities were mapped in %s'%
                (len(identifier2objs), len(exp.identifiers),exp.name()))

        return mapped_exp


    def common_neighbors(self,with_entity_types:list,of_ids1:list,reltypes12:list,effect12:list,dir12:str,
                                and_ids3:list,reltypes23:list,effect23:list,dir23:str):
        """
            Input
            -----
            dir12,dir23 must be '<','>' or empty string

            Returns
            -------
            ResnetGraph with 1-->2-->3 relations, where nodes in position 1 are from "of_ids1" and nodes in position 3 are from "and_ids3"\n
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

        r_n = f'Find entities common between {len(of_ids1)} and {len(and_ids3)} entities' 
        neighbors_entity_graph = self.iterate_oql2(f'{find_neighbors_oql}',of_ids1,and_ids3,request_name=r_n)
        neighbors_ids = list(neighbors_entity_graph.nodes())
        print('Found %d common neighbors' % len(neighbors_ids))

        if dir12 == '>':
            sel_12_graph_oql = sel_rels12 + r' AND downstream NeighborOf (SELECT Entity WHERE id = ({ids1})) AND upstream NeighborOf (SELECT Entity WHERE id = ({ids2}))'
        elif dir12 == '<':
            sel_12_graph_oql = sel_rels12 + r' AND upstream NeighborOf (SELECT Entity WHERE id = ({ids1})) AND downstream NeighborOf (SELECT Entity WHERE id = ({ids2}))'
        else:
            sel_12_graph_oql = sel_rels12 + r' AND NeighborOf (SELECT Entity WHERE id = ({ids1})) AND NeighborOf (SELECT Entity WHERE id = ({ids2}))'
        rn12 = f'Find 1-2 relations between {len(of_ids1)} entities in position 1 and {len(neighbors_ids)} common neighbors'
        #graph12 = self.iterate_oql2(f'{sel_12_graph_oql}',of_ids1,neighbors_ids,request_name=rn12)

        if dir23 == '>':
            sel_23_graph_oql = sel_rels23 + r' AND upstream NeighborOf (SELECT Entity WHERE id = ({ids1})) AND downstream NeighborOf (SELECT Entity WHERE id = ({ids2}))'
        elif dir12 == '<':
            sel_23_graph_oql = sel_rels23 + r' AND downstream NeighborOf (SELECT Entity WHERE id = ({ids1})) AND upstream NeighborOf (SELECT Entity WHERE id = ({ids2}))'
        else:
            sel_23_graph_oql = sel_rels23 + r' AND NeighborOf (SELECT Entity WHERE id = ({ids1})) AND NeighborOf (SELECT Entity WHERE id = ({ids2}))'
        rn23 = f'Find 2-3 relations between {len(and_ids3)} entities in position 3 and {len(neighbors_ids)} common neighbors'
       # graph23 = self.iterate_oql2(f'{sel_23_graph_oql}',and_ids3,neighbors_ids,request_name=rn23)

        with ThreadPoolExecutor(max_workers=4, thread_name_prefix='FindCommonNeighbors') as executor:
            future12 = executor.submit(self.iterate_oql2,f'{sel_12_graph_oql}',of_ids1,neighbors_ids, True,rn12) 
            future23 = executor.submit(self.iterate_oql2,f'{sel_23_graph_oql}',and_ids3,neighbors_ids, True,rn23)
            
            graph12 = future12.result()
            graph23 = future23.result()

        graph12.add_graph(graph23)
        return graph12


    def common_neighbors_with_effect(self,with_entity_types:list,of_ids1:list,reltypes12:list,dir12:str,
                                and_ids3:list,reltypes23:list,dir23:str,id_type='id'):
        # written to circumvent bug in Patwhay Studio
        """
            Input
            -----
            dir12,dir23 must be '<','>' or empty string

            Returns
            -------
            ResnetGraph with 1-->2-->3 relations, where nodes in position 1 are from "of_ids1" and nodes in position 3 are from "and_ids3"\n
            if nodes in positon 1 and positon 3 are neighbors returns graph to complete 3-vertex cliques
        """
        sel_rels = 'SELECT Relation WHERE objectType = ({}) AND NOT Effect = unknown'

        sel_rels12 = sel_rels.format(','.join(reltypes12))
        sel_rels23 = sel_rels.format(','.join(reltypes23))

        find_neighbors_oql = r'SELECT Entity objectType = ('+','.join(with_entity_types)+')'

        find_neighbors_oql = find_neighbors_oql+' AND Connected by ('+sel_rels12+') to (SELECT Entity WHERE '+id_type+' = ({ids1}))'
        find_neighbors_oql += ' AND Connected by ('+sel_rels23+') to (SELECT Entity WHERE '+id_type+' = ({ids2}))'

        r_n = f'Find entities common between {len(of_ids1)} and {len(and_ids3)} entities' 
        neighbors_entity_graph = self.iterate_oql2(f'{find_neighbors_oql}',of_ids1,and_ids3,request_name=r_n)
        neighbors_ids = list(neighbors_entity_graph.nodes())
        print('Found %d common neighbors' % len(neighbors_ids))

        if dir12 == '>':
            sel_12_graph_oql = sel_rels12 + r' AND downstream NeighborOf (SELECT Entity WHERE '+id_type+' = ({ids1})) AND upstream NeighborOf (SELECT Entity WHERE '+id_type+' = ({ids2}))'
        elif dir12 == '<':
            sel_12_graph_oql = sel_rels12 + r' AND upstream NeighborOf (SELECT Entity WHERE '+id_type+' = ({ids1})) AND downstream NeighborOf (SELECT Entity WHERE '+id_type+' = ({ids2}))'
        else:
            sel_12_graph_oql = sel_rels12 + r' AND NeighborOf (SELECT Entity WHERE '+id_type+' = ({ids1})) AND NeighborOf (SELECT Entity WHERE '+id_type+' = ({ids2}))'
        rn12 = f'Find 1-2 relations between {len(of_ids1)} entities in position 1 and {len(neighbors_ids)} common neighbors'
        #graph12 = self.iterate_oql2(f'{sel_12_graph_oql}',of_ids1,neighbors_ids,request_name=rn12)

        if dir23 == '>':
            sel_23_graph_oql = sel_rels23 + r' AND upstream NeighborOf (SELECT Entity WHERE '+id_type+' = ({ids1})) AND downstream NeighborOf (SELECT Entity WHERE '+id_type+' = ({ids2}))'
        elif dir12 == '<':
            sel_23_graph_oql = sel_rels23 + r' AND downstream NeighborOf (SELECT Entity WHERE '+id_type+' = ({ids1})) AND upstream NeighborOf (SELECT Entity WHERE '+id_type+' = ({ids2}))'
        else:
            sel_23_graph_oql = sel_rels23 + r' AND NeighborOf (SELECT Entity WHERE '+id_type+' = ({ids1})) AND NeighborOf (SELECT Entity WHERE '+id_type+' = ({ids2}))'
        rn23 = f'Find 2-3 relations between {len(and_ids3)} entities in position 3 and {len(neighbors_ids)} common neighbors'

        iterate_ids = True if id_type == 'id' else False
        with ThreadPoolExecutor(max_workers=4, thread_name_prefix='FindCommonNeighbors') as executor:
            session12 = self._clone()
            future12 = executor.submit(session12.iterate_oql2,f'{sel_12_graph_oql}',of_ids1,neighbors_ids,True,rn12,iterate_ids) 
            session23 = self._clone()
            future23 = executor.submit(session23.iterate_oql2,f'{sel_23_graph_oql}',and_ids3,neighbors_ids,True,rn23,iterate_ids)
            
            graph12 = future12.result()
            session12.close_connection()
            graph23 = future23.result()
            session23.close_connection()

        graph123 = graph12.compose(graph23)
        if self.add2self:
            self.Graph = self.Graph.compose(graph123)
        return graph123


    def gv2gene(self,gv_ids:list):
        """
        Input
        -----
        list of GV ids

        Returns
        -------
        {gv_id:[gene_names]}
        """
        prot2gvs_graph = ResnetGraph()
        print ('Finding genes for %d genetic variants' % len(gv_ids))
        number_of_iterations = int(len(gv_ids)/1000)+1
        for i in range(0, len(gv_ids),1000):
            chunk = gv_ids[i: i+1000]
            oql_query = OQL.expand_entity(PropertyValues=chunk, SearchByProperties=['id'], 
                                expand_by_rel_types=['GeneticChange'],expand2neighbors=['Protein'])
            request_name = f'{str(int(i/1000)+1)} iteration out of {str(number_of_iterations)} to find genes linked to GVs'
            prot2gvs_graph.add_graph(self.process_oql(oql_query,request_name))

        # making gvid2genes for subsequent annotation
        gvid2genes = dict()
        for gv_id, protein_id, rel in prot2gvs_graph.edges.data('relation'):
            protein_node = prot2gvs_graph._get_node(protein_id)
            gene_name = protein_node['Name'][0]
            try:
                gvid2genes[gv_id].append(gene_name)
            except KeyError:
                gvid2genes[gv_id] = [gene_name]

        #[nx.set_node_attributes(disease2gvs, {}) for gvid, gene_names in gvid2genes.items()]
        return gvid2genes

    '''
    def unified_rel(self,merge2rel:PSRelation):
        unified_rel = merge2rel.copy()
        graph_between = ResnetGraph()
        for regulator_id, target_id in merge2rel.get_regulators_targets():
            graph_between.add_graph(self.connect_nodes([regulator_id],[target_id],in_direction='>'))

        positive_refs,negative_refs = graph_between._effect_counts__() 
        if len(positive_refs) > len(negative_refs):
            my_effect =  'activated'
        elif len(negative_refs) > len(positive_refs):
            my_effect = 'repressed'
        else:
            my_effect = 'unknown'


        def choose_type(rels:list):
            type_ranks = ['DirectRegulation','Binding','ProtModification','PromoterBinding','Expression','MolTransport','MolSynthesis','Regulation']
            for type in type_ranks:
                for r in rels:
                    if r.objtype() == type:
                        return type
            return NotImplemented

        rel_between = graph_between.get_relations()
        if merge2rel.objtype() != 'ChemicalReaction':
            my_type = choose_type(rel_between)
        else:
            my_type = 'ChemicalReaction'
        
        [unified_rel.merge_rel(r) for r in rel_between]
        unified_rel.set_property('ObjTypeName',my_type)
        unified_rel.set_property('Effect',my_effect)
        return unified_rel
        

    def repair_graph(self,graph=ResnetGraph()):
        my_graph = graph if graph else self.Graph
        for regulator_id,target_id,rel in my_graph.edges(data='relation'):
            if rel.get_reference_count() == 0:
                new_rel = self.unified_rel(rel)
                my_graph.remove_relation(rel)
                my_graph.add_rel(new_rel)
    '''

        
        




        



        



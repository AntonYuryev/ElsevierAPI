import time
import math
import networkx as nx
from datetime import timedelta
from ElsevierAPI.ResnetAPI.ZeepToNetworkx import PSNetworx
from ElsevierAPI.ResnetAPI.ResnetGraph import ResnetGraph

class APISession(PSNetworx):
    pass
    ResultRef = str()
    ResultPos = int()
    ResultSize = int()
    PageSize = 100
    GOQLquery = str()
    __IsOn1st_page = True
    relProps = ['Name', 'RelationNumberOfReferences']
    entProps = ['Name']
    DumpFiles = ['ResnetAPIsessionDump.tsv']
    csv_delimeter = '\t'

    def __init__(self,url, username, password):
        super().__init__(url, username, password)
        self.DumpFiles = ['ResnetAPIsessionDump.tsv']

    def __init_session(self):
        if self.__getLinks:
            zeep_data, (self.ResultRef, self.ResultSize, self.ResultPos) = self.init_session(self.GOQLquery,
                                                                                             self.PageSize,
                                                                                             self.relProps)
            if type(zeep_data) != type(None):
                obj_ids = list(set([x['EntityId'] for x in zeep_data.Links.Link]))
                zeep_objects = self.get_object_properties(obj_ids, self.entProps)
                return self._load_graph(zeep_data, zeep_objects)
            else:
                return None
        else:
            zeep_data, (self.ResultRef, self.ResultSize, self.ResultPos) = self.init_session(self.GOQLquery,
                                                                                             self.PageSize,
                                                                                             self.entProps,
                                                                                             getLinks=False)
            if type(zeep_data) != type(None):
                return self._load_graph(None, zeep_data)
            else:
                return None

    def __get_next_page(self, no_mess=True):
        if self.ResultPos < self.ResultSize:
            if not no_mess:
                print('fetched %d results out of %d' % (self.ResultPos, self.ResultSize))
            if self.__getLinks:
                zeep_data, self.ResultPos = self.get_next_session_page(self.ResultRef, self.ResultPos, self.PageSize,
                                                                       self.ResultSize, self.relProps)
                if type(zeep_data) != type(None):
                    obj_ids = list(set([x['EntityId'] for x in zeep_data.Links.Link]))
                    zeep_objects = self.get_object_properties(obj_ids, self.entProps)
                    return self._load_graph(zeep_data, zeep_objects)
                else:
                    return None
            else:
                zeep_data, self.ResultPos = self.get_next_session_page(self.ResultRef, self.ResultPos, self.PageSize,
                                                                       self.ResultSize, self.entProps, getLinks=False)
                if type(zeep_data) != type(None):
                    return self._load_graph(None, zeep_data)
                else:
                    return None
        else:
            return None

    def __set_get_links(self):
        return self.GOQLquery[7:15] == 'Relation'

    def __replace_goql(self, oql_query: str):
        self.GOQLquery = oql_query
        self.__getLinks = self.__set_get_links()

    def add_rel_props(self, add_props: list):
        self.relProps = list(set(self.relProps + add_props))

    def add_ent_props(self, add_props: list):
        self.entProps = list(set(self.entProps + add_props))

    def add_graph(self, new_graph: ResnetGraph):
        pass  # placeholder to derive child class from APISession (see DiseaseNetwork(APISession) as example)

    def flash_dump_files(self):
        for f in self.DumpFiles:
            open(f, 'w').close()
            print('File "%s" was cleared before processing' % f)

    def add_dump_file(self, dump_fname, replace_main_dump=False):
        if replace_main_dump:
            self.DumpFiles = []
        self.DumpFiles.append(dump_fname)

    @staticmethod
    def execution_time(execution_start):
        return "{}".format(str(timedelta(seconds=time.time() - execution_start)))

    @staticmethod
    def reopen(fname):
        open(fname, "w", encoding='utf-8').close()
        return open(fname, "a", encoding='utf-8')

    def to_csv(self, file_out, in_graph: ResnetGraph, access_mode='w', debug=False):
        self.Graph.print_references(file_out, self.relProps, self.entProps, access_mode, 
                                    self.__IsOn1st_page, col_sep=self.csv_delimeter,debug=debug)

    def get_result_size(self):
        zeep_relations, (self.ResultRef, self.ResultSize, self.ResultPos) = self.init_session(self.GOQLquery, PageSize=1,
                                                                                              property_names=[])
        return self.ResultSize

    def process_oql(self, oql_query, request_name='', flash_dump=False, debug=False, no_mess=True, iteration_limit=1):
        global_start = time.time()
        if flash_dump:
            self.flash_dump_files()
        self.__replace_goql(oql_query)

        entire_graph = ResnetGraph()
        reference_counter = 0
        
        start_time = time.time()
        return_type = 'relations' if self.__getLinks else 'entities'
        page_graph = self.__init_session()
        number_of_iterations = math.ceil(self.ResultSize / self.PageSize)
        if number_of_iterations > 1:
            print('\n\"%s\"\nrequest found %d %s. Begin retrieval in %d iterations' % (request_name,self.ResultSize,return_type, number_of_iterations))
            if debug:
                print("Processing GOQL query:\n\"%s\"\n" % (self.GOQLquery))
            if no_mess:
                print('Progress report is suppressed. Retrieval may take long time - be patient!')

        while isinstance(page_graph, ResnetGraph):
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
            
            self.add_graph(page_graph)
            
            if len(self.DumpFiles) > 0:
                self.to_csv(self.DumpFiles[0], in_graph=page_graph, access_mode='a',debug=True)
                self.__IsOn1st_page = False
                if number_of_iterations > 1:
                    print("With %d in %d iterations, %d relations in %d with %d references saved into \"%s\" file. Retrieval time: %s" % 
                        (iteration,number_of_iterations,self.ResultPos,self.ResultSize,reference_counter,self.DumpFiles[0],exec_time))

            entire_graph = nx.compose(page_graph, entire_graph)

            if debug and number_of_iterations >= iteration_limit: break
            page_graph = self.__get_next_page(no_mess)
            

        if debug: print("GOQL query:\n \"%s\"\n was executed in %s in %d iterations" % 
                 (self.GOQLquery, self.execution_time(global_start), number_of_iterations))
            

        if len(request_name) > 0:
            if number_of_iterations == 1: print('\n') # this is the only message about retreival with one ietartion. 
            print('\"%s\" was retrieved in %s by %d iterations' % 
                 (request_name, self.execution_time(global_start), number_of_iterations))

        self.ResultRef = ''
        self.ResultPos = 0
        self.ResultSize = 0
        self.__IsOn1st_page = True
        return entire_graph

    def get_ppi_graph(self, fout):
        # Build PPI network between proteins
        graph_proteins = self.Graph.get_entity_ids(['Protein'])
        print('Retrieving PPI network for %d proteins' % (len(graph_proteins)))
        start_time = time.time()
        ppi_graph = self.get_ppi(graph_proteins, self.relProps, self.entProps)
        self.to_csv(fout, ppi_graph)
        print("PPI network with %d relations was retrieved in %s ---" % (
            ppi_graph.size(), self.execution_time(start_time)))
        return ppi_graph


    def iterate_oql(self, oql_query:str, id_set:set):
        to_return = ResnetGraph()
        to_return = self._iterate_oql(oql_query,id_set, self.relProps,self.entProps)
        return to_return

    def iterate_oql2(self, oql_query:str, id_set1:set, id_set2:set):
        to_return = ResnetGraph()
        to_return = self._iterate_oql2(oql_query,id_set1,id_set2,self.relProps,self.entProps)
        return to_return

    def to_pandas (self, in_graph=None, RefNumPrintLimit=0):
        return self.Graph.ref2pandas(self.relProps,self.entProps,in_graph,RefNumPrintLimit)

    def connect_nodes(self,node_ids1:set,node_ids2:set,by_relation_type=None,with_effect=None,in_direction=None):
        oql_query = r'SELECT Relation WHERE '
        if isinstance(in_direction,str) and in_direction == '>':
            oql_query = oql_query + 'NeighborOf downstream (SELECT Entity WHERE id = ({ids1})) AND NeighborOf upstream (SELECT Entity WHERE id = ({ids2}))'
        elif isinstance(in_direction,str) and in_direction == '<':
            oql_query = oql_query + 'NeighborOf upstream (SELECT Entity WHERE id = ({ids1})) AND NeighborOf downstream (SELECT Entity WHERE id = ({ids2}))'
        else:
            oql_query = oql_query + 'NeighborOf (SELECT Entity WHERE id = ({ids1})) AND NeighborOf (SELECT Entity WHERE id = ({ids2}))'

        if isinstance(with_effect,list): 
            oql_query = oql_query + ' AND Effect = ('+','.join(with_effect)+')'

        if isinstance(by_relation_type,list):
            oql_query = oql_query + ' AND objectType = ('+','.join(by_relation_type)+')'

        return self.iterate_oql2(f'{oql_query}',node_ids1,node_ids2)


    

    

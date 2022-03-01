import time
import math
import networkx as nx
import pandas as pd
from datetime import timedelta
from .ZeepToNetworkx import PSNetworx
from .ResnetGraph import ResnetGraph
from .PathwayStudioGOQL import OQL
from ..ETM_API.references import SENTENCE_PROPS
from xml.dom import minidom

class APISession(PSNetworx):
    pass
    ResultRef = None
    ResultPos = 0
    ResultSize = None
    PageSize = 100
    GOQLquery = str()
    __IsOn1st_page = True
    relProps = ['Name', 'RelationNumberOfReferences']
    entProps = ['Name']
    DumpFiles = ['ResnetAPIsessionDump.tsv']
    csv_delimeter = '\t'
    print_rel21row = False
    clear_graph_cache = False

    def __init__(self,url, username, password):
        super().__init__(url, username, password)
        self.DumpFiles = ['ResnetAPIsessionDump.tsv']

    def __init_session(self):
        obj_props = self.relProps if self.__getLinks else self.entProps
        first_iteration = self.ResultPos
        zeep_data, (self.ResultRef, self.ResultSize, self.ResultPos) = self.init_session(self.GOQLquery,
                                                                                                self.PageSize,
                                                                                                obj_props,
                                                                                                getLinks = self.__getLinks)

        if first_iteration > 0:
            self.ResultPos = first_iteration
            return self.__get_next_page()                                                  

        if type(zeep_data) != type(None):
            if self.__getLinks:
                obj_ids = list(set([x['EntityId'] for x in zeep_data.Links.Link]))
                zeep_objects = self.get_object_properties(obj_ids, self.entProps)
                return self._load_graph(zeep_data, zeep_objects)
            else:
                return self._load_graph(None, zeep_data)
        else:
            return ResnetGraph()

    def clear(self):
        self.Graph.clear()
        self.IDtoRelation.clear()
        self.ID2Children.clear()
        
    def start_download_from(self,result_position:int):
        self.ResultPos = result_position

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
                    return self._load_graph(zeep_data, zeep_objects)
                else:
                    return self._load_graph(None, zeep_data)
            else:
                return ResnetGraph()
        else: return ResnetGraph()


    def get_saved_results(self, results_names:list, get_links=True):
        query_names = OQL.join_with_quotes(',',results_names)
        oql_query = 'SELECT Result WHERE Name = ({names})'
        oql_query = oql_query.format(names=query_names)
        result_graph = self.load_graph_from_oql(oql_query,list(self.relProps),list(self.entProps),get_links=False)
        #result_ids = self._obj_id_by_oql(oql_query)
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

    def __set_get_links(self):
        return self.GOQLquery[7:15] == 'Relation'

    def __replace_goql(self, oql_query: str):
        self.GOQLquery = oql_query
        self.__getLinks = self.__set_get_links()
        self.ResultRef = None

    def add_rel_props(self, add_props:list):
        self.relProps = self.relProps+[i for i in add_props if i not in self.relProps]
        if set(SENTENCE_PROPS).difference(add_props) and 'TextRef' not in self.relProps:
            self.relProps.append('TextRef')

    def add_ent_props(self, add_props: list):
        self.entProps = self.entProps+[i for i in add_props if i not in self.entProps]

    def add_graph(self, new_graph: ResnetGraph):
        pass  # placeholder to derive child class from APISession (see DiseaseNetwork(APISession) as example)

    def flush_dump_files(self, no_mess=True):
        for f in self.DumpFiles:
            open(f, 'w').close()
            if not no_mess:
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

    def to_csv(self, file_out, in_graph: ResnetGraph=None, access_mode='w', debug=False):
        if not isinstance(in_graph,ResnetGraph): in_graph = self.Graph
        in_graph.print_references(file_out, self.relProps, self.entProps, access_mode, self.__IsOn1st_page,
                                    col_sep=self.csv_delimeter,debug=debug,single_rel_row=self.print_rel21row)

    def get_result_size(self,oql_query):
        zeep_relations, (self.ResultRef, self.ResultSize, self.ResultPos) = self.init_session(oql_query, PageSize=1,
                                                                                              property_names=[])
        return self.ResultSize

    def process_oql(self, oql_query, request_name='', flush_dump=False, debug=False, no_mess=True, iteration_limit=1) -> ResnetGraph:
        global_start = time.time()
        if flush_dump and self.ResultPos == 0:
            self.flush_dump_files(no_mess)
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
            
            self.add_graph(page_graph)
            
            if len(self.DumpFiles) > 0:
                self.to_csv(self.DumpFiles[0], in_graph=page_graph, access_mode='a',debug=True)
                self.__IsOn1st_page = False
                if number_of_iterations > 1:
                    print("With %d in %d iterations, %d %s in %d with %d references saved into \"%s\" file. Retrieval time: %s" % 
                        (iteration,number_of_iterations,self.ResultPos,return_type,self.ResultSize,reference_counter,self.DumpFiles[0],exec_time))

            entire_graph = nx.compose(page_graph, entire_graph)

            if debug and number_of_iterations >= iteration_limit: break
            if self.clear_graph_cache: self.clear()
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
        self.to_csv(fout, in_graph=ppi_graph)
        print("PPI network with %d relations was retrieved in %s ---" % (
            ppi_graph.size(), self.execution_time(start_time)))
        return ppi_graph

    def get_network_graph(self, InteractorIdList:set, connect_by_rel_types:list=None):
        return self.get_network(InteractorIdList,connect_by_rel_types, self.relProps,self.entProps)


    def iterate_oql(self, oql_query:str, id_set:set):
        to_return = ResnetGraph()
        to_return = self._iterate_oql(oql_query,id_set, self.relProps,self.entProps)
        return to_return


    def iterate_oql2(self, oql_query:str, id_set1:set, id_set2:set):
        to_return = ResnetGraph()
        to_return = self._iterate_oql2(oql_query,id_set1,id_set2,self.relProps,self.entProps)
        return to_return
        

    def to_pandas (self, in_graph=None, RefNumPrintLimit=0)-> 'pd.DataFrame':
        if not isinstance(in_graph, ResnetGraph): in_graph = self.Graph
        return in_graph.ref2pandas(self.relProps,self.entProps,RefNumPrintLimit)

    def connect_nodes(self,node_ids1:set,node_ids2:set,by_relation_type:list=None,with_effect:list=None,in_direction:str=None):
        oql_query = r'SELECT Relation WHERE '
        if isinstance(in_direction,str):
            if in_direction == '>':
                oql_query = oql_query + 'NeighborOf downstream (SELECT Entity WHERE id = ({ids1})) AND NeighborOf upstream (SELECT Entity WHERE id = ({ids2}))'
            elif in_direction == '<':
                oql_query = oql_query + 'NeighborOf upstream (SELECT Entity WHERE id = ({ids1})) AND NeighborOf downstream (SELECT Entity WHERE id = ({ids2}))'
            else: print('Direction sign is not recorgnized!')
        else:
            oql_query = oql_query + 'NeighborOf (SELECT Entity WHERE id = ({ids1})) AND NeighborOf (SELECT Entity WHERE id = ({ids2}))'

        if isinstance(with_effect,list): 
            oql_query = oql_query + ' AND Effect = ('+','.join(with_effect)+')'

        if isinstance(by_relation_type,list):
            oql_query = oql_query + ' AND objectType = ('+','.join(by_relation_type)+')'

        return self.iterate_oql2(f'{oql_query}',node_ids1,node_ids2)

    def get_group_members(self, group_names:list):
        oql_query = 'SELECT Entity WHERE MemberOf (SELECT Group WHERE Name = ({group}))'
        groups = OQL.join_with_quotes(',', group_names)
        req_name = 'Find members of groups: ' + ','.join(group_names)
        graph2return = self.process_oql(oql_query.format(group=groups), request_name=req_name)
        if len(graph2return) == 0:
            print('%s groups are empty or do not exist in databse' % str(group_names))
        else:
            print('loaded %d members from %s' % (graph2return.number_of_nodes(),str(group_names)))
        return graph2return


    def map_props2objs(self, propValues:list, prop_names:list,map2types=None,case_insensitive=False):
        propval2objs,objid2propval = self.Graph.get_props2obj_dic(propValues, prop_names,case_insensitive)
        need_db_mapping = set(propValues).difference(propval2objs.keys())
        only_object_types = [] if map2types is None else map2types
        self.add_ent_props(prop_names)

        step = 1000
        iteration_counter = math.ceil(len(need_db_mapping) / step)
        print('Will use %d %s identifiers to find entities in %d iterations' % 
             (len(need_db_mapping), ','.join(prop_names), iteration_counter))

        for i in range(0, len(need_db_mapping), step):
            propval_chunk = propValues[i:i + step]
            query_node = OQL.get_entities_by_props(propval_chunk, prop_names, only_object_types)
            self.process_oql(query_node)

        p2obj,id2p = self.Graph.get_props2obj_dic(need_db_mapping, prop_names,case_insensitive)
        propval2objs.update(p2obj)
        for id,prop_vals in id2p.items():
            try:
                mapped_values = set(objid2propval[id])
                mapped_values.update(map(lambda x:str(x).lower(),prop_vals))
                objid2propval[id] = list(mapped_values)
            except KeyError:
                objid2propval[id] = prop_vals

        return propval2objs, objid2propval


    @staticmethod
    def pretty_xml(xml_string:str, no_declaration = False):
        pretty_xml = str(minidom.parseString(xml_string).toprettyxml(indent='   '))
        if no_declaration:
            pretty_xml = pretty_xml[pretty_xml.find('\n')+1:]
        
        return pretty_xml





        

    

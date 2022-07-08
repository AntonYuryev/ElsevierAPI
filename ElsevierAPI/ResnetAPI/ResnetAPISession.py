import time
import math
import networkx as nx
from datetime import timedelta
from .ZeepToNetworkx import PSNetworx
from .ResnetGraph import ResnetGraph,PSObject,PSRelation,CHILDS,REFCOUNT
from .PathwayStudioGOQL import OQL
from .Zeep2Experiment import Experiment,Sample
from ..ETM_API.references import SENTENCE_PROPS
from xml.dom import minidom
from ..pandas.panda_tricks import df
import glob
from pathlib import Path
import os


class APISession(PSNetworx):
    pass
    ResultRef = None
    ResultPos = 0
    ResultSize = None
    PageSize = 100
    GOQLquery = str()
    __IsOn1st_page = True
    relProps = ['Name','RelationNumberOfReferences','URN']
    entProps = ['Name']
    DumpFiles = ['ResnetAPIsessionDump.tsv']
    csv_delimeter = '\t'
    print_rel21row = False
    clear_graph_cache = False
    __getLinks = True
    APIconfig = dict()
    data_dir = ''

    def __init__(self,url, username, password):
        super().__init__(url, username, password)
        self.APIconfig['ResnetURL'] = url
        self.APIconfig['PSuserName'] = username
        self.APIconfig['PSpassword'] = password
        self.DumpFiles = []


######################################  CONFIGURATION ######################################
    @classmethod
    def fromAPIconfig(cls,APIconfig:dict):
        session =  cls(APIconfig['ResnetURL'], APIconfig['PSuserName'], APIconfig['PSpassword'])
        session.APIconfig = APIconfig
        return session

    def __set_get_links(self):
        return self.GOQLquery[7:15] == 'Relation'

    def __replace_goql(self, oql_query: str):
        self.GOQLquery = oql_query
        self.__getLinks = self.__set_get_links()
        self.ResultRef = None


    def start_download_from(self,result_position:int):
        self.ResultPos = result_position


    def add_rel_props(self, add_props:list):
        self.relProps = self.relProps+[i for i in add_props if i not in self.relProps]
        if set(SENTENCE_PROPS).difference(add_props) and 'TextRef' not in self.relProps:
            self.relProps.append('TextRef')


    def add_ent_props(self, add_props: list):
        self.entProps = self.entProps+[i for i in add_props if i not in self.entProps]

    
    def add_dump_file(self, dump_fname, replace_main_dump=False):
        if replace_main_dump:
            self.DumpFiles = []
        self.DumpFiles.append(dump_fname)

#########################  PAGE BY PAGE PROCESSING, RETRIEVAL   ##################################
    def __init_session(self):
        self.__set_get_links()
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
            
            self.Graph.add_graph(page_graph)
            
            if len(self.DumpFiles) > 0:
                self.to_csv(self.DumpFiles[0], in_graph=page_graph, access_mode='a',debug=True)
                self.__IsOn1st_page = False
                if number_of_iterations > 1:
                    print("With %d in %d iterations, %d %s in %d with %d references saved into \"%s\" file. Retrieval time: %s" % 
                        (iteration,number_of_iterations,self.ResultPos,return_type,self.ResultSize,reference_counter,self.DumpFiles[0],exec_time))
            else:
                if number_of_iterations > 1:
                    print("With %d in %d iterations, %d %s in %d with %d references retrieved in: %s" % 
                        (iteration,number_of_iterations,self.ResultPos,return_type,self.ResultSize,reference_counter,exec_time))


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

    
    def iterate_oql(self, oql_query:str, id_set:set):
        """
        # oql_query MUST contain string placeholder called {ids}
        """
        to_return = ResnetGraph()
        to_return = self._iterate_oql(oql_query,id_set, self.relProps,self.entProps)
        return to_return


    def iterate_oql2(self, oql_query:str, id_set1:set, id_set2:set):
        """
        # oql_query MUST contain 2 string placeholders called {ids1} and {ids2}
        """
        to_return = ResnetGraph()
        to_return = self._iterate_oql2(oql_query,id_set1,id_set2,self.relProps,self.entProps)
        return to_return
        

    def clear(self):
        self.Graph.clear()
        self.IDtoRelation.clear()
        self.ID2Children.clear()


    def flush_dump_files(self, no_mess=True):
        for f in self.DumpFiles:
            open(f, 'w').close()
            if not no_mess:
                print('File "%s" was cleared before processing' % f)

    @staticmethod
    def execution_time(execution_start):
        return "{}".format(str(timedelta(seconds=time.time() - execution_start)))

    @staticmethod
    def reopen(fname):
        open(fname, "w", encoding='utf-8').close()
        return open(fname, "a", encoding='utf-8')

    def get_result_size(self,oql_query):
        zeep_relations, (self.ResultRef, self.ResultSize, self.ResultPos) = self.init_session(oql_query, PageSize=1,
                                                                                              property_names=[])
        return self.ResultSize



###########################GET GET RETRIEVE GET RETRIEVE GET GET ################################################
    def get_saved_results(self, results_names:list, get_links=True):
        query_names = OQL.join_with_quotes(results_names)
        oql_query = 'SELECT Result WHERE Name = ({names})'
        oql_query = oql_query.format(names=query_names)
        result_graph = self.load_graph_from_oql(oql_query,list(self.relProps),list(self.entProps),get_links=False)
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


    def get_ppi_graph(self, fout):
        # Build PPI network between proteins
        graph_proteins = self.Graph.get_node_ids(['Protein'])
        print('Retrieving PPI network for %d proteins' % (len(graph_proteins)))
        start_time = time.time()
        ppi_graph = self.get_ppi(graph_proteins, self.relProps, self.entProps)
        self.to_csv(self.data_dir+fout, in_graph=ppi_graph)
        print("PPI network with %d relations was retrieved in %s ---" % (
            ppi_graph.size(), self.execution_time(start_time)))
        return ppi_graph


    def get_network_graph(self, for_interactor_ids:set, connect_by_rel_types:list=None):
        return self.get_network(for_interactor_ids,connect_by_rel_types, self.relProps,self.entProps)


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
        groups = OQL.join_with_quotes( group_names)
        req_name = 'Find members of groups: ' + ','.join(group_names)
        graph2return = self.process_oql(oql_query.format(group=groups), request_name=req_name)
        if len(graph2return) == 0:
            print('%s groups are empty or do not exist in databse' % str(group_names))
        else:
            print('loaded %d members from %s' % (graph2return.number_of_nodes(),str(group_names)))
        return graph2return


    def map_props2objs(self, using_values:list, in_properties:list,map2types=[],case_insensitive=False):
        """
        returns:
            propval2objs = {prop_value:[PSObject]}
            objid2propval = {node_id:[prop_values]}
        where prop_value is from 'using_values'
        """
        propval2objs,objid2propval = self.Graph.get_props2obj_dic(using_values, in_properties,case_insensitive)
        need_db_mapping = set(using_values).difference(propval2objs.keys())
        self.add_ent_props(in_properties)

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
        return propval2objs, objid2propval

################################# ONTOLOGY  ONTOLOGY ############################################
    def child_graph(self, propValues:list, search_by_properties=[],include_parents=True):
        if not search_by_properties: search_by_properties = ['Name','Alias']
        oql_query = OQL.get_childs(propValues,search_by_properties,include_parents=include_parents)
        request_name = 'Find ontology children for {}'.format(','.join(propValues))
        ontology_graph = self.process_oql(oql_query,request_name)
        print('Found %d ontology children' % len(ontology_graph))
        return ontology_graph


    def load_children(self,in_graph=ResnetGraph(),for_parent_ids=[],add_new_nodes=False):
        """
        slower than child_graph but annotates each parent with CHILDS property
        """
        graph_registry = [in_graph, self.Graph] #need to pass as pointers because graph is modified
        index = 0 if in_graph else 1

        annotate_parent_ids = set(for_parent_ids) if for_parent_ids else set(graph_registry[index])
        self.get_children(annotate_parent_ids)

        if add_new_nodes:
            [graph_registry[index].property2node(parent_id,CHILDS,self.ID2Children[parent_id]) for parent_id in annotate_parent_ids]
        else:
            for parent_id in annotate_parent_ids:
                all_child_ids = set(self.ID2Children[parent_id])
                graph_child_ids = all_child_ids.intersection(set(graph_registry[index]))
                graph_registry[index].property2node(parent_id,CHILDS,list(graph_child_ids))


    def add2ontology_graph(self, members:list, to_parent:PSObject):
        """
        members - list of PSObjects
        """
        ontology_graph = ResnetGraph()
        parent_id = to_parent['Id'][0]
        ontology_graph.add1node(to_parent)
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
        #self.add_graph(ontology_graph)
        return ontology_graph


    def add_parents(self, to_graph=ResnetGraph(), for_child_ids=[], depth=1):
        """
        Adds ontology parents to the graph
        """
        if not for_child_ids:
            if to_graph:
                for_child_ids = list(to_graph)
            else:
                for_child_ids = list(self.Graph)

        get_parent_query = 'SELECT Entity WHERE InOntology (SELECT Annotation WHERE Ontology=\'Pathway Studio Ontology\' AND Relationship=\'is-a\') inRange {steps} over (SELECT OntologicalNode WHERE id = ({ids}))'
        oql_query = get_parent_query.format(steps=str(depth),ids=','.join(map(str,for_child_ids)))
        request_name = 'Find parents of {count} nodes with depth {d}'.format(count = str(len(for_child_ids)), d=str(depth))
        parent_graph = self.process_oql(oql_query,request_name)

        if to_graph:
            to_graph.add_graph(parent_graph)
            self.load_children(to_graph,for_parent_ids=list(parent_graph))
            return to_graph._get_node(list(parent_graph))
        else:
            self.load_children(for_parent_ids=list(parent_graph))
            return self.Graph._get_nodes(list(parent_graph))



###################################   WRITE, DUMP   ##############################################
    def to_csv(self, file_out, in_graph: ResnetGraph=None, access_mode='w', debug=False):
        if not isinstance(in_graph,ResnetGraph): in_graph = self.Graph
        in_graph.print_references(file_out, self.relProps, self.entProps, access_mode, self.__IsOn1st_page,
                                    col_sep=self.csv_delimeter,debug=debug,single_rel_row=self.print_rel21row)


    def to_pandas (self, in_graph=None, RefNumPrintLimit=0)-> 'df':
        if not isinstance(in_graph, ResnetGraph): in_graph = self.Graph
        return df(in_graph.ref2pandas(self.relProps,self.entProps,RefNumPrintLimit))


    @staticmethod
    def pretty_xml(xml_string:str, no_declaration = False):
        pretty_xml = str(minidom.parseString(xml_string).toprettyxml(indent='   '))
        if no_declaration:
            pretty_xml = pretty_xml[pretty_xml.find('\n')+1:]
        return pretty_xml


    def graph2rnef(self,fname:str, in_graph=ResnetGraph(), resnet_size=1000):
        if not in_graph: in_graph = self.Graph
        relation_count = in_graph.number_of_edges()
        print('Dumping graph with %d edges and %d nodes into %s file' % 
                (relation_count,in_graph.number_of_nodes(),self.data_dir+fname) )
        
        resnets_cnt = math.ceil(relation_count/resnet_size)
        print('%s cache file will have %d <resnet> sections with %d relations' % 
                (self.data_dir+fname, resnets_cnt, resnet_size))

        with open(self.data_dir+fname, 'w', encoding='utf=8') as f:
            resnet = ResnetGraph()
            f.write('<?xml version="1.0" encoding="UTF-8" standalone="no" ?>\n<batch>\n')
            for regulatorID, targetID, e in in_graph.edges(data='relation'):
                resnet.copy_rel(e,in_graph)
                if resnet.number_of_edges() == resnet_size:
                    rnef_str = resnet.to_rnef(ent_props=self.entProps,rel_props=self.relProps)
                    rnef_str = self.pretty_xml(rnef_str,no_declaration=True)
                    f.write(rnef_str)
                    resnet.clear()

            rnef_str = resnet.to_rnef(ent_props=self.entProps,rel_props=self.relProps)
            rnef_str = self.pretty_xml(rnef_str,no_declaration=True)
            f.write(rnef_str)
            f.write('</batch>')


    @staticmethod
    def __make_file_name(folder_or_object_name: str):
        new_str = list(folder_or_object_name)
        for i in range(0, len(new_str)):
            if new_str[i] in {'>', '<', '|','/'}:
                new_str[i] = '-'
            if new_str[i] in {':'}:
                new_str[i] = '_'
        return "".join(new_str)


    def find_folder(self, folder_name:str):
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
        my_dir = self.find_folder(of_folder)
        if not my_dir:
            if in2parent_folder:
                parent_dir = self.find_folder(in2parent_folder)
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


    def count_files(self, in_folder:str,in2parent_folder='', with_extension='rnef'):
        folder_path = self.dump_path(in_folder,in2parent_folder)
        listing = glob.glob(folder_path+'*.' + with_extension)
        if listing:
            return folder_path,len(listing)
        else:
            return folder_path,1

    def dump_base(self,of_folder:str):
        return 'content of '+ APISession.__make_file_name(of_folder)+'_'

    def dumpfile(self,of_folder:str,in2parent_folder='',new=False):
        folder_path,file_count = self.count_files(of_folder,in2parent_folder)
        if new:
            file_count += 1
        return folder_path + self.dump_base(of_folder)+str(file_count)+'.rnef'

    def dump_rnef(self,rnef_xml:str,to_folder='',in2parent_folder='',max_rnef_size=100000000,can_close=True):
        write2 = self.dumpfile(to_folder,in2parent_folder)
        if Path(write2).exists():
            f = open(write2,'a',encoding='utf-8')
            f.write(rnef_xml)
            file_size = os.path.getsize(write2)
            if file_size > max_rnef_size and can_close:
                f.write('</batch>')
                f.close()

                new_write2 = self.dumpfile(to_folder,in2parent_folder,new=True)
                f = open(new_write2,'w',encoding='utf-8')
                f.write('<batch>\n')
            f.close()
        else:
            with open(write2,'w',encoding='utf-8') as f:
                f.write('<batch>\n'+rnef_xml)


    def close_rnef_dump(self,to_folder='',in2parent_folder=''):
        last_dump_file = self.dumpfile(to_folder,in2parent_folder)
        f = open(last_dump_file,'a',encoding='utf-8')
        f.write('</batch>')
        f.close()


    
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


            

             


    
        
        




        



        



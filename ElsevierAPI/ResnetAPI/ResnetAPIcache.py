import  os, glob, json
from shutil import copyfile
from .ResnetAPISession import REFERENCE_IDENTIFIERS, REFCOUNT
from .ResnetAPISession import APISession,ResnetGraph,PSObject,PSRelation
from .ResnetGraph import EFFECT
from contextlib import redirect_stdout

CACHE_DIR = 'D:/Python/ENTELLECT_API/ElsevierAPI/ResnetAPI/__pscache__/'

class APIcache(APISession):
    '''
    manages work with cache files.\n Reads cache into ResnetGraph self.network
    self.Graph is used for download cache and emptied after self.network loading
    '''
    pass
    reference_cache_size = 1000000 # max number of reference allowed in self.Graph. Clears self.Graph if exceeded
    resnet_size = 1000 # number of <node><control> sections in RNEF dump
    max_rnef_size = 100000000 # max size of RNEF XML dump file. If dump file exceeds max_file_size new file is opened with index++
    example_node = PSObject()

######################################  CONFIGURATION  ######################################
    def __init__(self,*args,**kwargs):
        '''
        Input
        -----
        APIconfig = args[0]\n
        what2retrieve - defaults REFERENCE_IDENTIFIERS, other options:\n
        [NO_REL_PROPERTIES,DATABASE_REFCOUNT_ONLY,REFERENCE_IDENTIFIERS,BIBLIO_PROPERTIES,SNIPPET_PROPERTIES,ONLY_REL_PROPERTIES,ALL_PROPERTIES]
        \nno_mess - default True, if False script becomes more verbose
        connect2server - default True, set to False to run script using data in __pscache__ files instead of database
        '''
        
        default_cache_name = 'Resnet subset'
        self.network = ResnetGraph()
        my_kwargs ={
            'what2retrieve' : REFERENCE_IDENTIFIERS,
            'data_dir' : CACHE_DIR,
            'cache_name': default_cache_name, # used for folder in self.data_dir and in RNEF dump name
            'oql_queries':[], # [(oql_query,request_name),...] list of GOQL queries tuples for fetching data from database
            'connect_nodes' : [], # [[entity types or dbids],[relation types]]
            'ent_props' : ['Name'], # additional entity props to fetch from database
            'rel_props' : ['URN',EFFECT], # additional relation props to fetch from database
            'make_simple': True, # if True makes simple graph from database graph
            'ranks4simplifying' : [], # parameter for making simple graph
            'refprop2rel':dict(), # {reference_propname:relation_propname} - moves refprops to relprops using PSRelation._refprop2rel()
            'refprop_minmax':0, # parameter for PSRelation._refprop2rel(refprop2rel,refprop_minmax)
            #if min_max < 0 assigns single value to relprop that is the minimum of all ref_prop values
            #if min_max > 0 assigns single value to relprop that is the maximum of all ref_prop values
            #min_max works only if ref_prop values are numerical
            'predict_effect4' : [], # [_4enttypes:list,_4reltypes:list] - parameters for ResnetGraph.predict_effect4()
            'no_id_version': True,
            'max_threads' : 25 # controls download speed.  Make it 10 if what2retrieve=ALL_PROPERTIES
        }

        ent_props = list(kwargs.pop('ent_props',[]))
        rel_props = list(kwargs.pop('rel_props',[]))
        my_kwargs.update(kwargs)

        self.cache_was_modified = False

        super().__init__(*args, **my_kwargs)
        # properties have to added after super().__init__ adds properies from what2retrieve
        super().add_ent_props(ent_props)
        super().add_rel_props(rel_props)
        
        my_kwargs['ent_props'] = ent_props
        my_kwargs['rel_props'] = rel_props # to pass for printing
        cache_name = my_kwargs.pop('cache_name')
        
        self.relprops2rnef = list(rel_props)
        self.relprops2rnef += [p for p in [EFFECT,'URN',REFCOUNT] if p not in rel_props]
        
        # loads cache_name from self.data_dir:
        load_cache = my_kwargs.pop('load_cache',True)
        if load_cache:
            self.network = self._load_cache(cache_name,**my_kwargs)
            self.save_description()
            super().clear() # clearing self.Graph to save RAM
  

    def __path2cache(self,cache_name:str,extension='.rnef'):
        return self.data_dir+cache_name+extension


    def clean_graph(self,database_g:ResnetGraph,**kwargs)->'ResnetGraph':
        # filter to load relations only for nodes with desired properties\n
        rank4simplifying = kwargs.pop('ranks4simplifying',[])
        predict_effect4 = kwargs.pop('predict_effect4',dict()) #kwargs for ResnetGraph.predict_effect4(): _4enttypes, _4reltypes
        refprop2rel = dict(kwargs.pop('refprop2rel',dict()))
        refprop_minmax = kwargs.pop('refprop_minmax',0)

        clean_graph = database_g.copy()

        if kwargs['make_simple']:
            clean_graph = clean_graph.make_simple(rank4simplifying)

        #converting sentence properties to relation properties before making a dump.  Used for pX dump
        for refprop, relprop in refprop2rel.items():
            clean_graph.refprop2rel(refprop,relprop,refprop_minmax)
            self.relprops2rnef.remove(refprop)
            self.relprops2rnef.append(relprop)
            
        if predict_effect4:
            clean_graph,modified_rels = clean_graph.predict_effect4(**predict_effect4)
            modified_rels_graph = clean_graph.subgraph_by_rels(modified_rels)
            modified_rels_graph.name = 'relations with predicted effect'
            modified_rels_graph.dump2rnef('Effect predicted4',self.entProps,self.relprops2rnef)
        
        if kwargs['no_id_version']:
            for ent_prop in self.entProps:
                if ent_prop[-3:] == ' ID':
                    clean_graph = clean_graph.clean_version_number(ent_prop)

        return clean_graph
    

    def save_network(self, example_node_name=''):
        if example_node_name:
            self.example_node = self.network._psobjs_with(example_node_name)[0]
        if self.cache_was_modified:
            self.replace_cache(self.network.name,self.network,self.entProps,self.relprops2rnef)
            self.save_description()
        return
    

    def save_description(self):
        network_description = {'Name':[self.network.name],'OQLs':self.my_oql_queries,
                                   'nodes':self.network.number_of_nodes(),
                                   'edges':self.network.number_of_edges(),
                                   'diseases': len(self.network._psobjs_with('Disease','ObjTypeName')),
                                   'chemicals': len(self.network._psobjs_with('SmallMol','ObjTypeName')),
                                   'proteins': len(self.network._psobjs_with('Protein','ObjTypeName')),
                                   'DirectRegulation': len(self.network.psrels_with(['DirectRegulation'])),
                                   'Binding': len(self.network.psrels_with(['Binding'])),
                                   'Regulation': len(self.network.psrels_with(['Regulation']))
                                   }
        
        with open(self.data_dir+self.network.name+"_description.json", "w") as f:
            json.dump(network_description, f, indent=2)


    def _load_cache(self,cache_name, **kwargs)->'ResnetGraph':
        """
        Input
        -----

        oql_queries - list of tuples (oql_query,request_name)\n
        prop2values={prop_name:[values]} - filter to load relations only for nodes with desired properties\n
        predict_effect4 - kwargs for ResnetGraph.predict_effect4()
        refprop2rel - {refPropName:NewRelPropName}
        "ranks4simplifying" is used to merge relations during cache network simplification\n
        relation types must be ordered from most important type to

        if refprop_minmax < 0 assigns single value to relprop that is the minimum of all ref_prop values
        if refprop_minmax > 0 assigns single value to relprop that is the maximum of all ref_prop values

        Examples:\n
        ['DirectRegulation','Binding','ProtModification','Regulation']\n
        ['PromoterBinding','Expression','Regulation']\n
        ['Regulation','Biomarker','StateChange','QuantitativeChange','FunctionalAssociation']\n
        ['Biomarker','StateChange','FunctionalAssociation']\n
        ['Biomarker','QuantitativeChange','FunctionalAssociation']\n
        [MolSynthesis','Regulation','FunctionalAssociation']\n
        [MolTransport','Regulation','FunctionalAssociation']\n
        [MolTransport','CellExpression']\n
        ['Regulation','FunctionalAssosiation']\n

        Return
        ------
        ResnetGraph with name = cache_name
        """
        print(f'Loading {cache_name} cache')
        prop2values = kwargs.pop('prop2values',dict()) #prop2values={prop_name:[values]}
        # filter to load relations only for nodes with desired properties\n

        my_cache_file = self.__path2cache(cache_name)
        try:
            cached_graph = ResnetGraph.fromRNEF(my_cache_file,prop2values=prop2values)
            print(f'Loaded "{cache_name}" cache with {len(cached_graph)} nodes and {cached_graph.number_of_edges()} edges')
            cached_graph.name = cache_name
            return cached_graph
        except FileNotFoundError:
            print(f'Cannot find {my_cache_file} cache file')
            self.dump_folder = cache_name + '_raw'
            dump_dir = self.data_dir+self.dump_folder+'/'
            database_graph = ResnetGraph.fromRNEFdir(dump_dir,merge=False)
            if not database_graph:
                print(f'Cache "{dump_dir}" was not found\nBegin caching network from database')
                my_session_kwargs = {
                    'connect2server':True,
                    'no_mess' : self.no_mess
                    }
                my_session = self._clone_session(**my_session_kwargs) #need identifiers to make graph simple

                log_name = f'Download of {cache_name}.log'
                log_path = self.data_dir+log_name
                with open(log_path, 'w', buffering=1) as fout:
                    print(f'Download log is in {log_path}')
                    with redirect_stdout(fout):
                        if kwargs['connect_nodes']:
                            node_ids = set(kwargs['connect_nodes'][0])
                            # node_ids can be either {str} or {int}
                            if isinstance(next(iter(node_ids), None),str):
                                node_types_str = ",".join(node_ids)
                                query = f'SELECT Entity WHERE objectType = ({node_types_str})'
                                nodes_graph = my_session.process_oql(query,'Loading nodes from database')
                                if isinstance(nodes_graph,ResnetGraph):
                                    node_ids = set(ResnetGraph.dbids(nodes_graph.psobjs_with()))
                                else:
                                    print(f'No entities with {node_types_str} exist in database')
                                    node_ids = set()
                            
                            rel_types = kwargs['connect_nodes'][1]
                            my_session.get_network(node_ids,rel_types,download=True,threads=kwargs['max_threads'])
                        else:
                            for i, query in enumerate(self.my_oql_queries):
                                my_session.download_oql(query[0],query[0],threads=kwargs['max_threads'])


                my_session.clear()
                database_graph = ResnetGraph.fromRNEFdir(dump_dir,merge=False)
            
            database_graph.name = cache_name
            database_graph = self.clean_graph(database_graph,**kwargs)
            
            database_graph.dump2rnef(my_cache_file,self.entProps,self.relprops2rnef,with_section_size=1000) #with_section_size=1000
            
            
            print('%s with %d edges and %d nodes was written into %s' 
                % (database_graph.name,database_graph.number_of_edges(),database_graph.number_of_nodes(),my_cache_file))
            
            return database_graph


    def replace_cache(self,cache_name:str,with_network:ResnetGraph,
                      ent_props=['Name'],rel_props=['URN'],
                      do_backup=True,replace_raw=False):
        if replace_raw:
            my_cache_dir = self.data_dir+cache_name
            listing = glob.glob(my_cache_dir+'*.rnef')
            [os.remove(file) for file in listing]
            self._dump2rnef(with_network,cache_name)

        my_cache_file = self.__path2cache(cache_name)
        if do_backup:
            backup_cache = self.__path2cache(cache_name,'_backup.rnef')
            copyfile(my_cache_file, backup_cache)
        
        with_network.dump2rnef(my_cache_file,ent_props,rel_props)
        if self.example_node:
            example_cache = with_network.neighborhood({self.example_node})
            example_file = self.__path2cache(cache_name,extension='_example.rnef')
            example_cache.name = self.example_node.name()+' example'
            example_cache.dump2rnef(example_file,ent_props,rel_props)


    def clear(self):
        super().clear()
        self.network.clear_resnetgraph()


    def dump_subgraph_by_relprops(self,_2file:str,search_values:list,in_properties:list=['ObjTypeName']):
        '''
        Dumps
        ------
        subgraph with relations having search_values in_properties identified by ResnetGraph.subgraph_by_relprops
        into self.data_dir/_2file,
        if "search_values" is empty dumps subgraph with all relations annotated by "in_properties"
        '''
        my_subgraph = self.network.subgraph_by_relprops(search_values,in_properties)
        my_subgraph.dump2rnef(self.data_dir+_2file,self.entProps,self.relprops2rnef)


    def add2cache(self,graph:ResnetGraph,**kwargs):
        print(f'Adding {graph.number_of_nodes()} and {graph.number_of_edges()} edges to {self.network.name} cache')
        add2raw = kwargs.pop('add2raw',False)
        if add2raw:
           self._dump2rnef(graph,self.network.name)
        else:
            cache_name = self.network.name
            combined_network = graph.compose(self.network)
            print(f'New cache has {combined_network.number_of_nodes()} nodes, {combined_network.number_of_edges()} edges')
            self.network = self.clean_graph(combined_network,**kwargs)
            self.network.name = cache_name
            self.cache_was_modified = True
            self.save_network() #with_section_size=1000

            
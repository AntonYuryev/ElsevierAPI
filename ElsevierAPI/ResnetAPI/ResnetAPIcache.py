import os,glob,zipfile
from pathlib import Path
from shutil import copyfile
from .ResnetAPISession import REFERENCE_IDENTIFIERS, REFCOUNT
from .ResnetAPISession import APISession,ResnetGraph,PSObject
from .ResnetGraph import EFFECT
from contextlib import redirect_stdout

CACHE_DIR = 'D:/Python/ENTELLECT_API/ElsevierAPI/ResnetAPI/__pscache__/'
DEFAULT_CACHE_NAME = 'Resnet subset'

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
        self.network = ResnetGraph()
        default_kwargs ={
            'what2retrieve' : REFERENCE_IDENTIFIERS,
            'data_dir' : CACHE_DIR,
            'cache_name': DEFAULT_CACHE_NAME, # used for folder in self.data_dir and in RNEF dump name
            'oql_queries':[], # [(oql_query,request_name,resume_from_page),] 
            # list of GOQL queries tuples for fetching data from database. resume_from_page is optional
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
            'max_threads' : 25, # controls download speed.  Make it 10 if what2retrieve=ALL_PROPERTIES
            'read_raw' : False
        }

        ent_props = list(kwargs.pop('ent_props',[]))
        rel_props = list(kwargs.pop('rel_props',[]))

        my_kwargs = dict(default_kwargs)
        my_kwargs.update(kwargs)

        self.cache_was_modified = False

        super().__init__(*args, **my_kwargs)
        # additional properties from kwargs have to be added after super().__init__ adds properties from what2retrieve
        super().add_ent_props(ent_props)
        super().add_rel_props(rel_props)
        
        my_kwargs['ent_props'] = ent_props
        my_kwargs['rel_props'] = rel_props # to pass for printing
        
        self.relprops2rnef = list(rel_props)
        self.relprops2rnef += [p for p in [EFFECT,'URN',REFCOUNT] if p not in rel_props]
        
        # loads cache_name from self.data_dir:
        load_cache = my_kwargs.pop('load_cache',True)
        if load_cache:
            self.network = self._load_cache(**my_kwargs)
            if self.network:
                self.save_description()
  

    def __path2cache(self,cache_name:str,extension='.rnef'):
        return self.data_dir+'simple/'+cache_name+extension
    

    @staticmethod
    def cache_filename(parameters:dict):
        return parameters['cache_name']+'.rnef'


    def simplify_graph(self,database_g:ResnetGraph,**kwargs)->'ResnetGraph':
        # filter to load relations only for nodes with desired properties\n
        rank4simplifying = kwargs.pop('ranks4simplifying',[])
        predict_effect4 = kwargs.pop('predict_effect4',dict()) #kwargs for ResnetGraph.predict_effect4(): _4enttypes, _4reltypes
        refprop2rel = dict(kwargs.pop('refprop2rel',dict()))
        refprop_minmax = kwargs.pop('refprop_minmax',0)

        simple_graph = database_g.copy()

        if kwargs['make_simple']:
            simple_graph = simple_graph.make_simple(rank4simplifying)

        #converting sentence properties to relation properties before making a dump.  Used for pX dump
        for refprop, relprop in refprop2rel.items():
            simple_graph.refprop2rel(refprop,relprop,refprop_minmax)
            self.relprops2rnef.remove(refprop)
            self.relprops2rnef.append(relprop)
            
        if predict_effect4:
            simple_graph,modified_rels = simple_graph.predict_effect4(**predict_effect4)
            modified_rels_graph = simple_graph.subgraph_by_rels(modified_rels)
            modified_rels_graph.name = 'relations with predicted effect'
            modified_rels_graph.dump2rnef('Effect predicted4',self.entProps,self.relprops2rnef)
        
        if kwargs['no_id_version']:
            for ent_prop in self.entProps:
                if ent_prop[-3:] == ' ID':
                    simple_graph = simple_graph.clean_version_number(ent_prop)

        return simple_graph
    

    def save_network(self, example_node_name=''):
        if example_node_name:
            self.example_node = self.network._psobjs_with(example_node_name)[0]
        if self.cache_was_modified:
            self.replace_cache(self.network.name,self.network,self.entProps,self.relprops2rnef)
            self.save_description()
        return
    

    def save_description(self,do_backup=True):
        stats_df = self.network.get_stats()
        descr_file = self.data_dir+'descriptions/'+self.network.name+"_description"
        if do_backup:
            backup_descr = descr_file + '_old.tsv'
            try:
                copyfile(descr_file, backup_descr)
            except FileNotFoundError:
                pass
        stats_df.to_csv(descr_file+'.tsv',sep='\t',index=False)


    def __dump_dir_name(self,**kwargs):
        self.dump_folder = kwargs['cache_name'] + '_raw'
        return self.data_dir+'raw/'+self.dump_folder+'/'


    def __download(self,**kwargs):
        cache_name = kwargs['cache_name']
        dump_dir = self.__dump_dir_name(**kwargs)
        print(f'Downloading {cache_name} into {dump_dir}')
        my_session_kwargs = {
            'connect2server':True,
            'no_mess' : self.no_mess,
            'data_dir' : kwargs['data_dir']
            }
        my_session = self._clone_session(**my_session_kwargs) #need identifiers to make graph simple

        log_name = f'Download of {cache_name}.log'
        log_path = self.data_dir+'log/'+log_name
        with open(log_path, 'w', buffering=1) as fout:
            print(f'Download log is in {log_path}')
            with redirect_stdout(fout):
                if kwargs['connect_nodes']:
                    node_ids = set(kwargs['connect_nodes'][0])
                    # node_ids can be either {str} or {int}
                    if isinstance(next(iter(node_ids), None),str):
                        node_types_str = ",".join(node_ids)
                        rel_types = kwargs['connect_nodes'][1]
                        rel_types_str = ",".join(rel_types)
                        query = f'SELECT Entity WHERE objectType = ({node_types_str}) AND Connectivity > 0 AND NeighborOf (SELECT Relation WHERE objectType = ({rel_types_str}))'
                        nodes_graph = my_session.process_oql(query,'Loading nodes from database')
                        assert(isinstance(nodes_graph,ResnetGraph))                  
                        node_ids = set(ResnetGraph.dbids(nodes_graph.psobjs_with())) if ResnetGraph else set()                  
                        print(f'Found {len(node_ids)} entities of type: {node_types_str} connected with relations {rel_types_str} exist in database')
                    if node_ids:
                        print(f'Will connect {len(node_ids)} with {rel_types} relations')
                        my_session.get_network(node_ids,rel_types,download=True)
                else:
                    for i, query in enumerate(self.my_oql_queries):
                        oql = query[0]
                        reqname = query[1]
                        resume_from = query[2] if len(query)>2 else 0
                        my_session.download_oql(oql,reqname,resume_page=resume_from,threads=kwargs['max_threads'])


    def __read_raw(self,**kwargs):
        dump_dir = self.__dump_dir_name(**kwargs)
        database_graph = ResnetGraph.fromRNEFdir(dump_dir,merge=False)
        if not database_graph and kwargs['connect2server']:
            print(f'Cache "{dump_dir}" was not found\nBegin graph download from database')
            self.__download(**kwargs)
            database_graph = ResnetGraph.fromRNEFdir(dump_dir,merge=False)

        return database_graph


    def _load_cache(self,**kwargs)->'ResnetGraph':
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
        cache_name = kwargs.get('cache_name',DEFAULT_CACHE_NAME)
        print(f'Loading {cache_name} cache')
        prop2values = kwargs.pop('prop2values',dict()) #prop2values={prop_name:[values]}
        # filter to load relations only for nodes with desired properties\n

        if kwargs.get('read_raw',False):
            database_graph = self.__read_raw(**kwargs)
            database_graph.name = cache_name
            return self.simplify_graph(database_graph,**kwargs)
        else:
            my_cache_file = self.__path2cache(cache_name)
            try:
                cached_graph = ResnetGraph.fromRNEF(my_cache_file,prop2values=prop2values)
                print(f'Loaded "{cache_name}" cache with {len(cached_graph)} nodes and {cached_graph.number_of_edges()} edges')
                cached_graph.name = cache_name
                return cached_graph
            except FileNotFoundError:
                print(f'Cannot find {my_cache_file} cache file')
                self.clear() #clearing self.Graph and self.network to releasae RAM

                if kwargs.get('resume_download',False):
                    dump_dir = self.__dump_dir_name(**kwargs)
                    print(f'Resume graph download from database into {dump_dir}')
                    self.__download(**kwargs)

                database_graph = self.__read_raw(**kwargs)
                
                if database_graph:
                    database_graph.name = cache_name
                    database_graph = self.simplify_graph(database_graph,**kwargs)
                    database_graph_nodup = database_graph.remove_undirected_duplicates()
                    database_graph_nodup.dump2rnef(my_cache_file,self.entProps,self.relprops2rnef,with_section_size=1000) #with_section_size=1000
                    
                    print('%s with %d edges and %d nodes was written into %s' 
                        % (database_graph.name,database_graph.number_of_edges(),database_graph.number_of_nodes(),my_cache_file))
                else:
                    print(f'No data is available in database for {cache_name}')
                return database_graph


    def replace_cache(self,cache_name:str,with_network:ResnetGraph,
                      ent_props:list[str]=['Name'],rel_props:list[str]=['URN'],
                      do_backup=True,replace_raw=False):
        if replace_raw:
            my_cache_dir = self.data_dir+cache_name
            listing = glob.glob(my_cache_dir+'*.rnef')
            [os.remove(file) for file in listing]
            self._dump2rnef(with_network,cache_name)
            print(f'{cache_name} raw cache files were replaced')

        my_cache_file = self.__path2cache(cache_name)
        if do_backup:
            backup_cache = self.__path2cache(cache_name,'_backup.zip')
            with zipfile.ZipFile(backup_cache,'w',zipfile.ZIP_DEFLATED) as z:
                z.write(my_cache_file)
        
        with_network.dump2rnef(my_cache_file,ent_props,rel_props)
        print(f'{cache_name} cache file was replaced')

        if self.example_node:
            example_cache = with_network.neighborhood({self.example_node})
            example_file = self.__path2cache(cache_name,extension='_example.rnef')
            example_cache.name = self.example_node.name()+' example'
            example_cache.dump2rnef(example_file,ent_props,rel_props)
            print(f'{example_file} example file was replaced')

        
    def clear(self):
        super().clear()
        self.network.clear_resnetgraph()


    def dump_subgraph_by_relprops(self,_2file:str,search_values:list[str|int],
                                  in_properties:list[str]=['ObjTypeName']):
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
            self.network = self.simplify_graph(combined_network,**kwargs)
            self.network.name = cache_name
            self.cache_was_modified = True
            self.save_network() #with_section_size=1000

    @staticmethod
    def list_cache(data_dir:str):
        exist_cache = glob.glob(f"{data_dir}*.rnef")
        exist_cache = [Path(f).stem for f in exist_cache]
        print(f'dump directory contains {len(exist_cache)} cache files.')
        return exist_cache
  

    def subgraph(self,with_objs:list[PSObject],new_cache_name:str):
        subg = self.network.subgraph(ResnetGraph.uids(with_objs))
        if subg:
            subg.name = new_cache_name
            subg = subg.remove_undirected_duplicates()
            subg.dump2rnef(self.data_dir+new_cache_name,self.entProps,self.relprops2rnef)
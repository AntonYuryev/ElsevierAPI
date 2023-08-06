import  os, glob, json
from shutil import copyfile
from .ResnetAPISession import TO_RETRIEVE,NO_REL_PROPERTIES,DATABASE_REFCOUNT_ONLY,REFERENCE_IDENTIFIERS,BIBLIO_PROPERTIES,SNIPPET_PROPERTIES,ONLY_REL_PROPERTIES,ALL_PROPERTIES
from .ResnetAPISession import APISession,ResnetGraph,PSObject


class APIcache(APISession):
    '''
    manages work with cache files
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
        APIconfig = args[0]\nwhat2retrieve - defaults NO_REL_PROPERTIES, other options -\n
        [DATABASE_REFCOUNT_ONLY,REFERENCE_IDENTIFIERS,BIBLIO_PROPERTIES,SNIPPET_PROPERTIES,ONLY_REL_PROPERTIES,ALL_PROPERTIES]
        no_mess - default True, if False your script becomes more verbose
        connect2server - default True, set to False to run script using data in __pscache__ files instead of database
        '''
        my_kwargs = dict(kwargs)
        super().__init__(*args, **my_kwargs)
        try:
            self.cache_dir = my_kwargs['cache_dir']
            if self.cache_dir[-1] != '/':
                self.cache_dir = self.cache_dir + '/' 
        except KeyError:
            self.cache_dir = 'ElsevierAPI/ResnetAPI/__pscache__/'


    def __path2cache(self,cache_name:str,extension='.rnef'):
        return self.cache_dir+cache_name+extension


    def load_cache(self,cache_name, **kwargs)->'ResnetGraph':
        """
        Input
        -----

        oql_queries - list of tuples (oql_query,request_name)\n
        prop2values={prop_name:[values]} - filter to load relations only for nodes with desired properties\n
        predict_effect4 - kwargs for ResnetGraph.predict_effect4()
        refprop2rel - {refPropName:NewRelPropName}
        "rank4simplifying" is used to merge relations during cache network simplification\n
        relation types must be ordered from most important type to

        if refprop_minmax < 0 assigns single value to relprop that is the minimum of all ref_prop values
        if refprop_minmax > 0 assigns single value to relprop that is the maximum of all ref_prop values

        Examples:\n
        ['DirectRegulation','Binding','ProtModification','Regulation']\n
        ['PromoterBinding','Expression','Regulation']\n
        ['Regulation','Biomarker','StateChange','QuantitativeChange']\n
        ['Biomarker','StateChange']\n
        ['Biomarker','QuantitativeChange']\n
        [MolSynthesis','Regulation']\n
        [MolTransport','Regulation']\n
        [MolTransport','CellExpression']\n
        ['Regulation','FunctionalAssosiation']\n
        """
        oql_queries = kwargs.pop('oql_queries',[])
        ent_props = kwargs.pop('ent_props',['Name'])
        rel_props = kwargs.pop('rel_props',['URN'])
        prop2values = kwargs.pop('prop2values',dict())
        rank4simplifying = kwargs.pop('rank4simplifying',[])
        predict_effect4 = kwargs.pop('predict_effect4',dict()) #kwargs for ResnetGraph.predict_effect4(): _4enttypes, _4reltypes
        refprop2rel = kwargs.pop('refprop2rel',dict())
        refprop_minmax = kwargs.pop('refprop_minmax',0)

        my_cache_file = self.__path2cache(cache_name)
        try:
            cached_graph = ResnetGraph.fromRNEF(my_cache_file,prop2values=prop2values)
            return cached_graph
        except FileNotFoundError:
            dump_folder = cache_name + '_raw'
            database_graph = ResnetGraph.fromRNEFdir(self.cache_dir+dump_folder+'/',merge=False)
            if not database_graph:
                print(f'Cache "{self.dump_folder}" was not found\nBegin caching network from database')
                my_session_kwargs = {
                    TO_RETRIEVE:REFERENCE_IDENTIFIERS, 
                    'connect2server':True,
                    'no_mess' : self.no_mess}
                
                my_session = self._clone_session(**my_session_kwargs) #need identifiers to make graph simple
                my_session.dump_folder = dump_folder
                my_session.add_ent_props(ent_props)
                my_session.add_rel_props(rel_props)
                my_session.set_dir(self.cache_dir)
                can_close = False
                for i, oql_query in enumerate(oql_queries):
                    query = oql_query[0]
                    req_name = oql_query[1]
                    if i == len(oql_query) - 1:
                        can_close = True
                    my_session.download_oql(query,req_name,can_close=can_close)
                my_session.clear()
                database_graph = ResnetGraph.fromRNEFdir(self.cache_dir+dump_folder+'/',merge=False)
            
            database_graph.name = cache_name
            simplified_network = database_graph.make_simple(rank4simplifying)
            
            rnefrelporps = list(rel_props)
            for refprop, relprop in refprop2rel.items():
                simplified_network.refprop2rel(refprop,relprop,refprop_minmax)
                rnefrelporps.remove(refprop)
                rnefrelporps.append(relprop)
                
            if predict_effect4:
                simplified_network,modified_rels = simplified_network.predict_effect(**predict_effect4)
                modified_rels_graph = simplified_network.subgraph_by_rels(modified_rels)
                modified_rels_graph.name = 'relations with predicted effect'
                modified_rels_graph.dump2rnef('Effect predicted4',ent_props,rnefrelporps)
                    
            for ent_prop in ent_props:
                if ent_prop[-3:] == ' ID':
                    simplified_network = simplified_network.clean_version_number(ent_prop)
            
            simplified_network.dump2rnef(my_cache_file,ent_props,rnefrelporps) #with_section_size=1000
            network_description = {'Name':[cache_name],'OQLs':oql_queries}
            with open(self.cache_dir+cache_name+"_description.json", "w") as f:
                json.dump(network_description, f,indent=2)
            
            print('%s with %d edges and %d nodes was written into %s' 
                % (cache_name,simplified_network.number_of_edges(),simplified_network.number_of_nodes(),my_cache_file))
            
            return simplified_network


    def replace_cache(self,cache_name:str,with_network:ResnetGraph,
                      ent_props=['Name'],rel_props=['URN'],
                      do_backup=True,replace_raw=False):
        if replace_raw:
            my_cache_dir = self.cache_dir+cache_name
            listing = glob.glob(my_cache_dir+'*.rnef')
            [os.remove(file) for file in listing]
            self.__dump2rnef(with_network,cache_name)

        my_cache_file = self.__path2cache(cache_name)
        if do_backup:
            backup_cache = self.__path2cache(cache_name,'_backup.rnef')
            copyfile(my_cache_file, backup_cache)
        
        with_network.dump2rnef(my_cache_file,ent_props,rel_props)
        if self.example_node:
            example_cache = with_network.neighborhood([self.example_node])
            example_file = self.__path2cache(cache_name,extension='_example.rnef')
            example_cache.name = self.example_node.name()+' example'
            example_cache.dump2rnef(example_file,ent_props,rel_props)

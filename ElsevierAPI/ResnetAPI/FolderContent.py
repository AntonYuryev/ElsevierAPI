from .ResnetAPISession import APISession, time, len,SNIPPET_PROPERTIES
from .ResnetGraph import ResnetGraph,RESNET
import glob
from .NetworkxObjects import PSObject
from .PSPathway import PSPathway
import xml.etree.ElementTree as et
from urllib.parse import quote
from concurrent.futures import ThreadPoolExecutor,as_completed

PATHWAY_ID = 'Pathway ID'
FOLDER_OBJECTS = ['Pathway','Group','attributesSearch']
FOBJECTS_PROPS = ['Name','Description','Notes']
MAX_SESSIONS = 20

class FolderContent(APISession): 
    pass
    def __init__(self,  *args,**kwargs):
        my_kwargs = {'preload_folder_tree':True,'what2retrieve':SNIPPET_PROPERTIES}
        my_kwargs.update(kwargs)
        super().__init__(*args,**my_kwargs)
        self.PageSize = 1000
        self.DumpFiles = []
        self.id2folder = dict()
        self.dbid2folder_obj = dict() # {dbid:PSObject} 
        # where PSObject corresponds to types in FOLDER_OBJECTS

        if my_kwargs['preload_folder_tree']:
            self.id2folder = self.load_folder_tree() # {folder_id:folder_zobj, folder_name:folder_zobj}


    def clone(self,**kwargs):
        my_kwargs = dict(kwargs)
        copy_graph = my_kwargs.pop('copy_graph',False)
        api_session = super()._clone_session(**my_kwargs)
        my_kwargs['preload_folder_tree'] = False
        new_session = FolderContent(api_session.APIconfig,**my_kwargs)
        if copy_graph:
            new_session.Graph = self.Graph
        new_session.entProps = api_session.entProps
        new_session.relProps = api_session.relProps
        new_session.data_dir = api_session.data_dir
        new_session.id2folder = dict(self.id2folder)
        return new_session


    def __folder_id(self,folder_name:str) -> str:
        if not hasattr(self, "id2folder"): 
            self.id2folder = self.load_folder_tree()

        for k,v in self.id2folder.items():
            if v['Name'] == folder_name:
                return str(k)
        return ''
    

    def __folder_name(self,folder_dbid:int):
        return  self.id2folder[int(folder_dbid)]['Name']
    

    @staticmethod
    def __urn4result(o:PSObject):
        return 'urn:agi-result:'+quote(o['Name'][0])


    @staticmethod
    def __objtype4result(result_graph:ResnetGraph):
        if result_graph.number_of_edges():
            return 'Pathway'
        else:
            return 'Group'


    def load_subfolders(self, subfolder_ids:list):
        # FolderIds - list of folder IDs
        if not hasattr(self, "id2folder"): 
            self.id2folder = self.load_folder_tree()

        subfolders_ids = set()
        for folder_id in subfolder_ids:
            folder = self.id2folder[int(folder_id)]
            if type(folder['SubFolders']) != type(None):
                subfolders_ids.update(folder['SubFolders']['long'])

        return subfolders_ids


    def subfolder_ids(self, folder_name:str, include_parent=True):
        """
        Returns {subfolders_ids+parent_folder_id}
        -------  
        """
        FolderId = self.__folder_id(folder_name)
        accumulate_subfolder_ids = {FolderId} if include_parent else set()
        subfolder_ids = set(self.load_subfolders([FolderId]))
        while len(subfolder_ids) > 0:
            accumulate_subfolder_ids.update(subfolder_ids)
            subs = set()
            for db_id in subfolder_ids:
                subs.update(self.load_subfolders([db_id]))
            subfolder_ids = subs

        return accumulate_subfolder_ids


    def get_subfolder_tree(self,folder_name):
        """
        Return
        ------
        child2parent = {folder_id:folder_id}
        parent2child = {parent_id:[child_ids]}
        """
        if not hasattr(self, "id2folder"): 
            self.id2folder = self.load_folder_tree()

        folder_id = self.__folder_id(folder_name)
        child2parent = {folder_id:folder_id}
        parent2child = PSObject(dict())

        if type (self.id2folder[folder_id]['SubFolders']) != type(None):
            children_ids = self.id2folder[folder_id]['SubFolders']['long']
            p2c = {folder_id : children_ids} 
        else: 
            return {folder_id:folder_id}, dict()

        while len(p2c) > 0:
            parent2child.update(p2c)
            p2c_iter = dict()
            for par, chi in p2c.items():
                for c in chi:
                    child2parent[c] = par
                    if type (self.id2folder[c]['SubFolders']) != type(None):
                        c_subs = self.id2folder[c]['SubFolders']['long']
                        p2c_iter[c] = c_subs
                    else: continue
            p2c = p2c_iter

        return child2parent, dict(parent2child)


    def get_objects_from_folders(self, FolderIds:list, with_layout=False):
        """
        Updates
        -----
        self.dbid2folder_obj = {id:PSObject}
        """
        for fid in FolderIds:
            folder_name = self.__folder_name(fid)
            zeep_objects = self.get_folder_objects_props(fid, FOBJECTS_PROPS)
            self.dbid2folder_obj = self._zeep2psobj(zeep_objects)
            # annotating all folder object with "Folders" property:
            [o.update_with_value('Folders', folder_name) for o in self.dbid2folder_obj.values()]
        
        if with_layout:
            [o.update_with_value('layout', self.get_layout(o.dbid())) for o in self.dbid2folder_obj.values() if o.objtype() == 'Pathway']
   
        return


    @staticmethod
    def _put2folder(named:str, folder_obj:PSObject,resnet:et.Element):
        xml_nodes = et.SubElement(resnet, 'nodes')
        folder_local_id = 'F0'
        xml_node_folder = et.SubElement(xml_nodes, 'node', {'local_id':folder_local_id, 'urn': 'urn:agi-folder:xxxxx_yyyyy_zzzzz'})
        et.SubElement(xml_node_folder, 'attr', {'name': 'NodeType', 'value': 'Folder'})
        et.SubElement(xml_node_folder, 'attr', {'name': 'Name', 'value': named})
        pathway_local_id = 'P0'
        xml_node_pathway = et.SubElement(xml_nodes, 'node', {'local_id':pathway_local_id, 'urn':folder_obj.urn()})
        et.SubElement(xml_node_pathway, 'attr', {'name': 'NodeType', 'value': 'Pathway'})
        xml_controls = et.SubElement(resnet, 'controls')
        xml_control = et.SubElement(xml_controls, 'control', {'local_id':'CFE1'})
        et.SubElement(xml_control, 'attr', {'name':'ControlType', 'value':'MemberOf'})
        et.SubElement(xml_control, 'link', {'type':'in', 'ref':pathway_local_id})
        et.SubElement(xml_control, 'link', {'type':'out', 'ref':folder_local_id})
 

    def load_folder_obj(self,fobj:PSObject,props4rel=dict(),props4pathway=dict(),prettify=True):
        '''
        Return
        ------
        tuple (fobj:PSObject,fobj_graph:ResnetGraph,fobj_xml:str)\n
        if fobj.objtype() == 'attributesSearch' adds URN and ObjTypeName
        '''
        if fobj.objtype() not in FOLDER_OBJECTS:
            return fobj, ResnetGraph(),str()
        
        my_fobj = PSObject(fobj)
        obj_graph = self.graph4obj(my_fobj)
        if my_fobj.objtype() == 'attributesSearch':
            my_fobj.set_property('URN',self.__urn4result(my_fobj))
            my_fobj.set_property('ObjTypeName',self.__objtype4result(obj_graph))

        fobj_props = dict(props4pathway)
        descr = my_fobj.descr()
        if descr:
            fobj_props.update({'Description':[descr]})
            
        notes = my_fobj.notes()
        if notes:
            fobj_props.update({'Notes':[notes]})

        rnef = et.fromstring(self._2rnefs(obj_graph,props4rel,fobj_props))
        rnef.set('name', my_fobj.name())
        rnef.set('urn', my_fobj.urn())
        rnef.set('type', my_fobj.objtype())

        if fobj.objtype() == 'Pathway': 
        # here original fobj must be 'Pathway'. attributesSearch do not have layout
            lay_out = et.Element('attachments')
            attachments = self.get_layout(fobj.dbid())
            if attachments:
                lay_out.append(et.fromstring(attachments))
                rnef.append(lay_out)
          
        xml_str = str(et.tostring(rnef,encoding='utf-8').decode("utf-8"))
        # xml_str here must have xml_declaration in order for pretty_xml to work
        if prettify:
            xml_str = self.pretty_xml(xml_str,remove_declaration=True)
 
        print(f'Downloaded "{fobj.name()}" {str(fobj.objtype()).lower()} with {len(obj_graph)} nodes and {obj_graph.number_of_edges()} relations')
        return my_fobj,obj_graph,xml_str


    def folder_content(self,folder_id_or_name:int or str,parent_folder_name:str,folderobjs2skip:set=None,
                            props4rel=dict(),props4pathway=dict()):
        """
        Input
        -----
        either folder_id or folder_name must be supplied
        if folder_id is supplied folder_name is retreived from database
        folderobjs2skip - {PSObject}

        Updates
        -----
        self.dbid2folder_obj = {id:PSObject}

        Dumps
        -----
        Objects from "folder_id_or_name" into 'folder_name' inside 'parent_folder_name' located in "self.data_dir"
        if size of dump file exceeds 100000000, "rnef_xml" is splitted into several RNEF files\n
        dump RNEF files are named as: 'content of folder_name#',
        where # - dump file number

        Returns
        -------
        tuple (new_pathway_counter:int,symlinks:{PSObject},folder_name:str,parent_folder_name:str)
        """

        folder_id = self.__folder_id(folder_id_or_name) if isinstance(folder_id_or_name, str) else folder_id_or_name
        folder_name = str(self.id2folder[folder_id]['Name'])
        self.get_objects_from_folders([folder_id],with_layout=True)
        if self.dbid2folder_obj:
            print('Start downloading %d pathways from \"%s\" folder' % (len(self.dbid2folder_obj),folder_name))
        else:
            print('Folder \"%s\" has no pathways' % folder_name)
            return int(0),set(),folder_name,parent_folder_name
        
        new_pathway_counter = 0
        folder_object_counter = 0
        symlinks = set()
        folder_download_start = time.time()
        write2folder = '' if parent_folder_name == folder_name else parent_folder_name

        thread_name = f'{folder_name} download'
        futures = list()
        with ThreadPoolExecutor(max_workers=1, thread_name_prefix=thread_name) as e:
            start_time = time.time()
            for dbid, folder_obj in self.dbid2folder_obj.items():
                if folder_obj['IsSymlink'][0]:
                    symlinks.add(folder_obj)
                else:
                    #printing pathways only symlinks should be printed at the end and only if they were not downloaded for another folder
                    if isinstance(folderobjs2skip,set) and folder_obj in folderobjs2skip:
                        print(f'{folder_obj.name()} was downloaded for another folder')
                        continue

                    if folder_obj.objtype() in FOLDER_OBJECTS:
                        futures.append(e.submit(self.load_folder_obj,folder_obj,props4rel,props4pathway,True))
                        new_pathway_counter += 1
                    else:
                        print ('%s folder has object with unknown type %s: id = %d' % (folder_name,folder_obj['ObjTypeName'][0],dbid))
                        continue

            for f in as_completed(futures):
                folder_obj, obj_graph, obj_xml = f.result()
                self.rnefs2dump(obj_xml,folder_name,write2folder,can_close=False)
                folder_object_counter += 1
                
                if isinstance(folderobjs2skip,set):
                    folderobjs2skip.add(folder_obj)

                #printing folder object membership
                folder_local_id = 'F0'
                folder_resnet = et.Element(RESNET)
                folder_nodes = et.SubElement(folder_resnet, 'nodes')
                member_controls = et.SubElement(folder_resnet, 'controls')

                xml_node_folder = et.SubElement(folder_nodes, 'node', {'local_id':folder_local_id, 'urn':'urn:agi-folder:'+str(folder_id)})
                et.SubElement(xml_node_folder, 'attr', {'name':'NodeType', 'value':'Folder'})
                et.SubElement(xml_node_folder, 'attr', {'name':'Name', 'value':folder_name})

                pathway_local_id = 'P'+str(folder_object_counter)
                folder_pathway_node = et.SubElement(folder_nodes, 'node', {'local_id':pathway_local_id, 'urn':folder_obj.urn()})
                et.SubElement(folder_pathway_node, 'attr', {'name': 'NodeType', 'value': 'Pathway'})
                et.SubElement(folder_pathway_node, 'attr', {'name': 'Name', 'value': folder_obj.name()})
                
                control_local_id = 'L'+str(folder_object_counter)
                member_control = et.SubElement(member_controls, 'control', {'local_id':control_local_id})
                et.SubElement(member_control, 'attr', {'name':'ControlType', 'value':'MemberOf'})
                if folder_obj['IsSymlink'][0] == True:
                    et.SubElement(member_control, 'attr', {'name':'Relationship', 'value':'symlink'})
                et.SubElement(member_control, 'link', {'type':'in', 'ref':pathway_local_id})
                et.SubElement(member_control, 'link', {'type':'out', 'ref':folder_local_id})

                folder_rnef = et.tostring(folder_resnet, encoding='utf-8',xml_declaration=False).decode("utf-8")
                folder_rnef = self.pretty_xml(folder_rnef,remove_declaration=True)
                self.rnefs2dump(folder_rnef,folder_name,write2folder)

        print('%d out of %d objects in folder "%s" was downloaded in %s' %
                    (folder_object_counter, len(self.dbid2folder_obj),folder_name, self.execution_time(start_time)))
        print('Total folder download time: %s' % self.execution_time(folder_download_start))
        return new_pathway_counter,symlinks,folder_name,parent_folder_name


    def __make_cache_dir(self, root_folder_id):
        """
        Creates
        -------
        directory structure in self.data_dir mirroring foler structure in database
        """
        subfolder_ids = self.load_subfolders([root_folder_id])
        subfolders_dict = {root_folder_id:subfolder_ids} if subfolder_ids else dict()

        while subfolders_dict:
            subs = dict()
            for parent_id, subfolder_ids in subfolders_dict.items():
                parent_folder_name = self.id2folder[parent_id]['Name']
                self.dump_path(parent_folder_name)
                for subfolder_id in subfolder_ids:
                    subsub_folder_ids = self.load_subfolders([subfolder_id])
                    if subsub_folder_ids:
                        subs.update({subfolder_id : self.load_subfolders([subfolder_id])})

                    subfolder_name = self.id2folder[subfolder_id]['Name']
                    self.dump_path(subfolder_name,parent_folder_name)
            subfolders_dict = subs


    def __foldtree2rnef(self, root_folder_name:str):
        """
        Loads
        -------
        self.subfold2fold = {folder_id:folder_id}

        Writes
        ------
        folder tree in RNEF format to self.data_dir/root_folder_name

        Return
        ------
        [Folder Objects]
        
        """
        self.subfold2fold_ids, parent2child = self.get_subfolder_tree(root_folder_name)
        subtree_xml = str()

        for parent_id, child_ids in parent2child.items():
            folder_resnet = et.Element(RESNET)
            folder_name = self.id2folder[parent_id]['Name']
            folder_nodes = et.SubElement(folder_resnet, 'nodes')
            xml_controls = et.SubElement(folder_resnet, 'controls')
            parent_local_id = 'F0'
            xml_node_folder = et.SubElement(folder_nodes, 'node', {'local_id':parent_local_id, 'urn':'urn:agi-folder:'+str(parent_id)})
            et.SubElement(xml_node_folder, 'attr', {'name':'NodeType', 'value':'Folder'})
            et.SubElement(xml_node_folder, 'attr', {'name':'Name', 'value':folder_name})
            for child_id in child_ids:
                child_local_id = str(child_id)
                subfolder_name = self.id2folder[child_id]['Name']
                xml_node_pathway = et.SubElement(folder_nodes, 'node', {'local_id':child_local_id, 'urn':'urn:agi-folder:'+child_local_id})
                et.SubElement(xml_node_pathway, 'attr', {'name':'NodeType', 'value':'Folder'})
                et.SubElement(xml_node_pathway, 'attr', {'name':'Name', 'value':subfolder_name})
                
                control_local_id = 'L'+str(child_id)
                xml_control = et.SubElement(xml_controls, 'control', {'local_id':control_local_id})
                et.SubElement(xml_control, 'attr', {'name':'ControlType', 'value':'MemberOf'})
                et.SubElement(xml_control, 'link', {'type':'in', 'ref':child_local_id})
                et.SubElement(xml_control, 'link', {'type':'out', 'ref':parent_local_id})

            unpretty_xml = et.tostring(folder_resnet,encoding='utf-8',xml_declaration=False).decode("utf-8")
            subtree_xml += self.pretty_xml(unpretty_xml,remove_declaration=True)

        self.__make_cache_dir(self.__folder_id(root_folder_name))
        self.rnefs2dump(subtree_xml,root_folder_name)
        return [folder for dbid,folder in self.id2folder.items() if dbid in self.subfold2fold_ids]


    def __download_folders(self,folders:list,job_name:str,printed_pathways=set(),props4rel=dict(),props4pathway=dict()):
        '''
        Input
        -----
        folder_objs = {PSObject}
        '''
        chunks = list()
        symlinks = set()
        start_time = time.time()
        for start in range(0,len(folders),MAX_SESSIONS):
            end = start+MAX_SESSIONS
            chunks.append(folders[start:end])

        for i,chunk in enumerate(chunks):
            thread_name = f'{job_name} download {i+1}of{len(chunks)}'
            print(f'Downloading {len(folders)} folders in {MAX_SESSIONS} threads')
            
            with ThreadPoolExecutor(max_workers=MAX_SESSIONS, thread_name_prefix=thread_name) as e:
                futures = list()
                sessions = list()   
                for folder in chunk:
                    folder_id = folder['Id']
                    parent_id = self.subfold2fold_ids[folder_id]
                    parent_name = self.id2folder[parent_id]['Name']
                    new_session = self.clone(copy_graph=True)
                    futures.append(e.submit(new_session.folder_content,folder_id,parent_name,printed_pathways,props4rel,props4pathway))
                    sessions.append(new_session)

                download_counter = 0
                for f in as_completed(futures):
                    pathway_counter,syms,subfolder_name,parent_folder_name = f.result()
                    self.close_rnef_dump(subfolder_name,parent_folder_name)
                    download_counter += pathway_counter
                    symlinks.update(syms)

                for s in sessions:
                    self.Graph = self.Graph.compose(s.Graph) # accumulate s.Graphs into cache
                    s.close_connection()

            if self.Graph.weight() > self.reference_cache_size:
                print('Clearing cache due to large size: %d' % self.Graph.weight())
                self.clear()
            
        print(f'{job_name} downloaded {download_counter} objects from {len(folders)} out of {len(self.id2folder)} folders in {self.execution_time(start_time)}')
        print('Graph cache has %d references\n' % self.Graph.weight())
        return symlinks.difference(printed_pathways)


    def folder2rnef(self,root_folder_name:str,include_subfolders=True,add_props2rel=dict(),add_pathway_props=dict()):
        """
        Dumps
        -----
        content of root_folder_name into RNEF file 'content of root_folder_name.rnef' located in 'self.data_dir'
        """
        download_start_time = time.time()
        if not include_subfolders:
            return self.folder_content(None,root_folder_name,add_props2rel=add_props2rel,add_pathway_props=add_pathway_props)
        else:
            subfolders = self.__foldtree2rnef(root_folder_name)
            thread_name = f'{root_folder_name} download'
            symlinks2print = self.__download_folders(subfolders,thread_name)
        
            if symlinks2print:
                print(f'Will print {len(symlinks2print)} pathways to support symlinks')
                # dumping loose symlinks into root folder
                for folder_obj in symlinks2print:
                    fobj, pathway_graph, pathway_xml = self.load_folder_obj(folder_obj)
                    self.rnefs2dump(pathway_xml,root_folder_name)
                self.close_rnef_dump(root_folder_name)
            else:
                print(f'No additional symlinks are required for {root_folder_name} folder %s')
                
        print("Complete download execution time: %s" % self.execution_time(download_start_time))


    @staticmethod
    def __downloaded(in_rnef_file:str):
        folder_objs = set()
        with open(in_rnef_file, "r", encoding='utf-8') as f:
            line = f.readline()
            while line:
                line = line.strip()
                if line[:8] == '<resnet ':
                    urn_pos = line.find('urn=',8)
                    if urn_pos > 0 : 
                        urn_start = urn_pos+5
                        urn_end = line.find('\"', urn_start)
                        urn = line[urn_start:urn_end]
                        o = PSObject({'URN':[urn]})
                        folder_objs.add(0)
                line = f.readline()
        print(f'Read \"{in_rnef_file}\" with {len(folder_objs)} pathways')
        return folder_objs


    def resume_download(self, root_folder_name:str, last_downloaded_folder:str):
        child2parent = self.get_subfolder_tree(root_folder_name)[0]
        subfolder_ids = [list(child2parent.keys())]
        global_start = time.time()
        try:
            last_downloaded_folder_id = self.__folder_id(last_downloaded_folder)
            start_folder_idx = subfolder_ids.index(last_downloaded_folder_id)
            folder_counter = start_folder_idx+1
            continue_from_dir =  self.data_dir+self.filename4(root_folder_name)
            listing = glob.glob(continue_from_dir+"/**/*.rnef", recursive=True)
            downloaded_pathways = set()
            [downloaded_pathways.update(self.__downloaded(rnef_file)) for rnef_file in listing]
            download_counter = len(downloaded_pathways)

            last_rnef = listing[-1]
            print (f'Resuming download of {root_folder_name} folder from {last_rnef}')
            folder_objs = [folder for dbid,folder in self.id2folder.items() if dbid in self.subfold2fold_ids]
            folder_objs2download = folder_objs[start_folder_idx+1:]
            job_name = f'Resume {root_folder_name} download'
            symlinks2print = self.__download_folders(folder_objs2download,job_name,downloaded_pathways)
           
            if len(symlinks2print) > 0:
                print(f'Will print {len(symlinks2print)} pathways to support symlinks')
                # dumping loose symlinks into root folder
                for symlink in symlinks2print:
                    xml_str = self.load_folder_obj(symlink)[2]
                    self.rnefs2dump(xml_str,root_folder_name)
                self.close_rnef_dump(root_folder_name)
            else:
                print('All pathways for symlinks were printed for other folders')
                print(f'Resumed download was finished in {self.execution_time(global_start)}')
        except ValueError:
            print ('Folder %s was not found in folder tree' % last_downloaded_folder)
            return
    

    def load_containers(self,with_props=['Name','Description','Notes'],from_folders=[]):
        """
        retreives entire folder tree stricture from database 
        including all folder objects
        if from_folders is not specidied works for ~40sec
        Updates
        -----
        self.dbid2folder_obj = {id:PSObject}
        Returns
        -------
        urn2pathway = {urn:PSObject}
        """
        print('Retrieving identifiers of all pathways from database may take couple minutes')

        if not hasattr(self,'id2folder'): self.id2folder = self.load_folder_tree()

        urn2fobj = dict()
        folders =  [PSObject.from_zeep(folder) for folder in self.id2folder.values()]
        if from_folders:
            folders = {o for o in folders if o.name() in from_folders}

        for folder in folders:
            zeep_objects = self.get_folder_objects_props(folder['Id'],with_props)
            folder_objs = self._zeep2psobj(zeep_objects).values()

            # annotating all folder object with "Folders" property:
            for fobj in folder_objs:
                urn2fobj[fobj.urn()] = fobj
                fobj.update_with_value('Folders', folder['Name'])
                if fobj.objtype() == 'Pathway':
                    fobj.update_with_value('layout', self.get_layout(fobj.dbid()))

                if fobj.objtype() == 'attributesSearch':
                    fobj.set_property('URN', self.__urn4result(fobj))
        
        self.dbid2folder_obj.update(folder_objs)

        print('Found %d pathways,%d groups, %d results in the database' % (len(self.dbid2folder_obj), len(self.id2group), len(self.id2result)))
        return urn2fobj

    
    def download_pathways_by_urn(self, pathway_urns:list, fout:str, format:str='RNEF',from_folders=[]):
        global_start = time.time()
        urn2pathway = self.load_containers(from_folders=from_folders) 
        # works 40 sec if "from_folders" is not specified
        print(f'List of all pathways identifiers was retrieved in {self.execution_time(global_start)}')
        
        missedURNs = list()
        download_counter = 0
        pathways2return = list()

        file_ext = '.txt'
        if format == 'RNEF':
            file_ext = '.rnef'
        elif format == 'SBGN':
            file_ext = '.sbgn'
        elif format == 'JSON-LD':
            file_ext = '.jsonld'
        elif format == 'json':
            file_ext = '.json'

        with open(fout+file_ext, 'w',encoding='utf-8') as f:
            if format == 'RNEF': f.write('<batch>')
            for urn in pathway_urns:
                try: 
                    pathway = urn2pathway[urn]
                    start_time = time.time()
                    pathway_graph, pathway_str = self.load_folder_obj(pathway)[1:2]
                    f.write(pathway_str)
                    pathway = PSPathway(urn2pathway[urn],pathway_graph)
                    pathway[format] = pathway_str
                    pathways2return.append(pathway)

                    download_counter += 1
                    exec_time = self.execution_time(start_time)
                    print('%d out of %d pathways downloaded in %s, %d not found' %
                        (download_counter, len(pathway_urns), exec_time, len(missedURNs)))                
                except KeyError:
                    missedURNs.append(urn)
                    continue
            if format == 'RNEF': f.write('</batch>')

        print('Pathways not found:\n%s' % '\n'.join(missedURNs))
        print("Total download time: %s" % self.execution_time(global_start))
        return list(map(PSPathway.from_pathway,pathways2return))


    def find_pathways(self, for_entities:list, in_folders:list):
        """
        Input [PSObject]
        -----

        Return {entity_id:PSObject},
        -------
        where PSObject has 'Pathway ID' property

        Loads self.dbid2folder_obj = {pathway_id:PSObject}
        -----
        where PSObject type is one from FOLDER_OBJECTS annotated with 'Folders' property containing folder names
        """
        ent_dbids = list(map(str,ResnetGraph.dbids(for_entities)))
        to_return = dict()
        for folder in in_folders:
            folder_id = self.__folder_id(folder)
            sub_folder_ids = self.subfolder_ids(folder_id)
            self.get_objects_from_folders(list(sub_folder_ids)) # loads dbid2pathway
        
        for pathway_id in self.dbid2folder_obj.keys():
            id2psobj = self.get_pathway_members([pathway_id],None,ent_dbids,['id'])
            for id, psobj in id2psobj.items():
                try:
                    to_return[id].update_with_value(PATHWAY_ID,pathway_id)
                except KeyError:
                    psobj.update_with_value(PATHWAY_ID,pathway_id)
                    to_return[id] = psobj

        return to_return


    def folder2pspathways(self,folder_id_or_name:int or str,with_layout=True):
        """
        Input
        -----
        Either "folder_id" or "folder_name" must be supplied. 
        If "folder_id" is supplied "folder_name" is retreived from database

        Returns
        -------
        list of PSPathway objects from folder_id_or_name annotated with 'resnet' and 'Folders' properties
        """
        folder_id = self.__folder_id(folder_id_or_name) if isinstance(folder_id_or_name, str) else folder_id_or_name
        folder_name = self.id2folder[int(folder_id)]['Name']
        sub_folder_ids = self.subfolder_ids(folder_name)
        my_folders_ids = list(sub_folder_ids) + [folder_id]
        self.get_objects_from_folders(my_folders_ids,with_layout)
        if self.dbid2folder_obj:
            print('Start downloading %d pathways from \"%s\" folder' % (len(self.dbid2folder_obj),folder_name))
        else:
            print('Folder \"%s\" has no pathways' % folder_name)
            return list()
        
        thread_name = f'Loading {len(self.dbid2folder_obj)} pathways'
        pspathways2return = list()
        with ThreadPoolExecutor(max_workers=MAX_SESSIONS, thread_name_prefix=thread_name) as e: 
            futures = list()
            for dbid, folder_obj in self.dbid2folder_obj.items():
                if folder_obj.objtype() == 'Pathway':
                    futures.append(e.submit(self.load_folder_obj,folder_obj,dict(),dict(),False))
                
            for f in as_completed(futures):
                pathway_obj, pathway_graph,pathway_xml = f.result()
                ps_pathway = PSPathway(dict(pathway_obj),pathway_graph)
                ps_pathway[RESNET] = pathway_xml
                ps_pathway.update_with_value('Folders',folder_name)
                pspathways2return.append(ps_pathway)
 
        return pspathways2return
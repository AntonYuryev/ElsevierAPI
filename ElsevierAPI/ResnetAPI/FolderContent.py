from .ResnetAPISession import APISession, time, len,SNIPPET_PROPERTIES,glob,Path
from .ResnetGraph import ResnetGraph,RESNET,execution_time,OBJECT_TYPE,PSObject,nx, PSRelation, REGULATORS, TARGETS
from collections import Counter,defaultdict
from .PSPathway import PSPathway
from lxml import etree as et
from urllib.parse import quote
from concurrent.futures import ThreadPoolExecutor,as_completed
from ..utils import run_tasks,pretty_xml
import os


PATHWAY_ID = 'Pathway ID'
SYMLINK = 'Symlink'
FOLDER_OBJECTS = ['Pathway','Group','attributesSearch']
FOBJECTS_PROPS = ['Name','Description','Notes']
MAX_FOLDER_SESSIONS = 8

'''
NOTAION:
folder - zeep_object for database folder,
fobj - folder object: PSObject of type Pathway,Group,attributesSearch
'''

class FolderContent(APISession): 
    pass
    def __init__(self,  *args,**kwargs):
        my_kwargs = {'preload_folder_tree':True,'what2retrieve':SNIPPET_PROPERTIES}
        my_kwargs.update(kwargs)
        self.root_folder = self.filename4(my_kwargs.pop('root_folder','')) # default root folder
        super().__init__(*args,**my_kwargs)
        self.SOAPclient.transport.load_timeout = 600 # pathways tend to have relations with large number of references
        self.PageSize = 1000
        self.DumpFiles = []
        self.id2folder = dict() # {int:dict} = {folder_id:folder}
        self.dbid2fobj = dict() # {dbid:PSObject}
        self.parent2childs = dict()
        self.subfold_id2fold_id = dict()
        # where PSObject corresponds to types in FOLDER_OBJECTS: Pathway, Group, attributesSearch
        self.fobj_counter = Counter()
        self.FolderGraph = ResnetGraph() # graph of folders, used to store folder tree

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
      new_session.FolderGraph = self.FolderGraph
      new_session.root_folder = self.root_folder
      return new_session


    def __folder_id(self,folder_name:str) -> int:
        '''
        loads:
          self.id2folder = {dbid:folder_zobj} if necessary
        output:
          folder dbid
        '''
        if not hasattr(self, "id2folder") or self.id2folder is None: 
          self.id2folder = self.load_folder_tree()
        return next((int(k) for k, v in self.id2folder.items() if v['Name'] == folder_name), 0)
    

    def __folder_name(self,folder_dbid:int):
        return  str(self.id2folder[folder_dbid]['Name'])
    

    @staticmethod
    def __urn4result(o:PSObject):
        return 'urn:agi-result:'+quote(o['Name'][0])


    @staticmethod
    def __objtype4result(result_graph:ResnetGraph):
        return 'Pathway' if result_graph.number_of_edges() else 'Group'


    def load_subfolders(self, with_dbids:list[int])->set[int]:
        '''
        Return:
            subfolder_ids
        '''
        # FolderIds - list of folder dbids
        if not hasattr(self, "id2folder"): 
            self.id2folder = self.load_folder_tree()

        subfolders_ids = set()
        for folder_id in with_dbids:
            folder = self.id2folder[folder_id]
            if type(folder['SubFolders']) != type(None):
                subfolders_ids.update(folder['SubFolders']['long']) # folder['SubFolders']['long'] = [int]

        return subfolders_ids


    def subfolder_idsOLD(self, folder_name:str, include_parent=True):
        """
        Returns {subfolders_ids+parent_folder_id}
        -------  
        """
        FolderId = self.__folder_id(folder_name)
        accumulate_subfolder_ids = {FolderId} if include_parent else set()
        subfolder_ids = set(self.load_subfolders([FolderId]))
        while subfolder_ids:
            accumulate_subfolder_ids.update(subfolder_ids)
            subs = set()
            [subs.update(self.load_subfolders([db_id])) for db_id in subfolder_ids]
            subfolder_ids = subs
        return accumulate_subfolder_ids
    

    def subfolder_ids(self, folder_name:str, include_parent:bool=True) -> set[int]:
        """
        Returns a set of subfolder IDs, optionally including the parent folder ID.
        """
        folder_id = self.__folder_id(folder_name)
        all_subfolder_ids = {folder_id} if include_parent else set()

        stack = [folder_id]
        while stack:
            current_id = stack.pop()
            children = self.load_subfolders([current_id])
            all_subfolder_ids.update(children)
            stack.extend(children)

        return all_subfolder_ids


    def get_subfolder_trees(self,folder_name:str)->tuple[dict[int,int],dict[int,list[int]]]:
        """
        loads:
          self.id2folder = {dbid:folder_zobj}
          self.child2parent = {folder_id:folder_id}, self.parent2childs = {parent_id:[child_ids]}
          self.FolderGraph = ResnetGraph() - graph of folders
        """
        if not hasattr(self, "id2folder"): self.id2folder = self.load_folder_tree()

        folder_id = self.__folder_id(folder_name)
        parent2childs = defaultdict(list)
        if self.id2folder[folder_id]['SubFolders'] is None:
            return {folder_id: folder_id}, dict()
            
        stack = [(folder_id, self.id2folder[folder_id]['SubFolders']['long'])]
        while stack:
            parent, children = stack.pop()
            parent2childs[parent].extend(children)  # Directly extend for brevity
            for child in children:
                if self.id2folder[child]['SubFolders'] is not None:
                    stack.append((child, self.id2folder[child]['SubFolders']['long']))
        
        child2parent = {child: parent for parent, children in parent2childs.items() for child in children}
        print(f'Loaded folder tree with {len(parent2childs)} parent folders and {len(child2parent)} subfolders')
        self.subfold_id2fold_id = dict(sorted(child2parent.items()))
        self.parent2childs = dict(sorted(parent2childs.items()))
        self.__foldtree2graph()  # Create the graph from the folder tree
        return self.subfold_id2fold_id, self.parent2childs


    def __objects_from_folders(self, FolderIds:list[int], with_layout=False)->dict[int,PSObject]:
        '''
        Return:
            {dbid:PSObject}, where PSObject.objtype() in ['Pathway','Group','attributesSearch'] 
            PSObject are annotated with attributes "Folder" and "layout"
        '''
        updated_dbid2fobj = dict()
        for fid in FolderIds:
            folder_name = self.__folder_name(fid)
            zeep_objects = self.get_folder_objects_props(fid, FOBJECTS_PROPS)
            self.dbid2fobj = self._zeep2psobj(zeep_objects)
            # annotating all folder object with "Folders" property:
            for k,o in self.dbid2fobj.items():
                o.update_with_value('Folders', folder_name)
                updated_dbid2fobj[k] = o

        if with_layout:
            for k,o in self.dbid2fobj.items():
                o.update_with_value('layout', self.get_layout(o.dbid()))
                updated_dbid2fobj[k] = o

        return updated_dbid2fobj
    

    def __make_cache_dir(self, root_folder_id):
        """
        Creates a directory structure in self.data_dir mirroring the folder structure in the database.
        """
        stack = [(root_folder_id, '')]  # Stack to manage folders (id, parent_name)
        folder_count = 1
        while stack:
          folder_id, parent_name = stack.pop()  # Get the next folder to process
          folder_name = self.id2folder[folder_id].Name
          self.dump_path(folder_name, parent_name,self.root_folder)  # Create the folder
          folder_count += 1
          subfolders_ids = self.load_subfolders([folder_id])  # Get subfolders
          [stack.append((sub_id, folder_name)) for sub_id in subfolders_ids] # Add subfolders to the stack
        print(f'Created cache directory in {self.data_dir} with {folder_count} subfolders')
        return
    
    
    @staticmethod
    def __folder2psobj(folder):
       '''
        Converts a folder zeep object to a PSObject with URN and Id properties.
       '''
       folder_psobj = PSObject({
          'Id':[folder.Id],
          'URN':['urn:agi-folder:'+str(folder.Id)],
          'Name':[folder.Name],
          OBJECT_TYPE:['Folder']
        })
       return folder_psobj
    

    def __foldtree2graph(self):
      '''
      input:
        self.parent2childs
      '''
      psrels = []
      for parent_id, child_ids in self.parent2childs.items():
        parent_psobj = self.__folder2psobj(self.id2folder[parent_id])
        child_psobjs = [self.__folder2psobj(self.id2folder[child_id]) for child_id in child_ids]
        for child in child_psobjs:
          rel = PSRelation({OBJECT_TYPE:['MemberOf']})
          rel.Nodes[REGULATORS] = [child]
          rel.Nodes[TARGETS] = [parent_psobj]
          psrels.append(rel)

      self.FolderGraph = ResnetGraph.from_rels(psrels)
      return self.FolderGraph


    def __foldtree2rnef(self)->str:
      assert(self.FolderGraph), 'FolderGraph is not initialized. Call get_subfolder_trees() first'
      subtree_xml = self.FolderGraph.to_rnefstr(['Name'],['Source'])
      subtree_xml = pretty_xml(subtree_xml,remove_declaration=True)
      # need fake property 'Source' to avoid printing 'URN'and ObjTypeName attribute
      return subtree_xml
        

    def make_cache_dir(self, root_folder_name:str)->list:
        """
        loads:
          self.subfold2fold = {folder_id:folder_id}
        writes:
          folder tree in RNEF format to self.data_dir/root_folder_name
        output:
          [Folder Zeep Objects] with all subfolders of root_folder_name
        """
        if not self.subfold_id2fold_id:
          self.get_subfolder_trees(root_folder_name)

        self.__make_cache_dir(self.__folder_id(root_folder_name))
        subtree_xml = self.__foldtree2rnef()
        self.rnefs2dump(subtree_xml,root_folder_name)
        self.close_rnef_dump(root_folder_name)
        return [folder for dbid,folder in self.id2folder.items() if dbid in self.subfold_id2fold_id]

    
    def update_dbid2fobj(self, dbids4folders:list[int], with_layout=False):
        """
        Updates
        -----
        self.dbid2fobj = {dbid:PSObject}, where PSObject.objtype() in ['Pathway','Group','attributesSearch']

        Return
        ------
        {dbid:PSObject} folder objects from dbids4folders
        """
        my_fobjs = self.__objects_from_folders(dbids4folders,with_layout)
        self.dbid2fobj.update(my_fobjs)
        return my_fobjs


    def __download_fobj(self,fobj:PSObject,props4rel=dict(),props4pathway=dict(),prettify=True)->tuple[PSObject,ResnetGraph,str]:
        '''
        Return:
            tuple (fobj:PSObject,fobj_graph:ResnetGraph,fobj_xml:str)\n
            if fobj.objtype() == 'attributesSearch' adds URN and ObjTypeName
        '''
        fobj_download_start = time.time()
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
            attachments = self.get_layout(fobj.dbid())
            if attachments:
                lay_out = et.Element('attachments')
                lay_out.append(et.fromstring(attachments))
                rnef.append(lay_out)
          
        xml_str = str(et.tostring(rnef,encoding='utf-8').decode("utf-8"))
        if prettify:
            # pretty_xml reuires xml_str to have xml_declaration:
            xml_str = pretty_xml(xml_str,remove_declaration=True)
        print(f'Downloaded "{fobj.name()}" {str(fobj.objtype()).lower()} with {len(obj_graph)} nodes \
and {obj_graph.number_of_edges()} relations in {execution_time(fobj_download_start)}')
        return my_fobj,obj_graph,xml_str


    def __fobj2folder_rnef(self,fobjs:list[PSObject],_2folder):
      folder_local_id = 'F0'
      folder_resnet = et.Element(RESNET)
      folder_nodes = et.SubElement(folder_resnet, 'nodes')
      member_controls = et.SubElement(folder_resnet, 'controls')
      xml_node_folder = et.SubElement(folder_nodes, 'node', {'local_id':folder_local_id, 'urn':'urn:agi-folder:'+str(_2folder.Id)})
      et.SubElement(xml_node_folder, 'attr', {'name':'NodeType', 'value':'Folder'})
      et.SubElement(xml_node_folder, 'attr', {'name':'Name', 'value':_2folder.Name})
      for fobj_count, fobj in enumerate(fobjs):
        pathway_local_count = str(fobj_count)
        pathway_local_id = 'P'+pathway_local_count
        folder_pathway_node = et.SubElement(folder_nodes, 'node', {'local_id':pathway_local_id, 'urn':fobj.urn()})
        et.SubElement(folder_pathway_node, 'attr', {'name': 'NodeType', 'value': fobj.objtype()})
        et.SubElement(folder_pathway_node, 'attr', {'name': 'Name', 'value': fobj.name()})
      
        control_local_id = 'L'+pathway_local_count
        member_control = et.SubElement(member_controls, 'control', {'local_id':control_local_id})
        et.SubElement(member_control, 'attr', {'name':'ControlType', 'value':'MemberOf'})
        if fobj['IsSymlink'][0]:
            et.SubElement(member_control, 'attr', {'name':'Relationship', 'value':'symlink'})
        et.SubElement(member_control, 'link', {'type':'in', 'ref':pathway_local_id})
        et.SubElement(member_control, 'link', {'type':'out', 'ref':folder_local_id})
        fobj_count += 1

      folder_rnef_s = et.tostring(folder_resnet, encoding='utf-8',xml_declaration=False).decode("utf-8")
      folder_rnef_s = pretty_xml(folder_rnef_s,remove_declaration=True)
      return folder_rnef_s


    def folderpath_rnef(self, folder_dbid:int)->str:
      """
      Returns RNEF XML for folder path
      """
      folder_urn = 'urn:agi-folder:'+str(folder_dbid)
      folder_psobj = self.FolderGraph._psobjs_with(folder_urn,'URN')[0]
      folder_uid = folder_psobj.uid()
      path = nx.dfs_tree(self.FolderGraph, source=folder_uid)

      path_rels = []
      for r,t in path.edges():
        rels = [v['relation'] for v in dict(self.FolderGraph[r][t]).values()]
        path_rels += rels

      pathG = ResnetGraph.from_rels(path_rels)
      path_rnef_s = pathG.to_rnefstr(['Name'],['Source'])
      # need fake property 'Source' to avoid printing 'URN'and ObjTypeName attribute
      path_rnef_s = pretty_xml(path_rnef_s,remove_declaration=True)
      return path_rnef_s
    

    def __download_folder(self,folder_id_or_name:int|str,parent_folder_name:str,fobjs2skip:set[PSObject]=set(),
      props4rel:dict[str,str]=dict(),props4pathway:dict[str,str]=dict())->tuple[set[PSObject],str,str,Counter]:
      """
      Input:
        either folder_id or folder_name must be supplied
        if folder_id is supplied folder_name is retreived from self.id2folder cache\n
        fobjs2skip - {PSObject} - list of fobjs that were downloaded previously. Generated by __inspect()
        used to avoid duplications by resolve resume_download() 
      Updates:
        self.dbid2fobj = {id:PSObject}\n
        fobjs2skip - to find unsupported symlinks that will be printed at the end of a download
      Dumps:
        Objects from "folder_id_or_name" into 'folder_name' inside 'parent_folder_name' located in "self.data_dir"
        if size of dump file exceeds 100000000, "rnef_xml" is splitted into several RNEF files\n
        dump RNEF files are named as: 'content of folder_name#',
        where # - dump file number
      output:
        new_pathway_counter,symlinks,folder_name,parent_folder_name,fobj_counter
      """
      folder_id = self.__folder_id(folder_id_or_name) if isinstance(folder_id_or_name, str) else folder_id_or_name
      my_folder = self.id2folder[folder_id] # zeep object for folder
      my_fobjs = set(self.update_dbid2fobj([folder_id],with_layout=True).values())
      if not my_fobjs:
        print(f'\nFolder "{my_folder.Name}" has no objects')
        return set(),my_folder.Name,parent_folder_name,Counter()
      else:
        all_fobjs_count = len(my_fobjs)
        print(f'\n"{my_folder.Name}" folder in "{parent_folder_name}" has {all_fobjs_count} objects')
    
      folderpath_rnef = self.folderpath_rnef(my_folder.Id)
      fobj_counter = Counter()
      symlinks = set()
      fobjs2download = list()
      for fobj in my_fobjs:
        if fobj['IsSymlink'][0]:
          symlinks.add(fobj)
        elif fobj in fobjs2skip:
          fobj_counter[fobj.objtype()] += 1
        else:
          fobjs2download.append(fobj)

      if symlinks:
        symlink_count = len(symlinks)
        # overrides existing symlink file:
        my_folder_name = self.filename4(my_folder.Name)
        my_folder_path = self.dump_path(my_folder.Name,parent_folder_name,self.root_folder)
        print(f'Dumping {symlink_count} symlinks for "{my_folder.Name}" folder into "{my_folder_path}"')
        with open(os.path.join(my_folder_path,my_folder_name+'_symlinks.rnef'), 'w',encoding='utf-8') as symf:
          symf.write('<batch>\n')
          symf.write(self.__fobj2folder_rnef(symlinks,my_folder))
          symf.write('\n'+folderpath_rnef)
          symf.write('\n</batch>')
          fobj_counter[SYMLINK] += symlink_count

      fobjs2dload = len(fobjs2download)
      if fobjs2dload:
        print(f'{all_fobjs_count-fobjs2dload} objects for {my_folder.Name} were downloaded previously')
        print(f'Begin download remaining {fobjs2dload} objects for {my_folder.Name}')
      else:
        print(f'All objects in {my_folder.Name} folder were downloaded previosly. Skipping download')
        return symlinks,my_folder.Name,parent_folder_name,fobj_counter       
        
      # multithreading download of folder objects is a bad idea 
      # because pathways in folder often share same relations and nodes
      # their simultaneous downlaod will be lock database
      start_time = time.time()
      added_fobjs = []
      for fobj in fobjs2download:
        fobj_type = fobj.objtype()
        if fobj_type in FOLDER_OBJECTS:
          added_fobjs.append(fobj)
          fobj_counter[fobj_type] += 1
          _,_,fobj_xml  = self.__download_fobj(fobj,props4rel,props4pathway,True)
          my_parentfolder_name = '' if parent_folder_name == my_folder.Name else parent_folder_name
          path2file = self.rnefs2dump(fobj_xml+'\n',my_folder.Name,my_parentfolder_name,self.root_folder,can_close=False)
          fobjs2skip.add(fobj)
        else:
          print (f'{my_folder.Name} folder has object with id={fobj.dbid()} of unknown type "{fobj.objtype()}"')
          continue

      fobjs2folder_rnef = self.__fobj2folder_rnef(added_fobjs,my_folder)
      folder_rnef = fobjs2folder_rnef + '\n' +folderpath_rnef
      path2file = self.rnefs2dump(folder_rnef,my_folder.Name,my_parentfolder_name,self.root_folder,can_close=False)

      print('%d out of %d objects in folder "%s" with %d symlinks were downloaded to %s in %s' %
  (sum(fobj_counter.values()), len(self.dbid2fobj),my_folder.Name,len(symlinks),path2file,execution_time(start_time)))
      return symlinks,my_folder.Name,parent_folder_name, fobj_counter


    def __download_folders(self,folders:list,job_name:str,downloaded_fobjs:set[PSObject]=set(),
                           props4rel=dict(),props4pathway=dict())->tuple[set[PSObject],Counter]:
        '''
        input:
          folders - list of zeep objects

        output:
            symlinks to be printed. Updates downloaded_fobjs
        '''
        chunks = list()
        symlinks = set()
        start_time = time.time()
        chunks = [folders[i:i+MAX_FOLDER_SESSIONS] for i in range(0, len(folders), MAX_FOLDER_SESSIONS)]
        fobj_counter = Counter()

        print(f'Will download {len(folders)} folders by {MAX_FOLDER_SESSIONS} threads in {len(chunks)} iterations')
        for i,chunk in enumerate(chunks):
            thread_name = f'{job_name} chunk #{i+1} of {len(chunks)}'
            with ThreadPoolExecutor(max_workers=MAX_FOLDER_SESSIONS, thread_name_prefix=thread_name) as e:
              futures = list()
              sessions = list()
              for folder in chunk:
                parent_id = self.subfold_id2fold_id[folder.Id]
                parent_folder_name = self.id2folder[parent_id].Name
                new_session = self.clone(copy_graph=True)
                futures.append(e.submit(new_session.__download_folder,folder.Id,parent_folder_name,downloaded_fobjs,props4rel,props4pathway))
                sessions.append(new_session)

              future_to_index = {future: i for i, future in enumerate(futures)}
              for f in as_completed(futures):
                  index = future_to_index[f]
                  sessions[index].close_connection() # close connection ASAP to avoid locking database
                  self.Graph = self.Graph.compose(sessions[index].Graph)

                  syms,subfolder_name,parent_folder_name,fobc = f.result()
                  fobj_counter.update(fobc)
                  fobj_counter['Folder'] += 1
                  self.close_rnef_dump(subfolder_name,parent_folder_name,self.root_folder,True)

                  print(f'\nFinished downloading "{subfolder_name}" subfolder in "{parent_folder_name}" folder: {fobj_counter['Folder']} out of {len(folders)}')
                  if fobj_counter:
                    print(f'Download {fobj_counter}')

                  symlinks.update(syms)                  
              e.shutdown()

            if self.Graph.weight() > self.reference_cache_size:
                print('Clearing cache due to large size: %d' % self.Graph.weight())
                self.clear()

        print(f'{job_name} downloaded {sum(fobj_counter.values())} objects from {len(folders)} out of {len(self.id2folder)} folders in {execution_time(start_time)}')
        print('Graph cache has %d references\n' % self.Graph.weight())
        return symlinks.difference(downloaded_fobjs),fobj_counter #only symlinks that were not printed 


########################################  RESUME DOWNLOAD #########################################
    def __inspect_file(self,rnef_file:str)->tuple[set[PSObject],Counter]:
        '''
        Return:
            set of fake pathway objects identified by URN present in rnef_file
        '''
        fobjs = set()  # Using a set directly to avoid duplicate URNs
        fobj_counter = Counter()
        with open(rnef_file, "r", encoding='utf-8') as f:
          content = f.read().strip()
        
        # Unfinished RNEF files will be missing </batch> tag
        if not content.endswith('</batch>'):
          content += '</batch>'  # Ensure the file ends with a closing tag. 
            
        for resnet in et.fromstring(content).findall('./resnet'):
          ps_pathway = PSPathway.from_resnet(resnet, ignore_graph=True)
          objtype = ps_pathway.objtype()
          if objtype in FOLDER_OBJECTS:
            fobjs.add(PSObject({'URN':[ps_pathway.urn()],OBJECT_TYPE:[objtype]}))
            fobj_counter[objtype] += 1
        print(f'\nFound {len(fobjs)} objects in "{rnef_file}"')
        return fobjs, fobj_counter
    

    def __inspect_symlinks(self,symlink_file:str):
        members_graph = ResnetGraph.fromRNEF(symlink_file)
        symlinks = members_graph.psobjs_with(only_with_values=FOLDER_OBJECTS)
        return symlinks


    def _inspect_dir(self,dirname:str)->tuple[bool,set[PSObject],set[PSObject],Counter]:
        '''
        output:
          path_exist,downloaded_fobjs,symlinks,fobj_counter
        '''
        downloaded_fobjs = set()
        dir_symlinks = set()
        fobj_counter = Counter()
        path_exist = Path(dirname).exists()
        if path_exist:
          listing = glob.glob(dirname+"/**/*.rnef", recursive=True)
          if listing:
            max_workers = min(32,len(listing),(os.cpu_count() or 1) + 4)
            with ThreadPoolExecutor(max_workers=max_workers, thread_name_prefix='InspectDir') as i:
              futures = list()
              for rnef_file in listing:
                if 'symlink' in rnef_file:
                  file_symlinks = self.__inspect_symlinks(rnef_file)
                  dir_symlinks.update(file_symlinks)
                else:
                  futures.append(i.submit(self.__inspect_file,rnef_file))

              for f in as_completed(futures):
                fobj, counter = f.result()
                downloaded_fobjs.update(fobj)
                fobj_counter.update(counter)
              i.shutdown()

        fobj_counter[SYMLINK] = len(dir_symlinks)
        return path_exist,downloaded_fobjs,dir_symlinks,fobj_counter


    def download_dir(self):
       return os.path.join(self.data_dir,self.root_folder)


    def resume_download(self, root_folder_name:str):
      '''
        Inspects download directory and resumes download of folders and objects that have not beed downloaded in the previous run
      '''
      global_start = time.time()
      assert(self.data_dir+self.root_folder), 'Data directory is not set. Use set_dir() first'
      # inspect download directory and load folder structure from DB:
      tasks = [(self.get_subfolder_trees,(root_folder_name,)), (self._inspect_dir,(self.download_dir(),))]
      results = run_tasks(tasks)
      self.subfold_id2fold_id,_ = results['get_subfolder_trees']
      path_exist,downloaded_fobjs,_,self.fobj_counter = results['_inspect_dir']
      print(f'Stats from previous download attempt: {self.fobj_counter}')
      
      if not path_exist:# case when nothing was downloaded yet
        self.make_cache_dir(root_folder_name)
    
      subfolders = [folder for dbid,folder in self.id2folder.items() if dbid in self.subfold_id2fold_id]
      # subfolders are zeep objects
      job = 'Resume download'
      symlinks4topfolder,fobj_counter = self.__download_folders(subfolders,job,downloaded_fobjs)
      self.fobj_counter.update(fobj_counter)

      if symlinks4topfolder:
        print(f'Will print {len(symlinks4topfolder)} pathways into top-level folder to support symlinks')
        # dumping loose symlinks into root folder
        for symlink in symlinks4topfolder:
          _,_,xml_str = self.__download_fobj(symlink) #downloading missing fobjs necessary to support symlinks
          self.rnefs2dump(xml_str,root_folder_name,'',self.root_folder)
        self.close_rnef_dump(root_folder_name,'',self.root_folder)
      else:
        print('All pathways for symlinks were downloaded into other folders')
          
      print(f'Total download stats: {self.fobj_counter}')
      print(f'Resumed download was finished in {execution_time(global_start)}')
      

############################  OTHER DOWNLOAD FUNCTIONS #############################################
    def folder2rnef(self,root_folder_name:str,include_subfolders=True,
                    add_props2rel=dict(),add_pathway_props=dict()):
        """
        Dumps
        -----
        content of root_folder_name into RNEF file 'content of root_folder_name.rnef' located in 'self.data_dir'
        """
        print(f'Downloading content of {root_folder_name} into {self.data_dir}')
        if include_subfolders:
            print(f'Including subfolders of {root_folder_name} into download')
        download_start_time = time.time()
        if not include_subfolders:
            return self.__download_folder(root_folder_name,'',props4rel=add_props2rel,props4pathway=add_pathway_props)
        else:
            subfolders = self.make_cache_dir(root_folder_name)
            thread_name = f'{root_folder_name} download'
            symlinks2print,fobj_counter = self.__download_folders(subfolders,thread_name)
        
            if symlinks2print:
                print(f'Will print {len(symlinks2print)} pathways to support symlinks')
                # dumping loose symlinks into root folder
                for fobj in symlinks2print:
                    fobj, pathway_graph, pathway_xml = self.__download_fobj(fobj)
                    self.rnefs2dump(pathway_xml,root_folder_name)
                self.close_rnef_dump(root_folder_name)
            else:
                print(f'No additional symlinks are required for {root_folder_name} folder %s')
                
        print(f"Complete download execution time: {execution_time(download_start_time)}")


    def load_containers(self,with_props=['Name','Description','Notes'],from_folders=[]):
        """
        retreives entire folder tree stricture from database including all folder objects
        if from_folders is not specidied works for ~40sec

        Updates
        -----
        self.dbid2fobj = {id:PSObject}

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

        pathway_counter = 0
        group_counter = 0
        result_counter = 0
        for folder in folders:
            zeep_objects = self.get_folder_objects_props(folder['Id'],with_props)
            folder_objs = self._zeep2psobj(zeep_objects)

            # annotating all folder object with "Folders" property:
            for dbid,fobj in folder_objs.items():
                fobj.update_with_value('Folders', folder['Name'])
                if fobj.objtype() == 'Pathway':
                    fobj.update_with_value('layout', self.get_layout(fobj.dbid()))
                    pathway_counter += 1
                elif fobj.objtype() == 'attributesSearch':
                    fobj.set_property('URN', self.__urn4result(fobj))
                    result_counter += 1
                else:
                    group_counter += 1

                urn2fobj[fobj.urn()] = fobj
        
        self.dbid2fobj.update(folder_objs)

        print('Found %d pathways,%d groups, %d results in database' % 
              pathway_counter,group_counter,result_counter)
        return urn2fobj

    
    def download_pathways_by_urn(self, pathway_urns:list, fout:str, format:str='RNEF',from_folders=[]):
        '''
        fout - dump filename 
        '''
        global_start = time.time()
        urn2pathway = self.load_containers(from_folders=from_folders) 
        # works 40 sec if "from_folders" is not specified
        print(f'List of all pathways identifiers was retrieved in {execution_time(global_start)}')
        
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
                    pathway_obj,pathway_graph,pathway_str = self.__download_fobj(pathway)
                    f.write(pathway_str)
                    pathway = PSPathway(urn2pathway[urn],pathway_graph)
                    pathway[format] = pathway_str
                    pathways2return.append(pathway)

                    download_counter += 1
                    exec_time = execution_time(start_time)
                    print('%d out of %d pathways downloaded in %s, %d not found' %
                        (download_counter, len(pathway_urns), exec_time, len(missedURNs)))                
                except KeyError:
                    missedURNs.append(urn)
                    continue
            if format == 'RNEF': f.write('</batch>')

        print('Pathways not found:\n%s' % '\n'.join(missedURNs))
        print(f"Total download time: {execution_time(global_start)}")
        return list(map(PSPathway.from_pathway,pathways2return))


    def find_pathways(self, for_entities:list[PSObject], in_folders:list[str]):
        """
        Return 
        -------
        {entity_id:PSObject}, where PSObject has 'Pathway ID' property

        Loads
        -----
        self.dbid2fobj = {pathway_id:PSObject}, where PSObject.objtype() in FOLDER_OBJECTS and annotated with 'Folders' property containing folder names
        """
        ent_dbids = list(map(str,ResnetGraph.dbids(for_entities)))
        to_return = dict()
        for folder_name in in_folders:
            #folder_id = self.__folder_id(folder)
            sub_folder_ids = self.subfolder_ids(folder_name)
            self.update_dbid2fobj(list(sub_folder_ids)) # loads dbid2pathway
        
        for pathway_id in self.dbid2fobj.keys():
            id2psobj = self.get_pathway_members([pathway_id],None,ent_dbids,['id'])
            for id, psobj in id2psobj.items():
                try:
                    to_return[id].update_with_value(PATHWAY_ID,pathway_id)
                except KeyError:
                    psobj.update_with_value(PATHWAY_ID,pathway_id)
                    to_return[id] = psobj

        return to_return


    def folder2pspathways(self,folder_id_or_name:int|str,with_layout=True)->list[PSPathway]:
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
        fobjs = self.update_dbid2fobj(my_folders_ids,with_layout)
        if fobjs:
            print('Start downloading %d pathways from \"%s\" folder' % (len(self.dbid2fobj),folder_name))
        else:
            print(f'Folder \"{folder_name}\" has no pathways')
            return list()
        
        thread_name = f'Loading {len(self.dbid2fobj)} pathways'
        pspathways2return = list()
        futures = list()
        with ThreadPoolExecutor(max_workers=MAX_FOLDER_SESSIONS, thread_name_prefix=thread_name) as e:   
          for fobj in self.dbid2fobj.values():
            if fobj.objtype() == 'Pathway':
              futures.append(e.submit(self.__download_fobj,fobj,dict(),dict(),False))
            
        for f in as_completed(futures):
          pathway_obj, pathway_graph,pathway_xml = f.result()
          ps_pathway = PSPathway(dict(pathway_obj),pathway_graph)
          ps_pathway[RESNET] = pathway_xml
          ps_pathway.update_with_value('Folders',folder_name)
          pspathways2return.append(ps_pathway)
        e.shutdown(wait=True)
 
        return pspathways2return
  
    
    @staticmethod
    def _put2folder(named:str, fobj:PSObject,resnet:et.Element):
        xml_nodes = et.SubElement(resnet, 'nodes')
        folder_local_id = 'F0'
        xml_node_folder = et.SubElement(xml_nodes, 'node', {'local_id':folder_local_id, 'urn': 'urn:agi-folder:xxxxx_yyyyy_zzzzz'})
        et.SubElement(xml_node_folder, 'attr', {'name': 'NodeType', 'value': 'Folder'})
        et.SubElement(xml_node_folder, 'attr', {'name': 'Name', 'value': named})
        pathway_local_id = 'P0'
        xml_node_pathway = et.SubElement(xml_nodes, 'node', {'local_id':pathway_local_id, 'urn':fobj.urn()})
        et.SubElement(xml_node_pathway, 'attr', {'name': 'NodeType', 'value': 'Pathway'})
        xml_controls = et.SubElement(resnet, 'controls')
        xml_control = et.SubElement(xml_controls, 'control', {'local_id':'CFE1'})
        et.SubElement(xml_control, 'attr', {'name':'ControlType', 'value':'MemberOf'})
        et.SubElement(xml_control, 'link', {'type':'in', 'ref':pathway_local_id})
        et.SubElement(xml_control, 'link', {'type':'out', 'ref':folder_local_id})
 
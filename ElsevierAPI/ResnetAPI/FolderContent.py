from ElsevierAPI.ResnetAPI.ResnetAPISession import APISession
from ElsevierAPI.ResnetAPI.NetworkxObjects import PSObject, REF_ID_TYPES,SENTENCE_PROPS,RELATION_PROPS
from ElsevierAPI.ResnetAPI.rnef2sbgn import rnef2sbgn_str
import xml.etree.ElementTree as et
from xml.dom import minidom
import time

class FolderContent (APISession): 

    def __init__(self, APIconfig, preload_folder_tree=True):
        super().__init__(APIconfig['ResnetURL'], APIconfig['PSuserName'], APIconfig['PSpassword'])
        self.PageSize = 1000
        self.DumpFiles = []
        self.id2folders = {}
        if preload_folder_tree:
            self.id2folders = self.load_folder_tree() # {folder_id:folder_zobj, folder_name:folder_zobj}
        self.reference_cache_size = 1000000
        #self.id2pathways = dict() # {id:PSObj} loaded on demand
        #self.id2groups = dict() # {id:PSObj} loaded on demand
       
    def get_folder_id(self,folder_name:str) -> list:
        if not hasattr(self, "id2folders"): 
            self.id2folders = self.load_folder_tree()

        return self.id2folders[folder_name][0]['Id']


    def get_subfolders(self, FolderIds: list):
        if not hasattr(self, "id2folders"): 
            self.id2folders = self.load_folder_tree()

        subfolders_ids = set()
        for db_id in FolderIds:
            folders = self.id2folders[db_id]
            for folder in folders:
                if type(folder['SubFolders']) != type(None):
                    subfolders_ids.update(folder['SubFolders']['long'])
        return subfolders_ids

    def get_subfolders_recursively(self, FolderId):
        accumulate_subfolders = {FolderId}
        subfolders = set(self.get_subfolders([FolderId]))
        while len(subfolders) > 0:
            accumulate_subfolders.update(subfolders)
            subs = set()
            for db_id in subfolders:
                subs.update(self.get_subfolders([db_id]))
            subfolders = subs

        return accumulate_subfolders

    def get_subfolder_tree(self,folder_name):
        if not hasattr(self, "id2folders"): 
            self.id2folders = self.load_folder_tree()

        folder_id = self.id2folders[folder_name][0]['Id']
        child2parent = {folder_id:folder_id}
        parent2child = PSObject(dict())

        if type (self.id2folders[folder_id][0]['SubFolders']) != type(None):
            children_ids = self.id2folders[folder_id][0]['SubFolders']['long']
            p2c = {folder_id : children_ids} 
        else: 
            return {folder_id:folder_id}, dict()

        while len(p2c) > 0:
            parent2child.update(p2c)
            p2c_iter = dict()
            for par, chi in p2c.items():
                for c in chi:
                    child2parent[c] = par
                    if type (self.id2folders[c][0]['SubFolders']) != type(None):
                        c_subs = self.id2folders[c][0]['SubFolders']['long']
                        p2c_iter[c] = c_subs
                    else: continue
            p2c = p2c_iter

        return child2parent, dict(parent2child)

    def load_containers(self, property_names=None):
        if property_names is None: property_names = ['Name']
        print('Retrieving identifiers of all pathways from database may take couple minutes')

        if not hasattr(self,'id2folders'): self.id2folders = self.load_folder_tree()

        self.id2pathways = dict()
        self.id2groups = dict()
        urn2pathway = dict()
        for folderList in self.id2folders.values():
            for folder in folderList:
                zeep_objects = self.get_folder_objects_props(folder['Id'], property_names)
                ps_objects = self._zeep2psobj(zeep_objects)
                for Id, psObj in ps_objects.items():
                    if psObj['ObjTypeName'][0] == 'Pathway':
                        try:
                            self.id2pathways[Id].add_unique_property('Folders', folder['Name'])
                        except KeyError:
                            psObj['Folders'] = [folder['Name']]
                            self.id2pathways[Id] = psObj
                            urn2pathway[psObj['URN'][0]] = psObj
                    if psObj['ObjTypeName'][0] == 'Group':
                        try:
                            self.id2groups[Id].add_unique_property('Folders', folder['Name'])
                        except KeyError:
                                psObj['Folders'] = [folder['Name']]
                                self.id2groups[Id] = psObj

        print('Found %d pathways in the database' % (len(self.id2pathways)))
        return urn2pathway


    @staticmethod
    def pretty_xml(xml_string:str, no_declaration = False):
        pretty_xml = str(minidom.parseString(xml_string).toprettyxml(indent='   '))
        if no_declaration:
            pretty_xml = pretty_xml[pretty_xml.find('\n')+1:]
        
        return pretty_xml

    def get_pathway(self, pathwayId,path_urn:str=None,path_name:str=None,rel_props:set=set(), ent_props:set=set(),
                    xml_format='RNEF',put2folder:str=None, add_props2rel:dict=None, add_props2pathway:dict=None, as_batch=True, prettify=True):
    # add_rel_props, add_pathway_props structure - {PropName:[PropValues]}

        get_rel_props = REF_ID_TYPES | rel_props | {'TextRef'} | RELATION_PROPS
        #REF_ID_TYPES are needed to load_references
        # 'TextRef' is needed for references with several supporting snippets
        if not rel_props:
            get_rel_props = get_rel_props | SENTENCE_PROPS
        

        if hasattr(self,'id2pathways'):
            if not isinstance(path_urn,str):
                try:
                    path_urn = self.id2pathways[pathwayId]['URN'][0]
                    path_name = str(self.id2pathways[pathwayId]['Name'][0])
                except KeyError:
                    print('Pathway collection does not have %s pathway with URN %s' % (path_name,path_urn))

        if not isinstance(path_urn,str):
            print('Pathway has no URN specifed!!!! ')
            path_urn = 'no_urn'
        
        if not isinstance(path_name,str):
            print('Pathway has no Name specifed!!!! ')
            path_name = 'no_name'


        pathway_graph = self.get_pathway_components([pathwayId],'id',retrieve_rel_properties=get_rel_props,
                                                    retrieve_ent_properties=ent_props)
        pathway_graph.count_references()

        graph_xml = self.to_rnef(pathway_graph,add_props2rel,add_props2pathway)
        import xml.etree.ElementTree as et
        rnef_xml = et.fromstring(graph_xml)
        rnef_xml.set('name', path_name)
        rnef_xml.set('urn', path_urn)
        rnef_xml.set('type', 'Pathway')

        lay_out = et.Element('attachments')
        attachments = self.get_layout(pathwayId)
        if attachments:
            lay_out.append(et.fromstring(attachments))
        rnef_xml.append(lay_out)
        
        batch_xml = et.Element('batch')
        batch_xml.insert(0,rnef_xml)
                   
        if xml_format == ['SBGN']:
            pathway_xml = et.tostring(batch_xml,encoding='utf-8',xml_declaration=True).decode("utf-8")
            pathway_xml = minidom.parseString(pathway_xml).toprettyxml(indent='   ')
            pathway_xml = rnef2sbgn_str(pathway_xml, classmapfile='ElsevierAPI/ResnetAPI/rnef2sbgn_map.xml')
        else:
            if isinstance(put2folder,str):
                resnet = et.Element('resnet')
                xml_nodes = et.SubElement(resnet, 'nodes')
                folder_local_id = 'F0'
                xml_node_folder = et.SubElement(xml_nodes, 'node', {'local_id':folder_local_id, 'urn': 'urn:agi-folder:xxxxx_yyyyy_zzzzz'})
                et.SubElement(xml_node_folder, 'attr', {'name': 'NodeType', 'value': 'Folder'})
                et.SubElement(xml_node_folder, 'attr', {'name': 'Name', 'value': put2folder})
                pathway_local_id = 'P0'
                xml_node_pathway = et.SubElement(xml_nodes, 'node', {'local_id':pathway_local_id, 'urn': path_urn})
                et.SubElement(xml_node_pathway, 'attr', {'name': 'NodeType', 'value': 'Pathway'})
                xml_controls = et.SubElement(resnet, 'controls')
                xml_control = et.SubElement(xml_controls, 'control', {'local_id':'CFE1'})
                et.SubElement(xml_control, 'attr', {'name':'ControlType', 'value':'MemberOf'})
                et.SubElement(xml_control, 'link', {'type':'in', 'ref':pathway_local_id})
                et.SubElement(xml_control, 'link', {'type':'out', 'ref':folder_local_id})
                batch_xml.append(resnet)
            
            if as_batch:
                pathway_xml = et.tostring(batch_xml,encoding='utf-8',xml_declaration=True).decode("utf-8")
                if prettify: pathway_xml = minidom.parseString(pathway_xml).toprettyxml(indent='   ')
            else:
                pathway_xml = et.tostring(rnef_xml,encoding='utf-8',xml_declaration=True).decode("utf-8")
                if prettify:
                    pathway_xml = str(minidom.parseString(pathway_xml).toprettyxml(indent='   '))
                    pathway_xml = pathway_xml[pathway_xml.find('\n')+1:]
                #minidom does not work without xml_declaration

        print('\"%s\" pathway downloaded: %d nodes, %d edges supported by %d references' % 
            (path_name, pathway_graph.number_of_nodes(),pathway_graph.number_of_edges(),pathway_graph.size(weight="weight")))

        return pathway_graph, str(pathway_xml)


    def get_group(self, group_id,group_urn:str=None,group_name:str=None, ent_props:set=None,put2folder:str=None,as_batch=True, prettify=True):
        if hasattr(self,'id2groups'):
            if not isinstance(group_urn,str):
                try:
                    group_urn = self.id2groups[group_id]['URN'][0]
                    group_name = self.id2groups[group_id]['Name'][0]
                except KeyError:
                    print('Pathway collection does not have %s group with URN %s' % (group_name,group_urn))

        if not isinstance(group_urn,str):
            print('Pathway has no URN specified!!!!')
            group_urn = 'no_urn'
        
        if not isinstance(group_name,str):
            print('Pathway has no Name specified!!!!')
            group_name = 'no_name'

        group_graph = self.load_graph_from_oql('SELECT Entity WHERE MemberOf (SELECT Group WHERE Id = {objectId})'.format(objectId = group_id),entity_props=ent_props,get_links=False)
        rnef_xml = et.fromstring(self.to_rnef(group_graph))
        rnef_xml.set('name', group_name)
        rnef_xml.set('urn', group_urn)
        rnef_xml.set('type', 'Group')
        
        batch_xml = et.Element('batch')
        batch_xml.insert(0,rnef_xml)
                   
        if isinstance(put2folder,str):
            folder_resnet = et.Element('resnet')
            xml_nodes = et.SubElement(folder_resnet, 'nodes')
            folder_local_id = 'F0'
            xml_node_folder = et.SubElement(xml_nodes, 'node', {'local_id':folder_local_id, 'urn': 'urn:agi-folder:xxxxx_yyyyy_zzzzz'})
            et.SubElement(xml_node_folder, 'attr', {'name': 'NodeType', 'value': 'Folder'})
            et.SubElement(xml_node_folder, 'attr', {'name': 'Name', 'value': put2folder})
            pathway_local_id = 'P0'
            xml_node_pathway = et.SubElement(xml_nodes, 'node', {'local_id':pathway_local_id, 'urn':group_urn})
            et.SubElement(xml_node_pathway, 'attr', {'name': 'NodeType', 'value': 'Group'})
            xml_controls = et.SubElement(folder_resnet, 'controls')
            xml_control = et.SubElement(xml_controls, 'control', {'local_id':'CFE1'})
            et.SubElement(xml_control, 'attr', {'name':'ControlType', 'value':'MemberOf'})
            et.SubElement(xml_control, 'link', {'type':'in', 'ref':pathway_local_id})
            et.SubElement(xml_control, 'link', {'type':'out', 'ref':folder_local_id})
            batch_xml.append(folder_resnet)
        
        if as_batch:
            group_xml = et.tostring(batch_xml,encoding='utf-8',xml_declaration=True).decode("utf-8")
            if prettify: group_xml = minidom.parseString(group_xml).toprettyxml(indent='   ')
        else:
            group_xml = et.tostring(rnef_xml,encoding='utf-8',xml_declaration=True).decode("utf-8")
            if prettify: group_xml = self.pretty_xml(group_xml,no_declaration=True)
            #minidom does not work without xml_declaration

        print('\"%s\" group downloaded: %d nodes' % (group_name, group_graph.number_of_nodes()))
        return group_graph, str(group_xml)


    def get_objects_from_folders(self, FolderIds: list, property_names=None, with_layout=False):
        if property_names is None: property_names = ['Name']
        if not hasattr(self,'id2folders'): self.id2folders = self.load_folder_tree()
        if not hasattr(self,'id2pathways'): self.id2pathways = dict()
        if not hasattr(self,'id2groups'): self.id2groups = dict()
        id2objects = dict()
        for f in FolderIds:
            folder_name = self.id2folders[f][0]['Name']
            zeep_objects = self.get_folder_objects_props(f, property_names)
            id2objs = self._zeep2psobj(zeep_objects)
            for Id, psObj in id2objs.items():
                if psObj['ObjTypeName'][0] == 'Pathway':
                    try:
                        self.id2pathways[Id].add_unique_property('Folders', folder_name)
                    except KeyError:
                            psObj['Folders'] = [folder_name]
                            self.id2pathways[Id] = psObj
                    if with_layout:
                        psObj['layout'] = self.get_layout(Id)
                if psObj['ObjTypeName'][0] == 'Group':
                    try:
                        self.id2groups[Id].add_unique_property('Folders', folder_name)
                    except KeyError:
                            psObj['Folders'] = [folder_name]
                            self.id2groups[Id] = psObj

            id2objects.update(id2objs)
        return id2objects

    def pathways_from_folder(self,folder_id:int=None,folder_name:str=None,skip_id:set=None,skip_urn:set=None,
                            add_props2rel:dict=None,add_props2pathway:dict=None,as_batch=True):
        #either folder_id or folder_name must be supplied
        folder_id = folder_id if isinstance(folder_id, int) else self.get_folder_id(folder_name)
        folder_name = self.id2folders[folder_id][0]['Name']
        id2pathways = self.get_objects_from_folders([folder_id],with_layout=True)
        if id2pathways:
            print('Start downloading %d pathways from \"%s\" folder' % (len(id2pathways),folder_name))
        else:
            print('Folder \"%s\" has no pathways' % folder_name)
            return '', 0, {}
        
        folder_resnet = et.Element('resnet')
        folder_nodes = et.SubElement(folder_resnet, 'nodes')
        member_controls = et.SubElement(folder_resnet, 'controls')

        folder_local_id = 'F0'
        xml_node_folder = et.SubElement(folder_nodes, 'node', {'local_id':folder_local_id, 'urn':'urn:agi-folder:'+str(folder_id)})
        et.SubElement(xml_node_folder, 'attr', {'name':'NodeType', 'value':'Folder'})
        et.SubElement(xml_node_folder, 'attr', {'name':'Name', 'value':folder_name})

        new_pathway_counter = 0
        folder_object_counter = 0
        pathways_xml = str()
        symlinks = dict()
        folder_download_start = time.time()
        for pathway_id, pathway in id2pathways.items():
            if not pathway['IsSymlink'][0]:
                #printing pathwys only symlinks should be printed at the end and only if they were not downloaded for another folder
                start_time = time.time()
                folder_object_counter += 1
                if isinstance(skip_id,set):
                    if pathway_id in skip_id:
                        print('%s pathway was downloaded for another folder' % pathway['Name'][0])
                        continue
                if isinstance(skip_urn,set):
                    if pathway['URN'][0] in skip_urn:
                        print('%s pathway was downloaded for another folder' % pathway['Name'][0])
                        continue

                if pathway['ObjTypeName'][0] == 'Pathway':
                    pathway_graph, pathway_xml = self.get_pathway(pathway_id,rel_props=set(self.relProps),add_props2rel=add_props2rel,
                                    ent_props=set(self.entProps),add_props2pathway=add_props2pathway, as_batch=False)
                elif pathway['ObjTypeName'][0] == 'Group':
                    pathway_graph, pathway_xml = self.get_group(pathway_id, ent_props=set(self.entProps),as_batch=False)
                else:
                    print ('%s folder has object with unknown type %s: id = %d' % (folder_name,pathway['ObjTypeName'][0],pathway_id))
                    continue

                pathways_xml = pathways_xml + pathway_xml #adding pathways to XML
                new_pathway_counter += 1
                if isinstance(skip_id,set):skip_id.add(pathway_id)
                if isinstance(skip_urn,set):skip_urn.add(pathway['URN'][0])
                
                print('%d out of %d objects in folder %s downloaded in %s' %
                    (folder_object_counter, len(id2pathways),folder_name, self.execution_time(start_time)))
            else:
                symlinks[pathway_id] = pathway['URN'][0]

            #printing folder object membership
            pathway_urn = pathway['URN'][0]
            pathway_local_id = 'P'+str(folder_object_counter)
            folder_pathway_node = et.SubElement(folder_nodes, 'node', {'local_id':pathway_local_id, 'urn':pathway_urn})
            et.SubElement(folder_pathway_node, 'attr', {'name': 'NodeType', 'value': 'Pathway'})
            et.SubElement(folder_pathway_node, 'attr', {'name': 'Name', 'value': pathway['Name'][0]})
            
            control_local_id = 'L'+str(folder_object_counter)
            member_control = et.SubElement(member_controls, 'control', {'local_id':control_local_id})
            et.SubElement(member_control, 'attr', {'name':'ControlType', 'value':'MemberOf'})
            if pathway['IsSymlink'][0] == True:
                et.SubElement(member_control, 'attr', {'name':'Relationship', 'value':'symlink'})
            et.SubElement(member_control, 'link', {'type':'in', 'ref':pathway_local_id})
            et.SubElement(member_control, 'link', {'type':'out', 'ref':folder_local_id})

        
        return_xml = str()
        if id2pathways:
            folder_rnef = et.tostring(folder_resnet, encoding='utf-8',xml_declaration=False).decode("utf-8")
            folder_rnef = self.pretty_xml(folder_rnef,no_declaration=True)
            return_xml = '<batch>\n'+ pathways_xml + folder_rnef + '</batch>' if as_batch else pathways_xml + folder_rnef
        else:
            return_xml = '<batch/>' if as_batch else ''

        print('Total folder download time: %s' % self.execution_time(folder_download_start))
        return return_xml, new_pathway_counter, symlinks

    @staticmethod
    def name_output(folder_name):
        return 'content of '+folder_name+'.rnef'
    
    def content2rnef(self, folder_name:str,include_subfolders=True, add_props2rel:dict=None,add_pathway_props:dict=None):

        if not include_subfolders:
            return self.pathways_from_folder(None,folder_name,add_props2rel=add_props2rel,add_pathway_props=add_pathway_props)
        else:
            child2parent, parent2child = self.get_subfolder_tree(folder_name)

            global_start = time.time()
            
            with open(self.name_output(folder_name), "w", encoding='utf-8') as f:
                f.write('<?xml version="1.0" ?>\n')
                f.write('<batch>\n')

                #printing folder tree
                for parent_id, child_ids in parent2child.items():
                    folder_resnet = et.Element('resnet')
                    folder_name = self.id2folders[parent_id][0]['Name']
                    folder_nodes = et.SubElement(folder_resnet, 'nodes')
                    xml_controls = et.SubElement(folder_resnet, 'controls')
                    parent_local_id = 'F0'
                    xml_node_folder = et.SubElement(folder_nodes, 'node', {'local_id':parent_local_id, 'urn':'urn:agi-folder:'+str(parent_id)})
                    et.SubElement(xml_node_folder, 'attr', {'name':'NodeType', 'value':'Folder'})
                    et.SubElement(xml_node_folder, 'attr', {'name':'Name', 'value':folder_name})
                    for child_id in child_ids:
                        child_local_id = str(child_id)
                        subfolder_name = self.id2folders[child_id][0]['Name']
                        xml_node_pathway = et.SubElement(folder_nodes, 'node', {'local_id':child_local_id, 'urn':'urn:agi-folder:'+child_local_id})
                        et.SubElement(xml_node_pathway, 'attr', {'name':'NodeType', 'value':'Folder'})
                        et.SubElement(xml_node_pathway, 'attr', {'name':'Name', 'value':subfolder_name})
                        
                        control_local_id = 'L'+str(child_id)
                        xml_control = et.SubElement(xml_controls, 'control', {'local_id':control_local_id})
                        et.SubElement(xml_control, 'attr', {'name':'ControlType', 'value':'MemberOf'})
                        et.SubElement(xml_control, 'link', {'type':'in', 'ref':child_local_id})
                        et.SubElement(xml_control, 'link', {'type':'out', 'ref':parent_local_id})

                    subtree_xml = et.tostring(folder_resnet,encoding='utf-8',xml_declaration=False).decode("utf-8")
                    subtree_xml = self.pretty_xml(subtree_xml,no_declaration=True)
                    f.write(subtree_xml)
                print('Folder tree with %d folders was downloaded in %s' % (len(child2parent)+1, self.execution_time(global_start)))

                #fetching pathways
                download_counter = 0
                folder_counter = 0
                printed_pathway_ids = set()
                symlinks_ids = set()
                for folder_id in child2parent.keys():
                    folder_xml, pathway_counter, symlinks = self.pathways_from_folder(folder_id,as_batch=False,skip_id=printed_pathway_ids)
                    f.write(folder_xml)
                    symlinks_ids.update(symlinks.keys())
                    download_counter += pathway_counter
                    folder_counter +=1
                    if pathway_counter > 0:
                        print('Downloaded %d pathways from %d out of %d folders in %s' % (download_counter,folder_counter,len(child2parent),self.execution_time(global_start)))
                        print('Relations cache has %d relations supported by %d references\n' % (self.Graph.number_of_edges(), self.Graph.weight()))
                    if self.Graph.weight() > self.reference_cache_size:
                        print ('Clearing cache due to size %d' % self.Graph.weight())
                        self.clear()

                symlinks2print = symlinks_ids.difference(printed_pathway_ids)
                if len(symlinks2print) > 0:
                    print('Will print %d pathways to support symlinks' % len(symlinks2print))
                    for pathway_id in symlinks2print:
                        if pathway_id in self.id2pathways.keys():
                            pathway_graph, pathway_xml = self.get_pathway(pathway_id,rel_props=set(self.relProps),ent_props=set(self.entProps), as_batch=False)
                            f.write(pathway_xml)
                        elif pathway_id in self.id2groups.keys():
                            pathway_graph, pathway_xml = self.get_group(pathway_id, ent_props=set(self.entProps),as_batch=False)
                            f.write(pathway_xml)
                        else:
                            continue
                else:
                    print('All pathways for symlinks were printed for other folders')

                f.write('</batch>')
            print("Total execution time: %s" % self.execution_time(global_start))

    @staticmethod
    def load_urns(rnef_file:str):
        urns = set()
        with open(rnef_file, "r", encoding='utf-8') as f:
            line = f.readline()
            while line:
                line = line.strip()
                if line[:8] == '<resnet ':
                    urn_pos = line.find('urn=',8)
                    if urn_pos > 0 : 
                        urn_start = urn_pos+5
                        urn_end = line.find('\"', urn_start)
                        urn = line[urn_start:urn_end]
                        urns.add(urn)
                line = f.readline()
        print('Read \"%s\" with %d pathways' % (rnef_file,len(urns)))
        return urns


    def resume_download(self, folder_name:str, last_downloaded_folder:str):
        child2parent, parent2child = self.get_subfolder_tree(folder_name)
        subfolders = list(child2parent.keys())
        global_start = time.time()
        try:
            last_downloaded_folder_id = self.get_folder_id(last_downloaded_folder)
            start_folder_idx = subfolders.index(last_downloaded_folder_id)
            folder_counter = start_folder_idx+1
            rnef_file_to_continue =  self.name_output(folder_name)
            urns_printed = self.load_urns(rnef_file_to_continue)

            download_counter = len(urns_printed)
            print ('Resuming download of %s folder from %s' % (folder_name,last_downloaded_folder))
            symlinks_dict = dict()
            with open(rnef_file_to_continue, "a", encoding='utf-8') as f:
                for i in range(start_folder_idx+1, len(subfolders)):
                    folder_id = subfolders[i]
                    folder_xml, pathway_counter, symlinks = self.pathways_from_folder(folder_id,as_batch=False,skip_urn=urns_printed)
                    f.write(folder_xml)
                    download_counter += pathway_counter
                    folder_counter +=1
                    symlinks_dict.update(symlinks)
                    if self.Graph.weight() > self.reference_cache_size:
                        print ('Clearing cache due to size %d' % self.Graph.weight())
                        self.clear()
                    print('Downloaded %d pathways from %d out of %d folders in %s' % (download_counter,folder_counter,len(child2parent),self.execution_time(global_start)))
                    print('Graph cache has %d references\n' % self.Graph.weight())

                symlinks_urns = set(symlinks_dict.values())
                symlinks2print = symlinks_urns.difference(urns_printed)
                if len(symlinks2print) > 0:
                    print('Will print %d pathways to support symlinks' % len(symlinks2print))
                    for pathway_urn in symlinks2print:
                        pathway_id = list(symlinks_dict.keys())[list(symlinks_dict.values()).index(pathway_urn)]
                        if pathway_id in self.id2pathways.keys():
                            pathway_graph, pathway_xml = self.get_pathway(pathway_id,rel_props=set(self.relProps),ent_props=set(self.entProps), as_batch=False)
                            f.write(pathway_xml)
                        elif pathway_id in self.id2groups.keys():
                            pathway_graph, pathway_xml = self.get_group(pathway_id, ent_props=set(self.entProps),as_batch=False)
                            f.write(pathway_xml)
                        else:
                            continue
                else:
                    print('All pathways for symlinks were printed for other folders')

                print('Resumed download was finished in %s' % self.execution_time(global_start))
                f.write('</batch>')
        except ValueError:
            print ('Folder %s was not found in folder tree' % folder_name)
            return
    
    def urns2rnef(self, pathway_urns:list, fout:str, xml_format:str='RNEF'):
        global_start = time.time()
        urn2pathway = self.load_containers()#works 40 sec - cache result to your application
        print('List of all pathways identifiers was retrieved in %s' % self.execution_time(global_start))
        
        missedURNs = list()
        download_counter = 0
        with open(fout, 'w',encoding='utf-8') as f:
            if xml_format == 'RNEF': f.write('<batch>')
            for urn in pathway_urns:
                try: 
                    pathwayId = urn2pathway[urn]['Id'][0]
                    start_time = time.time()
                    pathway_graph, pathway_xml = self.get_pathway(pathwayId,xml_format=xml_format,as_batch=False)
                    f.write(pathway_xml)
                    download_counter += 1
                    exec_time = self.execution_time(start_time)
                    print('%d out of %d pathways downloaded in %s, %d not found' %
                        (download_counter, len(pathway_urns), exec_time, len(missedURNs)))                
                except KeyError:
                    missedURNs.append(urn)
                    continue
            if xml_format == 'RNEF': f.write('</batch>')

        print('Pathways not found:\n%s' % '\n'.join(missedURNs))
        print("Total download time: %s" % self.execution_time(global_start))
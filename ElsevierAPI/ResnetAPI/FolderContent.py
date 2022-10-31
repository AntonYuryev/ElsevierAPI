from .ResnetAPISession import APISession, time, len, NO_REL_PROPERTIES
from .Resnet2rdf import ResnetRDF
from .ResnetGraph import ResnetGraph,PSPathway,PSObject,RESNET
from .rnef2sbgn import rnef2sbgn_str, minidom
import xml.etree.ElementTree as et
from urllib.parse import quote

PATHWAY_ID = 'Pathway ID'

class FolderContent (APISession): 
    pass
    def __init__(self, APIconfig:dict, preload_folder_tree=True,what2retrieve=NO_REL_PROPERTIES):
        """
        self.id2pathway = {id:PSObject} is loaded on demand 
        self.id2group = {id:PSObject} is loaded on demand
        self.id2result = {id:PSObject} is loaded on demand\n
        loading is done by self::get_objects_from_folders() or self::load_containers()
        """
        super().__init__(APIconfig,what2retrieve)
        self.PageSize = 1000
        self.DumpFiles = []
        self.id2folder = {}
        if preload_folder_tree:
            self.id2folder = self.load_folder_tree() # {folder_id:folder_zobj, folder_name:folder_zobj}
        #self.reference_cache_size = 1000000
          

    def get_folder_id(self,folder_name:str) -> list:
        if not hasattr(self, "id2folder"): 
            self.id2folder = self.load_folder_tree()

        for k,v in self.id2folder.items():
            if v['Name'] == folder_name:
                return k

        return None


    def get_subfolders(self, folder_ids: list):
        # FolderIds - list of folder IDs
        if not hasattr(self, "id2folder"): 
            self.id2folder = self.load_folder_tree()

        subfolders_ids = set()
        for folder_id in folder_ids:
            folder = self.id2folder[folder_id]
            if type(folder['SubFolders']) != type(None):
                subfolders_ids.update(folder['SubFolders']['long'])

        return subfolders_ids


    def get_subfolders_recursively(self, FolderId):
        """
        Returns
        -------
        {subfolders_ids+parent_folder_id}
        """
        accumulate_subfolder_ids = {FolderId}
        subfolder_ids = set(self.get_subfolders([FolderId]))
        while len(subfolder_ids) > 0:
            accumulate_subfolder_ids.update(subfolder_ids)
            subs = set()
            for db_id in subfolder_ids:
                subs.update(self.get_subfolders([db_id]))
            subfolder_ids = subs

        return accumulate_subfolder_ids


    def get_subfolder_tree(self,folder_name):
        """
        returns child2parent = {folder_id:folder_id}
        parent2child = {parent_id:[child_ids]}
        """
        if not hasattr(self, "id2folder"): 
            self.id2folder = self.load_folder_tree()

        folder_id = self.get_folder_id(folder_name)
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


    def load_containers(self, property_names=None, from_folders=[]):
        """
        retreives entire folder tree stricture from database 
        including all folder objects
        if from_folders is not specidied works for ~40sec

        Updates
        -----
        self.id2folder = {id:PSObject}
        self.id2pathway = {id:PSObject}
        self.id2group = {id:PSObject}
        self.id2result = {id:PSObject}

        Returns
        -------
        urn2pathway = {urn:PSObject}
        """
        if property_names is None: property_names = ['Name']
        print('Retrieving identifiers of all pathways from database may take couple minutes')

        if not hasattr(self,'id2folder'): self.id2folder = self.load_folder_tree()

        self.id2pathway = dict() # self.id2pathway = {id:PSObject}
        self.id2group = dict() # self.id2pathway = {id:PSObject}
        self.id2result = dict() # self.id2result = {id:PSObject}
        urn2pathway = dict()
        folder_objs =  [PSObject.from_zeep(folder) for folder in self.id2folder.values()]
        if from_folders:
            folder_objs = {o for o in folder_objs if o.name() in from_folders}

        for folder in folder_objs:
            zeep_objects = self.get_folder_objects_props(folder['Id'][0], property_names)
            ps_objects = self._zeep2psobj(zeep_objects)
            for Id, psObj in ps_objects.items():
                if psObj['ObjTypeName'][0] == 'Pathway':
                    try:
                        self.id2pathway[Id].update_with_value('Folders', folder['Name'])
                    except KeyError:
                        psObj['Folders'] = [folder['Name']]
                        self.id2pathway[Id] = psObj
                        urn2pathway[psObj['URN'][0]] = psObj
                elif psObj['ObjTypeName'][0] == 'Group':
                    try:
                        self.id2group[Id].update_with_value('Folders', folder['Name'])
                    except KeyError:
                            psObj['Folders'] = [folder['Name']]
                            self.id2group[Id] = psObj
                elif psObj['ObjTypeName'][0] == 'attributesSearch':
                    try:
                        self.id2result[Id].update_with_value('Folders', folder['Name'])
                    except KeyError:
                            psObj['Folders'] = [folder['Name']]
                            psObj['URN'] = ['urn:agi-pathway:'+quote(psObj['Name'][0])]
                            self.id2result[Id] = psObj

        print('Found %d pathways,%d groups, %d results in the database' % (len(self.id2pathway), len(self.id2group), len(self.id2result)))
        return urn2pathway


    def get_pathway(self, pathway_id,pathway_urn:str=None,pathway_name:str=None,
                    format='RNEF',put2folder='', add_props2rel=dict(), 
                    add_props2pathway=dict(), as_batch=True, prettify=True):
        """
        Input
        -----
        format = ['RNEF', 'SBGN', 'JSON-LD']
        add_rel_props = optional. {PropName:[PropValues]} - adds properties to all relations in pathway
        add_pathway_props = optional. {PropName:[PropValues]} - adds properties to all nodes in pathway
        put2folder - optional. specifies folder destination for pathway upon reimport into Pathway Studio
        'as_batch' must be True to output returned output string to single rnef file
        if 'as_batch'=False output strings can be concatenated into one file.  <batch> element must be added after concatenation
        prettify - uses minidom to prettify xml output  

        Returns
        -------
        tuple ResnetGraph, XML "format" string
        """

        pathway_graph = self.pathway_components([pathway_id],'id',self.relProps,self.entProps)
        
        if format == 'JSON-LD':
            return pathway_graph, ResnetRDF.fromResnetGraph(pathway_graph).to_jsons()
        else:
            return pathway_graph, self.pathway2xml(pathway_id,pathway_graph,pathway_urn,pathway_name,
                    format,put2folder, add_props2rel,add_props2pathway, as_batch, prettify)



    def pathway2xml(self, pathway_id,pathway_graph:ResnetGraph, pathway_urn:str=None,pathway_name:str=None,
                    format='RNEF',put2folder='', add_props2rel=dict(), 
                    add_props2pathway=dict(), as_batch=True, prettify=True):

        if hasattr(self,'id2pathway'):
            if not isinstance(pathway_urn,str):
                try:
                    pathway_urn = self.id2pathway[pathway_id]['URN'][0]
                    pathway_name = str(self.id2pathway[pathway_id]['Name'][0])
                except KeyError:
                    print('Pathway collection does not have %s pathway with URN %s' % (pathway_name,path_urn))

        if not isinstance(pathway_urn,str):
            print('Pathway has no URN specifed!!!!')
            path_urn = 'no_urn'
        
        if not isinstance(pathway_name,str):
            print('Pathway has no Name specifed!!!!')
            pathway_name = 'no_name'
        
        pathway_graph.load_references()
        graph_xml = self.to_rnef(pathway_graph,add_props2rel,add_props2pathway)

        import xml.etree.ElementTree as et
        rnef_xml = et.fromstring(graph_xml)
        rnef_xml.set('name', pathway_name)
        rnef_xml.set('urn', pathway_urn)
        rnef_xml.set('type', 'Pathway')

        lay_out = et.Element('attachments')
        attachments = self.get_layout(pathway_id)
        if attachments:
            lay_out.append(et.fromstring(attachments))
        rnef_xml.append(lay_out)
        
        batch_xml = et.Element('batch')
        batch_xml.insert(0,rnef_xml)
                   
        if format == 'SBGN':
            pathway_xml = et.tostring(batch_xml,encoding='utf-8',xml_declaration=True).decode("utf-8")
            pathway_xml = rnef2sbgn_str(pathway_xml)
        else:
            if put2folder:
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
            (pathway_name, pathway_graph.number_of_nodes(),pathway_graph.number_of_edges(),pathway_graph.size(weight="weight")))

        return str(pathway_xml)


    def get_group(self, group_id,group_urn:str=None,group_name:str=None, put2folder:str=None,as_batch=True, prettify=True):
        if hasattr(self,'id2group'):
            if not isinstance(group_urn,str):
                try:
                    group_urn = self.id2group[group_id]['URN'][0]
                    group_name = self.id2group[group_id]['Name'][0]
                except KeyError:
                    print('Pathway collection does not have %s group with URN %s' % (group_name,group_urn))

        if not isinstance(group_urn,str):
            print('Pathway has no URN specified!!!!')
            group_urn = 'no_urn'
        
        if not isinstance(group_name,str):
            print('Pathway has no Name specified!!!!')
            group_name = 'no_name'

        oql_query = 'SELECT Entity WHERE MemberOf (SELECT Group WHERE Id = {})'.format(group_id)
        group_graph = self.load_graph_from_oql(oql_query,entity_props=self.entProps,get_links=False)
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


    def get_result(self, result_id,put2folder:str=None,add_props2rel:dict=None, add_props2pathway:dict=None, as_batch=True, prettify=True):
        if not hasattr(self,'id2result'): return ResnetGraph(), ''
        try:
            result = self.id2result[result_id]
        except KeyError: 
            print('No results with %d id exists!' % result_id)

        result_name = result['Name'][0]
        result_graph = self.get_saved_results(result['Name'])
        result_graph.load_references()
        rnef_xml = et.fromstring(self.to_rnef(result_graph,add_props2rel,add_props2pathway))
        rnef_xml.set('name', result['Name'][0])
        rnef_xml.set('type', 'Pathway')
        pathway_urn = result['URN'][0]
        
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
            xml_node_pathway = et.SubElement(xml_nodes, 'node', {'local_id':pathway_local_id, 'urn':pathway_urn})
            et.SubElement(xml_node_pathway, 'attr', {'name': 'NodeType', 'value': 'Group'})
            xml_controls = et.SubElement(folder_resnet, 'controls')
            xml_control = et.SubElement(xml_controls, 'control', {'local_id':'CFE1'})
            et.SubElement(xml_control, 'attr', {'name':'ControlType', 'value':'MemberOf'})
            et.SubElement(xml_control, 'link', {'type':'in', 'ref':pathway_local_id})
            et.SubElement(xml_control, 'link', {'type':'out', 'ref':folder_local_id})
            batch_xml.append(folder_resnet)
        
        if as_batch:
            pathway_xml = et.tostring(batch_xml,encoding='utf-8',xml_declaration=True).decode("utf-8")
            if prettify: group_xml = minidom.parseString(pathway_xml).toprettyxml(indent='   ')
        else:
            pathway_xml = et.tostring(rnef_xml,encoding='utf-8',xml_declaration=True).decode("utf-8")
            if prettify: pathway_xml = self.pretty_xml(pathway_xml,no_declaration=True)
            #minidom does not work without xml_declaration

        print('\"%s\" search result downloaded as pathway: %d nodes, %d edges supported by %d references' % 
        (result_name, result_graph.number_of_nodes(),result_graph.number_of_edges(),result_graph.size(weight="weight")))
        return result_graph, str(pathway_xml)


    def get_objects_from_folders(self, FolderIds: list, property_names=None, with_layout=False):
        """
        Updates
        -----
        self.id2folder = {id:PSObject}
        self.id2pathway = {id:PSObject}
        self.id2group = {id:PSObject}
        self.id2result = {id:PSObject}

        Returns
        -------
        id2objects={id:PSObject}, PSObjects are annotated with 'Folders' attribute,\n
        where PSObject is either Pathway, Group, Result
        """
        if property_names is None: property_names = ['Name']
        if not hasattr(self,'id2folder'): self.id2folder = self.load_folder_tree()
        if not hasattr(self,'id2pathway'): self.id2pathway = dict() # self.id2pathway = {id:PSObject}
        if not hasattr(self,'id2group'): self.id2group = dict() # self.id2pathway = {id:PSObject}
        if not hasattr(self,'id2result'): self.id2result = dict() # self.id2result = {id:PSObject}
        id2objects = dict()
        for fid in FolderIds:
            folder_name = self.id2folder[fid]['Name']
            zeep_objects = self.get_folder_objects_props(fid, property_names)
            id2objs = self._zeep2psobj(zeep_objects)
            for Id, psObj in id2objs.items():
                if psObj['ObjTypeName'][0] == 'Pathway':
                    try:
                        self.id2pathway[Id].update_with_value('Folders', folder_name)
                    except KeyError:
                            psObj['Folders'] = [folder_name]
                            self.id2pathway[Id] = psObj
                    if with_layout:
                        psObj['layout'] = self.get_layout(Id)
                elif psObj['ObjTypeName'][0] == 'Group':
                    try:
                        self.id2group[Id].update_with_value('Folders', folder_name)
                    except KeyError:
                            psObj['Folders'] = [folder_name]
                            self.id2group[Id] = psObj
                elif psObj['ObjTypeName'][0] == 'attributesSearch':
                    try:
                        self.id2result[Id].update_with_value('Folders', folder_name)
                    except KeyError:
                            psObj['Folders'] = [folder_name]
                            psObj['URN'] = ['urn:agi-pathway:'+quote(psObj['Name'][0])]
                            self.id2result[Id] = psObj

            id2objects.update(id2objs)
        return id2objects


    def folder_content(self,folder_id_or_name:int or str,parent_folder_name:str,skip_id:set=None,skip_urn:set=None,
                            add_props2rel:dict=None,add_props2pathway:dict=None):
        """
        Input
        -----
        either folder_id or folder_name must be supplied
        if folder_id is supplied folder_name is retreived from database

        Updates
        -----
        self.id2folder = {id:PSObject}
        self.id2pathway = {id:PSObject}
        self.id2group = {id:PSObject}
        self.id2result = {id:PSObject}

        Dumps
        -----
        Objects from "folder_id_or_name" into 'folder_name' inside 'parent_folder_name' located in "self.data_dir"
        if size of dump file exceeds 100000000, "rnef_xml" is splitted into several RNEF files\n
        dump RNEF files are named as: 'content of folder_name#',
        where # - dump file number

        Returns
        -------
        tuple new_pathway_counter, symlinks
        symlinks = {pathway_id : pathway_obj['URN'][0]}
        """

        folder_id = self.get_folder_id(folder_id_or_name) if isinstance(folder_id_or_name, str) else folder_id_or_name
        folder_name = self.id2folder[folder_id]['Name']
        id2folder_obj = self.get_objects_from_folders([folder_id],self.entProps,with_layout=True)
        if id2folder_obj:
            print('Start downloading %d pathways from \"%s\" folder' % (len(id2folder_obj),folder_name))
        else:
            print('Folder \"%s\" has no pathways' % folder_name)
            return 0, {}
        
        folder_local_id = 'F0'
        new_pathway_counter = 0
        folder_object_counter = 0
        symlinks = dict()
        folder_download_start = time.time()
        write2folder = '' if parent_folder_name == folder_name else parent_folder_name
        for pathway_id, pathway_obj in id2folder_obj.items():
            if not pathway_obj['IsSymlink'][0]:
                #printing pathways only symlinks should be printed at the end and only if they were not downloaded for another folder
                start_time = time.time()
                folder_object_counter += 1
                if isinstance(skip_id,set):
                    if pathway_id in skip_id:
                        print('%s pathway was downloaded for another folder' % pathway_obj['Name'][0])
                        continue
                if isinstance(skip_urn,set):
                    if pathway_obj['URN'][0] in skip_urn:
                        print('%s pathway was downloaded for another folder' % pathway_obj['Name'][0])
                        continue

                if pathway_obj['ObjTypeName'][0] == 'Pathway':
                    pathway_graph, pathway_xml = self.get_pathway(pathway_id,add_props2rel=add_props2rel,add_props2pathway=add_props2pathway,as_batch=False)
                elif pathway_obj['ObjTypeName'][0] == 'Group':
                    pathway_graph, pathway_xml = self.get_group(pathway_id, as_batch=False)
                elif pathway_obj['ObjTypeName'][0] == 'attributesSearch':
                    pathway_graph, pathway_xml = self.get_result(pathway_id,add_props2rel=add_props2rel,add_props2pathway=add_props2pathway,as_batch=False)
                else:
                    print ('%s folder has object with unknown type %s: id = %d' % (folder_name,pathway_obj['ObjTypeName'][0],pathway_id))
                    continue

                self.dump_rnef(pathway_xml,folder_name,write2folder,can_close=False)
                new_pathway_counter += 1
                if isinstance(skip_id,set):skip_id.add(pathway_id)
                if isinstance(skip_urn,set):skip_urn.add(pathway_obj['URN'][0])
                
                print('%d out of %d objects in folder "%s" was downloaded in %s' %
                    (folder_object_counter, len(id2folder_obj),folder_name, self.execution_time(start_time)))
            else:
                symlinks[pathway_id] = pathway_obj['URN'][0]

            #printing folder object membership
            folder_resnet = et.Element('resnet')
            folder_nodes = et.SubElement(folder_resnet, 'nodes')
            member_controls = et.SubElement(folder_resnet, 'controls')

            xml_node_folder = et.SubElement(folder_nodes, 'node', {'local_id':folder_local_id, 'urn':'urn:agi-folder:'+str(folder_id)})
            et.SubElement(xml_node_folder, 'attr', {'name':'NodeType', 'value':'Folder'})
            et.SubElement(xml_node_folder, 'attr', {'name':'Name', 'value':folder_name})

            pathway_urn = pathway_obj['URN'][0]
            pathway_local_id = 'P'+str(folder_object_counter)
            folder_pathway_node = et.SubElement(folder_nodes, 'node', {'local_id':pathway_local_id, 'urn':pathway_urn})
            et.SubElement(folder_pathway_node, 'attr', {'name': 'NodeType', 'value': 'Pathway'})
            et.SubElement(folder_pathway_node, 'attr', {'name': 'Name', 'value': pathway_obj['Name'][0]})
            
            control_local_id = 'L'+str(folder_object_counter)
            member_control = et.SubElement(member_controls, 'control', {'local_id':control_local_id})
            et.SubElement(member_control, 'attr', {'name':'ControlType', 'value':'MemberOf'})
            if pathway_obj['IsSymlink'][0] == True:
                et.SubElement(member_control, 'attr', {'name':'Relationship', 'value':'symlink'})
            et.SubElement(member_control, 'link', {'type':'in', 'ref':pathway_local_id})
            et.SubElement(member_control, 'link', {'type':'out', 'ref':folder_local_id})

            folder_rnef = et.tostring(folder_resnet, encoding='utf-8',xml_declaration=False).decode("utf-8")
            folder_rnef = self.pretty_xml(folder_rnef,no_declaration=True)
            self.dump_rnef(folder_rnef,folder_name,write2folder)

        print('Total folder download time: %s' % self.execution_time(folder_download_start))
        return new_pathway_counter, symlinks


    def make_cache_dir(self, top_folder_id):
        """
        Creates
        -------
        directory structure in self.data_dir mirroring foler structur in database
        """
        subfolder_ids = self.get_subfolders([top_folder_id])
        subfolders_dict = {top_folder_id:subfolder_ids} if subfolder_ids else dict()

        while subfolders_dict:
            subs = dict()
            for parent_id, subfolder_ids in subfolders_dict.items():
                parent_folder_name = self.id2folder[parent_id]['Name']
                self.dump_path(parent_folder_name)
                for subfolder_id in subfolder_ids:
                    subsub_folder_ids = self.get_subfolders([subfolder_id])
                    if subsub_folder_ids:
                        subs.update({parent_id : self.get_subfolders([subfolder_id])})

                    subfolder_name = self.id2folder[subfolder_id]['Name']
                    self.dump_path(subfolder_name,parent_folder_name)
            subfolders_dict = subs


    def __subfolders_rnef(self, top_folder_name:str):
        """
        Returns
        -------
        child2parent = {folder_id:folder_id}

        Writes 
        ------
        folder tree rnef to self.data_dir/top_folder_name
        """
        child2parent, parent2child = self.get_subfolder_tree(top_folder_name)
        subtree_xml = str()

        for parent_id, child_ids in parent2child.items():
            folder_resnet = et.Element('resnet')
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

            subtree_xml = et.tostring(folder_resnet,encoding='utf-8',xml_declaration=False).decode("utf-8")
            subtree_xml = self.pretty_xml(subtree_xml,no_declaration=True)

            self.make_cache_dir(self.get_folder_id(top_folder_name))
            self.dump_rnef(subtree_xml,top_folder_name)

        return child2parent


    def content2rnef(self, top_folder_name:str,include_subfolders=True, add_props2rel:dict=None,add_pathway_props:dict=None):
        """
        Dumps
        -----
        content of top_folder_name into RNEF file 'content of top_folder_name.rnef' located in 'self.data_dir'
        """
        download_start_time = time.time()
        if not include_subfolders:
            return self.folder_content(None,top_folder_name,add_props2rel=add_props2rel,add_pathway_props=add_pathway_props)
        else:
            child2parent = self.__subfolders_rnef(top_folder_name)
            #fetching pathways by subfolders
            download_counter = 0
            folder_counter = 0
            printed_pathway_ids = set()
            symlinks_ids = set()
        
            for folder_id,parent_id in child2parent.items():
            # child2parent has {top_folder_id:top_folder_id} to process pathways in top_folder_name
                parent_folder_name = self.id2folder[parent_id]['Name']
                pathway_counter, symlinks = self.folder_content(folder_id,parent_folder_name,skip_id=printed_pathway_ids)
                symlinks_ids.update(symlinks.keys())
                download_counter += pathway_counter
                folder_counter +=1
                if pathway_counter > 0:
                    print('Downloaded %d pathways from %d out of %d folders in %s' % (download_counter,folder_counter,len(child2parent),self.execution_time(download_start_time)))
                    print('Relations cache has %d relations supported by %d references\n' % (self.Graph.number_of_edges(), self.Graph.weight()))
                if self.Graph.weight() > self.reference_cache_size:
                    print('Clearing cache due to size %d' % self.Graph.weight())
                    self.clear()

                symlinks2print = symlinks_ids.difference(printed_pathway_ids)
                if len(symlinks2print) > 0:
                    print('Will print %d pathways to support symlinks' % len(symlinks2print))
                    for pathway_id in symlinks2print:
                        if pathway_id in self.id2pathway.keys():
                            pathway_graph, pathway_xml = self.get_pathway(pathway_id, as_batch=False)
                            self.dump_rnef(pathway_xml,top_folder_name,parent_folder_name)
                        elif pathway_id in self.id2group.keys():
                            pathway_graph, group_xml = self.get_group(pathway_id, as_batch=False)
                            self.dump_rnef(group_xml,top_folder_name,parent_folder_name)
                        elif pathway_id in self.id2result.keys():
                            pathway_graph, result_xml = self.get_result(pathway_id, as_batch=False)
                            self.dump_rnef(result_xml,top_folder_name,parent_folder_name)
                        else:
                            continue
                else:
                    print('No symlinks to print for folder %s' % top_folder_name)
                self.close_rnef_dump(self.id2folder[folder_id]['Name'],parent_folder_name)

        print("Total execution time: %s" % self.execution_time(download_start_time))


    @staticmethod
    def __load_urns(rnef_file:str):
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


    def resume_download(self, top_folder_name:str, last_downloaded_folder:str):
        child2parent, parent2child = self.get_subfolder_tree(top_folder_name)
        subfolders = list(child2parent.keys())
        global_start = time.time()
        try:
            last_downloaded_folder_id = self.get_folder_id(last_downloaded_folder)
            start_folder_idx = subfolders.index(last_downloaded_folder_id)
            folder_counter = start_folder_idx+1
            rnef_file_to_continue =  self.name_output(top_folder_name)
            urns_printed = self.__load_urns(rnef_file_to_continue)

            download_counter = len(urns_printed)
            print ('Resuming download of %s folder from %s' % (top_folder_name,last_downloaded_folder))
            symlinks_dict = dict()
            with open(rnef_file_to_continue, "a", encoding='utf-8') as f:
                for i in range(start_folder_idx+1, len(subfolders)):
                    folder_id = subfolders[i]
                    folder_xml, pathway_counter, symlinks = self.folder_content(folder_id,as_batch=False,skip_urn=urns_printed)
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
                        if pathway_id in self.id2pathway.keys():
                            pathway_graph, pathway_xml = self.get_pathway(pathway_id, as_batch=False)
                            f.write(pathway_xml)
                        elif pathway_id in self.id2group.keys():
                            pathway_graph, pathway_xml = self.get_group(pathway_id, as_batch=False)
                            f.write(pathway_xml)
                        elif pathway_id in self.id2result.keys():
                            pathway_graph, pathway_xml = self.get_result(pathway_id, as_batch=False)
                            f.write(pathway_xml)
                        else:
                            continue
                else:
                    print('All pathways for symlinks were printed for other folders')

                print('Resumed download was finished in %s' % self.execution_time(global_start))
                f.write('</batch>')
        except ValueError:
            print ('Folder %s was not found in folder tree' % last_downloaded_folder)
            return
    
    
    def download_pathways_by_urn(self, pathway_urns:list, fout:str, format:str='RNEF',from_folders=[]):
        global_start = time.time()
        urn2pathway = self.load_containers(from_folders=from_folders) 
        # works 40 sec if "from_folders" is not specified
        print('List of all pathways identifiers was retrieved in %s' % self.execution_time(global_start))
        
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
                    pathway_id = urn2pathway[urn]['Id'][0]
                    start_time = time.time()
                    pathway_graph, pathway_str = self.get_pathway(pathway_id,format=format,as_batch=False)
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


    def find_pathways(self, for_entity_ids:list, in_folders:list):
        """
        Return
        -------
        {entity_id:PSObject},  where PSObject has 'Pathway ID' property

        Loads
        -----
        self.id2pathway = {pathway_id:PSObject}.
        where PSObject['ObjTypeName'] == 'Pathway' with 'Folders' property containing folder names
        """
        ent_ids = list(map(str,for_entity_ids))
        to_return = dict()
        for folder in in_folders:
            folder_id = self.get_folder_id(folder)
            sub_folder_ids = self.get_subfolders_recursively(folder_id)
            self.get_objects_from_folders(sub_folder_ids,self.entProps) # loads id2pathway
        
        for pathway_id in self.id2pathway.keys():
            id2psobj = self.get_pathway_members([pathway_id],None,ent_ids,['id'])
            for id, psobj in id2psobj.items():
                try:
                    to_return[id].update_with_value(PATHWAY_ID,pathway_id)
                except KeyError:
                    psobj.update_with_value(PATHWAY_ID,pathway_id)
                    to_return[id] = psobj

        return to_return


    def folder2pspathways(self,folder_id_or_name:int or str):
        """
        Input
        -----
        Either "folder_id" or "folder_name" must be supplied. 
        If "folder_id" is supplied "folder_name" is retreived from database

        Returns
        -------
        list of PSPathway objects from folder_id_or_name annotated with 'resnet' and 'Folders' properties
        """
        folder_id = self.get_folder_id(folder_id_or_name) if isinstance(folder_id_or_name, str) else folder_id_or_name
        folder_name = self.id2folder[folder_id]['Name']
        id2folder_obj = self.get_objects_from_folders([folder_id],self.entProps,with_layout=True)
        if id2folder_obj:
            print('Start downloading %d pathways from \"%s\" folder' % (len(id2folder_obj),folder_name))
        else:
            print('Folder \"%s\" has no pathways' % folder_name)
            return list()
        
        pspathways2return = list()
        for pathway_id, pathway_obj in id2folder_obj.items():
            if pathway_obj['ObjTypeName'][0] == 'Pathway':
                pathway_graph, pathway_xml = self.get_pathway(pathway_id,as_batch=False)
                ps_pathway = PSPathway(dict(pathway_obj),pathway_graph)
                ps_pathway[RESNET] = pathway_xml
                ps_pathway.update_with_value('Folders',folder_name)
                pspathways2return.append(ps_pathway)
                
        return pspathways2return

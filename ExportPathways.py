import sqlite3
import xml.etree.ElementTree as et
from ElsevierAPI import open_api_session
from xml.dom import minidom

def reverse_tree(tree):
    """This function reverses the tree (dict): {child: parent} => {parent: {child}}; {parent: {child}} => {child: parent}."""
    ret = {}
    for key in tree:
        if type(tree[key]) in [int, str]:
            if tree[key] not in ret:
                ret[tree[key]] = set()
            ret[tree[key]].add(key)
        elif type(tree[key]) == set:
            for item in tree[key]:
                if item in ret:
                    raise KeyError
                ret[item] = key
        else:
            raise TypeError
    return ret

def expand_rev(g, node, subtree=None):
    """This function traces paths to all leaves in a reversed tree {parent: {child}}."""
    if subtree == None:
        subtree = {}
    if node in g:
        subtree[node] = g[node]
        for child in g[node]:
            expand_rev(g, child, subtree)
    return subtree

def get_folder_tree(api):
    """This function downloads entire folder structure and returns dicts {id: name} and {child: parent}."""
    raw_tree = api.load_folder_tree()
    folder_names, folder_graph = {}, {}
    for item in raw_tree:
        if type(item) != int:
            continue
        folder_names[item] = raw_tree[item][0]['Name']
        if raw_tree[item][0]['SubFolders']:
            for subitem in raw_tree[item][0]['SubFolders']['long']:
                if subitem in folder_graph:
                    raise KeyError('Attempted assigning folder %d to parent %d when it is already assigned to another parent %d' % (subitem, item, folder_graph[subitem]))
                folder_graph[subitem] = item
    return folder_names, folder_graph

def get_folder_subtree(folder_graph, folder_names, folder_ids=None, folder_name=None):
    """This function takes in downloaded folder structure and folder names, and separates complete subtree containing specified {folder_ids} or `folder_name`."""
    assert folder_ids is not None or folder_name is not None, 'Neither folder ID nor folder name is specified.'
    folder_names_rev = reverse_tree(folder_names)
    if folder_ids is None:
        folder_ids = folder_names_rev[folder_name]
    _folder_graph = {child: {folder_graph[child]} for child in folder_graph}
    paths_to_root = {}
    for folder_id in folder_ids:
        paths_to_root.update(expand_rev(_folder_graph, folder_id, None)) 
    paths_to_root = {child: next(iter(paths_to_root[child])) for child in paths_to_root}
    subtree = {}
    subtree.update(paths_to_root)
    _folder_graph_rev = reverse_tree(folder_graph)
    paths_to_leaves_rev = {}
    for folder_id in folder_ids:
        paths_to_leaves_rev.update(expand_rev(_folder_graph_rev, folder_id))
    paths_to_leaves = reverse_tree(paths_to_leaves_rev)
    subtree.update(paths_to_leaves)
    subtree_names = {x: folder_names[x] for x in subtree}
    subtree_names.update({subtree[x]: folder_names[subtree[x]] for x in subtree})
    return subtree_names, subtree

def get_content_metadata_and_symlinks(api, folder_subtree, list_props):
    """This function downloads properties of objects found in specified folder subtree.
    Returns dict {item_id: {propname: propvalue}} and dict {item_id: {folder_id}} for symlinks."""
    flat_folders = set([key for key in folder_subtree])
    flat_folders = flat_folders.union(set([folder_subtree[key] for key in folder_subtree]))
    metadata = {}
    placement = {} # {folder_id: {(item_id, IsSymlink)}}
    for folder_id in flat_folders:
        raw_metadata = api.get_folder_objects_props(folder_id, list_props)
        if raw_metadata is None:
            continue
        assert raw_metadata['Completed'], 'Folder object properties came incomplete'
        object_refs = list(raw_metadata['Objects']['ObjectRef']) if 'Objects' in raw_metadata and 'ObjectRef' in raw_metadata['Objects'] else []
        for object_ref in object_refs:
            if object_ref['ObjTypeName'] not in ['Pathway', 'Group']:
                continue
            if object_ref['Id'] not in metadata:
                metadata[object_ref['Id']] = {
                    'URN': object_ref['URN'],
                    'ObjTypeName': object_ref['ObjTypeName']
                }
            if folder_id not in placement:
                placement[folder_id] = set()
            placement[folder_id].add((object_ref['Id'], object_ref['IsSymlink']))
        properties = list(raw_metadata['Properties']['ObjectProperty']) if 'Properties' in raw_metadata and 'ObjectProperty' in raw_metadata['Properties'] else []
        for property in properties:
            if 'ObjId' not in property or 'PropName' not in property:
                continue
            prop_name = property['PropName']
            if property['ObjId'] not in metadata:
                continue
            if prop_name not in list_props:
                continue
            if prop_name not in metadata[property['ObjId']]:
                metadata[property['ObjId']][prop_name] = set()
            values = set([property['PropValues'][key] for key in [key for key in property['PropValues']]][0])
            metadata[property['ObjId']][prop_name] = metadata[property['ObjId']][prop_name].union(values)
    return metadata, placement

def flat_pathways_to_rnef_storage(apiconfig_filename, metadata, db_destination):
    """This method downloads content of pathways and groups by looking at IDs in `metadata`,
    conbines downloaded content with properties stored in `metadata`,
    and submits (id, xml) to SQLite database."""
    conn = sqlite3.connect(db_destination)
    cur = conn.cursor()
    cur.execute('create table if not exists content (id number, rnef text);')
    conn.commit()
    for key in metadata:
        cur.execute('select id from content where id = %d;' % (key))
        if cur.fetchone() is not None:
            continue
        if metadata[key]['ObjTypeName'] == 'Pathway':
            # create new API object each time, otherwise memory leaks somewhere in API degrading performance
            api = open_api_session(apiconfig_filename)
            print('Retrieving pathway "%s" (id=%d, urn="%s")' % (next(iter(metadata[key]['Name'])), key, metadata[key]['URN']))
            _, raw_xml = api.get_pathway(pathwayId=key, path_urn=metadata[key]['URN'], path_name=next(iter(metadata[key]['Name'])), as_batch=False, prettify=False)
            del api
        elif metadata[key]['ObjTypeName'] == 'Group':
            api = open_api_session(apiconfig_filename)
            print('Retrieving group "%s" (id=%d, urn="%s")' % (next(iter(metadata[key]['Name'])), key, metadata[key]['URN']))
            _, raw_xml = api.get_group(group_id=key, group_urn=metadata[key]['URN'], group_name=next(iter(metadata[key]['Name'])), as_batch=False, prettify=False)
            del api
        x = et.fromstring(raw_xml)
        x.set('type', metadata[key]['ObjTypeName'])
        props = et.Element('properties')
        for prop_name in metadata[key]:
            if prop_name in ['URN', 'Name', 'ObjTypeName']:
                continue
            for prop_value in metadata[key][prop_name]:
                attr = et.Element('attr')
                attr.set('name', prop_name)
                attr.set('value', prop_value)
                props.append(attr)
        x.append(props)
        xml_str = minidom.parseString(et.tostring(x, encoding='utf8').decode('utf8')).toprettyxml(indent='  ')
        xml_str = xml_str[xml_str.find('\n')+1:]
        #cur.execute('delete from content where id = %d;' % (key))
        cur.execute('insert into content (id, rnef) select ?, ?;', (key, xml_str))
        conn.commit()

PS_API = open_api_session('./ElsevierAPI/.misc/APIconfig.json')
folder_names, folder_graph = get_folder_tree(PS_API)
#this_names, this_subtree = get_folder_subtree(folder_graph, folder_names, folder_name='Integrin Receptors')
#this_metadata, this_placement = get_content_metadata_and_symlinks(PS_API, this_subtree, ['Name', 'Notes', 'Description', 'PMID', 'CellType', 'Organ', 'Tissue', 'PathwayType', 'Organ System'])
this_metadata, this_placement = get_content_metadata_and_symlinks(PS_API, folder_graph, ['Name', 'Notes', 'Description', 'PMID', 'CellType', 'Organ', 'Tissue', 'PathwayType', 'Organ System'])
del PS_API
flat_pathways_to_rnef_storage('./ElsevierAPI/.misc/APIconfig.json', this_metadata, '../pathways-pse/rnef.db')

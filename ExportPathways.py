import argparse
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

def flat_pathways_to_rnef_storage(apiconfig_filename, metadata, db_destination, bypass=False):
    """This method downloads content of pathways and groups by looking at IDs in `metadata`,
    conbines downloaded content with properties stored in `metadata`,
    and submits (id, xml) to SQLite database."""
    conn = sqlite3.connect(db_destination)
    cur = conn.cursor()
    cur.execute('create table if not exists content (id number, rnef text);')
    conn.commit()
    if not bypass:
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

def folders_to_rnef_storage(subtree, names, metadata, placement, db_destination, table_name, bypass=False):
    """This method generates RNEF <resnet> elements with folder/folder and folder/object placements,
    and submits them to SQLite database.
    """
    conn = sqlite3.connect(db_destination)
    cur = conn.cursor()
    cur.execute('drop table if exists %s;' % (table_name))
    cur.execute('drop table if exists %s_h;' % (table_name))
    cur.execute('drop table if exists %s_f;' % (table_name))
    cur.execute('create table if not exists %s (id number, rnef text);' % (table_name))
    cur.execute('create table if not exists %s_h (folderid number, objectid number);' % (table_name))
    cur.execute('create table if not exists %s_f (parentid number, childid number);' % (table_name))
    conn.commit()
    if not bypass:
        for folder_id in placement:
            for object_id in placement[folder_id]:
                cur.execute('insert into %s_h (folderid, objectid) select ?, ?;' % (table_name), (folder_id, object_id[0]))
            conn.commit()
        hierarchy = {}
        subtree_rev = reverse_tree(subtree)
        for parent_folder_id in subtree_rev:
            urn_parent_folder = 'urn:agi-folder:%d' % (parent_folder_id)
            name_parent_folder = names[parent_folder_id]
            type_parent_folder = 'Folder'
            local_id_parent_folder = 'F%d' % (parent_folder_id)
            if parent_folder_id not in hierarchy:
                hierarchy[parent_folder_id] = []
            for child_folder_id in subtree_rev[parent_folder_id]:
                urn_child_folder = 'urn:agi-folder:%d' % (child_folder_id)
                name_child_folder = names[child_folder_id]
                type_child_folder = 'Folder'
                local_id_child_folder = 'F%d' % (child_folder_id)
                hierarchy[parent_folder_id].append((local_id_parent_folder, urn_parent_folder, name_parent_folder, type_parent_folder, local_id_child_folder, urn_child_folder, name_child_folder, type_child_folder, False))
                cur.execute('insert into %s_f (parentid, childid) select ?, ?;' % (table_name), (parent_folder_id, child_folder_id))
            conn.commit()
            if parent_folder_id in placement:
                for item in placement[parent_folder_id]:
                    urn_item = metadata[item[0]]['URN']
                    name_item = next(iter(metadata[item[0]]['Name']))
                    type_item = metadata[item[0]]['ObjTypeName']
                    symlink_item = item[1]
                    local_id_item = 'M%d' % (item[0])
                    hierarchy[parent_folder_id].append((local_id_parent_folder, urn_parent_folder, name_parent_folder, type_parent_folder, local_id_item, urn_item, name_item, type_item, symlink_item))
        for child_folder_id in subtree:
            urn_child_folder = 'urn:agi-folder:%d' % (child_folder_id)
            name_child_folder = names[child_folder_id]
            type_child_folder = 'Folder'
            local_id_child_folder = 'F%d' % (child_folder_id)
            if child_folder_id in placement:
                if child_folder_id not in hierarchy:
                    hierarchy[child_folder_id] = []
                for item in placement[child_folder_id]:
                    urn_item = metadata[item[0]]['URN']
                    name_item = next(iter(metadata[item[0]]['Name']))
                    type_item = metadata[item[0]]['ObjTypeName']
                    symlink_item = item[1]
                    local_id_item = 'M%d' % (item[0])
                    hierarchy[child_folder_id].append((local_id_child_folder, urn_child_folder, name_child_folder, type_child_folder, local_id_item, urn_item, name_item, type_item, symlink_item))
        for parent_id in hierarchy:
            resnet = et.Element('resnet')
            nodes = et.Element('nodes')
            controls = et.Element('controls')
            parent = et.Element('node')
            nodetype = et.Element('attr')
            name = et.Element('attr')
            parent.set('local_id', hierarchy[parent_id][0][0])
            parent.set('urn', hierarchy[parent_id][0][1])
            nodetype.set('name', 'NodeType')
            nodetype.set('value', hierarchy[parent_id][0][3])
            name.set('name', 'Name')
            name.set('value', hierarchy[parent_id][0][2])
            parent.append(nodetype)
            parent.append(name)
            nodes.append(parent)
            for i in range(len(hierarchy[parent_id])):
                child = et.Element('node')
                nodetype = et.Element('attr')
                name = et.Element('attr')
                child.set('local_id', hierarchy[parent_id][i][4])
                child.set('urn', hierarchy[parent_id][i][5])
                nodetype.set('name', 'NodeType')
                nodetype.set('value', hierarchy[parent_id][i][7])
                name.set('name', 'Name')
                name.set('value', hierarchy[parent_id][i][6])
                child.append(nodetype)
                child.append(name)
                nodes.append(child)
                control = et.Element('control')
                controltype = et.Element('attr')
                link_in = et.Element('link')
                link_out = et.Element('link')
                control.set('local_id', 'CSF%d' % (i+1))
                controltype.set('name', 'ControlType')
                controltype.set('value', 'MemberOf')
                link_in.set('type', 'in')
                link_in.set('ref', hierarchy[parent_id][i][4])
                link_out.set('type', 'out')
                link_out.set('ref', hierarchy[parent_id][i][0])
                control.append(controltype)
                if hierarchy[parent_id][i][8]:
                    relationship = et.Element('attr')
                    relationship.set('name', 'Relationship')
                    relationship.set('value', 'symlink')
                    control.append(relationship)
                control.append(link_in)
                control.append(link_out)
                controls.append(control)
            resnet.append(nodes)
            resnet.append(controls)
            xml_str = minidom.parseString(et.tostring(resnet, encoding='utf8').decode('utf8')).toprettyxml(indent='  ')
            xml_str = xml_str[xml_str.find('\n')+1:]
            cur.execute('insert into %s (id, rnef) select ?, ?;' % (table_name), (parent_id, xml_str))
            conn.commit()

def write_rnef(rnef_storage, placement_table_name, output_filepath, bypass=False):
    """This function writes content of intermediate SQLite database into RNEF file."""
    if not bypass:
        conn = sqlite3.connect(rnef_storage)
        cur = conn.cursor()
        conn2 = sqlite3.connect(rnef_storage)
        cur2 = conn2.cursor()
        ids = set()
        rows = cur.execute('select distinct objectid from %s_h;' % (placement_table_name)).fetchall()
        for row in rows:
            ids.add(row[0])
        with open(output_filepath, mode='w', encoding='utf8') as o:
            o.write('<batch>\n')
            with open(output_filepath, mode='r', encoding='utf8') as i:
                schema = ''
                for line in i:
                    schema += line
                o.write(schema)
                o.write('\n')
                rows = cur.execute('select rowid, id from content order by rowid asc;').fetchall()
                total = 0
                for row in rows:
                    rowid = row[0]
                    itemid = row[1]
                    if itemid not in ids:
                        continue
                    rnef = cur2.execute('select rnef from content where rowid = ?;', (rowid, )).fetchone()[0]
                    rnef = rnef.replace(' name="Sentence"', ' name="msrc"')
                    o.write(rnef)
                    total += 1
                    print('Wrote %d pathways and groups' % (total))
                rows = cur.execute('select rnef from %s order by rowid asc;' % (placement_table_name)).fetchall()
                total = 0
                for row in rows:
                    rnef = row[0]
                    o.write(rnef)
                    total += 1
                    print('Wrote %d folder placement RNEFs' % (total))
                o.write('</batch>\n')

def do_the_job(api_config_json, export_pathways_content=True, export_folder_structure=True, write_rnef_file=True, top_folder_name=None, placement_table_name='placement', output_filepath='out.rnef', rnef_storage='rnef.db'):
    PS_API = open_api_session(api_config_json)    
    _names, _graph = get_folder_tree(PS_API)
    folder_names, subtree = _names, _graph
    if top_folder_name:
        folder_names, subtree = get_folder_subtree(_graph, _names, folder_name=top_folder_name)
    metadata, placement = get_content_metadata_and_symlinks(PS_API, subtree, ['Name', 'Notes', 'Description', 'PMID', 'CellType', 'Organ', 'Tissue', 'PathwayType', 'Organ System'])
    del PS_API
    flat_pathways_to_rnef_storage(apiconfig_filename=api_config_json, metadata=metadata, db_destination=rnef_storage, bypass=not export_pathways_content)
    folders_to_rnef_storage(subtree=subtree, names=folder_names, metadata=metadata, placement=placement, db_destination=rnef_storage, table_name=placement_table_name, bypass=not export_folder_structure)
    write_rnef(rnef_storage=rnef_storage, placement_table_name=placement_table_name, output_filepath=output_filepath, bypass=not write_rnef_file)

if __name__ == '__main__':
    """
    Test:
    python ExportPathways.py --config .\ElsevierAPI\.misc\APIconfig.json --folder-name "C-Type Lectin/Selectin Receptors/CD Molecules"
    results in something like:
    > Pathways
      > Signal Processing
        > Receptor Signaling
          > C-Type Lectin/Selectin Receptors/CD Molecules
            - SELE -> ELK-SRF Signaling
            - Sialophorin -> CTNNB/MYC/TP53 Signaling
    """
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', type=str, required=True, help='Path and name of config file')
    ap.add_argument('--folder-name', type=str, help='Top folder name (will download everything if folder name is not specified)')
    ap.add_argument('--storage', type=str, default='rnef.db', help='Location for intermediate local data storage')
    ap.add_argument('--output', type=str, default='out.rnef', help='Location of RNEF file to write')
    ap.add_argument('--placement-table-name', type=str, default='placement', help='Table name for folder placement data in intermediate local data storage')
    ap.add_argument('--skip-pathways', type=bool, nargs='?', const=True, default=False, help='Do not download pathways and groups')
    ap.add_argument('--skip-folders', type=bool, nargs='?', const=True, default=False, help='Do not download folder structure')
    ap.add_argument('--skip-write-rnef', type=bool, nargs='?', const=True, default=False, help='Do not write RNEF file at the end')
    args = ap.parse_args()
    api_config_json, export_pathways_content, export_folder_structure, write_rnef_file, top_folder_name, placement_table_name, output_filepath, rnef_storage = args.config, not args.skip_pathways, not args.skip_folders, not args.skip_write_rnef, args.folder_name, args.placement_table_name, args.output, args.storage
    do_the_job(api_config_json=api_config_json, export_pathways_content=export_pathways_content, export_folder_structure=export_folder_structure, write_rnef_file=write_rnef_file, top_folder_name=top_folder_name, placement_table_name=args.placement_table_name, output_filepath=output_filepath, rnef_storage=rnef_storage)

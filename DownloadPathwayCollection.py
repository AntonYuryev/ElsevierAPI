import time
import argparse
import textwrap
from ElsevierAPI import open_api_session
import xml.etree.ElementTree as et
from xml.dom import minidom

ps_api = open_api_session()

if __name__ == "__main__":
    instructions = '''
    folder - name of the folder in Resnet with pathway collection
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog=textwrap.dedent(instructions))
    parser.add_argument('-f', '--folder', type=str, required=True)
    args = parser.parse_args()

    global_start = time.time()
    print('Attempting to download pathways from "%s"' % args.folder)

    start_time = time.time()
    child2parent, parent2child = ps_api.get_subfolder_tree(args.folder)

    printed_pathway_ids = set()

    with open('pathways from '+args.folder+'.rnef', "w", encoding='utf-8') as f:
        f.write('<?xml version="1.0" ?>\n')
        f.write('<batch>\n')
        #fetching pathways
        download_counter = 0
        for folder_id in child2parent.keys():
            pathway_counter = 0
            folder_name = ps_api.id2folders[folder_id][0]['Name']
            print('Downloading pathways from %s folder' % folder_name)
            folder_resnet = et.Element('resnet')
            folder_nodes = et.SubElement(folder_resnet, 'nodes')
            xml_controls = et.SubElement(folder_resnet, 'controls')
            folder_local_id = 'F0'
            xml_node_folder = et.SubElement(folder_nodes, 'node', {'local_id':folder_local_id, 'urn': 'urn:agi-folder:'+str(folder_id)})
            et.SubElement(xml_node_folder, 'attr', {'name': 'NodeType', 'value': 'Folder'})
            et.SubElement(xml_node_folder, 'attr', {'name': 'Name', 'value': folder_name})

            id2pathways = ps_api.get_objects_from_folders([folder_id],with_layout=True)

            for pathway_id, pathway in id2pathways.items():
                start_time = time.time()

                if pathway_id not in printed_pathway_ids:
                    if pathway['ObjTypeName'][0] == 'Pathway':
                        pathway_graph, pathway_xml = ps_api.get_pathway(pathway_id,as_batch=False)
                    elif pathway['ObjTypeName'][0] == 'Group':
                        pathway_graph, pathway_xml = ps_api.get_group(pathway_id,as_batch=False)
                    else:
                        continue

                    f.write(pathway_xml) #printing pathways 
                    printed_pathway_ids.add(pathway_id)

                pathway_counter += 1
                print('%d out of %d pathways in folder %s downloaded in %s' %
                    (pathway_counter, len(id2pathways),folder_name, ps_api.execution_time(start_time)))

                pathway_urn = pathway['URN'][0]
                pathway_local_id = 'P'+str(pathway_counter)
                xml_node_pathway = et.SubElement(folder_nodes, 'node', {'local_id':pathway_local_id, 'urn':pathway_urn})
                et.SubElement(xml_node_pathway, 'attr', {'name': 'NodeType', 'value': 'Pathway'})
                
                control_local_id = 'L'+str(pathway_counter)
                xml_control = et.SubElement(xml_controls, 'control', {'local_id':control_local_id})
                et.SubElement(xml_control, 'attr', {'name':'ControlType', 'value':'MemberOf'})
                if pathway['IsSymlink'][0] == True:
                    et.SubElement(xml_control, 'attr', {'name':'Relationship', 'value':'symlink'})
                et.SubElement(xml_control, 'link', {'type':'in', 'ref':pathway_local_id})
                et.SubElement(xml_control, 'link', {'type':'out', 'ref':folder_local_id})

            if len(id2pathways) > 0: #printing folder with pathways
                folder_xml = et.tostring(folder_resnet,encoding='utf-8',xml_declaration=False).decode("utf-8")
                folder_xml = ps_api.pretty_xml(folder_xml,no_declaration=True)
                f.write(folder_xml)

            download_counter += pathway_counter
            ps_api.Graph.clear() #release memory during big dumps

        print('Downloaded %d pathways' % download_counter)

        #printing folder tree
        for parent_id, child_ids in parent2child.items():
            folder_name = ps_api.id2folders[parent_id][0]['Name']
            folder_resnet = et.Element('resnet')
            folder_nodes = et.SubElement(folder_resnet, 'nodes')
            xml_controls = et.SubElement(folder_resnet, 'controls')
            parent_local_id = 'F0'
            xml_node_folder = et.SubElement(folder_nodes, 'node', {'local_id':parent_local_id, 'urn':'urn:agi-folder:'+str(parent_id)})
            et.SubElement(xml_node_folder, 'attr', {'name':'NodeType', 'value':'Folder'})
            et.SubElement(xml_node_folder, 'attr', {'name':'Name', 'value':folder_name})
            for child_id in child_ids:
                child_local_id = str(child_id)
                xml_node_pathway = et.SubElement(folder_nodes, 'node', {'local_id':child_local_id, 'urn':'urn:agi-folder:'+child_local_id})
                et.SubElement(xml_node_pathway, 'attr', {'name':'NodeType', 'value':'Folder'})
                
                control_local_id = 'L'+str(child_id)
                xml_control = et.SubElement(xml_controls, 'control', {'local_id':control_local_id})
                et.SubElement(xml_control, 'attr', {'name':'ControlType', 'value':'MemberOf'})
                et.SubElement(xml_control, 'link', {'type':'in', 'ref':child_local_id})
                et.SubElement(xml_control, 'link', {'type':'out', 'ref':parent_local_id})


        tree_xml = et.tostring(folder_resnet,encoding='utf-8',xml_declaration=False).decode("utf-8")
        tree_xml = ps_api.pretty_xml(tree_xml,no_declaration=True)
        f.write(tree_xml)
        f.write('</batch>')
    print("Total execution time: %s" % ps_api.execution_time(global_start))



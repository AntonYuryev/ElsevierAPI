import time
import argparse
import textwrap
import networkx as nx
from ElsevierAPI import APIconfig 
from ElsevierAPI.ResnetAPI.NetworkxObjects import REF_PROPS,REF_ID_TYPES
from ElsevierAPI.ResnetAPI.ResnetAPISession import APISession

ps_api = APISession(APIconfig['ResnetURL'], APIconfig['PSuserName'], APIconfig['PSpassword'])

def DownloadPathwayXML(PathwayURN:str,URNtoPathway:dict,get_sbgn=False,out_dir=''):
    try: pathwayId = URNtoPathway[PathwayURN]['Id'][0]
    except KeyError: return False

    layout = ps_api.get_layout(pathwayId)
    pathwayName = URNtoPathway[PathwayURN]['Name'][0]

    pathway_graph = ps_api.get_pathway_components([pathwayId],'id',retrieve_rel_properties=REF_PROPS + REF_ID_TYPES)
    ps_api.Graph.count_references(pathway_graph)

    import xml.etree.ElementTree as et
    rnef_xml = et.fromstring(ps_api.to_rnef(pathway_graph))
    rnef_xml.set('name', pathwayName)
    rnef_xml.set('urn', PathwayURN)

    lay_out = et.Element('attachments')
    lay_out.insert(0,et.fromstring(layout))
    rnef_xml.insert(2,lay_out)

    batch_xml = et.Element('batch')
    batch_xml.insert(0,rnef_xml)
    
    from xml.dom import minidom
    pathway_rnef = et.tostring(batch_xml,encoding='utf-8',xml_declaration=True).decode("utf-8")
    pathway_rnef = minidom.parseString(pathway_rnef).toprettyxml(indent='   ')
    from  ElsevierAPI.ResnetAPI.rnef2sbgn import make_file_name
    fout_name = out_dir +'/' + make_file_name(pathwayName)
    with open(fout_name+'.rnef', mode='w', encoding='utf-8') as f2: f2.write(pathway_rnef)
    
    if get_sbgn == True:
        from ElsevierAPI.ResnetAPI.rnef2sbgn import rnef2sbgn_str
        pathway_sbgn = rnef2sbgn_str(pathway_rnef, classmapfile='ElsevierAPI/ResnetAPI/rnef2sbgn_map.xml')
        with open(fout_name+'.sbgn.xml', mode='w', encoding='utf-8') as f2: f2.write(pathway_sbgn)

    print('\"%s\" pathway downloaded: %d nodes, %d edges supported by %d references' % 
        (pathwayName, pathway_graph.number_of_nodes(),pathway_graph.number_of_edges(),pathway_graph.size()))
    return True

if __name__ == "__main__":
    instructions = '''
    infile - single column file with pathway URNs
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog=textwrap.dedent(instructions))
    parser.add_argument('-i', '--infile', type=str, required=True)
    args = parser.parse_args()

    global_start = time.time()

    import csv
    f = open(args.infile,"r")
    PathwayURNs = csv.reader(f, delimiter="\t")
    urnSet = [u[0] for u in PathwayURNs]

    print('Attempting to download %s pathways from %s' %(len(urnSet),args.infile))

    out_dir = str(args.infile)[:str(args.infile).rfind('/')]
    start_time = time.time()
    IDtoPathway,URNtoPathway = ps_api.get_all_pathways()#works 40 sec - cache result to your application
    print('List of all pathways identifiers was retrieved in %s' % ps_api.execution_time(start_time))
    missedURNs = list()
    download_counter = 0
    for urn in urnSet:
        start_time = time.time()
        if DownloadPathwayXML(urn,URNtoPathway,get_sbgn=True,out_dir=out_dir):
            download_counter += 1
            exec_time = ps_api.execution_time(start_time)
            print('%d out of %d pathways downloaded in %s, %d not found' %
                  (download_counter, len(urnSet), exec_time, len(missedURNs)))
        else:
            missedURNs.append(urn)

    print('Pathways not found:\n%s' % '\n'.join(missedURNs))
    print("Total execution time: %s" % ps_api.execution_time(global_start))



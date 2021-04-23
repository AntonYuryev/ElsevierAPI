import time
import networkx as nx
import argparse
import textwrap
from ElsevierAPI import networx as PSnx,ExecutionTime
from ElsevierAPI.ResnetAPI.NetworkxObjects import REF_PROPS,REF_ID_TYPES

def DownloadPathwayXML(PathwayURN:str,URNtoPathway:dict,get_sbgn=False,out_dir=''):
    try: pathwayId = URNtoPathway[PathwayURN]['Id'][0]
    except KeyError: return False

    layout = PSnx.get_layout(pathwayId)
    pathwayName = URNtoPathway[PathwayURN]['Name'][0]
    
    OQLquery = 'SELECT Entity WHERE MemberOf (SELECT Network WHERE id = {pId})'
    pathway_nodes = PSnx.LoadGraphFromOQL(OQLquery.format(pId=pathwayId),getLinks=False)

    OQLquery = 'SELECT Relation WHERE MemberOf (SELECT Network WHERE id = {pId})'
    pathway_relations = PSnx.LoadGraphFromOQL(OQLquery.format(pId=pathwayId), REF_PROPS+REF_ID_TYPES)
    pathway_graph = nx.compose(pathway_relations, pathway_nodes)
    PSnx.count_references(pathway_graph)

    import xml.etree.ElementTree as et
    rnef_xml = et.fromstring(PSnx.toRNEF(pathway_graph))
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
    from  ElsevierAPI.ResnetAPI.rnef2sbgn import makeFileName
    fout_name = out_dir+'/'+makeFileName(pathwayName)
    with open(fout_name+'.rnef', mode='w', encoding='utf-8') as f2: f2.write(pathway_rnef)
    
    if get_sbgn == True:
        from ElsevierAPI.ResnetAPI.rnef2sbgn import rnef2sbgnStr
        pathway_sbgn = rnef2sbgnStr(pathway_rnef,'ElsevierAPI/ResnetAPI/rnef2sbgn_map.xml.xml')
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
    urnList = set()
    for urn in PathwayURNs:
        urnList.add(urn[0])

    print('Will attempt downloading %s pathways from %s' %(len(urnList),args.infile))
    out_dir = str(args.infile)[:str(args.infile).rfind('/')]

    IDtoPathway,URNtoPathway = PSnx.GetAllPathways()
    get_counter = 0
    miss_counter = 0
    missedURNs = list()
    for urn in urnList:
        start_time = time.time()
        if DownloadPathwayXML(urn,URNtoPathway,get_sbgn=True,out_dir=out_dir) == True:
            get_counter += 1
            print('%d out of %d pathways downloaded in %s, %d not found' % 
                (get_counter,len(urnList),ExecutionTime(start_time),miss_counter))
        else:
            miss_counter +=1
            missedURNs.append(urn)

    print('Pathways not found:\n%s' % '\n'.join(missedURNs))
    print("Total execution time: %s" % ExecutionTime(global_start))



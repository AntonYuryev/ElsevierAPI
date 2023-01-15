from .NetworkxObjects import PSObject
from .ResnetGraph import ResnetGraph,RESNET
import xml.etree.ElementTree as et
from .rnef2sbgn import rnef2sbgn_str,minidom

################################## PSPathway, PSPathway, PSPathway, PSPathway##########################################
class PSPathway(PSObject):
    pass
    def __init__(self,dic=dict(),g=ResnetGraph()):
        super().__init__(dic)
        self.graph = ResnetGraph()
        
        self.update(dic)
        self.graph.update(g)
        self.graph.urn2obj.update(g.urn2obj)
        self.graph.urn2rel.update(g.urn2rel)


    @classmethod
    def from_resnet(cls, resnet:et.Element, add_annotation=dict()):
        pathway = PSPathway()
        pathway['Name'] = [resnet.get('name')]
        pathway['URN'] = [resnet.get('urn')]
        pathway.update(add_annotation)
        [pathway.update_with_value(p.get('name'),p.get('value')) for p in resnet.findall('properties/attr')]

        if pathway.graph._parse_nodes_controls(resnet):
            pathway[RESNET] = et.tostring(resnet, encoding='unicode', method='xml')
            return pathway
        else:
            print('Invalid <resnet> section')
            return None


    @classmethod
    def from_pathway(cls,other:"PSPathway"):
        return PSPathway(other,other.graph)


    def number_of_nodes(self, obj_type=''):
        if obj_type:
            try:
                return self['#'+obj_type]
            except KeyError:
                obj_count = len([o for o in self.graph.urn2obj.values() if obj_type in o['ObjTypeName']])
                self['#'+obj_type] = obj_count
                return obj_count
        else:
            return len(self.graph.urn2obj)


    def merge_pathway(self, other:"PSPathway"):
        self.merge_obj(other)
        self.graph.add_graph(other.graph)
        self.graph.urn2obj.update(other.graph.urn2obj)
        self.graph.urn2rel.update(other.graph.urn2rel)


    def to_rnef(self, format='RNEF',add_props2rel=dict(),add_props2pathway=dict(), as_batch=True, prettify=True):
        '''
        Returns
        -------
        graph in RNEF XML string
        '''
        pathway_urn = self.urn()
        if not pathway_urn:
            print('Pathway has no URN!!!!')
            pathway_urn = 'no_urn'
        
        pathway_name = self.name()
        if not pathway_name:
            print('Pathway has no Name!!!!')
            pathway_name = 'no_name'
        
        self.graph.load_references()
        graph_xml = self.to_rnef(add_props2rel,add_props2pathway)

        import xml.etree.ElementTree as et
        rnef_xml = et.fromstring(graph_xml)
        rnef_xml.set('name', pathway_name)
        rnef_xml.set('urn', pathway_urn)
        rnef_xml.set('type', self.objtype())

        if self.objtype() == 'Pathway':
            lay_out = et.Element('attachments')
            try:
                lay_out.append(et.fromstring(self['layout']))
            except KeyError: pass
            rnef_xml.append(lay_out)
        
        batch_xml = et.Element('batch')
        batch_xml.insert(0,rnef_xml)
                   
        if format == 'SBGN':
            pathway_xml = et.tostring(batch_xml,encoding='utf-8',xml_declaration=True).decode("utf-8")
            pathway_xml = rnef2sbgn_str(pathway_xml)
        else:
            '''
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
            '''
            if as_batch:
                pathway_xml = et.tostring(batch_xml,encoding='utf-8',xml_declaration=True).decode("utf-8")
                if prettify: pathway_xml = minidom.parseString(pathway_xml).toprettyxml(indent='   ')
            else:
                pathway_xml = et.tostring(rnef_xml,encoding='utf-8',xml_declaration=True).decode("utf-8")
                if prettify:
                    pathway_xml = str(minidom.parseString(pathway_xml).toprettyxml(indent='   '))
                    pathway_xml = pathway_xml[pathway_xml.find('\n')+1:]
                #minidom does not work without xml_declaration

        return str(pathway_xml)


 
            


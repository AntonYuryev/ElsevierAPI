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
        self.graph.urn2rel.update(g.urn2rel)


    @classmethod
    def from_resnet(cls, resnet:et.Element, add_annotation=PSObject(),make_rnef=True):
        pathway = PSPathway()
        pathway['Name'] = [resnet.get('name')]
        pathway['URN'] = [resnet.get('urn')]
        pathway.update(add_annotation)
        [pathway.update_with_value(p.get('name'),p.get('value')) for p in resnet.findall('properties/attr')]

        if pathway.graph._parse_nodes_controls(resnet):
            if make_rnef:
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
                obj_count = len([obj_t for i,obj_t in self.graph.nodes.data('ObjTypeName') if obj_type in obj_t['ObjTypeName']])
                self['#'+obj_type] = obj_count
                return obj_count
        else:
            return len(self.graph)


    def merge_pathway(self, other:"PSPathway"):
        '''
        Merges pathway.graph\n
        Removes properties 'layout','resnet' as they are no longer applicable
        '''
        self.pop('layout','')
        self.pop('resnet','')
        for prop_name,values in other.items():
            if prop_name not in ['layout','resnet']:
                self.update_with_list(prop_name,values)

        self.graph = self.graph.compose(other.graph)


    def member_dbids(self,with_values:list=[],in_properties:list=['ObjTypeName']):
        return self.graph.dbids4nodes(with_values,in_properties)

    
    def get_members(self,with_values:list=[],in_properties=['ObjTypeName']):
        return self.graph.psobjs_with(in_properties,with_values)


    def props2rnef(self):
        props_element = et.Element('properties')
        for prop_name,prop_val in self.items():
            for val in prop_val:
                et.SubElement(props_element, 'attr', {'name':str(prop_name), 'value':str(val)})
        return str(et.tostring(props_element,encoding='utf-8',xml_declaration=True).decode("utf-8"))


    def to_xml_s(self,ent_props:list,rel_props:list,add_props2rel=dict(),add_props2pathway=dict(),format='RNEF',as_batch=True, prettify=True):
        '''
        Returns
        -------
        pathway XML string. supports RNEF and SBGN XML formats
        '''
        self.graph.load_references()
        graph_xml = self.graph.to_rnefstr(ent_props,rel_props,add_props2rel,add_props2pathway)
        pathway_xml = et.fromstring(graph_xml)

        pathway_props = et.SubElement(pathway_xml, 'properties')
        for prop_name,prop_val in self.items():
            for val in prop_val:
                et.SubElement(pathway_props, 'attr', {'name':str(prop_name), 'value':str(val)})

        try:
            my_layout = et.fromstring(self['layout'])
            lay_out = et.Element('attachments')
            lay_out.append(my_layout)
        except KeyError: pass
        pathway_xml.append(lay_out)

        if format == 'SBGN':
            pathway_xml_str = et.tostring(pathway_xml,encoding='utf-8',xml_declaration=True).decode("utf-8")
            pathway_xml_str = rnef2sbgn_str(pathway_xml)
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
                batch_xml = et.Element('batch')
                batch_xml.insert(0,pathway_xml)
                pathway_xml_str = et.tostring(batch_xml,encoding='utf-8',xml_declaration=True).decode("utf-8")
                if prettify: pathway_xml_str = minidom.parseString(pathway_xml_str).toprettyxml(indent='   ')
            else:
                pathway_xml_str = et.tostring(pathway_xml,encoding='utf-8',xml_declaration=True).decode("utf-8")
                if prettify:
                    pathway_xml_str = str(minidom.parseString(pathway_xml_str).toprettyxml(indent='   '))
                    pathway_xml_str = pathway_xml[pathway_xml_str.find('\n')+1:]
                #minidom does not work without xml_declaration

        return str(pathway_xml_str)


    def rnef2file(self,fname, ent_props:list,rel_props:list,add_props2rel=dict(),add_props2pathway=dict()):
        self.update(add_props2pathway)
        with open(fname,'w',encoding='utf-8') as f:
            f.write('<batch>\n') 
            props_str = self.props2rnef()
            f.write(props_str+'\n')
            self.graph.dump2rnef(f,ent_props,rel_props,add_props2rel)
            f.write('</batch>') 

 
            


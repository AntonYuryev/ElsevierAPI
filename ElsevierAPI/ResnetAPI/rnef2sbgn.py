
from xml.dom import minidom
from xml.etree.ElementTree import Element, tostring, parse, fromstring

def count_occurrences(d, _id, tag):
    ret = 0
    for k in d:
        if _id in d[k] and 'endpoint' in d[k][_id] and d[k][_id]['endpoint'] == tag:
            ret += 1
    return ret


def get_properties(elements, names):
    ret = {}
    for element in elements:
        _id = element.attrib['local_id']
        _properties = {}
        for attr in element:
            if 'name' not in attr.attrib:
                continue
            if attr.attrib['name'] in names:
                if attr.attrib['name'] not in _properties:
                    _properties[attr.attrib['name']] = set()
                _properties[attr.attrib['name']].add(attr.attrib['value'])
        if _id not in ret:
            ret[_id] = {}
        ret[_id] = {
            k: ret[_id][k].union(_properties[k]) if k in ret[_id] and k in _properties \
                else ret[_id][k] if k in ret[_id] \
                else _properties[k] \
            for k in set(ret[_id].keys()).union(set(_properties.keys()))
        }
    return ret


def read_scene(rnefStr, plot_scale):
    ret = {
        'glyphs': {},
        'links': {},
        'arcs': {},
        'g_map': {},
        'g_map_rev': {},
        'a_map': {},
        'a_map_rev': {}
    }
    rnef = fromstring(rnefStr)
    nodes = rnef.find('./resnet/nodes')
    controls = rnef.find('./resnet/controls')
    vobjs = rnef.find('./resnet/attachments/layout/scene/vobjs')
    vlinks = rnef.find('./resnet/attachments/layout/scene/vlinks')
    nodeprops = get_properties(nodes, ['NodeType', 'Name'])
    controlprops = get_properties(controls, ['ControlType', 'Effect'])
    # any other property of a node or control from RNEF (<attr> element) can be retrieved the same way
    for vobj in vobjs:
        if vobj.attrib['type'] not in ['Node', 'Control']:
            continue
        _id, _type, _ref = [vobj.attrib[x] for x in ['local_id', 'type', 'ref']]
        _x, _y, _w, _h = None, None, None, None
        for attr in vobj:
            if attr.attrib['name'] == 'Position':
                _x, _y = [float(x) * plot_scale for x in attr.attrib['value'].split(' ')]
            if attr.attrib['name'] == 'Size':
                _w, _h = [float(x) * plot_scale for x in attr.attrib['value'].split(' ')]
        if _id not in ret['glyphs']:
            for x in [_id, _type, _ref, _x, _y, _w, _h]:
                assert x is not None
            ret['glyphs'][_id] = {
                'type': _type,
                'ref': _ref,
                'x': _x,
                'y': _y,
                'w': _w,
                'h': _h,
                # put properties retrieved earlier using get_properties() function into ret['glyphs'][_id] structure
                'class': next(iter(nodeprops[_ref]['NodeType'])) if _ref in nodeprops else next(
                    iter(controlprops[_ref]['ControlType'])),
                'label': next(iter(nodeprops[_ref]['Name'])) if _ref in nodeprops else next(
                    iter(controlprops[_ref]['Effect'])).lower() if 'Effect' in controlprops[_ref] else 'default'
            }
    for vlink in vlinks:
        _src_ref, _dst_ref = [vlink.attrib[x] for x in ['src_ref', 'dst_ref']]
        _arrow, _points = '', []
        for attr in vlink:
            if attr.attrib['name'] == 'ArrowheadShape':
                _arrow = attr.attrib['value']
            if attr.attrib['name'] == 'Points':
                _points = [float(x) * plot_scale for x in attr.attrib['value'].split(' ')]
        if (_src_ref, _dst_ref) not in ret['links']:
            ret['links'][(_src_ref, _dst_ref)] = {
                'arrow': _arrow,
                'points': _points
            }
    for link in ret['links']:
        # control is source, node is target
        if ret['glyphs'][link[0]]['type'] == 'Control':
            arc_id = link[0]
            if arc_id not in ret['arcs']:
                ret['arcs'][arc_id] = {}
            ret['arcs'][arc_id][link[1]] = ret['links'][link]
            ret['arcs'][arc_id][link[1]]['endpoint'] = 'target'
            if ret['arcs'][arc_id][link[1]]['arrow'] == 'Undirected':
                ret['arcs'][arc_id][link[1]]['endpoint'] = 'source'
        # node is source, control is target
        if ret['glyphs'][link[1]]['type'] == 'Control':
            arc_id = link[1]
            if arc_id not in ret['arcs']:
                ret['arcs'][arc_id] = {}
            ret['arcs'][arc_id][link[0]] = ret['links'][link]
            ret['arcs'][arc_id][link[0]]['endpoint'] = 'source'

    extra_arcs = {}
    for _id in ret['arcs']:
        if len(ret['arcs'][_id]) > 2:
            # more than 2 participants
            extra_arcs_count = 0
            for participant_id in ret['arcs'][_id]:
                participant = ret['arcs'][_id][participant_id]
                extra_arcs_count += 1
                extra_arcs['%s:%d' % (_id, extra_arcs_count)] = {}
                extra_arcs['%s:%d' % (_id, extra_arcs_count)][participant_id] = participant
                extra_arcs['%s:%d' % (_id, extra_arcs_count)][_id] = {
                    'arrow': 'Undirected' if participant['arrow'] != 'Undirected' else 'Arrow',
                    'points': [],
                    'endpoint': 'source' if participant['endpoint'] != 'source' else 'target'
                }
                ret['glyphs']['%s:%d' % (_id, extra_arcs_count)] = {
                    'type': 'Control',
                    'ref': ret['glyphs'][_id]['ref'],
                    'x': ret['glyphs'][_id]['x'],
                    'y': ret['glyphs'][_id]['y'],
                    'w': ret['glyphs'][_id]['w'],
                    'h': ret['glyphs'][_id]['h'],
                    'class': ret['glyphs'][_id]['class'],
                    'label': ret['glyphs'][_id]['label']
                }
            ret['glyphs']['%s' % (_id)] = {
                'type': 'Node',
                'ref': ret['glyphs'][_id]['ref'],
                'x': ret['glyphs'][_id]['x'],
                'y': ret['glyphs'][_id]['y'],
                'w': ret['glyphs'][_id]['w'],
                'h': ret['glyphs'][_id]['h'],
                'class': 'and',
                'label': ret['glyphs'][_id]['label']
            }
            ret['arcs'][_id] = {}
    ret['arcs'].update(extra_arcs)

    # This will count occurrences of all nodes:
    # sources_cnt, targets_cnt = {}, {}
    # for _id in ret['arcs']:
    #     if not ret['arcs'][_id]:
    #         continue
    #     for participant_id in ret['arcs'][_id]:
    #         sources_cnt[participant_id] = count_occurrences(ret['arcs'], participant_id, 'source')
    #         targets_cnt[participant_id] = count_occurrences(ret['arcs'], participant_id, 'target')

    for _id in ret['arcs']:
        if not ret['arcs'][_id]:
            continue
        if not ('target' in [ret['arcs'][_id][link_id]['endpoint'] for link_id in ret['arcs'][_id]]):
            # 2 participants, but both are sources
            # artificial rule:
            # target will be the one with first type in alphabetical order (Small Molecule --> Complex)
            # If both participants have same type, target will be the one with first local id in alphabetical order
            # (completely meaningless biology-wise)
            participants = [(participant_id, ret['glyphs'][participant_id]['class']) for participant_id in
                            ret['arcs'][_id]]
            if participants[0][1] != participants[1][1]:
                participants.sort(key=lambda x: x[1])
            else:
                participants.sort(key=lambda x: x[0])
            ret['arcs'][_id][participants[0][0]]['endpoint'] = 'target'

    g_count = 0
    for _id in ret['glyphs']:
        g_count += 1
        glyph_id = 'glyph%d' % (g_count)
        ret['g_map'][_id] = glyph_id
        ret['g_map_rev'][glyph_id] = [_id]

    a_count = 0
    for _id in ret['arcs']:
        if not ret['arcs'][_id]:
            continue
        a_count += 1
        arc_id = 'arc%d' % (a_count)
        ret['a_map'][_id] = arc_id
        ret['a_map_rev'][arc_id] = [_id]

    return ret


def read_classmap(filename, language):
    ret = {'language': language, 'nodes': {}, 'controls': {}}
    root = parse(filename)
    for classmap in root.find('.'):
        if classmap.attrib['language'] != language:
            continue
        for element in classmap:
            if element.tag == 'entity':
                ret['nodes'][element.attrib['rnef']] = element.attrib['sbgn']
            if element.tag == 'relation':
                rnef = element.attrib['rnef']
                effect = element.attrib['effect']
                sbgn = element.attrib['sbgn']
                if rnef not in ret['controls']:
                    ret['controls'][rnef] = {}
                ret['controls'][rnef][effect] = sbgn
    return ret


def compile_sbgn(scene, classmap):
    ret = Element('sbgn')
    ret.attrib['xmlns'] = 'http://sbgn.org/libsbgn/0.3'
    _map = Element('map')
    _map.attrib['language'] = classmap['language']
    for _id in scene['glyphs']:
        if scene['glyphs'][_id]['type'] == 'Node':
            glyph = Element('glyph')
            glyph.attrib['id'] = scene['g_map'][_id]
            glyph.attrib['class'] = classmap['nodes'][scene['glyphs'][_id]['class']] if scene['glyphs'][_id]['class'] in \
                                                                                        classmap['nodes'] else \
                scene['glyphs'][_id]['class']
            label = Element('label')
            label.attrib['text'] = scene['glyphs'][_id]['label']
            glyph.append(label)
            bbox = Element('bbox')
            for k in ['x', 'y', 'w', 'h']:
                bbox.attrib[k] = str(scene['glyphs'][_id][k])
                bbox.attrib['y'] = str(-scene['glyphs'][_id]['y'])  # Y coordinate in SBGN is supposed to be inverted
            glyph.append(bbox)
            _map.append(glyph)
    for _id in scene['arcs']:
        if not scene['arcs'][_id]:
            continue
        arc = Element('arc')
        arc.attrib['id'] = scene['a_map'][_id]
        effect = scene['glyphs'][_id]['label']
        if effect not in ['negative', 'positive']:
            effect = 'default'
        sc_class = scene['glyphs'][_id]['class']
        try:
            arc.attrib['class'] = classmap['controls'][sc_class][effect]
        except KeyError:
            print('!!!Relation "%s" with effect %s is not specified in rnef2sbgn_map.xml!!!' % (sc_class,effect))
            continue

        participants = [(participant_id, scene['arcs'][_id][participant_id]['endpoint']) for participant_id in
                        scene['arcs'][_id]]
        participants.sort(key=lambda x: x[1])
        [arc.attrib['source'], arc.attrib['target']] = [scene['g_map'][participant[0]] for participant in participants]
        arc_start_x = str(scene['glyphs'][scene['g_map_rev'][arc.attrib['source']][0]]['x'])
        arc_start_y = str(-scene['glyphs'][scene['g_map_rev'][arc.attrib['source']][0]][
            'y'])  # Y coordinate in SBGN is supposed to be inverted
        arc_end_x = str(scene['glyphs'][scene['g_map_rev'][arc.attrib['target']][0]]['x'])
        arc_end_y = str(-scene['glyphs'][scene['g_map_rev'][arc.attrib['target']][0]][
            'y'])  # Y coordinate in SBGN is supposed to be inverted
        arc_points_source, arc_points_target, arc_points = [], [], []

        # The following block attempts to transfer the inflections of arcs, but the result looks ugly
        # for participant_id in scene['arcs'][_id]:
        #     participant = scene['arcs'][_id][participant_id]
        #     for p in range(0, len(participant['points']), 2):
        #         if participant['endpoint'] == 'target':
        #             arc_points_target.append((str(participant['points'][p]), str(-participant['points'][p+1]))) # Y coordinate in SBGN is supposed to be inverted
        #         else:
        #             arc_points_source.append((str(participant['points'][p]), str(-participant['points'][p+1]))) # Y coordinate in SBGN is supposed to be inverted

        arc_points = arc_points_source + arc_points_target
        arc_start = Element('start')
        arc_start.attrib['x'], arc_start.attrib['y'] = arc_start_x, arc_start_y
        arc.append(arc_start)
        for point in arc_points:
            arc_next = Element('next')
            arc_next.attrib['x'], arc_next.attrib['y'] = point[0], point[1]
            arc.append(arc_next)
        arc_end = Element('end')
        arc_end.attrib['x'], arc_end.attrib['y'] = arc_end_x, arc_end_y
        arc.append(arc_end)
        _map.append(arc)
    ret.append(_map)
    return ret


def make_file_name(pathwayName: str):
    new_str = list(pathwayName)
    for i in range(0, len(new_str)):
        if new_str[i] in {'>', '<', '|','/'}:
            new_str[i] = '-'
        if new_str[i] in {':'}:
            new_str[i] = '_'
    return "".join(new_str)


def rnef2sbgn_str(rnef:str, classmapfile='ElsevierAPI/ResnetAPI/sbgn/rnef2sbgn_map.xml', language='activity flow', plot_scale=30):
    #pretty_rnef = minidom.parseString(rnef).toprettyxml(indent='   ')
    scene = read_scene(rnef, plot_scale)
    classmap = read_classmap(classmapfile, language)
    sbgn = compile_sbgn(scene, classmap)
    return minidom.parseString(tostring(sbgn)).toprettyxml(indent='   ')


def to_sbgn_file(rnef:str, fname:str, outdir:str):
    sbgn_xml = rnef2sbgn_str(rnef)
    fout_name = make_file_name(fname)
    with open(outdir+fout_name+'.sbgn', 'w', encoding='utf-8') as f:
        f.write(sbgn_xml)


def do_the_job(filename_in, dir_out, language, plot_scale, classmap='sbgn/rnef2sbgn_map.xml'):
    rnef = parse(filename_in).getroot()
    for resnet in rnef.findall('./resnet'):
        try:
            pathway_name = str(resnet.attrib['name'])
        except KeyError:
            print('Pathway must have a name')
            continue
        pathway_name = make_file_name(pathway_name)

        resnet_str = tostring(resnet, encoding='utf-8').decode("utf-8")
        resnet_str = "<?xml version='1.0' encoding='UTF-8'?><batch>" + resnet_str + '</batch>'
        rnef_out = dir_out + '\\' + pathway_name + '.rnef'
        with open(rnef_out, mode='w', encoding='utf8') as f1:
            f1.write(resnet_str)

        sbgn_str = rnef2sbgn_str(resnet_str, classmap, language, plot_scale)
        sbgn_out = dir_out + '\\' + pathway_name + '.sbgn.xml'
        with open(sbgn_out, mode='w', encoding='utf8') as f2:
            f2.write(sbgn_str)


if __name__ == '__main__':
    # Usage:
    # rnef2sbgn.py FILENAME_IN.RNEF FILENAME_OUT SBGN_TYPE PLOT_SCALE
    [fin, dir_out, lang, s_plot_scale] = ['COVID19 models backup.rnef', 'COVID19 models', "activity flow", 30]
    plot_scale = int(s_plot_scale)
    do_the_job(fin, dir_out, lang, plot_scale)

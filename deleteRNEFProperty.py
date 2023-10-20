from ElsevierAPI.ResnetAPI.ResnetGraph import ResnetGraph, PSRelation
import argparse, textwrap

if __name__ == "__main__":
    instructions = '''
    -i, --infile: full path to RNEF XML file
    -p, --properties: comma-separated list of properties to delete. 
    Example: -p "c,bad,wrong" will remove properties named "c", "bad", "wrong"
    '''

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog=textwrap.dedent(instructions))
    parser.add_argument('-i', '--infile', type=str)
    parser.add_argument('-p', '--properties', type=str)
    
    args = parser.parse_args()
    props2delete = str(args.properties).split(',')

    my_rnef = ResnetGraph.fromRNEF(str(args.infile))
    print(f'Read {str(args.infile)} file with {my_rnef.number_of_nodes()} nodes, {my_rnef.number_of_edges()} relations and {len(my_rnef.load_references())} references')
    for r,t,urn,rel in my_rnef.edges.data(True,keys=True):
        clean_rel = rel['relation'].remove_props(props2delete)
        my_rnef[r][t][urn]['relation'] = clean_rel

    fout = str(args.infile)[:str(args.infile).rfind('.')]+"_clean.rnef"
    my_rnef.dump2rnef(fout,with_section_size=1000)
    print(f'Created {fout} file with {my_rnef.number_of_nodes()} nodes, {my_rnef.number_of_edges()} relations and {len(my_rnef.load_references())} references')

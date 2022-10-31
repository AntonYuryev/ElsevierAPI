from ElsevierAPI import open_api_session
import pandas as pd
import argparse
import textwrap


if __name__ == "__main__":
    instructions = '''
                infile - tab-delimted file with mapping values in the 1st column. Specify mapping attribute in the column header.
                all other columns - new annotation values. Use column header to name the property
                '''

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog=textwrap.dedent(instructions))
    parser.add_argument('-i', '--infile', type=str,required=True)
    args = parser.parse_args()

    annotation_pd = pd.read_csv(args.infile, delimiter='\t', header=0)
    mapping_property = annotation_pd.columns[0]
    map2attributes = list(annotation_pd[mapping_property])

    api_config = 'D:/Python/ENTELLECT_API/ElsevierAPI/APIconfigLanzaTech.json'
    api_config = ''
    ps_api = open_api_session(api_config)
    ps_api.map_props2objs(map2attributes,[mapping_property]) # have to fetch objects for mapping first

    annotation_cols = annotation_pd.columns.to_list()[1:]
    for col_name in annotation_cols:
        map_dict = pd.Series(annotation_pd[col_name].values,index=annotation_pd[mapping_property]).to_dict()
        map_dict = {str(k):[str(v)] for k,v in map_dict.items()}
        ps_api.Graph.add_node_annotation(with_new_prop=col_name,map2prop=mapping_property,using_map=map_dict)

    ps_api.add_ent_props(annotation_cols)
    rnef_file = str(args.infile)[:-3]+'rnef'
    ps_api.graph2rnef(rnef_file)
    print('Entity annotation is in %s file' % rnef_file)
#import rnef_file into Pathway Studio to annotate proteins in the database



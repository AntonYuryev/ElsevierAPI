import time
import pandas as pd
import argparse
import textwrap
from ElsevierAPI.ResnetAPI.PathwayStudioGOQL import OQL
from ElsevierAPI import load_api_config
from ElsevierAPI.ResnetAPI.SemanticSearch import SemanticSearch

        
if __name__ == "__main__":
    EntityListFile = str()
    link2concepts = pd.DataFrame()
    instructions = '''
    infile - tab-delimted file with entity idnetifiers. Specify identifier type in each column header. 
    Identifiers will be used in the order of the columns: if entity is not found by identifier from the first column, identfier from the second column will be used.
    
    infile_has_header - specifies if infile has a header.  Defaults to no header. Use "Y" if header is present. 
    If header is not specified only identifiers in the 1st column will be used for mapping by Name+Alias to database entities

    entity_types - comma-separated objectTypeNames for entities in infile

    targets_file - tab-delimted file, alternative input for --targets option. Header: SearchPropertyName<>Effect<>RelationType<>Direction
        SearchPropertyName - Required. concepts that must be semantically linked to entities from infile. Header value defines what properties must be used to retreive concepts. Defaults to 'Name,Alias'
        Effect - optional.  The effect sign   ebetween entity in infile and concept in targets_file. Defaults to any.
        RelationType - optional. The type of relation between entity in infile and concept in targets_file. Defaults to any.
        Direction - optional. The direction of relation between entity in infile and concept in targets_file. Allowed values: '>', '<',''. Defaults to any.

    pathways - comma-separated list of pathway names which must be semantically linked to entities from infile.
    Semantic reference count to a pathway is calculated as sum of reference count of semantic links to every pathway component

    Output: tab-delimted file with rows from infile annotated with semantic reference count between entity (row) and semantic concept (column)  

    dump_references - create additional dump file containing all references supporitng relations found by the script
    resnet_retreive_props - list of entity properties to include into dump file
    use_cache_to_resume - must be 'Y' if you need to continue interrupted download. Will look for map_file.tsv,reference_cache.tsv,counts_cache.tsv to load pre-mapped database identifiers - saves time when infile is re-used. Defaults to start from scratch.
   
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog=textwrap.dedent(instructions))
    parser.add_argument('-i', '--infile', type=str,required=True)
    parser.add_argument('-H', '--infile_has_header',action='store_true')
    parser.add_argument('-E', '--entity_types', type=str,default='')
    parser.add_argument('-f', '--targets_file', type=str,required=True)
    parser.add_argument('-d', '--dump_references', action='store_true')
    parser.add_argument('-a', '--resnet_retreive_props', type=str)
    parser.add_argument('-p', '--pathways', type=str, default='')
    parser.add_argument('-c', '--use_cache_to_resume', action='store_true')
    parser.add_argument('--debug', action="store_true")

    args = parser.parse_args()

    if len(args.infile) > 0:
        EntityListFile = args.infile #full path to tab-delimted file with drugs in the first column
    else:
        print('No input file with entities was specified!!!')

    concepts_in_file = 0
    link2concepts = pd.read_csv(args.targets_file, delimiter='\t', header=0, index_col=0)
    concepts_in_file = len(link2concepts)
   
    #Script will work faster if you specify what list of object types in your input Entity list
    entity_types = str(args.entity_types)
    entity_types = entity_types.replace("Small molecule",'SmallMol')
    entity_types = str(entity_types).split(',') if len(args.entity_types) > 0 else [] 

    #If you want to use identifiers other than names enter appropriate Propeprty types into SearchByProperties list
    #consult with "Resnet Entities&Properties.txt" for the list of available identifier for nodes in the knowledge graph
    global_start_time = time.time()
    header_pos = 0 if args.infile_has_header else None

    search = SemanticSearch(load_api_config())
    AnnotateNodesWith = str(args.resnet_retreive_props).split(',') if len(args.resnet_retreive_props) > 0 else []
    search.add_ent_props(AnnotateNodesWith)
    if args.dump_references:
        search.add_rel_props(['Name','Sentence','PubYear','Title'])

    if args.use_cache_to_resume:
        search.load_pandas(None,None,use_cache=True)
    else:
        EntityPandas = pd.read_csv(EntityListFile,delimiter='\t', header=header_pos)
        search.load_pandas(EntityPandas,args.infile_has_header, use_cache=False, map2type=entity_types)
         
    if len(args.pathways) > 0:
        print("Begin linking entities mapped from infile to pathways")
        LinkToPathways = str(args.pathways).split(",")
        start_time = time.time()
        for PathwayName in LinkToPathways:
            PathwayMembersId2Entity = search.get_pathway_member_ids(
                [PathwayName], search_pathways_by=['Name'], only_entities=['Protein', 'FunctionalClass', 'Complex'], with_properties=['objectType'])
            pathway_components = set(PathwayMembersId2Entity.keys())
            QueryOntology = OQL.get_childs(list(pathway_components), ['id'])
            pathway_components.update(search._obj_id_by_oql(QueryOntology))

            if len(pathway_components) == 0:
                print ('No entity for %s found in the database' % PathwayName)
            else:
                search.link2concept(PathwayName, list(pathway_components))

        exec_time = search.execution_time(start_time)
        print("Entities in file %s were linked to %s pathway in %s" % (EntityListFile, PathwayName, exec_time))
    else:
        print('No pathways were specified for semantic linking with entities from \"%s\"' % (EntityListFile))


    if len(link2concepts) > 0:
        search_concept_by = str(link2concepts.index.name).split(',')
        print ('Will search concepts in %s by %s' % (args.targets_file,link2concepts.index.name))
        
        print("\nBegin linking entities mapped from infile to %d concepts from %s" % (concepts_in_file,args.targets_file))
        for concept_idx in range(0, len(link2concepts.index)):
            link2concept = link2concepts.index[concept_idx]
            concept_ids = search._get_obj_ids_by_props(
                propValues=[link2concept], search_by_properties=search_concept_by)
 
            if len(concept_ids) == 0: print ('No entity %s found in database' % link2concept)
            else:
                connect_by_rels = list()
                rel_effect = list()
                rel_dir = ''
                
                if 'RelationType' in link2concepts.columns:
                    rel_t = link2concepts.at[link2concept,'RelationType']
                    if pd.notna(rel_t): connect_by_rels= str(rel_t).split(',')

                if 'Effect' in link2concepts.columns:
                    ef = link2concepts.at[link2concept,'Effect']
                    if pd.notna(ef): rel_effect = str(ef).split(',')

                if 'Direction' in link2concepts.columns:
                    dr = link2concepts.at[link2concept,'Direction']
                    if pd.notna(dr): rel_dir = str(dr)

                search.set_how2connect(connect_by_rels,rel_effect,rel_dir)

                linked_entities_count= search.link2concept(link2concept,list(concept_ids))
                print("Total %d out of %d concepts interrogated in %s" % 
                     (concept_idx + 1, len(link2concepts), search.execution_time(global_start_time)))

    exec_time = search.execution_time(global_start_time)
    print("%d entities in file %s were mapped linked to %d concepts from %s in %s" % 
         (len(EntityPandas), EntityListFile, concepts_in_file, args.targets_file, exec_time))

    countsOutFile = EntityListFile[:len(EntityListFile)-4]+'+SemanticRefcount.tsv'
    refOutFile = EntityListFile[:len(EntityListFile)-4]+'+SemanticReferences.tsv' if args.dump_references else ''
    search.print_ref_count(countsOutFile,refOutFile)


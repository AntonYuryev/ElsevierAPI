import time
import pandas as pd
import argparse
import textwrap
from ElsevierAPI.ResnetAPI.PathwayStudioGOQL import OQL
from ElsevierAPI import load_api_config
from ElsevierAPI.ResnetAPI.SemanticSearch import SemanticSearch,len

        
if __name__ == "__main__":
    EntityListFile = str()
    link2concepts = pd.DataFrame()
    instructions = '''
    infile - tab-delimted file with entity idnetifiers. Specify identifier type in each column header. 
    Identifiers will be used in the order of the columns: if entity is not found by identifier from the first column, identfier from the second column will be used.
    
    infile_has_header - specifies if infile has a header.  Defaults to no header. Use "Y" if header is present. 
    If header is not specified only identifiers in the 1st column will be used for mapping by Name+Alias to database entities

    entity_types - comma-separated objectTypeNames for entities in infile

    targets_file - tab-delimted file. Header: SearchPropertyName<>Effect<>RelationType<>Direction
        SearchPropertyName - Required. concepts that must be semantically linked to entities from infile. Header value defines what properties must be used to retreive concepts. Defaults to 'Name,Alias'
        Effect - optional.  The effect sign   ebetween entity in infile and concept in targets_file. Defaults to any.
        RelationType - optional. The type of relation between entity in infile and concept in targets_file. Defaults to any.
        Direction - optional. The direction of relation between entity in infile and concept in targets_file. Allowed values: '>', '<',''. Defaults to any.

    pathways - comma-separated list of pathway names which must be semantically linked to entities from infile.
    Semantic reference count to a pathway is calculated as sum of reference count of semantic links to every pathway component

    Output: tab-delimted file with rows from infile annotated with semantic reference count between entity (row) and semantic concept (column)  

    add_references - create additional dump file containing all references supporitng relations found by the script
    retreive_entity_props - list of entity properties to include into dump file
    use_cache_to_resume - must be 'Y' if you need to continue interrupted download. Will look for map_file.tsv,reference_cache.tsv,counts_cache.tsv to load pre-mapped database identifiers - saves time when infile is re-used. Defaults to start from scratch.
    '''

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog=textwrap.dedent(instructions))
    parser.add_argument('-i', '--infile', type=str)
    parser.add_argument('-H', '--infile_has_header',action='store_true')
    parser.add_argument('-E', '--entity_types', type=str,default='')
    parser.add_argument('-f', '--targets_file', type=str)
    parser.add_argument('-d', '--add_references', action='store_true')
    parser.add_argument('-a', '--retreive_entity_props', type=str)
    parser.add_argument('-p', '--pathways', type=str, default='')
    parser.add_argument('-c', '--use_cache_to_resume', action='store_true')
    parser.add_argument('--debug', action="store_true")

    args = parser.parse_args()

    if args.use_cache_to_resume:
        search = SemanticSearch(load_api_config(),need_references=args.add_references,need_snippets=args.add_references)
        search.load_pandas(None,None,use_cache=True)
    elif args.infile:
        #Script will work faster if you specify what list of object types in your input Entity list
        search = SemanticSearch(load_api_config(),need_references=args.add_references,need_snippets=args.add_references)
        header_pos = 0 if args.infile_has_header else None
        entity_types = str(args.entity_types)
        entity_types = entity_types.replace("Small molecule",'SmallMol')
        entity_types = str(entity_types).split(',')
        EntityListFile = args.infile #full path to tab-delimted file with drugs in the first column
        link2concepts = pd.read_csv(args.targets_file, delimiter='\t', header=0, index_col=0)
        EntityPandas = pd.read_csv(EntityListFile,delimiter='\t', header=header_pos)
        search.RefCountPandas = search.RefCountPandas = search.load_pandas(EntityPandas,args.infile_has_header,map2type=entity_types)
        AnnotateNodesWith = str(args.resnet_retreive_props).split(',') if args.resnet_retreive_props else []
        report_fname = EntityListFile[:len(EntityListFile)-4]+'+SemanticRefcount'
    else:
        # add your code here to select entities from the database
        search = SemanticSearch(load_api_config(),need_references=True,need_snippets=True)
        entity_types = ['Protein']
        AnnotateNodesWith = ['Class']
        concept_name = 'cardiovascular disease' #'chronic obstructive pulmonary disease' #
        oql_query = "SELECT Entity WHERE objectType=({t}) AND Organism = 'Homo sapiens' AND Connected by ({r}) to ({d})"
        disease = "SELECT Entity WHERE InOntology (SELECT Annotation WHERE Ontology='Pathway Studio Ontology' AND Relationship='is-a') under (SELECT OntologicalNode WHERE (Name,Alias)=('{c}')) OR Name = '{c}'".format(c=concept_name)
        rel = 'SELECT Relation WHERE objectType =(Regulation,QuantitativeChange,Biomarker,GeneticChange,StateChange)'
        entity_types_str = ','.join(entity_types)
        oql_query = oql_query.format(t=entity_types_str,r=rel,d=disease)
        entity_graph = search.process_oql(oql_query, 'Find disease targets')
        id2names = entity_graph.get_properties('Name')
        protein_names = [v[0] for v in id2names.values()]
        link2concepts = pd.DataFrame.from_dict({'Name':[concept_name]})
        link2concepts.set_index('Name',inplace=True)
        EntityPandas = pd.DataFrame.from_dict({'Name':protein_names})
        search.RefCountPandas = search.load_pandas(EntityPandas,prop_names_in_header=True, map2type=entity_types)
        report_fname = 'SemanticRefcount4entities linked to {}'.format(concept_name)
        
    concepts_count = len(link2concepts)
    search.add_ent_props(AnnotateNodesWith)
    
    #If you want to use identifiers other than names enter appropriate Propeprty types into SearchByProperties list
    #consult with "Resnet Entities&Properties.txt" for the list of available identifier for nodes in the knowledge graph
    global_start_time = time.time()
        
    if len(args.pathways) > 0:
        print("Begin linking entities mapped from infile to pathways")
        LinkToPathways = str(args.pathways).split(",")
        start_time = time.time()
        for PathwayName in LinkToPathways:
            PathwayMembersId2Entity = search.get_pathway_members(
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
        print ('Will search %d concepts by %s' % (len(link2concepts),link2concepts.index.name))
        
        print("\nBegin linking entities mapped from infile to %d concepts" % (concepts_count))
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
         (len(EntityPandas), EntityListFile, concepts_count, args.targets_file, exec_time))

    search.add2report(search.make_count_df())
    if args.add_references:
        search.add_snippets()
    search.print_report(report_fname+'.xlsx')

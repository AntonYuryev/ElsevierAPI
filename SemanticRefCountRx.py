from ElsevierAPI import APIconfig
import time
import pandas as pd
import argparse
import textwrap
from ElsevierAPI.ResnetAPI.SemanticSearch import SemanticSearch
import ElsevierAPI.ResnetAPI.PathwayStudioGOQL as OQL
import ElsevierAPI.ReaxysAPI.Reaxys_API as RxAPI


LazyIdDict = dict()
def AddSemanticRefCount(EntityPandas:pd.DataFrame, LinkToEntityName:str,LinkToEntityIDs:list, ps_api: SemanticSearch):
    Entities = list(EntityPandas.index)
    NewColumnName = "RefCount to "+ LinkToEntityName
    EntityPandas.insert(len(EntityPandas.columns), NewColumnName, 0)
    
    for Entity in Entities:
        if Entity in LazyIdDict.keys():
            protids = LazyIdDict[Entity]
        else:
            protids = ps_api._get_obj_ids_by_props(propValues=[Entity],search_by_properties=['Name'])                     
            if len(protids) == 0 :
                protids = ps_api._get_obj_ids_by_props(propValues=[Entity],search_by_properties=['Alias']) 
                            
        LazyIdDict[Entity] = protids
        if len(protids) == 0:
            print('No protein with name or alias %s found in the database' % Entity)
                    
        refCount = set()
        if len(protids) > 0 : refCount = ps_api.semantic_refcount_by_ids(LinkToEntityIDs, protids)
                        
        if len(refCount) > 0 : print("pair %s:%s has %d semantic references" % (Entity,LinkToEntityName,len(refCount)))          
        else : print("No relations found for pair %s:%s" % (Entity,LinkToEntityName))

        EntityPandas.at[Entity,NewColumnName] = len(refCount)


if __name__ == "__main__":
    EntityListFile = str()
    LinkToEntities = list()
    instructions = '''
    infile - single column file with entity names. 
    If you want to use identifiers other than names enter appropriate Propeprty types into SearchByProperties list.

    targets - comma-separated list of concept that must be semantically linked to entities from infile
    target_type - objectTypeName for targets

    pathways - comma-separated list of pathway names which must be semantically linked to entities from infile.
    Semantic reference count to a pathway is calculated as sum of reference count of semantic links to every pathay component
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog=textwrap.dedent(instructions))
    parser.add_argument('-i', '--infile', type=str, required=True)
    parser.add_argument('-t', '--targets', type=str)
    parser.add_argument('-o', '--target_type', type=str,default="")
    parser.add_argument('-p', '--pathways', type=str, default="")
    parser.add_argument('--debug', action="store_true")
    args = parser.parse_args()

    if len(args.infile) > 0:
        EntityListFile = args.infile #full path to tab-delimted file with drugs in the first column
    else:
        print('No entity list file was specified')

    if len(args.targets) > 0:
        LinkToEntities = args.targets.split(",") #Enter comma-separated list of concepts that should be linked to your input Entity list as first argument 
    
    if len(args.target_type) > 0: #Script will work faster if you specify what list of object types in your input Entity list
        QueryObjType = args.target_type
        QueryObjType = QueryObjType.split(',')
    else:
        QueryObjType = []
        print('No object type for targets was specified. Will seacrh for entities in %s using all object types' % EntityListFile)

    if len(args.pathways) == 0:
        print('No pathways for linking with entities in \"%s\" were specified' % (EntityListFile))
    

    #If you want to use identifiers other than names enter appropriate Propeprty types into SearchByProperties list
    #consult with "Resnet Entities&Properties.txt" for the list of available identifier for nodes in the knowledge graph
    EntityPandas = pd.read_csv(EntityListFile, delimiter='\t', header=0, index_col=0)
    ps_api = SemanticSearch(APIconfig)

    if len(args.pathways) > 0:
        LinkToPathways = args.pathways
        LinkToPathways = LinkToPathways.split(",")
        start_time = time.time()
        for PathwayName in LinkToPathways:
            PathwayMembersIdToEntity = ps_api.get_pathway_member_ids([PathwayName], search_pathways_by=['Name'], FilterBy=['Protein', 'FunctionalClass', 'Complex'], InProperty='objectType')
            EntityIDs = set(PathwayMembersIdToEntity.keys())
            QueryOntology = OQL.get_childs(list(EntityIDs), ['id'])
            EntityIDs.update(ps_api.GetObjIds(QueryOntology))

            if len(EntityIDs) == 0:
                print ('No entity for %s found in the database' % PathwayName)
            else:
                AddSemanticRefCount(EntityPandas, PathwayName, EntityIDs)

        print("Entities in file %s were linked to %s pathway in %s" % (EntityListFile, PathwayName, PsAPI.ExecutionTime(start_time)))


    if len(args.targets) > 0:
        start_time = time.time()
        for LinkToEntityName in LinkToEntities:
            QueryNode = OQL.get_entities_by_props(PropertyValues=[LinkToEntityName], SearchByProperties=['Name', 'Alias'], only_object_types=QueryObjType)
            QueryOntology = OQL.get_childs(PropertyValues=[LinkToEntityName], SearchByProperties=['Name', 'Alias'], only_object_types=QueryObjType)
            EntityIDs = ps_api.GetObjIds(QueryNode)
            EntityIDs.update(ps_api.GetObjIds(QueryOntology))

            if len(EntityIDs) == 0: print ('No entity for %s found in the database' % LinkToEntityName)
            else: AddSemanticRefCount(EntityPandas,LinkToEntityName,EntityIDs)

        print("Entities in file %s were linked to %s in %s" % (EntityListFile, LinkToEntityName, ps_api.execution_time()(start_time) ))


    foutName = EntityListFile[:len(EntityListFile)-4]+'+SemanticRefcount.tsv'
    EntityPandas.to_csv(foutName, sep='\t', index=True)

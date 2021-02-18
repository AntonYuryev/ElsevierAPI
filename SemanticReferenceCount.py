import time
import pandas as pd
import argparse
import textwrap
import ElsevierAPI
import ElsevierAPI.ResnetAPI.PathwayStudioGOQL as OQL
from ElsevierAPI import networx as PSnx


def GenerateEntityIDdict(EntityPandas:pd.DataFrame, SearchByProp=['Name','Alias']):
    PandaEntityToID = dict()
    Entities = set(EntityPandas.index)
    for Entity in Entities:
        for p in SearchByProp:
            protids = PSnx.GetObjIdsByProps(PropertyValues=[Entity],SearchByProperties=[p])
            if len(protids) > 0: 
                PandaEntityToID[Entity] = protids
                break
        
        if Entity not in PandaEntityToID.keys():
            print('Cannot find %s entity in the database' % Entity)
    
    print ('Found %d from %d entities in the database', (len(PandaEntityToID), len(Entities)))
    return PandaEntityToID


def AddSemanticRefCount(EntityPandas:pd.DataFrame,PandaEntityToID:dict,LinkToConceptName:str,LinkToConceptIDs:list,ConnectByRelTypes=[],RelEffect=[],RelDirection='',REL_PROPS=[], ENT_PROPS=[]):
    Entities = set(EntityPandas.index)
    TotalEntotyCount = len(Entities)
    NewColumnName = "RefCount to "+ LinkToConceptName
    EntityPandas.insert(len(EntityPandas.columns),NewColumnName,0)

    effecStr = ','.join(RelEffect)
    if len(effecStr) == 0: effecStr = 'all'
    relTypeStr = ','.join(ConnectByRelTypes)
    if len(relTypeStr) == 0: relTypeStr = 'all'
    dirStr = ':'
    if len(RelDirection) > 0: dirStr = RelDirection

    entityCounter = 0
    for Entity in Entities:
        entityCounter += 1
        try: protids = PandaEntityToID[Entity]
        except KeyError: continue

        refCount = PSnx.SemanticRefCountByIds(LinkToConceptIDs,list(protids),ConnectByRelTypes,RelEffect,RelDirection,REL_PROPS=REL_PROPS, ENTITY_PROPS=ENT_PROPS)
        if len(refCount) > 0 : 
            print("%d out of %d pair %s:%s has %d semantic references for %s relations with %s effects" % (entityCounter, TotalEntotyCount, Entity,LinkToConceptName,len(refCount),relTypeStr,effecStr))
        else: print("No relations found for pair %s:%s" % (Entity,LinkToConceptName))

        EntityPandas.at[Entity,NewColumnName] = len(refCount)


if __name__ == "__main__":
    EntityListFile = str()
    LinkToConcepts = pd.DataFrame()
    instructions = '''
    infile - single column file with entity names. 
    If you want to use identifiers other than names enter appropriate Propeprty types into SearchByProperties list.
    infile_header - specifies if infile has a header.  Defaults to no header. Use "Y" if headeer is present 

    targets - comma-separated list of concept that must be semantically linked to entities from infile
    target_type - objectTypeName for targets
    targets_file - tab-delimted file. Header: SearchPropertyName<>Effect<>RelationType<>Direction
        SearchPropertyName - Required. concepts that must be semantically linked to entities from infile. Header value defines what properties must be used to retreive concepts. Defaults to 'Name,Alias'
        Effect - optional.  The effect sign   ebetween entity in infile and concept in targets_file. Defaults to any.
        RelationType - optional. The type of relation between entity in infile and concept in targets_file. Defaults to any.
        Direction - optional. The direction of relation between entity in infile and concept in targets_file. Allowed values: '>', '<',''. Defaults to any.

    pathways - comma-separated list of pathway names which must be semantically linked to entities from infile.
    Semantic reference count to a pathway is calculated as sum of reference count of semantic links to every pathway component
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog=textwrap.dedent(instructions))
    parser.add_argument('-i', '--infile', type=str, required=True)
    parser.add_argument('-H', '--infile_header',  default='')
    parser.add_argument('-t', '--targets', default='')
    parser.add_argument('-o', '--target_type', type=str,default='')
    parser.add_argument('-f', '--targets_file', type=str,default='')
    parser.add_argument('-d', '--dump_references', default='')
    parser.add_argument('-a', '--resnet_retreive_props', type=str)
    parser.add_argument('-p', '--pathways', type=str, default='')
    parser.add_argument('--debug', action="store_true")
    args = parser.parse_args()

    if len(args.infile) > 0:
        EntityListFile = args.infile #full path to tab-delimted file with drugs in the first column

    if len(args.targets_file) > 0:
        LinkToConcepts = pd.read_csv(args.targets_file, delimiter='\t', header=0, index_col=0)
    
    if len(args.targets) > 0:
        LinkToConcepts.append(args.targets.split(','))
    
    if len(args.target_type) > 0: #Script will work faster if you specify what list of object types in your input Entity list
        QueryObjType = args.target_type
        QueryObjType = QueryObjType.split(',')
    else:
        QueryObjType = []
        print('No object type for targets was specified. Will seacrh for entities in %s using all object types' % EntityListFile)

    if len(args.pathways) == 0:
        print('No pathways for linking with entities in \"%s\" were specified' % (EntityListFile))

    if args.dump_references in ['True', 'true','yes','y','Y']:
        REL_PROP = ['Name','Sentence','PubYear','Title', 'PMID', 'DOI']
        ENT_PROPS = args.resnet_retreive_props.split(',')
    else:
        REL_PROP = []
    

    #If you want to use identifiers other than names enter appropriate Propeprty types into SearchByProperties list
    #consult with "Resnet Entities&Properties.txt" for the list of available identifier for nodes in the knowledge graph
    has_header = None
    if args.infile_header in ['True','true','yes','y','Y']: has_header = 0
    EntityPandas = pd.read_csv(EntityListFile,delimiter='\t',header=has_header,index_col=0)
    print ('Findings %d entities in %d rows of %s in the database' % (len(set(EntityPandas.index)),len(EntityPandas),EntityListFile))
    PandaEntityToID = GenerateEntityIDdict(EntityPandas)

    if len(args.pathways) > 0:
        LinkToPathways = args.pathways
        LinkToPathways = LinkToPathways.split(",")
        start_time = time.time()
        for PathwayName in LinkToPathways:
            PathwayMembersIdToEntity = PSnx.GetPathwayMemberIds([PathwayName],SearchPathwaysBy=['Name'],FilterBy=['Protein','FunctionalClass','Complex'], InProperty='objectType')
            EntityIDs = set(PathwayMembersIdToEntity.keys())
            QueryOntology = OQL.GetChildEntities(list(EntityIDs), ['id'])
            EntityIDs.update(PSnx.GetObjIds(QueryOntology))

            if len(EntityIDs) == 0:
                print ('No entity for %s found in the database' % PathwayName)
            else:
                AddSemanticRefCount(EntityPandas, PathwayName, EntityIDs)

        print("Entities in file %s were linked to %s pathway in %s" % (EntityListFile, PathwayName, ElsevierAPI.ExecutionTime(start_time)))


    if len(LinkToConcepts) > 0:
        start_time = time.time()
        if type(LinkToConcepts.index.name) == type(None):
            SearchConceptBy=['Name','Alias']
        else:
            SearchConceptBy=LinkToConcepts.index.name.split(',')

        EntityPandas.index.name = ','.join(SearchConceptBy)
        
            
        for LinkToConcept in LinkToConcepts.index:
            ConceptIDs = PSnx.GetObjIdsByProps(PropertyValues=[LinkToConcept],SearchByProperties=SearchConceptBy,OnlyObjectTypes=QueryObjType)
 
            if len(ConceptIDs) == 0: print ('No entity %s found in database' % LinkToConcept)
            else:
                effect = []
                reltypes = []
                dir = ''
                if 'Effect' in LinkToConcepts.columns:
                    ef = LinkToConcepts.at[LinkToConcept,'Effect']
                    if pd.notna(ef): effect = str(ef).split(',')

                if 'RelationType' in LinkToConcepts.columns:
                    reltypes = str(LinkToConcepts.at[LinkToConcept,'RelationType']).split(',')

                if 'Direction' in LinkToConcepts.columns:
                    dir = str(LinkToConcepts.at[LinkToConcept,'RelationType'])

                AddSemanticRefCount(EntityPandas,PandaEntityToID, LinkToConcept,list(ConceptIDs),reltypes,effect,dir,REL_PROP,ENT_PROPS)

        print("Entities in file %s were linked to %s in %s" % (EntityListFile, LinkToConcept, ElsevierAPI.ExecutionTime(start_time) ))

    foutName = EntityListFile[:len(EntityListFile)-4]+'+SemanticRefcount.tsv'
    EntityPandas.to_csv(foutName, sep='\t', index=True)
    if args.dump_references in ['True', 'true','yes','y','Y']:
        fOut = EntityListFile[:len(EntityListFile)-4]+'+SemanticReferences.tsv'
        PSnx.PrintReferenceView(fOut,REL_PROP,entPropNames=REL_PROP)

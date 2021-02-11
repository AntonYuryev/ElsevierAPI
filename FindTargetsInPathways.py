import time
import pandas as pd
import argparse
import textwrap
from ElsevierAPI import networx as PSnx

def MapEntityToPathways(IdToFolderPathwayDict:dict(),IDtoPathways:dict(),folderName,FilterBy:list,SearchByProperties=['Name','Alias']):
    IdtoMembersCollector = dict()
    for pthwy in IDtoPathways.values():
        if pthwy['ObjTypeName'][0] == 'Pathway':
            PathwayName = pthwy['Name'][0]
            PathwayIds = pthwy['Id']
            #FolderPathwayTuple = (folderName, PathwayName)
            IdtoMembersInPathway = PSnx.GetPathwayMemberIds(PathwayIds=PathwayIds, OnlyEntities=FilterBy, WithProperties=SearchByProperties)
            for EntityId in IdtoMembersInPathway.keys():
                try:
                    IdToFolderPathwayDict[EntityId][folderName].append(PathwayName)
                except KeyError:
                    try: 
                        IdToFolderPathwayDict[EntityId][folderName] = [PathwayName]
                    except KeyError:
                        IdToFolderPathwayDict[EntityId]= {folderName:[PathwayName]}

            IdtoMembersCollector.update(IdtoMembersInPathway)

    return IdtoMembersCollector

def FindPathways(FoldersWithPathways:list, EntityPandas, SearchByProperty=['Name','Alias']):
    Entities = list(EntityPandas.index)
    IdtoMembers = dict()

    for PathwayFolder in FoldersWithPathways:
        subFolders = PSnx.GetSubfolders([PathwayFolder])
        if len(subFolders) == 0:
            print ('Input folder has no subfolders. Will use only pathways from %s' % PathwayFolder)
        else:
            print('Found %d subfolders in %s' % (len(subFolders), PathwayFolder))

        IdToFolderPathwayDict = dict()
        ParentFolder = PSnx.IdToFolders[PathwayFolder][0]
        IDtoPathways = PSnx.GetObjectsFromFolders([ParentFolder.Id])
        IdtoMembers.update(MapEntityToPathways(IdToFolderPathwayDict,IDtoPathways,PathwayFolder,Entities,SearchByProperty))

        for FolderId in subFolders:
            subFolder = PSnx.IdToFolders[FolderId][0]
            subFolderName = subFolder['Name']
            allSubSubFolders = PSnx.GetSubfolderTree(FolderId)
                
            IDtoPathways = PSnx.GetObjectsFromFolders(allSubSubFolders)
            IdtoMembers.update(MapEntityToPathways(IdToFolderPathwayDict,IDtoPathways,subFolderName,Entities,SearchByProperty))

    return IdToFolderPathwayDict, IdtoMembers
    

if __name__ == "__main__":
    EntityListFile =str()
    instructions = '''
    infile - single column file with entity names. 
    search_property - property to find entities from infile. Default search property: Name,Alias
    pathway_folder - Comma-separated list of Folders in Pathway Studio server with pathway collection. 
    Pathways in all subfolders will be also used. Default folder: 'Hallmarks of Cancer'
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog=textwrap.dedent(instructions))
    parser.add_argument('-i', '--infile', type=str, required=True)
    parser.add_argument('-p', '--search_by_properties', default='Name,Alias')
    parser.add_argument('-f', '--pathway_folder', default='Hallmarks of Cancer')
    args = parser.parse_args()

    if len(args.infile) > 0:
        EntityListFile = args.infile #full path to single column entity list file
    else:
        print('No entity list file was specified')

    if len(args.search_by_properties) > 0:
        SearchByProperty = args.search_by_properties.split(",")


    SearchPathwaysInFolders = args.pathway_folder.split(",")
    
    start_time = time.time()
    EntityPandas = pd.read_csv(EntityListFile, delimiter='\t', header=0, index_col=0)
    IdToFolderPathwayDict, IdtoMembers = FindPathways(SearchPathwaysInFolders, EntityPandas)


    HallmarkColumnName = 'Hallmarks'
    PathwayColumnName = 'Pathways'
    EntityPandas.insert(len(EntityPandas.columns), HallmarkColumnName, '')
    EntityPandas.insert(len(EntityPandas.columns), PathwayColumnName, '')

    for EntityId, FolderToPathwayDict in IdToFolderPathwayDict.items():
        TargetName = IdtoMembers[EntityId]['Name'][0]
        if TargetName in EntityPandas.index:
            for hallmark, patwhays in FolderToPathwayDict.items():
                EntityPandas.at[TargetName,HallmarkColumnName] = EntityPandas.at[TargetName,HallmarkColumnName]+'\''+ hallmark+'\''+ ';'
                allPathwaysinHallmark = '\';\''.join(patwhays)
                EntityPandas.at[TargetName,PathwayColumnName] = EntityPandas.at[TargetName,PathwayColumnName]+'\''+allPathwaysinHallmark+'\';'


    foutName = EntityListFile[:len(EntityListFile)-4]+' Pathways from '+','.join(SearchPathwaysInFolders)+'.tsv'
    EntityPandas.to_csv(foutName, sep='\t', index=True)

    PathwayCounter = set([z for z in [x for x in [v.values() for v in IdToFolderPathwayDict.values()]]])
    print('Found %d pathways in \"%s\" for %d entities from \"%s\"' % (len(PathwayCounter), ','.join(SearchPathwaysInFolders), len(IdToFolderPathwayDict), EntityListFile))
    print('Results were printed to \"%s\" in %s seconds' % (foutName,time.time() - start_time))


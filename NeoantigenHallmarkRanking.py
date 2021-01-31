import time
import pandas as pd
import argparse
import textwrap
from ElsevierAPI import networx as PSnetworkx

def MapEntityToPathways(IdToFolderPathwayDict:dict(),IDtoPathways:dict(),folderName,FilterBy:list,SearchByProperties=['Name','Alias']):
    IdtoMembersCollector = dict()
    for pthwy in IDtoPathways.values():
        if pthwy['ObjTypeName'][0] == 'Pathway':
            PathwayName = pthwy['Name'][0]
            PathwayIds = pthwy['Id']
            #FolderPathwayTuple = (folderName, PathwayName)
            IdtoMembersInPathway = PSnetworkx.GetPathwayMemberIds(PathwayIds=PathwayIds, OnlyEntities=FilterBy, WithProperties=SearchByProperties)
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
        subFolders = PSnetworkx.GetSubfolders([PathwayFolder])
        if len(subFolders) == 0:
            print ('Input folder has no subfolders. Will use only pathways from %s' % PathwayFolder)
        else:
            print('Found %d subfolders in %s' %  (len(subFolders), PathwayFolder))

        IdToFolderPathwayDict = dict()
        ParentFolder = PSnetworkx.IdToFolders[PathwayFolder][0]
        IDtoPathways = PSnetworkx.GetObjectsFromFolders([ParentFolder.Id])
        IdtoMembers.update(MapEntityToPathways(IdToFolderPathwayDict,IDtoPathways,PathwayFolder,Entities,SearchByProperty))

        for FolderId in subFolders:
            subFolder = PSnetworkx.IdToFolders[FolderId][0]
            subFolderName = subFolder['Name']
            allSubSubFolders = PSnetworkx.GetSubfolderTree(FolderId)
                
            IDtoPathways = PSnetworkx.GetObjectsFromFolders(allSubSubFolders)
            IdtoMembers.update(MapEntityToPathways(IdToFolderPathwayDict,IDtoPathways,subFolderName,Entities,SearchByProperty))

    return IdToFolderPathwayDict, IdtoMembers
    

if __name__ == "__main__":
    EntityListFile =str()
    instructions = '''
    infile - single column file with entity names. 
    
    search_property - property to find entities from infile. Default search property: Name
    
    Semantic reference count to a pathway is calculated as sum of reference count of semantic links to every pathay component
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog=textwrap.dedent(instructions))
    parser.add_argument('-i', '--infile', type=str, required=True)
    parser.add_argument('-p', '--search_by_properties', default='Name,Alias')
    args = parser.parse_args()

    if len(args.infile) > 0:
        EntityListFile = args.infile #full path to single column entity list file
    else:
        print('No entity list file was specified')

    if len(args.search_by_properties) > 0:
        SearchByProperty = args.search_by_properties.split(",")


    SearchPathwaysInFolders = ['Hallmarks of Cancer']
    Neo7HallmarkRanking = {
        'Hallmarks of Cancer (1): Sustaining Proliferative Signaling':14, 
        'Hallmarks of Cancer (6): Activating Invasion and Metastasis':13,
        'Hallmarks of Cancer (5): Inducing Angiogenesis':12,
        'Hallmarks of Cancer (2): Evading Growth Suppressors':11,
        'Hallmarks of Cancer (10): Tumor-Promoting Inflammation':10,
        'Hallmarks of Cancer (8): Evading Immune Destruction':9,
        'Hallmarks of Cancer (3): Resisting Cell Death':8,
        'Hallmarks of Cancer (7): Deregulated Metabolism':7,
        'Hallmarks of Cancer (11): Mitochondria Instability':6,
        'Hallmarks of Cancer (12): Tumor Microenvironment Physics Adaptation':5,
        'Hallmarks of Cancer (4): Enabling Replicative Immortality':4,
        'Hallmarks of Cancer (9): Genome Instability':3,
        'Hallmarks of Cancer (13): Drug resistance':2,
        'Hallmarks of Cancer (14): Cancer Epigenetics':1
    }

    start_time = time.time()
    EntityPandas = pd.read_csv(EntityListFile, delimiter='\t', header=0, index_col=0)
    IdToFolderPathwayDict, IdtoMembers = FindPathways(SearchPathwaysInFolders, EntityPandas)

    HallmarkRankColumnName = 'Hallmark rank'
    HallmarkColumnName = 'Hallmark'
    PathwayColumnName = 'Pathways'
    EntityPandas.insert(len(EntityPandas.columns), HallmarkRankColumnName, 0)
    EntityPandas.insert(len(EntityPandas.columns), HallmarkColumnName, '')
    EntityPandas.insert(len(EntityPandas.columns), PathwayColumnName, '')

    for EntityId, FolderToPathwayDict in IdToFolderPathwayDict.items():
        TargetName = IdtoMembers[EntityId]['Name'][0]
        for hallmark, patwhays in FolderToPathwayDict.items():
            rank = Neo7HallmarkRanking[hallmark]
            RankinPandas = EntityPandas.at[TargetName,HallmarkRankColumnName]
            if rank > RankinPandas:
                EntityPandas.at[TargetName,HallmarkRankColumnName] = rank
                EntityPandas.at[TargetName,HallmarkColumnName] = hallmark
                allPathwaysinHallmark = '\',\''.join(patwhays)
                EntityPandas.at[TargetName,PathwayColumnName] = '\''+allPathwaysinHallmark+'\''
        

    foutName = EntityListFile[:len(EntityListFile)-4]+' Pathways from '+','.join(SearchPathwaysInFolders)+'.tsv'
    EntityPandas.to_csv(foutName, sep='\t', index=True)

    PathwayCounter = set([z for z in [x for x in [v.values() for v in IdToFolderPathwayDict.values()]]])
    print('Found %d pathways in \"%s\" for %d entities from \"%s\"' % (len(PathwayCounter), ','.join(SearchPathwaysInFolders), len(IdToFolderPathwayDict), EntityListFile))
    print('Results were printed to \"%s\" in %s seconds' % (foutName,time.time() - start_time))


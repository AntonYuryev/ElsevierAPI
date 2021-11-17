import sys
import time
from ElsevierAPI import open_api_session

ps_api = open_api_session()

def MapEntityToPathways(EntityDict:dict(), IDtoPathways:dict(),FilterBy:list, folderName):
    for id, pthwy in IDtoPathways.items():
        if pthwy['ObjTypeName'][0] == 'Pathway':
            PathwayName = pthwy['Name'][0]
            PathwayIds = pthwy['Id']
            EntityDictValue = folderName + '\t' + PathwayName
            IdtoMembers = ps_api.get_pathway_member_ids(PathwayIds, FilterBy, InProperties=SearchByProperty)
            for ent in IdtoMembers.values():
                EntityName = ent['Name'][0]
                EntityURN = ent['URN'][0]
                EntityDictKey = EntityName+'\t'+EntityURN
                if EntityDictKey in EntityDict.keys():
                    EntityDict[EntityDictKey].update([EntityDictValue])
                else:
                    EntityDict[EntityDictKey] = set([EntityDictValue])

def FindPathways(FoldersWithPathways:list, EntityListFile):
    Entities = list()
    with open(EntityListFile) as file:
        Entities = [line.strip() for line in file]

    for PathwayFolder in FoldersWithPathways:
        subFolders = ps_api.get_subfolders([PathwayFolder])
        if len(subFolders) == 0:
            print ('Input folder has no subfolders. Will use only pathways from %s' % PathwayFolder)
        else:
            print('Found %d subfolders in %s' %  (len(subFolders), PathwayFolder))

        EntityDict = dict()
        ParentFolder = ps_api.IdToFolders[PathwayFolder][0]
        IDtoPathways = ps_api.get_objects_from_folders([ParentFolder.Id])
        MapEntityToPathways(EntityDict, IDtoPathways, Entities, PathwayFolder)

        for FolderId in subFolders:
            subFolder = ps_api.IdToFolders[FolderId][0]
            subFolderName = subFolder['Name']
            allSubSubFolders = ps_api.get_subfolders_tree(FolderId)
                
            IDtoPathways = ps_api.get_objects_from_folders(allSubSubFolders)
            MapEntityToPathways(EntityDict, IDtoPathways, Entities, subFolderName)

    return dict(sorted(EntityDict.items()))
    
if __name__ == "__main__":
    if len(sys.argv) > 1:
        EntityListFile = sys.argv[1] #full path to single column entity list file
    else:
        print('No entity list file was specified')
        EntityListFile = 'D:\\Python\\NEO7\\NEO7US208CNMLN008\\NEO7US208CNMLN008_Neoantigens.txt'

    if len(sys.argv) > 2:
        SearchByProperty = sys.argv[2] 
        SearchByProperty = SearchByProperty.split(",")
    else:
        SearchByProperty = ['Name', 'Alias']
        print('No search property to find entities from \"%s\" was specified. Will use %s' % (EntityListFile,','.join(SearchByProperty)))

    if len(sys.argv) > 3:
        SearchPathwaysInFolders = sys.argv[3] 
        SearchPathwaysInFolders = SearchPathwaysInFolders.split(",")
    else:
        SearchPathwaysInFolders = ['Hallmarks of Cancer']
        print('No folder list was specified. Will use pathways from \"%s\" folder ' % ','.join(SearchPathwaysInFolders))


    start_time = time.time()
    sorted_dict = FindPathways(SearchPathwaysInFolders, EntityListFile)
    foutName = EntityListFile[:len(EntityListFile)-4]+' Pathways from '+','.join(SearchPathwaysInFolders)+'.tsv'
    with open(foutName, 'w', encoding='utf-8') as fout:
        for k,v in sorted_dict.items():
            for p in v:
                fout.write(k+'\t'+p+'\n')

    PathwayCounter = {x for x in v for v in sorted_dict.values()}
    print('Found %d pathways in \"%s\" for %d entities from \"%s\"' % (len(PathwayCounter), ','.join(SearchPathwaysInFolders), len(sorted_dict), EntityListFile))
    print('Results were printed to \"%s\" in %s seconds' % (foutName,time.time() - start_time))



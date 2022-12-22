import time
import pandas as pd
import argparse
import textwrap
from ElsevierAPI import load_api_config
from ElsevierAPI.ResnetAPI.FolderContent import FolderContent

ps_api = FolderContent(load_api_config())

def map2pathways(entity2folder2pathway: dict, id2pathways: dict, folderName, FilterBy: list,
                 SearchByProperties=None):
    SearchByProperties = ['Name', 'Alias'] if SearchByProperties is None else SearchByProperties
    IdtoMembersCollector = dict()
    for pathway in id2pathways.values():
        if pathway['ObjTypeName'][0] == 'Pathway':
            PathwayName = pathway['Name'][0]
            PathwayIds = pathway['Id']
            IdtoMembersInPathway = ps_api.get_pathway_members(PathwayIds,FilterBy,SearchByProperties)
            for entity_Id in IdtoMembersInPathway.keys():
                try:
                    entity2folder2pathway[entity_Id][folderName].append(PathwayName)
                except KeyError:
                    try:
                        entity2folder2pathway[entity_Id][folderName] = [PathwayName]
                    except KeyError:
                        entity2folder2pathway[entity_Id] = {folderName: [PathwayName]}

            IdtoMembersCollector.update(IdtoMembersInPathway)

    return IdtoMembersCollector


def find_pathways(FoldersWithPathways: list, entity_pandas: pd.DataFrame, search_by_property=None):
    search_by_property = ['Name', 'Alias'] if search_by_property is None else search_by_property

    Entities = list(entity_pandas.index)
    id2members = dict()

    for PathwayFolder in FoldersWithPathways:
        subFolders = ps_api.get_subfolders([PathwayFolder])
        if len(subFolders) == 0:
            print('Input folder has no subfolders. Will use only pathways from %s' % PathwayFolder)
        else:
            print('Found %d subfolders in %s' % (len(subFolders), PathwayFolder))

        id2folder_pathway = dict()
        ParentFolder = ps_api.id2folder[PathwayFolder][0]
        id2pathways = ps_api.get_objects_from_folders([ParentFolder.Id])
        id2members.update(
            map2pathways(id2folder_pathway, id2pathways, PathwayFolder, Entities, search_by_property))

        for FolderId in subFolders:
            subFolder = ps_api.id2folders[FolderId][0]
            subFolderName = subFolder['Name']
            allSubSubFolders = ps_api.get_subfolder_tree(subFolderName)

            # {PathwayId:psObject} pathway ID to properties dict
            id2pathways = ps_api.get_objects_from_folders(allSubSubFolders)
            id2members.update(
                map2pathways(id2folder_pathway, id2pathways, subFolderName, Entities, search_by_property))

    return id2folder_pathway, id2members


if __name__ == "__main__":
    EntityListFile = str()
    instructions = '''
    infile - single column file with entity names. 
    search_property - property to find entities from infile. Default search property: Name,Alias
    pathway_folder - Comma-separated list of Folders in Pathway Studio server with pathway collection. 
    Pathways in all subfolders will be also used. Default folder: 'Hallmarks of Cancer'
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent(instructions))
    parser.add_argument('-i', '--infile', type=str, required=True)
    parser.add_argument('-p', '--search_by_properties', default='Name,Alias')
    parser.add_argument('-f', '--pathway_folder', default='Hallmarks of Cancer')
    args = parser.parse_args()

    if len(args.infile) > 0:
        EntityListFile = args.infile  # full path to single column entity list file
    else:
        print('No entity list file was specified')

    if len(args.search_by_properties) > 0:
        SearchByProperty = args.search_by_properties.split(",")

    SearchPathwaysInFolders = args.pathway_folder.split(",")

    start_time = time.time()
    EntityPandas = pd.read_csv(EntityListFile, delimiter='\t', header=0, index_col=0)
    IdToFolderPathwayDict, IdtoMembers = find_pathways(SearchPathwaysInFolders, EntityPandas)

    HallmarkColumnName = 'Hallmarks'
    PathwayColumnName = 'Pathways'
    EntityPandas.insert(len(EntityPandas.columns), HallmarkColumnName, '')
    EntityPandas.insert(len(EntityPandas.columns), PathwayColumnName, '')

    for EntityId, FolderToPathwayDict in IdToFolderPathwayDict.items():
        TargetName = IdtoMembers[EntityId]['Name'][0]
        if TargetName in EntityPandas.index:
            for hallmark, pathways in FolderToPathwayDict.items():
                EntityPandas.at[TargetName, HallmarkColumnName] = EntityPandas.at[
                                                                      TargetName, HallmarkColumnName] + '\'' + hallmark + '\'' + ';'
                allPathwaysInHallmark = '\';\''.join(pathways)
                EntityPandas.at[TargetName, PathwayColumnName] = EntityPandas.at[
                                                                     TargetName, PathwayColumnName] + '\'' + allPathwaysInHallmark + '\';'

    foutName = EntityListFile[:len(EntityListFile) - 4] + ' Pathways from ' + ','.join(SearchPathwaysInFolders) + '.tsv'
    EntityPandas.to_csv(foutName, sep='\t', index=True)

    PathwayCounter = {z for z in [x for x in [v.values() for v in IdToFolderPathwayDict.values()]]}
    print('Found %d pathways in \"%s\" for %d entities from \"%s\"' %
         (len(PathwayCounter), ','.join(SearchPathwaysInFolders), len(IdToFolderPathwayDict), EntityListFile))
    print('Results were printed to \"%s\" in %s seconds' % (foutName, time.time() - start_time))

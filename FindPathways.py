import sys
import time
from ElsevierAPI import load_api_config
from ElsevierAPI.ResnetAPI.FolderContent import FolderContent

ps_api = FolderContent(load_api_config())

def MapEntityToPathways(FilterBy:list):
    EntityDict = dict()
    for id, pthwy in ps_api.id2pathways.items():
        folder_name = ps_api.id2folders[id][0]['Name']
        if pthwy['ObjTypeName'][0] == 'Pathway':
            pathway_name = pthwy['Name'][0]
            PathwayIds = pthwy['Id']
            EntityDictValue = folder_name + '\t' + pathway_name
            IdtoMembers = ps_api.get_pathway_member_ids(PathwayIds, FilterBy, InProperties=SearchByProperty)
            for ent in IdtoMembers.values():
                EntityName = ent['Name'][0]
                EntityURN = ent['URN'][0]
                EntityDictKey = EntityName+'\t'+EntityURN
                try:
                    EntityDict[EntityDictKey].update([EntityDictValue])
                except KeyError:
                    EntityDict[EntityDictKey] = set([EntityDictValue])

    return EntityDict


def FindPathways(FoldersWithPathways:list, EntityListFile:str):
    # loading id2pathways, id2folders into ps_api
    for folder in FoldersWithPathways:
        folder_id = ps_api.get_folder_id(folder)
        sub_folder_ids = ps_api.get_subfolders_recursively([folder_id])
        if len(sub_folder_ids) == 0:
            print ('Input folder has no subfolders. Will use only pathways from %s' % folder)
        else:
            print('Found %d subfolders in %s' %  (len(sub_folder_ids), folder))

        ps_api.get_objects_from_folders([sub_folder_ids])

    # making {EntityName+'\t'+EntityURN : folderName + '\t' + PathwayName}
    Entities = [line.strip() for line in open(EntityListFile)]
    EntityDict= MapEntityToPathways(EntityDict, Entities)

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
    entities = [line.strip() for line in open(EntityListFile)]
    entity_ids = ps_api._get_obj_ids_by_props(entities,SearchByProperty)
    id2psobj = ps_api.find_pathways(entity_ids,SearchPathwaysInFolders)
    pathway_counter = set()
    # making {EntityName+'\t'+EntityURN : folderName + '\t' + PathwayName}
    foutName = EntityListFile[:len(EntityListFile)-4]+' Pathways from '+','.join(SearchPathwaysInFolders)+'.tsv'
    rows = list()
    for psobj in id2psobj.values():
        entity_name = psobj['Name'][0]
        entity_urn = psobj['URN'][0]
        for pathway_id in psobj['Pathway ID']:
            patwhay_obj = ps_api.id2pathways[pathway_id]
            pathway_name = patwhay_obj['Name'][0]
            pathway_counter.add(pathway_name)
            for folder_name in patwhay_obj['Folders']:
                rows.append(entity_name+'\t'+entity_urn+'\t'+folder_name+'\t'+pathway_name+'\n')

    rows.sort()
    with open(foutName, 'w', encoding='utf-8') as fout:
        [fout.write(row) for row in rows]

    print('Found %d pathways in \"%s\" for %d entities from \"%s\"' % (len(pathway_counter), ','.join(SearchPathwaysInFolders), len(entities), EntityListFile))
    print('Results were printed to \"%s\" in %s seconds' % (foutName,time.time() - start_time))



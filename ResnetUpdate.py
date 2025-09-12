from ENTELLECT_API.ElsevierAPI.ResnetAPI.ResnetGraph import ResnetGraph,PSRelation
from ENTELLECT_API.ElsevierAPI.ResnetAPI.NetworkxObjects import ALL_PSREL_PROPS
from ENTELLECT_API.ElsevierAPI.utils import dir2flist
from os import path

UPDATE_YEAR = 2024
DATA_DIR = 'C:/ResnetDump/BiomarkerTest'
UPDATE_DIR = 'C:/ResnetDump'
#RESNET_DUMP = path.join(DATA_DIR,'resnet18_dump_05062025.rnef')
RESNET_UPDATE =  path.join(UPDATE_DIR,f'resnet_since{UPDATE_YEAR}_')
MAX_REL_COUNT = 40000

def update_resnet(resnet_dump:str,year=UPDATE_YEAR,file_count=0):
    """
    Generates the ResNet graph with the latest data from resnet_dump
    :param year: The year of the update.
    :param resnet_dump: The path to the ResNet dump file.
    """
    update_rels = list()
    dump_files_count = file_count
    for nodes,rels in ResnetGraph.read_rnef(resnet_dump):
      for rel in rels:
        #assert(isinstance(rel,PSRelation))
        refs = rel.refs()
        update_refs = [r for r in refs if r.pubyear() >= year]
        if update_refs:
          rel.references = update_refs
          update_rels.append(rel)
        
      if len(update_rels) >= MAX_REL_COUNT:
        dump_files_count += 1
        updateG = ResnetGraph.from_rels(update_rels)
        updateG = updateG.remove_undirected_duplicates()
        updateG.dump2rnef(f'{RESNET_UPDATE}{dump_files_count}.rnef',rel_prop2print=ALL_PSREL_PROPS,with_section_size=1000)
        update_rels.clear()
        updateG.clear_resnetgraph()

    if update_rels:
      updateG = ResnetGraph.from_rels(update_rels)
      updateG = updateG.remove_undirected_duplicates()
      dump_files_count += 1
      updateG.dump2rnef(f'{RESNET_UPDATE}{dump_files_count}.rnef',rel_prop2print=ALL_PSREL_PROPS,with_section_size=1000)
      updateG.clear_resnetgraph()
    return dump_files_count, len(update_rels) + MAX_REL_COUNT*(dump_files_count-1)


if __name__ == "__main__":
    dump_chunks = dir2flist(DATA_DIR, 'rnef')
    total_edges = 0
    dump_files_count = 0
    print(f"Starting ResNet update for year {UPDATE_YEAR} with {len(dump_chunks)} chunks in {DATA_DIR}.")
    for chunk in dump_chunks:
      dump_files_count, chunk_edges = update_resnet(chunk, UPDATE_YEAR,dump_files_count)
      total_edges += chunk_edges
      print(f"Processed {chunk} with {chunk_edges} edges, total edges so far: {total_edges}")
    print(f"ResNet update completed for year {UPDATE_YEAR}.")
    print(f"{total_edges} relations dumped to {RESNET_UPDATE}")
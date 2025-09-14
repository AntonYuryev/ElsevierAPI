from ENTELLECT_API.ElsevierAPI.ResnetAPI.ResnetGraph import ResnetGraph
from ENTELLECT_API.ElsevierAPI.ResnetAPI.NetworkxObjects import ALL_PSREL_PROPS
from ENTELLECT_API.ElsevierAPI.utils import dir2flist, Tee, time, execution_time
from os import path

UPDATE_YEAR = 2025
DATA_DIR = 'C:/ResnetDump/OneDump'
UPDATE_DIR = 'C:/ResnetDump/ResnetClean/Update2025'
RESNET_UPDATE =  path.join(UPDATE_DIR,f'resnet_since{UPDATE_YEAR}')
MAX_REL_COUNT = 40000

my_relprops = ALL_PSREL_PROPS + ["TextMods"]
my_node_props = [
        "Name",
        "Alias",
        "MeSH ID",
        "MeSH Heading",
        "MedScan ID",
        "Hugo ID",
        "Description",
        "Notes",
        "Primary Cell Localization",
        "Organism",
        "Organ"
    ]


def update_resnet(resnet_dump:str,year=UPDATE_YEAR,file_count=0):
    """
    Generates the ResNet graph with the latest data from resnet_dump
    :param year: The year of the update.
    :param resnet_dump: The path to the ResNet dump file.
    """
    update_rels = list()
    dump_files_count = file_count
    all_nodes = dict()
    total_rels_processed = 0
    resnet_counter = 0
    for nodes,rels in ResnetGraph.read_rnef(resnet_dump,only4objs=all_nodes):
      total_rels_processed += len(rels)
      resnet_counter += 1
      was_updated = False
      for rel in rels:
        refs = rel.refs()
        update_refs = [r for r in refs if r.pubyear() >= year]
        if update_refs:
          rel.references = update_refs
          update_rels.append(rel)
          was_updated = True

      if resnet_counter % 1000 == 0:
        print(f"Processed {resnet_counter} resnet sections with {len(all_nodes)} nodes and {total_rels_processed} relations")
      
      if update_rels:
        rels4update_count = len(update_rels)
        if rels4update_count % 1000 == 0 and was_updated:
          print(f"Collected {rels4update_count} relations out of {total_rels_processed} processed for {year} update from {resnet_dump}")

      if len(update_rels) >= MAX_REL_COUNT:
        dump_files_count += 1
        updateG = ResnetGraph.from_rels(update_rels)
        updateG = updateG.remove_undirected_duplicates()
        updateG.dump2rnef(f'{RESNET_UPDATE}{dump_files_count}.rnef',ent_prop2print=my_node_props,rel_prop2print=my_relprops,with_section_size=1000)
        update_rels.clear()
        updateG.clear_resnetgraph()

    if update_rels:
      updateG = ResnetGraph.from_rels(update_rels)
      updateG = updateG.remove_undirected_duplicates()
      dump_files_count += 1
      updateG.dump2rnef(f'{RESNET_UPDATE}{dump_files_count}.rnef',ent_prop2print=my_node_props,rel_prop2print=my_relprops,with_section_size=1000)
      updateG.clear_resnetgraph()
    return dump_files_count, len(update_rels) + MAX_REL_COUNT*(dump_files_count-1)


if __name__ == "__main__":
    dump_files = dir2flist(DATA_DIR, 'rnef')
    total_edges = 0
    dump_files_count = 0
    start = time.time()
    with Tee(path.join(UPDATE_DIR,f'Resnet{UPDATE_YEAR}Update.log'), 'w') as log:
      print(f'Update generation started at {time.ctime()}')
      print(f"Starting ResNet update for year {UPDATE_YEAR} with {len(dump_files)} files in {DATA_DIR}.")
      for file in dump_files:
        dump_files_count, chunk_edges = update_resnet(file, UPDATE_YEAR,dump_files_count)
        total_edges += chunk_edges
        print(f"Processed {file} with {chunk_edges} edges, total edges so far: {total_edges}")
      print(f"ResNet update completed for year {UPDATE_YEAR}.")
      print(f"{total_edges} relations dumped to {RESNET_UPDATE}")
      print(f'Update generation finished at {time.ctime()}')
      print(f"Update was generated in: {execution_time(start)}")
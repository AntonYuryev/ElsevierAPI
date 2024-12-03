from ENTELLECT_API.ElsevierAPI.ResnetAPI.ResnetGraph import ResnetGraph,df
from ENTELLECT_API.ElsevierAPI.utils import dir2flist,str2str,Tee
import os,csv,argparse,textwrap

DATA_DIR = 'C:/ResnetDump/ResnetClean/'

if __name__ == "__main__":
  instructions = '''
    indir - full path to directory with .rnef files. Script will also use all .rnef files from all subdirectories
    outdir - full path to directory for writing output
    db - specifies format for the output.  Options: [neo4j, neptune]
    '''

  parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog=textwrap.dedent(instructions))
  parser.add_argument('-i', '--indir', type=str)
  parser.add_argument('-o', '--outdir', type=str)
  parser.add_argument('-d', "--db", choices=['neo4j', 'neptune'], default='neo4j',help="Select output database")

  args = parser.parse_args()
  resnet_dump = args.indir
  dump_name,_ = os.path.splitext(os.path.basename(resnet_dump))

  ext = '.txt'
  csv_kwargs = {'sep':'|','index':False}
  
  if args.db == 'neptune':
    csv_kwargs['sep'] = ','
    csv_kwargs.update({'quotechar':'"','quoting':csv.QUOTE_NONNUMERIC})
    ext = '.csv'

  refset_file = os.path.join(args.outdir,dump_name+f'.RefSets.{args.db}'+ext)
  open(refset_file,'w').close() # flash
  rel_file =  os.path.join(args.outdir,dump_name+f'.Relations.{args.db}'+ext)
  open(rel_file,'w').close() # flash
  nodes_file = os.path.join(args.outdir,dump_name+f'.Nodes.{args.db}'+ext)
  refs_file = os.path.join(dump_name+f'.References.{args.db}'+ext)

  all_nodes = set()
  all_refs = dict()
  dir_files = dir2flist(args.indir,file_ext='.rnef')
  with Tee(os.path.join(args.outdir,resnet_dump+f'2{args.db}.log')):
    for i,file in enumerate(dir_files):
      print(f'Reading {os.path.basename(file)} file {i+1} out of {len(dir_files)}')
      G = ResnetGraph.fromRNEF(file)
      G = G.remove_undirected_duplicates()
      nodes, refset_df,rel_df = G.neo4j_df()
      G.merge_refs(all_refs)

      len_before = len(all_nodes)
      all_nodes.update(nodes)
      all_nodes_len = len(all_nodes)
      print(f'Added {all_nodes_len-len_before} new nodes from {len(nodes)}')
      print(f'all_nodes has {all_nodes_len} nodes in total')

      rel_df.to_csv(**csv_kwargs, mode='a',path_or_buf=rel_file)
      refset_df.to_csv(**csv_kwargs, mode='a',path_or_buf=refset_file)
      G.clear_resnetgraph()
   
  nodes_df = df([str2str(o) for o in all_nodes])
  nodes_df.to_csv(**csv_kwargs, mode='w',path_or_buf=nodes_file)

  all_refs = list(set(all_refs.values()))
  refs_df = df([str2str(o) for o in all_refs])
  refs_df.to_csv(**csv_kwargs, mode='w',path_or_buf=refs_file)
  
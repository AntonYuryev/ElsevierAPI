import argparse
import textwrap
from ElsevierAPI.ResnetAPI.FolderContent import FolderContent
from ElsevierAPI.ResnetAPI.FolderContent import SNIPPET_PROPERTIES, ONLY_REL_PROPERTIES,ENTITY_IDENTIFIERS
from ElsevierAPI import load_api_config 
import csv
from contextlib import redirect_stdout


if __name__ == "__main__":
    instructions = '''
    infile - single column file with pathway URNs
    folder - folder in Pathway Studio server to download pathways from. Output is in RNEF
    dir - working directory to download folder content. Defaults to current
    format - output xml format. Choices: RNEF, SBGN. Defaults to RNEF.
    resume_from - last downloaded folder to resume download
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog=textwrap.dedent(instructions))
    parser.add_argument('-i', '--infile', type=str, default='')
    parser.add_argument('-f', '--folder', type=str, default='')
    parser.add_argument('-d', '--dir', type=str, default='')
    parser.add_argument('-x', '--format', type=str, default='RNEF', choices=['RNEF','SBGN','JSON-LD','json'])
    parser.add_argument('-r', '--resume_from', type=str, default='')
    args = parser.parse_args()

    api_cofig_file = 'D:/Python/ENTELLECT_API/ElsevierAPI/APIconfigMDACC.json'
    ps_api = FolderContent(load_api_config(api_cofig_file),what2retrieve=SNIPPET_PROPERTIES)
    ps_api.add_ent_props(ENTITY_IDENTIFIERS+['Alias'])
    # set what2retrieve=ONLY_REL_PROPERTIES to exclude supporting sentecnes from download
    ps_api.set_dir(args.dir)

    if args.infile:
        urnList = [u[0] for u in csv.reader(open(args.infile,"r"), delimiter="\t")]
        print('Attempting to download %s pathways from %s' %(len(urnList),args.infile))
        urnlistfile = args.infile[:len(args.infile)-4]
        fout_name = '{urnlist} pathways'.format(urnlist=urnlistfile)
        from_folders = [args.folder] if args.folder else []

        with open('download_by_pathway_urns.log', 'w') as fout:
            with redirect_stdout(fout):
                ps_api.download_pathways_by_urn(urnList,fout_name,args.format,args.folder)
    elif args.folder:
            ps_api.reference_cache_size = 800000
            #ps_api.Graph is used to cache relations for speed but can be too big when a lot of pathways are downloaded
            # reference_cache_size controls the memory use by specifying total number of references allowed in ps_api.Graph
            # The smaller max_size the slower is download for large pathway collections 
            logfile = ps_api.data_dir+args.folder+'_download.log' 
            if args.resume_from:
                with open(logfile, 'a',encoding='utf-8') as f:
                    with redirect_stdout(f):
                        ps_api.resume_download(args.folder, args.resume_from)
            else:
                with open(logfile, 'w',encoding='utf-8') as f:
                    with redirect_stdout(f):
                        ps_api.folder2rnef(args.folder)
                        #results will be written into the RNEF XML file called "content of args.folder" 
                        #results will contain pathways and groups from the specified folder
                        #command performing same operation is available in Pathway UI from folder context menu accesible by right-clicking on the folder 



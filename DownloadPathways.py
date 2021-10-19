import time
import argparse
import textwrap
from ElsevierAPI import load_api_config 
from ElsevierAPI.ResnetAPI.FolderContent import FolderContent
import csv
from contextlib import redirect_stdout

if __name__ == "__main__":
    instructions = '''
    infile - single column file with pathway URNs
    format - output xml format. Choices: RNEF, SBGN. Defaults to RNEF.
    folder - folder in Pathway Studio server to download pathways from. Output is in RNEF
    resume_from - last downloaded folder to resume download
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog=textwrap.dedent(instructions))
    parser.add_argument('-i', '--infile', type=str, default='')
    parser.add_argument('-x', '--format', type=str, default='RNEF', choices=['RNEF', 'SBGN'])
    parser.add_argument('-f', '--folder', type=str, default='')
    parser.add_argument('-r', '--resume_from', type=str, default='')
    args = parser.parse_args()

    ps_api = FolderContent(load_api_config())

    if args.infile:
        urnList = [u[0] for u in csv.reader(open(args.infile,"r"), delimiter="\t")]
        print('Attempting to download %s pathways from %s' %(len(urnList),args.infile))
        urnlistfile = args.infile[:len(args.infile)-4]
        fout_name = '{urnlist} pathways.rnef'.format(urnlist=urnlistfile)

        with open('download_by_pathway_urns.log', 'w') as fout:
            with redirect_stdout(fout):
                urn2pathway = ps_api.load_containers()#works 40 sec - cache result to your application
                ps_api.urns2rnef(urnList,fout_name,args.format)


    if args.folder:
        ps_api.reference_cache_size = 800000
        logfile = args.folder+'_download.log'
        #ps_api.Graph is used to cache relations for speed but can be too big when a lot of pathways are downloaded
        # reference_cache_size controls the memory use. The smaller max_size the slower is download for large pathway collections  
        if args.resume_from:
            with open(logfile, 'a',encoding='utf-8') as f:
                with redirect_stdout(f):
                    ps_api.resume_download(args.folder, args.resume_from)
        else:
            with open(logfile, 'w',encoding='utf-8') as f:
                with redirect_stdout(f):
                        # by default relations will have only attribute from:
                        # {'TexRef'}|NetworkObjects.RELATION_PROPS|NetworkObjects.REF_ID_TYPES
                        # Uncomment line below to download pathways with Relations annotated with Sentence
                        #ps_api.add_rel_props(['Sentence'])
                        # add names of other attributes to input list for addtional Relation annotation
                        #max number of references in the ps_api.Graph until clearing
                        ps_api.content2rnef(args.folder)
                        #results will be written into the RNEF XML file called "content of args.folder" 
                        #results will contain pathways and groups from the specified folder
                        #command performing same operation is available in Pathway UI from folder context menu accesible by right-clicking on the folder 

                   






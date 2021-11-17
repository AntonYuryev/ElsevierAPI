import time
import ElsevierAPI.ResnetAPI.PathwayStudioGOQL as OQL
from ElsevierAPI import open_api_session,load_api_config
import ElsevierAPI.ReaxysAPI.Reaxys_API as RxAPI
import argparse
import textwrap


start_time = time.time()
def GOQLtoFindDrugs(TargetIds:list, TargetType = 'Protein', drugEffect=['negative']):
    if TargetType == 'Protein':
        return OQL.get_drugs(for_targets_with_ids=TargetIds)
    elif TargetType == 'Small Molecule':
        REL_TYPES = ['Regulation', 'MolSynthesis']
        OQLquery = OQL.expand_entity(PropertyValues=TargetIds, SearchByProperties=['id'], expand_by_rel_types=REL_TYPES, expand2neighbors=['Small Molecule'], direction='upstream')
        OQLquery += ' AND Effect = (' + ','.join(drugEffect)+')'
        return OQLquery
    else: 
        REL_TYPES = ['Regulation']
        OQLquery = OQL.expand_entity(PropertyValues=TargetIds, SearchByProperties=['id'], expand_by_rel_types=REL_TYPES, expand2neighbors=['Small Molecule'], direction='upstream')
        OQLquery += ' AND Effect = (' + ','.join(drugEffect)+')'
        return OQLquery


if __name__ == "__main__":   
    instructions = '''
    infile - single column file with entity names that must be modulated by drugs. 
    If you want to use identifiers other than names enter appropriate Propeprty types into SearchByProperties list.

    target_type - indicates type of entities are in infile: Protein, Small Molecule, CellProcess, Disease, ClinicalParameter, etc. Only one type of entities is allowed in infile

    resnet_props - comma-separated list of property names used for searching entities from infile in Resnet. Default: Name+Alias
    effect - Find Agonists (effect=positive) or Antagonits (effect=negative). Default: negative

    reaxys_prop - comma-separated list of Reaxys properties to annotate drugs found in Resnet.
    Script finds drugs affecting entities in infile and then annotates them with Reaxys properties for output
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog=textwrap.dedent(instructions))
    parser.add_argument('-i', '--infile', type=str, required=True)
    parser.add_argument('-t', '--target_type', type=str, required=True, default='Disease')
    parser.add_argument('-p', '--resnet_search_props', type=str, default='Name,Alias')
    parser.add_argument('-a', '--resnet_retreive_props', type=str)
    parser.add_argument('-e', '--effect', type=str, default='negative')
    parser.add_argument('-r', '--reaxys_prop', type=str, default='')
    parser.add_argument('--debug', action="store_true")
    args = parser.parse_args()


    TargetType = args.target_type 
    SearchPsProps = args.resnet_search_props.split(',')
    drugEffect = [args.effect]
  
    if args.reaxys_prop:
        ReaxysFields = args.reaxys_prop.split(',')
    else:
        print ('No Reaxys properties specified!')
        ReaxysFields = []

    with open(str(args.infile)) as f:
        EntitiesToExpand = [line.rstrip('\n') for line in f]

    ps_api = open_api_session()
    TargetIDs = ps_api._get_obj_ids_by_props(PropertyValues=EntitiesToExpand, SearchByProperties=SearchPsProps)
    PSdumpFile = str(args.infile)[:len(str(args.infile))-4]+'_psdump.tsv'
    ps_api.add_dump_file(PSdumpFile, replace_main_dump=True)
    ps_api.entProps = args.resnet_retreive_props.split(',')
    ps_api.add_ent_props('Reaxys ID', 'InChIKey')
  
    ps_api.relProps = ['Name','Sentence','PMID','DOI','PubYear','RelationNumberOfReferences']
    #Data dump columns will be ordered according to the order in this list
    ps_api.process_oql(GOQLtoFindDrugs(TargetIDs, TargetType=TargetType, drugEffect=drugEffect))

    if len(ReaxysFields) > 0:
        FoundDrugs = [y for x,y in ps_api.Graph.nodes(data=True) if ((ps_api.Graph.out_degree(x)>0) & (y['ObjTypeName'][0] in ['Small Molecule', 'SmallMol']))]
        print('Found %d drugs in Resnet' % len(FoundDrugs))
        ReaxysAPI = RxAPI.Reaxys_API()
        ReaxysAPI.OpenSession(load_api_config())
        foundRxProps = 0
        print("Start looking for Reaxys properties")
        for drug in FoundDrugs:
            try: inchikeys = drug['InChIKey']
            except KeyError:
                try: RXNIds = drug['Reaxys ID']
                except KeyError: continue
            ReaxysProps = ReaxysAPI.GetCompoundProps(inchikeys, 'IDE.INCHI', ReaxysFields)
            if len(ReaxysProps) == 0:
                ReaxysProps = ReaxysAPI.GetCompoundProps(RXNIds, 'IDE.XRN', ReaxysFields)
            if len(ReaxysProps)> 0:
                drug.update(ReaxysProps)
                foundRxProps +=1
        ReaxysAPI.disconnect()
        print('Found Reaxys properties for %d out of %d Resnet drugs' % (foundRxProps, len(FoundDrugs)))

        fileRx = str(args.infile)[:len(str(args.infile))-4]+ '+Rx.tsv'
        EntityProps = ps_api.entProps+ReaxysFields
        ps_api.Graph.print_references(fileRx, ps_api.relProps, EntityProps)

        print("%d relations supported by %d references and annotated by Reaxys fields are in file: %s" % (ps_api.Graph.number_of_edges(),ps_api.Graph.size(weight='weight'),fileRx))
        print("Time to fetch drugs linked to %s found in %s ---" % (str(args.infile),ps_api.execution_time(start_time)))

import time
import ElsevierAPI.ResnetAPI as PsAPI
import ElsevierAPI.ResnetAPI.PathwayStudioGOQL as OQL
import ElsevierAPI.ResnetAPI.ResnetAPISession as PsAPIsn
import ReaxysAPI as RxAPI
import ElsevierAPI.ReaxysAPI.Reaxys_API as rxAPI
import argparse
import textwrap

#initiating Pathway Studio WSDL access
#RxAPI = ReaxysAPI.RxAPI
PSnx = PsAPI.networx

start_time = time.time()
#to view allowed Property and ObjectType names for use in OQL queries:
#ZeepPSAPI.DBcaller.DumpPropNames("D:\\Python\\PS_API\\DBPropertyTypes.txt")
#ZeepPSAPI.DBcaller.DumpObjNames("D:\\Python\\PS_API\\DBObjectTypes.txt")
def GOQLtoFindDrugs(TargetIds:list, TargetType = 'Protein', drugEffect=['negative']):
    if TargetType == 'Protein':
        return OQL.GetDrugs(ForTargetsIDlist=TargetIds)
    elif TargetType == 'Small Molecule':
        REL_TYPES = ['Regulation', 'MolSynthesis']
        OQLquery = OQL.ExpandEntity(PropertyValues=TargetIds,SearchByProperties=['id'],ExpandWithRelationTypes=REL_TYPES,ExpandToNeighborTypes=['Small Molecule'],direction='upstream')
        OQLquery += ' AND Effect = (' + ','.join(drugEffect)+')'
        return OQLquery
    else: 
        REL_TYPES = ['Regulation']
        OQLquery = OQL.ExpandEntity(PropertyValues=TargetIds,SearchByProperties=['id'],ExpandWithRelationTypes=REL_TYPES,ExpandToNeighborTypes=['Small Molecule'],direction='upstream')
        OQLquery += ' AND Effect = (' + ','.join(drugEffect)+')'
        return OQLquery



#fileIn = 'targets for anti-aging.txt'
#TargetType = 'Protein'
fileIn = 'aging.txt'
TargetType = ''
#fileIn = 'oxidizing agent.txt'
#TargetType = 'Small Molecule'

if __name__ == "__main__":
    EntityListFile =str()
    
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
    parser.add_argument('-p', '--resnet_props', type=str, default='Name,Alias')
    parser.add_argument('-e', '--effect', type=str, default='negative')
    parser.add_argument('-r', '--reaxys_prop', type=str)
    parser.add_argument('--debug', action="store_true")
    args = parser.parse_args()

    if type(args.infile) != type(None):
        EntityListFile = args.infile
    else: print('No entity list file was specified')

    TargetType = args.target_type 
    RetreivePsProps = args.resnet_props.split(',')
    drugEffect = [args.effect]
  
    if type(args.reaxys_prop) != type(None):
        ReaxysFields = args.reaxys_prop.split(',')
    else:
        print ('No Reaxys properties specified!')
        ReaxysFields = []


    with open(EntityListFile) as f:
        EntitiesToExpand = [line.rstrip('\n') for line in f]

    TargetIDs = PSnx.GetObjIdsByProps(PropertyValues=EntitiesToExpand,SearchByProperties=RetreivePsProps)
    OQLquery = GOQLtoFindDrugs(TargetIDs,TargetType=TargetType, drugEffect=drugEffect)
    PSdumpFile = EntityListFile[:len(EntityListFile)-4]+'_psdump.tsv'
    sn = PsAPIsn.APISession(OQLquery, PSnx)
    sn.AddDumpFile(PSdumpFile,replace_main_dump=True)
    sn.entProps = ['Name','Reaxys ID', 'Class']#Data dump columns will be ordered according to the order in this list
    sn.relProps = ['Name','Sentence','PMID','DOI','PubYear','RelationNumberOfReferences']#Data dump columns will be ordered according to the order in this list
    sn.GOQLquery = OQLquery
    sn.ProcessOQL()

    if len(ReaxysFields) > 0:
        FoundDrugs = [y for x,y in sn.Graph.nodes(data=True) if ((sn.Graph.out_degree(x)>0) & (y['ObjTypeName'][0] in ['Small Molecule', 'SmallMol']))]
        print('Found %d drugs in Resnet' % len(FoundDrugs))
        ReaxysAPI = rxAPI.Reaxys_API()
        ReaxysAPI.Load(RxAPI.APIconfig['ReaxysURL'],'',RxAPI.APIconfig['ReaxysUsername'],RxAPI.APIconfig['ReaxysPassowrd'],RxAPI.APIconfig['ReaxysAPIkey'])
        foundRxProps = 0
        print("Start looking for Reaxys properties")
        for drug in FoundDrugs:
            try: RXNIds = drug['Reaxys ID']
            except KeyError: continue
            ReaxysProps = ReaxysAPI.GetCompoundProps(RXNIds, 'IDE.XRN', ReaxysFields, RxAPI.APIconfig)
            if len(ReaxysProps)> 0:
                drug.update(ReaxysProps)
                foundRxProps +=1
        ReaxysAPI.disconnect()
        print('Found Reaxys properties for %d out of %d Resnet drugs' % (foundRxProps, len(FoundDrugs)))

        fileRx = EntityListFile[:len(EntityListFile)-4]+ '+Rx.tsv'
        EntityProps = sn.entProps+ReaxysFields
        sn.PrintReferenceView(fileRx, sn.relProps, EntityProps)

        print("%d relations supported by %d references and annotated by Reaxys fields are in file: %s" % (sn.Graph.number_of_edges(),sn.Graph.size(weight='weight'),fileRx))
        print("Time to fetch proteins linked to %s found in %s ---" % (fileIn,PsAPI.ExecutionTime(start_time)))

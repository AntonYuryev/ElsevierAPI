import time
import ElsevierAPI.ResnetAPI.ZeepToNetworkx as ZNx
import pandas as pd
import argparse
import textwrap
import ElsevierAPI
import ElsevierAPI.ResnetAPI.PathwayStudioGOQL as OQL
from ElsevierAPI import networx as PSnx


ReferenceViewNewColumn1 = PSnx.mappedBy_propName
ReferenceViewNewColumn2 = 'Concept name query'
temp_id_col = 'entity_IDs'
col_name_prefix = "RefCount to "

class SearchParam:
    MappedEntityPandas = pd.DataFrame()
    def __init__(self, EntityPandas:pd.DataFrame,prop_names_in_header,use_cache=False,need_references=False, AnnotateNodesWith=['Name']):
        self.REL_PROPS=[]
        self.ENT_PROPS=AnnotateNodesWith
        self.__cntCache='counts_cache.tsv'
        self.__refCache='reference_cache.tsv'
        self.__mapfile='map_file.tsv'
        self.need_references=need_references
        if self.need_references == True:
            self.REL_PROPS = ['Name','Sentence','PubYear','Title']

        if use_cache == False:
            self.all_entity_ids = self.MapEntities(EntityPandas,prop_names_in_header)
        else:
            try: self.ReadCache()
            except FileNotFoundError:
                self.all_entity_ids = self.MapEntities(EntityPandas,prop_names_in_header)


    def MapEntities(self,EntityPandas:pd.DataFrame, prop_names_in_header = False):
        start_mapping_time  = time.time()
        if prop_names_in_header == False:
            print('Entity infile does not have header with property names. Will use 1st column for mapping as Name and then as Alias')
            EntityPandas.rename(columns={0:'Name'}, inplace=True)
            EntityPandas['Alias'] = EntityPandas['Name']
        
        PropNameToPropToEntityID = dict()
        PropNameToPropToResnetName = dict()
        RemainToMapPandas = pd.DataFrame(EntityPandas)
        mapped_count = 0
        for propName in EntityPandas.columns:
            identifiers = list(RemainToMapPandas[propName])
            identifiers = list(map(str,identifiers))
            propValToPSobj = PSnx.MapPropToEntities(identifiers,propName,set(self.ENT_PROPS),GetChilds=True)
            PropNameToPropToEntityID[propName] = propValToPSobj
            PropNameToPropToResnetName[propName] = propValToPSobj
            mapped_count = mapped_count + len(propValToPSobj)
            RemainToMapPandas = RemainToMapPandas[~RemainToMapPandas[propName].isin(list(propValToPSobj.keys()))]

        print('%d rows out of %d were mapped to %d entities that have at least one connection in Resnet using %s' % (mapped_count,len(EntityPandas),PSnx.Graph.number_of_nodes(),','.join(EntityPandas.columns)))

        def get_entIDs(col_name,cell_value):
            propValToPSobj = PropNameToPropToEntityID[col_name]
            try:
                mappedPSobj = propValToPSobj[str(cell_value)].values()
                obj_names = [o['Name'] for o in mappedPSobj]
                name_union = set().union(*obj_names)
                resnet_names = ','.join(name_union)

                obj_ids = [o['Id'] for o in mappedPSobj]
                all_obj_ids = set().union(*obj_ids)
                try:
                    child_ids = [o[PSnx.child_ids_propName] for o in mappedPSobj]
                    all_obj_ids.update(set().union(*child_ids))
                    return resnet_names,str(cell_value),tuple(all_obj_ids)
                except KeyError:
                    return resnet_names,str(cell_value),tuple(all_obj_ids)
            except KeyError:
                return None,None,None
        
        def apply_and_concat(dataframe, field, func, column_names):
            return pd.concat((dataframe,dataframe[field].apply(lambda cell: pd.Series(func(field,cell),index=column_names))),axis=1)

        for propName in EntityPandas.columns:
            MappedEntitiesByProp = pd.DataFrame(EntityPandas)
            MappedEntitiesByProp= apply_and_concat(MappedEntitiesByProp,propName,get_entIDs,['Resnet name',PSnx.mappedBy_propName,temp_id_col])
            MappedEntitiesByProp = MappedEntitiesByProp[MappedEntitiesByProp[temp_id_col].notnull()]
            self.MappedEntityPandas = self.MappedEntityPandas.append(MappedEntitiesByProp)

        print ('Mapped %d out of %d identifiers to entities in database in %s\n' % (len(self.MappedEntityPandas), len(EntityPandas.index),ElsevierAPI.ExecutionTime(start_mapping_time)))
        entity_id_tuples = self.MappedEntityPandas[temp_id_col].tolist()
        set_list = [set(x) for x in entity_id_tuples]
        return list(set.union(*set_list))

    def LinkToConcept(self,ConceptName,ConceptIDs:list,ConnectByRelTypes=[],RelEffect=[],RelDirection=''):
        print('Linking input entities to %s concept with %d ontology childs' % (ConceptName,len(ConceptIDs)))
        if (len(ConceptIDs) > 500 and len(self.MappedEntityPandas) > 500):
            print('%s concept has %d ontology childs! Linking may take a while, be patient' % (ConceptName,len(ConceptIDs)-1))
        NewColumnName = col_name_prefix+ ConceptName
        self.MappedEntityPandas.insert(len(self.MappedEntityPandas.columns),NewColumnName,0)

        effecStr = ','.join(RelEffect)
        if len(effecStr) == 0: effecStr = 'all'

        relTypeStr = ','.join(ConnectByRelTypes)
        if len(relTypeStr) == 0: relTypeStr = 'all'

        linked_entities_counter = 0
        start_time  = time.time()
        FoundRelations = PSnx.SemanticRefCountByIds(ConceptIDs,self.all_entity_ids,ConnectByRelTypes,RelEffect,RelDirection,REL_PROPS=self.REL_PROPS, ENTITY_PROPS=self.ENT_PROPS)
        if FoundRelations.size() > 0:
            for regulatorID, targetID, rel in FoundRelations.edges.data('relation'):
                try:
                    entity_search_attr = ','.join(PSnx.Graph.nodes[regulatorID][PSnx.mappedBy_propName])
                except KeyError:
                    entity_search_attr = ','.join(PSnx.Graph.nodes[targetID][PSnx.mappedBy_propName])

                rel[ReferenceViewNewColumn1] = [entity_search_attr]
                rel[ReferenceViewNewColumn2] = [ConceptName]
      
            AccumulateReference = set()
            for idx in self.MappedEntityPandas.index:
                idx_entity_ids = list(self.MappedEntityPandas.at[idx,temp_id_col])
                reference_count = PSnx.CountReferences(idx_entity_ids,ConceptIDs,inGraph=FoundRelations)
                self.MappedEntityPandas.at[idx,NewColumnName] = len(reference_count)
                if len(reference_count) > 0:
                    AccumulateReference.update(reference_count)
                    linked_entities_counter += 1
            
            self.Print() #prints cache files for temp storage to handle network interruptions
            print("Concept \"%s\" is linked to %d rows from infile by %s relations of type \"%s\" supported by %d references with effect \"%s\" in %s" % 
                (ConceptName,linked_entities_counter,FoundRelations.number_of_edges(),relTypeStr,len(AccumulateReference),effecStr,ElsevierAPI.ExecutionTime(start_time)))

        else: print("Concept \"%s\" has no links to entities in infile" % (ConceptName))
        return linked_entities_counter
    

    def Print (self, refCountsOut='', referencesOut=''):  
        PandasToPrint = pd.DataFrame(self.MappedEntityPandas)    
  
        if len(refCountsOut)>0: 
            countFile = refCountsOut
            PandasToPrint = PandasToPrint.drop(temp_id_col,axis=1)
            PandasToPrint = PandasToPrint.loc[(PandasToPrint!=0).any(axis=1)]
        else:
            print('Data is saved data to %s for protection or re-use' % self.__cntCache)
            countFile = self.__cntCache

        PandasToPrint.to_csv(countFile, sep='\t', index=False)

        if self.need_references==True:
            if len(referencesOut) > 0:
                refFile = referencesOut
                access_mode = 'w'
            else:
                refFile = self.__refCache
                access_mode = 'a'
                print('Supporting semantic triples are saved to %s for protection or re-use' % self.__refCache)
            
            RelPropsToPrint = [ReferenceViewNewColumn1] + self.REL_PROPS + ZNx.REF_ID_TYPES + [ReferenceViewNewColumn2]
            PSnx.PrintReferenceView(refFile,relPropNames=RelPropsToPrint,entPropNames=self.ENT_PROPS,access_mode=access_mode)

    def ReadCache(self):#returns last concept linked in cache
        try:
            self.MappedEntityPandas.to_csv(self.__cntCache,'\t')
            last_col_name = self.MappedEntityPandas.columns[-1]
            return last_col_name[len(col_name_prefix):]
        except FileNotFoundError:
            print('Cache was not found! Will start processing all input files from scratch')
            return FileNotFoundError


        
if __name__ == "__main__":
    EntityListFile = str()
    LinkToConcepts = pd.DataFrame()
    instructions = '''
    infile - tab-delimted file with entity idnetifiers. Specify identifier type in each column header. 
    Identifiers will be used in the order of the columns: if entity is not found by identifier from the first column, identfier from the second column will be used.
    
    infile_has_header - specifies if infile has a header.  Defaults to no header. Use "Y" if header is present. 
    If header is not specified only identifiers in the 1st column will be used for mapping by Name+Alias to database entities

    targets - comma-separated list of concept that must be semantically linked to entities from infile
    target_type - objectTypeName for targets
    targets_file - tab-delimted file, alternative input for --targets option. Header: SearchPropertyName<>Effect<>RelationType<>Direction
        SearchPropertyName - Required. concepts that must be semantically linked to entities from infile. Header value defines what properties must be used to retreive concepts. Defaults to 'Name,Alias'
        Effect - optional.  The effect sign   ebetween entity in infile and concept in targets_file. Defaults to any.
        RelationType - optional. The type of relation between entity in infile and concept in targets_file. Defaults to any.
        Direction - optional. The direction of relation between entity in infile and concept in targets_file. Allowed values: '>', '<',''. Defaults to any.

    pathways - comma-separated list of pathway names which must be semantically linked to entities from infile.
    Semantic reference count to a pathway is calculated as sum of reference count of semantic links to every pathway component

    Output: tab-delimted file with rows from infile annotated with semantic reference count between entity (row) and semantic concept (column)  

    dump_references - create additional dump file containing all references supporitng relations found by the script
    resnet_retreive_props - list of entity properties to include into dump file
    use_cache_to_resume - must be 'Y' if you need to continue interrupted download. Will look for map_file.tsv,reference_cache.tsv,counts_cache.tsv to load pre-mapped database identifiers - saves time when infile is re-used. Defaults to start from scratch.
   
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog=textwrap.dedent(instructions))
    parser.add_argument('-i', '--infile', type=str, required=True)
    parser.add_argument('-H', '--infile_has_header', action='store_true')
    parser.add_argument('-t', '--targets', default='')
    parser.add_argument('-T', '--target_type', type=str,default='')
    parser.add_argument('-f', '--targets_file', type=str,default='')
    parser.add_argument('-d', '--dump_references', action='store_true')
    parser.add_argument('-a', '--resnet_retreive_props', type=str)
    parser.add_argument('-p', '--pathways', type=str, default='')
    parser.add_argument('-c', '--use_cache_to_resume', action='store_true')
    parser.add_argument('--debug', action="store_true")

    args = parser.parse_args()

    if len(args.infile) > 0:
        EntityListFile = args.infile #full path to tab-delimted file with drugs in the first column
    else:
        print('No input file with entities was specified!!!')

    concepts_in_file = 0
    if len(args.targets_file) > 0:
        LinkToConcepts = pd.read_csv(args.targets_file, delimiter='\t', header=0, index_col=0)
        concepts_in_file = len(LinkToConcepts)
    
    concepts_in_cmd = 0
    if len(args.targets) > 0:
        additional_concetps = args.targets.split(',')
        LinkToConcepts.set_index(additional_concetps,append=True)
        concepts_in_cmd = len(additional_concetps)
    
    if len(args.target_type) > 0: #Script will work faster if you specify what list of object types in your input Entity list
        QueryObjType = args.target_type
        QueryObjType = QueryObjType.split(',')
    else:
        QueryObjType = []
        print('No object type for targets was specified. Will seacrh for entities in %s using all object types' % EntityListFile)

 
    #If you want to use identifiers other than names enter appropriate Propeprty types into SearchByProperties list
    #consult with "Resnet Entities&Properties.txt" for the list of available identifier for nodes in the knowledge graph
    global_start_time = time.time()
    header_pos = None
    if args.infile_has_header == True: header_pos = 0

    AnnotateNodesWith = []
    if len(args.resnet_retreive_props) > 0: AnnotateNodesWith = args.resnet_retreive_props.split(',')

    if args.use_cache_to_resume == False:
        EntityPandas = pd.read_csv(EntityListFile,delimiter='\t', header=header_pos)
        search_param = SearchParam(EntityPandas,args.infile_has_header,False,args.dump_references,AnnotateNodesWith)
    else:
        search_param = SearchParam(None,None,True,args.dump_references,AnnotateNodesWith)
         
    if len(args.pathways) > 0:
        print("Begin linking entities mapped from infile to pathways")
        LinkToPathways = args.pathways
        LinkToPathways = LinkToPathways.split(",")
        start_time = time.time()
        for PathwayName in LinkToPathways:
            PathwayMembersIdToEntity = PSnx.GetPathwayMemberIds([PathwayName],SearchPathwaysBy=['Name'],FilterBy=['Protein','FunctionalClass','Complex'], InProperty='objectType')
            pathway_components = set(PathwayMembersIdToEntity.keys())
            QueryOntology = OQL.GetChildEntities(list(pathway_components), ['id'])
            pathway_components.update(PSnx.GetObjIds(QueryOntology))

            if len(pathway_components) == 0:
                print ('No entity for %s found in the database' % PathwayName)
            else:
                search_param.LinkToConcept(PathwayName, pathway_components)

        print("Entities in file %s were linked to %s pathway in %s" % (EntityListFile, PathwayName, ElsevierAPI.ExecutionTime(start_time)))
    else:
        print('No pathways were specified for semantic linking with entities from \"%s\"' % (EntityListFile))


    if len(LinkToConcepts) > 0:
        SearchConceptBy=['Name','Alias']
        if type(LinkToConcepts.index.name) != type(None):
            SearchConceptBy=LinkToConcepts.index.name.split(',')
        
        print("\nBegin linking entities mapped from infile to %d concepts from %s and %d concept from command line" % (concepts_in_file,args.targets_file,concepts_in_cmd))
        for conceptIdx in range(0, len(LinkToConcepts.index)):
            LinkToConcept = LinkToConcepts.index[conceptIdx]
            ConceptIDs = PSnx.GetObjIdsByProps(PropertyValues=[LinkToConcept],SearchByProperties=SearchConceptBy,OnlyObjectTypes=QueryObjType)
 
            if len(ConceptIDs) == 0: print ('No entity %s found in database' % LinkToConcept)
            else:
                ConnectByRelTypes = list()
                RelEffect = list()
                RelDirection = ''
                
                if 'RelationType' in LinkToConcepts.columns:
                    rel_t = LinkToConcepts.at[LinkToConcept,'RelationType']
                    if pd.notna(rel_t): ConnectByRelTypes = str(rel_t).split(',')

                if 'Effect' in LinkToConcepts.columns:
                    ef = LinkToConcepts.at[LinkToConcept,'Effect']
                    if pd.notna(ef): RelEffect = str(ef).split(',')

                if 'Direction' in LinkToConcepts.columns:
                    dr = LinkToConcepts.at[LinkToConcept,'Direction']
                    if pd.notna(dr): RelDirection = str(dr)

                linked_entities_count= search_param.LinkToConcept(LinkToConcept,list(ConceptIDs),ConnectByRelTypes,RelEffect,RelDirection)
                print("Total %d out of %d concepts interrogated in %s" % (conceptIdx+1, len(LinkToConcepts),ElsevierAPI.ExecutionTime(global_start_time)))

    print("%d entities in file %s were mapped linked to %d concepts from %s and %d concepts on cmd in %s" % (len(EntityPandas),EntityListFile,concepts_in_file,args.targets_file,concepts_in_cmd,ElsevierAPI.ExecutionTime(global_start_time)))

    countsOutFile = EntityListFile[:len(EntityListFile)-4]+'+SemanticRefcount.tsv'
    refOutFile = ''
    if search_param.need_references == True:
        refOutFile = EntityListFile[:len(EntityListFile)-4]+'+SemanticReferences.tsv'
    search_param.Print(countsOutFile,refOutFile)

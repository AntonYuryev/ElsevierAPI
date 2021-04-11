import ElsevierAPI
from ElsevierAPI.ResnetAPI import PathwayStudioZeepAPI as PSAPI
from ElsevierAPI.ResnetAPI import PathwayStudioGOQL as OQL
import json
import networkx as nx
import time
import math

class PSObject(dict): #{PropId:[values], PropName:[values]}
    def __init__(self, ZeepObjectRef):
        zeepIter = iter(ZeepObjectRef)
        while True:
            try: item = next(zeepIter)
            except StopIteration: break  #Iterator exhausted: stop the loop
            else: self[item]= [ZeepObjectRef[item]]
                
    def __hash__(self):
        return self['Id'][0]
    
    def AddSingleProperty(self, PropId, PropValue:str):
        self[PropId] = [PropValue]

    def AddProperty(self, PropId, PropValue:str):
        try: self[PropId].append(PropValue)
        except KeyError: self[PropId] = [PropValue]

    def AddUniqueProperty(self, PropId, PropValue:str):
        try:
            uniqPropList = set(self[PropId])
            uniqPropList.update([PropValue])
            self[PropId] = list(uniqPropList)
        except KeyError: self[PropId] = [PropValue]

    def AddUniquePropertyList(self, PropId, PropValues:list):
        try:
            uniqPropList = set(self[PropId])
            uniqPropList.update(PropValues)
            self[PropId] = list(uniqPropList)
        except KeyError: self[PropId] = list(set(PropValues))

    def PropValuesToStr(self, propID, cell_sep=';'):
        try: return cell_sep.join(self[propID])
        except KeyError: return ''
    
    def PropValuesToList(self, propID):
        try: return self[propID]
        except KeyError: return []

    def DataToStr (self, columnPropNames:list, col_sep='\t', cell_sep=';', endOfline='\n'):
        #assumes all properties in columnPropNames were fetched from Database otherwise will crash
        table_row = str()       
        for propName in columnPropNames:
            try:
                values = self[propName]
                propVal = cell_sep.join(values)
            except KeyError:
                propVal = ''
                
            table_row = table_row + propVal + col_sep

        return table_row[0:len(table_row)-1] + endOfline


REF_ID_TYPES = ['PMID','DOI','PII','PUI','EMBASE']
REF_PROPS = ['Sentence','PubYear','Authors','Journal','MedlineTA','CellType','CellLineName','Organ','Tissue','Organism']

class Reference(PSObject): #Identifiers{REF_ID_TYPES[i]:identifier}; self{REF_PROPS[i]:value}
    pass
    def __init__(self, idType:str, ID:str):
        self.Identifiers = {idType:ID}

    def __key(self):
        for i in range(0,len(REF_ID_TYPES)):
            try: 
                self_identifier = self.Identifiers[REF_ID_TYPES[i]]
                return self_identifier
            except KeyError: continue
        return NotImplemented

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self,other):
        if isinstance(other, Reference):
            for i in range(0,len(REF_ID_TYPES)):
                try: 
                    self_identifier = self.Identifiers[REF_ID_TYPES[i]]
                    try:
                        other_identifier = other.Identifiers[REF_ID_TYPES[i]]
                        return self_identifier == other_identifier
                    except KeyError: continue
                except KeyError: continue


    def ToString(self,IdType:str, sep='\t'):
        to_return = self.Identifiers[IdType]
        for propId, propValues in self.items():
            to_return = to_return+sep+propId+':'+';'.join(propValues)
        return to_return

    def MergeReference(self, ref_to_merge):
        self.Identifiers.update(ref_to_merge.Identifiers)



class PSRelation(PSObject):
    pass
    def __init__(self, ZeepObjectRef):
        zeepIter = iter(ZeepObjectRef)
        while True:
            try:
                item = next(zeepIter)
            except StopIteration:
                break  # Iterator exhausted: stop the loop
            else:
                val = ZeepObjectRef[item]
                self[item]= [val]
        self.PropSetToProps = dict() #{PropSetID:{PropID:[values]}}
        self.Nodes = dict() #{"Regulators':[(entityID, Dir, effect)], "Targets':[(entityID, Dir, effect)]}
        self.References = dict() #{refIdentifier (PMID, DOI, EMBASE, PUI, LUI, Title):Reference}
       
    def __hash__(self):
        return self['Id'][0]
        
    def IsDirectional(self, Links):
        propId = self['Id'][0]
        if len(Links[propId]) > 1:
            return True
        else:
            return False


    def PropValuesToStr(self, propID, cell_sep=';'):
        for id, values in self.items():
            if id == propID:
                return cell_sep.join(map(str,values))

        propSetVals = []
        for prop in self.PropSetToProps.values():            
            for id, values in prop.items():
                if id == propID:
                    propSetVals.append(cell_sep.join(map(str,values)))
                    break
             
        return cell_sep.join(propSetVals)


    def PropValuesToList(self, propID, cell_sep=';'):
        if propID in self.keys():
            return self[propID]
        else:
            to_return = []
            for prop in self.PropSetToProps.values():
                propSetVal = str()
                for id, values in prop.items():
                    if id == propID:
                        propSetVal = cell_sep.join(values)
                        break
                to_return.append(propSetVal)
            return to_return


    def load_references(self):
        for propSet in self.PropSetToProps.values():
            propset_references = list()
            for i in range(0,len(REF_ID_TYPES)):
                ref_id_type = REF_ID_TYPES[i]
                try:
                    ref_id = propSet[ref_id_type][0]
                    propset_references.append((ref_id_type,ref_id))
                    #propset_article_identifiers.append(ref_id_value)
                except KeyError: continue

            if len(propset_references) == 0:
                try: 
                    propsetTtitle = propSet['Title'][0] #trying id reference by title as a last resort since it does not have valid identifiers 
                    try:
                        Ref =  self.References[propsetTtitle]
                    except KeyError: 
                        Ref = Reference('Title', propsetTtitle)
                        self.References[propsetTtitle] = Ref
                except KeyError: continue
            else:
                #case when reference have valid identifiers
                propset_article_identifiers = {x[1] for x in propset_references}
                existing_ref = {r for i,r in self.References.items() if i in propset_article_identifiers}
                if len(existing_ref) == 0:#case when reference is new
                    id_type = propset_references[0][0]
                    id_value = propset_references[0][1]
                    Ref = Reference(id_type,id_value)
                    self.References[id_value] = Ref
                    for Id in range (1, len(propset_references)):
                        id_type = propset_references[Id][0]
                        id_value = propset_references[Id][1]
                        Ref.Identifiers[id_type] = id_value
                        self.References[id_value] = Ref
                elif len(existing_ref) > 1:#identifiers from one propset point to several references
                    #will merge all references from the propset with the first one 
                    conflictRefs = list(existing_ref)
                    anchorRef = conflictRefs[0]
                    for i in range(1,len(conflictRefs)):
                        ref_to_merge = conflictRefs[i]
                        anchorRef.MergeReference(ref_to_merge)
                        for id_value in ref_to_merge.Identifiers.values():
                            self.References[id_value] = anchorRef
                else:
                    Ref = next(iter(existing_ref))#when existing_ref == 1: there is nothing to do

            for propId, propValues in propSet.items():#adding all other valid properties to Ref
                if propId in REF_PROPS: 
                    Ref.AddUniquePropertyList(propId, propValues)


    def TripleToStr(self,columnPropNames:list,return_dict = False,col_sep='\t',cell_sep=';',endOfline='\n',RefNumPrintLimit=0):
        #assumes all properties in columnPropNames were fetched from Database otherwise will crash
        #initializing table
        colCount = len(columnPropNames)+2
        RelationNumberOfReferences = int(self['RelationNumberOfReferences'][0])

        rowCount = 1
        if RelationNumberOfReferences >= RefNumPrintLimit:
            rowCount = max(1, len(self.PropSetToProps))

        regulatorIDs = str()
        targetIDs = str()
        for k,v in self.Nodes.items():
            if k == 'Regulators':
                regIDs = [x[0] for x in v]
                #regulatorIDs = ','.join([str(int) for int in regIDs])
                regulatorIDs = ','.join(map(str,regIDs))
            else:
                trgtIDs = [x[0] for x in v]
                #targetIDs = ','.join([str(int) for int in trgtIDs])
                targetIDs = ','.join(map(str,trgtIDs))


        first_row = [''] * colCount
        first_row[colCount-2] = regulatorIDs 
        first_row[colCount-1] = targetIDs 
        table = {0: first_row}

        for r in range(1,rowCount):
            new_row = [''] * colCount
            new_row[colCount-2] = regulatorIDs 
            new_row[colCount-1] = targetIDs 
            table[r] = new_row

        
        for col in range (len(columnPropNames)):
            propId = columnPropNames[col]
            if propId in self.keys():
                for row in range(0, rowCount):
                    propValue = self.PropValuesToStr(propId)
                    table[row][col] = propValue
            elif RelationNumberOfReferences >= RefNumPrintLimit:
                row = 0
                for propList in self.PropSetToProps.values():
                    if propId in propList.keys():
                        propValues = propList[propId]
                        cellValue = cell_sep.join(propValues)
                        table[row][col] = cellValue
                        row +=1

        if return_dict == True: return table
        else:
            tableStr = str()
            for row in table.values():
                tableStr = tableStr + col_sep.join(row) + endOfline
            return tableStr

    def GetRegulatorTargetPairs(self):
        #relEntities = self.Links[relId]
        if len(self.Nodes) > 1:
            RegTargetPairs = []
            for regTriple in self.Nodes['Regulators']:
                for targetTriple in self.Nodes['Targets']:
                    pairTuple = (regTriple[0],targetTriple[0])
                    RegTargetPairs.append(pairTuple)
            return RegTargetPairs
        else:
            import itertools
            try:
                objIdList = [x[0] for x in self.Nodes['Regulators']]
                return itertools.combinations(objIdList, 2)
            except KeyError: 
                objIdList = [x[0] for x in self.Nodes['Targets']]
                return itertools.combinations(objIdList, 2)

    def GetEntitiesIDs(self):
        nodeIds = [x[0] for x in self.Nodes['Regulators']] + [x[0] for x in self.Nodes['Targets']]
        return list(set(nodeIds))

    def toJSON(self):
        str1 = '{"Relation Properties": ' + json.dumps(self) + '}'
        strP = '{"Relation References": ' + json.dumps(self.PropSetToProps) + '}'
        strR = json.dumps(self.Nodes)
        return str1 +'\n'+ strP +'\n'+strR +'\n'


class PSNetworx(PSAPI.DataModel):
    def __init__(self, DBModel:PSAPI.DataModel):
        self.SOAPclient = DBModel.SOAPclient
        self.IdtoObjectType = DBModel.IdtoObjectType #from DB.ObjectTypes
        self.IdToPropType = DBModel.IdToPropType #from DB.PropTypes 
        self.PropIdToDict = DBModel.PropIdToDict #for dictionary properties from DB.Dicts
        self.IDtoRelation = dict() #{relID:PSRelation} needs to be - Resnet relations may not be binary
        self.Graph = nx.MultiDiGraph()
        self.IdToFolders = DBModel.IdToFolders
        self.mappedBy_propName = 'mapped_by'
        self.child_ids_propName = 'Child Ids'

    def __ZeepToPSObjects(self, ZeepObjects):
        IDtoEntity = dict()
        if type(ZeepObjects) == type(None): return IDtoEntity
        for o in ZeepObjects.Objects.ObjectRef:
            psObj = PSObject(o)
            objId = o.Id
            IDtoEntity[objId] = psObj  

        for prop in ZeepObjects.Properties.ObjectProperty:
            objId = prop.ObjId
            propId = prop.PropId
            propName = prop.PropName
            vals = prop.PropValues.string
            IDtoEntity[objId][propId] = vals
            IDtoEntity[objId][propName] = vals

            try:
                propDisplayName = prop.PropDisplayName
                IDtoEntity[objId][propDisplayName] = vals
            except AttributeError: continue
                   
        return IDtoEntity


    def _load_graph(self, ZeepRelations, ZeepObjects):
        newGraph = nx.MultiDiGraph()

        #loading entities and their properties
        IDtoEntity = self.__ZeepToPSObjects(ZeepObjects)
        newGraph.add_nodes_from([(k,v.items()) for k,v in IDtoEntity.items()])
       
        newRelations = dict()
        for rel in ZeepRelations.Objects.ObjectRef:
            psRel = PSRelation(rel)
            relId = rel['Id']
            newRelations[relId] = psRel
        
        #loading relations and their properties
        for prop in ZeepRelations.Properties.ObjectProperty:
            relId = prop['ObjId']
            propId = prop['PropId']
            PropSetId = prop['PropSet']
            propName = prop['PropName']
            propDisplayName = prop['PropDisplayName']
            vals = prop['PropValues']['string']

            if self.IdToPropType[propId]['IsMultiple'] == False:
                newRelations[relId][propId] = vals
                newRelations[relId][propName] = vals
                newRelations[relId][propDisplayName] = vals
            elif PropSetId in newRelations[relId].PropSetToProps.keys():
                newRelations[relId].PropSetToProps[PropSetId][propId] = vals
                newRelations[relId].PropSetToProps[PropSetId][propName]= vals
                newRelations[relId].PropSetToProps[PropSetId][propDisplayName] = vals
            else:                  
                newRelations[relId].PropSetToProps[PropSetId] = {propId:vals}
                newRelations[relId].PropSetToProps[PropSetId] = {propName:vals}
                newRelations[relId].PropSetToProps[PropSetId] = {propDisplayName:vals}
        
        #loading connected entities from Links
        for l in ZeepRelations.Links.Link:
            relId = l['RelationId']
            Dir = l['Dir']
            link = (l['EntityId'], Dir, l['Effect'])

            if Dir == 1:
                if len(newRelations[relId].Nodes) < 2:
                    newRelations[relId].Nodes['Targets'] = [link]
                else:
                    newRelations[relId].Nodes['Targets'].append(link)
            else:
                if len(newRelations[relId].Nodes) < 1:
                    newRelations[relId].Nodes['Regulators'] = [link]
                else:
                    newRelations[relId].Nodes['Regulators'].append(link)

        
        for rel in newRelations.values():
            RegTargPairs = rel.GetRegulatorTargetPairs()
            for pair in RegTargPairs:
                refCount = rel['RelationNumberOfReferences'][0]
                
                newGraph.add_edge(pair[0], pair[1], relation=rel, weight=float(refCount))
                #print (newGraph.get_edge_data(pair[0], pair[1]))
                   
        self.IDtoRelation.update(newRelations)#must be kept since Resnet relation may not be binary
        self.Graph = nx.compose(newGraph,self.Graph)
        
        return newGraph
    
    def LoadGraphFromOQL(self,OQLquery:str,REL_PROPS=[],ENTITY_PROPS=[],getLinks=True):
        relPropsSet = set(REL_PROPS)
        relPropsSet.update(['Name','RelationNumberOfReferences'])
        relPropsSet = list(relPropsSet)

        entPropSet = set(ENTITY_PROPS)
        entPropSet.update(['Name'])
        entPropSet = list(entPropSet)

        if getLinks == True:
            ZeepRelations = self.GetData(OQLquery,relPropsSet,getLinks=True)
            if type(ZeepRelations) != type(None):
                objIdlist = list(set([x['EntityId'] for x in ZeepRelations.Links.Link]))
                ZeepObjects = self.GetObjProperties(objIdlist, entPropSet)
                return self._load_graph(ZeepRelations, ZeepObjects)
        else:
            return self.GetData(OQLquery,entPropSet,getLinks=False)


    def __objId_byOQL(self, OQLquery:str):
        ZeepEntities = self.GetData(OQLquery,RetreiveProperties=['Name'],getLinks=False)
        if type(ZeepEntities) != type(None):
            objIDs = set([x['Id'] for x in ZeepEntities.Objects.ObjectRef])
            return objIDs
        else: return set()

    def GetObjIdsByProps(self, PropertyValues:list,SearchByProperties=['Name','Alias'], GetChilds=True, OnlyObjectTypes=[]):
        QueryNode = OQL.GetEntitiesByProps(PropertyValues,SearchByProperties,OnlyObjectTypes)
        TargetIDs = self.__objId_byOQL(QueryNode)
        if GetChilds == True:
            QueryOntology = OQL.GetChildEntities(PropertyValues,SearchByProperties,OnlyObjectTypes)
            ChildIDs = self.__objId_byOQL(QueryOntology)
            TargetIDs.update(ChildIDs)
        return TargetIDs


    def MapPropToEntities(self,PropertyValues:list,propName:str,RetreiveNodeProperties:set,OnlyObjectTypes=[],GetChilds=False,MinConnectivity=1):
        entProps = set(RetreiveNodeProperties)
        entProps.update([propName,'Name'])
        entProps = list(entProps)

        step = 1000
        iterCount = math.ceil(len(PropertyValues)/step)
        print('Will use %d %s identifiers to find Resnet entities in %d iterations' % (len(PropertyValues),propName,iterCount))

        IDtoEntity = dict()
        for i in range (0,len(PropertyValues),step):
            start_time = time.time()
            propval_chunk = PropertyValues[i:i+step]
            QueryNode = OQL.GetEntitiesByProps(propval_chunk,[propName],OnlyObjectTypes,MinConnectivity)
            ZeepEntities = self.LoadGraphFromOQL(QueryNode,REL_PROPS=[],ENTITY_PROPS=entProps,getLinks=False)
            if type(ZeepEntities) != type(None):
                IDtoEntitychunk = self.__ZeepToPSObjects(ZeepEntities)
                IDtoEntity.update(IDtoEntitychunk)
                print("%d from %d iteration: %d entities were found for %d %s identifiers in %s" % (i/step+1, iterCount, len(IDtoEntitychunk), len(propval_chunk), propName, ElsevierAPI.ExecutionTime(start_time) ))
        
        LazyChildDict = dict()
        ChildIDtoPSobj = dict()
        has_childs = 0
        PropToPSobj = dict()
        for psobj in IDtoEntity.values():
            psobj_id = psobj['Id']
            prop_vals = [x for x in psobj[propName] if x in PropertyValues]
            mappedBy_propvalue = propName+':'+','.join(prop_vals)
            psobj.AddUniqueProperty(self.mappedBy_propName,mappedBy_propvalue)
            
            if GetChilds == True:
                lazy_key = tuple(psobj_id)
                try: 
                    ChildIDs = LazyChildDict[lazy_key]
                except KeyError:
                    QueryOntology = OQL.GetChildEntities(psobj_id,['Id'],OnlyObjectTypes)
                    ZeepEntities = self.LoadGraphFromOQL(QueryOntology,REL_PROPS=[],ENTITY_PROPS=entProps,getLinks=False)
                    if type(ZeepEntities) != type(None):
                        has_childs += 1
                        childIDstoEntities = self.__ZeepToPSObjects(ZeepEntities)
                        for child in childIDstoEntities.values():
                            child.AddUniqueProperty(self.mappedBy_propName,mappedBy_propvalue)
                        ChildIDtoPSobj.update(childIDstoEntities)
                        ChildIDs = list(childIDstoEntities.keys())
                    else: ChildIDs = []

                    LazyChildDict[lazy_key] = ChildIDs
                    psobj.AddUniquePropertyList(self.child_ids_propName,ChildIDs)
                
            for prop_val in prop_vals:
                try:
                    PropToPSobj[prop_val][psobj['Id'][0]] = psobj
                except KeyError: 
                    PropToPSobj[prop_val] = {psobj['Id'][0]:psobj}

        if GetChilds == True:
            print ('Childs for %d entities were found in database' % has_childs)
            IDtoEntity.update(ChildIDtoPSobj)

        self.Graph.add_nodes_from([(k,v.items()) for k,v in IDtoEntity.items()])
        print('%d out of %d %s identifiers were mapped on entities in the database' % (len(PropToPSobj),len(PropertyValues),propName))
        return PropToPSobj

    def GetNeighbors(self, EntityIDs:set, Graph=None):
        if type(Graph) == type(None):
            Graph = self.Graph
 
        IdToNeighbors = dict()
        for Id in EntityIDs:
            Neighbors =  set([x for x in nx.all_neighbors(Graph, Id)])                
            IdToNeighbors[Id] = Neighbors
        return IdToNeighbors

    def GetSubGraph(self,between_nodeids:list,and_nodeids:list,inGraph:nx.MultiDiGraph=None):
        G = inGraph if type(inGraph) != type(None) else self.Graph
        sub_graph = nx.MultiDiGraph()
        for n1 in between_nodeids:
            for n2 in and_nodeids:
                if G.has_edge(n1,n2) == True:
                    for i in range(0,len(G[n1][n2])):
                        sub_graph.add_edge(n1, n2, relation=G[n1][n2][i]['relation'], weight=G[n1][n2][i]['weight'])
                if G.has_edge(n2,n1) == True:
                    for i in range(0,len(G[n2][n1])):
                        sub_graph.add_edge(n2, n1, relation=G[n2][n1][i]['relation'], weight=G[n2][n1][i]['weight'])

        return sub_graph

    def __count_references(self, inGraph:nx.MultiDiGraph=None):
        G = inGraph if type(inGraph) != type(None) else self.Graph
        ReferenceSet = set()
        for regulatorID, targetID, rel in G.edges.data('relation'):
            rel.load_references()
            ReferenceSet.update(rel.References.values())
        return ReferenceSet

    def CountReferences(self,between_nodeids:list,and_nodeids:list,inGraph:nx.MultiDiGraph=None):
        G = inGraph if type(inGraph) != type(None) else self.Graph
        sub_graph = self.GetSubGraph(between_nodeids,and_nodeids,G)
        return self.__count_references(sub_graph)

        
    def SemanticRefCountByIds(self,node1ids:list,node2ids:list,ConnectByRelTypes=[],RelEffect=[],RelDirection='',REL_PROPS=[],ENTITY_PROPS=[]):
        relProps = set(REL_PROPS)
        relProps.update(REF_ID_TYPES)#required fields for loading references to count
        relProps.update(['Name','RelationNumberOfReferences','Source']) #just in case

        entProps = set(ENTITY_PROPS)
        entProps.update(['Name'])

        step_size = 1000
        AccumulateRelation = nx.MultiDiGraph()
        for n1 in range (0, len(node1ids), step_size):
            n1end = min(n1+step_size,len(node1ids))
            n1ids = node1ids[n1:n1end]
            for n2 in range (0, len(node2ids), step_size):
                n2end = min(n2+step_size,len(node2ids))
                n2ids = node2ids[n2:n2end]
                OQLConnectquery = OQL.ConnectEntitiesIds(n1ids,n2ids,ConnectByRelTypes,RelEffect,RelDirection)   
                FoundRelations = self.LoadGraphFromOQL(OQLConnectquery,REL_PROPS=list(relProps),ENTITY_PROPS=list(entProps))
                if type(FoundRelations) != type(None): 
                    AccumulateRelation = nx.compose(FoundRelations,AccumulateRelation)

        return AccumulateRelation
    
    def SemanticRefCountByIdsOLD(self,node1ids:list,node2ids:list,ConnectByRelTypes=[],RelEffect=[],RelDirection='',REL_PROPS=[],ENTITY_PROPS=[]):
        relProps = set(REL_PROPS)
        relProps.update(['Name','RelationNumberOfReferences','Source'])
        relProps.update(REF_ID_TYPES)

        entProps = set(ENTITY_PROPS)
        entProps.update(['Name'])

        step_size = 1000
        AccumulateRelation = nx.MultiDiGraph()
        AccumulateReference = set()
        for n1 in range (0, len(node1ids), step_size):
            n1end = min(n1+step_size,len(node1ids))
            n1ids = node1ids[n1:n1end]
            for n2 in range (0, len(node2ids), step_size):
                n2end = min(n2+step_size,len(node2ids))
                n2ids = node2ids[n2:n2end]
                OQLConnectquery = OQL.ConnectEntitiesIds(n1ids,n2ids,ConnectByRelTypes,RelEffect,RelDirection)   
                FoundRelations = self.LoadGraphFromOQL(OQLConnectquery,REL_PROPS=list(relProps),ENTITY_PROPS=list(entProps))
                if type(FoundRelations) != type(None):
                    newRefs = self.__count_references(FoundRelations)   
                    AccumulateReference.update(newRefs)
                    AccumulateRelation = nx.compose(FoundRelations,AccumulateRelation)

        return AccumulateReference, AccumulateRelation

    def SemanticRefCount(self, node1PropValues:list, node1PropTypes:list, node1ObjTypes:list, node2PropValues:list, node2PropTypes:list, node2ObjTypes:list):
        node1ids = self.GetObjIdsByProps(node1PropValues,node1PropTypes,node1ObjTypes)
        AccumulateReference = set()
        AccumulateRelation = nx.MultiDiGraph()
        if len(node1ids) > 0:
            node2ids = self.GetObjIdsByProps(node2PropValues,node2PropTypes,node2ObjTypes)
            if len(node2ids) > 0:
                AccumulateReference, AccumulateRelation = self.SemanticRefCountByIds(node1ids, node2ids)

        return AccumulateReference, AccumulateRelation


    def SetEdgeProperty(self, nodeId1, nodeId2, PropertyName, PropertyValues:list, bothDirs=True):
        if self.Graph.has_edge(nodeId1,nodeId2) == True:
            for i in range(0,len(self.Graph[nodeId1][nodeId2])):
                self.Graph[nodeId1][nodeId2][i]['relation'][PropertyName]=PropertyValues
        if bothDirs == True:
            if self.Graph.has_edge(nodeId2,nodeId1) == True:
                for i in range(0,len(self.Graph[nodeId2][nodeId1])):
                    self.Graph[nodeId2][nodeId1][i]['relation'][PropertyName]=PropertyValues


    def GetGraphEntityIds(self, OnlyRelationTypes:list, OnlyInRelations=[]):
        if len(OnlyInRelations) == 0:
            AllIDs = [x for x,y in self.Graph.nodes(data=True) if y['ObjTypeName'][0] in OnlyRelationTypes]
            return AllIDs
        else:
            LinkedEntitiesIds = set([rel.GetEntitiesIDs() for rel in OnlyInRelations])
            AllowedIDs = [x for x,y in self.Graph.nodes(data=True) if y['ObjTypeName'][0] in OnlyRelationTypes]
            FilteredRelationNeighbors = LinkedEntitiesIds.intersection(AllowedIDs)
            return list(FilteredRelationNeighbors)
        
        
    def GetProperties(self, IDList:set, PropertyName):
            IdToProps = {x:y[PropertyName] for x,y in self.Graph.nodes(data=True) if x in IDList}
            return IdToProps

    def PrintTriples(self, fileOut, relPropNames, relList=[], access_mode='w', printHeader=True):
        with open(fileOut, access_mode, encoding='utf-8') as f:
            if printHeader == True:
                header = '\t'.join(relPropNames)+'\t'+"Regulators Id"+'\t'+"Targets Id"
                f.write(header + '\n')
            
            relToPrint = relList
            if len(relList) == 0: relToPrint = self.IDtoRelation.values() #dump everything
                
            for rel in relToPrint:
                f.write(rel.TripleToStr(relPropNames))


    def __get_node(self, nodeId, G:nx.MultiDiGraph=None):
        if type(G) == type(None): G = self.Graph
        dic = {k:v for k,v in G.nodes[nodeId].items()}
        PSobj = PSObject([])
        PSobj.update(dic)
        return PSobj

    def PrintReferenceView(self, fileOut, relPropNames, entPropNames=[], G:nx.MultiDiGraph=None, access_mode='w', printHeader=True, RefNumPrintLimit=0):
        rel_props = relPropNames
        rel_props.append('Name') if 'Name' not in rel_props else rel_props

        if type(G) == type(None): G = self.Graph
        with open(fileOut, access_mode, encoding='utf-8') as f:
            if printHeader == True:
                header = '\t'.join(rel_props)+'\t'+"Regulators Id"+'\t'+"Targets Id"
                targetPropheader = [''] * len(entPropNames)
                regPropheader = [''] * len(entPropNames) 

                for i in range(0,len(entPropNames)):
                    regPropheader[i] = 'Regulator:'+entPropNames[i]
                    targetPropheader[i] = 'Target:'+entPropNames[i]

                if len(regPropheader) > 0:
                    header = '\t'.join(regPropheader)+'\t'+header+'\t'+'\t'.join(targetPropheader)
                    
                f.write(header + '\n')

            if len(entPropNames) == 0:
                for regulatorID, targetID, rel in G.edges.data('relation'):
                    ReferenceViewTriple = rel.TripleToStr(rel_props)
                    f.write(ReferenceViewTriple)
            else:
                for regulatorID, targetID, rel in G.edges.data('relation'):
                    reg = self.__get_node(regulatorID, G)
                    targ = self.__get_node(targetID, G)
                    regPropsStr = reg.DataToStr(entPropNames)
                    targPropsStr = targ.DataToStr(entPropNames)
                    relPropsStrList = rel.TripleToStr(rel_props, return_dict=True, RefNumPrintLimit=RefNumPrintLimit)

                    ReferenceTableView = str()
                    for row in relPropsStrList.values():
                        ReferenceTableView = ReferenceTableView+regPropsStr[0:len(regPropsStr)-1]+'\t'+'\t'.join(row)+'\t'+targPropsStr
                    f.write(ReferenceTableView)

    
    def DumpEntities (self, fileOut, PropNames, EntityIDlist=[]):
        header = '\t'.join(PropNames)
        NodeDict = dict([(i,v) for i,v in self.Graph.nodes(data=True)])
        with open(fileOut, 'w', encoding='utf-8') as f:
            f.write(header + '\n')
            if len(EntityIDlist) == 0:
                for ent in NodeDict.values():
                    f.write(ent.DataToStr(PropNames))
            else:
                for ent in NodeDict.values():
                    if ent['Id'] in EntityIDlist:
                        f.write(ent.DataToStr(PropNames))
        
    def FindDrugs(self, ForTargetsIDlist:list,REL_PROPS:list, ENTITY_PROPS:list):
        OQLquery = OQL.GetDrugs(ForTargetsIDlist)
        ZeepRelations = self.GetData(OQLquery, REL_PROPS)
        if type(ZeepRelations) != type(None):
            objIdlist = list(set([x['EntityId'] for x in ZeepRelations.Links.Link]))
            ZeepObjects = self.GetObjProperties(objIdlist, ENTITY_PROPS)
            newPSRelations = self._load_graph(ZeepRelations, ZeepObjects)
            return newPSRelations

    def FindReaxysSubstances(self, ForTargetsIDlist:list,REL_PROPS:list, ENTITY_PROPS:list):
        OQLquery = OQL.GetReaxysSubstances(ForTargetsIDlist)
        ZeepRelations = self.GetData(OQLquery, REL_PROPS)
        if type(ZeepRelations) != type(None):
            objIdlist = list(set([x['EntityId'] for x in ZeepRelations.Links.Link]))
            ZeepObjects = self.GetObjProperties(objIdlist, ENTITY_PROPS)
            newPSRelations = self._load_graph(ZeepRelations, ZeepObjects)
            return newPSRelations

    def ConnectEntities(self,PropertyValues1:list,SearchByProperties1:list,EntityTypes1:list,PropertyValues2:list,SearchByProperties2:list,EntityTypes2:list, REL_PROPS:list, ConnectByRelationTypes=[],ENTITY_PROPS=[]):
        OQLquery = OQL.ConnectEntities(PropertyValues1,SearchByProperties1,EntityTypes1,PropertyValues2,SearchByProperties2,EntityTypes2, ConnectByRelationTypes)
        ZeepRelations = self.GetData(OQLquery, REL_PROPS)
        if type(ZeepRelations) != type(None):
            objIdlist = list(set([x['EntityId'] for x in ZeepRelations.Links.Link]))
            ZeepObjects = self.GetObjProperties(objIdlist, ENTITY_PROPS)
            newPSRelations = self._load_graph(ZeepRelations, ZeepObjects)
            return newPSRelations

    def GetPPIs(self, InteractorIdList:set, REL_PROPS:list, ENTITY_PROPS:list):
        splitter = []
        splitter.append(list(InteractorIdList))
        import math
        number_of_splits = int(math.log2(len(InteractorIdList)))
        PPIskeeper = nx.MultiDiGraph()
        for s in range (1, number_of_splits):
            new_splitter = []
            half = int(len(splitter[0])/2)
            for split in splitter:
                uqList1 = split[0:half]
                uqList2 = split[half:]
                OQLquery = OQL.GetPPIs(set(uqList1),set(uqList2))
                ZeepRelations = self.GetData(OQLquery, REL_PROPS)
                if type(ZeepRelations) != type(None):
                    objIdlist = list(set([x['EntityId'] for x in ZeepRelations.Links.Link]))
                    ZeepObjects = self.GetObjProperties(objIdlist, ENTITY_PROPS)
                    newPSRelations = self._load_graph(ZeepRelations, ZeepObjects)
                    PPIskeeper = nx.compose(PPIskeeper,newPSRelations)
                
                new_splitter.append(uqList1)
                new_splitter.append(uqList2)

            splitter = new_splitter
            s += 1

        return PPIskeeper

    def GetObjectsFromFolders(self, FolderIds:list, PropertyNames=['Name']):
        IDtoEntity = dict()
        for id in FolderIds:
            ZeepObjects = self.GetFolderObjectsProps(id, PropertyNames)
            IDtoEntity.update(self.__ZeepToPSObjects(ZeepObjects))
        
        return IDtoEntity

    def GetPathwayMemberIds(self, PathwayIds:list, SearchPathwaysBy:list=['id'],OnlyEntities:list=[],WithProperties:list=['Name','Alias']):      
        if SearchPathwaysBy[0] in ('id', 'Id', 'ID'):
            GOQLquery = 'SELECT Entity WHERE MemberOf (SELECT Network WHERE id = ('+','.join([str(int) for int in PathwayIds])+'))'
        else:
            PropertyNames, Values = OQL.GetSearchStrings(SearchPathwaysBy,PathwayIds)
            GOQLquery = 'SELECT Entity WHERE MemberOf (SELECT Network WHERE ('+PropertyNames+') = ('+Values+'))'

        if len(OnlyEntities) > 0:
            FilterPropName, FilterValues = OQL.GetSearchStrings(WithProperties,OnlyEntities)
            GOQLquery = GOQLquery + ' AND ('+ FilterPropName +') = (' + FilterValues + ')'
        
        ZeepObjects = self.GetData(GOQLquery, RetreiveProperties=['Name'], getLinks=False)
        IDtoEntity = dict()
        IDtoEntity = self.__ZeepToPSObjects(ZeepObjects)
        return IDtoEntity

    def FindTargetsInPathways(self, DrugProps:list, DrugSearchPropertyNames:list, PathwayNames:list, RelationTypes=['DirectRegulation'], TargetTypes=['Protein']):
        REL_PROPS = ['Name','RelationNumberOfReferences','DOI','PMID','Source']
        REL_TYPES = ','.join(RelationTypes)
        TARGET_TYPES = ','.join(TargetTypes)
        DrugPropNames, DrugPropValues = OQL.GetSearchStrings(DrugSearchPropertyNames,DrugProps)
        DrugQuery = 'Select Entity WHERE ('+DrugPropNames +') = ('+ DrugPropValues +')'
        PathwayQuery =  'SELECT Network WHERE Name = ('+OQL.joinWithQuotes(',',PathwayNames)+')'
        OQLquery = 'SELECT Relation WHERE objectType = ({RelTypes}) AND NeighborOf upstream (SELECT Entity WHERE MemberOf ({pathways}) AND objectType = ({Target_Types})) AND NeighborOf downstream ({drug})'
        
        DrugsToTargetsDRGraph = self.LoadGraphFromOQL(OQLquery.format(RelTypes=REL_TYPES,pathways=PathwayQuery,Target_Types=TARGET_TYPES,drug=DrugQuery), REL_PROPS = REL_PROPS)
        OQLquery = 'SELECT Relation WHERE objectType = Binding AND NeighborOf (SELECT Entity WHERE MemberOf ({pathways}) AND objectType = ({Target_Types})) AND NeighborOf ({drug})'
        DrugsToTargetsBINDGraph = self.LoadGraphFromOQL(OQLquery.format(pathways=PathwayQuery,Target_Types=TARGET_TYPES,drug=DrugQuery),REL_PROPS = REL_PROPS)


        if type(DrugsToTargetsDRGraph) != type(None):
            if type(DrugsToTargetsBINDGraph) != type(None):
                return nx.compose(DrugsToTargetsDRGraph, DrugsToTargetsBINDGraph)
            else:
                return DrugsToTargetsDRGraph
        else:
            if type(DrugsToTargetsBINDGraph) != type(None):
                return DrugsToTargetsBINDGraph
            else:
                return None


    def FindDrugToxicities(self, DrugIds:list, DrugSearchPropertyNames:list, MinRefNumber=0, REL_PROPS=['Name','RelationNumberOfReferences'], ENTITY_PROPS=['Name']):
        DrugPropNames, DrugPropValues = OQL.GetSearchStrings(DrugSearchPropertyNames,DrugIds)
        DrugQuery = 'Select Entity WHERE ('+DrugPropNames +') = ('+ DrugPropValues +')'
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = positive AND RelationNumberOfReferences >= '+str(MinRefNumber) +' AND NeighborOf upstream (SELECT Entity WHERE objectType = Disease) AND NeighborOf downstream ({drug})'
        OQLquery = OQLquery.format(drug=DrugQuery)
        DrugsToDiseaseGraph = self.LoadGraphFromOQL(OQLquery, REL_PROPS, ENTITY_PROPS)
        ClinParamOntologyBranch = 'SELECT Entity WHERE InOntology (SELECT Annotation WHERE Ontology=\'Pathway Studio Ontology\' AND Relationship=\'is-a\') under (SELECT OntologicalNode WHERE Name=\'disease-related parameters\')'
        OQLquery = 'SELECT Relation WHERE objectType = Regulation AND Effect = positive AND RelationNumberOfReferences >= '+str(MinRefNumber) +' AND NeighborOf upstream ({ClinicalParameters}) AND NeighborOf downstream ({drug})'
        OQLquery = OQLquery.format(drug=DrugQuery, ClinicalParameters=ClinParamOntologyBranch)
        DrugsToClinParamGraph = self.LoadGraphFromOQL(OQLquery, REL_PROPS, ENTITY_PROPS)

        if type(DrugsToDiseaseGraph) != type(None):
            if type(DrugsToClinParamGraph) != type(None):
                return nx.compose(DrugsToDiseaseGraph, DrugsToClinParamGraph)
            else:
                return DrugsToDiseaseGraph
        else:
            if type(DrugsToClinParamGraph) != type(None):
                return DrugsToClinParamGraph
            else:
                return None

    def RenameRelationProperty(self, oldPropertyName='MedlineTA', newPropertyName='Journal'):
        for regulatorID, targetID, rel in self.Graph.edges.data('relation'):
            try: rel[newPropertyName] = rel.pop(oldPropertyName)
            except KeyError:
                for prop in rel.PropSetToProps.values():
                    try: prop[newPropertyName] = prop.pop(oldPropertyName)
                    except KeyError: continue
                continue

    def GetRelationsInPathways (self,PathwayNames:list, REL_PROPS=['Name','RelationNumberOfReferences'], ENTITY_PROPS=['Name']):
        OQLquery = 'SELECT Relation WHERE MemberOf (SELECT Network WHERE Name = {pathway})'
        PatwhayRelations = nx.multidigraph()
        for pthwy in PathwayNames:
            OQLquery = OQLquery.format(pathway = pthwy)
            pathwayGraph = self.LoadGraphFromOQL(OQLquery, REL_PROPS, ENTITY_PROPS)
            PatwhayRelations = nx.compose(PatwhayRelations,pathwayGraph)

        return PatwhayRelations

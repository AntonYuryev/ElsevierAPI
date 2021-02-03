from ElsevierAPI.ResnetAPI import PathwayStudioZeepAPI as PSAPI
from ElsevierAPI.ResnetAPI import PathwayStudioGOQL as OQL
import json
import itertools
import sys
import networkx as nx

class PSObject(dict): # {PropId:[values], PropName:[values]}
    def __init__(self, ZeepObjectRef):
        zeepIter = iter(ZeepObjectRef)
        while True:
            try:
                item = next(zeepIter)
            except StopIteration:
                break  #Iterator exhausted: stop the loop
            else:
                self[item]= [ZeepObjectRef[item]]

                
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


class Reference(PSObject):
    pass
    def __init__(self, idType:str, ID:str):
        self.Identifiers = {idType:ID}

    def __hash__(self):
        try: return self.Identifiers['PMID']
        except KeyError: 
            try: return self.Identifiers['DOI']
            except KeyError: 
                try: return self.Identifiers['PUI']
                except KeyError: 
                    try: return self.Identifiers['EMBASE']
                    except KeyError: return self.Identifiers['Title']

    def ToString(self,IdType:str, sep='\t'):
        to_return = self.Identifiers[IdType]
        for propId, propValues in self.items():
            to_return = to_return +sep+propId+':'+';'.join(propValues)
        return to_return

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
       
    def __hash__(self):
        return self['Id'][0]
        
    def IsDirectional(self, Links):
        propId = self['Id'][0]
        if len(Links[propId]) > 1:
            return True
        else:
            return False

    def AddSetProperty(self, propSet:int, propID:int, propValue):
        if propSet not in self.PropSetToProps.keys():
            prop = {propID:[propValue]}
            self.PropSetToProps[propSet] = prop
        elif propID not in self.PropSetToProps[propSet].keys():
            self.PropSetToProps[propSet][propID] = [propValue]
        else:
            self.PropSetToProps[propSet][propID].append(propValue)


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

    def ReferencesCounter(self):
        articleIDs = set()
        for prop in self.PropSetToProps.values():
            doi = ''
            pmid =''
            source = ''
            if 'DOI' in prop.keys():
                doi = prop['DOI'][0]
            if  'PMID' in prop.keys():
                pmid = prop['PMID'][0]
            if 'Source' in prop.keys():
                source = prop['Source'][0]

            ref_tuple = doi+'='+pmid+'='+source
            articleIDs.update([ref_tuple])

        return articleIDs


    def TripleToStr(self, columnPropNames:list, return_dict = False, col_sep='\t', cell_sep=';', endOfline='\n', RefNumPrintLimit=0):
        #assumes all properties in columnPropNames were fetched from Database otherwise will crash
        #initializing table
        colCount = len(columnPropNames)+2
        RelationNumberOfReferences = int(self['RelationNumberOfReferences'][0])
        if RelationNumberOfReferences >= RefNumPrintLimit:
            rowCount = max(1, len(self.PropSetToProps))
        else:
            rowCount = 1

        regulatorIDs = str()
        targetIDs = str()
        for k,v in self.Nodes.items():
            if k == 'Regulators':
                regIDs = [x[0] for x in v]
                regulatorIDs = ','.join([str(int) for int in regIDs])
            else:
                trgtIDs = [x[0] for x in v]
                targetIDs = ','.join([str(int) for int in trgtIDs])


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

        if return_dict == True:
            return table
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
            objIdList = [x[0] for x in self.Nodes['Regulators']]
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

    def ZeepToPSObjects(self, ZeepObjects):
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


    def LoadGraph(self, ZeepRelations, ZeepObjects):
        newGraph = nx.MultiDiGraph()
        #loading entities and their properties
          
        IDtoEntity = self.ZeepToPSObjects(ZeepObjects)
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
        self.Graph = nx.compose(self.Graph, newGraph)
           
        return newGraph
    
    def LoadGraphFromOQL(self, OQLquery:str, REL_PROPS=[], ENTITY_PROPS=[]):
        relPropsSet = set(REL_PROPS)
        relPropsSet.update(['Name','RelationNumberOfReferences'])
        entPropSet = set(ENTITY_PROPS)
        entPropSet.update(['Name'])
        ZeepRelations = self.GetData(OQLquery, relPropsSet)
        if type(ZeepRelations) != type(None):
            objIdlist = list(set([x['EntityId'] for x in ZeepRelations.Links.Link]))
            ZeepObjects = self.GetObjProperties(objIdlist, entPropSet)
            return self.LoadGraph(ZeepRelations, ZeepObjects)

    def __GetObjIdsByOQL(self, OQLquery:str):
        ZeepEntities = self.GetData(OQLquery,RetreiveProperties=['Name'],getLinks=False)
        if type(ZeepEntities) != type(None):
            objIDs = set([x['Id'] for x in ZeepEntities.Objects.ObjectRef])
            return objIDs
        else: return set()

    def GetObjIdsByProps(self, PropertyValues:list,SearchByProperties=['Name','Alias'], GetChilds=True, OnlyObjectTypes=[]):
        QueryNode = OQL.GetEntitiesByProps(PropertyValues,SearchByProperties,OnlyObjectTypes)
        TargetIDs = self.__GetObjIdsByOQL(QueryNode)
        if GetChilds == True:
            QueryOntology = OQL.GetChildEntities(PropertyValues,SearchByProperties,OnlyObjectTypes)
            ChildIDs = self.__GetObjIdsByOQL(QueryOntology)
            TargetIDs.update(ChildIDs)
        return TargetIDs
        
    def __GetPSObjectsByOQL(self, OQLquery:str, RetreiveProperties:list=['Name']):
        if 'Name' not in RetreiveProperties: RetreiveProperties.append('Name')
        ZeepEntities = self.GetData(OQLquery,RetreiveProperties=RetreiveProperties, getLinks=False)
        if type(ZeepEntities) != type(None):
            IDtoEntity = self.ZeepToPSObjects(ZeepEntities)
            return list(IDtoEntity.values())
        else: return list()
    
    def GetPSObjectsByProps(self, PropertyValues:list,SearchByProperties=['Name','Alias'],OnlyObjectTypes=[],RetreiveProperties=['Name'],GetChilds=True):
        TargetIDs = self.GetObjIdsByProps(PropertyValues,SearchByProperties,GetChilds,OnlyObjectTypes)
        OQLquery = OQL.GetEntitiesByProps(TargetIDs, ['Id'])

        if 'Name' not in RetreiveProperties: RetreiveProperties.append('Name')
        ZeepEntities = self.GetData(OQLquery,RetreiveProperties, getLinks=False)
        if type(ZeepEntities) != type(None):
            IDtoEntity = self.ZeepToPSObjects(ZeepEntities)
            return list(IDtoEntity.values())
        else: return []

    def CountReferences(self, nodeId1, nodeId2, bothDirs=True):
        articleIDcounter = set()
        if self.Graph.has_edge(nodeId1,nodeId2) == True:
            relations = [self.Graph[nodeId1][nodeId2][x]['relation'] for x in self.Graph[nodeId1][nodeId2]]
            for rel in relations:
                articleIDcounter.update(rel.ReferencesCounter())
        if bothDirs == True:
            if self.Graph.has_edge(nodeId2,nodeId1) == True:
                relations = [self.Graph[nodeId2][nodeId1][x]['relation'] for x in self.Graph[nodeId2][nodeId1]]
            for rel in relations:
                 articleIDcounter.update(rel.ReferencesCounter())

        return articleIDcounter

    def SemanticRefCountByIds(self, node1ids:list, node2ids:list):
        articleIDscounter =set()
        OQLConnectquery = OQL.ConnectEntitiesIds(node1ids, node2ids)   
        FoundRelations = self.LoadGraphFromOQL(OQLConnectquery,REL_PROPS=['Name','RelationNumberOfReferences','DOI','PMID','Source'], ENTITY_PROPS=['Name'])
        if type(FoundRelations) != type(None):
            for regulatorID, targetID, rel in FoundRelations.edges.data('relation'):
                articleIDscounter.update(rel.ReferencesCounter())
        return articleIDscounter

    def SemanticRefCount(self, node1PropValues:list, node1PropTypes:list, node1ObjTypes:list, node2PropValues:list, node2PropTypes:list, node2ObjTypes:list):
        node1ids = self.GetObjIdsByProps(node1PropValues,node1PropTypes,node1ObjTypes)
        articleIDscounter = set()
        if len(node1ids) > 0:
            node2ids = self.GetObjIdsByProps(node2PropValues,node2PropTypes,node2ObjTypes)
            if len(node2ids) > 0:
                articleIDscounter = self.SemanticRefCountByIds(node1ids, node2ids)

        return articleIDscounter


    def SetEdgeProperty(self, nodeId1, nodeId2, PropertyName, PropertyValues:list, bothDirs=True):
        if self.Graph.has_edge(nodeId1,nodeId2) == True:
            for i in range(0,len(self.Graph[nodeId1][nodeId2])):
                self.Graph[nodeId1][nodeId2][i]['relation'][PropertyName]=PropertyValues
        if bothDirs == True:
            if self.Graph.has_edge(nodeId2,nodeId1) == True:
                for i in range(0,len(self.Graph[nodeId2][nodeId1])):
                    self.Graph[nodeId2][nodeId1][i]['relation'][PropertyName]=PropertyValues


    def GetNode(self, nodeId, G:nx.MultiDiGraph=None):
        if type(G) == type(None):
            dic = {k:v for k,v in self.Graph.nodes[nodeId].items()}
        else:
            dic = {k:v for k,v in G.nodes[nodeId].items()}
        PSobj = PSObject([])
        PSobj.update(dic)
        return PSobj

    def GetGraphEntityIds(self, ObjTypeNames:list, Relations=[]):
        if len(Relations) == 0:
            AllIDs = [x for x,y in self.Graph.nodes(data=True) if y['ObjTypeName'][0] in ObjTypeNames]
            return AllIDs
        else:
            LinkedEntitiesIds = set([rel.GetEntitiesIDs() for rel in Relations])
            AllowedIDs = [x for x,y in self.Graph.nodes(data=True) if y['ObjTypeName'][0] in ObjTypeNames]
            FilteredRelationNeighbors = LinkedEntitiesIds.intersection(AllowedIDs)
            return list(FilteredRelationNeighbors)
        
    def NodePropertyValues(self, nodeId:int, PropertyName:list):
        return self.Graph.nodes[nodeId][PropertyName]
        
    def GetProperties(self, IDList:set, PropertyName):
            IdToProps = {x:y[PropertyName] for x,y in self.Graph.nodes(data=True) if x in IDList}
            return IdToProps
        

    def PrintTriples(self, fileOut, PropNames, relList=[], access_mode='w', printHeader=True):
        with open(fileOut, access_mode, encoding='utf-8') as f:
            if printHeader == True:
                header = '\t'.join(PropNames)
                header = header +'\t' + "Regulators Id"
                header = header +'\t' + "Targets Id"
                f.write(header + '\n')
            if len(relList) == 0:#dump everything
                for rel in self.IDtoRelation.values():
                    f.write(rel.TripleToStr(PropNames))
            else:
                for rel in relList:#dump by relList
                    f.write(rel.TripleToStr(PropNames))

    def PrintReferenceView(self, fileOut, relPropNames, entPropNames, G:nx.MultiDiGraph=None, access_mode='w', printHeader=True, RefNumPrintLimit=0):
        if type(G) == type(None): G = self.Graph
        with open(fileOut, access_mode, encoding='utf-8') as f:
            if printHeader == True:
                header = '\t'.join(relPropNames)
                header = header +'\t' + "Regulators Id"
                header = header +'\t' + "Targets Id"
                targetPropheader = [''] * len(entPropNames)
                regPropheader = [''] * len(entPropNames) 

                for i in range(0,len(entPropNames)):
                    regPropheader[i] = 'Regulator:'+entPropNames[i]
                    targetPropheader[i] = 'Target:'+entPropNames[i]

                header = '\t'.join(regPropheader) + '\t' + header 
                header = header +'\t' + '\t'.join(targetPropheader)
                f.write(header + '\n')

            if len(entPropNames) == 0:
                for regulatorID, targetID, rel in G.edges.data('relation'):
                    ReferenceViewTriple = rel.TripleToStr(relPropNames)
                    f.write(ReferenceViewTriple)
            else:
                for regulatorID, targetID, rel in G.edges.data('relation'):
                    reg = self.GetNode(regulatorID, G)
                    targ = self.GetNode(targetID, G)
                    regProps = reg.DataToStr(entPropNames)
                    targProps = targ.DataToStr(entPropNames)
                    relProps = rel.TripleToStr(relPropNames, return_dict=True, RefNumPrintLimit=RefNumPrintLimit)

                    ReferenceTableView = str()
                    for row in relProps.values():
                        rowStr = '\t'.join(row)
                        ReferenceTableView = ReferenceTableView+regProps[0:len(regProps)-1]+'\t'+rowStr+'\t'+targProps
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
            newPSRelations = self.LoadGraph(ZeepRelations, ZeepObjects)
            return newPSRelations

    def FindReaxysSubstances(self, ForTargetsIDlist:list,REL_PROPS:list, ENTITY_PROPS:list):
        OQLquery = OQL.GetReaxysSubstances(ForTargetsIDlist)
        ZeepRelations = self.GetData(OQLquery, REL_PROPS)
        if type(ZeepRelations) != type(None):
            objIdlist = list(set([x['EntityId'] for x in ZeepRelations.Links.Link]))
            ZeepObjects = self.GetObjProperties(objIdlist, ENTITY_PROPS)
            newPSRelations = self.LoadGraph(ZeepRelations, ZeepObjects)
            return newPSRelations

    def ConnectEntities(self,PropertyValues1:list,SearchByProperties1:list,EntityTypes1:list,PropertyValues2:list,SearchByProperties2:list,EntityTypes2:list, REL_PROPS:list, ConnectByRelationTypes=[],ENTITY_PROPS=[]):
        OQLquery = OQL.ConnectEntities(PropertyValues1,SearchByProperties1,EntityTypes1,PropertyValues2,SearchByProperties2,EntityTypes2, ConnectByRelationTypes)
        ZeepRelations = self.GetData(OQLquery, REL_PROPS)
        if type(ZeepRelations) != type(None):
            objIdlist = list(set([x['EntityId'] for x in ZeepRelations.Links.Link]))
            ZeepObjects = self.GetObjProperties(objIdlist, ENTITY_PROPS)
            newPSRelations = self.LoadGraph(ZeepRelations, ZeepObjects)
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
                    newPSRelations = self.LoadGraph(ZeepRelations, ZeepObjects)
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
            IDtoEntity.update(self.ZeepToPSObjects(ZeepObjects))
        
        return IDtoEntity

    def GetPathwayMemberIds(self, PathwayIds:list, SearchPathwaysBy:list=['id'],OnlyEntities:list=[],WithProperties:list=['Name','Alias']):      
        if SearchPathwaysBy[0] in ('id', 'Id', 'ID'):
            GOQLquery = 'SELECT Entity WHERE MemberOf (SELECT Network WHERE id = ('+','.join([str(int) for int in PathwayIds])+'))'
        else:
            PropertyNames, Values = OQL.GetSearchStrings(SearchPathwaysBy,PathwayIds)
            GOQLquery = 'SELECT Entity WHERE MemberOf (SELECT Network WHERE ('+PropertyNames+') = ('+Values+'))'

        if len(OnlyEntities) > 0:
            FilterPropName, FilterValues = OQL.GetSearchStrings(WithProperties,OnlyEntities)
            GOQLquery = GOQLquery + 'AND ('+ FilterPropName +') = (' + FilterValues + ')'
        
        ZeepObjects = self.GetData(GOQLquery, RetreiveProperties=['Name'], getLinks=False)
        IDtoEntity = dict()
        IDtoEntity = self.ZeepToPSObjects(ZeepObjects)
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


    def GetNeighbors(self, EntityIDs:set, Graph=None):
        if type(Graph) == type(None):
            Graph = self.Graph
 
        IdToNeighbors = dict()
        for Id in EntityIDs:
            Neighbors =  set([x for x in nx.all_neighbors(Graph, Id)])                
            IdToNeighbors[Id] = Neighbors
        return IdToNeighbors

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

from ZeepResnetAPI import PathwayStudioZeepAPI as PSAPI
from ZeepResnetAPI import PathwayStudioGOQL as GOQL
import json
import itertools

class PSObject(dict): # {PropId:[values], PropName:[values]}
    def __init__(self, ZeepObjectRef):
        zeepIter = iter(ZeepObjectRef)
        while True:
            try:
                item = next(zeepIter)
            except StopIteration:
                break  #Iterator exhausted: stop the loop
            else:
                val = ZeepObjectRef[item]
                self[item]= [val]
    

    def AddProperty(self, PropId, PropValue:str):
        if PropId not in self.keys():
            self[PropId] = [PropValue]
        else:
            self[PropId].append(PropValue)
    

    def PropValuesToStr(self, propID, cell_sep=';'):
        for id, values in self.items():
            if id == propID:
                return cell_sep.join(values)
        return ''
    
    def PropValuesToList(self, propID):
        for id, values in self.items():
            if id == propID:
                return values      
        return []

    def DataToStr (self, columnPropNames:list, col_sep='\t', cell_sep=';', endOfline='\n'):
        #assumes all properties in columnPropNames were fetched from Database otherwise will crash
        table_row = str()       
        for propName in columnPropNames:
            propVal = ''
            if propName in self.keys():
                values = self[propName]
                propVal = cell_sep.join(values)
            table_row = table_row + propVal + col_sep

        return table_row + endOfline


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
                return cell_sep.join(values)

        propSetVals = []
        for prop in self.PropSetToProps.values():            
            for id, values in prop.items():
                if id == propID:
                    propSetVals.append(cell_sep.join(values))
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


    def DataToStr (self, columnPropNames:list, col_sep='\t', cell_sep=';', endOfline='\n'):
        #assumes all properties in columnPropNames were fetched from Database otherwise will crash
        tableStr = str()
        #initializing table
        colCount = len(columnPropNames)
        first_row = [''] * colCount
        table = {0: first_row}
        for r in range(1,len(self.PropSetToProps)):
            new_row = [''] * colCount
            table[r] = new_row

        if len(self.PropSetToProps) == 0:
            assert len(table) == 1
            for col in range (colCount):
                propId = columnPropNames[col]
                if propId in self.keys():
                    propValue = self.PropValuesToStr(propId)
                    table[0][col] = propValue
        else:
            for col in range (colCount):
                propId = columnPropNames[col]
                if propId in self.keys():
                    #setting all rows in table to non-PropertySet value
                    propValue = self.PropValuesToStr(propId)
                    for row in range (len(self.PropSetToProps)):
                        table[row][col] = propValue
                else:
                    row = 0
                    for propList in self.PropSetToProps.values():
                        if propId in propList.keys():
                            propValues = propList[propId]
                            cellValue = cell_sep.join(propValues)
                            tablerow = table[row]
                            tablerow[col] = cellValue
                            table[row] = tablerow
                        row = row+1

        for row in table.values():
            tableStr = tableStr + col_sep.join(row) + endOfline

        return tableStr


    def TripleToStr(self, columnPropNames:list, relLinks:dict, col_sep='\t', cell_sep=';', endOfline='\n'):
        #assumes all properties in columnPropNames were fetched from Database otherwise will crash
        tableStr = str()
        #initializing table
        colCount = len(columnPropNames)+2
        rowCount = max(1, len(self.PropSetToProps))
        regulatorIDs = str()
        targetIDs = str()
        for k,v in relLinks.items():
            if k == 'Regulators':
                regIDs = [x[0] for x in v]
                regulatorIDs = ','.join([str(int) for int in regIDs])
            else:
                trgtIDs = [x[0] for x in v]
                targetIDs = ','.join([str(int) for int in trgtIDs])

        first_row = [''] * colCount
        table = {0: first_row}
        for r in range(1,len(self.PropSetToProps)):
            new_row = [''] * colCount
            table[r] = new_row

        for col in range (len(columnPropNames)):
            propId = columnPropNames[col]
            for row in range(0, rowCount):
                table[row][colCount-2] = regulatorIDs
                table[row][colCount-1] = targetIDs
                if propId in self.keys():
                    propValue = self.PropValuesToStr(propId)
                    table[row][col] = propValue
                else:
                    for propList in self.PropSetToProps.values():
                        if propId in propList.keys():
                            propValues = propList[propId]
                            cellValue = cell_sep.join(propValues)
                            tablerow = table[row]
                            tablerow[col] = cellValue
                            table[row] = tablerow

        for row in table.values():
            tableStr = tableStr + col_sep.join(row) + endOfline

        return tableStr

    def TripleToStrN(self, columnPropNames:list, col_sep='\t', cell_sep=';', endOfline='\n'):
        #assumes all properties in columnPropNames were fetched from Database otherwise will crash
        tableStr = str()
        #initializing table
        colCount = len(columnPropNames)+2
        rowCount = max(1, len(self.PropSetToProps))
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
        table = {0: first_row}
        for r in range(1,len(self.PropSetToProps)):
            new_row = [''] * colCount
            table[r] = new_row

        for col in range (len(columnPropNames)):
            propId = columnPropNames[col]
            for row in range(0, rowCount):
                table[row][colCount-2] = regulatorIDs
                table[row][colCount-1] = targetIDs
                if propId in self.keys():
                    propValue = self.PropValuesToStr(propId)
                    table[row][col] = propValue
                else:
                    for propList in self.PropSetToProps.values():
                        if propId in propList.keys():
                            propValues = propList[propId]
                            cellValue = cell_sep.join(propValues)
                            tablerow = table[row]
                            tablerow[col] = cellValue
                            table[row] = tablerow

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


class Network(PSAPI.DataModel):
    def __init__(self, DBModel:PSAPI.DataModel):
        self.SOAPclient = DBModel.SOAPclient
        self.IdtoObjectType = DBModel.IdtoObjectType
        self.IdToPropType = DBModel.IdToPropType
        self.PropIdToDict = DBModel.PropIdToDict #for dictionary properties
        self.IDtoEntity = dict() #{objID:PSObject}
        self.Links = dict() #stores {relID:{'Targets':[(entityID, Dir, effect)]}, relID:{'Regulators':[(entityID, Dir, effect)]}}
        self.IDtoRelation = dict() #{relID:PSRelation}
        
        
    def Load (self, ZeepRelations, ZeepObjects):
        #loading entities and their properties
        for o in ZeepObjects.Objects.ObjectRef:
            psObj = PSObject(o)
            self.IDtoEntity[o['Id']] = psObj

        for prop in ZeepObjects.Properties.ObjectProperty:
            objId = prop['ObjId']
            propId = prop['PropId']
            propName = prop['PropName']
            propDisplayName = prop['PropDisplayName']
            vals = prop['PropValues']['string']
            self.IDtoEntity[objId][propId] = vals
            self.IDtoEntity[objId][propName] = vals
            self.IDtoEntity[objId][propDisplayName] = vals

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
            elif len(self.IDtoRelation) > 0 and PropSetId in self.IDtoRelation[relId].PropSetToProps.keys():
                newRelations[relId].PropSetToProps[PropSetId][propId] = vals
                newRelations[relId].PropSetToProps[PropSetId][propName]= vals
                newRelations[relId].PropSetToProps[PropSetId][propDisplayName] = vals
            else:                  
                newRelations[relId].PropSetToProps[PropSetId] = {propId:vals}
                newRelations[relId].PropSetToProps[PropSetId] = {propName:vals}
                newRelations[relId].PropSetToProps[PropSetId] = {propDisplayName:vals}
        
        #loading Links
        for l in ZeepRelations.Links.Link:
            relId = l['RelationId']
            Dir = l['Dir']
            link = (l['EntityId'], Dir, l['Effect'])
            if relId not in self.Links.keys():
                if Dir == 1:
                    self.Links[relId]= {'Targets': [link]}
                else:
                    self.Links[relId]= {'Regulators' : [link]}
            else:
                if Dir == 1:
                    if 'Targets' not in self.Links[relId].keys():
                        self.Links[relId]['Targets'] = [link]
                    else:
                        self.Links[relId]['Targets'].append(link)
                else:
                    if 'Regulators' not in self.Links[relId].keys():
                        self.Links[relId]['Regulators'] = [link]
                    else:
                        self.Links[relId]['Regulators'].append(link)

        self.IDtoRelation.update(newRelations)
        return list(newRelations.values())
    
    def LoadFromOQL(self, OQLquery:str, REL_PROPS:list, ENTITY_PROPS:list):
        ZeepRelations = self.GetRelations(OQLquery, REL_PROPS)
        if type(ZeepRelations) != type(None):
            objIdlist = list(set([x['EntityId'] for x in ZeepRelations.Links.Link]))
            ZeepObjects = self.GetObjProperties(objIdlist, ENTITY_PROPS)
            return self.Load(ZeepRelations, ZeepObjects)

    def LoadRel (self, ZeepRelations, ZeepObjects):
        #loading entities and their properties
        for o in ZeepObjects.Objects.ObjectRef:
            psObj = PSObject(o)
            self.IDtoEntity[o['Id']] = psObj

        for prop in ZeepObjects.Properties.ObjectProperty:
            objId = prop['ObjId']
            propId = prop['PropId']
            propName = prop['PropName']
            propDisplayName = prop['PropDisplayName']
            vals = prop['PropValues']['string']
            self.IDtoEntity[objId][propId] = vals
            self.IDtoEntity[objId][propName] = vals
            self.IDtoEntity[objId][propDisplayName] = vals

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
            elif len(self.IDtoRelation) > 0 and PropSetId in newRelations[relId].PropSetToProps.keys():
                newRelations[relId].PropSetToProps[PropSetId][propId] = vals
                newRelations[relId].PropSetToProps[PropSetId][propName]= vals
                newRelations[relId].PropSetToProps[PropSetId][propDisplayName] = vals
            else:                  
                newRelations[relId].PropSetToProps[PropSetId] = {propId:vals}
                newRelations[relId].PropSetToProps[PropSetId] = {propName:vals}
                newRelations[relId].PropSetToProps[PropSetId] = {propDisplayName:vals}
        
        #loading regulators and targets
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

        self.IDtoRelation.update(newRelations)
        return list(newRelations.values())

    def LoadRelFromOQL(self, OQLquery:str, REL_PROPS:list, ENTITY_PROPS:list):
        ZeepRelations = self.GetRelations(OQLquery, REL_PROPS)
        if type(ZeepRelations) != type(None):
            objIdlist = list(set([x['EntityId'] for x in ZeepRelations.Links.Link]))
            ZeepObjects = self.GetObjProperties(objIdlist, ENTITY_PROPS)
            return self.LoadRel(ZeepRelations, ZeepObjects)

    def GetEntityIdList(self, ObjTypeNames:list):
        ObjList = {k:v for k,v in self.IDtoEntity.items() if v['ObjTypeName'][0] in ObjTypeNames}
        return list(ObjList.keys())

    def GetLinkedEntityIds(self, ObjTypeNames:list, Relations:list):
        LinkedEntitiesIds = []
        for rel in Relations:
            LinkedEntitiesIds += rel.GetEntitiesIDs()

        LinkedEntitiesIds = set(LinkedEntitiesIds)

        AllowedObjList = {k:v for k,v in self.IDtoEntity.items() if v['ObjTypeName'][0] in ObjTypeNames}
        AllowedIDs = set(AllowedObjList.keys())
        AllowedNeighbors = LinkedEntitiesIds.intersection(AllowedIDs)
        return list(AllowedNeighbors)
        

    def DumpRelations (self, fileOut, PropNames, relIDlist=[]):
        header = '\t'.join(PropNames)
        with open(fileOut, 'w', encoding='utf-8') as f:
            f.write(header + '\n')
            if len(relIDlist) == 0:
                for rel in self.IDtoRelation.values():
                    f.write(rel.DataToStr(PropNames))
            else:
                for rel in self.IDtoRelation.values():
                    if rel['Id'] in relIDlist:
                        f.write(rel.DataToStr(PropNames))

    def PrintRelations (self, PropNames, relIDlist=[]):
        header = '\t'.join(PropNames)
        print(header+'\n')
        if len(relIDlist) == 0:
            for rel in self.IDtoRelation.values():
                print(rel.DataToStr(PropNames))
        else:
            for rel in self.IDtoRelation.values():
                if rel['Id'] in relIDlist:
                    print(rel.DataToStr(PropNames))
    
    def DumpTriples(self, fileOut, PropNames, relIdList=[], mode='w'):
        header = '\t'.join(PropNames)
        header = header +'\t' + "Regulators Id"
        header = header +'\t' + "Targets Id"
        with open(fileOut, mode, encoding='utf-8') as f:
            f.write(header + '\n')
            if len(relIdList) == 0:#dump everything
                for relId, rel in self.IDtoRelation.items():
                    relLinks = self.Links[relId]
                    f.write(rel.TripleToStr(PropNames, relLinks))
            else:
                for relId in relIdList:
                    rel = self.IDtoRelation[relId['Id'][0]]
                    relLinks = self.Links[relId['Id'][0]]
                    f.write(rel.TripleToStr(PropNames, relLinks))

    def DumpTriplesN(self, fileOut, PropNames, relList=[], access_mode='w',printHeader=True):
        with open(fileOut, access_mode, encoding='utf-8') as f:
            if printHeader == True:
                header = '\t'.join(PropNames)
                header = header +'\t' + "Regulators Id"
                header = header +'\t' + "Targets Id"
                f.write(header + '\n')
            if len(relList) == 0:#dump everything
                for rel in self.IDtoRelation.values():
                    f.write(rel.TripleToStrN(PropNames))
            else:
                for rel in relList:#dump by relList
                    f.write(rel.TripleToStrN(PropNames))
    

    def DumpEntities (self, fileOut, PropNames, EntityIDlist=[]):
        header = '\t'.join(PropNames)
        with open(fileOut, 'w', encoding='utf-8') as f:
            f.write(header + '\n')
            if len(EntityIDlist) == 0:
                for ent in self.IDtoEntity.values():
                    f.write(ent.DataToStr(PropNames))
            else:
                for ent in self.IDtoEntity.values():
                    if ent['Id'] in EntityIDlist:
                        f.write(ent.DataToStr(PropNames))


    def GetRegulatorTargetPairs(self, relId):
        relEntities = self.Links[relId]
        if len(relEntities) > 1:
            RegTargetPairs = []
            for regTriple in relEntities['Regulators']:
                for targetTriple in relEntities['Targets']:
                    pairTuple = (regTriple[0],targetTriple[0])
                    RegTargetPairs.append(pairTuple)
            return RegTargetPairs
        else:
            objIdList = [x[0] for x in relEntities['Regulators']]
            return itertools.combinations(objIdList, 2)
        

    def FindDrugs(self, ForTargetsIDlist:list,REL_PROPS:list, ENTITY_PROPS:list):
        OQLquery = GOQL.GetDrugs(ForTargetsIDlist)
        ZeepRelations = self.GetRelations(OQLquery, REL_PROPS)
        if type(ZeepRelations) != type(None):
            objIdlist = list(set([x['EntityId'] for x in ZeepRelations.Links.Link]))
            ZeepObjects = self.GetObjProperties(objIdlist, ENTITY_PROPS)
            newPSRelations = self.LoadRel(ZeepRelations, ZeepObjects)
            return newPSRelations

    def FindReaxysSubstances(self, ForTargetsIDlist:list,REL_PROPS:list, ENTITY_PROPS:list):
        OQLquery = GOQL.GetReaxysSubstances(ForTargetsIDlist)
        ZeepRelations = self.GetRelations(OQLquery, REL_PROPS)
        if type(ZeepRelations) != type(None):
            objIdlist = list(set([x['EntityId'] for x in ZeepRelations.Links.Link]))
            ZeepObjects = self.GetObjProperties(objIdlist, ENTITY_PROPS)
            newPSRelations = self.Load(ZeepRelations, ZeepObjects)
            return newPSRelations

    def ConnectEntities(self,PropertyValues1:list,SearchByProperties1:list,EntityTypes1:list,PropertyValues2:list,SearchByProperties2:list,EntityTypes2:list, ConnectByRelationTypes=[],REL_PROPS=[], ENTITY_PROPS=[]):
        OQLquery = GOQL.ConnectEntities(PropertyValues1,SearchByProperties1,EntityTypes1,PropertyValues2,SearchByProperties2,EntityTypes2,ConnectByRelationTypes)
        ZeepRelations = self.GetRelations(OQLquery, REL_PROPS)
        if type(ZeepRelations) != type(None):
            objIdlist = list(set([x['EntityId'] for x in ZeepRelations.Links.Link]))
            ZeepObjects = self.GetObjProperties(objIdlist, ENTITY_PROPS)
            newPSRelations = self.LoadRel(ZeepRelations, ZeepObjects)
            return newPSRelations

    def GetPPIs(self, InteractorIdList:set, REL_PROPS:list, ENTITY_PROPS:list):
        splitter = []
        splitter.append(list(InteractorIdList))
        import math
        number_of_splits = int(math.log2(len(InteractorIdList)))
        PPIskeeper = []
        for s in range (1, number_of_splits):
            new_splitter = []
            half = int(len(splitter[0])/2)
            for split in splitter:
                uqList1 = split[0:half]
                uqList2 = split[half:]
                OQLquery = GOQL.GetPPIs(set(uqList1),set(uqList2))
                ZeepRelations = self.GetRelations(OQLquery, REL_PROPS)
                if type(ZeepRelations) != type(None):
                    objIdlist = list(set([x['EntityId'] for x in ZeepRelations.Links.Link]))
                    ZeepObjects = self.GetObjProperties(objIdlist, ENTITY_PROPS)
                    newPSRelations = self.LoadRel(ZeepRelations, ZeepObjects)
                    PPIskeeper = PPIskeeper + newPSRelations
                
                new_splitter.append(uqList1)
                new_splitter.append(uqList2)

            splitter = new_splitter
            s += 1

        return PPIskeeper

    

class APISession(Network):
    ResultRef = str()
    ResultPos = int()
    ResultSize = int()
    def __init__(self, net:Network):
        self.SOAPclient = net.SOAPclient
        self.IdtoObjectType = net.IdtoObjectType
        self.IdToPropType = net.IdToPropType
        self.PropIdToDict = net.PropIdToDict #for dictionary properties
        self.IDtoEntity = net.IDtoEntity #{objID:PSObject}
        self.IDtoRelation = net.IDtoRelation #{relID:PSRelation}
    
    def InitAPISession (self, OQLrequest, PageSize=0, REL_PROPS=[], ENTITY_PROPS=[]):
        ZeepRelations,(self.ResultRef, self.ResultSize, self.ResultPos)  = self.InitSession(OQLrequest, PageSize, REL_PROPS)
        if type(ZeepRelations) != type(None):
            objIdlist = list(set([x['EntityId'] for x in ZeepRelations.Links.Link]))
            ZeepObjects = self.GetObjProperties(objIdlist, ENTITY_PROPS)
            return self.LoadRel(ZeepRelations, ZeepObjects)
    
    def GetNextPage(self, PageSize=100, REL_PROPS=[], ENTITY_PROPS=[]):
        if self.ResultPos < self.ResultSize:
            ZeepRelations, self.ResultPos  = self.GetNextSessionPage(self.ResultRef, self.ResultPos, PageSize, self.ResultSize, REL_PROPS)
            if type(ZeepRelations) != type(None):
                objIdlist = list(set([x['EntityId'] for x in ZeepRelations.Links.Link]))
                ZeepObjects = self.GetObjProperties(objIdlist, ENTITY_PROPS)
                return self.LoadRel(ZeepRelations, ZeepObjects)
        else:
            return []


import requests
import pandas as pd
from ElsevierAPI.ResnetAPI import PathwayStudioGOQL as OQL

class DataModel():

    def __init__(self, url, username, password):
        self.IdToPropType = dict() 
        self.IdtoObjectType = dict()
        self.PropIdToDict = dict()
        self.IdToFolders = dict()
        from requests.auth import HTTPBasicAuth
        from requests import Session
        from zeep.cache import SqliteCache
        from zeep.transports import Transport
        session = Session()
        session.auth = HTTPBasicAuth(username, password)
        transport = Transport(cache=SqliteCache(), session = session)
        from zeep import Client 
        self.SOAPclient = Client(url, transport=transport)
        self.LoadModel()

    def LoadModel(self):
        ObjectTypes = self.SOAPclient.service.GetObjectTypes()
        PropertyTypes = self.SOAPclient.service.GetPropertyDefinitions()
        #Folders = self.SOAPclient.service.GetFoldersTree(0)

        for i in range(0, len(ObjectTypes)):
            id = ObjectTypes[i]['Id']
            ObjTypeName =  ObjectTypes[i]['Name']
            ObjTypeDisplayname =  ObjectTypes[i]['DisplayName']
            self.IdtoObjectType[ObjTypeName] = ObjectTypes[i]
            self.IdtoObjectType[ObjTypeDisplayname] = ObjectTypes[i]
            self.IdtoObjectType[id] = ObjectTypes[i]
            for s in ObjectTypes[i]['Synonyms']:
                self.IdtoObjectType[s] = ObjectTypes[i]

        for i in range(0, len(PropertyTypes)):
            id = PropertyTypes[i]['Id']
            PropTypeName =  PropertyTypes[i]['Name']
            PropTypeDisplayName =  PropertyTypes[i]['DisplayName']
            self.IdToPropType[PropTypeName] = PropertyTypes[i]
            self.IdToPropType[PropTypeDisplayName] = PropertyTypes[i]
            self.IdToPropType[id] = PropertyTypes[i]

    def __objtypes_by_classid(self, classID):
        objtype_list = list()
        ObjectTypes = self.SOAPclient.service.GetObjectTypes()
        for i in range(0, len(ObjectTypes)):
            ObjTypeName =  ObjectTypes[i]['Name']
            ObjTyepeClassID = ObjectTypes[i]['ObjClassId']
            if ObjTyepeClassID == classID:
                objtype_list.append(ObjTypeName)

        return objtype_list

    def GetRelationTypes(self):
        return self.__objtypes_by_classid(3)

    def GetEntityTypes(self):
        return self.__objtypes_by_classid(1)

    def LoadFolderTree(self):
        Folders = self.SOAPclient.service.GetFoldersTree(0)
        for i in range(0, len(Folders)):
            Folderid = Folders[i]['Id']
            FolderName =  Folders[i]['Name']
            self.IdToFolders[Folderid] = [Folders[i]]
            #several folders may have the same name:
            if FolderName in self.IdToFolders.keys():
                self.IdToFolders[FolderName].append(Folders[i])
            else:
                self.IdToFolders[FolderName] = [Folders[i]]

            
            
    def DumpPropNames (self, fileOut):
        IDPropNamesList = [  (p['Id'], p['DisplayName'], p['Name']) for p in self.IdToPropType.values() ]
        IDPropNamesList = list(set(IDPropNamesList))
        import json
        propNamesStr = json.dumps(IDPropNamesList)
        import codecs
        import sys
        UTF8Writer = codecs.getwriter('utf8')
        sys.stdout = UTF8Writer(sys.stdout)
        file_result=open(fileOut, "w", encoding='utf-8')
        file_result.write("Format: [PropertyType Id, Name in PS UI , Search string for GOQL] = \n" + propNamesStr)
        file_result.close()

    def DumpObjNames (self, fileOut):
        IDObjTypeList = [  (p['Id'], p['DisplayName'], p['Name']) for p in self.IdtoObjectType.values() ]
        IDObjTypeList = list(set(IDObjTypeList))
        import json
        ObjTypeListStr = json.dumps(IDObjTypeList)
        import codecs
        import sys
        UTF8Writer = codecs.getwriter('utf8')
        sys.stdout = UTF8Writer(sys.stdout)
        file_result=open(fileOut, "w", encoding='utf-8')
        file_result.write("Format: [ObjectType Id, Name in PS UI , Search string for GOQL] = \n" + ObjTypeListStr)
        file_result.close()

    def MapPropertyNamestoID(self, PropertyNames:list):
        id_list = [self.IdToPropType[x]['Id'] for x in PropertyNames]
        if 'Name' not in PropertyNames:
            id_list.append(self.IdToPropType['Name']['Id'])
        return id_list

    def GetDictionary(self, idProperty, idDictFolder:int):
        if idProperty not in self.PropIdToDict.keys():
            dictFolder = self.SOAPclient.service.GetDictFolder(idDictFolder)
            IdValuesToStr = dict()
            for idval in dictFolder.Values.DictValue:
                id = idval['Id']
                val = idval['Value']
                IdValuesToStr[id] = val
            self.PropIdToDict[idProperty] = IdValuesToStr
            self.PropIdToDict[dictFolder['Name']] = IdValuesToStr 
        return self.PropIdToDict[idProperty]
        
    def GetDictValue(self, IdProperty:int, DictIdValue:int):
        if IdProperty not in self.PropIdToDict[IdProperty].keys():
            idDictionary = self.IdToPropType[IdProperty]['DictFolderID']
            assert idDictionary > 0
            self.GetDictionary(IdProperty, idDictionary)
        
        return self.PropIdToDict[IdProperty][DictIdValue]

    def GetSubfolders(self, FolderIds:list):
        if len(self.IdToFolders) == 0:
            self.LoadFolderTree()
        
        SubfoldersIds = set()
        for id in FolderIds:
            Folders = self.IdToFolders[id]
            for folder in Folders:
                if type(folder['SubFolders']) != type(None):
                    subIds = folder['SubFolders']['long']
                    SubfoldersIds.update(subIds)
        return SubfoldersIds

    def GetSubfolderTree(self, FolderId):
        AccumulateSubfolders = set([FolderId])
        Subfolders = set(self.GetSubfolders([FolderId]))
        while len(Subfolders) > 0:
            AccumulateSubfolders.update(Subfolders)
            subs = set()
            for id in Subfolders:
                subs.update(self.GetSubfolders([id]))
            Subfolders = subs

        return AccumulateSubfolders

    def GetFolderObjects(self, FolderId,result_param):
        result = self.SOAPclient.service.FolderGetObjects(FolderId, result_param)
        return result

    def OQLresponse (self, OQLquery, result_param):
        result = self.SOAPclient.service.OQLSearch(OQLquery, result_param)
        return result
    
    def ResultGetData (self, result_param):
        result = self.SOAPclient.service.ResultGetData(result_param)
        return result

    def CreateResultParam (self, PropertyNames=['Name']):
        IdPropertyList = self.MapPropertyNamestoID(PropertyNames)
        PropRefs = []
        propRef = self.SOAPclient.get_type('ns0:PropertyRef')
        for id in IdPropertyList:
            pL = propRef(Id = id)
            PropRefs.append(pL)
        result_param = self.SOAPclient.get_type('ns0:ResultParam')
        rp = result_param(
            CreateResult = False,
            ResultRef = '?',
            Objects = [],
            Objects_size = 0,
            ResultPos = 0,
            MaxPageSize= 0,
            GetObjects = False,
            GetProperties = False,
            PropertyList = {'PropertyRef': PropRefs} ,
            PropertyList_Size = len(PropertyNames),
            GetLinks = False,
            GetParents = False,
            GetMembers = False,
            GetAddlCol = False,
            RefLimit = 0,
            SortColumn = 0,
            SortDescending = False,
            SortPropId = 0,
            AddlAttrs = [],
            AddlAttrs_Size = 0,
            ApplySourceFilter = False)
        return rp


    def GetObjProperties(self, objIDList:list, PropertyNames=['Name']):
        rp = self.CreateResultParam(PropertyNames)
        rp.GetObjects = True
        rp.GetProperties = True
        OQLrequest = OQL.GetObjects(objIDList)
        ObjProps = self.OQLresponse(OQLrequest, rp)
        #setting objectType name
        for obj in ObjProps.Objects.ObjectRef:
            objTypeID = obj['ObjTypeId']
            objTypeName = self.IdtoObjectType[objTypeID]['Name']
            obj['ObjTypeName'] = objTypeName

        #setting dict property values and property names
        for prop in ObjProps.Properties.ObjectProperty:
            idProperty = prop['PropId']
            PropName = self.IdToPropType[idProperty]['Name']
            PropDisplayName = self.IdToPropType[idProperty]['DisplayName']
            prop['PropName'] = PropName
            prop['PropDisplayName'] = PropDisplayName
            DictFolderId = self.IdToPropType[idProperty]['DictFolderId']
            if DictFolderId > 0:
                dictFolder = self.GetDictionary(idProperty, DictFolderId)
                for i in range(0, len(prop['PropValues']['string'])):
                    idDictPropValue = int(prop['PropValues']['string'][i])
                    newDictValue = dictFolder[idDictPropValue]
                    prop['PropValues']['string'][i] = newDictValue


        return ObjProps

    def GetFolderObjectsProps(self, FolderId, PropertyNames=['Name']):
        rp = self.CreateResultParam(PropertyNames)
        rp.GetObjects = True
        rp.GetProperties = True
        ObjProps= self.GetFolderObjects(FolderId,rp)
        if type(ObjProps.Objects) == type(None):
            return None
        #setting objectType name
        for obj in ObjProps.Objects.ObjectRef:
            objTypeID = obj['ObjTypeId']
            objTypeName = self.IdtoObjectType[objTypeID]['Name']
            obj['ObjTypeName'] = objTypeName

        #setting dict property values and property names
        for prop in ObjProps.Properties.ObjectProperty:
            idProperty = prop['PropId']
            PropName = self.IdToPropType[idProperty]['Name']
            prop['PropName'] = PropName
            DictFolderId = self.IdToPropType[idProperty]['DictFolderId']
            if DictFolderId > 0:
                dictFolder = self.GetDictionary(idProperty, DictFolderId)
                for i in range(0, len(prop['PropValues']['string'])):
                    idDictPropValue = int(prop['PropValues']['string'][i])
                    newDictValue = dictFolder[idDictPropValue]
                    prop['PropValues']['string'][i] = newDictValue

        return ObjProps

    
    def GetData(self, OQLrequest, RetreiveProperties=['Name', 'RelationNumberOfReferences'], getLinks=True):
        rp = self.CreateResultParam(RetreiveProperties)
        rp.GetObjects = True
        rp.GetProperties = True
        rp.GetLinks = getLinks
        #setting objectType name
        ObjProps = self.OQLresponse(OQLrequest, rp)
        if type(ObjProps.Objects) == type(None):
            #print('Your SOAP response is empty! Check your OQL query and try again\n')
            return None

        for obj in ObjProps.Objects.ObjectRef:
            objTypeID = obj['ObjTypeId']
            objTypeName = self.IdtoObjectType[objTypeID]['Name']
            obj['ObjTypeName'] = objTypeName
        #setting dict property values
        for prop in ObjProps.Properties.ObjectProperty:
            idProperty = prop['PropId']
            PropName = self.IdToPropType[idProperty]['Name']
            PropDisplayName = self.IdToPropType[idProperty]['DisplayName']
            prop['PropName'] = PropName
            prop['PropDisplayName'] = PropDisplayName
            DictFolderId = self.IdToPropType[idProperty]['DictFolderId']
            if DictFolderId > 0:
                dictFolder = self.GetDictionary(idProperty, DictFolderId)
                for i in range(0, len(prop['PropValues']['string'])):
                    idDictPropValue = int(prop['PropValues']['string'][i])
                    newDictValue = dictFolder[idDictPropValue]
                    prop['PropValues']['string'][i] = newDictValue

        return ObjProps


    def InitSession(self, OQLrequest, PageSize:int, PropertyNames=['Name']):
        rp = self.CreateResultParam(PropertyNames)
        rp.GetObjects = True
        rp.GetProperties = True
        rp.GetLinks = True
        rp.CreateResult = True
        rp.MaxPageSize = PageSize
        ObjProps = self.OQLresponse(OQLrequest, rp)
        if type(ObjProps.Objects) == type(None):
            #print('Your SOAP response is empty! Check your OQL query and try again\n')
            return None, (ObjProps.ResultRef, ObjProps.ResultSize, ObjProps.ResultPos)

        for obj in ObjProps.Objects.ObjectRef:
            objTypeID = obj['ObjTypeId']
            objTypeName = self.IdtoObjectType[objTypeID]['Name']
            obj['ObjTypeName'] = objTypeName
        #setting dict property values
        for prop in ObjProps.Properties.ObjectProperty:
            idProperty = prop['PropId']
            PropName = self.IdToPropType[idProperty]['Name']
            PropDisplayName = self.IdToPropType[idProperty]['DisplayName']
            prop['PropName'] = PropName
            prop['PropDisplayName'] = PropDisplayName
            DictFolderId = self.IdToPropType[idProperty]['DictFolderId']
            if DictFolderId > 0:
                dictFolder = self.GetDictionary(idProperty, DictFolderId)
                for i in range(0, len(prop['PropValues']['string'])):
                    idDictPropValue = int(prop['PropValues']['string'][i])
                    newDictValue = dictFolder[idDictPropValue]
                    prop['PropValues']['string'][i] = newDictValue

        return ObjProps, (ObjProps.ResultRef, ObjProps.ResultSize, ObjProps.ResultPos)


    def GetNextSessionPage(self, ResultRef, ResultPos, PageSize, ResultSize, PropertyNames=['Name']):
        rp = self.CreateResultParam(PropertyNames)
        rp.GetObjects = True
        rp.GetProperties = True
        rp.GetLinks = True

        rp.ResultSize = ResultSize
        rp.ResultRef = ResultRef
        rp.MaxPageSize = PageSize
        rp.ResultPos = ResultPos
        ObjProps = self.ResultGetData(rp)

        #setting objectType name
        if type(ObjProps.Objects) == type(None):
            #print('Your SOAP response is empty! Check your OQL query and try again\n')
            return

        for obj in ObjProps.Objects.ObjectRef:
            objTypeID = obj['ObjTypeId']
            objTypeName = self.IdtoObjectType[objTypeID]['Name']
            obj['ObjTypeName'] = objTypeName
        #setting dict property values
        for prop in ObjProps.Properties.ObjectProperty:
            idProperty = prop['PropId']
            PropName = self.IdToPropType[idProperty]['Name']
            PropDisplayName = self.IdToPropType[idProperty]['DisplayName']
            prop['PropName'] = PropName
            prop['PropDisplayName'] = PropDisplayName
            DictFolderId = self.IdToPropType[idProperty]['DictFolderId']
            if DictFolderId > 0:
                dictFolder = self.GetDictionary(idProperty, DictFolderId)
                for i in range(0, len(prop['PropValues']['string'])):
                    idDictPropValue = int(prop['PropValues']['string'][i])
                    newDictValue = dictFolder[idDictPropValue]
                    prop['PropValues']['string'][i] = newDictValue


        if rp.ResultPos >= rp.ResultSize:
            self.SOAPclient.service.ResultRelease(rp.ResultRef)

        return ObjProps, ObjProps.ResultPos

    def PutExperiment(self, dataframe:pd, expName, expType, entityType, MapEntitiestBy, hasPvalue=True, descriptn=''):
        #PSexp = self.SOAPclient.service('ns0:GetExperiment')
        rowCount = dataframe.shape[0]
        colCount = dataframe.shape[1]
        minFoldChange = dataframe['log2FoldChange'].min()
        maxFoldChange = dataframe['log2FoldChange'].max()
        sampleCnt = int(colCount -1)
        if hasPvalue == True:
            sampleCnt = int(sampleCnt/2)
        header = dataframe.columns

        zExperiment = self.SOAPclient.get_type('ns0:Experiment')
        zSample = self.SOAPclient.get_type('ns0:SampleDefinition')
    
        sample = zSample(
            Name = header[1],
            Type = 1,
            Subtype = 2,
            hasPValue=hasPvalue,
            Calculated=False,
            Attributes = {'string':[header[1]]},
            Attributes_size = 1 
            )
        
        #ns0:Experiment( Owner: xsd:string, DateCreated: xsd:string, )
        psExperiment = zExperiment(
            Id = 0,
            Name = expName, 
            Description = descriptn,
            SampleCount = sampleCnt,
            RowsMappedCount = 0,
            ExperimentType = 1, #ratio or logratio
            ExperimentSubType = 1, #logarithmic
            Organism = 'Homo sapiens',
            ExperimentTypeName = expType,
            EntityTypeName = entityType,
            RowsCount = rowCount,
            GeneAttributeNames = {'string':[header[0]]}, #{string: xsd:string[]}
            GeneAttributeNames_Size = 1,
            OriginalGeneID_index = 0,
            SampleDefinitions = {'SampleDefinition':[sample]}, #{SampleDefinition: ns0:SampleDefinition[]}, 
            SampleDefinitions_Size = 1,
            SampleAttributeNames = {'string':['phenotype']},
            SampleAttributeNames_Size = 1,
            MinValueSignal = minFoldChange,
            MinValueRatio = minFoldChange,
            MaxValueSignal = maxFoldChange,
            MaxValueRatio = maxFoldChange,
            ReadOnly = True,
            MaskUsage = False
        )

        expID = self.SOAPclient.service.PutExperiment(psExperiment)
        
        rowAttrs = [None] * rowCount
        sampleVals = [None] * rowCount
        for i in range(0,rowCount):
            gene_id = dataframe.iat[i, 0]
            rowAttr = self.SOAPclient.get_type('ns0:ExperimentRowAttributes')
            rt = rowAttr(
                EntityURN = '?',
                EntityAttributes = {'string':[]},
                EntityAttributes_Size =0,
                OriginalGeneID = '?',
                GeneAttributes = {'string':[gene_id]}, #GeneAttributes: {string: xsd:string[]}
                GeneAttributes_Size = 1,
                EntityId = 0
            )
            rowAttrs[i] = rt

            sampleVal = self.SOAPclient.get_type('ns0:SampleValue')
            sv = sampleVal(
                Value = dataframe.iat[i, 1],
                PValue = dataframe.iat[i, 2]
            )
            sampleVals[i] = sv


        self.SOAPclient.service.PutExperimentRowAttributes(expID, 0, 0, {'ExperimentRowAttributes':rowAttrs})
        self.SOAPclient.service.PutExperimentValues(expID, 0, 0, 0, {'SampleValue':sampleVals}) 
        
        
        zMapParam = self.SOAPclient.get_type('ns0:ExperimentMappingToolParameters')
        zm = zMapParam(
            ExperimentID = expID,
            DbType = entityType,
            DbField = MapEntitiestBy,
            IsChip = False,
            ChipID = 0,
            ChipName = '?'
        )

        #ns0:OperationStatus(OperationId: xsd:long, State: xsd:int, Error: xsd:int, Status: xsd:string, ResultId: xsd:long)
        opStatus = self.SOAPclient.service.StartExperimentMappingTool(zm)
        result =  self.SOAPclient.service.GetMappingToolResult(opStatus.OperationId)
        self.SOAPclient.service.ResultRelease(result.ResultRef)
        return expID, result


def WriteSOAPresponse(fout, SOAPresponse):
    import codecs
    import sys
    UTF8Writer = codecs.getwriter('utf8')
    sys.stdout = UTF8Writer(sys.stdout)
    file_result=open(fout, "w", encoding='utf-8')
    file_result.write(SOAPresponse)
    file_result.close()

import time
import networkx as nx
import ElsevierAPI
import ElsevierAPI.ResnetAPI.PathwayStudioGOQL as OQL
from ElsevierAPI.ResnetAPI.PathwayStudioZeepAPI import DataModel
from ElsevierAPI.ResnetAPI.NetworkxObjects import PSObject,PSRelation,REF_ID_TYPES

RNEF_EXCLUDE_NODE_PROPS = ['Id','URN','ObjClassId','ObjTypeId','ObjTypeName','OwnerId','DateCreated','DateModified']
RNEF_EXCLUDE_REL_PROPS =  RNEF_EXCLUDE_NODE_PROPS+['RelationNumberOfReferences','# of Total References','Name']

class PSNetworx(DataModel):
    def __init__(self, DBModel:DataModel):
        self.SOAPclient = DBModel.SOAPclient
        self.IdtoObjectType = DBModel.IdtoObjectType #from DB.ObjectTypes
        self.IdToPropType = DBModel.IdToPropType #from DB.PropTypes 
        self.PropIdToDict = DBModel.PropIdToDict #for dictionary properties from DB.Dicts
        self.IDtoRelation = dict() #{relID:PSRelation} needs to be - Resnet relations may not be binary
        self.Graph = nx.MultiDiGraph()
        self.IdToFolders = DBModel.IdToFolders
        self.RNEFnameToPropType = DBModel.RNEFnameToPropType
        self.mappedBy_propName = 'mapped_by'
        self.child_ids_propName = 'Child Ids'

    @staticmethod
    def __ZeepToPSObjects(ZeepObjects):
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
       
        if type(ZeepRelations) != type(None):
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
        relPropsSet = list(set(REL_PROPS+['Name','RelationNumberOfReferences']))
        entPropSet = list(set(ENTITY_PROPS+['Name']))

        if getLinks == True:
            ZeepRelations = self.GetData(OQLquery,relPropsSet,getLinks=True)
            if type(ZeepRelations) != type(None):
                objIdlist = list(set([x['EntityId'] for x in ZeepRelations.Links.Link]))
                ZeepObjects = self.GetObjProperties(objIdlist, entPropSet)
                return self._load_graph(ZeepRelations, ZeepObjects)
        else:
            ZeepObjects = self.GetData(OQLquery,entPropSet,getLinks=False)
            return self._load_graph(None, ZeepObjects)


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

        import math
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

    def GetNeighbors(self, EntityIDs:set, InGraph=None):
        G = self.Graph if type(InGraph) == type(None) else InGraph
        IdToNeighbors = dict()
        for Id in EntityIDs:
            Neighbors =  set([x for x in nx.all_neighbors(G, Id)])                
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

    def count_references(self, inGraph:nx.MultiDiGraph=None):
        G = inGraph if type(inGraph) != type(None) else self.Graph
        ReferenceSet = set()
        for regulatorID, targetID, rel in G.edges.data('relation'):
            rel.load_references()
            ReferenceSet.update(rel.References.values())
        return ReferenceSet

    def CountReferences(self,between_nodeids:list,and_nodeids:list,inGraph:nx.MultiDiGraph=None):
        G = inGraph if type(inGraph) != type(None) else self.Graph
        sub_graph = self.GetSubGraph(between_nodeids,and_nodeids,G)
        return self.count_references(sub_graph)

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


    def GetGraphEntityIds(self, OnlyForNodeTypes:list, inGraph:nx.MultiDiGraph=None):
        G = inGraph if type(inGraph) != type(None) else self.Graph
        AllIDs = [x for x,y in G.nodes(data=True) if y['ObjTypeName'][0] in OnlyForNodeTypes]
        return AllIDs
        
    def GetProperties(self, IDList:set, PropertyName):
            IdToProps = {x:y[PropertyName] for x,y in self.Graph.nodes(data=True) if x in IDList}
            return IdToProps

    def PrintTriples(self, fileOut, relPropNames, inGraph:nx.MultiDiGraph=None, relList=[], access_mode='w', printHeader=True):
        G = inGraph if type(inGraph) != type(None) else self.Graph
        with open(fileOut, access_mode, encoding='utf-8') as f:
            if printHeader == True:
                header = '\t'.join(relPropNames)+'\t'+"Regulators Id"+'\t'+"Targets Id"
                f.write(header + '\n')
            
        for regulatorID, targetID, rel in G.edges.data('relation'):
            f.write(rel.TripleToStr(relPropNames))


    def __get_node(self, nodeId, inGraph:nx.MultiDiGraph=None):
        G = inGraph if type(inGraph) != type(None) else self.Graph
        dic = {k:v for k,v in G.nodes[nodeId].items()}
        PSobj = PSObject([])
        PSobj.update(dic)
        return PSobj


    def __get_relations(self,inGraph:nx.MultiDiGraph=None):
        G = inGraph if type(inGraph) != type(None) else self.Graph
        graph_relations = set()
        for regulatorID, targetID, rel in G.edges.data('relation'):
            graph_relations.add(rel)

        return graph_relations


    def PrintReferenceView(self, fileOut, relPropNames, entPropNames=[], inGraph:nx.MultiDiGraph=None, access_mode='w', printHeader=True, RefNumPrintLimit=0,col_sep='\t'):
        rel_props = list(set(relPropNames+['Name']))
        G = inGraph if type(inGraph) != type(None) else self.Graph
        with open(fileOut, access_mode, encoding='utf-8') as f:
            if printHeader == True:
                header = col_sep.join(rel_props)+col_sep+"Regulators Id"+col_sep+"Targets Id"
                targetPropheader = [''] * len(entPropNames)
                regPropheader = [''] * len(entPropNames) 

                for i in range(0,len(entPropNames)):
                    regPropheader[i] = 'Regulator:'+entPropNames[i]
                    targetPropheader[i] = 'Target:'+entPropNames[i]

                if len(regPropheader) > 0:
                    header = col_sep.join(regPropheader)+col_sep+header+col_sep+col_sep.join(targetPropheader)
                    
                f.write(header + '\n')

            if G.number_of_edges() > 0:
                if len(entPropNames) == 0:
                    for regulatorID, targetID, rel in G.edges.data('relation'):
                        ReferenceViewTriple = str(rel.TripleToStr(rel_props))
                        f.write(ReferenceViewTriple)
                else:
                    for regulatorID, targetID, rel in G.edges.data('relation'):
                        reg = self.__get_node(regulatorID, G)
                        targ = self.__get_node(targetID, G)
                        regPropsStr = reg.DataToStr(entPropNames,col_sep=col_sep)
                        targPropsStr = targ.DataToStr(entPropNames,col_sep=col_sep)
                        relPropsStrList = dict(rel.TripleToStr(rel_props,return_dict=True,RefNumPrintLimit=RefNumPrintLimit,col_sep=col_sep))

                        ReferenceTableView = str()
                        for row in relPropsStrList.values():
                            ReferenceTableView = ReferenceTableView+regPropsStr[0:len(regPropsStr)-1]+col_sep+col_sep.join(row)+col_sep+targPropsStr
                        f.write(ReferenceTableView)
            else:
                for node_id in G.nodes:
                    n = self.__get_node(node_id, G)
                    nodePropStr = n.DataToStr(entPropNames,col_sep=col_sep)
                    f.write(nodePropStr)

    
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

    def ConnectEntities(self,PropertyValues1:list,SearchByProperties1:list,EntityTypes1:list,PropertyValues2:list,SearchByProperties2:list,EntityTypes2:list, REL_PROPS=['Name','RelationNumberOfReferences'], ConnectByRelationTypes=[],ENTITY_PROPS=['Name']):
        rel_props = set(REL_PROPS+['Name','RelationNumberOfReferences'])
        ent_props = set(ENTITY_PROPS+['Name'])
        OQLquery = OQL.ConnectEntities(PropertyValues1,SearchByProperties1,EntityTypes1,PropertyValues2,SearchByProperties2,EntityTypes2, ConnectByRelationTypes)
        ZeepRelations = self.GetData(OQLquery, list(rel_props))
        if type(ZeepRelations) != type(None):
            objIdlist = list(set([x['EntityId'] for x in ZeepRelations.Links.Link]))
            ZeepObjects = self.GetObjProperties(objIdlist, list(ent_props))
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

    def GetObjectsFromFolders(self, FolderIds:list, PropertyNames=['Name'], with_layout=False):
        IDtoEntity = dict()
        for f in FolderIds:
            ZeepObjects = self.GetFolderObjectsProps(f, PropertyNames)
            IDtoEntity.update(self.__ZeepToPSObjects(ZeepObjects))
        
        if with_layout == True:
            for i, psObj in IDtoEntity.items():
                if psObj['ObjTypeName'][0] == 'Pathway':
                    layout = self.get_layout(i)
                    psObj['layout'] = str(layout['Attachment'].decode('utf-8'))

        return IDtoEntity

    def GetAllPathways(self, PropertyNames=['Name']):
        start_time = time.time()
        print('retreiving all pathways from the database')

        if (len(self.IdToFolders)) == 0: self.LoadFolderTree()

        IDtoPathway = dict()
        URNtoPathway = dict()
        for folderList in self.IdToFolders.values():
            for folder in folderList:
                ZeepObjects = self.GetFolderObjectsProps(folder['Id'],PropertyNames)
                psObjects = self.__ZeepToPSObjects(ZeepObjects)
                for Id, psObj in psObjects.items():
                    if psObj['ObjTypeName'][0] == 'Pathway':
                        try:
                            IDtoPathway[Id].AddUniqueProperty('Folders', folder['Name'])
                        except KeyError:
                            psObj['Folders'] = [folder['Name']]
                            IDtoPathway[Id] = psObj
                            URNtoPathway[psObj['URN'][0]] = psObj

        print('Found %d pathways in the database in %s' % (len(IDtoPathway),ElsevierAPI.ExecutionTime(start_time)))
        return IDtoPathway, URNtoPathway

    def GetPathwayMemberIds(self, PathwayIds:list, SearchPathwaysBy:list=['id'],OnlyEntities:list=[],WithProperties:list=['Name','Alias']):      
        if SearchPathwaysBy[0] in ['id','Id','ID']:
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


    def toRNEF(self, inGraph:nx.MultiDiGraph=None):
        G = inGraph if type(inGraph) != type(None) else self.Graph
        import xml.etree.ElementTree as et
        resnet = et.Element('resnet')
        xml_nodes = et.SubElement(resnet,'nodes')
        localIDcounter = 0
        for nodeId,n in G.nodes(data=True):
            localID = n['URN'][0]
            xml_node = et.SubElement(xml_nodes,'node',{'local_id':localID,'urn':n['URN'][0]})
            et.SubElement(xml_node,'attr',{'name':str('NodeType'),'value':str(n['ObjTypeName'][0])})
            for prop_name, prop_values in n.items():
                if prop_name in self.RNEFnameToPropType.keys():
                    if prop_name not in RNEF_EXCLUDE_NODE_PROPS:
                        for prop_value in prop_values:
                            et.SubElement(xml_node,'attr',{'name':str(prop_name),'value':str(prop_value)})

            localIDcounter +=1

        xml_controls = et.SubElement(resnet, 'controls')
        
        graph_relations = self.__get_relations()
        for rel in graph_relations:
            control_id = rel['URN'][0]
            xml_control = et.SubElement(xml_controls,'control', {'local_id':control_id})
            et.SubElement(xml_control,'attr',{'name':str('ControlType'),'value':str(rel['ObjTypeName'][0])})
            #adding links
            regulators = rel.Nodes['Regulators']
            try:
                targets = rel.Nodes['Targets']
                for r in regulators:
                    regulator_localID = G.nodes[r[0]]['URN'][0]
                    et.SubElement(xml_control,'link', {'type':'in','ref':regulator_localID})

                for t in targets:
                    target_localID = G.nodes[t[0]]['URN'][0]
                    et.SubElement(xml_control,'link', {'type':'out','ref':target_localID})
            except KeyError:
                    for r in regulators:
                        regulator_localID = G.nodes[r[0]]['URN'][0]
                        et.SubElement(xml_control,'link', {'type':'in-out','ref':regulator_localID})

        #adding non-reference properties
            for prop_name, prop_values in rel.items():
                if prop_name in self.RNEFnameToPropType.keys():
                    if prop_name not in RNEF_EXCLUDE_REL_PROPS:
                        for prop_value in prop_values:
                            et.SubElement(xml_control,'attr',{'name':str(prop_name),'value':str(prop_value)})

        #adding references
            references = list(set(rel.References.values()))
            for i in range(0,len(references)):
                for prop_name, prop_values in references[i].items():
                    for prop_value in prop_values:
                        if prop_name in self.RNEFnameToPropType.keys():
                            et.SubElement(xml_control,'attr',{'name':str(prop_name),'value':str(prop_value), 'index':str(i)})

        return et.tostring(resnet,encoding='utf-8').decode("utf-8")

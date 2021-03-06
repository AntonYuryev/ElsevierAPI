from ElsevierAPI.ResnetAPI import ZeepToNetworkx as zNX
from ElsevierAPI.ResnetAPI import PathwayStudioGOQL as GOQL
import networkx as nx
import time
import math

class APISession(zNX.PSNetworx):
    ResultRef = str()
    ResultPos = int()
    ResultSize = int()
    PageSize = 100
    GOQLquery = str()
    IsOn1stpage = True
    relProps = ['Name', 'RelationNumberOfReferences']
    entProps = ['Name']
    DumpFiles = ['ResnetAPIsessionDump.tsv']

    def __init__(self, GOQLquery, net:zNX.PSNetworx):
        self.SOAPclient = net.SOAPclient
        self.IdtoObjectType = net.IdtoObjectType #from DB.ObjectTypes
        self.IdToPropType = net.IdToPropType #from DB.PropTypes
        self.PropIdToDict = net.PropIdToDict #for DB,Dicts
        self.IDtoRelation = net.IDtoRelation #{relID:PSRelation}
        self.Graph = net.Graph
        self.GOQLquery = GOQLquery
        self.DumpFiles = [] #array of filenames used for dumping the graph data. First element is used for dumping data obtained ProcessOQL
        #other files are used to dump additional data obtained by overridden AddGraph
    
    def __InitAPISession (self):
        ZeepRelations,(self.ResultRef, self.ResultSize, self.ResultPos)  = self.InitSession(self.GOQLquery, self.PageSize, self.relProps)
        if type(ZeepRelations) != type(None):
            objIdlist = list(set([x['EntityId'] for x in ZeepRelations.Links.Link]))
            ZeepObjects = self.GetObjProperties(objIdlist, self.entProps)
            return self.LoadGraph(ZeepRelations, ZeepObjects)
    
    def __GetNextPage(self):
        if self.ResultPos < self.ResultSize:
            ZeepRelations, self.ResultPos  = self.GetNextSessionPage(self.ResultRef, self.ResultPos, self.PageSize, self.ResultSize, self.relProps)
            if type(ZeepRelations) != type(None):
                objIdlist = list(set([x['EntityId'] for x in ZeepRelations.Links.Link]))
                ZeepObjects = self.GetObjProperties(objIdlist, self.entProps)
                return self.LoadGraph(ZeepRelations, ZeepObjects)
        else:
            return None

    def AddGraph(self, newGraph:nx.MultiDiGraph):
        pass #placeholder to derive child class from APISession (see DiseaseNetwork(APISession))

    def FlashDumpFiles(self):
        for fname in self.DumpFiles:
            open(fname,'w').close()
            print('File "%s" was cleared before processing' % fname)

    def AddDumpFile(self, dump_fname, replace_main_dump=False):
        if replace_main_dump == True:
            self.DumpFiles = []
        self.DumpFiles.append(dump_fname)

    def ExecutionTime(self, execution_start):
        from datetime import timedelta
        return "{}".format(str(timedelta(seconds=time.time()-execution_start)))

    def TabData(self, fileOut, G:nx.MultiDiGraph, access_mode='w'):
        PrintHeader = self.IsOn1stpage
        self.PrintReferenceView( fileOut, self.relProps, self.entProps, G, access_mode, PrintHeader)
    

    def ProcessOQL(self, flash_dump=False):
        from datetime import timedelta
        global_start = time.time()
        start_time = time.time()
        if flash_dump == True:
            self.FlashDumpFiles()

        OQLGraph = nx.MultiDiGraph()
        pageGraph = self.__InitAPISession()
        iterCount = math.ceil(self.ResultSize/self.PageSize)
        refCount = 0
        print("GOQL query:\n%s\n returns %d relations.\n It will be processed in %d iterations" % (self.GOQLquery, self.ResultSize, iterCount))
        for pos in range(0, self.ResultSize, self.PageSize):    
            execution_time = self.ExecutionTime(start_time)
            iteration = math.ceil(self.ResultPos/self.PageSize)
            PagerefCount = pageGraph.size(weight='weight')
            print("Iteration %s retreived %d relations supported by %d references in %s seconds from Resnet" % (iteration,pageGraph.size(),PagerefCount,execution_time))
            refCount += PagerefCount
            self.AddGraph(pageGraph)
            if len(self.DumpFiles) > 0:
                self.TabData(self.DumpFiles[0],pageGraph,'a')
                self.IsOn1stpage = False
                
            start_time = time.time()
            OQLGraph = nx.compose(OQLGraph, pageGraph)
            pageGraph = self.__GetNextPage()

        print("%d relations supported by %d references are in file: %s" % (self.ResultSize,refCount,self.DumpFiles[0]))
        self.ResultRef = ''
        self.ResultPos = 0
        self.ResultSize = 0
        self.DumpFIles = []
        self.IsOn1stpage = True
        execution_time = self.ExecutionTime(global_start)
        print("GOQL query was executed in %s in %d iterations" % (execution_time,iterCount))
        
        return OQLGraph

    def GetPPIgraph(self, foutName):
        #Build PPI network between proteins
        NetworkProteins = self.GetGraphEntityIds(ObjTypeNames=['Protein'])
        print('Retreiving PPI network for %d proteins' % (len(NetworkProteins)))
        start_time = time.time()
        PPIrelationsGraph = self.GetPPIs(NetworkProteins,self.relProps, self.entProps)
        execution_time = self.ExecutionTime(start_time)
        self.TabData(foutName, PPIrelationsGraph)
        print("PPI network with %d relations was retreived in %s ---" % (PPIrelationsGraph.size(), execution_time))
        return PPIrelationsGraph



class DiseaseNetwork(APISession):
    pass
#APISession retreives records from Resnet by iterations equal APISession.PageSize using input GOQL query
#each iteration creates newGraph object. By default, newGraph is simply added to APISession.Graph at each iteration using nx.compose()
#you can re-write APISession.AddNewGraph as it is done here function to expand newGraph with additional GOQL queries and to add more data to APISession.Graph

    ProteinsPreviousPage = set()
    InputDiseaseNames = str()
    
#add here names for entity properties you want to fecth from knowledge graph. 
    entProps = ['Name']# Data dump columns will be ordered according to the order in this list
#add here names for relation properties you want to fecth from knowledge graph
    relProps = ['Name','Mechanism','RelationNumberOfReferences', 'Source', 'DOI']#Data dump columns will be ordered according to the order in this list

    def AddGraph(self, newGraph:nx.MultiDiGraph):
        #self.Graph = nx.compose(self.Graph, newGraph)#obsolete NewGraph is added by PSNetworkx

        #pageCount = int(self.ResultSize/self.ResultPos)
        start_time = time.time()
        ProteinsCurrentPage = set(self.GetGraphEntityIds(ObjTypeNames=['Protein']))
        NewProteins = ProteinsCurrentPage.difference(self.ProteinsPreviousPage)#to avoid duplications from proteins connected to disease with multiple relations
        self.ProteinsPreviousPage = ProteinsCurrentPage
        NewProteinsCount = len(NewProteins)
        execution_time = self.ExecutionTime(start_time)
        print('This iteration retreived %d proteins' % (NewProteinsCount))
        if NewProteinsCount == 0:
            return

    #Find GVs in proteins linked to input diseases
        start_time = time.time()
        GOQLquery = GOQL.ExpandEntity(NewProteins, ['id'],['GeneticChange'], ['GeneticVariant'])
        print('Searching for GeneticVariant linked to %d proteins' % NewProteinsCount)
        RelGVgenesGraph = self.LoadGraphFromOQL(GOQLquery, self.relProps, self.entProps)
        execution_time = self.ExecutionTime(start_time)
        print("%d genetic variants for %d proteins linked to %s were retreived in %s ---" % (RelGVgenesGraph.number_of_edges(),NewProteinsCount, self.InputDiseaseNames, execution_time))
        self.TabData(self.DumpFiles[1],RelGVgenesGraph, 'a')
        self.Graph = nx.compose(self.Graph, RelGVgenesGraph)

    #Find druggable targets
        start_time = time.time()
        print('Searching for drugs binding to %d proteins' % NewProteinsCount)
        drugTOtargetsGraph = self.FindDrugs(NewProteins, self.relProps, self.entProps)
        execution_time = self.ExecutionTime(start_time)
        self.TabData(self.DumpFiles[2], drugTOtargetsGraph, 'a')
        self.Graph = nx.compose(self.Graph, drugTOtargetsGraph)

        DrugCount = set([x for x in drugTOtargetsGraph.nodes() if drugTOtargetsGraph.out_degree(x)>0])
        FoundTargets = set([x for x in drugTOtargetsGraph.nodes() if drugTOtargetsGraph.in_degree(x)>0])
        print("%d drugs for %d proteins linked to %s were retreived in %s ---" % (len(DrugCount),len(FoundTargets), self.InputDiseaseNames, execution_time))


    #find RMC compounds for nondruggable targets
        ProteinNoDrugs = NewProteins.difference(FoundTargets)
        if len(ProteinNoDrugs) > 0:
            start_time = time.time()
            print('Searching for RMC compounds binding to %d proteins' % NewProteinsCount)
            RMCtoTargetsGraph = self.FindReaxysSubstances(ProteinNoDrugs,self.relProps, self.entProps)
            if type(RMCtoTargetsGraph) != type(None):
                execution_time = self.ExecutionTime(start_time)
                CompoundCount = set([x for x in RMCtoTargetsGraph.nodes() if RMCtoTargetsGraph.out_degree(x)>0])
                TargetCount = set([x for x in RMCtoTargetsGraph.nodes() if RMCtoTargetsGraph.in_degree(x)>0])
                print("%d lead compounds for %d undruggable proteins linked to %s were retreived in %s ---" % (len(CompoundCount), len(TargetCount), self.InputDiseaseNames, execution_time))
                self.Graph = nx.compose(self.Graph, RMCtoTargetsGraph)
                self.TabData(self.DumpFiles[2], RMCtoTargetsGraph, 'a')
            else:
                print('No RMC compounds found on this iteration')

from ElsevierAPI.ResnetAPI.PathwayStudioGOQL import expand_entity
from ElsevierAPI.ResnetAPI.ResnetAPISession import APISession
import networkx as nx
import time


class DiseaseNetwork(APISession):
    pass
    # APISession retrieves records from Resnet by iterations equal APISession.PageSize using input GOQL query each
    # iteration creates newGraph object. By default, newGraph is simply added to APISession.Graph at each iteration
    # using nx.compose() you can re-write APISession.AddNewGraph as it is done here function to expand newGraph with
    # additional GOQL queries and to add more data to APISession.Graph

    ProteinsPreviousPage = set()
    InputDiseaseNames = str()

    # add here names for entity properties you want to fecth from knowledge graph.
    entProps = ['Name']  # Data dump columns will be ordered according to the order in this list
    # add here names for relation properties you want to fecth from knowledge graph
    relProps = ['Name', 'Mechanism', 'RelationNumberOfReferences', 'Source',
                'DOI']  # Data dump columns will be ordered according to the order in this list

    def add_graph(self, newGraph: nx.MultiDiGraph):
        # self.Graph = nx.compose(self.Graph, newGraph)#obsolete NewGraph is added by PSNetworkx
        # pageCount = int(self.ResultSize/self.ResultPos)
        start_time = time.time()
        ProteinsCurrentPage = set(self.get_graph_entity_ids(['Protein']))
        NewProteins = ProteinsCurrentPage.difference(
            self.ProteinsPreviousPage)  # to avoid duplications from proteins connected to disease with multiple relations
        self.ProteinsPreviousPage = ProteinsCurrentPage
        NewProteinsCount = len(NewProteins)
        execution_time = self.execution_time(start_time)
        print('This iteration retreived %d proteins' % (NewProteinsCount))
        if NewProteinsCount == 0:
            return

        # Find GVs in proteins linked to input diseases
        start_time = time.time()
        GOQLquery = expand_entity(NewProteins, ['id'], ['GeneticChange'], ['GeneticVariant'])
        print('Searching for GeneticVariant linked to %d proteins' % NewProteinsCount)
        RelGVgenesGraph = self.load_graph_from_oql(GOQLquery, self.relProps, self.entProps)
        execution_time = self.execution_time(start_time)
        print("%d genetic variants for %d proteins linked to %s were retreived in %s ---" % (
            RelGVgenesGraph.number_of_edges(), NewProteinsCount, self.InputDiseaseNames, execution_time))
        self.to_csv(self.DumpFiles[1], RelGVgenesGraph, 'a')
        self.Graph = nx.compose(self.Graph, RelGVgenesGraph)

        # Find druggable targets
        start_time = time.time()
        print('Searching for drugs binding to %d proteins' % NewProteinsCount)
        drugTOtargetsGraph = self.find_drugs(NewProteins, self.relProps, self.entProps)
        execution_time = self.execution_time(start_time)
        self.to_csv(self.DumpFiles[2], drugTOtargetsGraph, 'a')
        self.Graph = nx.compose(self.Graph, drugTOtargetsGraph)

        DrugCount = set([x for x in drugTOtargetsGraph.nodes() if drugTOtargetsGraph.out_degree(x) > 0])
        FoundTargets = set([x for x in drugTOtargetsGraph.nodes() if drugTOtargetsGraph.in_degree(x) > 0])
        print("%d drugs for %d proteins linked to %s were retreived in %s ---" % (
            len(DrugCount), len(FoundTargets), self.InputDiseaseNames, execution_time))

        # find RMC compounds for nondruggable targets
        ProteinNoDrugs = NewProteins.difference(FoundTargets)
        if len(ProteinNoDrugs) > 0:
            start_time = time.time()
            print('Searching for RMC compounds binding to %d proteins' % NewProteinsCount)
            RMCtoTargetsGraph = self.find_reaxys_substances(ProteinNoDrugs, self.relProps, self.entProps)
            if type(RMCtoTargetsGraph) != type(None):
                execution_time = self.execution_time(start_time)
                CompoundCount = set([x for x in RMCtoTargetsGraph.nodes() if RMCtoTargetsGraph.out_degree(x) > 0])
                TargetCount = set([x for x in RMCtoTargetsGraph.nodes() if RMCtoTargetsGraph.in_degree(x) > 0])
                print("%d lead compounds for %d undruggable proteins linked to %s were retreived in %s ---" % (
                    len(CompoundCount), len(TargetCount), self.InputDiseaseNames, execution_time))
                self.Graph = nx.compose(RMCtoTargetsGraph, self.Graph)
                self.to_csv(self.DumpFiles[2], RMCtoTargetsGraph, 'a')
            else:
                print('No RMC compounds found on this iteration')

from .ResnetAPISession import APISession, time, math,REFERENCE_IDENTIFIERS,NO_REL_PROPERTIES
from .ResnetGraph import REFCOUNT, ResnetGraph, nx, EFFECT,NUMBER_OF_TARGETS
from .PathwayStudioGOQL import OQL
from ..pandas.panda_tricks import df, ExcelWriter
from .Drugs4Disease import Drugs4Targets
from .Drugs4Disease import ANTAGONIST_TARGETS_WS,AGONIST_TARGETS_WS,RANK,DRUG2TARGET_REGULATOR_SCORE,PHARMAPENDIUM_ID
from .Zeep2Experiment import Experiment, Sample
import scipy.stats as stats
import numpy as np
from concurrent.futures import ThreadPoolExecutor

PPERNET = 'protein_expression_regulatory_network'
DPERNET = 'drug-protein expression regulatory network'
ACTIVATION_IN = 'activation in '


class SNEA(APISession):
    def __init__(self, APIconfig,experiment_name:str,sample_names=[],sample_ids=[]):
        super().__init__(APIconfig)
        self.PageSize = 1000
        self.DumpFiles = []
        self.APIconfig = APIconfig
        with ThreadPoolExecutor(max_workers=4, thread_name_prefix='Initializing SNEA') as executor:
            future1 = executor.submit(Experiment.download,experiment_name,APIconfig)
            future2 = executor.submit(self.__load_network,['Protein','FunctionalClass','Complex'],PPERNET)

            self.experiment = future1.result()
            expression_network = future2.result()

        self.Graph.clear_resnetgraph() # in case graph was downloaded from database
        self.Graph = self.network4experiment(expression_network)

        my_samples = self.experiment.get_samples(sample_names,sample_ids)
        self.__my_sample_names__ = [s.name() for s in my_samples]
        self.experiment.annotate(self.Graph,self.__my_sample_names__)

        self.subnetworks = dict() # {regulator_id:[target_ids]}
        self.activity_matrix = df()
        self.drugs_df = df()
        self.subnet_pvalue_cutoff = 0.05
        self.min_subnet_size = 2
        

    def __my_samples__(self):
        return self.experiment.get_samples(self.__my_sample_names__)


    @staticmethod
    def __sample_annotation(sample:Sample):
        return ACTIVATION_IN+sample.name()


    def __load_network(self,regulator_types:list,network_name:str)->'ResnetGraph':
        cache_path = 'ElsevierAPI/ResnetAPI/__pscache__/'
        cache_file = cache_path+network_name+'.rnef'
        start = time.time()
        try:
            expression_network = ResnetGraph.fromRNEF(cache_file)
        except FileNotFoundError:
            print('%s was not found\nBegin downloading relations from database' % cache_file)
            self.relProps = [EFFECT,REFCOUNT,'URN']
            self.PageSize = 1000
            self.DumpFiles = []
            reg_types_str = ','.join(regulator_types)
            if network_name == DPERNET:
                get_regulators = OQL.select_drugs()
                self.add_ent_props(['PharmaPendium ID'])
            else:
                get_regulators = f'SELECT Entity WHERE objectType = ({reg_types_str})'

            get_targets = 'SELECT Entity WHERE objectType = (Protein,FunctionalClass,Complex)'
            oql_query = f'SELECT Relation WHERE objectType = (Expression,PromoterBinding) \
    AND NeighborOf upstream ({get_targets}) AND NeighborOf downstream ({get_regulators})'
            #oql_query = 'SELECT Relation WHERE objectType = Metabolization' # for testing
            expression_network = self.process_oql(oql_query, f'Fetching {network_name}')
            expression_network = expression_network.make_simple()
            expression_network.rnef2file(cache_path+network_name+'.rnef',
            ent_prop2print=['Name','PharmaPendium ID'],rel_prop2print=['URN',EFFECT,REFCOUNT])
            print('%s with %d edges and %d nodes was downloaded from database in %s' 
            % (network_name,expression_network.number_of_edges(),expression_network.number_of_nodes(),self.execution_time(start)))
            
        expression_network.name = network_name
        return expression_network


    def network4experiment(self, regulatory_network:ResnetGraph):
        input_net_name = regulatory_network.name
        exp_name = self.experiment.name()
        urns = self.experiment.urns()
        print('"%s" experiment has %d mapped entities' % (exp_name,len(urns)))

        net_name = input_net_name+f' for "{exp_name}"'
        experiment_network = regulatory_network.regulatory_network(urns,net_name)

        print('%s experiment has %d nodes connected with %d edges' 
            %(net_name,experiment_network.number_of_edges(),experiment_network.number_of_nodes()))

        return experiment_network


    def regulators4sample(self,sample:Sample,regulomes:dict,from_graph=ResnetGraph(),target_annotation_prefix=''):
        '''
        Input
        -----
        regulomes = {regulator_id : [targets]} made by ResnetGraph.regulome_dict()
        '''
        my_graph = from_graph if from_graph else self.Graph

        sample_start = time.time()
        abs_sample_distribution = sample.data['value'].abs()
        sample_regulator_ids = set()
        annotation = target_annotation_prefix+self.experiment.name4annotation(sample)
        for regulator_id, targets in regulomes.items():
            subnetwork_values = [t[annotation][0][0] for t in targets]
            sub_pd = df.from_dict({'value':subnetwork_values})
            abs_sub_pd = sub_pd['value'].abs()
            mw_results = stats.mannwhitneyu(x=abs_sample_distribution, y=abs_sub_pd, alternative = 'less') 
            """
            'two-sided': the distributions are not equal, i.e. F(u) ≠ G(u) for at least one u.
            'less': the distribution underlying x is stochastically less than the distribution underlying y, i.e. F(u) > G(u) for all u.
            'greater': the distribution underlying x is stochastically greater than the distribution underlying y, i.e. F(u) < G(u) for all u.
            One-sided alternative hypothesis tests median from one group can be greater or lesser than other group.
            """
            mv_pvalue = mw_results[1]

            if mv_pvalue <= self.subnet_pvalue_cutoff:
                regulator_activation_score = 0.0
                effect_target_counter = 0
                for target in targets:
                    target_exp_value = target[annotation][0][0]
                    target_exp_pvalue = target[annotation][0][1]
                    if target_exp_pvalue < 0.05 or (str(target_exp_pvalue) == str(np.nan) and abs(target_exp_value) >= 1.0):
                        reg2target_rels = my_graph._relation4(regulator_id,target['Id'][0])
                        rel = reg2target_rels[0]
                        sign = rel.effect_sign()
                        if sign != 0:
                            regulator_activation_score = regulator_activation_score + sign*target_exp_value
                            effect_target_counter += 1

                subnet_size = len(abs_sub_pd)
                if subnet_size >= self.min_subnet_size and effect_target_counter > 0:
                    regulator_activation_score = regulator_activation_score/math.sqrt(effect_target_counter)
                else:
                    regulator_activation_score = np.nan

                new_prop_name = self.__sample_annotation(sample)
                prp_value = list((regulator_activation_score,mv_pvalue,effect_target_counter))
                nx.set_node_attributes(my_graph, {regulator_id:{new_prop_name:prp_value}})
                sample_regulator_ids.add(regulator_id)

        sample_time = self.execution_time(sample_start)
        print('Found %d regulators with pvalue < %.2f in %s'
                % (len(sample_regulator_ids),self.subnet_pvalue_cutoff,sample_time))

        return sample_regulator_ids


    def expression_regulators(self,regulomes:dict=dict(),from_graph=ResnetGraph(),prefix4target_annotation=''):
        """
        # all samples from self.__my_sample_names__ are used for annotation
        Returns
        -------
        [PSObject] with expression regulators annotated with tuple (activity, pvalue, # valid targets) for each sample\n
        Tuple is stored in PSObject['activation in ' + sample_name] annotation field
        """
        start_time = time.time()
        my_graph = from_graph if from_graph else self.Graph
        samples = self.experiment.get_samples(self.__my_sample_names__)
        my_regulomes = regulomes if regulomes else my_graph.regulome_dict(['Protein','FunctionalClass','Complex'],min_size=2)

        regulator_id_list = set()
        for counter, sample in enumerate(samples):
            print('Finding regulators for %s sample (%d out of %d)' % (sample['Name'][0],counter+1,len(samples)))
            sample_regulator_ids = self.regulators4sample(sample,my_regulomes,my_graph,prefix4target_annotation)
            regulator_id_list.update(sample_regulator_ids)

        print('SNEA execution time: %s' % self.execution_time(start_time))
        return my_graph._get_nodes(list(regulator_id_list))


    def __activity_columns(self):
        return [c for c in self.activity_matrix.columns if c[:13] == 'activation in']


    def activity_df(self, annotated_regulators:list):
        node_annotation_names = list()
        column_list = ['Regulator','URN','Average activation score','# samples',NUMBER_OF_TARGETS]
        for s in self.__my_samples__():
            annotation_name = 'activation in '+s['Name'][0]
            node_annotation_names.append(annotation_name)
            column_list += [annotation_name,'activation pvalue in '+s['Name'][0],'# valid targets in '+s['Name'][0]]

        self.activity_matrix = df(columns=column_list)
        row = 0
        for regulator in annotated_regulators:
            self.activity_matrix.loc[row,'URN'] = regulator['URN'][0]
            self.activity_matrix.loc[row,NUMBER_OF_TARGETS] = regulator[NUMBER_OF_TARGETS]
            for a in node_annotation_names:
                try:
                    annotation = regulator[a]
                    self.activity_matrix.loc[row,a] = annotation[0]
                    sampl_name = str(a)[14:]

                    pval_col = 'activation pvalue in '+sampl_name
                    pval = '{:12.5e}'.format(annotation[1]) if annotation[1] != np.nan else ''
                    self.activity_matrix.loc[row,pval_col] = pval

                    measured_targets_col = '# valid targets in ' + sampl_name
                    measured_targets_val = str(annotation[2]) if annotation[2] > 0 else ''
                    self.activity_matrix.loc[row,measured_targets_col] = measured_targets_val
                except KeyError:
                    continue
            row += 1

        activation_scores_columns = self.__activity_columns()
        activation_scores_only = self.activity_matrix[activation_scores_columns]
        averages = activation_scores_only.mean(axis=1,skipna=True)
        self.activity_matrix['Average activation score'] = averages

        pvalue_columns = [c for c in self.activity_matrix.columns if c[:17] == 'activation pvalue']
        pvalues_only = self.activity_matrix[pvalue_columns]
        self.activity_matrix['# samples'] = pvalues_only.count(axis = 1)

        self.activity_matrix.dropna(how='all',subset=['Average activation score'],inplace=True)
        self.activity_matrix.sort_values(['# samples','Average activation score'],ascending=False,inplace=True,ignore_index=True)
        self.activity_matrix['Regulator'] = self.activity_matrix['URN'].apply(lambda x: self.Graph.urn2obj[x]['Name'][0])

        print('Generated report table with %d regulators for %d samples in "%s" experiment' %
                            (len(annotated_regulators), len(self.__my_sample_names__), self.experiment.name()))

        return self.activity_matrix
      

    def find_drugs(self):
        '''
        Output
        ------
        self.sample2drugs

        Returns
        -------
        urn2drugs - {urn:[PSObject]}
        drugid2urns - {database Id:[URNs]}
        '''
        print('Finding drugs inhibiting differential expression for each sample')
        drug_target_expression_network = self.__load_network(['SmallMol'],DPERNET)
        DPERNET4experiment = self.network4experiment(drug_target_expression_network)
        self.experiment.annotate(DPERNET4experiment,self.__my_sample_names__)

        self.sample2drugs = dict()
        all_drugs = set()
        drug2proteins_regulomes = DPERNET4experiment.regulome_dict(['SmallMol'],min_size=2)
        annotated_regulator_drugs = self.expression_regulators(drug2proteins_regulomes,DPERNET4experiment)
        for sample in self.__my_samples__():
            sample_annotation_prop = self.__sample_annotation(sample)
            inhibitors4sample = [drug for drug in annotated_regulator_drugs if drug.get_prop(sample_annotation_prop,0,0)<0]
            #sample_drugs = self.__drugs4sample(s,DPERNET4experiment,drug2proteins_regulomes)
            self.sample2drugs[sample.name()] = inhibitors4sample
            all_drugs.update(inhibitors4sample)

        drug_urns = [d.urn() for d in all_drugs]
        oql_query = 'SELECT Entity WHERE URN = ({props})'
        req_name = 'Finding drugs in the database'
        drugs_entity_graph = self.iterate_oql(oql_query,drug_urns,use_cache=False,request_name=req_name,iterate_ids=False)
        drugs_in_database = drugs_entity_graph._get_nodes()

        drug_ids = [drug.id() for drug in drugs_in_database]
        return drug_ids


    def make_drugs_df(self):
        drug_ids = self.find_drugs() 
        df2map = df.copy_df(self.activity_matrix,only_columns=['Regulator'],rename2={'Regulator':'Name'})
        regulator_urns = self.activity_matrix['URN'].to_list()
        regulators_objtypes = {self.Graph.urn2obj[urn]['ObjTypeName'][0] for urn in regulator_urns}
        input_parameters = {'experiment':self.experiment.name(),
                            'add_bibliography' : False,
                            'consistency_correction4target_rank' : False,
                            'target_types' : list(regulators_objtypes)
                            }

        dt = Drugs4Targets(self.APIconfig,input_parameters,clone=False,what2retrieve=NO_REL_PROPERTIES)
        dt.add_rel_props([EFFECT])
        print('Loading targets')
        mapped_targets_df = dt.load_pandas(from_entity_df=df2map, max_children_count=11)
        all_target_ids  = dt._all_ids(mapped_targets_df)

        dt.load_drug_graph(all_target_ids,drug_ids)
        dt.columns2drop = [dt.__temp_id_col__,dt.__resnet_name__,dt.__mapped_by__]
        dt.add_targets2drugs_ws = False

        drugs_df = df(columns=['Name',PHARMAPENDIUM_ID])
        activation_scores_columns = self.__activity_columns()
        for i, col in enumerate(activation_scores_columns):
            sample_name = col[14:]
            sample_drugs = self.sample2drugs[sample_name]
            dt.init_drug_df(sample_drugs)
            print('Ranking %d drugs for %s sample' % (len(sample_drugs),sample_name))
            
            dt.params['sample'] = sample_name
            rank_df = df(self.activity_matrix[['Regulator',col]])
            rank_df.rename(columns={col:RANK,'Regulator':'Name'},inplace=True)
            ranked_df = mapped_targets_df.merge_df(rank_df,how='left',on='Name')
            
            targets4antagonist_df = df(ranked_df[ranked_df[RANK] > 0.0])
            targets4antagonist_df._name_ = ANTAGONIST_TARGETS_WS
            
            targets4agonist_df = df(ranked_df[ranked_df[RANK] < 0.0])
            targets4agonist_df[RANK] = targets4agonist_df[RANK].abs()
            targets4agonist_df._name_ = AGONIST_TARGETS_WS
            
            dt.add2report(targets4antagonist_df)
            dt.add2report(targets4agonist_df)
            
            drugs_df4sample = dt.score_drugs(normalize=False)

            name2rank = drugs_df4sample.to_dict('Name',DRUG2TARGET_REGULATOR_SCORE)
            #name2pharmapendium = drugs_df4sample.to_dict('Name',PHARMAPENDIUM_ID)
            new_rank_col = DRUG2TARGET_REGULATOR_SCORE +' in '+sample_name
            drugs_df = drugs_df.merge_dict(name2rank, new_rank_col, 'Name', add_all=True)
            drugs_df = drugs_df.add_values(PHARMAPENDIUM_ID,drugs_df4sample,PHARMAPENDIUM_ID,'Name')
            print(f'Added ranks for regulators of {sample_name} sample')

        rank_columns = [c for c in drugs_df.columns.to_list() if c[:len(DRUG2TARGET_REGULATOR_SCORE)]==DRUG2TARGET_REGULATOR_SCORE]
        drugs_df.not_null_counts(rank_columns)
        drugs_df['Average regulator score'] = self.drugs_df[rank_columns].mean(axis=1)
        drugs_df.sort_values(by=['Row count','Average regulator score'],ascending=False,inplace=True)
        self.drugs_df = df(drugs_df[[PHARMAPENDIUM_ID,'Name','Row count','Average regulator score']+rank_columns])
        self.drugs_df.copy_format(drugs_df)
        self.drugs_df._name_ = 'Drugs'
        return drugs_df


    def report_path(self,extension:str):
        return self.data_dir+self.experiment.name()+' regulators.'+extension


    def report(self):
        self.activity_matrix.fillna('', inplace=True)
        fout = self.report_path('xlsx')
        writer = ExcelWriter(fout, engine='xlsxwriter')
        abs_scores = self.activity_matrix['Average activation score'].abs()
        max_abs = abs_scores.max()
        cond_format = {
                        'type': '3_color_scale',
                        'min_value': -max_abs,
                        'max_value': max_abs,
                        'mid_value': 0,
                        'min_color': 'blue',
                        'max_color': 'red',
                        'mid_color': 'white',
                        'mid_type': "num",
                        'min_type': "num",
                        'max_type': "num"
                        }
        self.activity_matrix.set_conditional_frmt(cond_format,'Average activation score','Average activation score')
        self.activity_matrix.make_header_vertical()
        self.activity_matrix.df2excel(writer,'activity')

        if not self.drugs_df.empty:
            self.drugs_df.make_header_vertical()
            self.drugs_df.df2excel(writer,'Drugs')

        writer.save()
        print('Scored regulators are in %s file' % fout)
        

def test_run(APIconfig:dict,fast=True, find_drugs=True):
    no_outliers_samples = list()
    for i in range(1,33):
        no_outliers_samples.append(i-1)
    no_outliers_samples.pop(14)
    no_outliers_samples.pop(8)
    no_outliers_samples.pop(0)

    test_samples = [2,3]
    input = test_samples if fast else no_outliers_samples
    snea = SNEA(APIconfig,'PNOC003vsGSE120046',sample_ids=input)
    snea.set_dir('D:/Python/PBTA/')
    annotated_regulators = snea.expression_regulators()
    snea.activity_df(annotated_regulators)
    if find_drugs:
        snea.make_drugs_df()
    snea.report()

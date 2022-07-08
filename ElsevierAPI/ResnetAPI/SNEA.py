from .ResnetAPISession import APISession, time, ResnetGraph, PSObject, nx, math
from ..pandas.panda_tricks import df, ExcelWriter
from .Zeep2Experiment import Experiment
import scipy.stats as stats
import numpy as np


class SNEA(APISession):
    def __init__(self, APIconfig, experiment_name:str):
        super().__init__(APIconfig['ResnetURL'], APIconfig['PSuserName'], APIconfig['PSpassword'])
        self.PageSize = 1000
        self.DumpFiles = []
        self.APIconfig = APIconfig
        self.experiment = Experiment.download(experiment_name,APIconfig)
        self.__load_network() #expression_network is in self.Graph
        self.subnetworks = dict() # {regulator_id:[target_ids]}
        self.activation_pd = df()
        self.subnet_pvalue_cutoff = 0.05
        self.min_subnet_size = 2
        

    def __load_network(self)->'ResnetGraph':
        network_name = 'protein_expression_network'
        cache_path = 'ElsevierAPI/ResnetAPI/__pscache__/'
        cache_file = cache_path+network_name+'.rnef'
        start = time.time()
        try:
           # print('Loading %s from %s file' % (network_name,cache_file))
            self.Graph = ResnetGraph.fromRNEF(cache_file)
            #self.Graph = self.Graph.make_simple()
            self.Graph.name = network_name
            print('%s with %d edges and %d nodes was loaded in %s' 
                % (network_name,self.Graph.number_of_edges(),self.Graph.number_of_nodes(),self.execution_time(start)))
            return self.Graph
        except FileNotFoundError:
            pass

        print('%s was not found.  Begin downloading relations from database' % cache_file)
        self.relProps = ['Effect','RelationNumberOfReferences','URN']
        self.PageSize = 1000
        self.DumpFiles = []
        get_proteins = 'SELECT Entity WHERE objectType = (Protein,FunctionalClass,Complex)'
        oql_query = 'SELECT Relation WHERE objectType = (Expression,PromoterBinding) AND NeighborOf upstream ({n1}) AND NeighborOf downstream ({n2})'
        oql_query = oql_query.format(n1=get_proteins, n2=get_proteins)
        #oql_query = 'SELECT Relation WHERE objectType = Metabolization' # for testing
        self.process_oql(oql_query, 'Fetching protein expression network')
        self.Graph = self.Graph.make_simple()
        self.graph2rnef(cache_path+network_name+'.rnef')
        self.Graph.name = network_name
        print('%s with %d edges and %d nodes was downloaloaded from database in %s' 
                % (network_name,self.Graph.number_of_edges(),self.Graph.number_of_nodes(),self.execution_time(start)))
        return self.Graph


    def expression_regulators(self, sample_names=[],sample_ids=[]):
        """
        returns list of PSObjects for expression regulators annotated with tuple (activity, pvalue, # valid targets) for each sample\n
        annotation field is called 'activity in + sample name'
        """
        start_time = time.time()
        expname = self.experiment.name()
        urns = self.experiment.urns()
        print('"%s" experiment has %d mapped entities' % (expname,len(urns)))

        network_name = 'expression regulatory network for "{}"'.format(expname)
        expression_network = self.Graph.regulatory_network(urns,network_name)
        print('%s experiment has %d edges and %d nodes' 
                %(network_name,expression_network.number_of_edges(),expression_network.number_of_nodes()))
        
        self.experiment.annotate_graph(expression_network,sample_names,sample_ids)
        regulomes =  expression_network.regulome_dict(self.min_subnet_size)

        sample_counter = 0
        samples = self.experiment.get_samples(sample_names,sample_ids)

        regulator_id_list = set()
        for sample in samples:
            sample_counter += 1
            sample_name = sample['Name'][0]
            attribute_name = self.experiment.name4annotation(sample_name)
            print('Finding regulators for %s sample (%d out of %d)' 
                % (sample_name,sample_counter,len(samples)))
            regulator_counter = 0
            sample_start = time.time()
            abs_sample_distribution = sample.data['value'].abs()
            for regulator_id, targets in regulomes.items():
                subnetwork_values = [t[attribute_name][0][0] for t in targets]
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
                    regulator_counter += 1
                    regulator_activation_score = 0.0
                    effect_target_counter = 0
                    regulator = expression_network._get_node(regulator_id)
                    if regulator['Name'][0] == 'SPRY2':
                        print('')
                    for target in targets:
                        target_exp_value = target[attribute_name][0][0]
                        target_exp_pvalue = target[attribute_name][0][1]
                        if target_exp_pvalue < 0.05 or (str(target_exp_pvalue) == str(np.nan) and abs(target_exp_value) >= 1.0):
                            reg2target_rels = expression_network._relation4(regulator_id,target['Id'][0])
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

                    new_prop_name = 'activation in ' + sample_name
                    prp_value = list((regulator_activation_score,mv_pvalue,effect_target_counter))
                    nx.set_node_attributes(expression_network, {regulator_id:{new_prop_name:prp_value}})
                    regulator_id_list.add(regulator_id)

            sample_time = self.execution_time(sample_start)
            print('Found %d regulators with pvalue < %.2f in %s'
                    % (regulator_counter,self.subnet_pvalue_cutoff,sample_time))

        print('SNEA execution time: %s' % self.execution_time(start_time))
        return expression_network._get_nodes(list(regulator_id_list))


    def make_report_pd(self, annotated_regulators:list, sample_names=[],sample_ids=[]):
        node_annotation_names = list()
        column_list = ['Regulator','URN','Average activation score','# samples','# targets']
        samples = self.experiment.get_samples(sample_names,sample_ids)
        for s in samples:
            annotation_name = 'activation in '+s['Name'][0]
            node_annotation_names.append(annotation_name)
            column_list += [annotation_name,'activation pvalue in '+s['Name'][0],'# valid targets in '+s['Name'][0]]

        self.activation_pd = df(columns=column_list)
        row = 0
        for regulator in annotated_regulators:
            self.activation_pd.loc[row,'URN'] = regulator['URN'][0]
            self.activation_pd.loc[row,'# targets'] = regulator['# targets'][0]
            for a in node_annotation_names:
                try:
                    annotation = regulator[a]
                    self.activation_pd.loc[row,a] = annotation[0]
                    sampl_name = str(a)[14:]

                    pval_col = 'activation pvalue in '+sampl_name
                    pval = '{:12.5e}'.format(annotation[1]) if annotation[1] != np.nan else ''
                    self.activation_pd.loc[row,pval_col] = pval

                    measured_targets_col = '# valid targets in ' + sampl_name
                    measured_targets_val = str(annotation[2]) if annotation[2] > 0 else ''
                    self.activation_pd.loc[row,measured_targets_col] = measured_targets_val
                except KeyError:
                    continue
            row += 1

        activation_scores_columns = [c for c in self.activation_pd.columns if c[:13] == 'activation in']
        activation_scores_only = self.activation_pd[activation_scores_columns]
        averages = activation_scores_only.mean(axis=1,skipna=True)
        self.activation_pd['Average activation score'] = averages

        pvalue_columns = [c for c in self.activation_pd.columns if c[:17] == 'activation pvalue']
        pvalues_only = self.activation_pd[pvalue_columns]
        self.activation_pd['# samples'] = pvalues_only.count(axis = 1)

        self.activation_pd.dropna(how='all',subset=['Average activation score'],inplace=True)
        self.activation_pd.sort_values(['# samples','Average activation score'],ascending=False,inplace=True,ignore_index=True)
        self.activation_pd['Regulator'] = self.activation_pd['URN'].apply(lambda x: self.Graph.urn2obj[x]['Name'][0])

        print('Report table for %d regulators in %d samples from "%s" experiment was generated' %
                            (len(annotated_regulators), len(samples), self.experiment.name()))

        return self.activation_pd


    def report(self):
        self.activation_pd.fillna('', inplace=True)
        fout = self.__expname()+' regulators.xlsx'
        writer = ExcelWriter(fout, engine='xlsxwriter')
        abs_scores = self.activation_pd['Average activation score'].abs()
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
        self.activation_pd.df2excel(writer,'activity',vertical_header=True, 
                        cond_format=cond_format,for_areas=['C:C'])
        writer.save()
        print('Scored regulators are in %s file' % fout)


    def test_run(self, fast=True):
        #from ElsevierAPI import load_api_config
        #expname = 'PNOC003vsGSE120046'
        #snea = SNEA(load_api_config(),expname)

        no_outliers_samples = list()
        for i in range(1,33):
            no_outliers_samples.append(i-1)
        no_outliers_samples.pop(14)
        no_outliers_samples.pop(8)
        no_outliers_samples.pop(0)

        test_samples = [2,3]
        input = test_samples if fast else no_outliers_samples
        self.expression_regulators(sample_ids=input)
        self.report()

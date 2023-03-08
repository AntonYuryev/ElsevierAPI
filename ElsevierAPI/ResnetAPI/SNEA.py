from .ResnetAPISession import APISession, time, math, NO_REL_PROPERTIES
from .ResnetGraph import REFCOUNT, ResnetGraph, nx, EFFECT,NUMBER_OF_TARGETS,PROTEIN_TYPES
from ..pandas.panda_tricks import df, ExcelWriter
from .Drugs4Disease import Drugs4Targets
from .Drugs4Disease import ANTAGONIST_TARGETS_WS,AGONIST_TARGETS_WS,RANK,DRUG2TARGET_REGULATOR_SCORE,PHARMAPENDIUM_ID
from .Zeep2Experiment import Experiment, Sample
import scipy.stats as stats
import numpy as np
from concurrent.futures import ThreadPoolExecutor
from contextlib import redirect_stdout

PPERNET = 'protein_expression_regulatory_network'
DPERNET = 'drug-protein expression regulatory network'
ACTIVATION_IN = 'activation in '
EXPRESSION = ['Expression','PromoterBinding']


class SNEA(APISession):
    __regulators__ = dict() # = {sample_name:[regulator_uids]}}
    __sample2drugs__ = dict() # = {sample_name:[PSObject]}} stores only drugs having negative activation score
    __my_psobjs__ = list()

    def __init__(self,*args,**kwargs):
        '''
        Input
        -----
        experiment_name_or_path2file=args[0]\n
        APIconfig=args[1] - optional

        sample_names - deafult []\n
        sample_ids - default []\n
        find_drugs - default False\n
        connect2server - default True, set to False to run script using data in __pscache__ files instead of database
        '''
        my_kwargs = {
            'data_dir':'',
            'no_mess':True,
            'sample_names':[],
            'sample_ids':[],
            'find_drugs':False,
            'connect2server':True
            }
        my_kwargs.update(kwargs)

        APIconfig = dict(args[1]) if my_kwargs['connect2server'] else dict()
        self.set_dir(my_kwargs['data_dir'])

        super().__init__(APIconfig,**my_kwargs)
        self.PageSize = 1000
        self.DumpFiles = []
        self.experiment = Experiment({'Name':[args[0]]})
        if my_kwargs['no_mess']:
            log_name = self.report_path('log')
            print('Runtime messages are written to "%s"'%log_name)
            with open(log_name, 'w') as fout:
                with redirect_stdout(fout):
                    self.__load_data__(**my_kwargs)
        else:
            self.__load_data__(**my_kwargs)
        

    def __load_data__(self,**my_kwargs):
        experiment_name = self.experiment.name()
        if my_kwargs['connect2server']:
            with ThreadPoolExecutor(max_workers=4, thread_name_prefix='Initializing SNEA') as executor:
                future1 = executor.submit(Experiment.download,experiment_name,self.APIconfig)
                future2 = executor.submit(self.__load_network,PPERNET)
                self.experiment = future1.result()
                expression_network = future2.result()
        else:
            exp = Experiment.from_file(experiment_name,**my_kwargs)
            expression_network = self.__load_network(PPERNET)
            mapping_id_name = exp.identifier_name()
            self.experiment = exp.map2(expression_network,mapping_id_name)
            
        if my_kwargs['find_drugs']:
            dt_kwargs = {'experiment':experiment_name,
                    'add_bibliography' : False,
                    'consistency_correction4target_rank' : False,
                    'target_types' : PROTEIN_TYPES,
                    'strict_mode' : False,
                    'clone': False,
                    'what2retrieve' : NO_REL_PROPERTIES,
                    }
            
            dt_kwargs['connect2server'] = my_kwargs.get('connect2server',True)
            if dt_kwargs['connect2server']:
                self.dt = Drugs4Targets(self.APIconfig,dt_kwargs)
            else:
                self.dt = Drugs4Targets.from_files(**dt_kwargs)
            
            e = ThreadPoolExecutor(max_workers=4, thread_name_prefix='Reading drug-target network')
            self.dt_future = e.submit(self.dt.dt_consist.load_cache)

        self.Graph.clear_resnetgraph() # in case graph was downloaded from database
        self.diffexp_pvalue_cutoff = 0.05
        self.experiment = self.experiment.mask_by_pval(self.diffexp_pvalue_cutoff)
        self.Graph = self.experiment.annotated_subnetwork(expression_network)

        my_samples = self.experiment.get_samples(my_kwargs["sample_names"],my_kwargs['sample_ids'])
        self.__my_sample_names__ = [s.name() for s in my_samples]

        self.subnetworks = dict() # {regulator_id:[target_ids]}
        self.activity_matrix = df()
        self.drugs_df = df()
        
        self.subnet_pvalue_cutoff = 0.05
        self.min_subnet_size = 2
          
    
    @classmethod
    def from_files(cls,path2experiment:str,**kwargs):
        '''
        Input
        -----
        path2experiment - args[0]\n
        sample_names - deafult []\n
        sample_ids - default []\n
        find_drugs - default False\n
        data_type from ['LogFC', 'Intensity', 'LogIntensity', 'SignedFC']\n
        identifier_column - default 0\n
        phenotype_rows - default []\n
        header - default 0\n
        has_pvalue - default True\n
        last_sample - column with last sample, default None. If experiment has p-value columns "last_sample" must point to the last p-value column
        '''
        my_kwargs = dict(kwargs)
        my_kwargs['connect2server'] = False
        snea = cls(path2experiment,**my_kwargs)
        return snea


    def __my_samples(self):
        return self.experiment.get_samples(self.__my_sample_names__)


    @staticmethod
    def __sample_annotation(sample:Sample):
        return ACTIVATION_IN+sample.name()


    def __load_network(self,network_name:str)->'ResnetGraph':
        get_regulators = f'SELECT Entity WHERE objectType = ({PROTEIN_TYPES})'
        ent_props = ['Name','Ensembl ID'] # Ensembl ID to support DESeq2 output

        get_targets = f'SELECT Entity WHERE objectType = ({PROTEIN_TYPES})'
        oql_query = f'SELECT Relation WHERE objectType = ({EXPRESSION}) \
            AND NeighborOf upstream ({get_targets}) AND NeighborOf downstream ({get_regulators})'

        rel_props = ['URN',EFFECT,REFCOUNT]
        return self.load_cache(network_name,[oql_query],ent_props,rel_props)


    def __regulators4sample(self,sample:Sample,regulomes:dict,from_graph=ResnetGraph(),target_annotation_prefix=''):
        '''
        Input
        -----
        regulomes = {regulator_id : [targets]} made by ResnetGraph.regulome_dict()

        Return
        ------
        [uids] - uids of regulators 
        regulator objects in self.Graph are annotated with SNEA activation score and p-value calculated from "sample" 
        '''
        my_graph = from_graph if from_graph else self.Graph

        sample_start = time.time()
        abs_sample_distribution = sample.data['value'].abs()
        abs_sample_distribution = abs_sample_distribution.dropna(how='all')
        sample_regulator_uids = set()
        annotation = target_annotation_prefix+self.experiment.name4annotation(sample)
        for regulator_uid, targets in regulomes.items():
            subnetwork_values = list()
            annotated_targets = list()
            for t in targets:
                try:
                    subnetwork_values.append(t[annotation][0][0])
                    annotated_targets.append(t)
                except KeyError:
                    continue
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
                for target in annotated_targets:
                    target_exp_value = target[annotation][0][0]
                    target_exp_pvalue = target[annotation][0][1]
                    if target_exp_pvalue < 0.05 or (str(target_exp_pvalue) == str(np.nan) and abs(target_exp_value) >= 1.0):
                        reg2target_rels = my_graph._psrels4(regulator_uid,target.uid())
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
                nx.set_node_attributes(my_graph, {regulator_uid:{new_prop_name:prp_value}})
                sample_regulator_uids.add(regulator_uid)

        sample_time = self.execution_time(sample_start)
        print('Found %d regulators with pvalue < %.2f in %s'
                % (len(sample_regulator_uids),self.subnet_pvalue_cutoff,sample_time),flush=True)

        return sample_regulator_uids # must return uids and not PSObjects here 


    def expression_regulators(self,regulomes=dict(),from_graph=ResnetGraph(),prefix4target_annotation=''):
        """
        Input
        -----
        all samples from self.__my_sample_names__ are used for annotation

        Returns
        -------
        self.__regulators__ = {sample_name:[regulator_uids]} with expression regulators annotated with tuple (activity, pvalue, # valid targets) for each sample\n
        Tuple is stored in PSObject['activation in ' + sample_name] annotation field in self.Graph
        self.__regulators__ is used in self.make_drug_df
        self.activity_df uses all annotated regulators in self.__regulators__
        """
        start_time = time.time()
        my_graph = from_graph if from_graph else self.Graph
        samples = self.experiment.get_samples(self.__my_sample_names__)
        my_regulomes = regulomes if regulomes else my_graph.regulome_dict(PROTEIN_TYPES,min_size=2)

        for counter, sample in enumerate(samples):
            print('Finding regulators for %s sample (%d out of %d)' % (sample['Name'][0],counter+1,len(samples)))
            sample_regulators_uids = self.__regulators4sample(sample,my_regulomes,my_graph,prefix4target_annotation)
            self.__regulators__[sample.name()] = sample_regulators_uids

        print('SNEA execution time: %s' % self.execution_time(start_time))
        return
    

    def __regulators(self, sample_name='', from_graph=ResnetGraph()):
        my_graph = from_graph if from_graph else self.Graph
        all_regulators_uids = set()
        if sample_name:
            all_regulators_uids.update(self.__regulators__[sample_name])
        else:
            [all_regulators_uids.update(reg_uids) for reg_uids in self.__regulators__.values()]
        return  my_graph._get_nodes(all_regulators_uids)


    def activity_df(self):
        node_annotation_names = list()
        column_list = ['Regulator','URN','Average activation score','# samples',NUMBER_OF_TARGETS]
        for sample in self.__my_samples():
            s_n = sample.name()
            annotation_name = 'activation in '+ s_n
            node_annotation_names.append(annotation_name)
            column_list += [annotation_name,'activation pvalue in '+s_n,'# valid targets in '+s_n]

        self.activity_matrix = df(columns=column_list)
        row = 0
        annotated_regulators = self.__regulators()
        for regulator in annotated_regulators:
            self.activity_matrix.loc[row,'URN'] = regulator.urn()
            self.activity_matrix.loc[row,'Regulator'] = regulator.name()
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
  
        print('Generated report table with %d regulators for %d samples in "%s" experiment' %
                            (len(annotated_regulators), len(self.__my_sample_names__), self.experiment.name()),flush=True)

        return self.activity_matrix
      
    
    def __activity_columns(self):
        return [c for c in self.activity_matrix.columns if c[:13] == 'activation in']


    def find_drugs(self):
        '''
        Output
        ------
        self.__regulators__

        Returns
        -------
        [PSObject]
        '''
        print('Finding drugs inhibiting differential expression for each sample')
        drug_target_network = self.dt_future.result()
        drug_target_expression_network = drug_target_network.subgraph_by_relprops(['Expression'])
        #DPERNET4experiment = self.network4experiment(drug_target_expression_network)
        DPERNET4experiment = self.experiment.annotated_subnetwork(drug_target_expression_network,self.__my_sample_names__)

        all_drugs = set()
        drug2proteins_regulomes = DPERNET4experiment.regulome_dict(['SmallMol'],min_size=2)
        self.expression_regulators(drug2proteins_regulomes,DPERNET4experiment)
        for sample in self.__my_samples():
            drugs4sample = self.__regulators(sample.name(),from_graph=DPERNET4experiment)
            sample_annotation_prop = self.__sample_annotation(sample)
            drugs_inhibiting_sample = [drug for drug in drugs4sample if drug.get_prop(sample_annotation_prop,0,0)<0]
            self.__sample2drugs__[sample.name()] = drugs_inhibiting_sample
            all_drugs.update(drugs_inhibiting_sample)

        self.Graph = self.Graph.compose(DPERNET4experiment) 
        # to make sure all ranked targets are in self.Graph which is required by Drugs4Targets::load_target_ranks() 
        return all_drugs


    def make_drugs_df(self):
        self.dt.__targets__ = self.__regulators()
        self.find_drugs()
        # SNEA is run using RNEF files only no dbids required
        mapped_targets_df = self.dt.load_df(self.dt.__targets__)
        self.dt.columns2drop = [self.dt.__temp_id_col__,self.dt.__resnet_name__,self.dt.__mapped_by__]
            
        self.dt_future.result()
        drugs_df = df(columns=['Name',PHARMAPENDIUM_ID])
        activation_scores_columns = self.__activity_columns()
        if len(activation_scores_columns) > 1:
            self.dt.add_targets2drugs_ws = False
        else:
            self.dt.add_targets2drugs_ws = True

        for i, col in enumerate(activation_scores_columns):
            sample_name = col[14:]
            drugs_inhibit_sample = self.__sample2drugs__[sample_name]
            self.dt.load_drugs(limit2drugs=drugs_inhibit_sample)
            self.dt.init_drug_df(drugs_inhibit_sample)
            print('Ranking %d drugs for %s sample' % (len(drugs_inhibit_sample),sample_name))
            
            self.dt.params['sample'] = sample_name
            target_df = df(self.activity_matrix[['Regulator',col]])
            target_df.rename(columns={col:RANK,'Regulator':'Name'},inplace=True)
            target_df = mapped_targets_df.merge_df(target_df,how='left',on='Name')
            
            targets4antagonist_df = df(target_df[target_df[RANK] > 0.0])
            targets4antagonist_df._name_ = ANTAGONIST_TARGETS_WS
            
            targets4agonist_df = df(target_df[target_df[RANK] < 0.0])
            targets4agonist_df[RANK] = targets4agonist_df[RANK].abs()
            targets4agonist_df._name_ = AGONIST_TARGETS_WS
            
            self.dt.add2report(targets4antagonist_df)
            self.dt.add2report(targets4agonist_df)
            
            drugs_df4sample = self.dt.score_drugs(normalize=False)
            name2rank = drugs_df4sample.to_dict('Name',DRUG2TARGET_REGULATOR_SCORE)

            new_rank_col = DRUG2TARGET_REGULATOR_SCORE +' in '+sample_name
            drugs_df = drugs_df.merge_dict(name2rank, new_rank_col, 'Name', add_all=True)
            drugs_df = drugs_df.add_values(from_column=PHARMAPENDIUM_ID,
                                            in_df=drugs_df4sample,
                                            map_by_my_col='Name')
            if self.dt.add_targets2drugs_ws:
                drugs_df = drugs_df.add_values(from_column='Directly inhibited targets',
                                                in_df=drugs_df4sample,
                                                map_by_my_col='Name')
                drugs_df = drugs_df.add_values(from_column='Directly activated targets',
                                                in_df=drugs_df4sample,
                                                map_by_my_col='Name')
                drugs_df = drugs_df.add_values(from_column='Indirectly inhibited targets',
                                                in_df=drugs_df4sample,
                                                map_by_my_col='Name')
                drugs_df = drugs_df.add_values(from_column='Indirectly activated targets',
                                                in_df=drugs_df4sample,
                                                map_by_my_col='Name')


            self.dt.target_uid2rank.clear()
            print(f'Added ranks for regulators of {sample_name} sample')

        rank_columns = [c for c in drugs_df.columns.to_list() if c[:len(DRUG2TARGET_REGULATOR_SCORE)]==DRUG2TARGET_REGULATOR_SCORE]
        drugs_df.not_nulls(rank_columns)
        drugs_df['Average regulator score'] = drugs_df[rank_columns].mean(axis=1)
        drugs_df[RANK] = drugs_df['Average regulator score'] * drugs_df['Row count']/len(rank_columns)
        drugs_df.sort_values(by=[RANK],ascending=False,inplace=True)
        column_order = [PHARMAPENDIUM_ID,'Name',RANK,'Row count','Average regulator score']+rank_columns
        if self.dt.add_targets2drugs_ws:
            column_order += ['Directly inhibited targets','Directly activated targets','Indirectly inhibited targets','Indirectly activated targets']
        self.drugs_df = df(drugs_df[column_order])
        self.drugs_df.copy_format(drugs_df)
        self.drugs_df._name_ = 'Drugs'
        return self.drugs_df


    def report_path(self,extension:str):
        return self.data_dir+self.experiment.name()+' SNEA.'+extension


    def report(self,to_file=''):
        self.activity_matrix.fillna('', inplace=True)
        fout = self.data_dir+to_file+'.xlsx' if to_file else self.report_path('xlsx')
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

    test_samples = [1,2]
    input = test_samples if fast else no_outliers_samples
    experiment_name = 'PNOC003vsGSE120046'
    snea = SNEA(APIconfig,experiment_name,sample_ids=input,find_drugs=find_drugs)
    snea.set_dir('D:/Python/PBTA/')
    snea.expression_regulators()
    snea.activity_df()
    if find_drugs:
        snea.make_drugs_df()
    snea.report(experiment_name+'test')

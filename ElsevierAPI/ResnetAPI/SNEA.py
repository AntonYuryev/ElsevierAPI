from .ResnetAPISession import time, math, NO_REL_PROPERTIES
from .ResnetAPIcache import APIcache
from .ResnetGraph import REFCOUNT, ResnetGraph, nx, EFFECT,NUMBER_OF_TARGETS,PROTEIN_TYPES,execution_time
from ..pandas.panda_tricks import df, ExcelWriter
from .Drugs4Disease import Drugs4Targets
from .Drugs4Disease import ANTAGONIST_TARGETS_WS,AGONIST_TARGETS_WS,RANK,DRUG2TARGET_REGULATOR_SCORE,PHARMAPENDIUM_ID
from .Zeep2Experiment import Experiment, Sample
from scipy.stats._mannwhitneyu import mannwhitneyu
import numpy as np
from concurrent.futures import ThreadPoolExecutor
from contextlib import redirect_stdout

PPERNET = 'protein_expression_regulatory_network'
DPERNET = 'drug-protein expression regulatory network'
ACTIVATION_IN = 'activation in '
EXPRESSION = ['Expression','PromoterBinding']


class SNEA(APIcache):

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
            'no_mess':True,
            'sample_names':[],
            'sample_ids':[],
            'find_drugs':False,
            'connect2server':True,
            'load_cache' : False
            #'limitDrugs2expression_regulators':True
            }
        my_kwargs.update(kwargs)
        session_kwargs, self.params = self.split_kwargs(my_kwargs)

        #APIconfig = dict(args[1]) if my_kwargs['connect2server'] else dict()
        self.diffexp_pvalue_cutoff = self.params.pop('max_diffexp_pval',0.05)
        #self.limit2expr_regs = self.params.pop('limitDrugs2expression_regulators',True)

        self.__regulators__ = dict() # = {sample_name:[regulator_uids]}}
        self.__sample2drugs__ = dict() # = {sample_name:[PSObject]}} stores only drugs having negative activation score
        self.__my_psobjs__ = list()

        self.experiment_dir = str(self.params.pop('experiment_dir',""))
        if self.experiment_dir[-1:] != '/':
            self.experiment_dir += '/'

        self.experiment = Experiment({'Name':[args[0]]})
        if my_kwargs['no_mess']:
            log_name = self.data_dir+self.report_path('log')
            print(f'Runtime messages for SNEA initialization will be added to {log_name}')
            with open(log_name, 'a') as fout:
                with redirect_stdout(fout):
                   super().__init__(**my_kwargs)
                   self.__load_data__(**my_kwargs)
        else:
            super().__init__(**my_kwargs)
            self.__load_data__(**my_kwargs)
        

    def __load_data__(self,**my_kwargs):
        experiment_name = self.experiment.name()
        if my_kwargs['connect2server']:
            with ThreadPoolExecutor(max_workers=4, thread_name_prefix='Initializing SNEA') as executor:
                future1 = executor.submit(Experiment.download,experiment_name,self.APIconfig)
                future2 = executor.submit(self.__load_network,PPERNET,**my_kwargs)
                self.experiment = future1.result()
                expression_network = future2.result()
        else:
            exp = Experiment.from_file(experiment_name,**my_kwargs)
            expression_network = self.__load_network(PPERNET,**my_kwargs)
            mapping_id_name = str(exp.identifier_name())
            self.experiment = exp.map2(expression_network,mapping_id_name)

        if my_kwargs['find_drugs']:
            dt_kwargs = {'experiment':experiment_name,
                    'add_bibliography' : False,
                    'consistency_correction4target_rank' : False,
                    'target_types' : PROTEIN_TYPES,
                    'strict_mode' : False,
                    'clone': False,
                    'what2retrieve' : NO_REL_PROPERTIES,
                    'data_dir' : 'D:/Python/ENTELLECT_API/ElsevierAPI/ResnetAPI/__pscache__/'
                    }      
            
            dt_kwargs['connect2server'] = my_kwargs.get('connect2server',True)
            dt_kwargs['load_cache'] = False
            if dt_kwargs['connect2server']:
                self.dt = Drugs4Targets(**dt_kwargs)
            else:
                self.dt = Drugs4Targets.from_files(**dt_kwargs)
            
            e = ThreadPoolExecutor(max_workers=4, thread_name_prefix='Reading drug-target network')
            dt_kwargs['cache_name'] = 'drug2target'
            self.dt_future = e.submit(self.dt.dt_consist._load_cache,**dt_kwargs)

        self.Graph.clear_resnetgraph() # in case graph was downloaded from database
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


    def __load_network(self,network_name:str, **kwargs)->'ResnetGraph':
        get_regulators = f'SELECT Entity WHERE objectType = ({PROTEIN_TYPES})'
        ent_props = ['Name','Ensembl ID'] # Ensembl ID to support DESeq2 output

        get_targets = f'SELECT Entity WHERE objectType = ({PROTEIN_TYPES})'
        oql_query = f'SELECT Relation WHERE objectType = ({EXPRESSION}) \
            AND NeighborOf upstream ({get_targets}) AND NeighborOf downstream ({get_regulators})'

        get_regulators = f'SELECT Entity WHERE objectType = ({PROTEIN_TYPES})'
        get_targets = f'SELECT Entity WHERE objectType = ({PROTEIN_TYPES})'
        oql_query = f'SELECT Relation WHERE objectType = ({EXPRESSION}) \
            AND NeighborOf upstream ({get_targets}) AND NeighborOf downstream ({get_regulators})'
        
        my_kwargs = dict(kwargs)
        my_kwargs['oql_queries'] = [(oql_query,f'downloading {network_name}')]
        my_kwargs['ent_props'] = ['Name','Ensembl ID'] # Ensembl ID to support DESeq2 output
        my_kwargs['rel_props'] = ['URN',EFFECT,REFCOUNT]
        my_kwargs['cache_name'] = network_name
        return self._load_cache(**my_kwargs)


    def activity(self,_in:Sample,_4reg_uid:int,with_targets:list,according2:ResnetGraph):
            '''
            Input
            -----
            with_targets - [PSObject]"annotated_with_sample"

            Return
            ------
            regulator_activation_score,regulome_values,effect_target_counter
            '''
            annotated_with_sample = self.experiment.name4annotation(_in)
            regulator_activation_score = 0.0
            effect_target_counter = 0
            regulome_values = list()
            for target in with_targets:
                try:
                    target_exp_value, target_exp_pvalue = target[annotated_with_sample][0]
                    regulome_values.append(target_exp_value)
                    if target_exp_pvalue < 0.05 or (str(target_exp_pvalue) == str(np.nan) and abs(target_exp_value) >= 1.0):
                        reg2target_rels = according2._psrels4(_4reg_uid,target.uid())
                        rel = reg2target_rels[0]
                        sign = rel.effect_sign()
                        if sign != 0:
                            regulator_activation_score = regulator_activation_score + sign*target_exp_value
                            effect_target_counter += 1
                except KeyError:
                    continue

            if len(regulome_values) >= self.min_subnet_size and effect_target_counter > 0:
                regulator_activation_score = regulator_activation_score/math.sqrt(effect_target_counter)
            else:
                regulator_activation_score = np.nan
            return regulator_activation_score, regulome_values, effect_target_counter


    def __regulators4sample(self,sample:Sample,regulomes:dict,from_graph=ResnetGraph(),target_annotation_prefix=''):
        '''
        Input
        -----
        regulomes = {regulator_id:[PSObject]} made by ResnetGraph.regulome_dict(),\n
        where [PSObject] are targets of regulator_id

        Return
        ------
        [uids] - uids of regulators 
        regulators in self.Graph are annotated with SNEA activation score and p-value calculated from "sample" 
        '''
        my_graph = from_graph if from_graph else self.Graph
        sample_annotation_name = target_annotation_prefix+self.experiment.name4annotation(sample)

        sample_start = time.time()
        abs_sample_distribution = sample.data['value'].abs()
        abs_sample_distribution = abs_sample_distribution.dropna(how='all')
        sample_regulator_uids = set()
        
        for regulator_uid, targets in regulomes.items():
            regulator_activation_score,subnetwork_values,effect_target_counter = self.activity(sample,
                                            regulator_uid,targets,my_graph)

            if not subnetwork_values: continue
            sub_pd = df.from_dict({'value':subnetwork_values})
            abs_sub_pd = sub_pd['value'].abs()

            if regulator_activation_score > 0 or regulator_activation_score == np.nan:
                mw_stats,mv_pvalue = mannwhitneyu(x=abs_sub_pd, y=abs_sample_distribution, alternative = 'greater')
            else:
                mw_stats,mv_pvalue = mannwhitneyu(x=abs_sub_pd, y=abs_sample_distribution, alternative = 'less')
            """
            # alternative hypothesis: abs_sample_distribution is less than abs_sub_pd
            'two-sided': the distributions are not equal, i.e. F(u) ≠ G(u) for at least one u.
            'less': the distribution underlying x is stochastically less than the distribution underlying y, i.e. F(u) > G(u) for all u.
            'greater': the distribution underlying x is stochastically greater than the distribution underlying y, i.e. F(u) < G(u) for all u.
            One-sided alternative hypothesis tests median from one group can be greater or lesser than other group.
            """
        #    if from_graph._get_node(regulator_uid).name() == 'FOXM1':
        #        print('')
            if mv_pvalue <= self.subnet_pvalue_cutoff:
                new_prop_name = self.__sample_annotation(sample)
                prp_value = [float(regulator_activation_score),float(mv_pvalue),int(effect_target_counter)]
                nx.set_node_attributes(my_graph, {regulator_uid:{new_prop_name:prp_value}})
                sample_regulator_uids.add(regulator_uid)

        sample_time = execution_time(sample_start)
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

        print('SNEA execution time: %s' % execution_time(start_time))
        return
    

    def __regulators(self, sample_name='', from_graph=ResnetGraph()):
        '''
        Return
        ------
        if sample_name provided returns [PSobject] containing regulators for sample in self.__regulators__\n
        otherwise returns [PSobject] with regulators from all samples in self.__regulators__
        '''
        my_graph = from_graph if from_graph else self.Graph
        regulators_uids = set()
        if sample_name:
            regulators_uids.update(self.__regulators__[sample_name])
        else:
            [regulators_uids.update(reg_uids) for reg_uids in self.__regulators__.values()]
        return  my_graph._get_nodes(regulators_uids)


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


    def loadDPERNET(self):
        '''
        Return
        ------
        'drug-protein expression regulatory network'(DPERNET) ResnetGraph annotated with self.experiment values
        '''
        drug_target_network = self.dt_future.result()
        drug_target_expression_network = drug_target_network.subgraph_by_relprops(['Expression'])
        DPERNET4experiment = self.experiment.annotated_subnetwork(drug_target_expression_network,self.__my_sample_names__)
        self.Graph = self.Graph.compose(DPERNET4experiment)
        return DPERNET4experiment


    def find_drugs(self,DPERNET4experiment:ResnetGraph):
        '''
        Output
        ------
        self.__regulators__

        Returns
        -------
        [PSObject]
        '''
        print('Finding drugs inhibiting differential expression for each sample')
        all_drugs = set()
        drug2proteins_regulomes = DPERNET4experiment.regulome_dict(['SmallMol'],min_size=2)
        for sample in self.__my_samples():
            drugs4sample = list() 
            for drug_uid,targets in drug2proteins_regulomes.items():
                drug_activation_score = self.activity(sample,drug_uid,targets,DPERNET4experiment)[0]
                if drug_activation_score < 0:
                    drugs4sample.append(DPERNET4experiment._get_node(drug_uid))
            
            self.__sample2drugs__[sample.name()] = drugs4sample
            all_drugs.update(drugs4sample)
        # to make sure all ranked targets are in self.Graph which is required by Drugs4Targets::load_target_ranks() 
        return all_drugs


    def make_drugs_df(self):
        self.dt.__targets__ = set(self.__regulators())
        remove_targets_names = {'cytokine','inflammatory cytokine','anti-inflammatory cytokine','protein tyrosine kinase',
                          'mitogen-activated protein kinase kinase','mitogen-activated protein kinase','mixed lineage kinases',
                          'F-box domain'}
        
        DPERNET4experiment = self.loadDPERNET()
        self.find_drugs(DPERNET4experiment)

        remove_targets = {o for o in self.dt.__targets__ if o.name() in remove_targets_names}
        before_len = len(self.dt.__targets__)
        [self.dt.__targets__.remove(o) for o in remove_targets]
        print(f'{before_len-len(self.dt.__targets__)} unwanted regulators were removed before finding drugs')
        # SNEA is run using RNEF files only no dbids required
        mapped_targets_df = self.dt.load_df(list(self.dt.__targets__))
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
        if len(rank_columns) > 1:
            drugs_df['Average regulator score'] = drugs_df[rank_columns].mean(axis=1)
            drugs_df[RANK] = drugs_df['Average regulator score'] * drugs_df['Row count']/len(rank_columns)
            drugs_df.sort_values(by=[RANK],ascending=False,inplace=True)
            column_order = [PHARMAPENDIUM_ID,'Name',RANK,'Row count','Average regulator score']+rank_columns
        else:
            column_order = [PHARMAPENDIUM_ID,'Name']+rank_columns
        
        if self.dt.add_targets2drugs_ws:
            column_order += ['Directly inhibited targets','Directly activated targets','Indirectly inhibited targets','Indirectly activated targets']
        
        self.drugs_df = df(drugs_df[column_order])
        self.drugs_df.copy_format(drugs_df)
        self.drugs_df._name_ = 'Drugs'
        return self.drugs_df


    def report_path(self,extension:str):
        return self.experiment_dir+self.experiment.name()+' SNEA.'+extension


    def report(self,to_file=''):
        self.activity_matrix.fillna('', inplace=True)
        fout = self.experiment_dir +to_file+'.xlsx' if to_file else self.report_path('xlsx')
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
        self.activity_matrix.set_conditional_frmt(cond_format,
                'Average activation score','Average activation score')
        self.activity_matrix.make_header_vertical()
        self.activity_matrix.df2excel(writer,'activity')

        if not self.drugs_df.empty:
            self.drugs_df.make_header_vertical()
            self.drugs_df.df2excel(writer,'Drugs')

        writer.close()
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
    snea = SNEA(experiment_name,APIconfig,sample_ids=input,find_drugs=find_drugs)
    snea.set_dir('D:/Python/PBTA/')
    snea.expression_regulators()
    snea.activity_df()
    if find_drugs:
        snea.make_drugs_df()
    snea.report(experiment_name+'test')

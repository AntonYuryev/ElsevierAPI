from math import log2, isnan
import time
from .ResnetGraph import ResnetGraph,PSObject
from ..pandas.panda_tricks import df, pd, NaN
from  .PathwayStudioZeepAPI import DataModel
import scipy.stats as stats
from datetime import timedelta

HAS_PVALUE = 'hasPvalue'

class Sample(PSObject):

    def __init__(self,dic=dict()):
        super().__init__(dic)
        data_pd = pd.DataFrame(columns=['value','pvalue'])
        data_pd = data_pd.astype({'value': 'float64','pvalue': 'float64'}).dtypes
        self.data = df(data_pd)


    def copy(self, with_data=True):
        my_copy = Sample(self)
        if with_data:
            my_copy.data = df.copy_df(self.data)
        return my_copy


    def has_pvalue(self):
        try:
            return self[HAS_PVALUE]
        except KeyError:
            haspval = not self.data['pvalue'].dropna().empty
            self[HAS_PVALUE] = [haspval]
            return haspval


    def mapped_data(self):
        mapped_data = df(self.data.copy())
        mapped_data = mapped_data.dropna(subset=['URN'])
        return mapped_data

    
    def name(self,prefix=''):
        return prefix+':'+self['Name'][0] if prefix else super().name()


    def fisher_exact(self, urn2psobj:dict, difexp_pval_cutoff=0.05, using_prefix=''):
        """
        urn2psobj = {urn:PSObject}, where PSObject were annotated by self.nae4annotation()
        returns (oddsratio, ft_pvalue) of differential expression of PSObjects in urn2psobj\n
        compared with number of differential expressed genes in sample
        """
        #annotated_objs = self.annotate(urn2psobj,using_prefix) if psobjs_need_annotation else urn2psobj
        sample_annotation = self.name(using_prefix)

        de_input_count = len([x for x in urn2psobj.values() if x[sample_annotation][0][1] <= difexp_pval_cutoff])
        nonde_input_count = len(urn2psobj) - de_input_count

        try:
            de_sample_count = self['DEcount']
        except KeyError:
            de_sample_count = len(self.data.loc[(self.data['pvalue'] <= difexp_pval_cutoff)])
            self['DEcount'] = de_sample_count

        nonde_sample_count = len(self.data) - de_sample_count

        oddsratio, ft_pvalue = stats.fisher_exact([[de_input_count, nonde_input_count], [de_sample_count, nonde_sample_count]])
        return (oddsratio, ft_pvalue)


class Experiment(PSObject):
    pass
    
    def _add_sample(self, sample:Sample):
        self.samples[sample['Name'][0]] = sample


    def __init__(self, dic=dict(), identifier_name='',identifiers=[],samples=[]):
        super().__init__(dic)
        self.samples = dict() # {sample_name:Sample}
        if identifiers:
            urns = ['']*len(identifiers)
            self.identifiers = df.from_dict({identifier_name:identifiers,'URN':urns})
        else:
            identifiers_pd = pd.DataFrame(columns=['OriginalGeneID','URN'])
            identifiers_pd = identifiers_pd.astype({'OriginalGeneID': 'object','URN': 'object'}).dtypes
            self.identifiers = df(identifiers_pd)

        [self._add_sample(s) for s in samples]


    def size(self): return(len(self.identifiers))

    @staticmethod
    def execution_time(execution_start):
        return "{}".format(str(timedelta(seconds=time.time() - execution_start)))
    
    def urns(self):
        try:
            return self.identifiers['URN'].to_list()
        except KeyError:
            return list()


    def identifier_name(self):
        '''
        Return
        ------
        name of the first column in self.identifiers\n
        This name must be identical to property name of PSObject in ResnetGraph for mapping self.identifiers[URN] column
        '''
        try:
            return self.identifiers.columns[0]
        except KeyError:
            return str()


    def list_identifiers(self):
        identifier_col = self.identifier_name()
        return self.identifiers.loc[:,identifier_col].tolist()


    def original_identifiers(self, for_urns=set()):
        if for_urns:
            subset = self.identifiers[self.identifiers['URN'].isin(for_urns)]
            return subset.loc[:,'OriginalGeneID'].tolist()
        else:
            return self.identifiers.loc[:,'OriginalGeneID'].tolist()


    def get_sample_names(self,sample_ids):
        sample_name_list = list(self.samples.keys())
        return [str(k) for i,k in enumerate(sample_name_list) if i in sample_ids]


    def get_samples(self, sample_names=[],sample_ids=[]):
        s_names = sample_names
        if sample_ids:
            s_names += self.get_sample_names(sample_ids)

        if s_names:
            return [v for k,v in self.samples.items() if k in s_names]
        else:
            return [v for k,v in self.samples.items()]
        
    
    def mask_by_pval(self, max_pval=0.05):
        masked_experiment = Experiment(self)
        masked_experiment.identifiers = df.copy_df(self.identifiers)
        for sample in self.get_samples():
            masked_counter = 0
            sample_copy = sample.copy()
            if sample_copy.has_pvalue():
                for idx in sample_copy.data.index:
                    if sample_copy.data.at[idx,'pvalue'] > max_pval:
                        sample_copy.data.loc[idx,'value'] = NaN
                        masked_counter += 1
            masked_experiment.samples[sample['Name'][0]] = sample_copy
            print(f'{masked_counter} rows out of {len(sample.data)} total in sample {sample.name()} were masked because their p-value was greater than {max_pval}')
        return masked_experiment

#################  CONSTRUCTORS CONSTRUCTORS CONSTRUCTORS  ###################################
    @classmethod
    def download(cls,experiment_name:str, APIconfig:dict, only_log_ratio_samples=True):
        """
        returns Experiment object
        """
        print('\nDownloading "%s" experiment' % experiment_name,flush=True)
        start = time.time()
        ps_api = DataModel(APIconfig)
        experiment_zobj,experiment_identifiers,sample_values = ps_api._fetch_experiment_data(experiment_name,only_log_ratio_samples)
        
        psobj = PSObject.from_zeep_objects(experiment_zobj)
        experiment = cls(psobj)
        experiment.identifiers = df(experiment_identifiers)

        experiment.samples = dict()
        assert(len(experiment['SampleDefinitions'][0]['SampleDefinition']) == len(sample_values))
        assert(len(sample_values[0]) == len(experiment_identifiers))
        
        sample_id = 0
        samples_zobjects = experiment['SampleDefinitions'][0]['SampleDefinition']
        for s in samples_zobjects:
            sample_attributes = PSObject.from_zeep(s)
            sample = Sample(sample_attributes)
            sample['Type'] = ['Signal'] if sample['Type'] == 0 else ['Ratio']
            if sample['Subtype'] == 0: sample['Subtype'] = ['Unknown']
            elif sample['Subtype'] == 1: sample['Subtype'] = ['Bare']
            else: sample['Subtype'] = ['Log']

            results = ps_api.SOAPclient.service.GetExperimentValues(experiment['Id'][0],sample_id, 0, len(experiment_identifiers))
            values =  [row['Value'] for row in results]
            pvalues = [row['PValue'] for row in results]
            # <PValue>NaN</PValue> in experiments with no p-value columns

            sample.data = df.from_dict({'value':values,'pvalue':pvalues})
            experiment.samples[sample['Name'][0]] = sample
            sample_id += 1

        print('Experiment %s was downloaded in %s' % (experiment_name,Experiment.execution_time(start)))
        return experiment


    @classmethod
    def from_file(cls,*args,**kwargs):
        """
        Input
        -----
        filename = args[0], reads files with .xlsx, .tsv, .txt extensions\n
        data_type from ['LogFC', 'Intensity', 'LogIntensity', 'SignedFC']\n
        identifier_column - default 0\n
        phenotype_rows - default []\n
        header - default 0\n
        has_pvalue - default True\n
        sample_ids - defaults [], overrides "last_sample"
        last_sample - column with last sample, default None. If experiment has p-value columns "last_sample" must point to the last p-value column\n
        """
        experiment_file_name = str(args[0])
        identifier_column=int(kwargs.get('identifier_column',0))
        phenotype_rows=list(kwargs.get('phenotype_rows',[]))
        header=int(kwargs.get('header',0))
        data_type=str(kwargs.get('data_type','LogFC'))
        has_pvalue=bool(kwargs.get('has_pvalue',True))
        last_sample=kwargs.get('last_sample',None)
        sample_ids = list(kwargs.get('sample_ids',[]))
        if sample_ids:
            last_sample = max(sample_ids)
            if has_pvalue:
                last_sample += 1

        if phenotype_rows: 
            exp_df = df.read(experiment_file_name,skiprows=phenotype_rows,header=header)
            phenotypes_pd = df.read(experiment_file_name,header=header, nrows=len(phenotype_rows),index_col=False)
        else:
            exp_df = df.read(experiment_file_name,header=header,index_col=False)

        experiment_name = experiment_file_name[:experiment_file_name.rfind('.')]
        experiment_name = experiment_name[experiment_name.rfind('/')+1:]
        range_end = len(exp_df.columns) if type(last_sample) == type(None) else identifier_column+last_sample 

        samples = list()
        for i in range(identifier_column+1,range_end):
            column_name = exp_df.columns[i]
            sample_name = column_name
            sample = Sample({'Name':[sample_name]})
            if data_type in [ 'Intensity', 'LogIntensity']:
                sample['Type'] = ['Signal']
            else:
                sample['Type'] = ['Ratio']
            
            if data_type in ['LogFC', 'LogIntensity']:
                sample['Subtype'] = ['Log']
            else:
                sample['Subtype'] = ['Bare']

            if phenotype_rows:
                sample_phenotypes = list(phenotypes_pd[column_name])
                sample['phenotype'] = sample_phenotypes

            if has_pvalue:
                values = exp_df[exp_df.columns[i]].to_list()
                pvalues = exp_df[exp_df.columns[i+1]].to_list()
                sample.data = df.from_dict({'value':values,'pvalue':pvalues})
                sample['hasPValue'] = [0] # tp make consistent with database annotations
                samples.append(sample)
                i += 1
            else:
                values = exp_df[exp_df.columns[i]].to_list()
                sample.data = df.from_dict({'value':values})
                samples.append(sample)
        
        identifier_name = exp_df.columns[identifier_column]
        identifiers = exp_df[identifier_name].to_list()
        experiment = cls({'Name':[experiment_name]},identifier_name,identifiers,samples)
        return experiment


    def map(self, using_id2objs:dict):
        '''
        Input
        -----
        identifier2objs = {identifier:[PSObject]} - dictionary mapping self.identifier_name() values to corresponding [PSObject]
        '''
        identifier_name = self.identifier_name()
        identifiers = self.list_identifiers()

        unmapped_identifiers = set(identifiers).difference(set(using_id2objs.keys()))
        unmapped_str = '\n'.join(unmapped_identifiers)
        print(f'Following {len(unmapped_identifiers)} identifiers out of {len(identifiers)} in {self.name()} experiment were not mapped:\n{unmapped_str}')

        identifiers_columns = self.identifiers.columns.to_list()
        #mapped_entities_df = df(columns = identifiers_columns)
        identifiers_lst = list()
        urn_lst = list()
        for i in identifiers:
            try:
                psobjs = using_id2objs[i]
                urns = [o.urn() for o in psobjs]
                identifiers_lst += [i]*len(urns)
                urn_lst += urns
            except KeyError:
                continue     
        assert(len(identifiers_lst) == len(urn_lst))
                
        mapped_entities_df = df.from_dict({identifier_name:identifiers_lst,'URN':urn_lst})        
        samples_pd = self.samples2pd()
        mapped_data = mapped_entities_df.merge(samples_pd,'outer', on=identifier_name)
        mapped_data = mapped_data.drop('URN_y', axis=1)
        mapped_data.rename(columns={'URN_x':'URN'},inplace=True)

        mapped_exp = Experiment(self)
        mapped_exp.identifiers = df(mapped_data[identifiers_columns])

        column_counter = len(identifiers_columns)
        for idx, sample in enumerate(self.get_samples()):
            mapped_sample = Sample(sample)
            values = mapped_data.iloc[:,idx+column_counter].to_list()
            if sample.has_pvalue():
                pvalues = mapped_data.iloc[:,idx+column_counter+1].to_list()
                mapped_sample.data = df.from_dict({'value':values,'pvalue':pvalues})
            else:
                mapped_sample.data = df.from_dict({'value':values})
            mapped_exp._add_sample(mapped_sample)

        print('%d out of %d entities were mapped to URNs in %s experiment'%
                (len(mapped_exp.identifiers),len(self.identifiers),self.name()))

        return mapped_exp


    def map2(self,network:ResnetGraph,using_original_id2:str):
        original_id2psobjs, uid2original_id = network.get_prop2obj_dic(using_original_id2)
        if len(original_id2psobjs) < 2:
            print(f'Mapping using {using_original_id2} has failed!!!')
        return self.map(original_id2psobjs)


    def ttest(self,control_phenotype:str,case_phenotype:str):
        """
        samples must have 'phenotype' annotation\n
        returns logFC sample and adds it to self
        """
        control_cols = pd.concat([s.data['value'] for s in self.samples.values() if control_phenotype in s['phenotype']], axis=1)
        case_cols = pd.concat([s.data['value'] for s in self.samples.values() if case_phenotype in s['phenotype']], axis=1)

        sample_ave = case_cols.mean(axis=1)
        control_ave = control_cols.mean(axis=1)
        ratio = sample_ave/control_ave
        log2FCs  = df(ratio).apply(lambda x: log2(x),axis=1)
        ttest,ttest_pvalue = stats.ttest_ind(control_cols,case_cols,axis=1)
        identifier_name = self.identifier_name()

        sample_name = case_phenotype+'/'+control_phenotype
        logfc_sample = Sample({'Name':[sample_name],HAS_PVALUE:[True]})
        logfc_sample.data['value'] = log2FCs
        logfc_sample.data['pvalue'] = ttest_pvalue

        exp_attributes = dict(self)
        ttest_exp_name = self.name()+'_ttest'
        exp_attributes.update({'Type':['Ratio'], 'Subtype':['Log'], 'Name':[ttest_exp_name]})
        ttest_exp = Experiment(exp_attributes,identifier_name, self.list_identifiers(), [logfc_sample])
        return ttest_exp


############################ PSObject ANNOTATION PSObject ANNOTATION PSObject ANNOTATION ################################
    def prefix4annotation(self):
        '''
        Return
        ------
        Experiment name
        '''
        return self.name()


    def name4annotation(self, s:Sample):
        '''
        Return
        ------
        Experiment name : sample name
        '''
        return s.name(self.prefix4annotation())


    def samples2pd(self,sample_names=[],sample_ids=[]):
        samples = self.get_samples(sample_names, sample_ids)
        dfs2concatenate = [self.identifiers] + [s.data for s in samples]
        return pd.concat(dfs2concatenate,axis=1)


    def __sample2pd(self, sample:Sample, deduplicate=True):
        #sample = self.samples[sample_name]
        sample_pd = pd.concat([self.identifiers, sample.data],axis=1)

        if deduplicate:
            if sample.has_pvalue():
                sample_pd.sort_values(['URN','pvalue'], ascending=True, inplace=True)
                sample_pd.drop_duplicates(subset=['URN'], inplace=True, ignore_index=True)
            else:
                sample_pd.sort_values(['URN','value'], ascending=False, inplace=True)
                sample_pd.drop_duplicates(subset=['URN'], inplace=True, ignore_index=True)

        return sample_pd


    def __make_dict4annotation(self, sample:Sample):
        sample_pd = self.__sample2pd(sample)
        urn2value = dict()
        for idx in sample_pd.index:
            urn = sample_pd.loc[idx,'URN']
            value = sample_pd.loc[idx,'value']
            if not isnan(value):
                pvalue = sample_pd.loc[idx,'pvalue']
                urn2value[urn] = [(value,pvalue)]

        return self.name4annotation(sample), urn2value


    def annotated_subnetwork(self,graph:ResnetGraph,with_values_from_sample_names:list=[],with_values_from_sample_ids:list=[]):
        '''
        Returns
        -------
        ResnetGraph with nodes annotated with "self.name4annotation(sample)" properties for all samples
        Nodes that do not have value in all samples and are not regulators in input graph are removed
        '''
        graph_copy = graph.copy()
        urns_with_values = set()
        for sample in self.get_samples(with_values_from_sample_names,with_values_from_sample_ids):
            node_annotation_name, urn2expvalue = self.__make_dict4annotation(sample)
            graph_copy.set_node_annotation(urn2expvalue,node_annotation_name)
            urns_with_values.update(urn2expvalue.keys())

        subnetwork_name = graph.name+f' for "{self.name()}"'
        return graph_copy.regulatory_network_urn(urns_with_values,subnetwork_name)


    def annotate_objs(self,urn2obj:dict,sample_names=[],sample_ids=[]):
        for sample in self.get_samples(sample_names,sample_ids):
            node_annotation_name, urn2value = self.__make_dict4annotation(sample)
            for u,values in urn2value.items():
                try:
                    urn2obj[u][node_annotation_name] = values
                except KeyError:
                    continue
                
    




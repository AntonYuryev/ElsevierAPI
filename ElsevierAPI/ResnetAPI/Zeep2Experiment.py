from .ResnetGraph import ResnetGraph,PSObject
from ..pandas.panda_tricks import df, np, pd
from  .PathwayStudioZeepAPI import DataModel
import scipy.stats as stats
from math import log2


class Sample(PSObject):
    # self['Type'] = ['Signal' or 'Ratio']
    # self['Subtype'] = ['Unknown' or 'Bare' or 'Log']
    # self['phenotype'] = [phenotypes]
    #def __init__(self,columns=['value','pvalue'], dic=dict()):


    def __init__(self,dic=dict()):
        super().__init__(dic)
        self.data = df(columns=['value','pvalue'])
        self.data.astype({'value': 'float64'}).dtypes
        self.data.astype({'pvalue': 'float64'}).dtypes


    def copy(self, with_data=True):
        my_copy = Sample(self)
        if with_data:
            my_copy.data = df.make_df(self.data)
        return my_copy


    def has_pvalue(self):
        try:
            return self['hasPvalue']
        except KeyError:
            haspval = not self.data['pvalue'].dropna().empty
            self['hasPvalue'] = haspval
            return haspval


    def mapped_data(self):
        mapped_data = df(self.data.copy())
        mapped_data = mapped_data.dropna(subset=['URN'])
        return mapped_data

    
    def name(self,prefix:str):
        return prefix+':'+self['Name'][0]


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
            self.identifiers = df(columns=['OriginalGeneID','URN'])
            self.identifiers.astype({'OriginalGeneID': 'object'}).dtypes
            self.identifiers.astype({'URN': 'object'}).dtypes

        [self._add_sample(s) for s in samples]

    
    def urns(self):
        try:
            return self.identifiers['URN'].to_list()
        except KeyError:
            return list()


    def identifier_name(self):
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
        return [k for i,k in enumerate(sample_name_list) if i in sample_ids]


    def get_samples(self, sample_names=[],sample_ids=[]):
        s_names = sample_names
        if sample_ids:
            s_names += self.get_sample_names(sample_ids)

        if s_names:
            return [v for k,v in self.samples.items() if k in s_names]
        else:
            return [v for k,v in self.samples.items()]


#################    CONSTRUCTORS CONSTRUCTORS CONSTRUCTORS  ###################################
    @classmethod
    def download(cls,experiment_name:str, APIconfig:dict, only_log_ratio_samples=True):
        """
        returns Experiment object
        """
        print('Downloading "%s" experiment' % experiment_name)
        ps_api = DataModel(APIconfig['ResnetURL'], APIconfig['PSuserName'], APIconfig['PSpassword'])
        experiment_zobj,experiment_identifiers,sample_values = ps_api._fetch_experiment_data(experiment_name,only_log_ratio_samples)
        
        psobj = PSObject.from_zeep_objects(experiment_zobj)
        experiment = Experiment(psobj)
        experiment.identifiers = experiment_identifiers

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

            sample.data  = df.from_dict({'value':values,'pvalue':pvalues})
            experiment.samples[sample['Name'][0]] = sample
            sample_id += 1

        return experiment


    @classmethod
    def from_file(cls,experiment_file_name:str,identifier_column=0,phenotype_rows=list(),
                    header=0,data_type='LogFC',has_pvalue=True):
        """
        reads files with .xlsx, .tsv, .txt extensions
        data_type in ['LogFC', 'Intensity', 'LogIntensity', 'SignedFC']
        """

        if phenotype_rows: 
            exp_df = df.read(experiment_file_name,skiprows=phenotype_rows,header=header)
            phenotypes_pd = df.read(experiment_file_name,header=header, nrows=len(phenotype_rows))
        else:
            exp_df = df.read(experiment_file_name,header=header)

        experiment_name = experiment_file_name[:experiment_file_name.rfind('.')]
        experiment_name = experiment_name[experiment_name.rfind('/')+1:]

        experiment = Experiment({'Name':[experiment_name]})
        experiment.identifiers['OriginalGeneID'] = exp_df.iloc[:,identifier_column]

        for i in range(identifier_column+1,len(exp_df.columns)):
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
                sample.data['value'] = exp_df[exp_df.columns[i]]
                sample.data['pvalue'] = exp_df[exp_df.columns[i+1]]
                sample['hasPValue'] = [0]
                experiment.__add_sample(sample)
                i += 1
            else:
                sample.data['value'] = exp_df[exp_df.columns[i]]
                experiment._add_sample(sample)
        
        return experiment


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
        logfc_sample = Sample({'Name':[sample_name],'hasPvalue':[True]})
        logfc_sample.data['value'] = log2FCs
        logfc_sample.data['pvalue'] = ttest_pvalue

        exp_attributes = dict(self)
        ttest_exp_name = self.name()+'_ttest'
        exp_attributes.update({'Type':['Ratio'], 'Subtype':['Log'], 'Name':[ttest_exp_name]})
        ttest_exp = Experiment(exp_attributes,identifier_name, self.list_identifiers(), [logfc_sample])
        return ttest_exp


############################ PSObject ANNOTATION PSObject ANNOTATION PSObject ANNOTATION ################################
    def prefix4annotation(self):
        return self['Name'][0]


    def name4annotation(self, s:Sample):
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
            pvalue = sample_pd.loc[idx,'pvalue']
            urn2value[urn] = [(value,pvalue)]

        return self.name4annotation(sample), urn2value


    def annotate(self, graph:ResnetGraph, sample_names=[],sample_ids=[]):
        for sample in self.get_samples(sample_names,sample_ids):
            node_annotation_name, urn2value = self.__make_dict4annotation(sample)
            if graph:
                graph.set_annotation(urn2value,node_annotation_name)
                graph.urn2obj = {o['URN'][0]:o for i,o in graph.nodes(data=True)}
            
            
    def annotate_objs(self,urn2obj:dict,sample_names=[],sample_ids=[]):
        for sample in self.get_samples(sample_names,sample_ids):
            node_annotation_name, urn2value = self.__make_dict4annotation(sample)
            for u,values in urn2value.items():
                try:
                    urn2obj[u][node_annotation_name] = values
                except KeyError:
                    continue
                
    




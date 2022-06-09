from .ResnetGraph import ResnetGraph,PSObject,PSRelation
from ..pandas.panda_tricks import df,np
from  .PathwayStudioZeepAPI import DataModel

class Sample(PSObject):
    name = str()
    type = str()  # 'Signal' or 'Ratio'
    subtype = str() # 'Unknown', 'Bare', 'Log'
    phenotype = str()
    data = df()  # df(columns=['URN','value','pvalue'])

    def has_pvalue(self): return type(self['hasPValue'][0]) != type(None)
    def urns(self): return self.data['URN'].to_list()


class Experiment(PSObject):
    pass
    samples = dict()
        
    def __init__(self,experiment_name:str, APIconfig:dict, only_log_ratio_samples=True):
        """
        returns Experiment object
        """
        print('Downloading "%s" experiment' % experiment_name)
        self.samples = dict()
        ps_api = DataModel(APIconfig['ResnetURL'], APIconfig['PSuserName'], APIconfig['PSpassword'])
        experiment,experiment_urns,sample_values = ps_api._fetch_experiment_data(experiment_name,only_log_ratio_samples)

        super().__init__(PSObject.from_zeep_objects(experiment))
        self.samples = dict()
        assert(len(self['SampleDefinitions'][0]['SampleDefinition']) == len(sample_values))
        assert(len(sample_values[0]) == len(experiment_urns))

        sample_id = 0
        samples_zobj = self['SampleDefinitions'][0]['SampleDefinition']
        for s in samples_zobj:
            sample = Sample.from_zeep(s)
            sample['Type'] = ['Signal'] if sample['Type'] == 0 else ['Ratio']
            if sample['Subtype'] == 0: sample['Subtype'] = ['Unknown']
            elif sample['Subtype'] == 1: sample['Subtype'] = ['Bare']
            else: sample['Subtype'] = ['Log']

            results = ps_api.SOAPclient.service.GetExperimentValues(self['Id'][0],sample_id, 0, len(experiment_urns))
            values =  [row['Value'] for row in results]
            if sample.has_pvalue():
                pvalues = [row['PValue'] for row in results]
                sample.data = df.from_dict({'URN':experiment_urns,'value':values,'pvalue':pvalues})
                sample.data = sample.data.sort_values(['URN','pvalue'],ascending=True,inplace=True)
                sample.data.drop_duplicates(subset=['URN'],inplace=True, ignore_index=True)
            else:
                sample.data = df.from_dict({'URN':experiment_urns,'value':values})
                sample.data.sort_values(['URN','value'],ascending=False,inplace=True)
                sample.data.drop_duplicates(subset=['URN'],inplace=True, ignore_index=True)

            self.samples[sample['Name'][0]] = sample
            sample_id += 1

    def urns(self):
        smpls = list(self.samples.values())
        return smpls[0].urns()

    def start_relevant_networks(self,experiment_id:int, sample_id=0, 
                            expand_by_reltypes=[], expand2types=[], find_regulators=True):

        grn_prams = self.create_GenRelNetsParameters(experiment_id,sample_id,expand_by_reltypes,expand2types,find_regulators)
        result = self.SOAPclient.service.StartGenerateRelevantNetworks(grn_prams)
        """
        result['State']:
                0 - waiting
                1 - running
                2 - finished
                3 - killed
        """
        
        return result

    def get_samples(self,sample_ids=[]):
        all_samples = list(self.samples.values())
        if sample_ids:
            selected_samples = [all_samples[s] for s in sample_ids] 
            return selected_samples
        else: 
            return all_samples


    def annotate_graph(self, g:ResnetGraph,sample_ids=[]):
        experiment_name = self['Name'][0]
        for sample in self.get_samples(sample_ids):
            sample_name = sample['Name'][0]
            urn2value = dict()
            for idx in sample.data.index:
                urn = sample.data.loc[idx,'URN']
                value = sample.data.loc[idx,'value']
                if sample.has_pvalue():
                    pvalue = sample.data.loc[idx,'pvalue']
                else:
                    pvalue = np.nan

                urn2value[urn] = [(value,pvalue)]
            g.set_annotation(urn2value,experiment_name+':'+sample_name)

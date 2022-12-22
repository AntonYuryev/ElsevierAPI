import time
import glob
import os
import scipy.stats as stats
import xml.etree.ElementTree as et
from .rnef2sbgn import make_file_name,to_sbgn_file
from .FolderContent import FolderContent,PSPathway
from .Zeep2Experiment import Experiment
from ..pandas.panda_tricks import df,ExcelWriter
from ..ETM_API.references import RELATION_PROPS

MEASURED_COUNT = '# measured entities'
MEASURED_ENTITIES = 'measured entities'

def get_file_listing(dirname:str):
    listing = glob.glob(dirname+'/*.rnef')
    return list(map(os.path.basename, listing))


class GSEA(FolderContent):
    def __init__(self, APIconfig):
        super().__init__(APIconfig)
        self.mv_pvalue_cutoff = 0.05
        self.diffexp_pvalue_cuoff= 0.05
        self.id2pathway = dict() # self.id2pathway = {id:PSPathway}
        self.data_dir = 'ElsevierAPI/ResnetAPI/__pscache__/'
        self.report_dir = ''
        wsdl_url = str(self.SOAPclient.wsdl.location)
        self.base_url = wsdl_url[:wsdl_url.find('/services/')]+'/app/sd?urn='
        

    def __add_pathway(self,p:PSPathway):
        try:
            p_id = p.id()
            self.id2pathway[p_id] = p
            print('Pathway %s already exist' % p.name())
        except KeyError:
            p_id = len(self.id2pathway)
            p['Id'] = [p_id] # makes fake IDs
            self.id2pathway[p_id] = p
            self.Graph.urn2obj.update(p.graph.urn2obj)


    @staticmethod
    def __find_valid_resnets(rnef_file:str, must_have_urns:list, min_overlap:int):
        bottom = min_overlap-1
        try:
            with open(rnef_file, 'r', encoding='utf-8') as f:
                resnet = str()
                line = ' '
                valid_resnets = list()
                resnet_count = 0
                while True:
                    while line and line[:7] != '<resnet':
                        line = f.readline().strip()

                    if not line:
                        break

                    resnet += line
                    node_count = 0
                    while line[:8] != '</resnet':
                        line = f.readline().strip()
                        if line[:6] =='<node ':
                            if node_count < min_overlap:
                                urn_start = line.find('urn=',20)
                                assert(urn_start > 0)
                                urn_end = line.rfind('\"')
                                urn = line[urn_start+5:urn_end]
                                if urn in must_have_urns:
                                    node_count += 1
                                
                        resnet  += line

                    if node_count > bottom:
                        valid_resnets.append(resnet)

                    resnet_count += 1
                    resnet = ''
            
                print ('Found %d out of %d pathways with > %d valid entities in "%s"' % 
                        (len(valid_resnets), resnet_count, bottom, rnef_file[rnef_file.rfind('\\')+1:]))
                return valid_resnets
        except KeyError:
            print('Cannot open %s file' % rnef_file)
            return []


    def read_rnef(self, rnef_file:str, must_have_urns:list, min_overlap=3, add_annotation=dict()):
        valid_resnets = self.__find_valid_resnets(rnef_file,must_have_urns,min_overlap)
 
        for resnet in valid_resnets:
            resnet_et = et.fromstring(resnet)
            if str(resnet_et.get('type')).lower() in ['pathway','group']:
                pathway = PSPathway.from_resnet(resnet_et,add_annotation)
            if type(pathway != type(None)):
                self.__add_pathway(pathway)


    def __load_cache(self, must_have_urns:list, folder_name='',min_overlap=3):
        dump_path = self.dump_path(folder_name)
        listing = glob.glob(dump_path+'**/*.rnef',recursive = True)
        if not listing:
            raise FileNotFoundError
        file_count = 0
        for fname in listing:
            dirpath = os.path.dirname(fname)
            file_folder = os.path.basename(dirpath)
            file_name = os.path.basename(fname)
            add_annotation = {'Folder':[file_folder],'File':[file_name]}
            self.read_rnef(fname,must_have_urns,add_annotation=add_annotation,min_overlap=min_overlap)
            file_count += 1
            print('%d out of %d files loaded' % (file_count,len(listing)))


    def load_pathways(self, pathway_folders:list, must_have_urns:list,min_overlap=3):
        missing_in_cache = list()
        for folder_name in pathway_folders:
            try:
                self.__load_cache(must_have_urns,folder_name,min_overlap)
            except FileNotFoundError:
                print('%s folder is not in cache. It will be downloaded from database\nDownload may take long time'
                        % folder_name)
                missing_in_cache.append(folder_name)

        self.relProps = list(RELATION_PROPS) + ['URN']
        if missing_in_cache:
            for folder_name in missing_in_cache:
                self.folder2rnef(folder_name)
                self.__load_cache(must_have_urns,folder_name)


    def gsea(self, experiment:Experiment, sample_names=[], sample_ids=[], calculate_activity=False):
        start_time = time.time()
        print('Performing GSEA and Fisher Exact test')

        samples = experiment.get_samples(sample_names,sample_ids)
        annotated_pathway_counter = set()
        experiment.annotate_objs(self.Graph.urn2obj,sample_names=[],sample_ids=[0])

        for sample in samples:
            sample_annotation = experiment.name4annotation(sample)
            abs_sample_distribution = sample.data['value'].abs()

            for ps_pathway in self.id2pathway.values():
                annotated_pathway_obj = {u:o for u,o in self.Graph.urn2obj.items() if u in ps_pathway.graph.urn2obj and sample_annotation in o}
                # self.Graph.urn2obj can have duplicate 'Id' because it collects objects from multiple pathways

                pathway_values = [v[sample_annotation][0][0] for v in annotated_pathway_obj.values()]
                ps_pathway[MEASURED_COUNT]= [len(pathway_values)]
                measured_entities_urns = {n.urn() for n in annotated_pathway_obj.values()}
                ps_pathway[MEASURED_ENTITIES] = experiment.original_identifiers(measured_entities_urns)
                
                sub_pd = df.from_dict({'value':pathway_values})
                abs_sub_pd = sub_pd['value'].abs()
                mw_results = stats.mannwhitneyu(x=abs_sample_distribution, y=abs_sub_pd, alternative = 'less') 
                mv_pvalue = mw_results[1]

                ft_result = sample.fisher_exact(annotated_pathway_obj,using_prefix=experiment.prefix4annotation())
                if ft_result[1] < 0.05:
                    ps_pathway[sample_annotation+':FisherExact'] = [ft_result]
                                
                if mv_pvalue <= self.mv_pvalue_cutoff:
                    logfc_sum = sum([n[sample_annotation][0][0] for n in annotated_pathway_obj.values()])
                    enity_count = ps_pathway.graph.number_of_nodes(experiment.objtype())
                    ps_pathway[sample_annotation+':GSEA'] = [(logfc_sum/enity_count, mv_pvalue)]
                    
                    if calculate_activity:
                        pathway_activity = 0.0
                        for urn, node in annotated_pathway_obj.items():
                            node_id = ps_pathway.urn2obj[urn]['Id'][0]
                            node_downstrem_rels = ps_pathway.graph.downstream_relations(node_id)
                            node_impact = 0
                            node_impact += sum([r.effect_sign() for r in node_downstrem_rels])
                            pathway_activity += node_impact*node[sample_annotation][0][0]
                    
                        ps_pathway[sample_annotation+':Activity'] = [pathway_activity]

                    annotated_pathway_counter.add(ps_pathway)

        print('GSEA was done in %s' % self.execution_time(start_time))
        print('%d pathways were annotated with %d samples from %s' % 
                (len(annotated_pathway_counter), len(experiment.samples),experiment.name() ))
        

    def make_report(self,experiment:Experiment)->tuple([df,list]):
        report_pd = df(columns=['Pathway','Folder','File',MEASURED_COUNT, MEASURED_ENTITIES])
        significant_pathway_ids = list()
        row = 0
        for pathway_obj in self.id2pathway.values():
            is_significant = False
            for sample in experiment.samples.values():
                node_annotation = experiment.name4annotation(sample)
                try:
                    ft_result = pathway_obj[node_annotation+':FisherExact'][0]
                    report_pd.loc[row,node_annotation+':Odds Ratio'] = ft_result[0]
                    report_pd.loc[row,node_annotation+':FisherExact pvalue'] = ft_result[1]
                    is_significant = True
                except KeyError:
                    pass

                try:
                    gsea_result = pathway_obj[node_annotation+':GSEA'][0]
                    report_pd.loc[row,node_annotation+':logFC mean'] = gsea_result[0]
                    pvalue_column = node_annotation+':GSEA pvalue'
                    report_pd.loc[row,pvalue_column] = gsea_result[1]
                    is_significant = True
                except KeyError:
                    pass

                try:
                    activity = pathway_obj[node_annotation+':Activity'][0]
                    report_pd.loc[row,node_annotation+':Activity'] = activity
                except KeyError:
                    pass

                if is_significant:
                    report_pd.loc[row,'Pathway'] = '=HYPERLINK("'+self.base_url+pathway_obj.urn()+'",\"'+pathway_obj.name()+'")'
                    report_pd.loc[row,'Folder'] = pathway_obj['Folder'][0]
                    report_pd.loc[row,'File'] = pathway_obj['File'][0]
                    report_pd.loc[row,MEASURED_COUNT] = pathway_obj[MEASURED_COUNT][0]
                    report_pd.loc[row,MEASURED_ENTITIES] = '; '.join(pathway_obj[MEASURED_ENTITIES])
                    significant_pathway_ids.append(pathway_obj.id())
                    row+=1

        pvalue_columns = [c for c in report_pd.columns if 'GSEA pvalue' in c]
        if pvalue_columns:
            report_pd = report_pd.sort_values(pvalue_columns,ascending=True)
            return report_pd,significant_pathway_ids
        else:
            print('No pathways with pvalue < %g were found' % self.mv_pvalue_cutoff)
            return df(), list()


    def report(self, folders_with_target_pathways:list,experiment:Experiment, format='SBGN'):
        report_pd, significant_pathway_ids = self.make_report(experiment)

        folders_name = ','.join(folders_with_target_pathways)

        fout = experiment.name()+' GSEA of '+ make_file_name(folders_name)
        report_fout = self.report_dir+fout

        clean_df = df(report_pd).clean()
        writer = ExcelWriter(report_fout+'.xlsx', engine='xlsxwriter')
        clean_df.make_header_vertical()
        clean_df.df2excel(writer,'GSEA')
        writer.save()
        print('Ranked pathways are in %s file' % fout)

        file_ext = '.sbgn' if format =='SBGN' else '.rnef'
        for i in significant_pathway_ids:
            pathway = self.id2pathway[i]
            pathway_xml = pathway['resnet']
            foutname = make_file_name(pathway.name())
            fout_name = foutname + file_ext
            if format == 'SBGN':
                to_sbgn_file('<batch>'+pathway_xml+'</batch>', foutname,self.report_dir)
            else:
                with open(fout_name, 'w', encoding='utf-8') as f:
                    f.write(pathway_xml)


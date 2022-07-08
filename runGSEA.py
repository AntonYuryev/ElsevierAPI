from ElsevierAPI.ResnetAPI.GSEA import GSEA
from ElsevierAPI.ResnetAPI.Zeep2Experiment import Experiment
from ElsevierAPI import load_api_config, execution_time, time

def make_ps_name(exp:Experiment):
    def convert2ps(lipid_name:str):
        head_group = lipid_name[:lipid_name.find('_')]
        chains = lipid_name[lipid_name.rfind('_')+1:]
        return head_group+'('+chains+')'

    exp.identifiers.insert(0,'Name,Alias', exp.identifiers.iloc[:,0].apply(lambda x: convert2ps(x)))


if __name__ == "__main__":
    start = time.time()
    exp_fname = 'D:/Python/MDACC/Patwhay-studio_demo_data_Iqbal-MDACC-BCB.csv'
    my_experimet = Experiment.from_file(exp_fname,header=0,phenotype_rows=[1],data_type='Intensity',has_pvalue=False)

    '''
    0 - header row
    1 - phenotype rows must start immediately after header rows
    '''
    # calculate differential expression:
    gsea_exp = my_experimet.ttest(control_phenotype='WT',case_phenotype='KO')
    # adds 'Name' columns with canonical lipid identifiers:
    make_ps_name(gsea_exp)
    gsea_exp['ObjTypeName'] = ['SmallMol']

    APIcofig = 'D:/Python/ENTELLECT_API/ElsevierAPI/APIconfigMDACC.json'
    gsea = GSEA(load_api_config(APIcofig))
    mapped_exp = gsea.map_experiment(gsea_exp)

    folders_with_target_pathways = ['TG triglyceride metabolism']
    min_overlap = 3
    #['Globosides metabolism','Sphingomyelins metabolism']
    #['Sphingolipids metabolism']
    #['PE-PC-PS cycle']
    #['Plasmalogen biosynthesis']
    #['PE metabolism']
    #['CL cardiolipin metabolism']
    #['PAF biosynthesis']#['PS metabolism (FA->DG->PC)'] 
    #['TG triglyceride metabolism']
    # #['Test'] # #['Neutral and Phospholipids metabolism'] # 
    gsea.entProps=['Name'] # for download speed
    gsea.load_pathways(folders_with_target_pathways,mapped_exp.urns(),min_overlap=min_overlap)
    
    gsea.gsea(mapped_exp)
    gsea.report_dir = 'D:/Python/MDACC/' 
    #change out_dir to directory on your computer where you want to see the report Excel file and pathway hits in SBGN files
    fout_name = mapped_exp.name()+' GSEA of '+','.join(folders_with_target_pathways)+' results'

    gsea.base_url = 'https://mdacc.pathwaystudio.com/app/sd?urn='
    report_df = gsea.report(fout_name,mapped_exp)
    print('GSEA analysis was done in %s' % execution_time(start))

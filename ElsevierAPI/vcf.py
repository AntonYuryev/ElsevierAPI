import myvariant, time, os, ast
from .utils import execution_time,fname,dirList
from .pandas.panda_tricks import df,pd
from statistics import mean

vcf_format_definitions = {
    "GT": "Genotype. Encodes the alleles found in the sample (e.g., 0/1 for heterozygous, 1/1 for homozygous alternate).",
    "AD": "Allele Depth. A list of integers representing the number of reads supporting each allele (e.g., 'Ref_Count,Alt_Count').",
    "AF": "Allele Fraction. The fraction of reads in this sample that support the alternate allele (calculated as Alt_Count / Total_Depth).",
    "DP": "Depth. The total number of filtered, high-quality reads covering this position used for calling.",
    "F1R2": "Forward 1st / Reverse 2nd. Count of reads supporting the ref/alt alleles in the F1R2 orientation pairs. Used to detect sequencing artifacts (like OxoG).",
    "F2R1": "Forward 2nd / Reverse 1st. Count of reads supporting the ref/alt alleles in the F2R1 orientation pairs. The counterpart to F1R2.",
    "GQ": "Genotype Quality. A Phred-scaled score representing the confidence that the Genotype (GT) call is correct.",
    "PL": "Phred-scaled Likelihoods. The normalized likelihoods for each potential genotype (e.g., 0/0, 0/1, 1/1). The lowest value is always 0 (the most likely genotype).",
    "GP": "Genotype Posterior Probabilities. The probability of each possible genotype, ranging from 0 to 1 (often used in imputation).",
    "PRI": "Prior Probability. The prior probability of the variant existing, often used in Bayesian somatic calling models.",
    "SB": "Strand Bias. Counts of reads supporting alleles on specific strands (usually: Ref_Forward, Ref_Reverse, Alt_Forward, Alt_Reverse).",
    "MB": "Mismatch Bias. Statistics related to mismatches within the reads mapping to this site, used to filter false positives caused by poor mapping."
}

def read_vcf(vcf_path: str) -> tuple[str, pd.DataFrame]:
  """
  args:
    vcf_path: Path to the VCF file.
  output:
    tuple (VCF header string, pandas DataFrame with the variant data)
  """
  vcf_header = []
  data_start_line = 0

  with open(vcf_path, 'r') as f:
    for i, line in enumerate(f):
      if line.startswith('##'):
        vcf_header.append(line)
      elif line.startswith('#'):
        data_start_line = i + 1 # Column header line just before the data
        break
      else:
        break # Handle cases where the VCF file might be malformed or only has data
  vcf_header_str = "".join(vcf_header)

  # 2. Second pass: Efficiently read the data into a pandas DataFrame
  vcfname = fname(vcf_path)
  try:
    vcf_pd = pd.read_csv(
      vcf_path,
      sep='\t',
      skiprows=data_start_line, # Skip all header lines in one go
      comment=None,             # The comment parameter is handled by skiprows
      header=0                  # Use the first row after skipping as the column names
    )
    # VCF column headers are prepended with '#CHROM'. pd.read_csv treats the first
    # column header (the one at skiprows) as the column names. We need to rename 
    # the first column from '#CHROM' to 'CHROM'.
    if vcf_pd.columns[0].startswith('#'):
      vcf_pd.rename(columns={vcf_pd.columns[0]: vcf_pd.columns[0][1:]}, inplace=True)

  except pd.errors.EmptyDataError: # Handle case where the file only contains header lines
    print(f"File {vcfname} contains no variant data.")
    vcf_pd = pd.DataFrame() # Return an empty DataFrame

  print(f'Read {vcfname} vcf file with {len(df)} variants and columns: {list(df.columns)}')
  return vcf_header_str, df.from_pd(vcf_pd,dfname=vcfname)


def HGVS2vcf(vcf_header:str,vcf_df:df,sample_col:str):
  '''
  output:
    modified VCF header
    vcf_df - dataframe with VCF data
    genotype_df -  dataframe with columns: HGVS, two alleles genotype infered from VCF data
  '''
  genotype_df_rows = []
  for i in vcf_df.index:
    sample_fields = vcf_df.at[i,sample_col].split(':')
    gt_str = '1/1' if sample_fields[0] =='1' else sample_fields[0].replace('|','/')
    GTcalls = list(map(int,gt_str.split('/')))
    ref_allele = vcf_df.at[i,'REF']
    alleles = [ref_allele] + vcf_df.at[i,'ALT'].split(',')
    genotype = [alleles[GTcalls[0]],alleles[GTcalls[1]]]
    hgvs = []
    for alt_allele in alleles[1:]:
      if alt_allele in genotype:
        hgv = f"{vcf_df.at[i,'CHROM']}:g.{vcf_df.at[i,'POS']}{ref_allele}>{alt_allele}"
        genotype_df_rows.append([hgv,''.join(genotype)])
        hgvs.append(hgv)

    frmt = vcf_df.at[i,'FORMAT'].split(':')
    frmt = ['HGVS'] + [frmt[0]] + ['GB']+ frmt[1:]
    vcf_df.loc[i,'FORMAT'] = ':'.join(frmt)
    sample_vals = vcf_df.at[i,sample_col].split(':')
    sample_vals = [','.join(hgvs)] + [sample_vals[0]]+['/'.join(genotype)]+sample_vals[1:]
    vcf_df.loc[i,sample_col] = ':'.join(sample_vals)

  vcf_header += '##FORMAT=<ID=HGVS,Number=A,Type=String,Description="HGVS nomenclature for the alternate allele (ALT)">\n'
  vcf_header += '##FORMAT=<ID=GB,Number=1,Type=String,Description="Genotype expressed as nucleotide bases (e.g., A/T or C|C)">\n'
  genotype_df = df.from_rows(genotype_df_rows,header=['HGVS',f'{sample_col}'],dfname=sample_col)
  return vcf_header, vcf_df, genotype_df


def mergegenotypedfs(df1:df,df2:df):
  df2cols = df2.columns.to_list()
  keys = df2['HGVS'].values
  df2genotype_col = df2cols[1]
  values = df2[df2genotype_col].values
  hgvs2genotype = dict(zip(keys, values))
  df1[df2genotype_col] = df1['HGVS'].map(hgvs2genotype)
  return df1


def trio_analysis(data_dir:str,out_dir:str,patient_id='WRFZO3066'):
  '''
  output:
    path to .tsv file containing genotypes found only in patient genome
  '''
  trio_genotype_dfs = []
  for i,vcf_path in dirList(data_dir,file_ext='vcf',include_subdirs=False):
    vcf_fname = fname(vcf_path)
    sample_name = vcf_fname[vcf_fname.find('_')+1: vcf_fname.find('.')]

    vcf_header, vcf_df = read_vcf(vcf_path)
    print(f'Processing file {vcf_fname} in {data_dir} with sample {sample_name}')
    vcf_header,vcf_df,genotype_df = HGVS2vcf(vcf_header,vcf_df,sample_name)

    annot_fname = os.path.join(out_dir,vcf_fname+'_annotated.vcf')
    open(annot_fname,'w').close() #flash
    with open(annot_fname, mode='a',encoding='utf-8') as f:
      f.write(vcf_header)
      f.write('#')
      vcf_df.to_csv(f,sep='\t',index=False,lineterminator='\n')
    trio_genotype_dfs.append(genotype_df)

  parent_dfs = []
  patient_df = None
  for i in range(0,len(trio_genotype_dfs)):
    if trio_genotype_dfs[i]._name_ == patient_id:
      patient_df = trio_genotype_dfs[i]
    else:
      parent_dfs.append(trio_genotype_dfs[i])

  trio_df = mergegenotypedfs(patient_df,parent_dfs[0])
  trio_df = mergegenotypedfs(trio_df,parent_dfs[1])
  path2trio = os.path.join(out_dir,'Trio genotypes.tsv')
  trio_df.to_csv(open(path2trio,mode='w',encoding='utf-8'),sep='\t',index=False,lineterminator='\n')

  parent_columns = parent_dfs[0].columns.to_list()+parent_dfs[1].columns.to_list()
  parent_columns = [c for c in parent_columns if c != 'HGVS']
  assert(len(parent_columns) == 2), f'Only two parents are possible, but you have: {parent_columns}'
  transmitted_genotypes = (trio_df[patient_id] == trio_df[parent_columns[0]]) | (trio_df[patient_id] == trio_df[parent_columns[1]])
  patient_uniqueGTs = trio_df[~transmitted_genotypes]
  path2unique_patient_genotypes = os.path.join(out_dir,'Only in patient.tsv')
  print(f'Found {len(patient_uniqueGTs)} genotypes usinque in patient')
  patient_uniqueGTs.to_csv(open(path2unique_patient_genotypes,mode='w',encoding='utf-8'),sep='\t',index=False,lineterminator='\n')
  return path2unique_patient_genotypes


def add_population_freq(df_path:str,patientid_column='WRFZO3066'):
  to_df = df.read(df_path,header=0)
  for i in to_df.index:
    genotype = to_df.at[i,patientid_column]
    if genotype[0] == genotype[1]: # is homozygous:
      allele_freqs = to_df.at[i,'dbsnp.alleles']
      if not pd.isna(allele_freqs):
        for d in ast.literal_eval(allele_freqs):
          if d['allele'] == genotype[0]:
            allele_freq = mean(d['freq'].values())
            to_df.loc[i,'Allele population frequence'] = allele_freq
            to_df.loc[i,'Genotype population frequence'] = allele_freq * allele_freq
  newdf_path = df_path[:df_path.rfind('.')] + 'freqs.tsv'
  to_df.to_csv(newdf_path,header=True,index=False,sep='\t')
  return to_df
  

MYVARIANT_FIELDS = (
    'dbsnp.rsid,'
    'clinvar.rcv.accession,'
    'dbnsfp.genename,'  # The Gene Symbol
    'snpeff.ann.effect,' # The consequence of the variant
    'dbsnp' # all dbSNP annotation fields 
)
dbsnp_columns = ['dbsnp.chrom','dbsnp.alleles','dbsnp.ref','dbsnp.alt','dbsnp.gene.geneid','dbsnp.gene.is_pseudo']
annotated_dfcols = ['HGVS','rsID','ClinVar ID','Gene','dbsnp.gene.symbol','dbsnp.gene.name','Predicted_Effect']+dbsnp_columns

def myvariant_annotate(genotype_df_path:str, output_dir:str):
  mv = myvariant.MyVariantInfo()
  print(f"Parsing {genotype_df_path}...")
  genotype_name = fname(genotype_df_path)
  genotype_df = df.read(genotype_df_path, header=0,name=genotype_name)
  hgvs_ids = genotype_df['HGVS'].to_list()
  print(f"Found {len(hgvs_ids)} variants. Querying MyVariant.info (hg38)...")
  start = time.time()
  annotated_pd = mv.getvariants(hgvs_ids, 
                      fields=MYVARIANT_FIELDS,
                      assembly='hg38',
                      as_dataframe=True)
  print(f'VCF was annoated by myvariant API in {execution_time(start)}')
  annotated_df = df.from_pd(annotated_pd,f'Annotated_{genotype_df._name_}')
  annotated_df = annotated_df.dfcopy(rename2={'dbsnp.rsid': 'rsID',
                                            'clinvar.rcv.accession': 'ClinVar ID',
                                            'dbnsfp.genename': 'Gene',
                                            'snpeff.ann.effect' : 'Predicted_Effect',
                                            '_id' :'HGVS'
                                            })
  annotated_df = annotated_df.merge_df(genotype_df,on='HGVS')
  annotated_df = add_population_freq(annotated_df)
  annotated_df = annotated_df.dfcopy(annotated_dfcols)
  output_csv_path = os.path.join(output_dir,genotype_name+'.myvariant.tsv')
  annotated_df.to_csv(output_csv_path, index=False,sep='\t',lineterminator='\n')
  print(f"Success! Annotated {len(annotated_df)} variants.\nData saved to {output_csv_path}")

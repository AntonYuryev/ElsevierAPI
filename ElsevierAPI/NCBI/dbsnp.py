import urllib.request,urllib.parse,os,time
from lxml import etree as ET
from urllib.error import HTTPError
from ..utils import dir2flist,execution_time,pretty_xml,next_tag,dir2flist,PCT
from ..ETM_API.references import Reference,TITLE,SENTENCE,AUTHORS,PUBYEAR,JOURNAL
from ..ResnetAPI.NetworkxObjects import PSRelation,OBJECT_TYPE
from ..ResnetAPI.ResnetAPIcache import APIcache, PSObject, ResnetGraph
import pandas as pd

BASE_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
CACHE_DIR = os.path.join(os.getcwd(),'ENTELLECT_API/ElsevierAPI/NCBI/__snpcache__/')
DBSNP_CACHE_BASE = 'dbsnp'

class SNP(PSObject):
  def __init__(self,name:str,id:str,database='dbSNP'):
    super().__init__()
    self.set_property('Name',name)
    self.Idendifires = {database:id}
    urn = 'urn:agi-gv-'+database.lower()+':'+id
    self.set_property('URN',urn)
    self.MAFs = dict() # {allele:(allele_count, pop_size)}

  def to_psobj(self):
    psobj = PSObject(self)
    [psobj.set_property(db+' ID',id) for db,id in self.Idendifires.items()]
    alleles = str()
    for allele, allele_count,pop_size in self.MAFs.items():
      alleles += allele + f'({round(float(allele_count/pop_size),2)}'+PCT+')'


def cache_path(basename=DBSNP_CACHE_BASE):
    cache_files = dir2flist(CACHE_DIR,fnames_has=basename,file_ext='.xml')
    return os.path.join(CACHE_DIR,DBSNP_CACHE_BASE+str(len(cache_files))+'.xml')


def gv2SNP(gv:PSObject,rsid2SNP:dict[str,SNP]):
    my_rsids = [gv.name()]+list(gv.get_props('Alias'))
    for rsid in my_rsids:
      if rsid.startswith('rs'):
          return rsid2SNP[rsid]
    return SNP()


def downloadSNP(rsids:list,mapdic:dict[str, dict[str, dict[str, PSObject]]])->tuple[dict[str,SNP],ResnetGraph]:
    """
    input:
        list of rs identifiers
    output:
        {snp_id:{allele:(allele_count, pop_size)}}, some {allele:(allele_count, pop_size)} may be empty
    """
    stepSize = 200
    params = {'db':'snp'}
    start = time.time()
    id2snp,gvs2genes = xmldir2SNP(CACHE_DIR,mapdic,rsids)
    rsids2download = list(set(rsids).difference(id2snp.keys()))
    print(f'{len(rsids)-len(rsids2download)} SNPs were found in cache.')
    rsids_len = len(rsids2download)
    if rsids_len:
        print(f'{len(rsids2download)} SNPs will be downloaded from dbSNP')
        with open(cache_path(),'w',encoding='utf-8') as f:
            f.write('<batch>\n')
            for i in range(0, rsids_len, stepSize):
                ids = ','.join(s[2:] for s in rsids2download[i:i+stepSize])
                
                params.update({'id':ids})
                req = urllib.request.Request(url=BASE_URL+'efetch.fcgi?'+urllib.parse.urlencode(params))
                for attempt in range(1,11):
                    try:
                        response = urllib.request.urlopen(req).read()
                        break
                    except HTTPError as e:
                        print(f'Fetch EUtils attempt {attempt} failed with error {e}.  Will try again in 5 seconds')
                        time.sleep(5)
                        continue

                snps = ET.fromstring('<documents>'+response.decode().strip()+'</documents>')
                f.write(pretty_xml(ET.tostring(snps),True))
                id2snps, gvs2genes_rels = xml2SNP(snps,mapdic)
                id2snp.update(id2snps)
                gvs2genes += gvs2genes_rels
                print(f'Downloaded {i+stepSize} SNPs out of {rsids_len}')
                time.sleep(1)
            f.write('</batch>')
    print(f'Download was done in {execution_time(start)}')
    return id2snp, ResnetGraph.from_rels(gvs2genes)
    

def dbsnp_hyperlink(rs_ids:list, as_count=True):
    """
    Input - list of rs identifiers\n
    Output - '=HYPERLINK(to search dbSNP with Input)
    """
    base_url = 'https://www.ncbi.nlm.nih.gov/snp/?'
    params = {'term':' OR '.join(rs_ids)}
    data = urllib.parse.urlencode(params, quote_via=urllib.parse.quote)
    if as_count:
        return '=HYPERLINK("'+base_url+data+'",\"{}\")'.format(str(len(rs_ids)))
    else:
        return '=HYPERLINK("'+base_url+data+'",\"{}\")'.format( ','.join(rs_ids))
    

def parseMAF(snp_record:ET._Element):
    '''
    output:
      {allele:(allele_count, pop_size)}
    '''
    alleles = dict() # {allele:(allele_count, pop_size)}
    for maf in snp_record.findall('./GLOBAL_MAFS/MAF'):
        #study = maf.find('STUDY').text
        allele_freq_pop = maf.find('FREQ').text
        eq_pos = allele_freq_pop.find('=')
        slash_pos = allele_freq_pop.find('/', eq_pos)
        allele = allele_freq_pop[:eq_pos]
        population_size = int(allele_freq_pop[slash_pos+1:])
        freq = float(allele_freq_pop[eq_pos+1:slash_pos])
        try:
            allele_count, pop_size = alleles[allele]
            alleles[allele] = (allele_count+freq*population_size, pop_size+population_size)
        except KeyError:
            alleles[allele] = (freq*population_size,population_size)

    to_return = dict()
    for a,f in alleles.items():
      frequency, population_size = f[0],f[1]
      if population_size:
        to_return[a] = str('%.4g' % float(f[0]/f[1]))
      else:
        to_return[a] = str(frequency)

    return to_return


def parseGene(snp_record:ET._Element)->list[tuple[str,str]]:
    '''
    output:
      [(gene_name,gene_id)...]
    '''
    genes = list() #[(gene_name, gene_id)]
    genes_elem = snp_record.find('./GENES')
    for gene in genes_elem.findall('GENE_E'):
        gene_name = gene.find('NAME').text
        geneid_element = gene.find('GENE_ID')
        gene_id = '' if geneid_element is None else geneid_element.text
        genes.append((gene_name,gene_id))

    return genes


def parseFunction(snp_record:ET._Element):
    fxns = list() 
    for fxn in snp_record.findall('./FXN_CLASS'):
        fn = str(fxn.text)
        if fn != 'None':
          fxns += fn.split(',')
    return fxns


def parseSNP(snpdocsum:ET._Element,mapdic:dict[str,dict[str,dict[str,PSObject]]])->tuple[SNP,list[PSRelation]]:
  gv2genes = list() #[(gene_name, gene_id)]
  snp_id = 'rs'+snpdocsum.find('./SNP_ID').text
  snp = SNP(snp_id,snp_id)
  snp.MAFs = parseMAF(snpdocsum)
  fxns = parseFunction(snpdocsum)
  if fxns:
    snp.update_with_list('Functional impact',fxns)
  genes = parseGene(snpdocsum)
  if genes:
    gv_psobj = snp.to_psobj()
    snp.update_with_list('Genes',genes)
    for gene_name,geneid in genes:
      gene = PSObject({'URN':'urn:agi-llid:'+geneid,'Name':gene_name})
      ref = Reference('PMID','11125122')
      ref[TITLE] = ['dbSNP: the NCBI database of genetic variation']
      ref[AUTHORS] = ['Sherry,S.T.;Ward,M.H.;Kholodov,M.;Baker,J.;Phan,L.;Smigielski,E.M.;Sirotkin,K']
      ref[PUBYEAR] = ['2001']
      ref[JOURNAL] = ['Nucleic Acids Res']
      sentence = 'Record from dbSNP links this Genetic Variant (SNV) to the indicated gene.'
      if fxns:
        sentence = ','.join(fxns)+' '+sentence
  
      ref.add_sentence_prop(ref._make_textref,SENTENCE,sentence)
      gv2gene = PSRelation.make_rel(gv_psobj,gene,{OBJECT_TYPE:'GeneticChange'},[ref])
      gv2genes.append(gv2gene)

  return snp, gv2genes


def xml2SNP(doc_summaries:ET._Element,mapdic: dict[str, dict[str, dict[str, PSObject]]],_4rsids:list=[])->tuple[dict[str,SNP],list[PSRelation]]:
    id2snp = dict()
    gvs2genes = list()
    if _4rsids:
      for snp_record in doc_summaries.findall('./DocumentSummary'):
          snp_id = snp_record.find('./SNP_ID')
          if snp_id is not None:
            snp_id = 'rs'+snp_record.find('./SNP_ID').text
            if snp_id in _4rsids:
              snp, gv2g = parseSNP(snp_record,mapdic)
              gvs2genes += gv2g
              id2snp[snp_id] = snp
    else:
      for snp_record in doc_summaries.findall('./DocumentSummary'):
        snp_id = snp_record.find('./SNP_ID').text
        if snp_id is not None:
          snp, gv2g = parseSNP(snp_record,mapdic)
          gvs2genes += gv2g
          id2snp[snp_id] = snp
                
    return id2snp,gvs2genes


def xmlfile2SNP(fname:str,mapdic:dict[str, dict[str, dict[str, PSObject]]],_4rsids:list=[])->tuple[dict[dict[str,SNP]],list[PSRelation]]:
  id2snp = dict()
  gvs2genes = list()
  print(f'Reading {fname} dbSNP dump file')
  for snps in next_tag(fname,"documents"):
      snp_dic,gv2gs = xml2SNP(snps,mapdic,_4rsids)
      gvs2genes += gv2gs
      id2snp.update(snp_dic)

  print(f'Loaded frequencies for {len(id2snp)} SNPs')
  return id2snp,gvs2genes


def xmldir2SNP(dirname:str,mapdic:dict[str, dict[str, dict[str, PSObject]]],_4rsids:list=[])->tuple[dict[dict[str,SNP]],list[PSRelation]]:
    id2snp = dict()
    listing = dir2flist(dirname,file_ext='xml')
    gvs2genes = list()
    for fname in listing:
        snp_dic,gv2gs = xmlfile2SNP(fname,mapdic,_4rsids)
        id2snp.update(snp_dic)
        gvs2genes += gv2gs
    return id2snp,gvs2genes


def __minor_allele(allele2freq:dict):
    sorted_dic = list(sorted(allele2freq.items(),key=lambda item: float(item[1])))
    return str(sorted_dic[0][0]),float(sorted_dic[0][1])


def genotype_frequency(genotype:str, allele2freq:dict)->tuple[float,str,float]:
    assert(len(genotype) == 2)
    minor_allele, minor_freq = __minor_allele(allele2freq)
    if minor_allele in genotype:
        if genotype[0] == genotype[1]:
            genotype_freq = minor_freq*minor_freq
        else:
            genotype_freq = 2*minor_freq*(1-minor_freq)
    else:
        major_freq = (1.0-minor_freq)
        genotype_freq = major_freq*major_freq
    
    return genotype_freq, minor_allele, minor_freq


def add_frequency_column(to_df:pd.DataFrame,rsid_colname:str,genotype_col:str):
  rsids = to_df[rsid_colname].to_list()
  rsids = [x for x in rsids if str(x).startswith('rs')]
  id2snp = downloadSNP(rsids)
  #assert isinstance(id2snp,dict[str,SNP])

  for idx in to_df.index:
    rsid = str(to_df.at[idx,rsid_colname])
    if rsid.startswith('rs'):
      genotype = to_df.at[idx,genotype_col]
      if len(genotype) == 2:
        if rsid in id2snp:
          snp = id2snp[rsid]
          assert isinstance(snp,SNP), 'Expect SNP object in id2snp dictionary'
          if snp.MAFs:
            genotype_freq, minor_allele, minor_freq = genotype_frequency(genotype,snp.MAFs)
            to_df.loc[idx,'Minor allele'] = minor_allele 
            to_df.loc[idx,'Minor allele frequency'] = float(minor_freq)
            to_df.loc[idx,'Genotype frequency'] = genotype_freq
          to_df.loc[idx,'Functional impact'] = ','.join(snp.get_props(['Functional impact']))
          genes = snp.get_props(['Genes'])
          gene_names = [g[0] for g in genes]
          to_df.loc[idx,'Genes'] = ','.join(gene_names)

  return to_df

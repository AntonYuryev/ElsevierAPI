import urllib.request, urllib.parse, os, time
from lxml import etree as ET
from urllib.error import HTTPError
from time import sleep
from ..utils import dir2flist, execution_time,pretty_xml,next_tag,dir2flist
from ..ResnetAPI.NetworkxObjects import PSObject

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
       self.MAFs = dict() # 


def cache_path(basename=DBSNP_CACHE_BASE):
    cache_files = dir2flist(CACHE_DIR,fnames_has=basename,file_ext='.xml')
    return os.path.join(CACHE_DIR,DBSNP_CACHE_BASE+str(len(cache_files))+'.xml')


def gv2SNP(gv:PSObject,rsid2SNP:dict[str,SNP]):
    my_rsids = [gv.name()]+list(gv.get_props(['Alias']))
    for rsid in my_rsids:
      if rsid.startswith('rs'):
          return rsid2SNP[rsid]
    return SNP()


def downloadSNP(rsids:list)->dict[str,SNP]:
    """
    input:
        list of rs identifiers
    output:
        {snp_id:{allele:(allele_count, pop_size)}}, some {allele:(allele_count, pop_size)} may be empty
    """
    stepSize = 200
    params = {'db':'snp'}
    start = time.time()
    
    id2snp = xmldir2SNP(CACHE_DIR,rsids)
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
                        sleep(5)
                        continue

                snps = ET.fromstring('<documents>'+response.decode().strip()+'</documents>')
                f.write(pretty_xml(ET.tostring(snps),True))
                id2snp.update(xml2SNP(snps))
                print(f'Downloaded {i+stepSize} SNPs out of {rsids_len}')
                sleep(1)
            f.write('</batch>')
    print(f'Download was done in {execution_time(start)}')
    return id2snp
    

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


def parseGene(snp_record:ET._Element):
    genes = list() #[(gene_name, gene_id)]
    genes_elem = snp_record.find('./GENES')
    for gene in genes_elem.findall('GENE_E'):
        gene_name = gene.find('NAME').text
        geneid_element = gene.find('GENE_ID')
        gene_id = '' if geneid_element is None else geneid_element.text
        genes.append((gene_name,gene_id))

    return genes


def parseFunction(snp_record:ET._Element):
    fxns = list() #[(gene_name, gene_id)]
    for fxn in snp_record.findall('./FXN_CLASS'):
        fxns += str(fxn.text).split(',')
    return fxns


def parseSNP(snp_record:ET._Element):
    snp_id = 'rs'+snp_record.find('./SNP_ID').text
    snp = SNP(snp_id,snp_id)
    snp.MAFs = parseMAF(snp_record)
    fxns = parseFunction(snp_record)
    if fxns:
        snp.update_with_list('Functional impact',fxns)
    genes = parseGene(snp_record)
    if genes:
        snp.update_with_list('Genes',genes)
    return snp


def xml2SNP(snps:ET._Element,_4rsids:list=[])->dict[dict[str,str]]:
    id2snp = dict()
    if _4rsids:
        for snp_record in snps.findall('./DocumentSummary'):
            snp_id = snp_record.find('./SNP_ID')
            if snp_id is not None:
                snp_id = 'rs'+snp_record.find('./SNP_ID').text
                if snp_id in _4rsids:
                    snp = parseSNP(snp_record)
                    id2snp[snp_id] = snp
    else:
        for snp_record in snps.findall('./DocumentSummary'):
            snp_id = snp_record.find('./SNP_ID').text
            if snp_id is not None:
                snp = parseSNP(snp_record)
                id2snp[snp_id] = snp
                
    return id2snp


def xmlfile2SNP(fname:str,_4rsids:list=[]):
    id2snp = dict()
    print(f'Reading {fname} dbSNP dump file')
    for snps in next_tag(fname,"documents"):
        snp_dic = xml2SNP(snps,_4rsids)
        id2snp.update(snp_dic)

    print(f'Loaded frequencies for {len(id2snp)} SNPs')
    return id2snp


def xmldir2SNP(dirname:str,_4rsids:list=[]):
    id2snp = dict()
    listing = dir2flist(dirname,file_ext='xml')
    for fname in listing:
        snp_dic = xmlfile2SNP(fname,_4rsids)
        id2snp.update(snp_dic)
    return id2snp


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

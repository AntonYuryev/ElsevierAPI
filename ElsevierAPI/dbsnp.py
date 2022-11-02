import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET

BASE_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'

def minor_allele(rsids:list):
    """
    Input - list of rs identifiers
    Output - {snp_id:{allele:(allele_count, pop_size)}}
    """
    
    stepSize = 200
    params = {'db':'snp'}
    
    snp2allele2freq = dict() # {snp_id:{allele:(allele_count, pop_size)}}
    for i in range(0, len(rsids), stepSize):
        ids = ','.join(s[2:] for s in rsids[i:i+stepSize])
        
        params.update({'id':ids})
        req = urllib.request.Request(url=BASE_URL+'efetch.fcgi?'+urllib.parse.urlencode(params))
        
        response = urllib.request.urlopen(req).read()
        snps = ET.fromstring('<documents>'+response.decode()+'</documents>')
        for snp_record in snps.findall('./DocumentSummary'):
            alleles = dict() # {allele:(allele_count, pop_size)}
            snp_id = snp_record.find('./SNP_ID')
            if type(snp_id) == type(None): continue
            snp_id = snp_record.find('./SNP_ID').text
            for maf in snp_record.findall('./GLOBAL_MAFS/MAF'):
                study = maf.find('STUDY').text
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

            alele_freqs = {a: str('%.4g' % float(f[0]/f[1])) for a,f in alleles.items() if f[1]>0}
            if alele_freqs:
                snp2allele2freq['rs'+snp_id] = alele_freqs

    return snp2allele2freq
    

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
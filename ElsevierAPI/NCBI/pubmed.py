
#C:Windows> py -m pip install entrezpy --user
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
import json

PUBMED_URL = 'https://pubmed.ncbi.nlm.nih.gov/?'
DOI_URL = 'http://dx.doi.org/'
baseURL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
#query = 'Aquilegia OR Arabidopsis OR Hordeum OR Brassica OR Theobroma OR Phaseolus OR Gossypium OR "Vitis vinifera" OR Mesembryanthemum OR Lactuca OR "Zea mays" OR Medicago OR "Nicotiana benthamiana" OR Allium OR Capsicum OR "Petunia hybrida" OR Pinus OR Strobus OR Solanum OR Oryza OR Secale OR "Sorghum bicolor" OR Picea OR Saccharum OR Helianthus OR Festuca OR Lycopersicon OR Triphysaria OR Triticum OR Populus OR anthocyanin OR "abscisic acid" OR brassinolide OR brassinosteroid OR brassinosteroids OR cytokinin OR cytokinins OR gibberellin OR gibberellins OR "gibberellic acid" OR "jasmonic acid" OR jasmonate OR jasmonates OR cultivar OR cultivars OR endosperm AND ((2021/05/01[PDat]:3019/12/31[PDat])'
#query = '(rice OR oryza OR sativa) AND (2021/05/01[PDat]:3019/12/31[PDat])'
#query= '(corn OR maize OR zea OR mays) AND (2021/05/01[PDat]:3019/12/31[PDat])'
#query= '(COVID19 OR SARS OR coronavirus) AND ("2020/09/17"[Date - Create] : "3000"[Date - Create])'
#query= "cytolethal distending toxin"
#query = '"inflammaging" OR "Inflamm-aging" OR "inflamm-ageing" OR "inflammageing" OR "age-increased inflammation" OR "age-associated inflammation" OR "age increased inflammation" OR "age associated inflammation"'
#query= '("lewy body" OR "lewy bodies" OR protofibrils OR protofibril)'
#query='(infection OR infectious OR pathogen OR parasite) AND (1989 [dp] OR 1994 [dp] OR 1999 [dp] OR 2004 [dp] OR 2009 [dp] OR 2014 [dp] OR 2019 [dp]) AND aug NOT (virus OR viral)'
#query = 'MAIT OR "Mucosa-associated invariant T" OR "Mucosa-associated invariant T cells" OR "Mucosa-associated invariant T cell"'
#query = 'bacillus OR subtilis AND ("2021/10/21"[Date - Create] : "3000"[Date - Create])'
#qname = 'Bsubtilis_since10212021'
#query = 'Clostridium AND ("2021/10/21"[Date - Create] : "3000"[Date - Create])'
#qname = 'Clostridium_since10212021'
#query = 'Cupriavidus OR necator AND ("2021/10/21"[Date - Create] : "3000"[Date - Create])'
#qname = 'Cnecator_since10212021'
query = 'escherichia OR coli AND ("2021/10/21"[Date - Create] : "3000"[Date - Create])'
qname = 'Ecoli_since10212021'
#query = 'H3K27M OR (K27M AND histone)'
#query = '(mutant OR mutated OR mutation OR mutations) AND (p53 OR tp53)'
#query = '("2019/11/29"[Date - Entry] : "3000"[Date - Entry])'
#query = '"Perineuronal nets" OR "Perineuronal net" OR PNNs'
#query = 'paxalisib OR GDC0084 OR "GDC-0084" OR P5DKZ70636 OR "EX-A1019"'
#query = '(sulfoquinovosyl OR monogalactosylmonoacylglycerol OR digalactosyldiacylglycerol OR monogalactosyldiacylglycerol OR digalactosylmonoacylglycerol OR monogalactosylmonoacylglycerol OR digalactosyldiacylglycerols OR monogalactosyldiacylglycerols OR digalactosylmonoacylglycerols) AND (human OR sapiens OR mammal OR rat OR mouse)'
#query = 'serine OR threonine OR tyrosine OR histidine'


def download_pubmed(query:str, query_name:str):
    """
    Input - query: NCBI query as used in UI
            query_name: used to make output files
    Output - query_name.json: has list of pmids found by query
           - query_name.xml: pubmed abstracts
    """

    params = {'db':'pubmed','term':query}

    count_param = dict(params)
    count_param.update({'rettype':'count'})
    req = urllib.request.Request(url=baseURL+'esearch.fcgi?'+urllib.parse.urlencode(count_param))
    response = ET.fromstring(urllib.request.urlopen(req).read())
    count = int(response.find('./Count').text)

    json_id_dum = query_name+'.json'
    try:
        all_pmids = json.load(open(json_id_dum,'r'))
        print('Loaded %d PMIDs from pmidlist.json' % len(all_pmids))
    except FileNotFoundError:
        all_pmids = list()
        params.update({'retmax':100000})
        for retstart in range(0, count, 100000):
            params.update({'retstart':retstart})
            req = urllib.request.Request(url=baseURL+'esearch.fcgi?'+urllib.parse.urlencode(params))
            response = ET.fromstring(urllib.request.urlopen(req).read())
            ids = response.findall('./IdList/Id')
            pmids = [x.text for x in ids]
            all_pmids += pmids
        json.dump(all_pmids, open(json_id_dum,'w'), indent = 2)

    #flushing content of "pumbed_results.xml"
    fname = query_name+".xml"
    file_result=open(fname, "w", encoding='utf-8')
    file_result.close()

    file_result=open(fname, "a", encoding='utf-8')
    file_result.write('<PubmedArticleSet>\n')

    params = {'db':'pubmed','rettype':'XML'}
    stepSize = 200
    for i in range(0, len(all_pmids), stepSize):
        ids = ','.join(str(s) for s in all_pmids[i:i+stepSize])
        params.update({'id':ids})
        req = urllib.request.Request(url=baseURL+'efetch.fcgi?'+urllib.parse.urlencode(params))
        abstractTree = ET.fromstring(urllib.request.urlopen(req).read())
        articles = abstractTree.findall('PubmedArticle')
        for article in articles:
            file_result.write(ET.tostring(article, encoding="unicode"))

    file_result.write('</PubmedArticleSet>\n')
    file_result.close()


def pubmed_hyperlink(pmids:list,display_str:int or str='',as_excel_formula=True):
    if as_excel_formula:
        ids4hyperlink = pmids[:20] 
        # hyperlink in Excel does not work with long URLs, 
        params = {'term':','.join(ids4hyperlink)}
        if display_str:
            to_display = str(display_str)
        elif len(pmids) == 1:
            to_display = str(pmids[0])
        else:
            to_display = str(len(pmids))

        data = urllib.parse.urlencode(params, quote_via=urllib.parse.quote)
        return '=HYPERLINK("'+PUBMED_URL+data+'",\"{}\")'.format(to_display)
    else:
        ids4hyperlink = pmids[:100]
        params = {'term':','.join(ids4hyperlink)}
        data = urllib.parse.urlencode(params, quote_via=urllib.parse.quote)
        return PUBMED_URL+data

if __name__ == "__main__":
    download_pubmed(query,qname)

'''
def minor_allele(rsids:list):
    """
    Input
    -----
    list of rsIDs
    
    Returns
    -------
    {rsID:Minor Allele:Frequency}
    """
    stepSize = 200
    params = {'db':'snp'}
    snp2allele2freq = dict() # {snp_id:{allele:(allele_count, pop_size)}}
    for i in range(0, len(rsids), stepSize):
        ids = ','.join(s[2:] for s in rsids[i:i+stepSize])

        baseURL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
        params.update({'id':ids})
        req = urllib.request.Request(url=baseURL+'efetch.fcgi?'+urllib.parse.urlencode(params))
        
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

            alele_freqs = {a:f[0]/f[1] for a,f in alleles.items() if f[1]>0}
            if alele_freqs:
                snp2allele2freq['rs'+snp_id] = alele_freqs

    return snp2allele2freq
'''
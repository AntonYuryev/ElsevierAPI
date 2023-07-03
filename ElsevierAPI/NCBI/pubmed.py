
#C:Windows> py -m pip install entrezpy --user
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
import json

PUBMED_URL = 'https://pubmed.ncbi.nlm.nih.gov/?'
DOI_URL = 'http://dx.doi.org/'
baseURL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
query = 'escherichia OR coli AND ("2021/10/21"[Date - Create] : "3000"[Date - Create])'
qname = 'Ecoli_since10212021'

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

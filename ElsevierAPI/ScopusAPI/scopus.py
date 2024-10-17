import urllib.request,urllib.parse,json
from xml.etree.ElementTree import fromstring
from time import sleep
from ..ETM_API.references import Reference, Author, DocMine
from ..ETM_API.references import SCOPUS_CI,ARTICLE_ID_TYPES,INSTITUTIONS,AUTHORS,PUBLISHER
from concurrent.futures import ThreadPoolExecutor, as_completed
from titlecase import titlecase


AUTHOR_SEARCH = 0
AUTHOR_RETRIEVAL = 1
SCOPUS_API_BASEURL = 'https://api.elsevier.com/content/'
SCOPUS_CACHE_DIR = 'D:/Python/ENTELLECT_API/ElsevierAPI/ScopusAPI/__scpcache__/'
# for Scopus query limits read https://dev.elsevier.com/guides/Scopus%20API%20Guide_V1_20230907.pdf
SCOPUS_AUTHORIDS = 'scopusAuthors'
SCOPUS_CITESCORE = 'CiteScore'
SCOPUS_SJR = 'Scientific Journal Rankings'
SCOPUS_SNIP = 'Source Normalized Impact per Paper'


class Scopus:
    base_url = 'https://api.elsevier.com/content/'
    page_size = 25

    @staticmethod
    def __journal_cache_name():
        return SCOPUS_CACHE_DIR+'JournalInfo.json'
    
    @staticmethod
    def __aff_cache_name():
        return SCOPUS_CACHE_DIR+'AffiliationInfo.json'
    
    
    def __init__(self,APIconfig:dict,add_param=dict()):
        self.params = {'apiKey':APIconfig['ELSapikey'], 'insttoken':APIconfig['insttoken'],'httpAccept':'application/json'}
        self.params.update(add_param)
        try:
            self.JournalInfo = json.load(open(self.__journal_cache_name(),'r',encoding='utf-8'))
        except FileNotFoundError:
            self.JournalInfo = dict()

        try:
            self.AffiliationInfo = json.load(open(self.__aff_cache_name(),'r',encoding='utf-8'))
        except FileNotFoundError:
            self.AffiliationInfo = dict()
        

    def _get_param_str(self):
        return urllib.parse.urlencode(self.params)


    def _url_request(self):
        return self.base_url+self._get_param_str()


    def _get_results(self,sleep_time = 0.3):
        for attempt in range(0,1):
            try:
                url_reg = self._url_request()
                response = urllib.request.urlopen(url_reg)
                scopus_view = response.read().decode('utf-8')
                result = json.loads(scopus_view)
                sleep(sleep_time)
                return result
            except Exception:
                sleep(sleep_time)
    
        return dict()
    

    def close(self):
        def write_journal_info():
            with open(self.__journal_cache_name(), 'w',encoding='utf-8') as f:
                json.dump(self.JournalInfo,f,indent=1)

        def write_instituion_info():
            with open(self.__aff_cache_name(), 'w',encoding='utf-8') as f1:
                json.dump(self.AffiliationInfo,f1,indent=1)        

        with ThreadPoolExecutor(2, thread_name_prefix='CloseScopusCache') as e:
                e.submit(write_journal_info)
                e.submit(write_instituion_info)
        e.shutdown()


    @staticmethod
    def oa_status(doi:str):
        url = 'https://api.elsevier.com/content/abstract/doi/'+doi
        try:
            result = urllib.request.urlopen(url).read().decode('utf-8')
            record = fromstring(result)
            ns = {
                "abstract": "http://www.elsevier.com/xml/svapi/abstract/dtd"
            }
            oa = record.find('abstract:coredata',ns).find('abstract:openaccess',ns).text
            return bool(oa)
        except Exception as e:
                return True
        
    
    def is_in_open_access(self,ref:Reference):
        doi = ref.identifier('DOI')
        return Scopus.oa_status(doi) if doi else True
    

    def __journal_info(self,issn:str,j_title=''):
        try:
            return self.JournalInfo[issn]
        except KeyError:
            self.base_url = SCOPUS_API_BASEURL+'serial/title/issn/'+issn+'?'
            self.params.update({'view':'CITESCORE'})
            #result_str = urllib.request.urlopen(url).read().decode('utf-8')
            #result = json.loads(result_str)
            result = self._get_results()
            if result:
                entry = result['serial-metadata-response']['entry'][0]
                assert( isinstance(entry,dict))
                j_title = entry['dc:title']
                publisher = str(entry['dc:publisher'])
                CiteScore = entry.get('citeScoreYearInfoList',{}).get('citeScoreCurrentMetric','')
                if isinstance(CiteScore,str):
                    CiteScore = float(CiteScore) if CiteScore else 0.0
                SJRscore = entry.get('SJRList',{}).get('SJR',[])
                if SJRscore:
                    SJRscore = SJRscore[0]
                    SJRscore = SJRscore['$'] + ' ('+SJRscore['@year']+')'
                SNIPscore = entry.get('SNIPList',{}).get('SNIP',[])
                if SNIPscore:
                    SNIPscore = SNIPscore[0]
                    SNIPscore = SNIPscore['$'] + ' ('+SNIPscore['@year']+')'
                record = [j_title,publisher,CiteScore,SJRscore,SNIPscore]
            else:
                record = [j_title,f'with ISSN {issn} has no Scopus record' ,'','','']
                print(f'{j_title} with ISSN {issn} has no Scopus record')

            self.JournalInfo[issn] = record
            self.params.pop('view','')
            return record


    def get_affiliation(self,institution:str):
        '''
        output:
            canonical institution name for input institution name.  In case of multiple hits returns name of the most relevant Scopus hit  
        '''
        titlecase_institution = titlecase(institution)
        try:
            return str(self.AffiliationInfo[titlecase_institution])
        except KeyError:
            self.base_url = SCOPUS_API_BASEURL+'search/affiliation?'
            query = f'affil({institution})'
            self.params.update({'query': query,'sort':'relevancy'}) # ensures the most relevant hit to be 1st
            result = self._get_results()
            if result:
                entry = result["search-results"]["entry"][0]
                if 'error' in entry: return ''
                affil_name = str(entry['affiliation-name'])
                affil_variants = entry['name-variant']
                for variant in affil_variants:
                    self.AffiliationInfo[titlecase(variant['$'])] = affil_name
                    self.AffiliationInfo[titlecase_institution] = affil_name
                return affil_name
            else:
                return ''


    def normalize_affiliations(self,ref:Reference):
        '''
        input:
            list of not normalized institution name
        output:
            list of institution names normalized by Scopus
        '''
        affils = ref.get(INSTITUTIONS,[])
        normalized_institutions = list()
        for aff in affils:
            norm_aff = self.get_affiliation(aff)
            if norm_aff:
                normalized_institutions.append(norm_aff)

        if normalized_institutions:
            ref[INSTITUTIONS] = normalized_institutions
        return normalized_institutions
    

    def scopus_stats4(self,ref:Reference)->tuple[str,str,str,str,str]:
        '''
        output:
            journal_title,publisher,CiteScore,SJRscore,SNIPscore
        '''
        for issn in ref.get_props('ISSN'):
            j_title,publisher,CiteScore,SJRscore,SNIPscore = self.__journal_info(issn,ref.journal())
            if publisher:
                ref[PUBLISHER] = [publisher]
                ref[SCOPUS_CITESCORE] = [CiteScore]
                ref[SCOPUS_SJR] = SJRscore
                ref[SCOPUS_SNIP] = SNIPscore

                return j_title,publisher,CiteScore,SJRscore,SNIPscore
        return '','','','',''
    

    def citation_overview(self,ref:Reference):
        id_type,identifier = ref.get_doc_id()
        id_name = 'pubmed_id' if id_type == 'PMID' else id_type
        self.url = 'https://api.elsevier.com/content/abstract/citations?'
        self.params.update({id_name:identifier})

        result = self._get_results()
        return result["abstract-citations-response"]
    

class AuthorRetreival(Scopus):
    def __init__(self,author_id:str, APIconfig:dict):
        self.base_url = super().base_url+'author/author_id/'+author_id+'?'
        super().__init__(APIconfig)


class AuthorSearch(Scopus):
    author_cache = dict()

    def __init__(self,APIconfig:dict,add_param=dict()):
        super().__init__(APIconfig,add_param)
        self.base_url = super().base_url+'search/author?'

    
    def __parse_results(self,author:dict):
        '''
        Return
        ------
        author_id,given_name,initial,surname,affiliation
        '''
        dc_identifier = str(author["dc:identifier"])
        author_id = dc_identifier[dc_identifier.find(':')+1:]
        author_name = author["preferred-name"]
        affiliation = str(author['affiliation-current']['affiliation-name'])
        given_name = author_name['given-name']
        initials = author_name['initials']
        surname = author_name['surname']
        self.author_cache[affiliation] = {surname:{given_name:author_id}}
        return int(author_id),str(given_name),str(initials),str(surname),affiliation


    def get_author_id(self, auLast:str, au1st='', institutions=[])->tuple[int,str,str,str]:
        '''
        Return
        ------
        [author_id,author['given-name'],author['initials'],author['surname'],institution]
        '''
        for institution in institutions:
            try:
                author_id = self.author_cache[institution][auLast][au1st]
                return author_id,au1st,auLast,institution
            except KeyError:
                continue
        
        query = 'authlast({})'.format(auLast)
        if au1st: query += ' and authfirst({})'.format(au1st)
        self.params.update({'query': query})
        authors = self._get_results()
        if not authors: return (0,'','','')
        if len(authors) == 1:
            au_info = self.__parse_results(authors[0])
            self.author_cache[auLast+' '+au1st] = au_info
            return au_info
        else:
            if institutions:
                for institution in institutions:
                    #sleep(1)
                    query_with_aff = query + f' and affil({institution})'
                    self.params['query'] = query_with_aff
                    authors = self._get_results()
                    if len(authors) == 1:
                        au_info = self.__parse_results(authors[0])
                        self.author_cache[auLast+' '+au1st] = au_info
                        return au_info
                    else:       
                        continue

                insts = ','.join(institutions)
                print(f'Author info {auLast},{au1st} from {insts} is ambigious')
                return (0,'','','')
            else:
                print('Author info %s,%s is ambigious' % (auLast,au1st))
                return (0,'','','')


    def get_authors(self,ref:DocMine):
        authors = ref[AUTHORS]
        instituts = ref[INSTITUTIONS]
        for author in authors:
            au_names = str(author).split(' ')
            surname = au_names[0]
            _1stname = au_names[1] if len(au_names) > 1 else ''
            scopus_author = self.get_author_id(surname,_1stname,instituts)
            ref.authors.append(Author(surname,_1stname,instituts[0]))
        return scopus_author


    @staticmethod
    def authorid(author_info:tuple):
        return author_info[0][10:]


class ScopusSearch(Scopus):

    def add_query(self,query:str):
        self.params.update({'query':query})

    def __init__(self,query:str,APIconfig:dict):
        add_param = {'query':query}
        super().__init__(APIconfig,add_param)
        self.base_url = super().base_url+'search/scopus?'

    def _get_results(self):
        result = super()._get_results()
        hit_count = int(result['search-results']['opensearch:totalResults'])
        articles = list(result["search-results"]["entry"])
        for s in range (self.page_size, hit_count, self.page_size):
            self.params.update({'start':str(s), 'count':str(self.page_size)})
            result = super()._get_results()
            articles += list(result["search-results"]["entry"])

        return articles
    

    @staticmethod
    def biblio_str(article:dict):
        title = article['dc:title']
        journal = article['prism:publicationName']
        volume = article['prism:volume']
        issue = article['prism:issueIdentifier']
        pages = article['prism:pageRange']
        pub_date = article['prism:coverDisplayDate']
        return str('"'+title+'" '+journal+' v.'+volume+':'+issue+'('+pages+'), '+pub_date)
    

    @staticmethod
    def get_doc_id(article:dict):
        try:
            return str('PMID'), str(article["pubmed-id"])
        except KeyError:
            try:
              return str('DOI'), str(article["prism:doi"])
            except KeyError:
                return str(),str() 


def g_index(APIconfig:dict,last_name:str,first_name:str,institution = ''):
    au = AuthorSearch(APIconfig)
    author = au.get_author_id(last_name,first_name,institution)
    query = f'AU-ID({au.authorid(author)}'
    search = ScopusSearch(query,APIconfig)
    author_articles = list(search._get_results())
    author_articles.sort(key=lambda x: int(x['citedby-count']), reverse=True)
    
    citation_sum = 0
    gindex = len(author_articles)
    index = 0
    for article in author_articles:
        citation_sum += int(article['citedby-count'])
        index += 1
        if index*index >= citation_sum:
            gindex = index+1
            break
    
    return gindex,len(author_articles),search.biblio_str(author_articles[0])


def loadCI(APIconfig:dict, references:set[Reference])->tuple[set[Reference], set[Reference]]:
    '''
    Return
    ------
    articles_with_ci - {Reference} annotated with [SCOPUS_CI]\n
    no_ci_articles - {Reference}
    '''
    idtype2id2ref = dict()
    for ref in references:
        id_type,refid = ref.get_doc_id()
        if id_type in ARTICLE_ID_TYPES:
            try:
                idtype2id2ref[id_type][refid] = ref
            except KeyError:
                idtype2id2ref[id_type] = {refid:ref}

    with ThreadPoolExecutor(100, thread_name_prefix='ScopusCI') as e:
        scopus_futures = list()
        for id_type, id2ref in idtype2id2ref.items():
            ids = list(id2ref.keys())
            for i in range(0,len(ids),100):
                query = '{id_type}(' + ' OR '.join(ids[i:i+100])+')'
                query = query.format(id_type=id_type)
                search = ScopusSearch(query,APIconfig)
                scopus_futures.append(e.submit(search._get_results))

        scopus_articles = [f.result() for f in as_completed(scopus_futures)]
        e.shutdown()

    for article in scopus_articles:
        id_type,id = ScopusSearch.get_doc_id(article)
        try:
            idtype2id2ref[id_type][id][SCOPUS_CI] = article['citedby-count']
        except KeyError:
            continue
    
    articles_with_ci = set()
    no_ci_articles = set()
    for id_type, id2ref in idtype2id2ref.items():
        for id, ref in id2ref.items():
            assert(isinstance(ref, Reference))
            [no_ci_articles,articles_with_ci][ref.has_property(SCOPUS_CI)].add(ref)

    return articles_with_ci, no_ci_articles



'''
if __name__ == "__main__":
    sc = CitationOverview('scopus_id','28773',load_api_config())
    co = sc._get_results()
    print('')
'''


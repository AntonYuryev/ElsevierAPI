import urllib.request
import urllib.parse
import urllib.error as http_error
import json
from time import sleep

AUTHOR_SEARCH = 0
AUTHOR_RETRIEVAL = 1

class Scopus:
    url = 'https://api.elsevier.com/content/'
    page_size = 25
    
    def __init__(self,APIconfig:dict,add_param=dict()):
        self.params = {'apiKey':APIconfig['ELSapikey'], 'insttoken':APIconfig['insttoken'],'httpAccept':'application/json'}
        self.params.update(add_param)
        

    def _get_param_str(self):
        return urllib.parse.urlencode(self.params)

    def _url_request(self):
        return self.url+self._get_param_str()
    
    def _get_results(self):
        try:
            url_reg = self._url_request()
            response = urllib.request.urlopen(url_reg)
        except http_error.HTTPError:
            sleep(30)
        scopus_view = response.read()
        result = json.loads(scopus_view.decode('utf-8'))
        sleep(0.34)
        return result

class AuthorRetreival(Scopus):
    def __init__(self,author_id:str, APIconfig:dict):
        self.url = super().url+'author/author_id/'+author_id+'?'
        super().__init__(APIconfig)

class AuthorSearch(Scopus):
    author_cache = dict()

    def __init__(self,APIconfig:dict,add_param=dict()):
        super().__init__(APIconfig,add_param)
        self.url = super().url+'search/author?'

    @staticmethod
    def __parse_results(author:dict):
        author_id = author["dc:identifier"]
        author_name = author["preferred-name"]
        return author_id,author_name['given-name'],author_name['initials'],author_name['surname']

    def _get_results(self):
        result = super()._get_results()
        return result["search-results"]["entry"]

    def get_author_id(self, auLast:str, au1st='', institutions=[]):
        query = 'authlast({})'.format(auLast)
        if au1st: query += ' and authfirst({})'.format(au1st)
        self.params.update({'query': query})
        try:
            return self.author_cache[auLast+' '+au1st]
        except KeyError:
            authors = self._get_results()
            if not authors: return []
            if len(authors) == 1:
                au_info = self.__parse_results(authors[0])
                self.author_cache[auLast+' '+au1st] = au_info
                return au_info
            else:
                if institutions:
                    for institution in institutions:
                        #sleep(1)
                        query += ' and affil({})'.format(institution)
                        authors = self._get_results()
                        if len(authors) == 1:
                            self.author_cache[auLast+' '+au1st] = au_info
                            return self.__parse_results[authors[0]]
                        else:       
                            continue
                    print('Author info %s,%s from *s is ambigious' % (auLast,au1st, ','.join(institutions)))
                    return []
                else:
                    print('Author info %s,%s is ambigious' % (auLast,au1st))
                    return []

    @staticmethod
    def authorid(author_info:tuple):
        return author_info[0][10:]

class ScopusSearch(Scopus):
    def __init__(self,query:str,APIconfig:dict):
        add_param = {'query':query}
        super().__init__(APIconfig,add_param)
        self.url = super().url+'search/scopus?'

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
    def biblio_str(article):
        title = article['dc:title']
        journal = article['prism:publicationName']
        volume = article['prism:volume']
        issue = article['prism:issueIdentifier']
        pages = article['prism:pageRange']
        pub_date = article['prism:coverDisplayDate']
        return '"'+title+'" '+journal+' v.'+volume+':'+issue+'('+pages+'), '+pub_date


class CitationOverview(Scopus):
    def __init__(self,identifier_type:str,identifier:str,APIconfig:dict):
        add_param = {identifier_type:identifier}
        super().__init__(APIconfig,add_param)
        self.url = 'https://api.elsevier.com/content/abstract/citations?'

    def _get_results(self):
        result = super()._get_results()
        return result["abstract-citations-response"]




def g_index(APIconfig:dict,last_name:str,first_name:str,institution = ''):
    au = AuthorSearch(APIconfig)
    author_id = au.get_author_id(last_name,first_name,institution)
    aid = au.authorid(author_id)
    query = 'AU-ID({})'.format(aid)
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



            



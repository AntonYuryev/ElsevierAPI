from urllib.error import HTTPError
import json, http.client, time, urllib.request, urllib.parse,ssl,re
from collections import defaultdict
from ..utils import  load_api_config
from time import sleep
from .references import AUTHORS,INSTITUTIONS,JOURNAL,SENTENCE,EMAIL,RELEVANCE,PUBLISHER,GRANT_APPLICATION
from .references import DocMine,Reference

DEFAULT_ETM = 'https://covid19-services.elseviertextmining.com/api'

def remove_the(t:str):
    return t[4:] if t.startswith('The ') else t


class ETMjson(DocMine):
    def __init__(self,article:dict):
        is_article, article_ids = self.__parse_ids(article)
        if is_article:
            id_type, id_value = article_ids.popitem()
            if id_type == 'ELSEVIER': id_type = 'PII'
            super().__init__(id_type, id_value)
            self.Identifiers.update(article_ids)
        else:
            return None # __init__ must return None not dict()
            
        try:
            title = article['article']['front']['article-meta']['title-group']['article-title']
            self._set_title(title)
        except KeyError: pass

        try:
            year = article['article']['front']['article-meta']['pub-date']['year']
            self.set_date(year)
        except KeyError: pass

        self.add2section('Abstract',' '.join(self.__parse_abstact(article)))
        authors,institutions,self.addresses = ETMjson.__parse_contributors(article)
        if authors: self[AUTHORS] = authors
        if institutions: self[INSTITUTIONS] = institutions

        try:
            journal = article['article']['front']['journal-meta']['journal-title-group']['journal-title']['content']
            #journal = self.__parse_journal(journal)
        except KeyError:
            journal = GRANT_APPLICATION
        self.append_property(JOURNAL, journal)

        issn = article.get('article', {}).get('front', {}).get('journal-meta', {}).get('issn', {}).get('x', '')
        if issn:
            self.append_property('ISSN', issn)

        publisher = article.get('article', {}).get('front', {}).get('journal-meta', {}).get('publisher',{}).get('publisherNameAndPublisherLoc',{}).get('publisher-name',{}).get('content', '')
        if publisher:
            self.append_property(PUBLISHER, publisher)

        self['categories'] = self.__parse_category(article)
        text_ref = self._make_textref()
        snippet, self.term_ids, self.keywords = self.__parse_snippet(article) #term_ids = {term_id+'\t'+term_name}
        self.add_sentence_prop(text_ref, SENTENCE, snippet)


    @staticmethod
    def dump_fname(search_name): return search_name + '.json'

    @staticmethod
    def __parse_ids(article:dict):
        is_article = True
        article_ids = dict() #{id_type:[id_values]}
        ids = article['article']['front']['article-meta']['article-id']
        if not isinstance(ids,list): ids = [ids]
        for article_id in ids:
            id_type = article_id['pub-id-type'].upper()
            id_value = str(article_id['value'])
            try:
                if id_type == 'ELSEVIER': continue
                exist_ids = article_ids[id_type]
                if isinstance (exist_ids, str) : 
                    exist_ids = [exist_ids]
                exist_ids.append(id_value)
                article_ids[id_type] = exist_ids
                is_article = False #record has multiple IDs of the same type.  
                # It is likely abstract from conference proceedings
                # Reference functions will not work.  Ref must be ignored in subsequent code
            except KeyError:
                article_ids[id_type] = id_value
        
        return is_article, article_ids


    @staticmethod
    def __parse_category(article:dict):
        try:
            categories = article['article']['front']['article-meta']['article-categories']['subj-group']['subject']
            if isinstance(categories,list):
                clean_categories = set()
                for category in categories:
                    category = str(category).title()
                    if category == 'Article':
                        category = 'Journal Article'
                    clean_categories.add(category)
                return list(clean_categories)
            else:
                categories = str(categories).title()
                return [categories]
        except KeyError:
            return list()


    @staticmethod
    def __parse_abstact(article:dict):
        try: abstract_secs = article['article']['front']['article-meta']['abstract']['sec']
        except KeyError: 
            try: abstract_secs = article['article']['body']
            except KeyError:  return [str()]

        if isinstance(abstract_secs, dict):
            try:
                abstract_secs = abstract_secs['sec']
                if isinstance(abstract_secs, dict): 
                    abstract_secs = [abstract_secs]
            except KeyError:
                try:
                    return [str(abstract_secs['addressOrAlternativesOrArray']['p'])]
                except KeyError:
                    print('Article abstract has no addressOrAlternativesOrArray', abstract_secs)
                    return [str()]
        
        #at this point abstract_secs can be only list of dict
        abstract_sentences = list()
        for paragraph in abstract_secs:
            try: abstract_sentences.append(str(paragraph['title'])) #some abstracts are secionalized and have titles
            except KeyError: pass

            try:
                abstract_sentences.append(str(paragraph['addressOrAlternativesOrArray']['p']))
            except KeyError:
                try:
                    scs = paragraph['sec']['sec'] if isinstance(paragraph['sec'], dict) else paragraph['sec']
                    if isinstance(scs,dict):
                        scs = [scs]
                    for s in scs:
                        try:
                            abstract_sentences.append(str(s['addressOrAlternativesOrArray']['p']))
                        except KeyError:
                            continue
                except KeyError: continue
                        
        return abstract_sentences


    @staticmethod
    def __parse_institution(correspondence:str)->tuple[list[str],str,str]:
        email_pattern = EMAIL.search(correspondence)
        if email_pattern:
            email_pos = email_pattern.start()
            email = correspondence[email_pos:email_pattern.end()]
            for i in range(3, email_pos-3):
                sub_start = email_pos-i
                sub = correspondence[sub_start : sub_start+3]
                if sub in {'. e', '. E'}: # . Electronic address
                    email_pos = sub_start
                    break
            full_address = correspondence[0:email_pos].strip(' .,')
        else:
            full_address = correspondence.strip(' .,')
            email = ''
        
        names = list()
        if full_address:
            # remove leading non-letters:
            full_address = re.sub(r'^[^a-zA-Z]+', '', full_address)
            names = full_address.split(', ')

        affiliations = list()
        if names:
            has_institution = False
            for n,name in enumerate(reversed(names)):
                if DocMine.has_institution_keyword(name):
                    has_institution = True
                    break

            if not has_institution:
                #print(f'article coresspondence "{correspondence}" has no institutional keyword!') # makes a lot of noise
                if len(names) > 2: #evidence that correspondence is not junk
                    affiliations.append(names[0].title())
            else:
                affiliations = names[:len(names)-n]
        return affiliations, full_address, email


    @staticmethod
    def __parse_contributors(article:dict)->tuple[list[str],list[str],list[str]]:
        '''
        Return
        ------
        [authors],[institutions], [addresses] where\n
        "institutions" - [global_affiliation]\n
        [addresses] - [detailed_affiliation]
        '''
        try:
            author_info = article['article']['front']['article-meta']['contribGroupOrAffOrAffAlternatives']
            if not author_info:
                return [],[],[]
        except KeyError:
            return [],[],[]
        if isinstance(author_info, dict): author_info = [author_info]

        authors = list()
        institutions = set()
        addresses = list()
        assert(len(author_info) <= 2) # either aff or contrib-group
        for item in author_info:
            # list of affiliations is not linked to list of contributors unfortunately in ETM
            try:
                affs = item['aff']['content']
                if not isinstance(affs, list): affs = [affs]
                for aff in affs:
                    affiliations, address, email = ETMjson.__parse_institution(aff['institution'])
                    addresses.append(address)
                    if affiliations:
                        institutions.add(affiliations[-1]) #use only last institution name defining global organization
            except KeyError:
                contribOrAddressOrAff = item['contrib-group']['contribOrAddressOrAff']
                if not isinstance(contribOrAddressOrAff, list): contribOrAddressOrAff = [contribOrAddressOrAff]

                au_name = ''
                for contributor in contribOrAddressOrAff:
                    try:
                        author = contributor['contrib']['contribIdOrAnonymousOrCollab']['name']
                    except KeyError:
                        continue
                    try:
                        given_names = author['given-names']
                        try:
                            au1stname = str(given_names['content'])
                        except KeyError:
                            try:
                                au1stname = str(given_names['initials'])+'.'
                            except KeyError:
                                au1stname = ''
                    except KeyError:
                        given_names = dict()
                        au1stname = ''

                    try:
                        au_name = author['surname']
                        if au1stname: au_name += ' ' + au1stname
                        au_name = str(au_name).title()
                        authors.append(au_name)
                        # list of affiliations is not linked to list of contributors unfortunately in ETM
                    except KeyError: 
                        continue

        return authors,list(institutions),addresses #authors: list(LastName 1stName), institutions: list([institution_name,address])


    @staticmethod
    def __parse_snippet(article:dict, markup_all=False):
        try:
            data = article['data']
            snippet = str(data['snippet']['text'])
            term_ids = set()
            keywords = set()
            pos_shift = 0
            highlights = article['data']['snippet']['highlight']
            if not isinstance(highlights, list): highlights = [highlights]
            last_markup_end = 0
            for markup in highlights:
                try:
                    query_end = pos_shift + markup['end']
                    if query_end <= last_markup_end: 
                        continue

                    query_term = markup['queryTerm']
                    query_start = pos_shift+markup['start']
                    snippet = snippet[:query_start] + '{'+query_term+'}='+snippet[query_start:]
                    pos_shift += len(query_term)+3
                    last_markup_end = query_end+pos_shift
                
                    try:
                        terms = markup['term']
                        if isinstance(terms,dict):
                            term_ids.add(terms['id']+'\t'+terms['value'])
                        else: 
                            [term_ids.add(term['id']+'\t'+term['value']) for term in terms]
                    except KeyError: continue
                except KeyError:
                    if markup_all:
                        keyword_end = pos_shift + markup['end']
                        if keyword_end <= last_markup_end: continue

                        keywords.add(markup['text'])
                        keyword_start = pos_shift + markup['start']
                        snippet = snippet[:keyword_start] + '{'+snippet[keyword_start:keyword_end]+'}'+snippet[keyword_end:]
                        pos_shift += 2
                        last_markup_end = keyword_end+pos_shift
                    else:
                        continue

            return snippet,term_ids,keywords
        except KeyError:
            try:
                article_body = article['article']
                sec = article_body['body']['sec']
                if isinstance(sec,dict):
                    return sec['addressOrAlternativesOrArray']['p'], set(), set()
                else:
                    snippets = str()
                    for s in sec:
                        snippets += s['addressOrAlternativesOrArray']['p']+'...'
                    return snippets, set(), set()
            except KeyError:
                return '', set(), set()


class ETMsearch:
    def __init__(self,APIconfig:dict=dict(), **kwargs):
        self.APIconfig = APIconfig if APIconfig else load_api_config()
        self.url = self.APIconfig.get('ETMURL',DEFAULT_ETM)
        self.params = {
            'searchTarget': 'full_index',
            'apikey':self.APIconfig['ETMapikey'],
            'limit': kwargs.get('limit',5),
            'insttoken':APIconfig['insttoken'] # ScopusAPI needs it
            #'so_c':'17~' # filter by Publication Source excluding NIH reporter; 
            # 'so_c':'17~' does not work on COVID server
            #'so_p':'ffb~' # filter by Publication Type
            }
        self.params.update(kwargs.get('add_param',dict()))
        self.min_relevance = kwargs.get('min_relevance',0.0) 
        self.hit_count = 0
        self.page_size = 100
        self.request_type = '/search/basic?'  # '/search/advanced?'


    def __base_url(self): 
        return self.url+self.request_type


    def clone(self, to_url=''):
        myApiconfig = dict(self.APIconfig)
        myApiconfig['ETMURL'] = to_url if to_url else self.url
        newEtMsearch =  ETMsearch(myApiconfig,self._limit(),self.params)
        newEtMsearch.page_size = self.page_size
        newEtMsearch.request_type = self.request_type
        newEtMsearch.hit_count = 0
        newEtMsearch.min_relevance = self.min_relevance
        return newEtMsearch


    def _limit(self): 
        """
        Returns
        -------
        self.params['limit'] - max number of reference to retreive
        """
        return self.params['limit']

    def references(self)->set[DocMine]:
        return set([x[0] for x in self.ref_counter.values()])
 
    def _set_query(self,query:str):
        self.params['query'] = query

    def _query(self):
        try:
            return self.params['query']
        except KeyError:
            print('ETMsearch class instance has no query specified')

    def __get_param_str(self):
        return urllib.parse.urlencode(self.params)

    def _url_request(self):
        return self.__base_url()+self.__get_param_str()

    def search_reviews_only(self):
        self.params.update({'so_p': '1~'})

    def use_advanced_search(self):
        self.request_type = '/search/advanced?'

    @staticmethod
    def article2ref(article:dict):
        return Reference(ETMjson(article))


    def _get_articles(self, page_start=0, need_snippets=True)->tuple[list,int]:
        """
        Returns tuple
        -------
        result['article-data'], result['total-hits=']
        """
        if page_start: self.params['start'] = page_start
        if need_snippets:
            self.params.update({'snip': '1.desc'})
        for attempt in range(1, 11):
            try:
                context = ssl._create_unverified_context()
                #context=ssl.create_default_context(cafile=certifi.core.where())
                the_page = urllib.request.urlopen(self._url_request(),context=context).read()         
                if the_page:
                    result = json.loads(the_page.decode('utf-8'))
                    if attempt > 1:
                        print(f'ETM connection was restored on the {attempt} attempt')
                    sleep(5)
                    return list(result['article-data']), int(result['total-hits='])
                else:
                    return list(),int(0)
            except HTTPError as err:
                raise err # usually means that ETM is not available
                timeout = 10
                print(f'ETMStat HTTPError. Will attempt to reconnect in {timeout}')
                time.sleep(timeout)
            except http.client.IncompleteRead:
                timeout = 10
                print(f'ETMStat IncompleteRead. Will attempt to reconnect in {timeout}')
                time.sleep(timeout)

        print(f'Cannot connect to {self.__base_url()} after 10 attempts')
        return list(),int(0)


    @staticmethod
    def basic_query(entity1:str,entity2:str,add2query:list):
        '''
        Returns
        -------
        ETM query for  basic serach where "terms" are joined with ";"
        '''
        return ';'.join([entity1,entity2]+add2query)
        #return '{'+'};{'.join([entity1,entity2]+add2query)+'}'


    def __search(self)->tuple[int,dict[str,list],list[Reference]]:
        """
        output:
        Tuple:
            [0] hit_count - TOTAL number of reference found by ETM basic search 
            [1] ref_ids = {id_type:[identifiers]}, where len(ref_ids) == RefStats.params['limit']\n
            id_type is from [PMID, DOI, 'PII', 'PUI', 'EMBASE','NCT ID']\n
            [2] references = [ref] list of Reference objects sorted by ETM relevance. len(references) == RefStats.params['limit'] 
            Relevance score is stored in ref['Relevance'] for every reference
        """
        articles, hit_count = self._get_articles(need_snippets=False)

        if self._limit() > 100:
            for page_start in range(len(articles), self._limit(), 100):
                more_articles,_ = self._get_articles(page_start=page_start,need_snippets=False)
                articles += more_articles

        references = list()
        for article in articles:
            relevance_score = float(article['score'])
            if relevance_score >= self.min_relevance:
                etm_ref = ETMjson(article)
                if hasattr(etm_ref,'Identifiers'): #sometime etm_ref is created empty 
                    etm_ref[RELEVANCE] = [relevance_score]
                    references.append(etm_ref)

        ref_ids = defaultdict(list)
        for ref in references:
            id_type, identifier = ref.get_doc_id()
            ref_ids[id_type].append(identifier)

        self.hit_count = hit_count #self.hit_count can be corrupted in parallel ETM requests
        return hit_count, dict(ref_ids), references


    def basic_search(self,search4concepts:list):
        '''
        output tuple:
                [0] hit_count - TOTAL number of reference found by ETM basic search 
                [1] ref_ids = {id_type:[identifiers]}, where len(ref_ids) == ETMsearch.params['limit']\n
                id_type is from [PMID, DOI, 'PII', 'PUI', 'EMBASE','NCT ID']\n
                [2] references = [ref] list of Reference objects sorted by ETM relevance. len(references) == ETMsearch.params['limit'] 
                Relevance score is stored in ref['Relevance'] for every reference
        '''  
        self._set_query(';'.join(search4concepts) )
        return self.__search()
    

    def advanced_query(self,entity1:str,entity2:str,add2query:list):
        '''
        Returns
        -------
        ETM query for advanced serach where "terms" are joined with AND operator
        '''
        self.use_advanced_search()
        query = '{'+'} AND {'.join([entity1,entity2]+add2query)+'}'
        self._set_query(query)
        return self.__search()


    def advanced_query_rel(self,entity1:str,entity2:str,add2query:list):
        '''
        Input
        -----
        Two terms.  Terms can have expansion function specified by /:
        /syn - find synonyms
        /exp - find child terms
        /inf - find inflected form
        /prt - find the term as part of composite term
        '''
        my_terms = list()
        for term in [entity1,entity2]+add2query:
            slash_pos = term.find('/')
            if slash_pos > 0:
                term2add = '{'+term[:slash_pos]+'}'+term[slash_pos:]
            else:
                term2add = '{'+term+'}'
            my_terms.append(term2add)

        query = 'rel({'+'} AND {'.join(my_terms)+'})'
        self.use_advanced_search()
        query = '{'+'} AND {'.join([entity1,entity2]+add2query)+'}'
        self._set_query(query)
        return self.__search()



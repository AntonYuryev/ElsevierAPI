from .references import DocMine, Reference, len, pubmed_hyperlink, make_hyperlink
from .references import AUTHORS,INSTITUTIONS,JOURNAL,PUBYEAR,SENTENCE,EMAIL,REF_ID_TYPES,RELEVANCE, ETM_CITATION_INDEX
import urllib.request
import urllib.parse
from urllib.error import HTTPError
import json, http.client, time, xlsxwriter, os, math
from datetime import datetime,timedelta
from ..ScopusAPI.scopus import Scopus, AuthorSearch
from ..pandas.panda_tricks import pd, df
from concurrent.futures import ThreadPoolExecutor


SCOPUS_AUTHORIDS = 'scopusAuthors'
ETM_REFS_COLUMN = 'Number of references. Link opens most relevant articles in PubMed'
IDENTIFIER_COLUMN = 'Document identifier: PMID or DOI'
MAX_ETM_SESSIONS = 5 #ETM performance deteriorates if number of concurrent sessions is too big
DEFAULT_ETM = 'https://covid19-services.elseviertextmining.com/api'
DEFAULT_APICONFIG = 'D:/Python/ENTELLECT_API/ElsevierAPI/APIconfig.json'


def execution_time(execution_start):
    return "{}".format(str(timedelta(seconds=time.time() - execution_start)))


def load_api_config(api_config_file=''):# file with your API keys and API URLs
    #cur_dir = os.getcwd()
    default_apiconfig = DEFAULT_APICONFIG
    if not api_config_file:
        print('No API config file was specified\nWill use default %s instead'% default_apiconfig)
        api_config_file = default_apiconfig
    try:
        return dict(json.load(open(api_config_file)))
    except FileNotFoundError:
        print("Cannot find API config file: %s" % api_config_file)
        if api_config_file != default_apiconfig:
            print('Cannot open %s config file\nWill use default %s instead'% (api_config_file, default_apiconfig))
            return dict(json.load(open(default_apiconfig)))
        else:
            print('No working API server was specified!!! Goodbye')
            return dict({})
        
class ETMjson(DocMine):

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
                exist_id = article_ids[id_type]
                if isinstance (exist_id, str) : exist_id = [exist_id]
                exist_id.append(id_value)
                article_ids[id_type] = exist_id
                is_article = False #record has multiple IDs of the same type.  It is likely conference proceedings
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
    def __parse_institution(correspondence:str):
        email_pattern = EMAIL.search(correspondence)
        address = correspondence
        if email_pattern:
            email_pos = email_pattern.start()
            email = correspondence[email_pattern.start():email_pattern.end()]
            for i in range(3, email_pattern.start()-3):
                sub_start = email_pos-i
                sub = correspondence[sub_start : sub_start+3]
                if sub in {'. e', '. E'}:
                    email_pos = sub_start
                    break

            address = correspondence[0:email_pos].strip(' .,')
        else:
            email = ''
        
        names = list()
        if address:
            if (not address[0].isalpha()) or address[0].islower():  address = address[1:].strip()
            names = address.split(', ')    

        last_name_index = 0
        has_institutional_keyword = False
        for n in range(1,len(names)):
            if DocMine.has_institution_keyword(names[n]):
                has_institutional_keyword = True
                last_name_index = n
                names[n] = names[n].title()
            else: break

        if not has_institutional_keyword and last_name_index > 0:
            print('article coresspondence %s has no institutional keyword!' % correspondence)

        institutional_names =  names[:last_name_index+1] if names else []

        return institutional_names, address, email


    @staticmethod
    def __parse_contributors(article:dict, scopus_api=None):
        try:
            author_info = article['article']['front']['article-meta']['contribGroupOrAffOrAffAlternatives']
            if not author_info:
                #print('Article has no author info') 
                return [],[],[]
        except KeyError:
            #print('Article has no author info')
            return [],[],[]
            
        if isinstance(author_info, dict): author_info = [author_info]

        authors = list()
        institutions = list()
        scopus_infos = list()
        org_names = list()

        for item in author_info:
            try:
                instituts = item['aff']['content']
                if not isinstance(instituts, list): instituts = [instituts]
                for inst in instituts:
                    org_names, address, email = ETMjson.__parse_institution(inst['institution'])
                    if org_names:
                        institutions.append((org_names[-1],address)) #use only last name of institution defining global organization
                continue
            except KeyError:
                contribOrAddressOrAff_list = item['contrib-group']['contribOrAddressOrAff']
                if not isinstance(contribOrAddressOrAff_list, list): contribOrAddressOrAff_list = [contribOrAddressOrAff_list]

                au_name = ''
                for author in contribOrAddressOrAff_list:
                    au1stname = str()
                    try:
                        au1stname = author['contrib']['contribIdOrAnonymousOrCollab']['name']['given-names']['content']
                    except KeyError:
                        try:
                            au1stname = author['contrib']['contribIdOrAnonymousOrCollab']['name']['given-names']['initials']+'.'
                        except KeyError: pass

                    try:
                        au_surname = author['contrib']['contribIdOrAnonymousOrCollab']['name']['surname']
                        au_name = au_surname + ' ' + au1stname if au1stname else au_surname
                        authors.append(str(au_name).title())

                        if isinstance(scopus_api,AuthorSearch):
                            scopus_info = scopus_api.get_author_id(au_surname,au1stname,org_names)
                        else:
                            scopus_info = []

                        if scopus_info:
                            scopus_infos.append(scopus_info)

                    except KeyError: continue

        return authors, institutions, scopus_infos #authors: list(str), institutions: list([orname,address])


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
                    if query_end <= last_markup_end: continue

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


    def __init__(self, article:dict,scopus_api=None): 
        is_article, article_ids = self.__parse_ids(article)
        if isinstance(scopus_api,Scopus): 
            self.scopusInfo = dict()
        if not is_article:
            return None # __init__ must return None not dict()
        else:
            id_type, id_value = article_ids.popitem()
            if id_type == 'ELSEVIER': id_type = 'PII'
            super().__init__(id_type, id_value)
            self.Identifiers.update(article_ids)

        try:
            title = article['article']['front']['article-meta']['title-group']['article-title']
            self._set_title(title)
        except KeyError: pass

        try:
            year = article['article']['front']['article-meta']['pub-date']['year']
            self.set_date(year)
        except KeyError: pass

        self.add2section('Abstract',' '.join(self.__parse_abstact(article)))
        authors,institutions,scopus_infos = self.__parse_contributors(article,scopus_api)
        if authors: self[AUTHORS] = authors
        self.addresses = {i[0]:i[1] for i in institutions}
        if scopus_infos: 
            self[SCOPUS_AUTHORIDS] = [i[0] for i in scopus_infos]
            self.scopusInfo.update({i[0]:i for i in scopus_infos})
        try:
            journal = article['article']['front']['journal-meta']['journal-title-group']['journal-title']
            journal = ','.join(journal.values())
        except KeyError:
            journal = 'Grant application'
        self.append_property(JOURNAL, journal)

        self['categories'] = self.__parse_category(article)

        text_ref = self._make_textref()
        snippet, self.term_ids, self.keywords = self.__parse_snippet(article) #term_ids = {term_id+'\t'+term_name}
        self.add_sentence_prop(text_ref, SENTENCE, snippet)


class ETMstat:
    """
    self.ref_counter = {str(id_type+':'+identifier):(ref,count)}
    """
    url = DEFAULT_ETM
    page_size = 100
    request_type = '/search/basic?'  # '/search/advanced?'
    def __base_url(self): return self.url+self.request_type
    min_relevance = float(0.0)
    

    def __init__(self,APIconfig:dict=dict({}), limit=5, add_param=dict()):
        """
        self.ref_counter = {str(id_type+':'+identifier):(ref,count)}
        """
        myAPIconfig = APIconfig if APIconfig else load_api_config()

        self.url = myAPIconfig['ETMURL']
        self.params = {
            'searchTarget': 'full_index',
            'apikey':myAPIconfig['ETMapikey'],
            'limit': limit,
            'so_c':'17~' #  filter by Publication Source excluding NIH reporter
            #'so_p':'ffb~' # filter by Publication Type
            }
        self.params.update(add_param)       
        self.ref_counter = dict() # {str(id_type+':'+identifier):(ref,count)}
        self.etm_ref_column_name = list()
        self.etm_doi_column_name = list()
        self.hit_count = 0


    def clone(self, to_url=''):
        if not to_url:
            to_url = self.url
        myApiconfig = {'ETMURL':to_url,'ETMapikey':self.params['apikey']}
        newEtMstat =  ETMstat(myApiconfig,self._limit(),self.params)
        newEtMstat.page_size = self.page_size
        newEtMstat.request_type = self.request_type
        newEtMstat.hit_count = 0
        newEtMstat.min_relevance = self.min_relevance
        return newEtMstat


    def _limit(self): 
        """
        Returns
        -------
        self.params['limit'] - max number of reference to retreive
        """
        return self.params['limit']

    def references(self):
        return set([x[0] for x in self.ref_counter.values()])
 
    def _set_query(self,query:str):
        self.params['query'] = query

    def _query(self):
        try:
            return self.params['query']
        except KeyError:
            print('ETMstat class instance has no query specified')

    def __get_param_str(self):
        return urllib.parse.urlencode(self.params)

    def _url_request(self):
        return self.__base_url()+self.__get_param_str()

    def search_reviews_only(self):
        self.params.update({'so_p': '1~'})

    def use_advanced_search(self):
        self.request_type = '/search/advanced?'


    def _get_articles(self, page_start=0, need_snippets=True):
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
                the_page = urllib.request.urlopen(self._url_request()).read()
                if the_page:
                    result = json.loads(the_page.decode('utf-8'))
                    if attempt > 1:
                        print(f'ETM connection was restored on the {attempt} attempt')
                    return list(result['article-data']), int(result['total-hits='])
                else:
                    return list(),int(0)
            except HTTPError:
                tiout = 10
                print(f'HTTPError. Will attempt to reconnect in {tiout}')
                time.sleep(tiout)
            except http.client.IncompleteRead:
                tiout = 10
                print(f'IncompleteRead. Will attempt to reconnect in {tiout}')
                time.sleep(tiout)

        print(f'Cannot connect to {self.__base_url()} after 10 attempts')
        return list(),int(0)


    def _add2counter(self, ref:Reference):
        for id_type in REF_ID_TYPES:
            try:
                identifier = ref.Identifiers[id_type]
                counter_key = id_type+':'+identifier
                try:
                    count_exist = self.ref_counter[counter_key][1]
                    self.ref_counter[counter_key] = (ref, count_exist+1)
                    self.ref_counter[counter_key][0][RELEVANCE][0] += ref.relevance()
                except KeyError:
                    self.ref_counter[counter_key] = (ref,1)
                return id_type,identifier
            except KeyError:
                continue
        
        return '',''


    def counter2df(self, use_relevance=True):
        """
        used for internal self.ref_counter = {identifier:(ref,ref_count)}
        used to count references from ETM
        """
        table_rows = set()
        for identifier, ref_count in self.ref_counter.items():
            ref = ref_count[0]
            refcount = ref_count[1]
            score = math.ceil(float(ref[RELEVANCE][0]) * refcount) if use_relevance else refcount
            biblio_str, id_type, identifier = ref._biblio_tuple()
            if id_type == 'PMID':
                identifier = pubmed_hyperlink([identifier],identifier)
            elif id_type == 'DOI':
                identifier = make_hyperlink(identifier,'http://dx.doi.org/')
            table_rows.add(tuple([score,biblio_str,id_type,identifier]))

        first_col = RELEVANCE if use_relevance else ETM_CITATION_INDEX
        header = [first_col,'Citation','Identifier type',IDENTIFIER_COLUMN]
        return_pd = df.from_rows(list(table_rows),header)
        if use_relevance:
            return_pd[RELEVANCE] = return_pd[RELEVANCE].astype(float).round(2)
        else:
            return_pd[ETM_CITATION_INDEX] = return_pd[ETM_CITATION_INDEX].astype(int)
        return_pd.sort_values(first_col,ascending=False,inplace=True)
        return_pd.add_column_format('Citation','width',150)
        return_pd.add_column_format('Citation','wrap_text',True)
        return_pd.make_header_horizontal()
        return_pd.set_hyperlink_color([IDENTIFIER_COLUMN])
        return return_pd


    @staticmethod
    def external_counter2pd(ref_counter:set, stat_prop=RELEVANCE):
        '''
        use for external ref_counter = {Reference} where Reference objects are annotated as Reference[stat_prop]
        used to count references in ResnetGraph()
        '''
        table_rows = set()
        for ref in ref_counter:
            biblio_str, id_type, identifier = ref._biblio_tuple()
            if id_type == 'PMID':
                identifier = pubmed_hyperlink([identifier],identifier)
            elif id_type == 'DOI':
                identifier = make_hyperlink(identifier,'http://dx.doi.org/')
            elif id_type == 'NCT ID':
                identifier = make_hyperlink(identifier,'https://clinicaltrials.gov/ct2/show/')
            table_rows.add(tuple([ref[stat_prop][0],biblio_str,id_type,identifier]))

        header = [stat_prop,'Citation','Identifier type',IDENTIFIER_COLUMN]
        return_pd = df.from_rows(list(table_rows),header)
        if stat_prop == RELEVANCE:
            return_pd[RELEVANCE] = return_pd[RELEVANCE].astype(float).round(2)
        else:
            return_pd[stat_prop] = return_pd[stat_prop].astype(int)
        return_pd.sort_values(stat_prop,ascending=False,inplace=True)
        return_pd.add_column_format('Citation','width',150)
        return_pd.add_column_format('Citation','wrap_text',True)
        return_pd.make_header_horizontal()
        return_pd.set_hyperlink_color([IDENTIFIER_COLUMN])
        return return_pd

   
    @staticmethod
    def basic_query(entity1:str,entity2:str,add2query:list):
        '''
        Returns
        -------
        ETM query for  basic serach where "terms" are joined with ";"
        '''
        return ';'.join([entity1,entity2]+add2query)
        #return '{'+'};{'.join([entity1,entity2]+add2query)+'}'


    @staticmethod
    def advanced_query(entity1:str,entity2:str,add2query:list):
        '''
        Returns
        -------
        ETM query for advanced serach where "terms" are joined with AND operator
        '''
        return '{'+'} AND {'.join([entity1,entity2]+add2query)+'}'


    @staticmethod
    def advanced_query_rel(entity1:str,entity2:str,add2query:list):
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

        return 'rel({'+'} AND {'.join(my_terms)+'})'
   

    def __get_stats(self):
        """
        Returns
        -------
        Tuple:
            [0] hit_count - TOTAL number of reference found by ETM basic search 
            [1] ref_ids = {id_type:[identifiers]}, where len(ref_ids) == ETMstat.params['limit']\n
            id_type is from [PMID, DOI, 'PII', 'PUI', 'EMBASE','NCT ID']\n
            [2] references = [ref] list of Reference objects sorted by ETM relevance. len(references) == ETMstat.params['limit'] 
            Relevance score is stored in ref['Relevance'] for every reference
        """
        articles, hit_count = self._get_articles(need_snippets=False)

        if self._limit() > 100:
            for page_start in range(len(articles), self._limit(), 100):
                more_articles, discard = self._get_articles(page_start=page_start,need_snippets=False)
                articles += more_articles

        references = list()
        for article in articles:
            relevance_score = float(article['score'])
            if relevance_score >= self.min_relevance:
                etm_ref = ETMjson(article)
                etm_ref[RELEVANCE] = [relevance_score]
                if hasattr(etm_ref,"Identifiers"):
                    references.append(etm_ref)

        ref_ids = dict()
        for ref in references:
            id_type, identifier = self._add2counter(ref)
            try:
                ref_ids[id_type].append(identifier)
            except KeyError:
                ref_ids[id_type] = [identifier]

        self.hit_count = hit_count #self.hit_count can be corrupted in parallel ETM requests
        return hit_count, ref_ids, references


    @staticmethod
    def count_refs(ref_counter:set, references:list):
        '''
        updates "ref_counter" with "references"
        updates 'Citation index' for reference in "ref_counter" by 1 if reference is in "references"
        '''
        ref_counter.update(references)
        counter_refs = ref_counter.intersection(set(references))
        for ref in counter_refs:
            try:
                count = ref['Citation index'][0]
                ref['Citation index'] = [count + 1]
            except KeyError:
                ref['Citation index'] = [1]
        return


    @staticmethod
    def refcounter2tsv(fname:str, ref_counter:set, use_relevance=False,include_idtype=False):
        to_sort = list(ref_counter)
        if use_relevance:
            for ref in to_sort:
                try:
                    relevance_index = float(ref[ETM_CITATION_INDEX][0])*float(ref[RELEVANCE][0])
                except KeyError: relevance_index = 0.0
                ref['Relevance index'] = [float(relevance_index)]
            to_sort.sort(key=lambda x: x['Relevance index'][0], reverse=True)
            with open(fname, 'w', encoding='utf-8') as f:
                f.write('Relevance\tCitation\tPMID or DOI\n')
                for ref in to_sort:
                    biblio_tup = ref._biblio_tuple()
                    if not biblio_tup[2]: 
                        continue
                    if include_idtype:
                        biblio_str = biblio_tup[0]+'\t'+biblio_tup[1]+':'+biblio_tup[2]
                    else:
                        biblio_str = biblio_tup[0]+'\t'+biblio_tup[2]

                    f.write(str(ref['Relevance index'][0])+'\t'+biblio_str+'\n')
        else:        
            to_sort.sort(key=lambda x: x['Citation index'][0], reverse=True)
            with open(fname, 'w', encoding='utf-8') as f:
                f.write('Citation index\tCitation\tPMID or DOI\n')
                for ref in to_sort:
                    biblio_tup = ref._biblio_tuple()
                    if not biblio_tup[2]: 
                        continue #to remove reference with no ID type
                    biblio_tup = ref._biblio_tuple()
                    if include_idtype:
                        biblio_str = biblio_tup[0]+'\t'+biblio_tup[1]+':'+biblio_tup[2]
                    else:
                        biblio_str = biblio_tup[0]+'\t'+biblio_tup[2]
                    f.write(str(ref['Citation index'][0])+'\t'+biblio_str+'\n')


    @staticmethod
    def _sort_dict(dic:dict, sort_by_value = True, reverse = True):
        item_idx = 1 if sort_by_value else 0
        return dict(sorted(dic.items(), key=lambda item: item[item_idx],reverse=reverse))  


    def basic_search(self,search4concepts:list):
        '''
        Return
        ------
        Tuple:
            [0] hit_count - TOTAL number of reference found by ETM basic search 
            [1] ref_ids = {id_type:[identifiers]}, where len(ref_ids) == ETMstat.params['limit']\n
            id_type is from [PMID, DOI, 'PII', 'PUI', 'EMBASE','NCT ID']\n
            [2] references = [ref] list of Reference objects sorted by ETM relevance. len(references) == ETMstat.params['limit'] 
            Relevance score is stored in ref['Relevance'] for every reference
        '''
        query = ';'.join(search4concepts)       
        self._set_query(query)
        return self.__get_stats()
    
    
    def __multiple_search(self,for_entity:str,and_concepts:list,my_query,add2query=[]):
        references = set()
        total_hits = 0
        for concept in and_concepts:
            query = my_query(for_entity,concept,add2query)
            self._set_query(query)
            hit_count,ref_ids,etm_refs = self.__get_stats()
            total_hits += hit_count
            references.update(etm_refs)

        return references, total_hits
    

    def __get_refs(self, entity_name:str, concepts2link:list,my_query,add2query=[]):
        references = set()
  #      if entity_name == 'siltuximab':
  #          print('')
        references, total_hits = self.__multiple_search(entity_name,concepts2link,my_query,add2query)

        references = list(references)
        references.sort(key=lambda x:x[RELEVANCE][0],reverse=True)
        best_refid_list = references[0:self._limit()]

        pmids=list()
        dois = list()
        for ref in best_refid_list:
            id_type,identifier = ref.get_doc_id()
            if id_type == 'PMID':
                pmids.append(identifier)
            elif id_type == 'DOI':
                dois.append(identifier)
            else:
                continue

        hyperlink2pubmed = pubmed_hyperlink(pmids,str(total_hits)) if pmids else '' 
        doi_str = make_hyperlink(dois[0],url='http://dx.doi.org/', display_str=';'.join(dois)) if dois else ''
        return [hyperlink2pubmed,doi_str]


    @staticmethod
    def _etm_ref_column_name(between_column:str, and_concepts:str|list):
        if isinstance(and_concepts,str):
            return ETM_REFS_COLUMN + ' between '+between_column+' and '+and_concepts
        else:
            return ETM_REFS_COLUMN + ' between '+between_column+' and '+','.join(and_concepts)

    @staticmethod
    def _etm_doi_column_name(between_column:str, and_concepts:str|list):
        if isinstance(and_concepts,str):
            return 'DOIs' + ' between '+between_column+' and '+and_concepts
        else:
            return 'DOIs' + ' between '+between_column+' and '+','.join(and_concepts)


    def __add_refs1(self,to_df:df,between_names_in_col:str,and_concepts:list,my_query,add2query=[]):
        my_df = df.copy_df(to_df)
        etm_ref_column_name = self._etm_ref_column_name(between_names_in_col,and_concepts)
        doi_ref_column_name = self._etm_doi_column_name(between_names_in_col,and_concepts)
        my_df[[etm_ref_column_name,doi_ref_column_name]] = my_df[between_names_in_col].apply(lambda row: self.__get_refs(row,and_concepts,my_query,add2query)).apply(pd.Series)
        return my_df,self.ref_counter # need to return self.ref_counter to concatenate ref_counters during cloning


    def add_etm_references(self,to_df:df,between_names_in_col:str,and_concepts:list,use_query,add2query=[],max_row=0):
        """
        Input
        -----
        my_query - function to generate query from "between_names_in_col" and each concept in "and_concepts"
        my_query must have 3 arguments: my_query(entity1:str, entity2:str, add2query:list)\n
        where add2query - list of additinal keywords used for all pairs "between_names_in_col" and "and_concepts"

        Returns
        ----
        copy of input to_df with added columns ETM_REFS_COLUMN,"DOIs"\n
        self.etm_ref_column_name - []\n
        self.etm_doi_column_name - []\n
        """
        if to_df.empty: return to_df
        start_time = time.time()
        if max_row:
            df2annotate = df.from_pd(to_df.iloc[:max_row])
            unannoated_rows = df.from_pd(to_df.iloc[max_row:])
        else:
            df2annotate = df.copy_df(to_df)
            unannoated_rows = df()
        #my_df = to_df.copy_df(to_df)
        row_count = len(to_df)
        
        if row_count > 9:
            thread_name = f'ETM refs 4 {len(df2annotate)} rows in {MAX_ETM_SESSIONS} threads'
            partition_size = 100
            with ThreadPoolExecutor(max_workers=MAX_ETM_SESSIONS, thread_name_prefix=thread_name) as e:
                futures = list()
                for i in range(0,len(df2annotate),partition_size):
                    df_part = df.from_pd(df2annotate.iloc[i:i+partition_size])
                    new_session = self.clone() # need to clone here to avoid self.ref_counter mutation
                    futures.append(e.submit(new_session.__add_refs1,df_part,between_names_in_col,and_concepts,use_query,add2query))

                dfs2concat = list()
                for f in futures: # cannot use as_completed here to ensure the same order for concatenation
                    session_df, session_refcounter = f.result()
                    dfs2concat.append(session_df)
                    [self._add2counter(t[0]) for t in session_refcounter.values()] # combine references from all futures
                
                dfs2concat.append(unannoated_rows)
                annotated_df = df.concat_df(dfs2concat,to_df._name_)
        else:
            annotated_df = self. __add_refs1(to_df,between_names_in_col,and_concepts,use_query,add2query)[0]
        
        annotated_df.copy_format(to_df)
        etm_ref_column_name = self._etm_ref_column_name(between_names_in_col,and_concepts)
        self.etm_ref_column_name.append(etm_ref_column_name)
        annotated_df.add_column_format(etm_ref_column_name,'align','center')
        etm_doi_column_name = self._etm_doi_column_name(between_names_in_col,and_concepts)
        self.etm_doi_column_name.append(etm_ref_column_name)
        annotated_df.set_hyperlink_color([etm_ref_column_name,etm_doi_column_name])
        
        rows_annotated = max_row if max_row else len(to_df)
        print('Annotated with ETM references %d rows from %s in %s' % 
                (rows_annotated,to_df._name_,execution_time(start_time)))
        return annotated_df

    '''
    def add_etm_refsOLD(self,to_df:df,between_names_in_col:str,and_concepts:list,use_query,add2query=[]):
        """
        Input
        -----
        my_query - function to generate query from "between_names_in_col" and each concept in "and_concepts"
        my_query must have 3 arguments: my_query(entity1:str, entity2:str, add2query:list)\n
        where add2query - list of additinal keywords used for all pairs "between_names_in_col" and "and_concepts"

        Returns
        ----
        copy of input to_df with added columns "ETM_REFS_COLUMN","DOIs"
        """
        start_time = time.time()
        etm1 = self.clone(DEFAULT_ETM)
        etm2 = self.clone('https://discover.elseviertextmining.com/api')
        etm3 = self.clone('https://research.elseviertextmining.com/api')
        my_df = to_df.copy_df(to_df)
        
        row_count = len(my_df)
        if row_count > 9:
            one3rd = int(row_count/3)
            df1 = df(my_df.iloc[:one3rd])
            df2 = df(my_df.iloc[one3rd : 2*one3rd])
            df3 = df(my_df.iloc[2*one3rd:])

            t1 = Thread(target=etm1. __add_refs1, args=(df1,between_names_in_col,and_concepts,use_query,add2query),name='etm_demo')
            t1.start()
            t2 = Thread(target=etm2. __add_refs1, args=(df2,between_names_in_col,and_concepts,use_query,add2query),name='etm_discover')
            t2.start()
            t3 = Thread(target=etm3. __add_refs1, args=(df3,between_names_in_col,and_concepts,use_query,add2query),name='etm_research')
            t3.start()

            t1.join()
            t2.join()
            t3.join()

            annotated_df = df(pd.concat([df1,df2,df3]),name=to_df._name_)
            [self._add2counter(ref) for ref in etm1.references()]
            [self._add2counter(ref) for ref in etm2.references()]
            [self._add2counter(ref) for ref in etm3.references()]
        else:
            self. __add_refs1(my_df,between_names_in_col,and_concepts,use_query,add2query)
            annotated_df = my_df

        annotated_df.copy_format(to_df)
        etm_ref_column_name = self._etm_ref_column_name(between_names_in_col,and_concepts)
        self.etm_ref_column_name.append(etm_ref_column_name)
        annotated_df.add_column_format(etm_ref_column_name,'align','center')
        etm_doi_column_name = self._etm_doi_column_name(between_names_in_col,and_concepts)
        self.etm_doi_column_name.append(etm_ref_column_name)
        annotated_df.set_hyperlink_color([etm_ref_column_name,etm_doi_column_name])
        print('Annotated %d rows from %s with ETM references in %s' % 
                (len(to_df),to_df._name_,execution_time(start_time)))
        return annotated_df
    '''


    def __etm42columns(self,in_df:df,between_col:str,and_col:str,my_query,add2query=[]):
        my_df = df.copy_df(in_df)
        etm_ref_column_name = self._etm_ref_column_name(between_col,and_col)
        doi_ref_column_name = self._etm_doi_column_name(between_col,and_col)
        for i in in_df.index:
            col1 = my_df.loc[i][between_col]
            col2 = my_df.loc[i][and_col]
            hyperlink2pubmed,hyperlink2doi = self.__get_refs(col1,[col2],my_query,add2query)
            my_df.at[i,etm_ref_column_name] = hyperlink2pubmed
            my_df.at[i,doi_ref_column_name] = hyperlink2doi

        return my_df,self.ref_counter


    def add_etm42columns(self,in_df:df,between_col:str,and_col:str,my_query,add2query=[],max_row=0):
        start_time = time.time()
        if max_row:
            df2annotate = df.from_pd(in_df.iloc[:max_row])
            unannoated_rows = df(in_df.iloc[max_row:])
        else:
            df2annotate = df.copy_df(in_df)
            unannoated_rows = df()

        partition_size = 100
        thread_name = f'Retrieve {partition_size} ETM references'
        dfs2concat = list()
        with ThreadPoolExecutor(max_workers=MAX_ETM_SESSIONS, thread_name_prefix=thread_name) as e:    
            futures = list()
            for i in range(0,len(df2annotate),partition_size):
                df_part = df(df2annotate.iloc[i:i+partition_size])
                new_session = self.clone() # need to clone here to avoid self.ref_counter mutation
                futures.append(e.submit(new_session.__etm42columns,df_part,between_col,and_col,my_query,add2query))

            for f in futures: # cannot use as_completed here to ensure the same order for concatenation
                session_df, session_refcounter = f.result()
                dfs2concat.append(session_df)
                [self._add2counter(t[0]) for t in session_refcounter.values()] # combine references from all futures 

            dfs2concat.append(unannoated_rows)

        annotated_df = df.concat_df(dfs2concat,in_df._name_)
        
        etm_ref_column_name = self._etm_ref_column_name(between_col,and_col)
        doi_ref_column_name = self._etm_doi_column_name(between_col,and_col)
        self.etm_ref_column_name.append(etm_ref_column_name)
        self.etm_doi_column_name.append(doi_ref_column_name)

        annotated_df.set_hyperlink_color([etm_ref_column_name,doi_ref_column_name])
        annotated_df.add_column_format(etm_ref_column_name,'align','center')
        annotated_df.set_hyperlink_color([etm_ref_column_name,doi_ref_column_name])
        print('Annotated %d rows from %s with ETM references in %s' % 
                (len(in_df),in_df._name_,execution_time(start_time)))
        return annotated_df


ETM_CACHE_DIR = 'ElsevierAPI/ETM_API/__etmcache__'
class ETMcache (ETMstat):
    """
    self.statistics = {prop_name:{prop_value:count}}
    """
    pass
    request_type = '/search/advanced?'
    etm_results_dir = ETM_CACHE_DIR
    etm_stat_dir = ETM_CACHE_DIR
    
    def __init__(self,query,search_name,APIconfig=dict(),etm_dump_dir='',etm_stat_dir='',add_params={}):
        super().__init__(APIconfig,add_param=add_params)
        self.params.pop('limit')
        self._set_query(query)
        self.hit_count = self._get_articles()[1]
        self.search_name = search_name
        self.etm_results_dir = etm_dump_dir if etm_dump_dir else ETM_CACHE_DIR
        self.etm_stat_dir = etm_stat_dir if etm_stat_dir else ETM_CACHE_DIR
        self.statistics = dict() # {prop_name:{prop_value:count}}
        self.scopusAPI = None


    def set_query(self,query:str,search_name:str):
        super()._set_query(query)
        self.search_name = search_name


    def set_stat_props(self,stat_props:list):
        for p in stat_props:
            if p == AUTHORS: p = SCOPUS_AUTHORIDS
            self.statistics[p] = dict()
        

    def __download(self,cached_articles=[],dump_date_create=''):
        print('\nPerforming ETM search "%s"' % self.search_name)
        articles, self.hit_count = self._get_articles(need_snippets=False)
        print (f'Query: "{self._query()}" found {self.hit_count} results' )
        if  self.hit_count > len(cached_articles): # self.hit_count 1-based index !!!
            print('Today ETM search finds %d results. %d more than in cache' %(self.hit_count, self.hit_count-len(articles)))
            start = time.time()
            if dump_date_create:
                self.params.update({'so_d':dump_date_create}) # format: 2003-03-19
            for page_start in range(0,self.hit_count,self.page_size):
                more_articles = self._get_articles(page_start=page_start)[0]
                articles += more_articles
                download_count = len(articles)
                print("Downloaded %d out of %d hits in %s" % (download_count,self.hit_count, execution_time(start)))
            self.params.pop('so_d','')
            return articles
        else:
            print('No new articles found. Will use articles from cache')
            return []


    def __dumps(self, articles:list, into_file:str):
        try:
            dump = open(into_file, "w", encoding='utf-8')
            dump.write(json.dumps(articles,indent=1))
        except FileNotFoundError:
            # path was specified incorrectly, dumps into the current dir
            dump = open(self.search_name+'.json', "w", encoding='utf-8')
            dump.write(json.dumps(articles,indent=1))
 

    def load_from_json(self,use_cache=True):
        if use_cache:
            dir_name = self.etm_results_dir
            if dir_name[-1] != '/':
                dir_name += '/'
            dumpfile_name = dir_name+self.search_name+'.json'
            try:
                # attepts to find json file with ETM results saved after ETM API call below
                f = open(dumpfile_name,'r')
                articles = json.load(f)
                f.close()

                path_end = dumpfile_name.rfind('/')
                if path_end < 0:
                    path_end = dumpfile_name.rfind('\\')
                fname = dumpfile_name[path_end+1:]
                cache_path = dumpfile_name[:path_end]
                print('Found "%s" file in "%s" cache with %d articles. Will use it to load results' % (fname,cache_path,len(articles))) 
                d =  datetime.strptime(time.ctime(os.path.getctime(dumpfile_name)), "%c")
                dump_create_date =  d.strftime('%Y-%m-%d')
                update_articles = self.__download(articles,dump_create_date)
                articles = articles + update_articles
                self.__dumps(articles,dumpfile_name)
            except FileNotFoundError:
                #if json dump file is not found new ETM search is initiated
                articles = self.__download()
                self.__dumps(articles,dumpfile_name)
        else:
            articles = self.__download()

        for article in articles:
            etm_ref = ETMjson(article,self.scopusAPI)
            if etm_ref: # etm_ref can be empty if it is conference proceedings
                relevance_score = float(article['score'])
                etm_ref[RELEVANCE] = [relevance_score]
                self._add2counter(etm_ref)

        
    def get_statistics(self, stat_prop_list):
        self.term2refs = dict() # {term:{ref}}
        self.keywords2ref = dict() #{keyword:{ref}}
        self.scopusid2name = dict()
        references = self.references()
        for ref in references:
            for p in stat_prop_list:
                ref.count_property(self.statistics[p], p)

            if hasattr(ref,'term_ids'):
                for term in ref.term_ids: #term_ids = {term_id+'\t'+term_name}
                    try:
                        self.term2refs[term].add(ref)
                    except KeyError:
                        self.term2refs[term] = {ref}
            
            if hasattr(ref,'keywords'):
                for k in ref.keywords:
                    try:
                        self.term2refs[k].add(ref)
                    except KeyError:
                        self.term2refs[k] = {ref}

            if hasattr(ref,'scopusInfo'):
                self.scopusid2name.update({k:v[1]+' '+v[2] for k,v in dict(ref.scopusInfo).items()})


    def _org2address_stats(self):
        org2addres = dict()
        references = self.references()
        for ref in references:
            for name, address in ref.addresses.items():
                try:
                    exist_address = org2addres[name]
                    if len(exist_address) < len(address):
                        org2addres[name] = address
                except KeyError:
                    org2addres[name] = address
        
        try:
            institution_counter = dict(self.statistics[INSTITUTIONS])
            org_address_couner = {org2addres[k]:v for k,v in institution_counter.items()}
            return dict(sorted(org_address_couner.items(), key=lambda item: item[1],reverse=True))
        except KeyError: return dict()


    def to_excel(self, stat_props:list):
        if self.scopusid2name:
            replaceid4name = dict()
            scopus_stats = dict(self.statistics[SCOPUS_AUTHORIDS])
            for k,v in scopus_stats.items():
                au_name = self.scopusid2name[k]
                replaceid4name[au_name] = v
            
            self.statistics[SCOPUS_AUTHORIDS] = replaceid4name
                
        workbook = xlsxwriter.Workbook(self.etm_stat_dir+self.search_name+'.xlsx')
        for p in stat_props:
            try:
                dic = dict(self.statistics[p])
                sort_by_value = True if p != PUBYEAR else False
                prop_stat = self._org2address_stats() if p == INSTITUTIONS else self._sort_dict(dic,sort_by_value)
                DocMine.dict2worksheet(workbook,p,prop_stat, [p, '#References'],100)
            except KeyError: continue

        worksheet_ref = workbook.add_worksheet('References')
        worksheet_ref.write_string(0,0,'Total number of articles: '+str(len(self.references())))
        row_counter = 1
        for ref in self.references():
            row_str = ref.to_str(['PMID','DOI','PUI'],col_sep='\n')+'\n\t\t\t\n'
            worksheet_ref.write_string(row_counter,0,row_str)
            row_counter += 1

        workbook.close()


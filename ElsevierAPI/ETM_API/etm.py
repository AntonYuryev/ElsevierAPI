from .references import DocMine, Reference
from .references import AUTHORS,INSTITUTIONS,JOURNAL,PUBYEAR,SENTENCE,EMAIL,REF_ID_TYPES
import urllib.request
import urllib.parse
import json
import time
from datetime import timedelta
import xlsxwriter
import os
import math

def execution_time(execution_start):
    return "{}".format(str(timedelta(seconds=time.time() - execution_start)))

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
                    return ''
        
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
                    for s in scs:
                        abstract_sentences.append(str(s['addressOrAlternativesOrArray']['p']))
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
    def __parse_contributors(article:dict):
        try:
            author_info = article['article']['front']['article-meta']['contribGroupOrAffOrAffAlternatives']
            if not author_info:
                #print('Article has no author info') 
                return [],[]
        except KeyError:
            #print('Article has no author info')
            return [],[]
            
        if isinstance(author_info, dict): author_info = [author_info]

        authors = list()
        institutions = list()
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
                    try:
                        au_name = author['contrib']['contribIdOrAnonymousOrCollab']['name']['given-names']['content']
                    except KeyError:
                        try:
                            au_name = author['contrib']['contribIdOrAnonymousOrCollab']['name']['given-names']['initials']+'.'
                        except KeyError: pass

                    try:
                        au_surname = author['contrib']['contribIdOrAnonymousOrCollab']['name']['surname']
                        au_name = au_surname + ' ' + au_name  
                        authors.append(str(au_name).title())
                    except KeyError: continue

        return authors, institutions #authors: list(str), institutions: list([orname,address])

    @staticmethod
    def __parse_snippet(article:dict, markup_all=False):
        try:
            snippet = str(article['data']['snippet']['text'])
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
                return article['article']['body']['sec']['addressOrAlternativesOrArray']['p'], set(), set()
            except KeyError:
                return '', set(), set()


    def __init__(self, article:dict): 
        is_article, article_ids = self.__parse_ids(article)
        if not is_article: 
            return None
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

        self.add2section('Abstract',self.__parse_abstact(article))
        authors,institutions = self.__parse_contributors(article)
        self[AUTHORS] = authors
        self.addresses = {i[0]:i[1] for i in institutions}
        
        try:
            journal = article['article']['front']['journal-meta']['journal-title-group']['journal-title']
            journal = ','.join(journal.values())
        except KeyError:
            journal = 'Grant application'
        self.append_property(JOURNAL, journal)

        text_ref = self._make_textref()
        snippet, self.term_ids, self.keywords = self.__parse_snippet(article) #term_ids = {term_id+'\t'+term_name}
        self.add_sentence_prop(text_ref, SENTENCE, snippet)

class ETMstat:
    url = 'https://demo.elseviertextmining.com/api'
    page_size = 100
    request_type = '/search/basic?'  # '/search/advanced?'
    def __base_url(self): return self.url+self.request_type
    searchTarget = 'full_index'
    number_snippets = 1
    hit_count = 0
    params = dict()
    ref_counter = dict() # {str(id_type+':'+identifier):(ref,count)}

    def __init__(self,APIconfig:dict, limit=5, add_param=dict()):
        self.url = APIconfig['ETMURL']
        self.params = {
            'searchTarget': 'full_index',
            'snip': '1.desc',
            'apikey':APIconfig['ETMapikey'],
            'limit': limit,
            }
        self.params.update(add_param)

    def __limit(self): return self.params['limit']
 
    def _set_query(self,query:str):
        self.params['query'] = query

    def __get_param_str(self):
        return urllib.parse.urlencode(self.params)

    def _url_request(self):
        return self.__base_url()+self.__get_param_str()

    def _get_articles(self, page_start=0):
        if page_start: self.params['start'] = page_start
        the_page = urllib.request.urlopen(self._url_request()).read()
        if the_page:
            result = json.loads(the_page.decode('utf-8'))
            return result['article-data'], result['total-hits=']
        else:
            return [], 0


    def _add2counter(self, ref:Reference):
        for id_type in REF_ID_TYPES:
            try:
                identifier = ref.Identifiers[id_type]
                counter_key = id_type+':'+identifier
                try:
                    count_exist = self.ref_counter[counter_key][1]
                    self.ref_counter[counter_key] = (ref, count_exist+1)  # {str(id_type+':'+identifier):(ref,count)}
                except KeyError:
                    self.ref_counter[counter_key] = (ref,1)
                return identifier
            except KeyError:
                continue


    def print_counter(self, fname, use_relevance=True):
        table_rows = list()
        for identifier, ref_count in self.ref_counter.items():
            ref = ref_count[0]
            refcount = ref_count[1]
            if use_relevance:
                score = math.ceil(float(ref['Relevance rank'][0]) * refcount)
                first_col = 'Relevance'
            else:
                score = refcount
                first_col = 'Citation index'

            biblio_str, identifiers = ref.get_biblio_str()
            table_rows.append((score,biblio_str,identifier))

        table_rows.sort(key=lambda x:x[0],reverse=True)
        with open(fname, 'w', encoding='utf-8') as f:
            f.write(first_col+'\tCitation\tPMID or DOI\n')
            for tup in table_rows:
                f.write(str(tup[0])+'\t'+tup[1]+'\t'+tup[2]+'\n')


    def relevant_articles(self, terms:list):
        # add_param controls number of best references to return. Defaults to 5
        query = '{'+'};{'.join(terms)+'}'
        self._set_query(query)
        articles = list()

        articles, self.hit_count = self._get_articles()

        if self.__limit() > 100:
            for page_start in range(len(articles), self.__limit(), 100):
                more_articles, discard = self._get_articles(page_start=page_start)
                articles += more_articles

        references = list()
        for article in articles:
            etm_ref = ETMjson(article)
            relevance_score = float(article['score'])
            etm_ref['Relevance rank'] = [relevance_score]
            if hasattr(etm_ref,"Identifiers"):
                references.append(etm_ref)

        ref_ids = list()
        for ref in references:
            ref_ids.append(self._add2counter(ref))
        
        #best_refs_str = ';'.join(ref_ids)
        return self.hit_count, ref_ids, references


class ETMcache (ETMstat):
    pass
    #url = 'https://demo.elseviertextmining.com/api'
    #api_key = str()
    request_type = '/search/advanced'
    #page_size = 100
    #def __base_url(self): return self.url+self.request_type+'?'
    #searchTarget = 'full_index'
    #number_snippets = 1
    #hit_count = 0
    search_name = str()
    #query = str()
    #url_request = str()
    #params = dict()
    etm_results_dir = ''
    etm_stat_dir = ''
    statistics = dict() # {prop_name:{prop_value:count}}
    
    def references(self):
        return set([x[0] for x in self.ref_counter.values()])

    def __download(self):
        print('\nPerforming ETM search "%s"' % self.search_name)
        print ('Query: %s found %d results' % (self.query,self.hit_count))
        start = time.time()
        articles = list()
        for page_start in range(0,self.hit_count,self.page_size):
            more_articles, discard = self._get_articles(page_start=page_start)
            articles += more_articles
            download_count =len(articles)
            print("Downloaded %d out of %d hits in %s" % (download_count,self.hit_count, execution_time(start)))
        return articles


    def __dumps(self, articles:list, into_file:str):
        try:
            dump = open(into_file, "w", encoding='utf-8')
            dump.write(json.dumps(articles,indent=1))
        except FileNotFoundError:
            # path was specified incorrectly, dumps into the current dir
            dump = open(self.search_name+'.json', "w", encoding='utf-8')
            dump.write(json.dumps(articles,indent=1))
 
    def __init__(self,query, search_name, APIconfig:dict, etm_dump_dir:str, etm_stat_dir, add_params={}):
        super().__init__(APIconfig,add_param=add_params)
        self.params.pop('limit')
        #self.url = APIconfig['ETMURL']  #ETMurl = 'https://discover.elseviertextmining.com/api/'
        #self.api_key = APIconfig['ETMapikey']
        #self.query = query
        #self.__set_params()
        #if add_params:
         #   self.params.update(add_params)
        #result = self.__get_result()
        articles, self.hit_count = self._get_articles()
        #result['total-hits='] if result else 0
        #print('\nETM query %s returns %d documents' % (self.query, self.hit_count))
        self.search_name = search_name
        self.etm_results_dir = etm_dump_dir
        self.etm_stat_dir = etm_stat_dir
        #self.references = set()


    def load_from_json(self,use_cache=True):
        if use_cache:
            dumpfile_name = self.etm_results_dir+self.search_name+'.json'
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
                
                if (self.hit_count > len(articles)): # hit_count is 1-based index!!!
                    print('Today ETM search finds %d results. %d more than in cache' %(self.hit_count, self.hit_count-len(articles)))
                    file_date_modified = os.path.getctime(dumpfile_name)
                    self.params.update({'so_d':str(file_date_modified)}) #2003-03-19
                    #result = self.__get_result()
                    #hit_count = result['total-hits='] if result else 0
                    update_articles = self.__download()
                    articles = articles + update_articles
                    self.__dumps(articles,dumpfile_name)
            except FileNotFoundError:
                #if json dump file is not found new ETM search is initiated
                articles = self.__download()
                self.__dumps(articles,dumpfile_name)            
        else:
            articles = self.__download()

        for article in articles:
            etm_ref = ETMjson(article)
            relevance_score = float(article['score'])
            etm_ref['Relevance rank'] = [relevance_score]
            if hasattr(etm_ref,"Identifiers"):
                #self.references.add(etm_ref)
                self._add2counter(etm_ref)

        
    def get_statistics(self, stat_prop_list):
       # stat_prop_list = [PUBYEAR,AUTHORS,INSTITUTIONS,JOURNAL]
        for p in stat_prop_list:
            self.statistics[p] = dict()

        self.term2refs = dict() # {term:{ref}}
        self.keywords2ref = dict() #{keyword:{ref}}
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


    @staticmethod
    def _sort_dict(dic:dict, sort_by_value = True, reverse = True):
        item_idx = 1 if sort_by_value else 0
        return dict(sorted(dic.items(), key=lambda item: item[item_idx],reverse=reverse))


    def to_excel(self, stat_props:list):
        workbook = xlsxwriter.Workbook(self.etm_stat_dir+self.search_name+'.xlsx')
        for p in stat_props:
            try:
                dic = dict(self.statistics[p])
                sort_by_value = True if p != PUBYEAR else False
                prop_stat = self._org2address_stats() if p == INSTITUTIONS else self._sort_dict(dic,sort_by_value)
                DocMine.dict2worksheet(workbook,p,prop_stat, [p, '#References'],100)
            except KeyError: continue

        worksheet_ref = workbook.add_worksheet('References')
        worksheet_ref.write_string(0,0,'Total number of articles: '+str(len(self.references)))
        row_counter = 1
        for ref in self.references:
            row_str = ref.to_str(['PMID','DOI','PUI'],col_sep='\n')+'\n\t\t\t\n'
            worksheet_ref.write_string(row_counter,0,row_str)
            row_counter += 1

        workbook.close()


   
    @staticmethod
    def count_refs(ref_counter:set, references:list):
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
    def print_citation_index(ref_counter:set,fname:str):
        to_sort = list(ref_counter)
        to_sort.sort(key=lambda x: x['Citation index'][0], reverse=True)
        with open(fname, 'w', encoding='utf-8') as f:
            for ref in to_sort:
                ref_str = ref.get_biblio_str()
                count = ref['Citation index'][0]
                f.write(str(count)+'\t'+ref_str+'\n')


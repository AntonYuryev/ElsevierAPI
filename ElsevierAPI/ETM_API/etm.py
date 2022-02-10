from ElsevierAPI.ETM_API.references import DocMine
from ElsevierAPI.ETM_API.references import AUTHORS,INSTITUTIONS,JOURNAL,PUBYEAR,SENTENCE,EMAIL,PS_ID_TYPES
from ElsevierAPI import execution_time
import urllib.request
import urllib.parse
import json
import time
import xlsxwriter

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
                        term_id = markup['term']['id']
                        term_name = markup['term']['value']
                        term_ids.add(term_id+'\t'+term_name)
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



class ETM:
    url = 'https://demo.elseviertextmining.com/api'
    api_key = str()
    request_type = '/search/advanced'
    page_size = 100
    def __base_url(self): return self.url+self.request_type+'?'
    searchTarget = 'full_index'
    number_snippets = 1
    hit_count = 0
    search_name = str()
    query = str()
    url_request = str()
    params = dict()
    etm_results_dir = ''
    etm_stat_dir = ''
    statistics = dict() # {prop_name:{prop_value:count}}
    references = set() #set of DocMine objects

    def __set_params(self):
        self.params = {
        'query':self.query,
        'searchTarget': self.searchTarget,
        'snip': str(self.number_snippets)+'.desc',
        'apikey':self.api_key,
        'limit':self.page_size
        }
            

    def __get_param_str(self):
        return urllib.parse.urlencode(self.params)

    def __get_result(self):
        param_str = self.__get_param_str()
        full_url = self.__base_url()+param_str
        the_page = urllib.request.urlopen(full_url).read()
        if the_page:
            return json.loads(the_page.decode('utf-8'))
        else:
            return []
 
    def __init__(self,query, search_name, APIconfig:dict, etm_dump_dir:str, etm_stat_dir, add_params={}):
        self.url = APIconfig['ETMURL']  #ETMurl = 'https://discover.elseviertextmining.com/api/'
        self.api_key = APIconfig['ETMapikey']
        self.query = query
        self.__set_params()
        if add_params:
            self.params.update(add_params)
        result = self.__get_result()
        self.hit_count = result['total-hits='] if result else 0
        #print('\nETM query %s returns %d documents' % (self.query, self.hit_count))
        self.search_name = search_name
        self.etm_results_dir = etm_dump_dir
        self.etm_stat_dir = etm_stat_dir
        self.references = set()


    def load_from_json(self, use_cache = True):
        dumpfile_name = self.etm_results_dir+self.search_name+'.json'
        search_again = False
        if use_cache:
            try:
                # attepts to find json file with ETM results saved after ETM API call below
                f = open(dumpfile_name,'r')
                articles = json.load(f)
                path_end = dumpfile_name.rfind('/')
                if path_end < 0:
                    path_end = dumpfile_name.rfind('\\')
                fname = dumpfile_name[path_end+1:]
                cache_path = dumpfile_name[:path_end]
                print('Found "%s" file in "%s" cache with %d articles. Will use it to load results' % (fname,cache_path,len(articles)))
                if (self.hit_count > len(articles)): # hit_count is 1-based index!!!
                    print('Current ETM search finds %d results. %d more than in cache' %(self.hit_count, self.hit_count-len(articles)))
                f.close()
            except FileNotFoundError:
                search_again = True
                pass

        if search_again or not use_cache:
            #if json dump file is not found new ETM search is initiated
            print('\nPerforming ETM search "%s"' % self.search_name)
            print ('Query: %s found %d results' % (self.query,self.hit_count))
            start = time.time()
            articles = list()
            for page in range(0,self.hit_count,self.page_size):
                self.params['start'] = page
                result = self.__get_result()
                articles += result['article-data']
                download_count =len(articles)
                print("Downloaded %d out of %d hits in %s" % (download_count,self.hit_count, execution_time(start)))
            
            if use_cache:
                with open(dumpfile_name, "w", encoding='utf-8') as dump:
                    dump.write(json.dumps(articles,indent=1))

        relevance_rank = 1
        for article in articles:
            etm_ref = ETMjson(article)
            etm_ref['Relevance rank'] = [relevance_rank]
            relevance_rank += 1
            if hasattr(etm_ref,"Identifiers"):
                self.references.add(etm_ref)

        
    def get_statistics(self, stat_prop_list):
       # stat_prop_list = [PUBYEAR,AUTHORS,INSTITUTIONS,JOURNAL]
        for p in stat_prop_list:
            self.statistics[p] = dict()

        self.term2refs = dict() # {term:{ref}}
        self.keywords2ref = dict() #{keyword:{ref}}
        for ref in self.references:
       #     if isinstance(ref, ETMjson):
       #         self.references.add(ref)
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
        for ref in self.references:
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
    def etm_relevance(terms:list, APIconfig, add_param:dict={}):
        # add_param controls number of best references to return. Defaults to 5
        try: 
            limit = add_param['limit']
        except KeyError:
            limit = 5

        query = '{'+'};{'.join(terms)+'}'
        params = {
        'query':query,
        'searchTarget': 'full_index',
        'snip': '1.desc',
        'apikey':APIconfig['ETMapikey'],
        'limit': limit,
        }

        params.update(add_param)
        baseurl = APIconfig['ETMURL']+'/search/basic?'
  
        articles = list()
        param_str = urllib.parse.urlencode(params)
        url_request = baseurl+param_str
        the_page = urllib.request.urlopen(url_request).read()
        if the_page:
            result = json.loads(the_page.decode('utf-8'))
            hit_count = result['total-hits=']
            articles += result['article-data']
        else:
            hit_count = 0

        if limit > 100:
            for page in range(len(articles), limit, 100):
                params['start'] = page
                param_str = urllib.parse.urlencode(params)
                page = urllib.request.urlopen(url_request).read()
                result = json.loads(the_page.decode('utf-8'))
                articles += result['article-data']

        references = list()
        relevance_rank = 1
        for article in articles:
            etm_ref = ETMjson(article)
            etm_ref['Relevance rank'] = [relevance_rank]
            relevance_rank += 1
            if hasattr(etm_ref,"Identifiers"):
                references.append(etm_ref)

        ref_ids = list()
        for ref in references:
            for id_type in PS_ID_TYPES:
                try:
                    identifier = ref.Identifiers[id_type]
                    ref_ids.append(identifier)
                    break
                except KeyError:
                    continue
        
        best_refs_str = ';'.join(ref_ids)
        return hit_count, best_refs_str, references
        

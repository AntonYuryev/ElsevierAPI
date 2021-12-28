from ElsevierAPI.ETM_API.references import DocMine
from ElsevierAPI.ETM_API.references import REF_ID_TYPES,AUTHORS,INSTITUTIONS,JOURNAL,PUBYEAR,SENTENCE, EMAIL
from ElsevierAPI import execution_time
import urllib.request
import urllib.parse
import json
import time
import xlsxwriter
import re


class ETMdoc:
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
    params = dict()
    articles = list() # list of DocMine objects
    etm_results_dir = 'D:/Python/ENTELLECT_API/ETM outputs/Aerosolization taxonomy/'
    statistics = dict() # {prop_name:{prop_value:count}}
    references = set()

    @staticmethod
    def dump_fname(search_name): return search_name + '.json'

    def __set_params(self):
        if not self.params:
            self.params = {
            'query':self.query,
            'searchTarget': self.searchTarget,
            'snip': str(self.number_snippets)+'.desc',
            'apikey':self.api_key,
            'limit':self.page_size
            # 'so_d' = '2021-06-19 -' # from date inclusive
            }

    def __get_param_str(self):
        return urllib.parse.urlencode(self.params)

    def __get_result(self):
        param_str = self.__get_param_str()
        req = urllib.request.Request(url=self.__base_url()+param_str)
        the_page = urllib.request.urlopen(req).read()
        return json.loads(the_page.decode('utf-8'))
 
    def __init__(self,query, search_name, APIconfig:dict):
        self.url = APIconfig['ETMURL']  #ETMurl = 'https://discover.elseviertextmining.com/api/'
        self.api_key = APIconfig['ETMapikey']
        self.query = query
        self.__set_params()
        result = self.__get_result()
        self.hit_count = result['total-hits=']
        print('ETM query %s returns %d documents' % (self.query, self.hit_count))
        self.search_name = search_name
        

    def load_from_json(self, use_cache = True):
        dumpfile_name = self.etm_results_dir+self.search_name+'.json'
        if use_cache:
            try:
                # attepts to find json file with ETM results saved after ETM API call below
                articles = json.load(open(dumpfile_name))
                print('Found cache %s file. Will use it to load results' % dumpfile_name)
                return articles
            except FileNotFoundError:
                pass

            print('Performing ETM search "%s"' % self.search_name)
            print ('Query: %s' % self.query)
            start = time.time()

            articles = list()
            for page in range(1,self.hit_count,self.page_size):
                self.params['start'] = page
                result = self.__get_result()
                articles = articles + result['article-data']
                download_count = min(page+self.page_size,self.hit_count)
                print("Downloaded %d out of %d hits in %s" % (download_count,self.hit_count, execution_time(start)))

            with open(dumpfile_name, "w", encoding='utf-8') as dump:
                dump.write(json.dumps(articles,indent=1))
        return articles


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
                return [str(abstract_secs['addressOrAlternativesOrArray']['p'])]
        
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

        return names[:last_name_index+1], address, email


    @staticmethod
    def __parse_contributors(article:dict):
        try:
            author_info = article['article']['front']['article-meta']['contribGroupOrAffOrAffAlternatives']
            if not author_info:
                print('Article has no author info') 
                return [],[]
        except KeyError:
            print('Article has no author info')
            return [],[]
            
        if isinstance(author_info, dict): author_info = [author_info]

        authors = list()
        institutions = list()
        for item in author_info:
            try:
                instituts = item['aff']['content']
                if not isinstance(instituts, list): instituts = [instituts]
                for inst in instituts:
                    org_names, address, email = ETMdoc.__parse_institution(inst['institution'])
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


    def make_ref(self, article:dict):
        is_article, article_ids = self.__parse_ids(article)
        if not is_article: return dict(), set(), set()
    
        for id_type in REF_ID_TYPES:
            try:
                ref = DocMine(id_type,article_ids[id_type])
                ref.Identifiers.update(article_ids)
                break
            except KeyError: continue

        try:
            title = article['article']['front']['article-meta']['title-group']['article-title']
            ref.append_property('Title', title)
        except KeyError: pass

        try:
            year = article['article']['front']['article-meta']['pub-date']['year']
            ref.set_date(year)
        except KeyError: pass

        ref.add2section('Abstract',self.__parse_abstact(article))
        authors,institutions = self.__parse_contributors(article)
        ref[AUTHORS] = authors
        #ref[INSTITUTIONS] = [i[0] for i in institutions]
        ref.addresses = {i[0]:i[1] for i in institutions}
        
        try:
            journal = article['article']['front']['journal-meta']['journal-title-group']['journal-title']
            journal = ','.join(journal.values())
        except KeyError:
            journal = 'Grant application'
        ref.append_property(JOURNAL, journal)

        text_ref = ref._make_textref()
        snippet, term_ids, keywords = self.__parse_snippet(article)
        ref.add_sentence_prop(text_ref, SENTENCE, snippet)

        return ref, term_ids, keywords #term_ids = {term_id+'\t'+term_name}


    def get_statistics(self, articles:list, stat_prop_list):
       # stat_prop_list = [PUBYEAR,AUTHORS,INSTITUTIONS,JOURNAL]
        for p in stat_prop_list:
            self.statistics[p] = dict()

        self.term2refs = dict() # {term:{ref}}
        self.keywords2ref = dict() #{keyword:{ref}}
        for article in articles:
            ref, term_ids, keywords = self.make_ref(article)
            if isinstance(ref, DocMine):
                self.references.add(ref)
                for p in stat_prop_list:
                    ref.count_str(self.statistics[p], p)

                for term in term_ids: #term_ids = {term_id+'\t'+term_name}
                    try:
                        self.term2refs[term].add(ref)
                    except KeyError:
                        self.term2refs[term] = {ref}
                
                for k in keywords:
                    try:
                        self.term2refs[k].add(ref)
                    except KeyError:
                        self.term2refs[k] = {ref}


    def __org2adress_stats(self):
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
    def __sort_dict(dic:dict(), sort_by_value = True, reverse = True):
        item_idx = 1 if sort_by_value else 0
        return dict(sorted(dic.items(), key=lambda item: item[item_idx],reverse=reverse))


    def to_excel(self, stat_props:list):
        workbook = xlsxwriter.Workbook(self.etm_results_dir+self.search_name+'.xlsx')
        for p in stat_props:
            try:
                dic = dict(self.statistics[p])
                sort_by_value = True if p != PUBYEAR else False
                prop_stat = self.__org2adress_stats() if p == INSTITUTIONS else self.__sort_dict(dic,sort_by_value)
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
        
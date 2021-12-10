import urllib.request
import urllib.parse
import json
from ElsevierAPI import load_api_config
from ElsevierAPI.ResnetAPI.NetworkxObjects import Reference, REF_ID_TYPES
from ElsevierAPI.ResnetAPI.ResnetAPISession import APISession
import time
import xlsxwriter

REFERENCE_TYPES = REF_ID_TYPES | {'REPORTER','GRANTNUMREPORTER'}
INSTITUTION_KEYWORDS = {'institute', 'institut', 'clinic', 'hospital', 'university', 'universitat', 'universiti', 'centre', 'center', 'inc', 'colleges',
                        'ltd','gmbh', 'school','politecnic','politecnico', 'college', 'department', 'division', 'council', 'academy'}

def OpenFile(fname):
    open(fname, "w", encoding='utf-8').close()
    return open(fname, "a", encoding='utf-8')

def make_json_dump_name(search_name):
    return search_name + '.json'

def count_str(counter:dict, name:str):
    try:
        current_count = counter[name]
        counter[name] = current_count +1
    except KeyError:
        counter[name] = 1

def has_institution_keyword(name:str):
    name_words = name.split(' ')
    for w in name_words:
        if w.lower() in INSTITUTION_KEYWORDS: return True
    return False


def download_json_etm(search_name, query, APIconfig:dict):
    ETMurl = APIconfig['ETMURL']  #ETMurl = 'https://discover.elseviertextmining.com/api/'
    request_type = '/search/advanced'
    baseURL = ETMurl+request_type+'?'
    PageSize = 100 # cannot be more than 100

    params = {'query':query,'searchTarget':'full_index','snip':'1.desc','apikey':APIconfig['ETMapikey'],'limit':PageSize}
    #params['so_d'] = '2021-06-19 -' #from date inclusive
    paramtURL = urllib.parse.urlencode(params)
    req = urllib.request.Request(url=baseURL+paramtURL)
    
    dump_fname = make_json_dump_name(search_name)
    dump_file = OpenFile(dump_fname)
    start = time.time()

    with urllib.request.urlopen(req) as response:
        the_page = response.read()
        result = json.loads(the_page.decode('utf-8'))
        HitCount = result['total-hits='] 
        print('Downloading %d hits from ETM' % HitCount)

    articles = list()
    for page in range(1,HitCount,PageSize):
        params['start'] = page
        paramtURL = urllib.parse.urlencode(params)
        req = urllib.request.Request(url=baseURL+paramtURL)
        response = urllib.request.urlopen(req)
        the_page = response.read()
        result = json.loads(the_page.decode('utf-8'))
        articles = articles + result['article-data']

        download_count = min(page+PageSize,HitCount)
        print("Downloaded %d out of %d hits in %s" % (download_count,HitCount, APISession.execution_time(start)))

    dump_file.write(json.dumps(articles,indent=1))
    dump_file.close()
    return HitCount, dump_fname

def dict2worksheet(workbook:xlsxwriter.Workbook, worksheet_name, dict2print:dict, header:list, key_col_width=50):
    worksheet = workbook.add_worksheet(worksheet_name)
    col = 0
    for h in header:
        worksheet.write(0, col, h)
        col +=1
        
    row = 1
    for k,v in dict2print.items():
        worksheet.write(row, 0, k)
        worksheet.write(row, 1, v)
        row += 1

    worksheet.set_column('A:A', key_col_width)

def parse_ids(article:dict, article_ids:dict):
    ids = article['article']['front']['article-meta']['article-id']
    if not isinstance(ids,list): ids = [ids]
    is_conference_proceeeding = False
    for article_id in ids:
        id_type = article_id['pub-id-type'].upper()
        Id = str(article_id['value'])
        if id_type not in article_ids.keys():
            article_ids[id_type] = Id
        else: 
            return  False #article has more than one identifier of the same type - it is conference proceedings
    
    return True

def parse_institution(institution:str):
    email_pos = institution.lower().find('; e')
    if email_pos < 0:  email_pos = institution.find('. e')

    address = institution[0:email_pos] if email_pos > 0  else institution
    if (not address[0].isalpha()) or address[0].islower(): 
        address = address[1:].strip()
    names = address.split(', ')
    #to_return = [names[0]]

    last_name_index = 0
    for n in range(1,len(names)):
        if has_institution_keyword(names[n]):
            last_name_index = n
        else: break

    return names[:last_name_index+1], address, institution[email_pos+2:]


def parse_abstacts(article:dict):
    try:
        abstract_secs = article['article']['front']['article-meta']['abstract']['sec']
    except KeyError: 
        abstract_secs = article['article']['body']

    if isinstance(abstract_secs, dict):
        try:
            abstract_secs = abstract_secs['sec']
            if isinstance(abstract_secs, dict): 
                abstract_secs = [abstract_secs]
        except KeyError:
            return [abstract_secs['addressOrAlternativesOrArray']['p']]
    
    #at this point abstract_secs can be only list of dict
    abstract_text = str() 
    for paragraph in abstract_secs:
        try: 
            abs_title = paragraph['title'] +'\n' #some abstracts are secionalized
        except KeyError: 
            abs_title = ''

        try:
            abs_text = paragraph['addressOrAlternativesOrArray']['p']
        except KeyError:
            abs_text = paragraph['sec']['addressOrAlternativesOrArray']['p']
        abstract_text = abstract_text + abs_title + abs_text+'\n'
    return [abstract_text[:-1]]


def parse_contributors(article:dict):
    try:
        author_info = article['article']['front']['article-meta']['contribGroupOrAffOrAffAlternatives']
    except KeyError:
        print('Article has no author info')
        return []
        
    if isinstance(author_info, dict): author_info = [author_info]

    if len(author_info) == 0:
        print('Article has no author info') 
        return []
 
    contributors = set()
    authors = list()
    institutions = list()
    for item in author_info:
        try:
            instituts = item['aff']['content']
            if not isinstance(instituts, list): instituts = [instituts]
            for inst in instituts:
                org_names, address, email = parse_institution(inst['institution'])
                institutions.append((org_names[-1],address)) #will take only last name of institution defining global organization
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
                    #count_str(auth_counter,au_surname)
                    authors.append(au_name)
                except KeyError: 
                    continue

    return authors, institutions

def parse_snippet(article:dict):
    term_ids = set()
    try:
        snippet = str(article['data']['snippet']['text'])
        pos_shift = 0
        highlights = article['data']['snippet']['highlight']
        if not isinstance(highlights, list): highlights = [highlights]
        for markup in highlights:
            try:
                query_term = markup['queryTerm']
                query_start = markup['start']
                snippet = snippet[:pos_shift+query_start] + '{'+query_term+'}='+snippet[pos_shift+query_start:]
                pos_shift += len(query_term)+3
            
                try:
                    term_id = markup['term']['id']
                    term_name = markup['term']['value']
                    term_ids.update([term_id+'\t'+term_name])
                except KeyError:
                    continue
            except KeyError: continue
    except KeyError:
        try:
            snippet = article['article']['body']['sec']['addressOrAlternativesOrArray']['p']
        except KeyError:
            snippet = ''

    return snippet, term_ids

if __name__ == "__main__":
    APIconfig = load_api_config()
    #query = '((inhaled OR inhalable OR inhalant OR volatile OR aerosol OR mist) NEXT/2 ({drugs}/exp OR {natural products and their synthetic derivatives}/exp)) OR (({drugs}/exp OR {natural products and their synthetic derivatives}/exp) NEXT/2 (inhaled OR inhalable OR inhalant OR aerosol OR mist))'
    search_name = 'PTH inhalation'
    query = 'parag({PTH} AND {inhalation}) OR title({PTH} AND {inhalation}) OR abs({PTH} AND {inhalation})'

    file_result_name = make_json_dump_name(search_name)

    #hit_count,file_result_name = download_json_etm(search_name, query, APIconfig)

    term2refs = dict()
    auth_counter = dict()
    institution_counter = dict()
    year_counter = dict()
    journal_counter = dict()
    orgname2address = dict()
    articles = json.load(open(file_result_name))
    dump_file_name = search_name+'.tsv'
    dump_file = OpenFile(dump_file_name)

    article_counter = 0
    for article in articles:
        try:
            year = article['article']['front']['article-meta']['pub-date']['year']
            count_str(year_counter,year)
        except KeyError:
            year = ''

        article_ids = dict()
        if not parse_ids(article,article_ids): continue
        article_counter += 1

        abstracts = parse_abstacts(article)
        authors,institutions = parse_contributors(article)

        #creating reference object using one of article identifiers
        for id_type in REFERENCE_TYPES:
            try:
                ref = Reference(id_type,article_ids[id_type])
                ref.Identifiers.update(article_ids)
                break
            except KeyError: continue

        try:
            article_title = article['article']['front']['article-meta']['title-group']['article-title']
        except KeyError: 
            article_title = ''

        try:
            journal = article['article']['front']['journal-meta']['journal-title-group']['journal-title']
            journal = ','.join(journal.values())
        except KeyError:
            journal = 'Grant application'               
        
        count_str(journal_counter,journal)
 

        uniq_instituts = set()
        for i in institutions:
            if i[0] not in uniq_instituts:
                count_str(institution_counter,i[0])
                uniq_instituts.add(i[0])
            orgname2address[i[0]] = i[1]

        for a in authors:
            author_last_name = a[:str(a).find(' ')]
            count_str(auth_counter,author_last_name)
        
        if len(authors) == 0 and len(institutions) == 0:
            print('Article \"%s\" has alternative author info' % article_title)

        ref.append_property('Title',article_title)
        ref.append_property('PubYear', year)
        ref.append_property('Journal', journal)
        ref.append_property('Abstract', abstracts[0])
        ref.update_with_list('Authors',authors)
        ref.update_with_list('Affiliations',[i[0] for i in institutions])

        snippet, term_ids = parse_snippet(article)

        text_ref = ref._make_textref()
        ref.add_sentence_prop(text_ref, 'Sentence', snippet)
        dump_file.write(ref.to_str(['PMID','DOI','PUI'],col_sep='\n')+'\n\t\t\t\n')

        for term in term_ids:
            try:
                term2refs[term].update([ref])
            except KeyError:
                term2refs[term] = {ref}

    dump_file.close()
    workbook = xlsxwriter.Workbook(search_name+'.xlsx')
    
    auth_counter = dict(sorted(auth_counter.items(),key=lambda item: item[1],reverse=True))
    dict2worksheet(workbook,'Authors',auth_counter, ['Author', '#References'], 20)
    
    org_address_couner = dict()
    for k,v in institution_counter.items():
        address = orgname2address[k]
        org_address_couner[address] = v
    institution_counter = dict(sorted(org_address_couner.items(), key=lambda item: item[1],reverse=True))
    dict2worksheet(workbook,'Institutions',institution_counter, ['Institution', '#References'],100)

    year_counter = dict(sorted(year_counter.items(), key=lambda item: item[0],reverse=True))
    dict2worksheet(workbook,'Timeline',year_counter, ['Year', '#References'],15)
    journal_counter = dict(sorted(journal_counter.items(), key=lambda item: item[1],reverse=True))
    dict2worksheet(workbook,'Journals',journal_counter, ['Journal', '#References'])


    worksheet = workbook.add_worksheet('References')
    worksheet.write_string(0,0,'Total number of articles: '+str(article_counter))
    #references_str = str()
    with open(dump_file_name, "r", encoding='utf-8') as f:
        line = '\n'
        excel_insert_row = 1
        record = str()
        while line:
            line = f.readline()
            if line != '\n':
                record = record + line
            if line == '\t\t\t\n':
                worksheet.write_string(excel_insert_row,0,record)
                record = ''
                excel_insert_row += 1
                

    workbook.close()

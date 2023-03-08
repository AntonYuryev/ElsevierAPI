import urllib.request
import urllib.parse
import json
from ElsevierAPI import load_api_config
from ElsevierAPI.ETM_API.references import Reference,REF_ID_TYPES
from ElsevierAPI.ResnetAPI.ResnetAPISession import APISession
from threading import Thread
import time


def get_children(terms:list = None):
    if isinstance(terms,list):
        for term in terms:
            term_id, term_name = term.split('\t')
            term_msid= str(term_id[1:])
            if term_msid.isdigit():
                return ps_api._props2psobj([term_msid],['MedScan ID'])
            else:
                return ps_api._props2psobj([term_name],['Name','Alias'])

if __name__ == "__main__":
    APIconfig = load_api_config()
    #ETMurl = 'https://demo.elseviertextmining.com/api/'
    ETMurl = 'https://discover.elseviertextmining.com/api/'
    request_type = 'search/advanced'
    baseURL = ETMurl+request_type+'?'
    PageSize = 100 # cannot be more than 100

    query = '((inhaled OR inhalable OR inhalant OR volatile OR aerosol OR mist) NEXT/2 ({drugs}/exp OR {natural products and their synthetic derivatives}/exp)) OR (({drugs}/exp OR {natural products and their synthetic derivatives}/exp) NEXT/2 (inhaled OR inhalable OR inhalant OR aerosol OR mist))'
    #query = '(({drugs}/exp OR {natural products and their synthetic derivatives}/exp) NEXT/2 (aerosol OR mist)) OR ((aerosol OR mist) NEXT/2 ({drugs}/exp OR {natural products and their synthetic derivatives}/exp))'
    #query = '{drugs}/exp NEXT inhalants' # returns 83 hits, used for testing
    #query = 'volatile NEXT/2 ({drugs}/exp OR {natural products and their synthetic derivatives}/exp)'

    params = {'query':query,'searchTarget':'full_index','snip':'1.desc','apikey':APIconfig['ETMapikey'],'limit':PageSize}
    #params['so_d'] = '2021-06-19 -' #from date inclusive

    paramtURL = urllib.parse.urlencode(params)
    req = urllib.request.Request(url=baseURL+paramtURL)
    
    ps_api = APISession(APIconfig['ResnetURL'], APIconfig['PSuserName'], APIconfig['PSpassword'])
    dump_file = ps_api.reopen("ETMdump.tsv")
    term2refs = dict()
    start = time.time()

    t = Thread(target=get_children, args=([None]))
    t.start() #initializing thread

    with urllib.request.urlopen(req) as response:
        the_page = response.read()
        result = json.loads(the_page.decode('utf-8'))
        HitCount = result['total-hits='] 
        print('Downloading %d hits from ETM' % HitCount)

        for page in range(1,HitCount,PageSize):
            articles = result['article-data']
            for article in result['article-data']:
                try:
                    journal = article['article']['front']['journal-meta']['journal-title-group']['journal-title']
                except KeyError:
                    journal = 'Grant application'

                try:
                    author_info = article['article']['front']['article-meta']['contribGroupOrAffOrAffAlternatives']
                except KeyError: author_info = []
                if isinstance(author_info, dict): author_info = [author_info]

                if len(author_info) == 0:
                    print('Article has no author info')
                else:
                    institutions = list()
                    authors = list()            
                    for item in author_info:
                        try:
                            contribOrAddressOrAff_list = item['contrib-group']['contribOrAddressOrAff']
                            if not isinstance(contribOrAddressOrAff_list, list): contribOrAddressOrAff_list = [contribOrAddressOrAff_list]

                            au_name = ''
                            for author in contribOrAddressOrAff_list:
                                try:
                                    au_name = author['contrib']['contribIdOrAnonymousOrCollab']['name']['given-names']['content']
                                    #if au_name == 'Tao':
                                    #    print('')

                                except KeyError:
                                    try:
                                        au_name = author['contrib']['contribIdOrAnonymousOrCollab']['name']['given-names']['initials']+'.'
                                    except KeyError: continue

                                try:
                                    au_name = au_name + ' ' + author['contrib']['contribIdOrAnonymousOrCollab']['name']['surname']
                                except KeyError: continue

                                authors.append(au_name)

                        except KeyError: 
                            try:
                                instituts = item['aff']['content']
                                if not isinstance(instituts, list): instituts = [instituts]
                                for inst in instituts:
                                    institutions.append(inst['institution'])
                            except KeyError: continue
                    
                    if len(authors) == 0 and len(institutions) == 0:
                        print('Article has alternative author info')



                if isinstance(journal,dict):
                    journal = ','.join(journal.values())

                try:
                    year = article['article']['front']['article-meta']['pub-date']['year']
                except KeyError:
                    year = ''

                article_ids = dict()
                ids = article['article']['front']['article-meta']['article-id']
                if not isinstance(ids,list): ids = [ids]
                for article_id in ids:
                    id_type = article_id['pub-id-type'].upper()
                    Id = str(article_id['value'])
                    article_ids[id_type] = Id

                for id_type in REF_ID_TYPES:
                    try:
                        ref = Reference(id_type,article_ids[id_type])
                        ref.Identifiers.update(article_ids)
                        break
                    except KeyError: continue

                try:
                    article_title = article['article']['front']['article-meta']['title-group']['article-title']
                except KeyError: 
                    article_title = ''

                ref.append_property('Title',article_title)
                ref.append_property('PubYear', year)
                ref.append_property('Journal', journal)
                ref.update_with_list('Authors',authors)
                ref.update_with_list('Affiliations',institutions)
                #dump_file.write(ref.to_string(['PMID','DOI'])+'\n')

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


                text_ref = ref._make_textref()
                ref.add_sentence_prop(text_ref, 'Sentence', snippet)
                dump_file.write(ref.to_str(['PMID','DOI'])+'\n')

                for term in term_ids:
                    try:
                        term2refs[term].update([ref])
                    except KeyError:
                        term2refs[term] = {ref}

                t.join()
                t = Thread(target=get_children, args=([term_ids]))
                t.start()

            if page+PageSize < HitCount:
                print("Downloaded %d out of %d hits" % (page+PageSize,HitCount))
                params['start'] = page+PageSize
                paramtURL = urllib.parse.urlencode(params)
                req = urllib.request.Request(url=baseURL+paramtURL)
                response = urllib.request.urlopen(req)
                the_page = response.read()
                result = json.loads(the_page.decode('utf-8'))
            else:
                t.join()
                fout_results = "ETMresults.tsv"
                file_result = ps_api.reopen(fout_results)
                file_result.write('term ID\tterm name\t#children\t#references\tPMID\tDOI\tTitle\tPubYear\tJournal\tSentence\n')
                
                fout_au = "ETMresults_authors.tsv"
                file_au = ps_api.reopen(fout_au)
                file_au.write('Authors\tterm ID\tterm name\t#children\t#references\tPMID\tDOI\tTitle\tPubYear\tJournal\tSentence\n') 

                fout_inst = "ETMresults_affiliations.tsv"
                file_inst = ps_api.reopen(fout_inst)
                file_inst.write('term ID\tterm name\t#children\t#references\tPMID\tDOI\tTitle\tPubYear\tJournal\tSentence\tInstitution\n')
                
                for term, refs in term2refs.items():
                    obj_child_ids = get_children([term])
                    number_of_children = len(obj_child_ids)
                    ref_count = str(len(refs))
                    for ref in refs:
                        row = term+'\t'+str(number_of_children)+'\t'+ref_count+'\t'+ref.to_string(['PMID','DOI'])
                        file_result.write(row+'\n')
                        for author in ref['Authors']:
                            row = author +'\t'+ term+'\t'+str(number_of_children)+'\t'+ref_count+'\t'+ref.to_string(['PMID','DOI'])
                            file_au.write(row+'\n')

                        for institution in ref['Affiliations']:
                            row = term+'\t'+str(number_of_children)+'\t'+ref_count+'\t'+ref.to_string(['PMID','DOI'])+'\t'+ institution
                            file_inst.write(row+'\n')


                print("Results are in %s file\nExecution time: %s" % (fout_results,ps_api.execution_time(start)))
                file_result.close()
                break

    dump_file.close()

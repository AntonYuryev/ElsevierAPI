import urllib.request
import urllib.parse
import json
from ElsevierAPI import APIconfig
from ElsevierAPI.ResnetAPI.NetworkxObjects import Reference, REF_ID_TYPES
from ElsevierAPI.ResnetAPI.ResnetAPISession import APISession
from threading import Thread
import time

REFERENCE_TYPES = REF_ID_TYPES + ['REPORTER','GRANTNUMREPORTER']

def get_children(terms = None):
    if terms is not None:
        for term in terms:
            term_id, term_name = term.split('\t')
            term_msid= str(term_id[1:])
            if term_msid.isdigit():
                return ps_api._get_obj_ids_by_props([term_msid],['MedScan ID'])
            else:
                return ps_api._get_obj_ids_by_props([term_name],['Name','Alias'])

if __name__ == "__main__":
    ETMurl = 'https://demo.elseviertextmining.com/api/'
    request_type = 'search/advanced'
    baseURL = ETMurl+request_type+'?'
    PageSize = 100 # cannot be more than 100

    #query = '((inhaled OR inhalable OR inhalant) NEXT/2 ({drugs}/exp OR {natural products and their synthetic derivatives}/exp)) OR (({drugs}/exp OR {natural products and their synthetic derivatives}/exp) NEXT/2 (inhaled OR inhalable OR inhalant))'
    query = '(({drugs}/exp OR {natural products and their synthetic derivatives}/exp) NEXT/2 (aerosol OR mist)) OR ((aerosol OR mist) NEXT/2 ({drugs}/exp OR {natural products and their synthetic derivatives}/exp))'
    #query = '{drugs}/exp NEXT inhalants' # returns 83 hits, used for testing

    params = {'query':query,'searchTarget':'full_index','snip':'1.desc','apikey':APIconfig['ETMapikey'],'limit':PageSize}
    paramtURL = urllib.parse.urlencode(params)
    req = urllib.request.Request(url=baseURL+paramtURL)

    ps_api = APISession(APIconfig['ResnetURL'], APIconfig['PSuserName'], APIconfig['PSpassword'])
    dump_file = ps_api.reopen("EMTdump.tsv")
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

                ref.add_property('Title',article_title)
                ref.add_property('PubYear', year)
                ref.add_property('Journal', journal) 
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


                ref.add_property('Sentence',snippet)
                dump_file.write(ref.to_string(['PMID','DOI'])+'\n')

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
                fout_name = "EMTresults.tsv"
                file_result = ps_api.reopen(fout_name)
                file_result.write('term ID\tterm name\t#children\t#references\tPMID\tDOI\tTitle\tPubYear\tJournal\tSentence')   
                for term, refs in term2refs.items():
                    obj_child_ids = get_children([term])
                    number_of_children = len(obj_child_ids)
                    ref_count = str(len(refs))
                    for ref in refs:
                        row = term+'\t'+str(number_of_children)+'\t'+ref_count+'\t'+ref.to_string(['PMID','DOI'])
                        file_result.write(row+'\n')

                print("Results are in %s file\nExecution time: %s" % (fout_name,ps_api.execution_time(start)))
                break

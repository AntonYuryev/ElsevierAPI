
import json

from datetime import timedelta
from ElsevierAPI import load_api_config
from ElsevierAPI.ETM_API.references import Reference, REF_ID_TYPES
from ElsevierAPI.ETM_API.medscan import MedScan
from textblob import TextBlob

import time
import xlsxwriter


if __name__ == "__main__":
    APIconfig = load_api_config()

    #this reads 2-column tab-delimeted file with ETM query in the 1st column, Search Name in the 2nd
    # Search Name is used as name of report files
    etm_queries = [x.strip().split('\t') for x in open('D:/Python/Patents/ETMqueries2.txt', 'r').readlines()]

    for search in etm_queries:
        search_name = search[1]
        query = search[0]

        dump_fname = make_json_dump_name(search_name)
        try:
            # attepts to find json file with ETM results saved after ETM API call below
            articles = json.load(open(dump_fname))
        except FileNotFoundError:
            # If a search is performed for the first time ETM API is called
            # search results are saved in the json file
            HitCount, dump_fname = download_json_etm(search_name, query, APIconfig)

        dump_file_name = search_name+'_cache.tsv'
        dump_file = OpenFile(dump_file_name)

        term2refs = dict()
        auth_counter = dict()
        institution_counter = dict()
        year_counter = dict()
        journal_counter = dict()
        orgname2address = dict()


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
            ref.append_property('Abstract', abstracts[0] if abstracts else '')
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
        workbook = xlsxwriter.Workbook(search_name+'_report.xlsx')
        
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

    #query = '((inhaled OR inhalable OR inhalant OR volatile OR aerosol OR mist) NEXT/2 ({drugs}/exp OR {natural products and their synthetic derivatives}/exp)) OR (({drugs}/exp OR {natural products and their synthetic derivatives}/exp) NEXT/2 (inhaled OR inhalable OR inhalant OR aerosol OR mist))'
    #search_name = 'PTH inhalation'
    #query = 'parag({PTH} AND {inhalation}) OR title({PTH} AND {inhalation}) OR abs({PTH} AND {inhalation})'
    #query = '(({drugs}/exp OR {natural products and their synthetic derivatives}/exp) NEXT/2 (aerosol OR mist)) OR ((aerosol OR mist) NEXT/2 ({drugs}/exp OR {natural products and their synthetic derivatives}/exp))'
    #query = '{drugs}/exp NEXT inhalants' # returns 83 hits, used for testing
    #query = 'volatile NEXT/2 ({drugs}/exp OR {natural products and their synthetic derivatives}/exp)'
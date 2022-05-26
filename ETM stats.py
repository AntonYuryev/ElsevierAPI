from ElsevierAPI import load_api_config
from ElsevierAPI.ETM_API.etm import ETMcache
from ElsevierAPI.ETM_API.references import AUTHORS, INSTITUTIONS, JOURNAL, PUBYEAR


if __name__ == "__main__":
    APIconfig = load_api_config()

    #reads 2-column tab-delimeted file with ETM query in the 1st column, Search Name in the 2nd
    # Search Name is used as name of report files
    file_with_searches = 'D:/Python/PMI/ETMqueries4cannabinoids.txt'
    etm_queries = [x.strip().split('\t') for x in open(file_with_searches, 'r').readlines()]

    hit_count_report_file = file_with_searches[:-4]+' request URLs.txt'
    fout = open(hit_count_report_file, 'w')
    stat_props = [AUTHORS,JOURNAL,INSTITUTIONS,PUBYEAR]
    etm_dump_dir = 'D:/Python/ENTELLECT_API/Data/ETM/Cannabinoids/'
    stats_dir = 'D:/Python/PMI/'
    for search in etm_queries:
        if len(search) < 2: continue
        search_name = search[1]
        query = search[0]
        etm_miner = ETMcache(query,search_name,APIconfig,etm_dump_dir,stats_dir)
        etm_miner.set_stat_props(stat_props)
        etm_miner.load_from_json()
        fout.write('%s\t%s\t%s\t%d\n' % (search_name,query,etm_miner._url_request(),etm_miner.hit_count))
        etm_miner.get_statistics(stat_props)
        etm_miner.to_excel(stat_props)

    fout.close()

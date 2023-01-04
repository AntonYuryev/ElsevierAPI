from ElsevierAPI import open_api_session
#from ElsevierAPI.ETM_API.references import PS_ID_TYPES,PS_BIBLIO_PROPS,SENTENCE_PROPS,CLINTRIAL_PROPS,RELATION_PROPS
from ElsevierAPI.ResnetAPI.ResnetAPISession import SNIPPET_PROPERTIES

# ps_api retreives data from the database and loads it into APISession.Graph derived from Networkx:MultiDiGraph 
ps_api = open_api_session(api_config_file='',what2retrieve=SNIPPET_PROPERTIES)#specify here path to your APIconfig file. 
#If api_config_file not specified the default APIConfig from __init__.py will be used

#all_relation_properties = list(PS_ID_TYPES)+list(PS_BIBLIO_PROPS)+list(SENTENCE_PROPS)+list(CLINTRIAL_PROPS)+list(RELATION_PROPS)

#ps_api.add_rel_props(['Name','Effect','Mechanism','ChangeType','BiomarkerType','QuantitativeType','Sentence','Title','PMID','DOI'])
#add_rel_props specifies what attributes to retreive for relations from the database. The list order defines the column order in the dump file
ps_api.add_ent_props(['Description','Alias'])
#add_ent_props specifies what attributes to retreive for nodes (entities) from the database.The list order defines the column order in the dump file
ps_api.clear_graph_cache = False #set it to True for large downloads
#by default ps_api.clear_graph_cache = False to keep all downloaded data in ps_api.Graph

pcnt = '%'
#my_goql_query = 'SELECT Relation WHERE objectType=StateChange AND CellType LIKE \'' + pcnt + 'hepatocyte' + pcnt + '\''
my_goql_query = 'SELECT Relation WHERE NeighborOf (SELECT Entity WHERE Alias = NSCLC) AND NeighborOf (SELECT Entity WHERE objectType = GeneticVariant)'
request_name = 'Find relations reported in hepatocytes'

#dafault print_rel21row = False to print 1 row per reference in every relation
ps_api.print_rel21row = False #if True ResnetAPIsessionDump.tsv will have only 1 row per each relation
# with reference properties concatenated into 1 string per property

if __name__ == "__main__":
    ps_api.start_download_from(0) #if download was interrupted change this paramater to resume download from certain position
#position must be specified as the number of relations (or entities) downloaded previously

    my_goql_query = 'SELECT Relation WHERE NeighborOf (SELECT Entity WHERE Name = PKC) AND NeighborOf (SELECT Entity WHERE Name = PDPK1)'
    my_graph = ps_api.process_oql(my_goql_query,request_name, debug=False, flush_dump=True)
    my_graph.name = ''

    ps_api.csv_delimeter = '\t'
    ps_api.to_csv('my_graph.tsv')
    ps_api.Graph.rnef2file('my_graph.rnef',ent_prop2print=ps_api.entProps, rel_prop2print=ps_api.relProps)
# process_oql retreives data by iterations. Iteration size is controled by ps_api.PageSize
# ps_api.PageSize defaults to 100 relations per iteration and cannot be bigger than 10,000
# debug=False retreives resuts only from the first iteration
# set debug=True to retreive all results from your GOQL query.
# during retreival ps_api caches data into dump file with default name 'ResnetAPIsessionDump.tsv' using APISession.to_csv()
# dump file delimeter is defined by APISession.csv_delimeter parameter
# all resultes are appended into dump file specified in ps_api.DumpFiles = ['ResnetAPIsessionDump.tsv']
# flush_dump=True deletes content of the old 'ResnetAPIsessionDump.tsv' before processing new GOQL query

# ps_api.process_oql rerutn NetworkX::my_graph object that can be writtent to file using variety of graph formats described here: 
# https://networkx.org/documentation/stable/reference/readwrite/index.html

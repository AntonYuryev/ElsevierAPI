from ElsevierAPI import open_api_session
from ElsevierAPI.ResnetAPI.NetworkxObjects import REF_PROPS,REF_ID_TYPES

ps_api = open_api_session(api_config_file=None)#specify here path to your APIconfig file. 
#If api_config_file not specified the default APIConfig from __init__.py will be used
# ps_api retreives data from the database and loads it into APISession.Graph derived from Networkx:MultiDiGraph 

ps_api.add_rel_props(REF_PROPS+REF_ID_TYPES)
#add_rel_props specifies what attributes to retreive for relations from the database
ps_api.add_ent_props(['Name','URN'])
#add_ent_props specifies what attributes to retreive for nodes (entities) from the database

pcnt = '%'
my_goql_query = 'select Relation where CellType LIKE \''+pcnt+'hepatocyte'+pcnt+'\''
my_goql_query = 'select Relation where objectType=Expression AND CellType LIKE \'' + pcnt + 'hepatocyte' + pcnt + '\''
my_goql_query = 'select Relation where objectType=StateChange AND CellType LIKE \'' + pcnt + 'hepatocyte' + pcnt + '\''
request_name = 'Find relations reported in hepatocytes'
my_graph = ps_api.process_oql(my_goql_query,request_name, debug=False, flush_dump=True)
# process_oql retreives data by iterations. Iteration size is controled by ps_api.PageSize
# ps_api.PageSize defaults to 100 relations per iteration and cannot be bigger than 10,000
# debug=False retreives resuts only from the first iteration
# set debug=True to retreive all results from your GOQL query. 
# during retreival ps_api caches data into dump file with default name 'ResnetAPIsessionDump.tsv' using APISession.to_csv()
# dump file delimeter is defined by APISession.csv_delimeter parameter
# all resultes are appended into dump file specified in ps_api.DumpFiles = ['ResnetAPIsessionDump.tsv']
# flush_dump=True deletes content of the old 'ResnetAPIsessionDump.tsv' before processing new GOQL query

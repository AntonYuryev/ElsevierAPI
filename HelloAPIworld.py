from ElsevierAPI.ResnetAPI.ResnetAPISession import APISession
from ElsevierAPI import APIconfig
from ElsevierAPI.ResnetAPI.NetworkxObjects import REF_PROPS,REF_ID_TYPES

ps_api = APISession(APIconfig['ResnetURL'],APIconfig['PSuserName'],APIconfig['PSpassword'])
# ps_api retreives data from the database and loads it into APISession.Graph derived from Networkx:MultiDiGraph 

ps_api.add_rel_props(REF_PROPS+REF_ID_TYPES)
#add_rel_props and add_end_props specify what attributes to retreive for relations and nodes (entities) from the database

pcnt = '%'
my_goql_query = 'select Relation where CellType LIKE \''+pcnt+'hepatocyte'+pcnt+'\''
my_graph = ps_api.process_oql(my_goql_query)
#process_oql retreives data by iterations. Iteration size is controled by ps_api.PageSize
# ps_api.PageSize defaults to 100 relations per iteration and cannot be bigger than 10,000
# during retreival ps_api generates cache dump file with default name 'ResnetAPIsessionDump.tsv' using to_csv()
# dump file delimeter is defined by ps_api.csv_delimeter
# output file name is in ps_api.DumpFiles = ['ResnetAPIsessionDump.tsv']






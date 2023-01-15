from ElsevierAPI import open_api_session
from ElsevierAPI.ResnetAPI.ResnetAPISession import SNIPPET_PROPERTIES,BIBLIO_PROPERTIES,REFERENCE_IDENTIFIERS,DATABASE_REFCOUNT_ONLY
from ElsevierAPI.ResnetAPI.ResnetGraph import REFCOUNT,EFFECT
from ElsevierAPI.ResnetAPI.PathwayStudioGOQL import OQL

ps_api = open_api_session(api_config_file='',what2retrieve=REFERENCE_IDENTIFIERS)
ps_api.add_rel_props([EFFECT])

request_name = 'Proteins linked to Diabetes'
select_disease = OQL.get_childs(PropertyValues=['Diabetes Mellitus'],SearchByProperties=['Name'],include_parents=True)
my_goql_query = f'SELECT Relation WHERE objectType = Regulation AND NeighborOf ({select_disease}) AND NeighborOf (SELECT Entity WHERE objectType = Protein)'

ps_api.print_rel21row = True

if __name__ == "__main__":
    my_graph = ps_api.process_oql(my_goql_query,request_name)

    fname = request_name+' relations' if ps_api.print_rel21row else request_name+' references'
    ps_api.entProps = ['ObjTypeName','Name']
    ps_api.relProps = ['ObjTypeName','Effect',REFCOUNT,'PMID','DOI']
    ps_api.to_csv(fname+'.tsv')

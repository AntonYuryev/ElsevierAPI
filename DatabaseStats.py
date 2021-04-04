import ElsevierAPI
from ElsevierAPI import networx as PSnx
import ElsevierAPI.ResnetAPI.ResnetAPISession as ssn
import time

NON_DIRECTIONAL = ['Binding','FunctionalAssociation']

def IsDirectional(relType:str):
    return relType not in NON_DIRECTIONAL


def CountRelations (EntityIN, EntityOUT, RelType, isDirectional=True):
    if isDirectional == True:
        OQLquery = 'SELECT Relation WHERE objectType = '+ RelType +' AND NeighborOf downstream (SELECT Entity WHERE objectType = '+EntityIN+') AND NeighborOf upstream (SELECT Entity WHERE objectType = '+EntityOUT+')'
    else:
        OQLquery = 'SELECT Relation WHERE objectType = '+ RelType +' AND NeighborOf downstream (SELECT Entity WHERE objectType = '+EntityIN+') AND NeighborOf downstream (SELECT Entity WHERE objectType = '+EntityOUT+')'

    sn = ssn.APISession(OQLquery, PSnx)
    relcnt = sn.GetResultSize()
    return relcnt


DBRelTypes = PSnx.GetRelationTypes()
DBEntities = PSnx.GetEntityTypes()
non_directional_counted = set()

start_time = time.time()
with open('DBstats.txt', 'w', encoding='utf-8') as f:
    f.write('Regulator\tTarget\tRelation\tRelation counts\n')
    for rel_type in DBRelTypes:
        is_directional = IsDirectional(rel_type)
        for EntityIN in DBEntities:
            for EntityOUT in DBEntities:
                if is_directional == False:
                    if EntityOUT+EntityIN in non_directional_counted: 
                        print('skipping: %s-%s-%s was already counted' % (EntityIN,rel_type,EntityOUT))
                        continue
                    else: non_directional_counted.update([EntityIN+EntityOUT])

                rel_count = CountRelations(EntityIN,EntityOUT,rel_type,is_directional)
                print('%s\t%s\t%s\t%d' % (EntityIN,EntityOUT,rel_type,rel_count))
                if rel_count > 0:
                    f.write('%s\t%s\t%s\t%d\n' % (EntityIN,EntityOUT,rel_type,rel_count))

print("Time to retreive database stats: %s" % (ElsevierAPI.ExecutionTime(start_time)))
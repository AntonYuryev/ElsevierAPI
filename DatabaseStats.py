import ElsevierAPI
from ElsevierAPI import open_api_session
import time
import argparse
import textwrap


ps_api = open_api_session()
NON_DIRECTIONAL = ['Binding','FunctionalAssociation','CellExpression','Metabolization','Paralog']

def IsDirectional(relType:str):
    return relType not in NON_DIRECTIONAL

def CountRelations (EntityIN, EntityOUT, RelType, isDirectional=True):
    if isDirectional == True:
        ps_api.GOQLquery = 'SELECT Relation WHERE objectType = '+ RelType +' AND NeighborOf downstream (SELECT Entity WHERE objectType = '+EntityIN+') AND NeighborOf upstream (SELECT Entity WHERE objectType = '+EntityOUT+')'
    else:
        ps_api.GOQLquery = 'SELECT Relation WHERE objectType = '+ RelType +' AND NeighborOf downstream (SELECT Entity WHERE objectType = '+EntityIN+') AND NeighborOf downstream (SELECT Entity WHERE objectType = '+EntityOUT+')'

    relcnt = ps_api.get_result_size()
    return relcnt


def GetStats(Triples:list):
    non_directional_counted = set()
    start_time = time.time()
    with open('DBstats.txt', 'w', encoding='utf-8') as f:
        f.write('Regulator\tTarget\tRelation\tRelation counts\n')
        for triple in Triples:
            EntityIN = triple[0]
            rel_type = triple[1]
            EntityOUT = triple[2]

            is_directional = IsDirectional(rel_type)
            if is_directional == False:
                if EntityOUT+EntityIN in non_directional_counted: 
                    print('skipping: %s-%s-%s was already counted' % (EntityIN,rel_type,EntityOUT))
                    continue
                else: non_directional_counted.update([EntityIN+EntityOUT])

            rel_count = CountRelations(EntityIN,EntityOUT,rel_type,is_directional)
            print('%s\t%s\t%s\t%d' % (EntityIN,EntityOUT,rel_type,rel_count))
            if rel_count > 0:
                f.write('%s\t%s\t%s\t%d\n' % (EntityIN,EntityOUT,rel_type,rel_count))

    print("Time to retreive database stats: %s" % (ps_api.execution_time(start_time)))


if __name__ == "__main__":
    instructions = '''
    infile - 3 column tab-delimted file with triplets description: Regulator<>Relation<tab>Target. 
    If infile is not supplied counts for all possible triplets for all entity types and all relation types will be outputed
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog=textwrap.dedent(instructions))
    parser.add_argument('-i', '--infile', default='')
    args = parser.parse_args()

    DBRelTypes = list()
    DBEntities = list()
    DBtriples_types = list()
    if len(args.infile) == 0:
        DBRelTypes = ps_api.get_relation_types()
        DBEntities = ps_api.get_entity_types()
        for rel_type in DBRelTypes:
            for EntityIN in DBEntities:
                for EntityOUT in DBEntities:
                    DBtriples_types.append([EntityIN,rel_type,EntityOUT])
    else:
        with open(args.infile, 'r', encoding='utf-8') as f:
            for line in f:
                triple = line.strip().split('\t')
                DBtriples_types.append(triple)

    GetStats(DBtriples_types)
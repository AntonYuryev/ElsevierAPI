from ElsevierAPI import open_api_session
import time
import argparse
import textwrap


ps_api = open_api_session()
NON_DIRECTIONAL = ['Binding','FunctionalAssociation','CellExpression','Metabolization','Paralog']

def IsDirectional(relType:str):
    return relType not in NON_DIRECTIONAL

def CountRelations (EntityIN, EntityOUT, RelType, isDirectional=True, TextRefs:list=None):
    if isDirectional == True:
        base_GOQLquery = 'SELECT Relation WHERE objectType = '+ RelType +' AND NeighborOf downstream (SELECT Entity WHERE objectType = '+EntityIN+') AND NeighborOf upstream (SELECT Entity WHERE objectType = '+EntityOUT+')'
    elif EntityIN != EntityOUT:
        base_GOQLquery = 'SELECT Relation WHERE objectType = '+ RelType +' AND NeighborOf (SELECT Entity WHERE objectType = '+EntityIN+') AND NeighborOf (SELECT Entity WHERE objectType = '+EntityOUT+')'
    else:
        base_GOQLquery = 'SELECT Relation WHERE objectType = '+ RelType +' AND NeighborOf downstream (SELECT Entity WHERE objectType = '+EntityIN+') AND NeighborOf downstream (SELECT Entity WHERE objectType = '+EntityOUT+')'

    textref_counts = {'Total':ps_api.get_result_size(base_GOQLquery)}
    if isinstance(TextRefs,list):
        for t in TextRefs:
            texref_clause = " AND TextRef LIKE '%{textref_type}%'".format(textref_type=t)
            textref_counts[texref_clause] = ps_api.get_result_size(base_GOQLquery + texref_clause)

    return textref_counts


def GetStats(Triples:list, TextRefs:list=None):
    non_directional_counted = set()
    start_time = time.time()
    with open('DBstats.txt', 'w', encoding='utf-8') as f:
        additional_columns = str()
        if isinstance(TextRefs,list):
            additional_columns = '\t'.join(TextRefs)
        f.write('Regulator\tTarget\tRelation\tRelation counts'+additional_columns+'\n')
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

            rel_counts_dict = CountRelations(EntityIN,EntityOUT,rel_type,is_directional,TextRefs)
            if len(rel_counts_dict) == 0:
                print('%s\t%s\t%s triple has no relations' % (EntityIN,EntityOUT,rel_type))
            else:
                for stat, rel_count in rel_counts_dict.items():
                    print('%s\t%s\t%s:' % (EntityIN,EntityOUT,rel_type))
                    print('%s - %d' % stat,rel_count)
                    f.write('%s\t%s\t%s' % (EntityIN,EntityOUT,rel_type))
                    f.write('\t%d' % rel_count)
                f.write('\n')

    print("Time to retreive database stats: %s" % (ps_api.execution_time(start_time)))


if __name__ == "__main__":
    instructions = '''
    infile - 3 column tab-delimted file with triplets description: Regulator<>Relation<tab>Target.
    TextRefs - comma-separated list of TextRef substrings to obtain document sections stats 
    If infile is not supplied counts for all possible triplets for all entity types and all relation types will be outputed

    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog=textwrap.dedent(instructions))
    parser.add_argument('-i', '--infile', default='')
    parser.add_argument('-t', '--TextRefs', default='')
    args = parser.parse_args()


    DBRelTypes = list()
    DBEntities = list()
    DBtriples_types = list()
    TextRefs  = None
    if args.TextRefs:
        TextRefs = args.TextRefs.split(',')

    if not args.infile:
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

    GetStats(DBtriples_types,TextRefs)
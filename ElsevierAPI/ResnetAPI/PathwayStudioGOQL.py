
from builtins import len
NEED_QUOTES = {' ','-','/','(',')','[',']','+','#',':','&'}
PERCNT = '%'

class OQL:
    @staticmethod
    def join_with_quotes(values:list,separator=','):
        to_return = str()
        for name in values:
            if name: #property value may not be empty
                clean_name = str(name).replace('\'',' ') # protection from properties with quote that will crash OQL request
                if 1 in [c in clean_name for c in NEED_QUOTES]:
                    to_return = to_return + '\'' + clean_name + '\'' + separator
                else:
                    to_return = to_return + clean_name + separator

        return to_return[:-1]


    @staticmethod
    def join_prop_names(prop_names:list):
        property_names = str()
        for prop_name in prop_names:
            if ' ' in prop_name:
                property_names += '\"' + prop_name + '\",'
            else:
                property_names += prop_name + ','
        return property_names[:-1]


    @staticmethod
    def get_search_strings(PropertyNameList:list, PropValuesList:list):
        some_values_have_quotes = False
        values = str()
        unique_values_list = set(PropValuesList)
        for v in unique_values_list:
            if v: #property value may not be empty
                val = str(v)
                val = val.replace('\'', '')
                if 1 in [c in val for c in NEED_QUOTES]:
                    val = '\'' + val + '\''
                    some_values_have_quotes = True
                values = values + val + ','
        values = values[:len(values) - 1]

        property_names = str()
        for n in range(0, len(PropertyNameList)):
            prop_name = str(PropertyNameList[n])
            if prop_name.find(' ') > 0:
                if not some_values_have_quotes:
                    prop_name = '\"' + prop_name + '\"'
                else:
                    message = "if you know how to concatenate string variable with double quotes and string variable with single quote in Python please let us know.\n \
                        Otherwise please remove values with white spaces from either your value list or property name list"
                    print(message)
                    raise ValueError(message)
            property_names = property_names + prop_name + ','
        property_names = property_names[:-1]

        return str(property_names), str(values)


    @staticmethod
    def get_entities_by_props(PropertyValues:list, SearchByProperties:list, only_object_types=[], MinConnectivity=0):
        if SearchByProperties[0] in ('id', 'Id', 'ID'):
            oql_query = "SELECT Entity WHERE id = " + '(' + ','.join([str(i) for i in PropertyValues]) + ')'
        else:
            prop_names, values = OQL.get_search_strings(SearchByProperties, PropertyValues)
            oql_query = "SELECT Entity WHERE (" + prop_names + ") = (" + values + ')'

        if only_object_types:
            object_types = OQL.join_with_quotes(only_object_types)
            oql_query += ' AND objectType = (' + object_types + ')'

        if MinConnectivity:
            oql_query += ' AND Connectivity >= ' + str(MinConnectivity)

        return oql_query

    @staticmethod
    def get_group_by_props(PropertyValues: list, SearchByProperties: list):
        if SearchByProperties[0] in ('id', 'Id', 'ID'):
            return "SELECT Group WHERE id = " + '(' + ','.join([str(int) for int in PropertyValues]) + ')'
        else:
            prop_names, values = OQL.get_search_strings(SearchByProperties, PropertyValues)
            return"SELECT Group WHERE (" + prop_names + ") = (" + values + ')'


    @staticmethod
    def id2str(ids:list):
        return ','.join(map(str,ids))


    @staticmethod
    def get_relations_by_props(PropertyValues:list, SearchByProperties:list, only_object_types=[], MinRef=0):
        if SearchByProperties[0] in ('id', 'Id', 'ID'):
            oql_query = "SELECT Relation WHERE id = " + '(' + OQL.id2str(PropertyValues) + ')'
        else:
            prop_names, values = OQL.get_search_strings(SearchByProperties, PropertyValues)
            oql_query = "SELECT Relation WHERE (" + prop_names + ") = (" + values + ')'

        if only_object_types:
            object_types = OQL.join_with_quotes(only_object_types)
            oql_query = oql_query + ' AND objectType = (' + object_types + ')'

        return oql_query if MinRef == 0 else oql_query + ' AND RelationNumberOfReferences >= ' + str(MinRef)


    @staticmethod
    def get_childs(PropertyValues:list,SearchByProperties:list,only_object_types=[],include_parents=False,depth=0):
        if depth:
            ontology_query = 'SELECT Entity WHERE InOntology (SELECT Annotation WHERE Ontology=\'Pathway Studio Ontology\' AND Relationship=\'is-a\') inRange '+str(depth)+' under (SELECT OntologicalNode WHERE {entities})'
        else:
            ontology_query = 'SELECT Entity WHERE InOntology (SELECT Annotation WHERE Ontology=\'Pathway Studio Ontology\' AND Relationship=\'is-a\') under (SELECT OntologicalNode WHERE {entities})'
        
        search_by_id = False
        if SearchByProperties[0] in ('id', 'Id', 'ID'):
            entity_query = "id = (" + OQL.id2str(PropertyValues) + ')'
            search_by_id = True
        else:
            prop_names, values = OQL.get_search_strings(SearchByProperties, PropertyValues)
            entity_query = '(' + prop_names + ') = (' + values + ')'

        search_query = ontology_query.format(entities=entity_query)
        if only_object_types:
            object_types = OQL.join_with_quotes(only_object_types)
            search_query = search_query + ' AND objectType = (' + object_types + ')'

        if include_parents:
          if not search_by_id:
            return str(search_query+' OR '+entity_query)
          else:
            print('Cannot include parents into ontology query which uses database Id for searching')
            return ''

        return search_query


    @staticmethod
    def expand_entity_by_id(IDlist: list, expand_by_rel_types=None, expand2neighbors=None, direction=''):
        expand2neighbors = [] if expand2neighbors is None else expand2neighbors
        expand_by_rel_types = [] if expand_by_rel_types is None else expand_by_rel_types
        expand_by_rel_types_str = OQL.join_with_quotes(expand_by_rel_types)
        expand2neighbors_str = OQL.join_with_quotes( expand2neighbors)

        values = ','.join([str(i) for i in IDlist])
        expand = 'Select Relation WHERE NeighborOf ' + direction + ' (SELECT Entity WHERE id = (' + values + '))'
        if direction == 'upstream':
            opposite_direction = 'downstream'
        elif direction == 'downstream':
            opposite_direction = 'upstream'
        else:
            opposite_direction = ''

        if expand_by_rel_types:
            if expand2neighbors:
                expand += " AND objectType = (" + expand_by_rel_types_str + ') AND NeighborOf ' + opposite_direction
                expand += ' (SELECT Entity WHERE objectType = (' + expand2neighbors_str + "))"
            else:
                expand += " AND objectType = (" + expand_by_rel_types_str + ")"
        else:
            if expand2neighbors:
                expand += ' AND NeighborOf ' + opposite_direction
                expand += ' (SELECT Entity WHERE objectType = (' + expand2neighbors_str + "))"

        return expand


    @staticmethod
    def expand_entity(PropertyValues:list, SearchByProperties:list, 
        expand_by_rel_types:list=[], expand2neighbors:list=[], direction=''):
        '''
        Input
        -----
        direction - upstream, downstream, ""
        '''
        expand2neighbors_str = OQL.join_with_quotes(expand2neighbors)
        expand_by_rel_types_str = OQL.join_with_quotes(expand_by_rel_types)

        if SearchByProperties[0] in ('id', 'Id', 'ID'):
            return OQL.expand_entity_by_id(PropertyValues, expand_by_rel_types, expand2neighbors, direction)

        property_names, values = OQL.get_search_strings(SearchByProperties, PropertyValues)
        expand = 'Select Relation WHERE NeighborOf ' + direction + ' (SELECT Entity WHERE (' + property_names + ') = (' + values + '))'

        oql_dir_hint = 'downstream' if direction == 'upstream' else 'upstream' if direction == 'downstream' else ''

        if expand_by_rel_types:
            if expand2neighbors:
                expand += " AND objectType = (" + expand_by_rel_types_str + ') AND NeighborOf ' + oql_dir_hint
                expand += ' (SELECT Entity WHERE objectType = (' + expand2neighbors_str + "))"
            else:
                expand += " AND objectType = (" + expand_by_rel_types_str + ")"
        else:
            if expand2neighbors:
                expand += ' AND NeighborOf ' + oql_dir_hint + ' (SELECT Entity WHERE objectType = ('
                expand += expand2neighbors_str + "))"

        return expand


    @staticmethod
    def get_neighbors(PropertyValues: list, SearchByProperties: list, expand_by_rel_types=None,
                    expand2neighbors=None):
        expand2neighbors = [] if expand2neighbors is None else expand2neighbors
        expand2neighbors_str = OQL.join_with_quotes( expand2neighbors)
        expand_by_rel_types = [] if expand_by_rel_types is None else expand_by_rel_types
        expand_by_rel_types_str = OQL.join_with_quotes(expand_by_rel_types)

        property_names, values = OQL.get_search_strings(SearchByProperties, PropertyValues)
        connect_to_str = " to (SELECT Entity WHERE (" + property_names + ") = (" + values + ")"
        if expand_by_rel_types:
            if expand2neighbors:
                return "SELECT Entity WHERE objectType = (" + expand2neighbors_str + ") AND Connected by (SELECT Relation WHERE objectType = " + expand_by_rel_types_str + ")" + connect_to_str + ")"
            else:
                return "SELECT Entity WHERE Connected by (SELECT Relation WHERE objectType= " + expand_by_rel_types_str + ")" + connect_to_str + ")"
        else:
            if expand2neighbors:
                return "SELECT Entity WHERE objectType = (" + expand2neighbors_str + ") AND Connected by (SELECT Relation WHERE NOT (URN = NULL))" + connect_to_str + ")"
            else:
                # no ExpandWithRelationTypes specified and no ExpandToNeighborTypes specified -> get ALL neigbors
                return "SELECT Entity WHERE Connected by (SELECT Relation WHERE NOT (URN = NULL))" + connect_to_str + ")"

    @staticmethod
    def get_objects(by_dbids:list):
        return "SELECT Entity WHERE id = (" + OQL.id2str(by_dbids) + ")"


    @staticmethod
    def select_drugs():
        return "SELECT Entity WHERE InOntology (SELECT Annotation WHERE Ontology='Pathway Studio Ontology' AND Relationship='is-a') under (SELECT OntologicalNode WHERE Name = (drugs,'plant medicinal product'))"


    @staticmethod
    def select_metabolites():
        return "SELECT Entity WHERE InOntology (SELECT Annotation WHERE Ontology='Pathway Studio Ontology' AND Relationship='is-a') under (SELECT OntologicalNode WHERE Name = ('mammal endogenous compounds and their derivatives'))"


    @staticmethod
    def drugs4(targets_with_dbids:list):
        strings = [str(integer) for integer in targets_with_dbids]
        db_ids = ",".join(strings)
        return "SELECT Relation WHERE objectType = (DirectRegulation,Binding) AND NeighborOf upstream (SELECT Entity WHERE id = (" + db_ids + ")) \
            AND NeighborOf downstream (SELECT Entity WHERE InOntology (SELECT Annotation WHERE Ontology='Pathway Studio Ontology' AND Relationship='is-a') under (SELECT OntologicalNode WHERE Name = drugs))"


    @staticmethod
    def get_reaxys_substances(ForTargetsIDlist: list):
        strings = [str(integer) for integer in ForTargetsIDlist]
        db_ids = ",".join(strings)
        return "SELECT Relation WHERE objectType = (DirectRegulation,Binding) AND NeighborOf upstream (SELECT Entity WHERE id = (" + db_ids + ")) AND Source = Reaxys"


    @staticmethod
    def connect_entities(PropertyValues1: list, SearchByProperties1: list, EntityTypes1: list, PropertyValues2: list,
                        SearchByProperties2: list, EntityTypes2: list, connect_by_rel_types:list=[]):

        prop_names1, prop_values1 = OQL.get_search_strings(PropertyNameList=SearchByProperties1, PropValuesList=PropertyValues1)
        prop_names2, prop_values2 = OQL.get_search_strings(PropertyNameList=SearchByProperties2, PropValuesList=PropertyValues2)
        object_type1 = OQL.join_with_quotes( EntityTypes1)
        object_type2 = OQL.join_with_quotes( EntityTypes2)

        entity_query = 'SELECT Entity WHERE {entity}'

        if len(SearchByProperties1) > 1:
            entity1 = '(' + prop_names1 + ') = (' + prop_values1 + ') AND objectType = (' + object_type1 + ')'
        else:
            entity1 =  prop_names1 + ' = (' + prop_values1 + ') AND objectType = (' + object_type1 + ')'
        
        if len(SearchByProperties2) > 1:
            entity2 = '(' + prop_names2 + ') = (' + prop_values2 + ') AND objectType = (' + object_type2 + ')'
        else:
            entity2 = prop_names2 + ' = (' + prop_values2 + ') AND objectType = (' + object_type2 + ')'

        entity1_query = entity_query.format(entity=entity1)
        entity2_query = entity_query.format(entity=entity2)

        oql_query = "SELECT Relation WHERE NeighborOf (" + entity1_query + ") AND NeighborOf (" + entity2_query + ")"
        if len(connect_by_rel_types) > 0:
            rel_type_list = OQL.join_with_quotes( connect_by_rel_types)
            oql_query += " AND objectType = (" + rel_type_list + ")"

        return oql_query


    @staticmethod
    def connect_ids (idlist1: list, idlist2: list, connect_by_rel_types=None, rel_effect=None, RelDirection=''):
        rel_effect = [] if rel_effect is None else rel_effect
        connect_by_rel_types = [] if connect_by_rel_types is None else connect_by_rel_types

        string_ids1 = [str(integer) for integer in idlist1]
        string_ids2 = [str(integer) for integer in idlist2]
        id1s = ",".join(string_ids1)
        id2s = ",".join(string_ids2)
        entity_query = 'SELECT Entity WHERE id = ({idlist})'
        entity1_query = entity_query.format(idlist=id1s)
        entity2_query = entity_query.format(idlist=id2s)

        oql_query = "SELECT Relation WHERE NeighborOf {dir1} (" + entity1_query + ") AND NeighborOf {dir2} (" + entity2_query + ")"
        if len(connect_by_rel_types) > 0:
            rel_type_list = OQL.join_with_quotes( connect_by_rel_types)
            oql_query = oql_query + ' AND objectType = (' + rel_type_list + ')'

        if len(rel_effect) > 0:
            effect_list = OQL.join_with_quotes( rel_effect)
            oql_query = oql_query + ' AND Effect = (' + effect_list + ')'

        if RelDirection == '<':
            oql_query = oql_query.format(dir1='upstream', dir2='downstream')
        elif RelDirection == '>':
            oql_query = oql_query.format(dir1='downstream', dir2='upstream')
        else:
            oql_query = oql_query.format(dir1='', dir2='')

        return oql_query


    @staticmethod
    def find_targets(RegulatorsIDs: list, TargetIDs:list=[], relation_types:list=[]):
        regids_str = ",".join([str(integer) for integer in RegulatorsIDs])
        oql_query = f"SELECT Relation WHERE NeighborOf downstream (SELECT Entity WHERE id = ({regids_str}))"

        if TargetIDs:
            targetids_str = ",".join([str(integer) for integer in TargetIDs])
            oql_query += f" AND NeighborOf upstream (SELECT Entity WHERE id = ({targetids_str}))"

        if relation_types:
            rel_type_list = OQL.join_with_quotes( relation_types)
            oql_query = f" AND objectType = ({rel_type_list})"
        
        return oql_query


    @staticmethod
    def get_ppi(between_dbids1:set, and_dbids2:set):
        '''
        PPI relation types: Binding, DirectRegulation, ProtModification
        '''
        regulators = [str(dbid) for dbid in between_dbids1]
        targets = [str(dbid) for dbid in and_dbids2]
        reg_ids_str = ",".join(regulators)
        target_ids_str = ",".join(targets)
        oql_query = "SELECT Relation WHERE objectType = (Binding, DirectRegulation, ProtModification) "
        oql_query += f"AND NeighborOf (SELECT Entity WHERE id = ({reg_ids_str})) "
        oql_query += f"AND NeighborOf (SELECT Entity WHERE id = ({target_ids_str}))"
        return oql_query

def joinWithQuotes(separator, NameList:list):
    to_return = str()
    for i in range(0, len(NameList)):
        Name = str(NameList[i])
        if Name.find(' ') > 0 or Name.find('-'):
            to_return = to_return + '\'' + Name + '\'' + separator
        else:
            to_return = to_return + Name + separator
        
    return to_return[:len(to_return)-1]


def GetSearchStrings(PropertyNameList:list, PropValuesList:list):
    needQuotes = set([' ','-','/','(',')','[',']','+', '#'])
    some_values_have_quotes = False
    Values = str()
    unique_values_list = set(PropValuesList)
    for v in unique_values_list:
        val = str(v)
        val = val.replace('\'', '')
        if 1 in [c in val for c in needQuotes]:
            val = '\'' + val + '\''
            some_values_have_quotes = True
        Values = Values + val + ','
    Values = Values[:len(Values)-1]

    PropertyNames  = str()
    for n in range(0, len(PropertyNameList)):
        propName =  PropertyNameList[n]
        if propName.find(' ') > 0:
            if some_values_have_quotes == False:
                propName = '\"'+propName+'\"'
            else:
                print ("if you know how to concatenate string variable with double quotes and string variable with single quote in Python please let us know.\n"
                "Otherwise please remove values with white spaces from either your value list or property name list")
                return
        PropertyNames = PropertyNames+propName+','
    PropertyNames = PropertyNames[:len(PropertyNames)-1]

    return PropertyNames, Values


def GetEntitiesByProps(PropertyValues:list, SearchByProperties:list, OnlyObjectTypes=[]):
    OQLquery = str()
    if SearchByProperties[0] in ('id', 'Id', 'ID'):
        OQLquery = "SELECT Entity WHERE id = " +  '('+','.join([str(int) for int in PropertyValues]) +')'
    else:
        PropNames, Values = GetSearchStrings(SearchByProperties, PropertyValues)
        OQLquery = "SELECT Entity WHERE " + '('+PropNames+')' + " = " +  '('+Values+')'
    
    if len(OnlyObjectTypes) > 0:
        objectTypes = joinWithQuotes(',',OnlyObjectTypes)
        OQLquery = OQLquery + ' AND objectType = ('+objectTypes+')'   
    return OQLquery


def GetChildEntities(PropertyValues:list, SearchByProperties:list, OnlyObjectTypes=[]):
    OntologyQuery = 'SELECT Entity WHERE InOntology (SELECT Annotation WHERE Ontology=\'Pathway Studio Ontology\' AND Relationship=\'is-a\') under (SELECT OntologicalNode WHERE {entities})'
    EntityQuery = str()
    if SearchByProperties[0] in ('id', 'Id', 'ID'):
        EntityQuery = "id = " +  '('+','.join([str(int) for int in PropertyValues]) +')'
    else:
        PropNames, Values = GetSearchStrings(SearchByProperties, PropertyValues)
        EntityQuery = '('+PropNames+')' + " = " +  '('+Values+')'

    SearchQuery = OntologyQuery.format(entities=EntityQuery)
    if len(OnlyObjectTypes) > 0:
        objectTypes = joinWithQuotes(',',OnlyObjectTypes)
        SearchQuery = SearchQuery + ' AND objectType = ('+objectTypes+')'

    return SearchQuery


def ExpandEntityById(IDlist:list, ExpandWithRelationTypes:list=[], ExpandToNeighborTypes:list=[], direction=''):
    Values =  ','.join([str(int) for int in IDlist])
    expandALL =  'Select Relation WHERE NeighborOf ' + direction +' (SELECT Entity WHERE id = (' + Values + '))'
    if direction == 'upstream': oppositeDirection = 'downstream'
    elif direction == 'downstream': oppositeDirection = 'upstream'
    else: oppositeDirection = ''

    if len(ExpandWithRelationTypes) > 0:
        if len(ExpandToNeighborTypes) > 0:
            return expandALL + " AND objectType = (" + joinWithQuotes(',', ExpandWithRelationTypes) + ') AND NeighborOf '+oppositeDirection + ' (SELECT Entity WHERE objectType = (' +  joinWithQuotes(',', ExpandToNeighborTypes) + "))"
        else:
            return expandALL + " AND objectType = (" + joinWithQuotes(',', ExpandWithRelationTypes) + ")" 
    else:
        if len(ExpandToNeighborTypes) > 0:
            return expandALL + ' AND NeighborOf '+oppositeDirection + ' (SELECT Entity WHERE objectType = (' +  joinWithQuotes(',', ExpandToNeighborTypes) + "))"
        else:
            return expandALL


def ExpandEntity (PropertyValues:list, SearchByProperties:list, ExpandWithRelationTypes=[], ExpandToNeighborTypes=[], direction=''):
    if SearchByProperties[0] in ('id','Id','ID'):
        return ExpandEntityById(PropertyValues, ExpandWithRelationTypes, ExpandToNeighborTypes, direction)

    PropertyNames,Values = GetSearchStrings(SearchByProperties, PropertyValues)
    expandALL = 'Select Relation WHERE NeighborOf ' + direction + ' (SELECT Entity WHERE (' + PropertyNames + ') = ('+ Values+ '))'
    
    if direction == 'upstream': oppositeDirection = 'downstream'
    elif direction == 'downstream': oppositeDirection = 'upstream'
    else: oppositeDirection = ''

    if len(ExpandWithRelationTypes) > 0:
        if len(ExpandToNeighborTypes) > 0:
            return expandALL + " AND objectType = (" + joinWithQuotes(',', ExpandWithRelationTypes) + ') AND NeighborOf '+oppositeDirection+' (SELECT Entity WHERE objectType = (' +  joinWithQuotes(',', ExpandToNeighborTypes) + "))"
        else:
            return expandALL + " AND objectType = (" + joinWithQuotes(',', ExpandWithRelationTypes) + ")" 
    else:
        if len(ExpandToNeighborTypes) > 0:
            return expandALL + ' AND NeighborOf '+oppositeDirection + ' (SELECT Entity WHERE objectType = (' +  joinWithQuotes(',', ExpandToNeighborTypes) + "))"
        else:
            return expandALL
        
def GetNeighbors(PropertyValues:list, SearchByProperties:list, ExpandWithRelationTypes:list=[], ExpandToNeighborTypes:list=[]):
    PropertyNames,Values = GetSearchStrings(SearchByProperties, PropertyValues)
    connect_to_str = "to (SELECT Entity WHERE (" + PropertyNames + ") = " + Values
    if len(ExpandWithRelationTypes) > 0:
        if len(ExpandToNeighborTypes)>0:
            return "SELECT Entity WHERE objectType = ("+joinWithQuotes(',', ExpandToNeighborTypes)+") AND Connected by (SELECT Relation WHERE objectType= "+ joinWithQuotes(',', ExpandWithRelationTypes)+")" +connect_to_str +")"
        else:
            return "SELECT Entity WHERE Connected by (SELECT Relation WHERE objectType= "+ joinWithQuotes(',', ExpandWithRelationTypes)+")" +connect_to_str +")"
    else:
        if len(ExpandToNeighborTypes)>0:
            return "SELECT Entity WHERE objectType = ("+joinWithQuotes(',', ExpandToNeighborTypes)+") AND Connected by (SELECT Relation WHERE NOT (URN = NULL))"+ connect_to_str +")"
        else:
            # no ExpandWithRelationTypes specified and no ExpandToNeighborTypes specified -> get ALL neigbors
            return "SELECT Entity WHERE Connected by (SELECT Relation WHERE NOT (URN = NULL))" + connect_to_str +")"

def GetObjects(objIDlist:list):
    strings = [str(integer) for integer in objIDlist]
    IDlist = ",".join(strings)
    return "SELECT Entity WHERE id = ("+ IDlist + ")"

def GetDrugs(ForTargetsIDlist:list):
    strings = [str(integer) for integer in ForTargetsIDlist]
    IDlist = ",".join(strings)
    return "SELECT Relation WHERE objectType = (DirectRegulation,Binding) AND NeighborOf upstream (SELECT Entity WHERE id = ("+ IDlist + ")) AND NeighborOf downstream (SELECT Entity WHERE InOntology (SELECT Annotation WHERE Ontology='Pathway Studio Ontology' AND Relationship='is-a') under (SELECT OntologicalNode WHERE Name = drugs))"
    
def GetReaxysSubstances(ForTargetsIDlist:list):
    strings = [str(integer) for integer in ForTargetsIDlist]
    IDlist = ",".join(strings)
    return "SELECT Relation WHERE objectType = (DirectRegulation,Binding) AND NeighborOf upstream (SELECT Entity WHERE id = ("+ IDlist + ")) AND Source = Reaxys"

def ConnectEntities(PropertyValues1:list,SearchByProperties1:list,EntityTypes1:list,PropertyValues2:list,SearchByProperties2:list,EntityTypes2:list,ConnectByRelationTypes:list=[]):
    PropNames1, PropValues1 = GetSearchStrings(PropertyNameList=SearchByProperties1,PropValuesList=PropertyValues1)
    PropNames2, PropValues2 = GetSearchStrings(PropertyNameList=SearchByProperties2,PropValuesList=PropertyValues2)
    objectType1 = joinWithQuotes(',',EntityTypes1)
    objectType2 = joinWithQuotes(',',EntityTypes2)

    EntityQuery = 'SELECT Entity WHERE {entity}'

    Entity1 = '('+PropNames1+') = ('+PropValues1+') AND objectType = ('+objectType1+')'
    Entity2 = '('+PropNames2+') = ('+PropValues2+') AND objectType = ('+objectType2+')'
    #OntologyQuery ="SELECT Entity WHERE ({entity}) OR InOntology (SELECT Annotation WHERE Ontology='Pathway Studio Ontology' AND Relationship='is-a') under (SELECT OntologicalNode WHERE {entity})"
 
    Entity1Query = EntityQuery.format(entity = Entity1)
    Entity2Query = EntityQuery.format(entity = Entity2)

    if len(ConnectByRelationTypes) > 0:
        relTypeList = joinWithQuotes(',', ConnectByRelationTypes)
        OQLquery = "SELECT Relation WHERE objectType = ("+relTypeList+") AND NeighborOf ("+Entity1Query+") AND NeighborOf ("+Entity2Query+")"
        return OQLquery
    else:
        return "SELECT Relation WHERE NeighborOf ("+Entity1Query+") AND NeighborOf ("+Entity2Query+")"
 

def ConnectEntitiesIds(idlist1:list,idlist2:list,ConnectByRelTypes=[],RelEffect=[],RelDirection=''):
    stringIdlist1 = [str(integer) for integer in idlist1]
    stringIdlist2 = [str(integer) for integer in idlist2]
    Id1s = ",".join(stringIdlist1)
    Id2s = ",".join(stringIdlist2)
    EntityQuery = 'SELECT Entity WHERE id = ({idlist})'  
    Entity1Query = EntityQuery.format(idlist=Id1s)
    Entity2Query = EntityQuery.format(idlist=Id2s)

    OQLquery = "SELECT Relation WHERE NeighborOf {dir1} ("+Entity1Query+") AND NeighborOf {dir2} ("+Entity2Query+")"
    if len(ConnectByRelTypes) > 0:
        relTypeList = joinWithQuotes(',', ConnectByRelTypes)
        OQLquery = OQLquery + ' AND objectType = ('+relTypeList+')'
        
    if len(RelEffect) > 0:
        effectList = joinWithQuotes(',', RelEffect)
        OQLquery = OQLquery + ' AND Effect = ('+effectList+')'

    if RelDirection == '<':
        OQLquery = OQLquery.format(dir1='upstream',dir2='downstream')
    elif RelDirection == '>':
        OQLquery = OQLquery.format(dir1='downstream',dir2='upstream')
    else: OQLquery = OQLquery.format(dir1='',dir2='')
    
    return OQLquery



def FindTargets(RegulatorsIDs:list, TargetIDs:list, RelationTypeList:list=[]):
    regstrIDs = [str(integer) for integer in RegulatorsIDs]
    trgtstrIDs = [str(integer) for integer in TargetIDs]
    regIDs = ",".join(regstrIDs)
    trgtIDs = ",".join(trgtstrIDs)
    relTypeList = joinWithQuotes(',', RelationTypeList)
    OQLquery = "SELECT Relation WHERE objectType = ("+relTypeList+") AND NeighborOf downstream (SELECT Entity WHERE id = ("+ regIDs + ")) AND NeighborOf upstream (SELECT Entity WHERE id = ("+ trgtIDs + "))"
    return OQLquery

def GetPPIs(uniqObjIdList1:set, uniqObjIdList2:set):
    regstrIDs =  [str(integer) for integer in uniqObjIdList1]
    trgtstrIDs = [str(integer) for integer in uniqObjIdList2]
    regIDs = ",".join(regstrIDs)
    trgtIDs = ",".join(trgtstrIDs)
    OQLquery = "SELECT Relation WHERE objectType = (Binding, DirectRegulation, ProtModification) AND NeighborOf (SELECT Entity WHERE id = ("+ regIDs + ")) AND NeighborOf (SELECT Entity WHERE id = ("+ trgtIDs + "))"
    return OQLquery

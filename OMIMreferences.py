#Initializing Pathway Studio data model
import ElsevierAPI.ResnetAPI.PathwayStudioGOQL as OQL
from ElsevierAPI import ps_api


import csv
OMIMpairsFile = "OMIM/OMIMDisease-Gene_Demo.txt"
tsv_file = open(OMIMpairsFile)
OMIMpairs = csv.reader(tsv_file, delimiter="\t")

OmimProp = "OMIM relation"
ENTITY_PROPS = ['Name']
REL_PROPS = ['Name','Mechanism','RelationNumberOfReferences']


for pair in OMIMpairs:
    diseaseAliases  = [pair[0]]
    nodash = pair[0].replace('-',' ')
    diseaseAliases.append(nodash)
    if pair[0][len(pair[0])-5:] in [" type", " form"]:
        diseaseAliases.append(pair[0][:len(pair[0])-5])
    noS = pair[0].replace('\'s','')
    diseaseAliases.append(noS)
    diseaseAliases = list(set(diseaseAliases))
    gene = pair[1]
    FoundRelations = ps_api.connect_entities([gene], ["Name"], ['Protein', 'Complex', 'FunctionalClass'], diseaseAliases, ["Name", "Alias"], ['Disease'], REL_PROPS=REL_PROPS)
    if type(FoundRelations) == type(None):
        FoundRelation = ps_api.connect_entities([gene], ["Alias"], ['Protein', 'Complex', 'FunctionalClass'], diseaseAliases, ["Name", "Alias"], ['Disease'], REL_PROPS=REL_PROPS)
    if type(FoundRelations) != type(None):  
        print('Found relations for %s-%s pair' % (pair[0],gene))
        for regulatorID, targetID, rel in FoundRelations.edges.data('relation'):
            rel[OmimProp] = [pair[0]+' ---> '+gene]


ps_api.Graph.count_references()

fileOut = OMIMpairsFile[:len(OMIMpairsFile)-4] + ' references.tsv'
REL_PROPS.append(OmimProp)
ps_api.Graph.print_references(fileOut, REL_PROPS, ['Name'])
class PSObject(dict): #{PropId:[values], PropName:[values]}
    def __init__(self, ZeepObjectRef):
        zeepIter = iter(ZeepObjectRef)
        while True:
            try: item = next(zeepIter)
            except StopIteration: break  #Iterator exhausted: stop the loop
            else: self[item]= [ZeepObjectRef[item]]
                
    def __hash__(self):
        return self['Id'][0]
    
    def AddSingleProperty(self, PropId, PropValue:str):
        self[PropId] = [PropValue]

    def AddProperty(self, PropId, PropValue:str):
        try: list(self[PropId]).append(PropValue)
        except KeyError: self[PropId] = [PropValue]

    def AddUniqueProperty(self, PropId, PropValue:str):
        try: self[PropId] = list(set(self[PropId]+[PropValue]))
        except KeyError: self[PropId] = [PropValue]

    def AddUniquePropertyList(self, PropId, PropValues:list):
        try: self[PropId] = list(set(self[PropId]+PropValues))
        except KeyError: self[PropId] = list(set(PropValues))

    def PropValuesToStr(self, propID, cell_sep=';'):
        try: return cell_sep.join(self[propID])
        except KeyError: return ''
    
    def PropValuesToList(self, propID):
        try: return self[propID]
        except KeyError: return []

    def DataToStr (self, columnPropNames:list, col_sep='\t', cell_sep=';', endOfline='\n'):
        #assumes all properties in columnPropNames were fetched from Database otherwise will crash
        table_row = str()       
        for propName in columnPropNames:
            try:
                values = self[propName]
                propVal = cell_sep.join(values)
            except KeyError: propVal = ''
                
            table_row = table_row + propVal + col_sep

        return table_row[0:len(table_row)-1] + endOfline


REF_ID_TYPES = ['PMID','DOI','PII','PUI','EMBASE','NCT ID']
REF_PROPS = ['Sentence','PubYear','Authors','Journal','MedlineTA','CellType','CellLineName','Organ','Tissue','Organism','Source',
            'TrialStatus','Phase','StudyType','Start','Intervention','Condition','Company','Collaborator','TextRef']

class Reference(PSObject): #Identifiers{REF_ID_TYPES[i]:identifier}; self{REF_PROPS[i]:value}
    pass
    def __init__(self, idType:str, ID:str):
        self.Identifiers = {idType:ID}

    def __key(self):
        for id_type in REF_ID_TYPES:
            try: return self.Identifiers[id_type]
            except KeyError: continue
        return NotImplemented

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self,other):
        if isinstance(other, Reference):
            for id_type in REF_ID_TYPES:
                try: 
                    self_identifier = self.Identifiers[id_type]
                    try:
                        other_identifier = other.Identifiers[id_type]
                        return self_identifier == other_identifier
                    except KeyError: continue
                except KeyError: continue


    def ToString(self,IdType:str, sep='\t'):
        to_return = self.Identifiers[IdType]
        for propId, propValues in self.items():
            to_return = to_return+sep+propId+':'+';'.join(propValues)
        return to_return

    def MergeReference(self, ref_to_merge):
        self.Identifiers.update(ref_to_merge.Identifiers)

    def IsFromAbstract(self):
        for textref in self['TextRef']:
            try:
                return bool(textref.rindex('#abs',len(textref)-8, len(textref)-3))
            except ValueError: continue
        return False



class PSRelation(PSObject):
    pass
    def __init__(self, ZeepObjectRef):
        zeepIter = iter(ZeepObjectRef)
        while True:
            try:
                item = next(zeepIter)
            except StopIteration: break  #Iterator exhausted: stop the loop
            else:
                val = ZeepObjectRef[item]
                self[item]= [val]
        self.PropSetToProps = dict() #{PropSetID:{PropID:[values]}}
        self.Nodes = dict() #{"Regulators':[(entityID, Dir, effect)], "Targets':[(entityID, Dir, effect)]}
        self.References = dict() #{refIdentifier (PMID, DOI, EMBASE, PUI, LUI, Title):Reference}
       
    def __hash__(self):
        return self['Id'][0]
        
    def IsDirectional(self, Links):
        if len(Links[self['Id'][0]]) > 1: return True
        else: return False


    def PropValuesToStr(self, propID, cell_sep=';'):
        for id, values in self.items():
            if id == propID: return cell_sep.join(map(str,values))

        propSetVals = []
        for prop in self.PropSetToProps.values():            
            for id, values in prop.items():
                if id == propID:
                    propSetVals.append(cell_sep.join(map(str,values)))
                    break
             
        return cell_sep.join(propSetVals)


    def PropValuesToList(self, propID, cell_sep=';'):
        if propID in self.keys(): return self[propID]
        else:
            to_return = []
            for prop in self.PropSetToProps.values():
                propSetVal = str()
                for id, values in prop.items():
                    if id == propID:
                        propSetVal = cell_sep.join(values)
                        break
                to_return.append(propSetVal)
            return to_return


    def load_references(self):
        for propSet in self.PropSetToProps.values():
            propset_references = list()
            for ref_id_type in REF_ID_TYPES:
                try:
                    ref_id = propSet[ref_id_type][0]
                    propset_references.append((ref_id_type,ref_id))
                    #propset_article_identifiers.append(ref_id_value)
                except KeyError: continue

            if len(propset_references) == 0:
                try: 
                    propsetTtitle = propSet['Title'][0] #trying id reference by title as a last resort since it does not have valid identifiers 
                    try:
                        Ref =  self.References[propsetTtitle]
                    except KeyError: 
                        Ref = Reference('Title', propsetTtitle)
                        self.References[propsetTtitle] = Ref
                except KeyError: continue
            else:
                #case when reference have valid identifiers
                propset_article_identifiers = {x[1] for x in propset_references}
                existing_ref = {r for i,r in self.References.items() if i in propset_article_identifiers}
                if len(existing_ref) == 0:#case when reference is new
                    id_type, id_value = propset_references[0][0], propset_references[0][1]
                    #id_value = propset_references[0][1]
                    Ref = Reference(id_type,id_value)
                    self.References[id_value] = Ref
                    for Id in range (1, len(propset_references)):
                        id_type, id_value = propset_references[Id][0], propset_references[Id][1]
                        #id_value = propset_references[Id][1]
                        Ref.Identifiers[id_type] = id_value
                        self.References[id_value] = Ref
                elif len(existing_ref) > 1:#identifiers from one propset point to different references
                    #will merge all references from propset with the first one 
                    conflictRefs = list(existing_ref)
                    anchorRef = conflictRefs[0]
                    for i in range(1,len(conflictRefs)):
                        ref_to_merge = conflictRefs[i]
                        anchorRef.MergeReference(ref_to_merge)
                        for id_value in ref_to_merge.Identifiers.values():
                            self.References[id_value] = anchorRef
                else:
                    Ref = next(iter(existing_ref))#when existing_ref == 1: there is nothing to do

            for propId, propValues in propSet.items():#adding all other valid properties to Ref
                if propId in REF_PROPS: 
                    Ref.AddUniquePropertyList(propId, propValues)

    def GetReferenceCount(self, count_abstracts=False): 
        if count_abstracts == True:
            ref_from_abstract = set([x for x in self.References.values() if x.IsFromAbstract() == True])
            return len(ref_from_abstract)
        else: return len(self.References)

    def TripleToStr(self,columnPropNames:list,return_dict = False,col_sep='\t',cell_sep=';',endOfline='\n',RefNumPrintLimit=0):
        #assumes all properties in columnPropNames were fetched from Database otherwise will crash
        #initializing table
        colCount = len(columnPropNames)+2
        RelationNumberOfReferences = int(self['RelationNumberOfReferences'][0])

        rowCount = 1
        if RelationNumberOfReferences >= RefNumPrintLimit:
            rowCount = max(1, len(self.PropSetToProps))

        regulatorIDs = str()
        targetIDs = str()
        for k,v in self.Nodes.items():
            if k == 'Regulators':
                regIDs = [x[0] for x in v]
                regulatorIDs = ','.join(map(str,regIDs))
            else:
                trgtIDs = [x[0] for x in v]
                targetIDs = ','.join(map(str,trgtIDs))


        first_row = [''] * colCount
        first_row[colCount-2] = regulatorIDs 
        first_row[colCount-1] = targetIDs 
        table = {0: first_row}

        for r in range(1,rowCount):
            new_row = [''] * colCount
            new_row[colCount-2] = regulatorIDs 
            new_row[colCount-1] = targetIDs 
            table[r] = new_row

        
        for col in range (len(columnPropNames)):
            propId = columnPropNames[col]
            if propId in self.keys():
                for row in range(0, rowCount):
                    propValue = self.PropValuesToStr(propId)
                    table[row][col] = propValue
            elif RelationNumberOfReferences >= RefNumPrintLimit:
                row = 0
                for propList in self.PropSetToProps.values():
                    if propId in propList.keys():
                        propValues = propList[propId]
                        cellValue = cell_sep.join(propValues)
                        table[row][col] = cellValue
                        row +=1

        if return_dict == True: return table
        else:
            tableStr = str()
            for row in table.values():
                tableStr = tableStr + col_sep.join(row) + endOfline
            return tableStr

    def GetRegulatorTargetPairs(self):
        #relEntities = self.Links[relId]
        if len(self.Nodes) > 1:
            RegTargetPairs = []
            for regTriple in self.Nodes['Regulators']:
                for targetTriple in self.Nodes['Targets']:
                    pairTuple = (regTriple[0],targetTriple[0])
                    RegTargetPairs.append(pairTuple)
            return RegTargetPairs
        else:
            import itertools
            try:
                objIdList = [x[0] for x in self.Nodes['Regulators']]
                return itertools.combinations(objIdList, 2)
            except KeyError: 
                objIdList = [x[0] for x in self.Nodes['Targets']]
                return itertools.combinations(objIdList, 2)

    def GetEntitiesIDs(self):
        nodeIds = [x[0] for x in self.Nodes['Regulators']] + [x[0] for x in self.Nodes['Targets']]
        return list(set(nodeIds))

    def toJSON(self):
        import json
        str1 = '{"Relation Properties": ' + json.dumps(self) + '}'
        strP = '{"Relation References": ' + json.dumps(self.PropSetToProps) + '}'
        strR = json.dumps(self.Nodes)
        return str1 +'\n'+ strP +'\n'+strR +'\n'

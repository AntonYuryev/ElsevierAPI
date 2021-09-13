import pandas as pd
import re

class PSObject(dict):  # {PropId:[values], PropName:[values]}
    def __init__(self, dic: dict):
        super().__init__(dic)

    @staticmethod
    def from_zeep_objects(ZeepObjectRef):
        return_dict = dict()
        zeep_iter = iter(ZeepObjectRef)
        while True:
            try:
                item = next(zeep_iter)
            except StopIteration:
                break  # Iterator exhausted: stop the loop
            else:
                return_dict[item] = [ZeepObjectRef[item]]
        return return_dict

    @classmethod
    def from_zeep(cls, ZeepObjectRef):
        return cls(cls.from_zeep_objects(ZeepObjectRef))

    def __hash__(self):
        return self['Id'][0]

    def add_single_property(self, PropId, PropValue: str):
        self[PropId] = [PropValue]

    def add_property(self, PropId, PropValue):
        try:
            self[PropId].append(PropValue)
        except KeyError:
            self[PropId] = [PropValue]

    def add_unique_property(self, PropId, PropValue: str):
        try:
            self[PropId] = list(set(self[PropId] + [PropValue]))
        except KeyError:
            self[PropId] = [PropValue]

    def add_properties(self, PropId, PropValues: list):
        try:
            self[PropId] = list(set(self[PropId] + PropValues))
        except KeyError:
            self[PropId] = list(set(PropValues))

    def prop_values2str(self, propID, cell_sep=';'):
        try:
            return cell_sep.join(self[propID])
        except KeyError:
            return ''

    def prop_values2list(self, propID):
        try:
            return self[propID]
        except KeyError:
            return []

    def data2str(self, columnPropNames: list, col_sep='\t', cell_sep=';', endOfline='\n'):
        # assumes all properties in columnPropNames were fetched from Database otherwise will crash
        table_row = str()
        for propName in columnPropNames:
            try:
                values = self[propName]
                prop_val = cell_sep.join(values)
            except KeyError:
                prop_val = ''

            table_row = table_row + prop_val + col_sep

        return table_row[0:len(table_row) - 1] + endOfline

    def to_pandas(self, columnPropNames: list):
        # assumes all properties in columnPropNames were fetched from Database otherwise will crash
        properties_pandas = pd.DataFrame(columns = columnPropNames)
        for propName in columnPropNames:
            try:
                # values = self[propName]
                # prop_val = ';'.join(values)
                properties_pandas.at[0,propName] = ';'.join(self[propName])
            except KeyError:
                continue

        return properties_pandas

    def has_str_property(self, PropName, PropValues: list, case_insensitive=True):
        if case_insensitive:
            search_set = {x.lower() for x in PropValues}
            try:
                props = [p.lower() for p in self[PropName]]
                return not search_set.isdisjoint(props)
            except KeyError: return False
        else:
            search_set = set(PropValues)
            try:
                return not search_set.isdisjoint(self[PropName])
            except KeyError: return False



REF_ID_TYPES = ['PMID', 'DOI', 'PII', 'PUI', 'EMBASE', 'NCT ID']
REF_PROPS = ['Sentence', 'PubYear', 'Authors', 'Journal', 'MedlineTA', 'CellType', 'CellLineName', 'Organ', 'Tissue',
             'Organism', 'Source',
             'TrialStatus', 'Phase', 'StudyType', 'Start', 'Intervention', 'Condition', 'Company', 'Collaborator',
             'TextRef', 'Title']

NOT_ALLOWED_IN_SENTENCE='[\t\r\n\v\f]' # regex to clean up special characters in sentences


class Reference(PSObject):  # Identifiers{REF_ID_TYPES[i]:identifier}; self{REF_PROPS[i]:value}
    pass

    def __init__(self, idType: str, ID: str):
        super().__init__(dict())
        self.Identifiers = {idType:ID}

    def __key(self):
        for id_type in REF_ID_TYPES:
            try:
                return self.Identifiers[id_type]
            except KeyError:
                continue
        return NotImplemented

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other):
        if isinstance(other, Reference):
            for id_type in REF_ID_TYPES:
                try:
                    return self.Identifiers[id_type] == other.Identifiers[id_type]
                except KeyError:
                    continue

    def to_string(self, id_types: list=None, sep='\t'):
        id_types = id_types if isinstance(id_types,list) else ['PMID']
        row = str()
        for t in id_types:
            try:
                row = row + t+':'+self.Identifiers[t]+sep
            except KeyError:
                row = row + sep

        row = row[:len(row)-1] #remove last separator
        for prop_id, prop_values in self.items():
            prop_val = ';'.join(prop_values)
            if prop_id in ['Sentence','Title']:
                prop_val = re.sub(NOT_ALLOWED_IN_SENTENCE, ' ', prop_val)
            row = row + sep + prop_id + ':' + prop_val
        return row


    def merge_reference(self, ref_to_merge):
        self.Identifiers.update(ref_to_merge.Identifiers)

    def is_from_abstract(self):
        for textref in self['TextRef']:
            try:
                return bool(textref.rindex('#abs', len(textref) - 8, len(textref) - 3))
            except ValueError:
                continue
        return False


class PSRelation(PSObject):
    pass
    PropSetToProps = dict()
    Nodes = dict()
    References = dict()

    def __init__(self, dic: dict):
        super().__init__(dic)
        self.PropSetToProps = dict()  # {PropSetID:{PropID:[values]}}
        self.Nodes = dict()  # {"Regulators':[(entityID, Dir, effect)], "Targets':[(entityID, Dir, effect)]}
        self.References = dict()  # {refIdentifier (PMID, DOI, EMBASE, PUI, LUI, Title):Reference}

    def copy(self, other):
        if isinstance(other, PSRelation):
            self.update(other)
            self.PropSetToProps.update(other.PropSetToProps)
            self.Nodes.update(other.Nodes)
            self.References.update(other.References)

    def __hash__(self):
        return self['Id'][0]

    def is_directional(self, Links):
        if len(Links[self['Id'][0]]) > 1:
            return True
        else:
            return False

    def prop_values2str(self, prop_id, cell_sep: str =';'):
        try:
            return cell_sep.join(map(str, self[prop_id]))
        except KeyError:
            prop_set_values = []
            for prop in self.PropSetToProps.values():
                try:
                    prop_set_values.append(cell_sep.join(map(str, prop[prop_id])))
                except KeyError:
                    continue
            
            to_return = cell_sep.join(prop_set_values)
            if prop_id in ['Sentence','Title']:
                to_return = re.sub(NOT_ALLOWED_IN_SENTENCE, ' ', to_return)
            return to_return

    def prop_values2list(self, propID, cell_sep=';'):
        try:
            return self[propID]
        except KeyError:
            to_return = list()
            for prop in self.PropSetToProps.values():
                prop_set_val = str()
                for prop_id, values in prop.items():
                    if prop_id == propID:
                        prop_set_val = cell_sep.join(values)
                        break
                to_return.append(prop_set_val)
            return to_return

    def load_references(self):
        if len(self.References): return
        for propSet in self.PropSetToProps.values():
            propset_references = list()
            for ref_id_type in REF_ID_TYPES:
                try:
                    ref_id = propSet[ref_id_type][0]
                    propset_references.append((ref_id_type, ref_id))
                except KeyError:
                    continue

            if len(propset_references) == 0:
                try:
                    # trying id reference by title as a last resort since it does not have valid identifiers
                    propset_title = propSet['Title'][0]
                    try:
                        ref = self.References[propset_title]
                    except KeyError:
                        ref = Reference('Title', propset_title)
                        self.References[propset_title] = ref
                except KeyError:
                    continue
            else:
                # case when reference have valid identifiers
                propset_article_identifiers = {x[1] for x in propset_references}
                existing_ref = {r for i, r in self.References.items() if i in propset_article_identifiers}
                if len(existing_ref) == 0:  # case when reference is new
                    id_type, id_value = propset_references[0][0], propset_references[0][1]
                    # id_value = propset_references[0][1]
                    ref = Reference(id_type, id_value)
                    self.References[id_value] = ref
                    for Id in range(1, len(propset_references)):
                        id_type, id_value = propset_references[Id][0], propset_references[Id][1]
                        ref.Identifiers[id_type] = id_value
                        self.References[id_value] = ref
                elif len(existing_ref) > 1:  # identifiers from one propset point to different references
                    # will merge all references from propset with the first one
                    conflict_refs = list(existing_ref)
                    anchor_ref = conflict_refs[0]
                    for i in range(1, len(conflict_refs)):
                        ref_to_merge = conflict_refs[i]
                        anchor_ref.merge_reference(ref_to_merge)
                        for id_value in ref_to_merge.Identifiers.values():
                            self.References[id_value] = anchor_ref
                        ref = anchor_ref
                else:
                    ref = next(iter(existing_ref))  # if len(existing_ref) == 1: there is nothing to do

            for propId, propValues in propSet.items():  # adding all other valid properties to Ref
                if propId in REF_PROPS:
                    ref.add_properties(propId, propValues)

    def get_reference_count(self, count_abstracts=False):
        if count_abstracts:
            ref_from_abstract = set([x for x in self.References.values() if x.is_from_abstract()])
            return len(ref_from_abstract)
        else:
            return len(self.References)

    def triple2str(self, columnPropNames: list, return_dict=False, col_sep='\t', cell_sep=';', endOfline='\n',
                   RefNumPrintLimit=0,add_entities=False):
        # assumes all properties in columnPropNames were fetched from Database otherwise will crash
        # initializing table
        col_count = len(columnPropNames) +2 if add_entities else len(columnPropNames)
        RelationNumberOfReferences = int(self['RelationNumberOfReferences'][0])

        rowCount = 1
        if RelationNumberOfReferences >= RefNumPrintLimit:
            rowCount = max(1, len(self.PropSetToProps))

        table = dict()
        for r in range(0, rowCount):
            table[r] = [''] * col_count

        if add_entities:
            regulatorIDs = str()
            targetIDs = str()
            for k, v in self.Nodes.items():
                if k == 'Regulators':
                    regulatorIDs = ','.join(map(str, [x[0] for x in v]))
                else:
                    targetIDs = ','.join(map(str, [x[0] for x in v]))

            for r in range(0, rowCount):
                table[r][col_count-2] = regulatorIDs
                table[r][col_count-1] = targetIDs

        for col in range(len(columnPropNames)):
            propId = columnPropNames[col]
            if propId in self.keys(): # filling columns with relation properties
                for row in range(0, rowCount):
                    propValue = self.prop_values2str(propId)
                    table[row][col] = propValue
            elif RelationNumberOfReferences >= RefNumPrintLimit: #filling columns with reference properties
                row = 0
                for propList in self.PropSetToProps.values():
                    if propId in propList.keys():
                        propValues = propList[propId]
                        cellValue = cell_sep.join(propValues)
                        table[row][col] = cellValue
                    row += 1

        if return_dict: return table
        else:
            tableStr = str()
            for row in table.values():
                tableStr = tableStr + col_sep.join(row) + endOfline
            return tableStr

    def to_pandas(self, columnPropNames: list, RefNumPrintLimit=0):
        # assumes all properties in columnPropNames were fetched from Database otherwise will crash
        # initializing table
        col_count = len(columnPropNames) + 2
        RelationNumberOfReferences = int(self['RelationNumberOfReferences'][0])

        rowCount = 1
        if RelationNumberOfReferences >= RefNumPrintLimit:
            rowCount = max(1, len(self.PropSetToProps))

        regulatorIDs = str()
        targetIDs = str()
        for k, v in self.Nodes.items():
            if k == 'Regulators':
                regulatorIDs = ','.join(map(str, [x[0] for x in v]))
            else:
                targetIDs = ','.join(map(str, [x[0] for x in v]))

        
        reference_pandas = pd.DataFrame(columns = columnPropNames+["Regulators Id","Targets Id"])
        for row in range(0, rowCount):
            reference_pandas.at[row,"Regulators Id"] = regulatorIDs
            reference_pandas.at[row,"Targets Id"] = targetIDs

        for col in columnPropNames:
            if col in self.keys():
                reference_pandas[col] = self.prop_values2str(col)
            elif RelationNumberOfReferences >= RefNumPrintLimit:
                row = 0
                for propList in self.PropSetToProps.values():
                    if col in propList.keys():
                        propValues = propList[col]
                        cellValue = ';'.join(propValues)
                        reference_pandas.at[row,col] = cellValue
                    row += 1
                
        return reference_pandas

    def get_regulators_targets(self):
        # relEntities = self.Links[relId]
        if len(self.Nodes) > 1:
            RegTargetPairs = []
            for regTriple in self.Nodes['Regulators']:
                for targetTriple in self.Nodes['Targets']:
                    pairTuple = (regTriple[0], targetTriple[0])
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

    def get_entities_ids(self):
        nodeIds = {[x[0] for x in self.Nodes['Regulators']] + [x[0] for x in self.Nodes['Targets']]}
        return list(nodeIds)

    def to_json(self):
        import json
        str1 = '{"Relation Properties": ' + json.dumps(self) + '}'
        strP = '{"Relation References": ' + json.dumps(self.PropSetToProps) + '}'
        strR = json.dumps(self.Nodes)
        return str1 + '\n' + strP + '\n' + strR + '\n'

from ElsevierAPI.ResnetAPI.NetworkxObjects import PSObject
import pandas as pd
import logging
import sys
import ElsevierAPI.ResnetAPI.PathwayStudioGOQL as OQL
import zeep


def configure_logging(logger):
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger

class DataModel:
    def __init__(self, url, username, password):
        self.IdToPropType = dict()
        self.IdtoObjectType = dict()
        self.PropIdToDict = dict()
        #self.id2folders = dict() # {id:folder_zobj, name:folder_zobj} loaded on demand
        #self.id2pathways = dict() # {id:PSObj} loaded on demand
        self.RNEFnameToPropType = dict()
        from requests.auth import HTTPBasicAuth
        from requests import Session
        from zeep.cache import SqliteCache
        from zeep.transports import Transport
        session = Session()
        session.auth = HTTPBasicAuth(username, password)
        transport = Transport(cache=SqliteCache(), session=session)
        from zeep import Client, Settings
        settings = Settings(strict=False, xml_huge_tree=True)
        #settings = zeep.Settings(extra_http_headers={'Authorization': 'Bearer ' + token})
        self.logger = configure_logging(logging.getLogger(__name__))
        try:
            self.SOAPclient = Client(wsdl=url, transport=transport, settings=settings)
            self.__load_model()
            print('Connected to Resnet API server:\n%s' % url)
        except Exception as error:
            self.logger.error("Pathway Studio server connection failed: {error}".format(error=error))
            raise ConnectionError(f"Server connection failed. Wrong or inaccessible url: {url}") from None

    def __load_model(self):
        object_types = self.SOAPclient.service.GetObjectTypes()
        property_types = self.SOAPclient.service.GetPropertyDefinitions()
    
        for i in range(0, len(object_types)):
            db_id = object_types[i]['Id']
            obj_type_name = object_types[i]['Name']
            obj_type_display_name = object_types[i]['DisplayName']
            self.IdtoObjectType[obj_type_name] = object_types[i]
            self.IdtoObjectType[obj_type_display_name] = object_types[i]
            self.IdtoObjectType[db_id] = object_types[i]
            for s in object_types[i]['Synonyms']:
                self.IdtoObjectType[s] = object_types[i]

        for i in range(0, len(property_types)):
            db_id = property_types[i]['Id']
            prop_type_name = property_types[i]['Name']
            prop_type_display_name = property_types[i]['DisplayName']
            self.IdToPropType[prop_type_name] = property_types[i]
            self.IdToPropType[prop_type_display_name] = property_types[i]
            self.IdToPropType[db_id] = property_types[i]
            self.RNEFnameToPropType[prop_type_name] = property_types[i]

    def __object_types_by_class_id(self, classID):
        object_types = list()
        zeep_object_types = self.SOAPclient.service.GetObjectTypes()
        for i in range(0, len(zeep_object_types)):
            obj_type_name = zeep_object_types[i]['Name']
            obj_type_class_id = zeep_object_types[i]['ObjClassId']
            if obj_type_class_id == classID:
                object_types.append(obj_type_name)

        return object_types

    def get_relation_types(self):
        return self.__object_types_by_class_id(3)

    def get_entity_types(self):
        return self.__object_types_by_class_id(1)

    def load_folder_tree(self):
        folders = self.SOAPclient.service.GetFoldersTree(0)
        IdToFolders = dict()
        for folder in folders:
            folder_id = folder['Id']
            folder_name = folder['Name']
            IdToFolders[folder_id] = [folder]
            # several folders may have the same name:
            if folder_name in IdToFolders.keys():
                IdToFolders[folder_name].append(folder)
            else:
                IdToFolders[folder_name] = [folder]
        
        return IdToFolders

    def dump_prop_names(self, fileOut):
        id_prop_names_list = list({(p['DisplayName'], p['Name'], p['Id']) for p in self.IdToPropType.values()})
        id_prop_names_list.sort(key=lambda x: x[0])
        import json
        with open(fileOut, 'w', encoding='utf-8') as f:
            f.write("Format - [\'Name in PS UI\' ,\'Search string for GOQL\',\'PropertyType Id\']: \n")
            for p in id_prop_names_list:
                f.write(json.dumps(p) + '\n')

    def dump_object_names(self, fileOut):
        id_obj_type_list = list({(p['DisplayName'], p['Name'], p['Id']) for p in self.IdtoObjectType.values()})
        id_obj_type_list.sort(key=lambda x: x[0])
        import json
        with open(fileOut, 'w', encoding='utf-8') as f:
            f.write("Format - [\'Name in PS UI\' ,\'Search string for GOQL\',\'ObjectType Id\']: \n")
            for o in id_obj_type_list:
                f.write(json.dumps(o) + '\n')

    def map_property_names_to_id(self, PropertyNames: list):
        id_list = [self.IdToPropType[x]['Id'] for x in PropertyNames]
        if 'Name' not in PropertyNames:
            id_list.append(self.IdToPropType['Name']['Id'])
        return id_list

    def get_dictionary(self, idProperty, idDictFolder: int):
        if idProperty not in self.PropIdToDict.keys():
            dict_folder = self.SOAPclient.service.GetDictFolder(idDictFolder)
            id_values_to_str = dict()
            for val in dict_folder.Values.DictValue:
                db_id = val['Id']
                val = val['Value']
                id_values_to_str[db_id] = val
            self.PropIdToDict[idProperty] = id_values_to_str
            self.PropIdToDict[dict_folder['Name']] = id_values_to_str
        return self.PropIdToDict[idProperty]

    def get_dict_value(self, IdProperty: int, DictIdValue: int):
        if IdProperty not in self.PropIdToDict[IdProperty].keys():
            id_dictionary = self.IdToPropType[IdProperty]['DictFolderID']
            assert id_dictionary > 0
            self.get_dictionary(IdProperty, id_dictionary)

        return self.PropIdToDict[IdProperty][DictIdValue]

    def get_subfolders(self, FolderIds: list):
        if not hasattr(self, "IdToFolders"): 
            self.id2folders = self.load_folder_tree()

        subfolders_ids = set()
        for db_id in FolderIds:
            folders = self.id2folders[db_id]
            for folder in folders:
                if type(folder['SubFolders']) != type(None):
                    subfolders_ids.update(folder['SubFolders']['long'])
        return subfolders_ids

    def get_subfolders_recursively(self, FolderId):
        accumulate_subfolders = {FolderId}
        subfolders = set(self.get_subfolders([FolderId]))
        while len(subfolders) > 0:
            accumulate_subfolders.update(subfolders)
            subs = set()
            for db_id in subfolders:
                subs.update(self.get_subfolders([db_id]))
            subfolders = subs

        return accumulate_subfolders

    def get_subfolder_tree(self,folder_name):
        if not hasattr(self, "id2folders"): 
            self.id2folders = self.load_folder_tree()

        folder_id = self.id2folders[folder_name][0]['Id']
        child2parent = {folder_id:folder_id}
        parent2child = PSObject(dict())
        subfolder_ids = self.get_subfolders([folder_id])
        while len(subfolder_ids) > 0:
            child2parent.update({i:folder_id for i in subfolder_ids})
            parent2child.add_properties(folder_id,list(subfolder_ids))
            subs = set()
            for db_id in subfolder_ids:
                subs = self.get_subfolders([db_id])
                if len(subs) > 0:
                    child2parent.update({i:db_id for i in subs})
                    parent2child.add_properties(db_id,subs)
            subfolder_ids = subs

        return child2parent, dict(parent2child)

    def get_folder_objects(self, FolderId, result_param):
        from zeep import exceptions
        for i in range (0,3):
            try:
                result = self.SOAPclient.service.FolderGetObjects(FolderId, result_param)
                return result
            except exceptions.Fault: continue

    def oql_response(self, OQLquery, result_param):
        from zeep import exceptions
        for i in range (0,3):
            try:
                result = self.SOAPclient.service.OQLSearch(OQLquery, result_param)
                return result
            except exceptions.Fault: continue


    def result_get_data(self, result_param):
        result = self.SOAPclient.service.ResultGetData(result_param)
        return result

    def create_result_param(self, property_names=None):
        property_names = ['Name'] if property_names is None else property_names

        id_property_list = self.map_property_names_to_id(property_names)
        prop_refs = list()
        prop_ref = self.SOAPclient.get_type('ns0:PropertyRef')
        for db_id in id_property_list:
            pl = prop_ref(Id=db_id)
            prop_refs.append(pl)
        result_param = self.SOAPclient.get_type('ns0:ResultParam')
        rp = result_param(
            CreateResult=False,
            ResultRef='?',
            Objects=[],
            Objects_size=0,
            ResultPos=0,
            MaxPageSize=0,
            GetObjects=False,
            GetProperties=False,
            PropertyList={'PropertyRef': prop_refs},
            PropertyList_Size=len(property_names),
            GetLinks=False,
            GetParents=False,
            GetMembers=False,
            GetAddlCol=False,
            RefLimit=0,
            SortColumn=0,
            SortDescending=False,
            SortPropId=0,
            AddlAttrs=[],
            AddlAttrs_Size=0,
            ApplySourceFilter=False)
        return rp

    def get_object_properties(self, obj_ids: list, property_names: list=None):
        property_names = ['Name'] if property_names is None else property_names

        rp = self.create_result_param(property_names)
        rp.GetObjects = True
        rp.GetProperties = True
        oql_request = OQL.get_objects(obj_ids)
        obj_props = self.oql_response(oql_request, rp)
        # setting objectType name
        for obj in obj_props.Objects.ObjectRef:
            obj_type_id = obj['ObjTypeId']
            obj_type_name = self.IdtoObjectType[obj_type_id]['Name']
            obj['ObjTypeName'] = obj_type_name

        # setting dict property values and property names
        for prop in obj_props.Properties.ObjectProperty:
            id_property = prop['PropId']
            prop_name = self.IdToPropType[id_property]['Name']
            prop_display_name = self.IdToPropType[id_property]['DisplayName']
            prop['PropName'] = prop_name
            prop['PropDisplayName'] = prop_display_name
            dict_folder_id = self.IdToPropType[id_property]['DictFolderId']
            if dict_folder_id > 0:
                dict_folder = self.get_dictionary(id_property, dict_folder_id)
                for i in range(0, len(prop['PropValues']['string'])):
                    id_dict_prop_value = int(prop['PropValues']['string'][i])
                    new_dict_value = dict_folder[id_dict_prop_value]
                    prop['PropValues']['string'][i] = new_dict_value

        return obj_props

    def get_folder_objects_props(self, FolderId, property_names=None):
        property_names = ['Name'] if property_names is None else property_names

        rp = self.create_result_param(property_names)
        rp.GetObjects = True
        rp.GetProperties = True
        obj_props = self.get_folder_objects(FolderId, rp)
        if type(obj_props.Objects) == type(None):return None

        #identifying synlinks 
        id2symlink = dict()
        for col in obj_props.AdditionalColumns.ResultAdditionalColumn:
            obj_id = col['ObjId']
            id2symlink[obj_id] = True if col['Value'] == '0' else False


        # setting objectType name
        for obj in obj_props.Objects.ObjectRef:
            obj_type_id = obj['ObjTypeId']
            obj_type_name = self.IdtoObjectType[obj_type_id]['Name']
            obj['ObjTypeName'] = obj_type_name
            obj['IsSymlink'] = id2symlink[obj['Id']]

        # setting dict property values and property names
        for prop in obj_props.Properties.ObjectProperty:
            id_property = prop['PropId']
            prop_name = self.IdToPropType[id_property]['Name']
            prop['PropName'] = prop_name
            dict_folder_id = self.IdToPropType[id_property]['DictFolderId']
            if dict_folder_id > 0:
                dict_folder = self.get_dictionary(id_property, dict_folder_id)
                for i in range(0, len(prop['PropValues']['string'])):
                    id_dict_prop_value = int(prop['PropValues']['string'][i])
                    new_dict_value = dict_folder[id_dict_prop_value]
                    prop['PropValues']['string'][i] = new_dict_value

        


        return obj_props

    def get_layout(self, PathwayId):
        # GetObjectAttachment = self.SOAPclient.get_type('ns0:GetObjectAttachment')
        result = self.SOAPclient.service.GetObjectAttachment(PathwayId, 1)
        return str(result['Attachment'].decode('utf-8'))

    def get_data(self, OQLrequest, retrieve_props: list=None, getLinks=True):
        retrieve_props = ['Name', 'RelationNumberOfReferences'] if retrieve_props is None else retrieve_props

        rp = self.create_result_param(retrieve_props)
        rp.GetObjects = True
        rp.GetProperties = True
        rp.GetLinks = getLinks
        # setting objectType name
        obj_props = self.oql_response(OQLrequest, rp)
        if type(obj_props.Objects) == type(None):
            # print('Your SOAP response is empty! Check your OQL query and try again\n')
            return None

        if len(obj_props.Objects.ObjectRef) == 0:
            return None

        for obj in obj_props.Objects.ObjectRef:
            obj_type_id = obj['ObjTypeId']
            obj_type_name = self.IdtoObjectType[obj_type_id]['Name']
            obj['ObjTypeName'] = obj_type_name
        # setting dict property values
        for prop in obj_props.Properties.ObjectProperty:
            id_property = prop['PropId']
            prop_name = self.IdToPropType[id_property]['Name']
            prop_display_name = self.IdToPropType[id_property]['DisplayName']
            prop['PropName'] = prop_name
            prop['PropDisplayName'] = prop_display_name
            dict_folder_id = self.IdToPropType[id_property]['DictFolderId']
            if dict_folder_id > 0:
                dict_folder = self.get_dictionary(id_property, dict_folder_id)
                for i in range(0, len(prop['PropValues']['string'])):
                    id_dict_prop_value = int(prop['PropValues']['string'][i])
                    new_dict_value = dict_folder[id_dict_prop_value]
                    prop['PropValues']['string'][i] = new_dict_value

        return obj_props

    def init_session(self, OQLrequest, PageSize: int, property_names=None, getLinks=True):
        property_names = ['Name'] if property_names is None else property_names

        rp = self.create_result_param(property_names)
        rp.GetObjects = True
        rp.GetProperties = True
        rp.GetLinks = getLinks
        rp.CreateResult = True
        rp.MaxPageSize = PageSize
        obj_props = self.oql_response(OQLrequest, rp)
        if type(obj_props.Objects) == type(None):
            # print('Your SOAP response is empty! Check your OQL query and try again\n')
            return None, (obj_props.ResultRef, obj_props.ResultSize, obj_props.ResultPos)

        for obj in obj_props.Objects.ObjectRef:
            obj_type_id = obj['ObjTypeId']
            obj_type_name = self.IdtoObjectType[obj_type_id]['Name']
            obj['ObjTypeName'] = obj_type_name
        # setting dict property values
        for prop in obj_props.Properties.ObjectProperty:
            id_property = prop['PropId']
            prop_name = self.IdToPropType[id_property]['Name']
            prop_display_name = self.IdToPropType[id_property]['DisplayName']
            prop['PropName'] = prop_name
            prop['PropDisplayName'] = prop_display_name
            dict_folder_id = self.IdToPropType[id_property]['DictFolderId']
            if dict_folder_id > 0:
                dict_folder = self.get_dictionary(id_property, dict_folder_id)
                for i in range(0, len(prop['PropValues']['string'])):
                    id_dict_prop_value = int(prop['PropValues']['string'][i])
                    new_dict_value = dict_folder[id_dict_prop_value]
                    prop['PropValues']['string'][i] = new_dict_value

        return obj_props, (obj_props.ResultRef, obj_props.ResultSize, obj_props.ResultPos)

    def get_session_page(self, ResultRef, ResultPos, PageSize, ResultSize, property_names=None, getLinks=True):
        property_names = ['Name'] if property_names is None else property_names

        rp = self.create_result_param(property_names)
        rp.GetObjects = True
        rp.GetProperties = True
        rp.GetLinks = getLinks

        rp.ResultSize = ResultSize
        rp.ResultRef = ResultRef
        rp.MaxPageSize = PageSize
        rp.ResultPos = ResultPos
        obj_props = self.result_get_data(rp)

        # setting objectType name
        if type(obj_props.Objects) == type(None):
            # print('Your SOAP response is empty! Check your OQL query and try again\n')
            return

        for obj in obj_props.Objects.ObjectRef:
            obj_type_id = obj['ObjTypeId']
            obj_type_name = self.IdtoObjectType[obj_type_id]['Name']
            obj['ObjTypeName'] = obj_type_name
        # setting dict property values
        for prop in obj_props.Properties.ObjectProperty:
            id_property = prop['PropId']
            prop_name = self.IdToPropType[id_property]['Name']
            prop_display_name = self.IdToPropType[id_property]['DisplayName']
            prop['PropName'] = prop_name
            prop['PropDisplayName'] = prop_display_name
            dict_folder_id = self.IdToPropType[id_property]['DictFolderId']
            if dict_folder_id > 0:
                dict_folder = self.get_dictionary(id_property, dict_folder_id)
                for i in range(0, len(prop['PropValues']['string'])):
                    id_dict_prop_value = int(prop['PropValues']['string'][i])
                    new_dict_value = dict_folder[id_dict_prop_value]
                    prop['PropValues']['string'][i] = new_dict_value

        if rp.ResultPos >= rp.ResultSize:
            self.SOAPclient.service.ResultRelease(rp.ResultRef)

        return obj_props,  obj_props.ResultSize, obj_props.ResultPos

    def put_experiment(self, dataframe: pd, expName, expType, entityType, map_entities_by, has_pvalue=True, description=''):
        # PSexp = self.SOAPclient.service('ns0:GetExperiment')
        row_count = dataframe.shape[0]
        col_count = dataframe.shape[1]
        min_fold_change = dataframe['log2FoldChange'].min()
        max_fold_change = dataframe['log2FoldChange'].max()
        sample_cnt = int(col_count - 1)
        if has_pvalue:
            sample_cnt = int(sample_cnt / 2)
        header = dataframe.columns

        z_experiment = self.SOAPclient.get_type('ns0:Experiment')
        z_sample = self.SOAPclient.get_type('ns0:SampleDefinition')

        sample = z_sample(
            Name=header[1],
            Type=1,
            Subtype=2,
            hasPValue=has_pvalue,
            Calculated=False,
            Attributes={'string': [header[1]]},
            Attributes_size=1
        )

        # ns0:Experiment( Owner: xsd:string, DateCreated: xsd:string, )
        ps_experiment = z_experiment(
            Id=0,
            Name=expName,
            Description=description,
            SampleCount=sample_cnt,
            RowsMappedCount=0,
            ExperimentType=1,  # ratio or logratio
            ExperimentSubType=1,  # logarithmic
            Organism='Homo sapiens',
            ExperimentTypeName=expType,
            EntityTypeName=entityType,
            RowsCount=row_count,
            GeneAttributeNames={'string': [header[0]]},  # {string: xsd:string[]}
            GeneAttributeNames_Size=1,
            OriginalGeneID_index=0,
            SampleDefinitions={'SampleDefinition': [sample]},  # {SampleDefinition: ns0:SampleDefinition[]},
            SampleDefinitions_Size=1,
            SampleAttributeNames={'string': ['phenotype']},
            SampleAttributeNames_Size=1,
            MinValueSignal=min_fold_change,
            MinValueRatio=min_fold_change,
            MaxValueSignal=max_fold_change,
            MaxValueRatio=max_fold_change,
            ReadOnly=True,
            MaskUsage=False
        )

        experiment_id = self.SOAPclient.service.PutExperiment(ps_experiment)

        row_attrs = list()
        sample_values = list()
        for i in range(0, row_count):
            gene_id = dataframe.iat[i, 0]
            row_attr = self.SOAPclient.get_type('ns0:ExperimentRowAttributes')
            rt = row_attr(
                EntityURN='?',
                EntityAttributes={'string': []},
                EntityAttributes_Size=0,
                OriginalGeneID='?',
                GeneAttributes={'string': [gene_id]},  # GeneAttributes: {string: xsd:string[]}
                GeneAttributes_Size=1,
                EntityId=0
            )
            row_attrs[i] = rt

            sample_val = self.SOAPclient.get_type('ns0:SampleValue')
            sv = sample_val(
                Value=dataframe.iat[i, 1],
                PValue=dataframe.iat[i, 2]
            )
            sample_values[i] = sv

        self.SOAPclient.service.PutExperimentRowAttributes(experiment_id, 0, 0, {'ExperimentRowAttributes': row_attrs})
        self.SOAPclient.service.PutExperimentValues(experiment_id, 0, 0, 0, {'SampleValue': sample_values})

        z_map_param = self.SOAPclient.get_type('ns0:ExperimentMappingToolParameters')
        zm = z_map_param(
            ExperimentID=experiment_id,
            DbType=entityType,
            DbField=map_entities_by,
            IsChip=False,
            ChipID=0,
            ChipName='?'
        )

        # ns0:OperationStatus(OperationId: xsd:long, State: xsd:int, Error: xsd:int, Status: xsd:string, ResultId: xsd:long)
        op_status = self.SOAPclient.service.StartExperimentMappingTool(zm)
        result = self.SOAPclient.service.GetMappingToolResult(op_status.OperationId)
        self.SOAPclient.service.ResultRelease(result.ResultRef)
        return experiment_id, result


def write_soap_response(fout, soap_response):
    with open(fout, "w", encoding='utf-8') as f:
        f.write(soap_response)

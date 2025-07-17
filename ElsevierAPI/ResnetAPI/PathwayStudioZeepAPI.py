import pandas as pd
import sys, time, logging, http.client, itertools
from time import sleep
from .PathwayStudioGOQL import OQL,len
from zeep import exceptions as zeep_exceptions
import requests.exceptions as req_exceptions
from zeep import helpers # converts zeep objects to dict
from ..utils import execution_time, load_api_config
from math import pow


CONNECTION_TIMEOUT = 20

def configure_logging(logger):
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger


class DataModel:
    def __init__(self, *args,**kwargs):
        '''
        APIconfig - args[0]\n
        no_mess - default True, if False your script becomes more verbose\n
        connect2server - default True, set to False to run script using data in __pscache__ files instead of database
        '''
        if args:
          self.APIconfig = load_api_config(args[0]) if isinstance(args[0], str) else dict(args[0])
        else:
          self.APIconfig = load_api_config()
        assert isinstance(self.APIconfig, dict), "APIconfig must be a dictionary or a path to a json file"
        my_kwargs = {'connect2server':True,'no_mess':True,'load_model':True}
        my_kwargs.update(kwargs)
        self.no_mess = my_kwargs.get('no_mess',True)
        
        self.IdToPropType = dict()
        self.IdtoObjectType = dict()
        self.propId2dict = dict()
        self.RNEFnameToPropType = dict()
           
        if my_kwargs['connect2server']:
            url = self.APIconfig[str('ResnetURL')]
            username = self.APIconfig['PSuserName']
            password = self.APIconfig['PSpassword']
        
            from requests.auth import HTTPBasicAuth
            from requests import Session
            from zeep.cache import SqliteCache
            from zeep.transports import Transport
            session = Session()
            #session.verify = False # not recommended. Uncomment if PS certificate has expired
            session.auth = HTTPBasicAuth(username, password)
            transport = Transport(cache=SqliteCache(), session=session)
            from zeep import Client, Settings
            settings = Settings(strict=False, xml_huge_tree=True)
            #settings = zeep.Settings(extra_http_headers={'Authorization': 'Bearer ' + token})
            self.logger = configure_logging(logging.getLogger(__name__))
            if my_kwargs['connect2server']:
                for attempt in range(1,11):
                    try:
                        self.SOAPclient = Client(wsdl=url, transport=transport, settings=settings)
                        if my_kwargs['load_model']:
                            self.__load_model()
                        if not self.no_mess:
                            print(f'New connection to Pathway Studio API server:\n{url} as {username}')
                        return
                    except Exception as error:
                        err_str = f"Pathway Studio server connection failed: {error}"
                        self.logger.error(err_str)
                        print(err_str)
                        attempt_timeout = CONNECTION_TIMEOUT*pow(attempt,2)
                        print(f'Pausing for {attempt_timeout} seconds on attmpt {attempt} due to "{error}"')
                        sleep(attempt_timeout)
                        continue
                raise req_exceptions.ConnectionError(f"Server connection failed. Wrong or inaccessible url: {url}") from None


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
          prop_type_name = str(property_types[i]['Name'])
          prop_type_display_name = str(property_types[i]['DisplayName'])
          self.IdToPropType[prop_type_name] = property_types[i]
          self.IdToPropType[prop_type_display_name] = property_types[i]
          self.IdToPropType[db_id] = property_types[i]
          
          
          self.RNEFnameToPropType[prop_type_name] = property_types[i]
      return
      


    def __id2objtype(self, classID):
        object_types = list()
        zeep_object_types = self.SOAPclient.service.GetObjectTypes()
        for i in range(0, len(zeep_object_types)):
            obj_type_name = zeep_object_types[i]['Name']
            obj_type_class_id = zeep_object_types[i]['ObjClassId']
            if obj_type_class_id == classID:
                object_types.append(obj_type_name)

        return object_types


    def get_relation_types(self):
        return self.__id2objtype(3)


    def get_entity_types(self):
        return self.__id2objtype(1)
    

    def db_triples(self,skip_nodetypes:list=[],skip_reltypes:list=[]):
        node_types = self.get_entity_types()
        relation_types = self.get_relation_types()
        node_types = [x for x in node_types if x not in skip_nodetypes]
        relation_types = [x for x in relation_types if x not in skip_reltypes]
        triples = list(set(itertools.product(node_types, relation_types,node_types)))
        return sorted(node_types), sorted(relation_types), sorted(triples)


    def load_folder_tree(self)->dict:
        """
        Returns
        -------
        {folder_id:folder}, where folder is zeep_object
        """
        print('Loading folders tree from database')
        folders = self.SOAPclient.service.GetFoldersTree(0)
        id2folders = dict()
        for folder in folders:
            id2folders[int(folder['Id'])] = folder
        
        return id2folders


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
        id_list = list()
        for x in PropertyNames:
            try:
                id_list.append(self.IdToPropType[x]['Id'])
            except KeyError:
                continue
        if 'Name' not in PropertyNames:
            id_list.append(self.IdToPropType['Name']['Id'])
        return id_list


    def __propid2name(self,idProperty:int,idDictFolder:int):
        try:
            return self.propId2dict[idProperty]
        except KeyError:
            for attempt in range(1,10):
                try:
                    dict_folder = self.SOAPclient.service.GetDictFolder(idDictFolder)
                    id_values_to_str = dict()
                    for val in dict_folder.Values.DictValue:
                        db_id = val['Id']
                        val = val['Value']
                        id_values_to_str[db_id] = val
                    self.propId2dict[idProperty] = id_values_to_str
                    self.propId2dict[dict_folder['Name']] = id_values_to_str
                    return self.propId2dict[idProperty]
                except Exception as error:
                    attempt_timeout = CONNECTION_TIMEOUT*pow(attempt,2)
                    print(f'Pausing for {attempt_timeout} seconds due to "{error}" on {attempt} attempt')
                    sleep(attempt_timeout)
                    continue
            raise req_exceptions.ConnectionError("Server connection failed after 10 attempts") 
             

    def get_folder_objects(self, FolderId, result_param):
        for i in range (0,3):
            try:
                result = self.SOAPclient.service.FolderGetObjects(FolderId, result_param)
                return result
            except zeep_exceptions.Fault: 
                continue


    def oql_response(self, OQLquery, result_param):
        timeout = 300 # 300 sec is default timeout in zeep
        max_iter = 10
        start = time.time()
        for attempt in range (1,max_iter):
            try:
                result = self.SOAPclient.service.OQLSearch(OQLquery, result_param)
                if attempt > 1:
                    print(f'{OQLquery[:100]} was succesful on {attempt} attempt')
                return result
            except (zeep_exceptions.Error,req_exceptions.RequestException) as err:
                print(f'\nSOAPclient session timed out after {execution_time(start)} on {attempt} \
out of {max_iter} attempt at GOQL query with length {len(OQLquery)}:\n{OQLquery[:100]}',flush=True)
                print(f'Timeout error message:\n{err}',flush=True)
                timeout = (attempt)*CONNECTION_TIMEOUT # None - default timeout in zeep
                print(f'\nWill make attempt #{attempt+1} with the same query after {timeout} seconds',flush=True)
                self.SOAPclient.transport.load_timeout = timeout
                sleep(timeout)
                continue
        
        tout = self.SOAPclient.transport.load_timeout
        if tout > 300:
            raise zeep_exceptions.TransportError(f'Timed out on GOQL query {OQLquery} after {tout} seconds',status_code=504)
        else:
            raise zeep_exceptions.Fault(f'Could not reconnect after 10 attempts while executing GOQL query "{OQLquery}"')


    def result_get_data(self, result_param):
        start = time.time()
        max_iter = 10
        for attempt in range (1,11):
            try:
                result = self.SOAPclient.service.ResultGetData(result_param)
                return result
            except (zeep_exceptions.Fault,req_exceptions.ChunkedEncodingError,http.client.RemoteDisconnected,req_exceptions.ConnectionError) as err:
                print(f'{err} error while retrieving results {result_param.ResultRef}')
                attempt_timeout = CONNECTION_TIMEOUT*attempt
                print(f'Attempt {attempt+1} to reconnect will be in {attempt_timeout} sec\n',flush=True)
                sleep(attempt_timeout)
                continue
            except zeep_exceptions.TransportError:
                print(f'\nSOAPclient session timed out after {execution_time(start)} on {attempt+1} attempt out of {max_iter}')
                timeout = (attempt+1)*300 # 300 - default timeout in zeep
                print(f'\nWill make attempt #{attempt+1} with the same query after {timeout} seconds',flush=True)
                self.SOAPclient.transport.load_timeout = timeout
                sleep(timeout)
                continue
        return


    def create_result_param(self, property_names=[]):
        property_names = property_names if property_names else ['Name'] 

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


    def __get_obj_props(self, obj_ids:list[int], property_names:list[str]=[]):
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
                dict_folder = self.__propid2name(id_property, dict_folder_id)
                for i in range(0, len(prop['PropValues']['string'])):
                    id_dict_prop_value = int(prop['PropValues']['string'][i])
                    new_dict_value = dict_folder[id_dict_prop_value]
                    prop['PropValues']['string'][i] = new_dict_value
        return obj_props
    

    def get_object_properties(self, obj_ids:list[int], property_names:list[str]=[]):
        max_chunk_size = 1000
        obj_props = self.__get_obj_props(obj_ids[:max_chunk_size],property_names)
        for chunk in range(max_chunk_size,len(obj_ids),max_chunk_size):
            chunk_props = self.__get_obj_props(obj_ids[chunk:chunk+max_chunk_size],property_names)
            obj_props['Objects']['ObjectRef'] += chunk_props.Objects.ObjectRef
            obj_props['Properties']['ObjectProperty'] += chunk_props.Properties.ObjectProperty
        return obj_props



    def get_folder_objects_props(self, FolderId:int, property_names=None):
        property_names = ['Name'] if property_names is None else property_names

        rp = self.create_result_param(property_names)
        rp.GetObjects = True
        rp.GetProperties = True
        obj_props = self.get_folder_objects(FolderId, rp)
        if obj_props is None: return None
        if obj_props.Objects is None: return None

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
            if prop['PropValues'] and dict_folder_id > 0:
                dict_folder = self.__propid2name(id_property, dict_folder_id)
                for i in range(0, len(prop['PropValues']['string'])):
                    id_dict_prop_value = prop['PropValues']['string'][i]
                    new_dict_value = dict_folder[id_dict_prop_value] if id_dict_prop_value in dict_folder else id_dict_prop_value
                    prop['PropValues']['string'][i] = new_dict_value

        return obj_props


    def get_layout(self, PathwayId,max_retries=10,retry_delay:float = 1.0):
        for attempt in range (max_retries):
            try:
                result = self.SOAPclient.service.GetObjectAttachment(PathwayId, 1)
                if 'Attachment' in result and result['Attachment'] is not None:
                    return str(result['Attachment'].decode('utf-8'))  
                else: 
                    return ''
            except (zeep_exceptions.Error,req_exceptions.RequestException):
                sleep(retry_delay)
                continue
        return ''
        

    def get_data(self, OQLrequest, retrieve_props=[], getLinks=True):

        retrieve_props = retrieve_props if retrieve_props  else ['Name', 'RelationNumberOfReferences']

        rp = self.create_result_param(retrieve_props)
        rp.GetObjects = True
        rp.GetProperties = True
        rp.GetLinks = getLinks
        # setting objectType name
        obj_props = self.oql_response(OQLrequest, rp)
        if obj_props is None:
            return None
        if obj_props.Objects is None:
            return None
        if not obj_props.Objects.ObjectRef:
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
                dict_folder = self.__propid2name(id_property, dict_folder_id)
                for i in range(0, len(prop['PropValues']['string'])):
                    id_dict_prop_value = int(prop['PropValues']['string'][i])
                    new_dict_value = dict_folder[id_dict_prop_value]
                    prop['PropValues']['string'][i] = new_dict_value
        return obj_props #helpers.serialize_object(result, target_cls=dict)


    def init_session(self, OQLrequest:str, PageSize:int, property_names=None, getLinks=True):
        rp = self.create_result_param(property_names)
        rp.GetObjects = True
        rp.GetProperties = True
        rp.GetLinks = getLinks
        rp.CreateResult = True
        rp.MaxPageSize = PageSize
        obj_props = self.oql_response(OQLrequest, rp)

        if type(obj_props) == type(None):
            return None, ('', 0, 0)
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
                dict_folder = self.__propid2name(id_property, dict_folder_id)
                for i in range(0, len(prop['PropValues']['string'])):
                    id_dict_prop_value = int(prop['PropValues']['string'][i])
                    new_dict_value = dict_folder[id_dict_prop_value]
                    prop['PropValues']['string'][i] = new_dict_value

        return obj_props, (obj_props.ResultRef, obj_props.ResultSize, obj_props.ResultPos)


    def get_session_page(self, ResultRef, ResultPos, PageSize, ResultSize, property_names=None, getLinks=True):
        '''
        Return
        ------
        tuple zeep_data, ResultSize, ResultPos
        '''
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

        if type(obj_props) == type(None):
            return None, int(0), int(0)
        if type(obj_props.Objects) == type(None):
            # print('Your SOAP response is empty! Check your OQL query and try again\n')
            return None,int(obj_props.ResultSize), int(obj_props.ResultPos)
        
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
                dict_folder = self.__propid2name(id_property, dict_folder_id)
                for i in range(0, len(prop['PropValues']['string'])):
                    id_dict_prop_value = int(prop['PropValues']['string'][i])
                    new_dict_value = dict_folder[id_dict_prop_value]
                    prop['PropValues']['string'][i] = new_dict_value

        if rp.ResultPos >= rp.ResultSize:
            self.SOAPclient.service.ResultRelease(rp.ResultRef)

        return obj_props, int(obj_props.ResultSize), int(obj_props.ResultPos)


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

    @staticmethod
    def write_soap_response(fout, soap_response):
        with open(fout, "w", encoding='utf-8') as f:
            f.write(soap_response)


    def create_GenRelNetsParameters(self, experiment_id:int, sample_id=0, 
                            expand_by_reltypes=[], expand2types=[], find_regulators=True):
            gnr_param = self.SOAPclient.get_type('ns0:GenRelNetsParameters')
            if find_regulators:
                fwd_reltypes = expand_by_reltypes
                fwd_targtypes = expand2types
                rvrs_reltypes = []
                rvrs_targtypes = []
            else:
                fwd_reltypes = []
                fwd_targtypes = []
                rvrs_reltypes = expand_by_reltypes
                rvrs_targtypes = expand2types

            gnr = gnr_param(
                    ExperimentID = experiment_id,
                    RatioSample = sample_id,
                    Mode = 1, # 1-SNEA; 3-BDEN; 4-GSEA
                    ForwardTypes = {'long':fwd_reltypes}, #for SNEA only
                    ForwardTypes_size = len(fwd_reltypes),
                    ReverseTypes = {'long':rvrs_reltypes}, #for SNEA only
                    ReverseTypes_size = len(rvrs_reltypes),
                    ForwardEntityTypes = {'long':fwd_targtypes}, #for SNEA only
                    ForwardEntityTypes_size = len(fwd_targtypes),
                    ReverseEntityTypes = {'long':rvrs_targtypes}, #for SNEA only
                    ReverseEntityTypes_size = len(rvrs_targtypes),
                    IsFindGroups = False, #obsolete
                    CutoffNumber = 100, #Cut-off by number of generated/found networks
                    MinPValue = 0.0,
                    MaxPValue = 0.5,
                    IsMetabolomicsMode = False,
                    IsIncludeOnlyMeasuredTargets = True, #for SNEA only
                    Folders = {'long':[0]}, #List of folders
                    Folders_size = 1,
                    Granularity = 0.1, # for BDEN
                    Targets = {'string':['']},   #AOQL queries specifying targets. for GSEA only.
                    Targets_size = 100,
                    Algorithm = 0, #for GSEA only: 0 - Mann-Whitney U-Test; 1 - Kolmogorov-Smirnov
                    Permute = 0, # for Kolmogorov-Smirnov GSEA only: 0 - permute genes, 1 - permute samples
                    NumberOfPermutations = 100, # for Kolmogorov-Smirnov GSEA only
                    EnrichmentStatistics = 0, # for Kolmogorov-Smirnov GSEA only: 0 - classic, 1 - weighted
                    ExpandContainers = False, # for GSEA only
                    UseMask = False,  #Use probe mask stored in the experiment
                    AnalyzePValues = False,  # calculate FDR?
                    UseConcordance  = True,
                    UnknownEffectBehavior = 0,
                    ActivationScoreUnknownEffectBehavior = 0, # for activation score in SNEA 
                    ChangeThreshold = 2.0 # for activation score in SNEA 
            )
            return gnr


    def __experiment_samples(self,experiment_id:int):
        for i in range (0,3):
            try:
                result = self.SOAPclient.service.GetExperiment(experiment_id)
                return result
            except zeep_exceptions.Fault: 
                continue


    def __get_experiment_identifiers(self,experiment_id:int, experiment_size:int):
        """
        returns dataframe of identifiers for Experiment
        """
        results = self.SOAPclient.service.GetExperimentRowAttributes(experiment_id, [], 0, experiment_size)
        urns = [a['EntityURN'] for a in results]
        ids = [a['OriginalGeneID'] for a in results]
        assert(len(urns) == len(ids))
        identifiers = pd.DataFrame.from_dict({'OriginalGeneID':ids,'URN':urns})
        return identifiers


    def _fetch_experiment_data(self, experiment_name:str, only_log_ratio_samples=True):
        """
        returns zeep objects for experiment, samples, urns, values
        """
        oql_query = 'SELECT Experiment WHERE Name = \'{}\''.format(experiment_name)
        experiment_zobj = self.get_data(oql_query)
        
        if type(experiment_zobj) == type(None):
            print('No experiment with name "%s" exist in database' % experiment_name)
            return None
        elif len(experiment_zobj.Objects.ObjectRef) > 1:
                print('!Experiment has duplicate names!  Only first experiment will be used')

        experiment_id = experiment_zobj.Objects.ObjectRef[0]['Id']
        experiment = self.__experiment_samples(experiment_id)
        experiment_identifiers = self.__get_experiment_identifiers(experiment_id, experiment.RowsCount)
        
        sample_values = list()
        sample_id = 0
        if only_log_ratio_samples:
            for sample in experiment.SampleDefinitions.SampleDefinition:
                if sample['Type'] == 1 and sample['Subtype'] == 2:
                    sample_values.append(self.SOAPclient.service.GetExperimentValues(experiment_id,sample_id, 0, len(experiment_identifiers)))
                sample_id += 1
        else:
            for sample in experiment.SampleDefinitions.SampleDefinition:
                sample_values.append(self.SOAPclient.service.GetExperimentValues(experiment_id,sample_id, 0, len(experiment_identifiers)))
                sample_id += 1
  
        return experiment, experiment_identifiers, sample_values


    def close_connection(self):
        self.SOAPclient.transport.session.close()

# -*- coding: utf-8 -*-
# Python wrapper for the Reaxys API
#
# Version: 1.1.0-beta.2
#
# Author:  Dr. Sebastian Radestock, Elsevier
# Author:  Dr. Alexander Riemer, Elsevier
# Author:  Dr. Markus Fischer, Elsevier
# Author:  Dr. Anton Yuryev, Elsevier
# Date:    July 26th, 2019
# Change Log 1.1.0-beta.1, July 26th, 2019
# A. Support for Python 3
# B. get_field_content modifications
# B.1. returns values for elements with highlights
# B.2. new method argument highlight_only. If True will return only a value if field contains highlights
#
# Change Log 1.1.0-beta.2, July 26th, 2019
# A. Method retrieve now supports clustering hitsets
# A.1. Added optional arguments dbname and context, that are required to formulate group by statements
#Change Log 1.2.0-beta.1, January 29th, 2021

import http.cookiejar, xml.dom.minidom, re
from urllib.request import Request, urlopen
from urllib.error import URLError
import http.client
from time import sleep
import json

DEFAULT_APICONFIG = 'D:/Python/ENTELLECT_API/ElsevierAPI/APIconfig.json'

def load_api_config(api_config_file=''):# file with your API keys and API URLs
    #cur_dir = os.getcwd()
    default_apiconfig = DEFAULT_APICONFIG
    if not api_config_file:
        print('No API config file was specified\nWill use default %s instead'% default_apiconfig)
        api_config_file = default_apiconfig
    try:
        return dict(json.load(open(api_config_file)))
    except FileNotFoundError:
        print("Cannot find API config file: %s" % api_config_file)
        if api_config_file != default_apiconfig:
            print('Cannot open %s config file\nWill use default %s instead'% (api_config_file, default_apiconfig))
            return dict(json.load(open(default_apiconfig)))
        else:
            print('No working API server was specified!!! Goodbye')
            return None


class Reaxys_API:

    def __init__(self, proxy=None, port=None):
        self.url = ""
        self.headers = {'Content-type' : 'text/xml; charset="UTF-8"'}
        self.callername = ""
        self.sessionid = ""
        self.resultname = ""
        self.resultsize = ""
        self.citationset = ""
        self.citationcount = ""
        self.proxy = proxy
        self.port = port
        #specify here the location of your API session file:
        self.ReaxysAPIjson = 'ElsevierAPI/ReaxysAPI/RxAPIsession.json'
        self.APIconfig = dict()

        # Set True for verbose output:
        self.debug = False

    def _get_resultname(self, response_xml):
        
        response_dom = xml.dom.minidom.parseString(response_xml)

        # Length of response_dom.getElementsByTagName("resultsname") should always be 1.
        # Node resultsname should not conatin subnodes.
        try:
            resultname = response_dom.getElementsByTagName("resultname")[0].childNodes[0].nodeValue
        except IndexError:
            resultname = None
        return resultname

    def _get_resultsize(self, response_xml):
        
        response_dom = xml.dom.minidom.parseString(response_xml)

        # Length of response_dom.getElementsByTagName("resultsize") should always be 1.
        # Node resultsize should not conatin subnodes.
        try:
            resultsize = response_dom.getElementsByTagName("resultsize")[0].childNodes[0].nodeValue
        except IndexError:
            resultsize = None

        return resultsize

    def _get_citationset(self, response_xml):
        
        response_dom = xml.dom.minidom.parseString(response_xml)

        # Length of response_dom.getElementsByTagName("citationset") should always be 1.
        # Node citationset should not conatin subnodes.          
        return response_dom.getElementsByTagName("citationset")[0].childNodes[0].nodeValue

    def _get_citationcount(self, response_xml):
        
        response_dom = xml.dom.minidom.parseString(response_xml)

        # Length of response_dom.getElementsByTagName("citationcount") should always be 1.
        # Node citationcount should not conatin subnodes.          
        return response_dom.getElementsByTagName("citationcount")[0].childNodes[0].nodeValue

    def get_facts_availability(self, response_xml, field):

        facts_availability = "0"
        
        response_dom = xml.dom.minidom.parseString(response_xml)

        facts = response_dom.getElementsByTagName("facts")[0]
        for fact in facts.childNodes:
            if 'name="' + field + '"' in fact.toxml():
                facts_availability = fact.childNodes[0].nodeValue.split("(")[0]

        return facts_availability

    def get_field_content(self, response_xml, field, highlight_only=False):
        
        field_content = []
        
        response_dom = xml.dom.minidom.parseString(response_xml)
        
        for element in response_dom.getElementsByTagName(field):

            # Concatenate text values if highlight is present
            if element.getAttribute('highlight') == 'true':
                field_content.append(
                    ''.join([e.data
                             if type(e) == xml.dom.minidom.Text
                             else e.childNodes[0].data for e in element.childNodes]))

            # If node contains further sub-nodes: return full xml.
            elif len(element.childNodes) > 1 and highlight_only is False:
                field_content.append(element.toxml())

            # If node does not conatin further sub-nodes: return node value.
            elif len(element.childNodes) == 1 and highlight_only is False:
                field_content.append(element.childNodes[0].nodeValue)
                
        return field_content

    def connect(self, url, username, password, callername):       
        self.url = url
        self.callername = callername
        cookies = http.cookiejar.CookieJar()
        
        connect_template = """<?xml version="1.0"?>
          <xf>
            <request caller="%s">
              <statement command="connect" username="%s" password="%s"/>
            </request>
          </xf>\n"""
        payload = connect_template % (callername, username, password)
        data = payload.encode()

        # Header reset.
        self.headers = {'Content-type' : 'text/xml; charset="UTF-8"'}

        # ELSAPI support
        self.headers['X-ELS-APIKey'] = callername
        self.headers['Accept'] = "*/*"
        request = Request(self.url, data=data, headers=self.headers)
        
        if self.debug:
            print('-----------------------\nQuery headers from connect:')
            print(self.headers)
            print('-----------------------\nQuery from connect:')
            print(payload)

        response = urlopen(request)
        response_xml = response.read()

        if self.debug:
            print('-----------------------\nResponse headers from connect:')
            print(response.info())
            print('-----------------------\nResponse from connect:')
            print(response_xml)

        # Get sessionid.
        response_dom = xml.dom.minidom.parseString(response_xml)
        element = response_dom.getElementsByTagName("sessionid")
        if len(element) > 0:
            self.sessionid = element[0].childNodes[0].nodeValue
        
        # Cookies are read from the response and stored in self.header
        #     which is used as a request header for subsequent requests.
        cookies.extract_cookies(response, request)
        
        # Cookie handling 3.0: Simply store and resend ALL cookies received from server
        self.headers['Cookie'] = "; ".join(re.findall(r"(?<=Cookie ).*?=\S*", str(cookies)))

    def disconnect(self):
        disconnect_template = """<?xml version="1.0"?>
          <xf>
            <request caller="%s">
              <statement command="disconnect" sessionid="%s"/>
            </request>
          </xf>\n"""
        payload = disconnect_template%(self.callername, self.sessionid)
        data = payload.encode()

        request = Request(self.url, data=data, headers=self.headers)

        if self.debug:
            print('-----------------------\nQuery headers from disconnect:')
            print(self.headers)
            print('-----------------------\nQuery from disconnect:')
            print(payload)

        response = urlopen(request)
        response_xml = response.read()
        
        if self.debug:
            print('-----------------------\nResponse headers from disconnect:')
            print(response.info())
            print('-----------------------\nResponse from disconnect:')
            print(response_xml)

    def select(self, dbname, context, where_clause, order_by, options):
        
        select_template = """<?xml version="1.0" encoding="UTF-8"?>
          <xf>
            <request caller="%s" sessionid="">
              <statement command="select"/>
              <select_list>
                <select_item/>
              </select_list>
              <from_clause dbname="%s" context="%s">
              </from_clause>
              <where_clause>%s</where_clause>
              <order_by_clause>%s</order_by_clause>
              <options>%s</options>
            </request>
          </xf>\n"""
        payload = select_template%(self.callername, dbname, context, where_clause, order_by, options)
        data = payload.encode()
        request = Request(self.url, data=data, headers=self.headers)

        if self.debug:
            print('-----------------------\nQuery headers from select:')
            print(self.headers)
            print('-----------------------\nQuery from select:')
            print(payload)

        try:
            response = urlopen(request)
            response_xml = response.read()
        except URLError:
            sleep(60)
            self.OpenSession()
            response = urlopen(request)
            response_xml = response.read()
        except http.client.RemoteDisconnected:
            sleep(60)
            self.OpenSession()
            response = urlopen(request)
            response_xml = response.read()
        
        if self.debug:
            print('-----------------------\nResponse headers from select:')
            print(response.info())
            print('-----------------------\nResponse from select:')
            print(response_xml)

        self.resultname = self._get_resultname(response_xml)
        self.resultsize = self._get_resultsize(response_xml)
        
        if ("NO_CORESULT" not in options) and ("C" not in context):
            self.citationset = self._get_citationset(response_xml)
            self.citationcount = self._get_citationcount(response_xml)

        return response_xml

    def expand(self, dbname, first_item, last_item, where_clause):
        
        select_template = """<?xml version="1.0" encoding="UTF-8"?>
          <xf>
            <request caller="%s" sessionid="%s">
              <statement command="expand"/>
              <from_clause dbname="%s" first_item="%s" last_item="%s">
              </from_clause>
              <where_clause>%s</where_clause>
            </request>
          </xf>\n"""
        payload = select_template%(self.callername, self.sessionid, dbname, first_item, last_item, where_clause)
        data = payload.encode()
        request = Request(self.url, data=data, headers=self.headers)

        if self.debug:
            print('-----------------------\nQuery headers from expand:')
            print(self.headers)
            print('-----------------------\nQuery from expand:')
            print(payload)
        
        response = urlopen(request)
        response_xml = response.read()
        
        if self.debug:
            print('-----------------------\nResponse headers from expand:')
            print(response.info())
            print('-----------------------\nResponse from expand:')
            print(response_xml)

        return response_xml

    def post(self, payload):

        data = payload.encode()
        request = Request(self.url, data=data, headers=self.headers)

        if self.debug:
            print('-----------------------\nQuery headers from post:')
            print(self.headers)
            print('-----------------------\nQuery from post:')
            print(payload)
        
        response = urlopen(request)
        response_xml = response.read()
        
        if self.debug:
            print('-----------------------\nResponse headers from post:')
            print(response.info())
            print('-----------------------\nResponse from post:')
            print(response_xml)

    def retrieve(self, resultname, select_items, first_item, last_item, order_by, group_by, group_item, options,
                 dbname=None, context=None):
        # if group_by is given, please provide group_item, e.g. "1" or "1,2"
        
        if group_by != '':
            grouplist = 'grouplist="' + group_item + '"'
        else:
            grouplist = ''

        db_template = ''
        if dbname is not None:
            db_template = 'dbname="%s"' % dbname

        context_template = ''
        if context is not None:
            context_template = 'context="%s"' % context


        select_item_template = """                <select_item>%s</select_item>\n"""
        select_template = """<?xml version="1.0" encoding="UTF-8"?>
          <xf>
            <request caller="%s" sessionid="%s">
              <statement command="select"/>
              <select_list>\n"""
        for index in range (0,len(select_items)):
            select_template = select_template + select_item_template%(select_items[index])
        select_template = select_template + """              </select_list>
              <from_clause %s %s resultname="%s" %s first_item="%s" last_item="%s">
              </from_clause>
              <order_by_clause>%s</order_by_clause>
              <group_by_clause>%s</group_by_clause>
              <options>%s</options>
            </request>
          </xf>\n"""
        payload = select_template % (
            self.callername, self.sessionid, db_template, context_template, resultname, grouplist,
            first_item, last_item, order_by, group_by, options)
        data = payload.encode()
        
        request = Request(self.url, data=data, headers=self.headers)

        if self.debug:
            print('-----------------------\nQuery headers from retrieve:')
            print(self.headers)
            print('-----------------------\nQuery from retrieve:')
            print(payload)
        
        response = urlopen(request)
        response_xml = response.read().decode()
        
        if self.debug:
            print('-----------------------\nResponse headers from retrieve:')
            print(response.info())
            print('-----------------------\nResponse from retrieve:')
            print(response_xml)

        return response_xml

    def __session_exists(self):
        hydrogenRXN = 3587189
        self.select(dbname="RX", context="S", where_clause="IDE.XRN = '"+str(hydrogenRXN)+"'",order_by="",options="WORKER,NO_CORESULT")
        if type(self.resultsize) == type(None):
            return False
        else: 
            return True

    def __read_from_json(self):
        try: 
            with open(self.ReaxysAPIjson) as f:
                data = json.load(f)
                self.url = data['url']
                #self.headers = {'Content-type' : 'text/xml; charset="UTF-8"'}
                self.headers = data['headers']
                self.callername = data['callername']
                self.sessionid = data['sessionid']
                self.resultname = data['resultname']
                self.resultsize = data['resultsize']
                self.citationset = data['citationset']
                self.citationcount = data['citationcount']
                self.proxy = data['proxy']
                self.port = data['port']
            #print('Read %s session file' % (self.ReaxysAPIjson))
        except FileNotFoundError:
            print('Cannot find %s session file' % self.ReaxysAPIjson)
            return FileNotFoundError

    def _start_new_session(self, url, username, password, callername):
        self.connect(url,username,password,callername)
        if self.__session_exists() == False:
            print('Cannot create new session. Too many current active logins')
            self.debug = True
            print('Attempting to disconnect session %s' % self.sessionid)
            self.disconnect()
            self.debug = False
            print('Please wait for 30 minutes to expire your session on Reaxys server')
            return
        cache = open(self.ReaxysAPIjson, "w")
        apiJson = json.dumps(self.__dict__,indent=1, sort_keys=True)+'\n'
        cache.write(apiJson)
        cache.close()
        print("New %s session file was created" % self.ReaxysAPIjson)


    def OpenSession(self,APIconfig:dict=dict()):
        if not APIconfig:
            APIconfig = load_api_config()

        if self.__read_from_json() == FileNotFoundError:
            self._start_new_session(APIconfig['ReaxysURL'],APIconfig['ReaxysUsername'],APIconfig['ReaxysPassowrd'],APIconfig['ReaxysAPIkey'])
    
        if self.__session_exists() == False:
            self._start_new_session(APIconfig['ReaxysURL'],APIconfig['ReaxysUsername'],APIconfig['ReaxysPassowrd'],APIconfig['ReaxysAPIkey'])
        

    def GetCompoundProps(self,CompoundIDs:list,CompIdType:str,Fields:list,retrive_options=''):
        FieldToValues = dict()
        #move next line to the function that calls this one for speed
        #self.OpenSession()
        ctx = 'S'
        for Id in CompoundIDs:
            self.select(dbname="RX", context=ctx, where_clause=CompIdType+" = '"+str(Id)+"'",order_by="",options="WORKER,NO_CORESULT")
            if type(self.resultsize) != type(None):
                selectItems = {f.split('.')[0] for f in Fields}
                response = self.retrieve(self.resultname, list(selectItems), str(1),str(self.resultsize), "", "", "", retrive_options)
                for field in Fields:
                    values = set(self.get_field_content(response, field))
                    try:
                        v = set(FieldToValues[field])
                        v.update(values)
                        FieldToValues[field] = list(v)
                    except KeyError:
                        if len(values) > 0:
                            FieldToValues[field] = list(values)

        #self.disconnect() #move this line to the function that calls this one for speed
        return FieldToValues
    

    def GetTargetProps(self,target_names:list,Fields:list,retrive_options=''):
        FieldToValues = dict()
        #move next line to the function that calls this one for speed
        #self.OpenSession()
        ctx = 'DPI'
        query = "DAT.TNAME='{}'"
        for name in target_names:
            q_t = query.format(name)
            self.select(dbname="RX", context=ctx, where_clause=q_t,order_by="",options="WORKER,NO_CORESULT")
            if type(self.resultsize) != type(None):
                selectItems = {f.split('.')[0] for f in Fields}
                response = self.retrieve(self.resultname,list(selectItems),str(1),str(self.resultsize),"", "", "", "")
                for field in Fields:
                    values = set(self.get_field_content(response, field))
                    try:
                        v = set(FieldToValues[field])
                        v.update(values)
                        FieldToValues[field] = list(v)
                    except KeyError:
                        if len(values) > 0:
                            FieldToValues[field] = list(values)

        #self.disconnect() #move this line to the function that calls this one for speed
        return FieldToValues
    


    def getpX(self,drug_name:str,target_name:str):
        '''
        Return
        ------
        max pX for drug-target pair
        '''

        #move next line to the function that calls this one for speed
        #self.OpenSession()
        ctx = 'DPI'
        #fields = ['DAT.PAUREUS']
        q = f"IDE.CN={drug_name} and DAT.TNAME='{target_name}' and DAT.CATEG = 'in vitro (efficacy)'"

        self.select(dbname="RX", context=ctx, where_clause=q,order_by="DAT.PAUREUS desc",options="WORKER,NO_CORESULT")
        if type(self.resultsize) != type(None):
            #selectItems = {f.split('.')[0] for f in fields}
            response = self.retrieve(self.resultname, ["DAT"], "1", str(self.resultsize), "", "", "", "")
            values = set(self.get_field_content(response, 'DAT.PAUREUS'))
            if values:
                max_pX = max(list(map(float,values)))
                return max_pX
            else:
                return 0

    


        

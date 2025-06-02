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

import http.cookiejar, xml.dom.minidom, re, os, json, http.client
from lxml import etree as et 
from urllib.request import Request, urlopen
from urllib.error import URLError
from time import sleep
from ..utils import load_api_config, most_frequent, multithread
from ..ETM_API.references import Reference,PUBYEAR,AUTHORS,PATENT_APP_NUM,JOURNAL,SENTENCE,MEASUREMENT
from ..ResnetAPI.NetworkxObjects import EFFECT,PSObject
from collections import defaultdict



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
      self.ReaxysAPIjson = os.path.join(os.getcwd(),'ENTELLECT_API/ElsevierAPI/ReaxysAPI/RxAPIsession.json')
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
        field_content.append( ''.join([e.data
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
      where_clause = f"IDE.XRN = '{hydrogenRXN}'"
      self.select(dbname="RX", context="S", where_clause=where_clause,order_by="",options="WORKER,NO_CORESULT")
      #response = self.retrieve(self.resultname, ['IDE.XRN'], str(1),str(self.resultsize), "", "", "", '')
      return bool(self.resultsize)


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
      if not self.__session_exists():
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


  def OpenSession(self,APIconfig=dict()):
      if not APIconfig:
        APIconfig = load_api_config()

      if self.__read_from_json() == FileNotFoundError:
        self._start_new_session(APIconfig['ReaxysURL'],APIconfig['ReaxysUsername'],APIconfig['ReaxysPassowrd'],APIconfig['ReaxysAPIkey'])
  
      if not self.__session_exists(): # case when old session file is expired
        self._start_new_session(APIconfig['ReaxysURL'],APIconfig['ReaxysUsername'],APIconfig['ReaxysPassowrd'],APIconfig['ReaxysAPIkey'])
  

  def GetCompoundProps(self,CompoundIDs:list[str],CompIdType='IDE.XRN',Fields=['IDE.MF'],
                        retrive_options='')->dict[str,PSObject]:
    '''
    output:
      {Id:PSObject}, where PSObject has values form Fields
    '''
    compounds = dict()
    # for speed open session - self.OpenSession() - prior to calling "GetCompoundProps" 
    # to avoid re-opening Reaxys session multiple sessions
    for Id in CompoundIDs:
      rx_compound = PSObject({CompIdType:[Id]})
      wh_clause = f"{CompIdType} = '{str(Id)}'"
      self.select(dbname="RX", context='S', where_clause=wh_clause,order_by="",options="WORKER,NO_CORESULT")
      if self.resultsize:
        selectItems = {f.split('.')[0] for f in Fields}
        response = self.retrieve(self.resultname, list(selectItems), str(1),str(self.resultsize), "", "", "", retrive_options)
        for field in Fields:
          values = set(self.get_field_content(response, field))
          rx_compound.update({field:list(values)})
        compounds[Id]=rx_compound
    # for speed disconnect Reaxys session - self.disconnect() - after calling "GetCompoundProps"
    return compounds
  

  def chemid2props(self,id:str,fields:list[str],id_code = 'IDE.XRN')->dict[str,list]:
    self.select(dbname="RX", context='S', where_clause=f"{id_code} = '{str(id)}'",order_by="",options="WORKER,NO_CORESULT")
    if self.resultsize:
      selectItems = {f.split('.')[0] for f in fields}
      response = self.retrieve(self.resultname, list(selectItems), str(1),str(self.resultsize), "", "", "", "")
      return {f:self.get_field_content(response, f) for f in fields}
    else:
      return dict()
  

  def drug2props(self,drug:PSObject,fields:dict[str,str])->dict[str,dict]:
    '''
    input:
      drugs: list of PSObject with Reaxys ID
      fields: {Reaxys field name : human readable field name}
    output:
      {human readable field name : most_frequent value}
    '''
   # if drug.name() in ['vactosertib','breviscapine']:
   #     print()
    f2props = defaultdict(list)
    for rxid in drug.get_props(['Reaxys ID']):
      f2prop = self.chemid2props(rxid,fields.keys())
      [f2props[f].extend(v) for f,v in f2prop.items()]
    
    if not f2props:
      f2props = self.chemid2props(drug.name(),fields.keys(),id_code='IDE.CN')
    
    if f2props:
      most_common_values = {fields[f]:most_frequent(v) for f,v in f2props.items()}
      return most_common_values
    else:
        print(f'No properties found in Reaxys for {drug.name()}')
        return dict()
  

  def GetTargetProps(self,target_names:list,Fields:list):
      FieldToValues = dict()
      #move next line to the function that calls this one for speed
      #self.OpenSession()
      for name in target_names:
          q_t =  f"DAT.TNAME='{name}'"
          self.select(dbname="RX", context='DPI', where_clause=q_t,order_by="",options="WORKER,NO_CORESULT")
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
      q = f"IDE.CN='{drug_name}' and DAT.TNAME='{target_name}' and DAT.CATEG = 'in vitro (efficacy)'"
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
  
  
  @staticmethod
  def _parse_citations(data:et._Element,props:dict=dict())->list[Reference]:
      refs = list()
      for citation in data.findall('./citations/citation/CIT'):
          pui = citation.find('./CIT.PUI').text
          ref = Reference('PUI',pui)
          ref[PUBYEAR] = [citation.find('./CIT.PREPY').text]
          cittype = citation.find('./CIT.DT').text
          if cittype == 'Patent':
              ref[AUTHORS] = [citation.find('./CIT.OPA').text]
              #country = citation.find('./CIT02/CIT.PCC').text
              pn = str(citation.find('./CIT02/CIT.PN').text)
              pn = pn.replace('/','')
              version = citation.find('./CIT02/CIT.PK').text
              patnum = pn+version
              ref.Identifiers[PATENT_APP_NUM] = patnum
          elif cittype in ['Article','Note']:
              ref[AUTHORS] = [citation.find('./CIT.AU').text]
              ref[JOURNAL] = [citation.find('./CIT01/CIT.JTS').text]
              ref['ISSN'] = [citation.find('./CIT01/CIT.ISSN').text]
              ref.Identifiers['DOI'] = citation.find('./CIT01/CIT.DOI').text
          else:
              print('Unknown reference type')
          
          textref = ref._make_textref()+'#curation'
          ref.add_sentence_prop(textref,SENTENCE,'Imported from Reaxys')
          [ref.add_sentence_prop(textref,propname,value) for propname,value in props.items()]
          refs.append(ref)
      return refs


  @staticmethod
  def __effect_str(ACTTRG:str):
      low = ACTTRG.lower()
      if low in {'agonist','activator','partial agonist'}: return 'positive'
      if low.startswith('stimulator'): return 'positive'
      if low.startswith(('radioligand','allosteric')): return 'unknown'
      return 'negative'


  def getTargets(self,drug_name:str)->list[tuple[list,dict,dict]]:
      '''
      Return
      ------
      [[refs],rel,target]
      '''
      #move next line to the function that calls this one for speed
      #self.OpenSession()
      ALLOWED_ORGANISM = {'human','mouse','rat','homo sapience','mus musculus','rattus norvegicus'}
      def organism_name(text:str):
          vocab = {'human':'Homo sapiens','mouse':'Mus musculus','rat':'Rattus norvegicus'}
          return vocab.get(text,'')

      q = f"IDE.CN='{drug_name}' and DAT.TNAME='*' and DAT.CATEG = 'in vitro (efficacy)'"
      rel_target = list()
      self.select(dbname="RX", context='DPI', where_clause=q,order_by="DAT.PAUREUS desc",options="WORKER,NO_CORESULT")
      if type(self.resultsize) != type(None):
          r = self.retrieve(self.resultname, ["DAT","DATIDS"], "1", str(self.resultsize), "", "", "", "")
          response = et.fromstring(r.encode())
          #response = json.loads(response)
          for dpitem in response.findall('./dpitems/dpitem'):
              data = dpitem.find('./DAT')
              _organism = data.find('./DAT02/DAT.TSPECIE')
              if _organism is None or _organism.text.lower() not in ALLOWED_ORGANISM:
                  continue
              
              _eff = data.find('./DAT.ACTTRG')
              effect = 'unknown' if _eff is None else self.__effect_str(_eff.text)
              rel = defaultdict(list)
              rel[EFFECT].append(effect)

              ref_props = dict()
              _px = data.find('./DAT.PAUREUS')
              if _px is not None:
                  ref_props['pX'] = _px.text

              _cell = data.find('./DAT05/DAT.BCELL')
              if _cell is not None:
                  ref_props['CellLineName'] = _cell.text

              _organ = data.find('./DAT05/DAT.BTISSUE')
              if _organ is not None:
                  ref_props['Organ'] = _organ.text

              measurement = ''
              _measurment = data.find('./DAT.VTYPE')
              if _measurment is not None:
                  _value = data.find('./DAT.VALUE')
                  if _value is not None:
                      _unit = data.find('./DAT.UNIT')
                      if _unit is not None:
                          measurement = _measurment.text+'='+_value.text+_unit.text
              if measurement:
                  ref_props[MEASUREMENT] = measurement

              ref_props['Organism'] = organism_name(_organism.text.lower()).capitalize()

              refs = self._parse_citations(data,ref_props)

              dataids = dpitem.find('./DATIDS')
              target = defaultdict(list)
              target['Name'].append(data.find('./DAT02/DAT.TSUBUNIT').text)
              l =  [e.text for e in dataids.findall('./DATIDS.TSUBSYN')]
              if l:
                  target['Alias'] += l
              g = [e.text for e in dataids.findall('./DATIDS.TUNIPROT')]
              if g:
                  target['GenBank ID'] += g
              rel_target.append((refs,dict(rel),dict(target)))

      return rel_target


def drugs2props(drugs:list[PSObject],fields:dict[str,str])->dict[str,dict]:
  '''
  input:
    drugs: list of PSObject with Reaxys ID
    fields: {Reaxys field name : human readable field name}
  output:
    {field:drug_name:value}
  '''
  rxapi = Reaxys_API()
  rxapi.OpenSession()
  f2d2props = {f:dict() for f in fields.values()}
  drug_props = multithread(drugs,rxapi.drug2props,fields=fields,max_workers=2)
  assert(len(drug_props) == len(drugs)), f'Error: {len(drug_props)} != {len(drugs)}'
  for i, drug in enumerate(drugs):
    for f,v in drug_props[i].items():
      if v:
        f2d2props[f][drug.name()] = v
  rxapi.disconnect()
  return f2d2props
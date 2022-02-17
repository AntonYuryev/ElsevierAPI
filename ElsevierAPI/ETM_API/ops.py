
import os
import xml.etree.ElementTree as et
from xml.dom import minidom
import epo_ops
import time
import json
import requests
import urllib.parse
from .references import DocMine,TITLE,AUTHORS,ABSTRACT,PATENT_APP_NUM,PATENT_GRANT_NUM,CLAIMS

from requests.models import HTTPError

inhaling_preparations = {'A61K009/12', 'A61K009/0073', 'A61K009/0075','A61K009/0078','A61K009/008','A61K49/1815','A61K8/046','A61K51/1206','A61K51/1231'}
inhaling_mixing_methods = {'B01F003/04992','B01F003/04007','B01F2003/0057'} #CPC/B01F3/04992 - to search in uspto
inhaling_devices= {'A61M015/00','A61M015/0001','A61M015/0003','A61M015/0005','A61M015/0006','A61M015/0008','A61M015/001','A61M015/0011','A61M015/0013','A61M015/0015','A61M015/0016','A61M015/0018','A61M015/002','A61M015/0021','A61M015/0023','A61M015/0025','A61M015/0026','A61M015/0028','A61M015/003','A61M015/0031','A61M015/0033','A61M015/0035','A61M015/0036','A61M015/0038','A61M015/004','A61M015/0041','A61M015/0043','A61M015/0045','A61M015/0046','A61M015/0048','A61M015/005','A61M015/0051','A61M015/0053','A61M015/0055','A61M015/0056','A61M015/0058','A61M015/006','A61M015/0061','A61M015/0063','A61M015/0065','A61M015/0066','A61M015/0068','A61M015/007','A61M015/0071','A61M015/0073','A61M015/0075','A61M015/0076','A61M015/0078','A61M015/008','A61M015/0081','A61M015/0083','A61M015/0085','A61M015/0086','A61M015/0088','A61M015/009','A61M015/0091','A61M015/0093','A61M015/0095','A61M015/0096','A61M015/0098','A61M015/02','A61M015/025','A61M015/06','A61M015/08','A61M015/085'}
smoking_devices= {'A24F040', 'A24F042', 'A24F047'}

allowed_ipc_codes = inhaling_preparations|inhaling_mixing_methods|inhaling_devices
allowed_ipc_codes = inhaling_mixing_methods

inhaling_preparations = ['A61K9/0073', 'A61K9/0075','A61K9/0078','A61K8/046','A61K51/1206','A61K51/1231','A61K49/1815']#,'A61K9/008','A61K9/12']

OPS_APIKEY = 'zMxluGSzvCvEtdQ1xP1QK7zA9Gt6fpZc'
OPS_SECRET = 'QyM48SYC04K7QSyh'

 
def cleanup_authors(names:list):
    to_return = set()
    for name in names:
        n = str(name)
        n = n.strip(' .,\n\t')
        if n[-1] == ']':
            n = n[:n.rfind('[')]
        if n[-1] == '\u2002':
            n = n[:-1]
        
        n = n.title().replace(',','')
        to_return.add(n)
    return to_return

def cleanup_institutions(names:list):
    to_return = list(cleanup_authors(names))
    for i in range(0, len(to_return)):
        n = str(to_return[i])
        n = n.replace(' Corporation', '')
        n = n.replace(' Corp', '')
        n = n.replace(' Incorporated', '')
        n = n.replace(' Limited', '')
        n = n.replace(' S.P.A', '')
        n = n.replace(' Gmbh', '')
        if n[:4] == 'The': 
            n = n[4:]

        if n[-3:] in {' Co',' Sa',' Bv',' Cv',' Ag'}:
            n = n[:-3]
        if n[-4:] in {' Spa',' Ltd',' Llc',' Inc',' B V',' B.V',' S.A',' C.V'}:
            n = n[:-4]
        to_return[i] = n
    return set(to_return)

def print_tags(tree:et.ElementTree):
    print ([e.tag for e in tree.iter()]) #to get list of tags


class OPSxml(DocMine):
    @staticmethod
    def tag(tag):
        return '{http://www.epo.org/exchange}'+tag


    def __init__(self, patent_in_opsxml:et.Element):
        self.document = patent_in_opsxml.find(self.tag('exchange-documents')+'/'+self.tag('exchange-document'))

        doc_id = self.document.find(self.tag('bibliographic-data')+'/'+self.tag('publication-reference')+'/'+self.tag('document-id'))
        pat_num = self.document.get('country') + self.document.get('doc-number')+ self.document.get('kind')
        super().__init__(PATENT_GRANT_NUM, pat_num)

        app_ref = self.document.find(self.tag('bibliographic-data')+'/'+self.tag('application-reference'))
        
        applications = list()
        doc = app_ref.find(self.tag('document-id'))
            #print_tags(doc)
        f = doc.find(self.tag('country'))
        app_country = f.text
        app_num = doc.find(self.tag('doc-number')).text
        app_kind = doc.find(self.tag('kind')).text
        self.Identifiers[PATENT_APP_NUM] = app_country+app_num+app_kind

        #self.patent = patent_in_opsxml
        self['file'] = ['']
        self.append_property('doc_number',self.document.get('doc-number'))
        self.append_property('country', self.document.get('country'))
        self.append_property('kind', self.document.get('kind'))
        self.append_property('family-id',self.document.get('family-id'))
        date = doc_id.find(self.tag('date')).text
        self.set_date(date[:4],date[4:6],date[6:8])

    @staticmethod
    def loads(xml_file:str):
        try:
            tree = et.parse(xml_file)
        except et.ParseError: 
            return None
        return tree.getroot()
    
    @classmethod
    def from_file(cls, xml_file:str):
        try:
            patent = cls.loads(xml_file)
            if isinstance(patent,et.Element):
                c= cls(patent)
                file = xml_file[xml_file.rfind('/')+1:]
                c['file'] = [file[file.rfind('\\')+1:]]
                return c
        except et.ParseError: 
            return None

    @staticmethod
    def parse_patnum(patnum:str):
        allowed_chars = ['A','B','C','U','T','L','E']
        country = patnum[:2]

        if patnum[-2] in allowed_chars:
            kind_pos = len(patnum)-2
        elif patnum[-1] in allowed_chars:
            kind_pos = len(patnum)-1
        else: kind_pos = len(patnum)

        doc_number = patnum[2:kind_pos]
        kind = patnum[kind_pos:]

        return country, doc_number, kind


    def get_cpc_codes(self): ##get document from patent using opsxml_get_docnum
        self.cpc_codes = set()
        patent_classifications = self.document.findall('./{ns1}bibliographic-data/{ns1}patent-classifications/{ns1}patent-classification'.format(ns1=self.ns1))
        for classification in patent_classifications:
            section = classification.find(self.tag('section')).text
            code_class = classification.find(self.tag('class')).text
            subclass = classification.find(self.tag('subclass')).text
            main_group = classification.find(self.tag('main-group')).text
            subgroup = classification.find(self.tag('subgroup')).text
            self.update_with_value('cpc_codes',section+code_class+subclass+main_group+'-'+subgroup)

    def has_cpc_codes(self, allowed_cpc_codes:set): ##get document from patent using opsxml_get_docnum
        patent_classifications = self.document.findall('./{ns1}bibliographic-data/{ns1}patent-classifications/{ns1}patent-classification'.format(ns1='{http://www.epo.org/exchange}'))
        for classification in patent_classifications:
            section = classification.find(self.tag('section')).text
            code_class = classification.find(self.tag('class')).text
            subclass = classification.find(self.tag('subclass')).text
            main_group = classification.find(self.tag('main-group')).text
            subgroup = classification.find(self.tag('subgroup')).text
            cpc_code = section+code_class+subclass+main_group+'-'+subgroup
            if cpc_code in allowed_cpc_codes: return cpc_code
        return ''

    def set_title(self):
        if not hasattr(self,'title'):
            tit = self.document.find('./{ns1}bibliographic-data/{ns1}invention-title[@lang=\'en\']'.format(ns1='{http://www.epo.org/exchange}'))
            title = tit.text.strip().title() if type(tit) != type(None) else ''
            self._set_title(title)

    def get_data(self):
        applicants = self.document.findall('./{ns1}bibliographic-data/{ns1}parties/{ns1}applicants/{ns1}applicant/{ns1}applicant-name/{ns1}name'.format(ns1='{http://www.epo.org/exchange}'))
        applicant_names = {x.text for x in applicants}
        applicant_names = cleanup_institutions(list(applicant_names))
        inventors = self.document.findall('./{ns1}bibliographic-data/{ns1}parties/{ns1}inventors/{ns1}inventor/{ns1}inventor-name/{ns1}name'.format(ns1='{http://www.epo.org/exchange}'))
        inventor_names = {x.text for x in inventors}
        self.update_with_list(AUTHORS,cleanup_authors(list(inventor_names)))

        #self.update_with_list(INSTITUTIONS,list(applicant_names.difference(inventor_names)))
        self.addresses = {i:i for i in list(applicant_names.difference(inventor_names))}

        self.set_title()
        abs = self.document.find('./{ns1}abstract[@lang=\'en\']/{ns1}p'.format(ns1='{http://www.epo.org/exchange}'))
        abstract = abs.text if type(abs) != type(None) else ''
        self.add2section(ABSTRACT,abstract.strip())


    def make_fname(self, cpc_code:str='', file_extension='xml'):
        try:
            title = str(self[TITLE][0])
        except KeyError:
            self.set_title()
            title = str(self[TITLE][0])
        
        title = title.replace('/','_')[:225]
        title = title.replace('\"','')
        title = title.replace('?','xxx')
        title = title.replace('\n', ' ')
        title = title.replace(',', '')

        code = cpc_code.replace('/','-')
        code = code.replace('_','-')

        return code+'_'+self['family-id']+'_('+self[PATENT_GRANT_NUM]+')_'+title+'.'+file_extension

    def write(self, patent_dir:str, allowed_cpc_codes:set):
        fname = self.make_fname(self.has_cpc_codes(allowed_cpc_codes))
        with open(patent_dir+fname,'w',encoding='utf-8') as f:
            patent_xml = minidom.parseString(et.tostring(self.patent)).toprettyxml(indent=' ')
            f.write(patent_xml)


class EPOxml(DocMine):
    def __init__(self, patent_in_epoxml:et.Element):
        self.document = self.patent = patent_in_epoxml
        pat_num = self.document.get('country') + self.document.get('doc-number')
        super().__init__(PATENT_GRANT_NUM, pat_num)

        app_ref = self.document.find('./SDOBI/B100')
        app_country = app_ref.find('./B190').text
        app_num = app_ref.find('./B110').text
        app_kind = app_ref.find('./B130').text
        self.Identifiers[PATENT_APP_NUM] = app_country+app_num+app_kind

        self.append_property('file','')
        self.append_property('doc_number',self.document.get('doc-number'))
        self.append_property('country', self.document.get('country'))
        self.append_property('kind', self.document.get('kind'))

        
    def get_data(self):
        language_names = self.patent.findall('./SDOBI/B500/B540/B541')
        titles = self.patent.findall('./SDOBI/B500/B540/B542')
        for i in range(0, len(language_names)):
            if language_names[i].text == 'en':
                self._set_title(titles[i].text)
                break

        date = self.patent.find('./SDOBI/B100/B140/date').text
        self.set_date(date[:4],date[4:6],date[6:])

        abstract = self.patent.find('./abstract[@lang=\'en\']/p')
        abstract = abstract.text if type(abstract) != type(None) else ''
        self.add2section(ABSTRACT,abstract)

        institutions = self.patent.findall('./SDOBI/B700/B710/B711/snm')
        institutions = {x.text for x in institutions}
        institutions = cleanup_institutions(institutions)
        self.addresses = {x:x for x in institutions}
        #self.update_with_list(list(cleanup_institutions(institutions)))

        authors = self.patent.findall('./SDOBI/B700/B720/B721/snm')
        authors = {x.text for x in authors}
        self.update_with_list(AUTHORS, list(cleanup_authors(list(authors))))

        claims = self.patent.findall('./claims/claim/claim-text')
        claims = {x.text for x in claims if type(x.text) != type(None)}
        for c in claims:
                self.add2section(CLAIMS,c)

    @staticmethod
    def loads(xml_file:str):
        try:
            tree = et.parse(xml_file)
        except et.ParseError: 
            return None
        return tree.getroot()
    
    @classmethod
    def from_file(cls, xml_file:str):
        try:
            c = cls(cls.loads(xml_file))
            file = xml_file[xml_file.rfind('/')+1:]
            c['file'] = [file[file.rfind('\\')+1:]]
            return c
        except et.ParseError: 
            return None


class USPTOjson (DocMine):
    def __init__(self, result_elements:list):
        self.patent = result_elements
        try:
            doc_number = self.patent['grantDocumentIdentifier']
            id_type = PATENT_GRANT_NUM
            date = self.patent['grantDate']            
            self.append_property('Application Date', self.patent['filingDate'])
        except KeyError:
            doc_number = self.patent['patentApplicationNumber']
            id_type = PATENT_APP_NUM
            date = self.patent['filingDate']
                   
        super().__init__(id_type,doc_number)
        self.Identifiers[PATENT_APP_NUM] = self.patent['patentApplicationNumber']

        self.append_property('doc_number', doc_number[2:])
        self.append_property('country', doc_number[:2])
        self.set_date(date[-4:],date[:2],date[3:5])

    def get_cpc_codes(self):
        main_cpc = str(self.patent['mainCPCSymbolText'])
        self.append_property('CPC',main_cpc)
        further_cpcs = list(self.patent['furtherCPCSymbolArrayText'])
        self.update_with_list('CPC',further_cpcs)


    def has_cpc_codes(self, allowed_cpc_codes:set): 
        main_cpc = str(self.patent['mainCPCSymbolText']).replace('/','-')
        if main_cpc in allowed_cpc_codes: return main_cpc

        cpcs = self.patent['furtherCPCSymbolArrayText']
        if type(cpcs) == type(None): return ''
        for c in cpcs:
            cpc_for_fname = str(c).replace('/','-')
            if cpc_for_fname in allowed_cpc_codes: return cpc_for_fname
        return ''

    def set_title(self):
        self._set_title(self.patent['inventionTitle'])

    def get_data(self):
        self.set_title()
        self.inventor_names = list(self.patent['inventorNameArrayText'])
        self.update_with_list(AUTHORS, list(cleanup_authors(self.inventor_names)))
        if type (self.patent['assigneeEntityName']) != type(None):
            institutions = [self.patent['assigneeEntityName']]
            institutions = cleanup_institutions(institutions)
            self.addresses = {x:x for x in institutions}
            #self.update_with_list(INSTITUTIONS, cleanup_institutions(institutions))

        abstract = self.patent['abstractText'][0]
        self.add2section(ABSTRACT,abstract)

        ct = self.patent['claimText']
        claims_text = str(self.patent['claimText'][0]) if type (ct) != type(None) else ''
        for i in range (1,9):
            claims_text = claims_text.replace(str(i)+'. ','|')
        claims = claims_text.split('|')
        
        for c in claims:
            if c:
                claims_pretty = list()
                ctext = c[:-1] if c[-1].isdigit() else c
                split_text = ctext.split(';')
                for s in split_text:
                    more_splits = s.split(', ') if len(s) >= 3500 else [s]
                    claims_pretty += more_splits
                
                claim_text = '. '.join(claims_pretty)
                self.add2section(CLAIMS,claim_text)

    def make_fname(self, cpc_code:str='', file_extension='xml'):
        try:
            title = str(self[TITLE][0])
        except KeyError:
            self.set_title()
            title = str(self[TITLE][0])
        
        title = title.replace('/','_')[:225]
        title = title.replace('\"','')
        title = title.replace('?','xxx')
        title = title.replace('\n', ' ')
        title = title.replace(',', '')

        code = cpc_code.replace('/','-')
        code = code.replace('_','-')

        pat_grant_num = str()
        try:
            pat_grant_num = self.Identifiers[PATENT_GRANT_NUM]
        except KeyError: pass
        pat_app_num = self.Identifiers[PATENT_APP_NUM]

        if pat_grant_num:
            fname = code+'_'+pat_grant_num+'_'+title+'.'+file_extension
        else:
            fname = code+'(_'+pat_app_num+')_'+title+'.'+file_extension

        return fname

    def write(self, patent_dir:str, allowed_cpc_codes:set):
        fname = self.make_fname(self.has_cpc_codes(allowed_cpc_codes),'json')
        with open(patent_dir+fname,'w',encoding='utf-8') as f:
            f.write(json.dumps(self.patent, indent=1))

    @staticmethod
    def loads(json_file:str):
        return json.load(open(json_file))
    
    @classmethod
    def from_file(cls, json_file:str):
        c = cls(cls.loads(json_file))
        file = json_file[json_file.rfind('/')+1:]
        c['file'] = [file[file.rfind('\\')+1:]]
        return c

def has_ipc(ipc_code:str):
    if ipc_code in allowed_ipc_codes: return True
    ipc_base = ipc_code[:ipc_code.find('/')]
    return ipc_base in smoking_devices


def is_patent_start(line:str):
    if line[:23] == '<us-patent-application ': return True 
    elif line == '<patent-application-publication>': return True
    else: return False

def is_patent_end(line:str):
    if line[:24] == '</us-patent-application>': return True
    elif line == '</patent-application-publication>': return True
    else: return False

def normalizeEPOipc(ipc_code:str, skip=1):
    ipc_code = ipc_code.strip()
    slash_pos = ipc_code.find('/')
    main_class = ipc_code[skip:slash_pos]
    main_class = main_class.split(' ')
    if len(main_class) > 2:
        classification_section = main_class[0]
        classification_class = main_class[1]
        classification_main_group = main_class[-1]
    elif len(main_class) == 2:
        classification_section = main_class[0]
        classification_class = main_class[1]
        classification_main_group = ''
    else:
        print('Unknown ipc code format: %s' % ipc_code)

    if len(classification_main_group) == 0: classification_main_group = '000'                                          
    if len(classification_main_group) == 1: classification_main_group = '00'+classification_main_group
    if len(classification_main_group) == 2: classification_main_group = '0'+classification_main_group

    subgroup = ipc_code[slash_pos+1 : ].strip()
    subgroup = subgroup.split(' ')
    subgroup = subgroup[0]

    normalized_ipc_code = classification_section+classification_class+classification_main_group+'/'+subgroup
    return normalized_ipc_code

def load_docnumbers(fname:str, col_number:int=0, has_header=False):
    doc_numbers = list()
    
    with open(fname, 'r') as listing:
        if has_header: line = listing.readline()
        line = listing.readline().strip()
        while line:
            columns = line.split('\t')
            #country, doc_number, kind = OPSxml.parse_patnum(columns[0])
            doc_numbers.append(columns[col_number])
            line = listing.readline().strip()
    return doc_numbers


def do_biblio_search(cql_query):
    client = epo_ops.Client(key=OPS_APIKEY, secret=OPS_SECRET)
    try:
        # getting search result counts
        search_response = client.published_data_search(cql=cql_query)
    except HTTPError: 
        print('no patents found for cql query %s' % cql_query)
        return []

    search_response_xml = str(search_response.content,encoding='utf-8')
    tree = et.fromstring(search_response_xml)
    biblio_search = tree.find('{http://ops.epo.org}biblio-search')
    found_patents = int(biblio_search.attrib['total-result-count'])
    print('%d patents found for %s query' % (found_patents,cql_query))
    refs = biblio_search.findall('./{http://ops.epo.org}search-result/{http://ops.epo.org}publication-reference/')

    doc_ids = set()
    for pat_id in refs:
        #parsing search results
        #id_type = pat_id.attrib['document-id-type']
        #print ([elem.tag for elem in pat_id.iter()]) #to get list of tags
        country = pat_id.find('./{http://www.epo.org/exchange}country').text
        doc_number = pat_id.find('./{http://www.epo.org/exchange}doc-number').text
        kind = pat_id.find('./{http://www.epo.org/exchange}kind').text
        doc_ids.add(country+doc_number+kind)

    return doc_ids

def cache_exist_patnum(patent_dir:str,allowed_countries:set=None,with_kind=True):
    exist_docnums = set()
    for root, subdirs, files in os.walk(patent_dir):
        if files:
            print('searching in %s directory' % root)
        for xml_file in files:
            file_path = os.path.join(root, xml_file)
            try:
                tree = et.parse(file_path)
            except et.ParseError or FileNotFoundError: continue
            patent = OPSxml(tree.getroot())
            patent.get_docnum()
            if isinstance(allowed_countries, set):
                if patent.country in allowed_countries:
                    docnum = patent.pat_number+patent.kind if with_kind else patent.pat_number
                    exist_docnums.add(docnum)
            else:
                docnum = patent.pat_number+patent.kind if with_kind else patent.pat_number
                exist_docnums.add(docnum)
    return exist_docnums

def download_patent(country, doc_number, kind):
    client = epo_ops.Client(key=OPS_APIKEY, secret=OPS_SECRET)
    return client.published_data( # Retrieve bibliography data
            reference_type = 'publication', # publication, application, priority
            input = epo_ops.models.Docdb(doc_number, country, kind),  # original, docdb, epodoc
            endpoint = 'biblio',  # biblio includes authors, title, abstract, document id, year, applicant-name:
            # https://worldwide.espacenet.com/help?locale=en_EP&method=handleHelpTopic&topic=bibliographic
            constituents = [] # optional, list of constituents
    )

# The doc number is a number given to a patent application when it is filed, published or or re-registred. 
# The number contains a two digit series code followed by a six digit serial number 
# assigned by the USPTO (Example: US9748552). Document kind must not be included into the doc_number
def download_patents_from_ops(patnum_file:str, patent_dir:str,allowed_cpc_codes:set):
    #patnum_file must have document IDs for download in the first column
    #middlewares = [mw.Dogpile(),mw.Throttler()] 
    exist_docnums = cache_exist_patnum(patent_dir)
    print('Patent directory "%s" has %d patents'%(patent_dir,len(exist_docnums)))

    download_count = 0
    exist_counter = 0
    line_counter = 0
    doc_numbers = load_docnumbers(patnum_file)
    #doc_numbers = do_biblio_search(cql_query)# not implelented need to pass cql_query
    start = time.time()
    for patnum in doc_numbers:
        if patnum in exist_docnums: 
            exist_counter +=1
            continue

        country, doc_number, kind = OPSxml.parse_patnum(patnum)
        if not kind: kind = ' '
        try:
            #downloading actual patent data
            patent_data = download_patent(country, doc_number, kind)
        except HTTPError: 
            print('publication %s is not available' % (patnum))
            continue
        
        patent = OPSxml(et.fromstring(patent_data.content))
        patent.get_docnum()
        #assert(patent.pat_number == patnum)

        patent.write(patent_dir,allowed_cpc_codes)
        download_count +=1

        if download_count % 100 == 0:
            print('Download status: %d patents downloaded, %d found in current directory out of %d in "%s"' % (
                download_count,exist_counter, len(doc_numbers), patnum_file))
            print('Download time for %d patents is %s' %(download_count,OPSxml.execution_time(start)))
            start = time.time()
            sleep_time = 30
            print('will sleep for %d sec to avoid IP address block by Espacenet' %sleep_time)
            time.sleep(sleep_time)

    print('Total download: %d patents out of %d' % (download_count,line_counter))

def cleanup_patent_dir(validpatentnums:str, dir2clean:str, allowed_cpc_codes:set):
    valid_doc_numbers = load_docnumbers(validpatentnums)
    found_patnum = set()
    duplicate = dict()
    to_delete = open('to_delete.cmd', 'w', encoding='utf-8')
    
    for root, subdirs, files in os.walk(dir2clean):
        if files:
            print('searching in %s directory' % root)
        for xml_file in files:
            file_path = os.path.join(root, xml_file)
            try:
                tree = et.parse(file_path)
            except et.ParseError: continue
            patent = OPSxml(tree.getroot())
            #patent.get_docnum()

            if patent.pat_number not in valid_doc_numbers:
                to_delete.write('del %s\n'% file_path)
                continue
            
            if patent.pat_number in found_patnum:
                try:
                    duplicate[patent.pat_number].append(file_path)
                except KeyError:
                    duplicate[patent.pat_number] = [file_path]
                continue
            else:
                found_patnum.add(patent.pat_number)
                patent.write(dir2clean,allowed_cpc_codes)
            
    to_delete.close()
    with open('dups_to_delete.cmd', 'w', encoding='utf-8') as dups:
        for k,v in duplicate.items():
            for f in v:
                dups.write(k+'\t'+'del '+f+'\n')


def download_patent_from_uspto(granted_patent_num, patent_dir:str,allowed_cpc_codes:set):
    params = {'patentNumber': str(granted_patent_num), 'start':0, 'rows':100, 'largeTextSearchFlag':'N'}
    data = urllib.parse.urlencode(params, quote_via=urllib.parse.quote)
    baseURL = 'https://developer.uspto.gov/ibd-api/v1/application/grants?'
    url = baseURL+data
    response = requests.get(url)
    result = json.loads(response.text)
    if result['results']:
        patent = USPTOjson(result['results'][0])
        patent.write(patent_dir,allowed_cpc_codes)

def download_patents_from_uspto(docnum_file,column_with_id,patent_dir:str,allowed_cpc_codes:set):
    docnums = load_docnumbers(docnum_file,column_with_id,has_header=True)
    for n in docnums:
        download_patents_from_uspto(n,patent_dir,allowed_cpc_codes)

if __name__ == "__main__":
    patent_dir = 'D:/Python/Patents/Espacenet/Englis2/'
    allowed_cpc_codes = {'A61K9-0073', 'A61K9-0075', 'A61K9-0078'}
    #cleanup_patent_dir('A61K9 patent numbers all.txt',patent_dir,allowed_cpc_codes)
    #download_patents_from_ops('A61K9 patent numbers all.txt',patent_dir,allowed_cpc_codes)
    #download_patents_from_uspto('USPTO formulation patents.txt',1,'D:/Python/Patents/full patents with claims/USPTO/',allowed_cpc_codes)
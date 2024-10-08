
import time,sys,os,json,urllib.request, requests
from datetime import timedelta
from xml.dom import minidom
from urllib.parse import quote
from requests.auth import HTTPBasicAuth
from lxml import etree as et
from base64 import b64encode
DEFAULT_CONFIG_DIR = os.path.join(os.getcwd(),'ENTELLECT_API/ElsevierAPI/')
DEFAULT_APICONFIG = os.path.join(DEFAULT_CONFIG_DIR,'APIconfig.json')


def execution_time(execution_start):
    return "{}".format(str(timedelta(seconds=time.time() - execution_start)))


def execution_time2(execution_start:float,current_iteration:int,number_of_iterations:int):
    '''
    input:
        if "number_of_iterations" is supplied assumes that "execution_start" is global start
        otherwise assumes "execution_start" is the start of the current iteration if "remaining_iterations" is supplied
        
    output:
        tuple: time passed from execution_start, remaining_time
    '''
    delta = time.time() - execution_start
    time_passed = "{}".format(str(timedelta(seconds=delta)))
    remaining_iterations = number_of_iterations - current_iteration
    remaining_time = delta*float(remaining_iterations/current_iteration)
    remaining_time_str = "{}".format(str(timedelta(seconds=remaining_time)))
    return time_passed, remaining_time_str



def load_api_config(api_config_file=''):# file with your API keys and API URLs
    if not api_config_file:
        print('No API config file was specified\nWill use default %s instead'% DEFAULT_APICONFIG)
        api_config_file = DEFAULT_APICONFIG
    else:
        if not os.path.isabs(api_config_file):
            print(f'APIconfig is specified only by file name {api_config_file}')
            print(f'Will look for {api_config_file} in default "{DEFAULT_CONFIG_DIR}" directory')
            api_config_file = os.path.join(DEFAULT_CONFIG_DIR,api_config_file)

    try:
        return dict(json.load(open(api_config_file,'r')))
    except FileNotFoundError:
        print("Cannot find API config file: %s" % api_config_file)
        if api_config_file != DEFAULT_APICONFIG:
            print('Cannot open %s config file\nWill use default %s instead'% (api_config_file, DEFAULT_APICONFIG))
            return dict(json.load(open(DEFAULT_APICONFIG,'r')))
        else:
            print('No working API server was specified!!! Goodbye')
            return dict()

    

def dir2flist(path2dir:str,include_subdirs=True,subdirs_only=False,file_ext='',fnames_has:list[str]=[]):
    my_path = os.path.join(path2dir, '')
    my_files = []
    for root, dirs, files in os.walk(my_path):
        if not include_subdirs and root != my_path: continue
        if subdirs_only and root == my_path: continue
        [my_files.append(os.path.join(root, f)) for f in files if f.lower().endswith(file_ext)]
        # The endswith() function in Python will return True if the input suffix is an empty string.
        # This is because an empty string is considered to be a suffix of any string, including itself.

    return [f for f in my_files if any(s in f for s in fnames_has)] if fnames_has else my_files


def normalize_filename(name:str) -> str:
    """
    Normalizes a filename by replacing illegal characters.
    """
    replacements = {'>': '-', '<': '-', '|': '-', '/': '-', ':': '_'}
    return "".join(replacements.get(char, char) for char in name)
    

def pretty_xml(xml_string:str, remove_declaration = False):
        '''
        xml_string must have xml declration
        '''
        pretty_xml = str(minidom.parseString(xml_string).toprettyxml(indent='   '))
        return pretty_xml[pretty_xml.find('\n')+1:] if remove_declaration else pretty_xml


def file_head(full_path:str,number_of_lines = 10000):
    path, filename = os.path.split(full_path)
    filename,ext = os.path.splitext(filename)
    with open(full_path,'r',encoding='utf-8') as f:
        path2test = os.path.join(path,filename+f'-first{number_of_lines}'+ext)
        with open(path2test,'w',encoding='utf-8') as t:
            for counter in range(0,number_of_lines):
                t.write(f.readline())
    return


def unpack(list_of_lists:list,make_unique=True):
    flat_list = [item for sublist in list_of_lists for item in sublist]
    return list(set(flat_list)) if make_unique else flat_list


def next_tag(in_xml_file:str,tag:str,namespace=''):
    if namespace:
        tag = '{'+namespace+'}'+tag

    context = et.iterparse(in_xml_file,tag=tag)
    for event, elem in context:
        assert(isinstance(elem,et._Element))
        yield elem
        #yield et.tostring(elem).decode()
        elem.clear()
    return


def all_tags(element:et._Element):
    tags = set()  # Use a set to store unique tag names
    tags.add((element.tag))
    [tags.update(all_tags(child)) for child in element]
    #print(tags)
    return tags


def all_child_parents(element:et._Element):
    child2parent = {c.tag:p.tag for p in element.iter() for c in p}
    print(child2parent)
    return child2parent


def urn_encode(string:str,prefix:str):
  """
    input:
        prefix - desired string after urn: prefix
    output:
        urn:prefix:URN-encoded(string) 
  """
  encoded_string = quote(string, safe='-_.!~')
  return f"urn:{prefix}:{encoded_string}"


def sortdict(indic:dict,by_key=True,reverse=False):
    i = 0 if by_key else 1
    return dict(sorted(indic.items(), key=lambda item: item[i],reverse=reverse))


def get_auth_token(**kwargs):
    """
    kwargs:
        token_url,
        client_id,
        secret,
        username,
        password,
    output:
        token, retreival time stamp
    """
    try:
        auth = HTTPBasicAuth(kwargs['client_id'], kwargs['secret'])
    except KeyError:
        pass

    data = {'username': kwargs['username'], 'password': kwargs['password']}
    response = requests.post(kwargs['token_url'], auth=auth, data=data)

    '''
    req = urllib.request.Request(kwargs['token_url'], method='POST')
    us_pas = '{}:{}'.format(kwargs['username'],kwargs['password'])
    us_pas_e = us_pas.encode()
    base64string = b64encode(us_pas_e)
    req.add_header("Authorization", "Basic %s" % base64string.decode())
    response = urllib.request.urlopen(req)
    '''
    body = json.loads(response.read())
    token = body['access_token']
    return token, time.time()


class Tee(object):
    def __init__(self, filename, mode="w"):
        self.file = open(filename, mode,encoding='utf-8')
        self.stdout = sys.stdout

    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)

    def flush(self):
        self.file.flush()
        self.stdout.flush()

    def __enter__(self):
        sys.stdout = self  # Redirect stdout to this object

    def __exit__(self, exc_type, exc_value, traceback):
        sys.stdout = self.stdout  # Restore original stdout
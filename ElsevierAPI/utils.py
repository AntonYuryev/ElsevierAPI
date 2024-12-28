
import time,sys,os,json, requests,re,traceback
from datetime import timedelta
from xml.dom import minidom
from urllib.parse import quote
from requests.auth import HTTPBasicAuth
from lxml import etree as et
from concurrent.futures import ThreadPoolExecutor
DEFAULT_CONFIG_DIR = os.path.join(os.getcwd(),'ENTELLECT_API/ElsevierAPI/')
DEFAULT_APICONFIG = os.path.join(DEFAULT_CONFIG_DIR,'APIconfig.json')

PCT = '%'


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
  pretty_xml = "\n".join([line for line in pretty_xml.splitlines() if line.strip()])
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


def str2str(dic:dict):
    new_dic = dict()
    for k,v in dic.items():
      new_dic[k] = ';'.join(map(str,v))
    return new_dic


GREEK2ENGLISH = {
      "α": "alpha",
      "β": "beta",
      "γ": "gamma",
      "δ": "delta",
      "ε": "epsilon",
      "ζ": "zeta",
      "η": "eta",
      "θ": "theta",
      "ι": "iota",
      "κ": "kappa",
      "λ": "lambda",
      "μ": "mu",
      "ν": "nu",
      "ξ": "xi",
      "ο": "omicron",
      "π": "pi",
      "ρ": "rho",
      "σ": "sigma",
      "τ": "tau",
      "υ": "upsilon",
      "φ": "phi",
      "χ": "chi",
      "ψ": "psi",
      "ω": "omega",
      "Α": "Alpha",
      "Β": "Beta",
      "Γ": "Gamma",
      "Δ": "Delta",
      "Ε": "Epsilon",
      "Ζ": "Zeta",
      "Η": "Eta",
      "Θ": "Theta",
      "Ι": "Iota",
      "Κ": "Kappa",
      "Λ": "Lambda",
      "Μ": "Mu",
      "Ν": "Nu",
      "Ξ": "Xi",
      "Ο": "Omicron",
      "Π": "Pi",
      "Ρ": "Rho",
      "Σ": "Sigma",
      "Τ": "Tau",
      "Υ": "Upsilon",
      "Φ": "Phi",
      "Χ": "Chi",
      "Ψ": "Psi",
      "Ω": "Omega"
  }
def greek2english(text:str):
  for symbol, spelling in GREEK2ENGLISH.items():
    text = text.replace(symbol, spelling)
  return text


def normalize(s:str):
  text = greek2english(s)
  text = re.sub(r'[^a-zA-Z0-9\s]', ' ', text)  # Remove non-alphanumeric characters
  text = text.replace('  ',' ')
  return text


def tokenize(s:str):
  text = normalize(s)
  tokens = text.lower().split()  # Split into words and convert to lowercase
  return tokens


def match_tokens(tokens1:list,tokens2:list):
  if len(tokens1) == len(tokens2):
    for i, token1 in enumerate(tokens1):
      if token1 != tokens2[1]:
          return False
    return True
  else:
      return False


def get_auth_token(**kwargs):
    """
    kwargs:
        token_url,
        client_id,
        client_secret,
        username,
        password
    output:
        authorization header, retreival time stamp
    """
    try:
        auth = HTTPBasicAuth(kwargs.pop('username'), kwargs.pop('password'))
        data = {"grant_type": 'password'}
        response = requests.post(kwargs['token_url'], auth=auth, data=data)
    except KeyError:
        try:
            data = {'client_id':kwargs['client_id'],'client_secret':kwargs['client_secret']}
            data.update({"grant_type": 'client_credentials'})
            response = requests.post(kwargs['token_url'],  data=data)
        except KeyError:
            print('No valid credetials are supplied to access SBS server')
            return None
        
    body = response.json()
    token = str(body['access_token'])

    '''
    using urllib.request for some websites:
    req = urllib.request.Request(kwargs['token_url'], method='POST')
    us_pas = '{}:{}'.format(kwargs['username'],kwargs['password'])
    us_pas_e = us_pas.encode()
    base64string = b64encode(us_pas_e)
    req.add_header("Authorization", "Basic %s" % base64string.decode())
    response = urllib.request.urlopen(req)
    body = json.loads(response.read())
    '''
    return {"Authorization": "Bearer " + token}, time.time()

def print_error_info(x:Exception,thread_name =''):
  exc_type, exc_value, exc_traceback = sys.exc_info()
  traceback_list = traceback.extract_tb(exc_traceback)
  error_message = f'{thread_name} thread has finished with error "{x}":'
  for tb_info in traceback_list:
    filename = tb_info.filename
    module_name = tb_info.name
    line_number = tb_info.lineno
    error_message += f"  - File: {filename}, Function: {module_name}, Line: {line_number}\n"
  print(error_message)


def run_tasks(tasks:list):
  '''
  Executes a list of tasks concurrently using ThreadPoolExecutor
  Args:
    tasks: A list of tuples, where each tuple contains a function and a tuple of its arguments. For example:
            [(func1, (arg1, arg2)), (func2, (arg1,))]
    if function argument is a single list convert it to tuple as (my_list,)
  '''
  future_dic = {}
  with ThreadPoolExecutor() as ex:
    for func, args in tasks:
      future_dic[func.__name__] = ex.submit(func, *args)

    result_dic = dict()
    for func_name, future in future_dic.items():
      try:
        result_dic[func_name] = future.result()
      except Exception as e:
        print_error_info(e,func_name)
    return result_dic


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
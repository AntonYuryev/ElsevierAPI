
import time,sys,os,json, requests,re,traceback,urllib.request,unicodedata
from urllib.parse import quote as urlencode
from collections import Counter
from itertools import chain as iterchain
from time import sleep
from math import ceil
from typing import Generator
from datetime import timedelta,datetime
from xml.dom import minidom
from requests.auth import HTTPBasicAuth
from lxml import etree as et
from concurrent.futures import ThreadPoolExecutor,as_completed
DEFAULT_CONFIG_DIR = os.path.join(os.getcwd(),'ENTELLECT_API/ElsevierAPI/')
DEFAULT_APICONFIG = os.path.join(DEFAULT_CONFIG_DIR,'APIconfig.json')
PCT = '%'


def current_time():
  """Prints the current date and time in a human-readable format."""
  now = datetime.now()
  return now.strftime('%Y-%m-%d %H:%M:%S')


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


def path2folderlist(path2dir:str, first_folder=''):
  '''
  output:
    list of folder until first_folder, no first_folder in the list
  '''
  drive, path_without_drive = os.path.splitdrive(path2dir)
  folders = []
  while True:
    path_without_drive, tail = os.path.split(path_without_drive)
    if tail == first_folder:
        break
    folders.insert(0, tail) # Insert at the beginning to maintain order
  return folders


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


def remove_duplicates(items:list):
  '''
    keeps uids order in uids
  '''
  seen = set()
  return [i for i in items if i not in seen and not seen.add(i)]


def unpack(list_of_lists:list[list]|list[tuple],make_unique=True):
  if make_unique:
    return list(dict.fromkeys(iterchain.from_iterable(list_of_lists)))
  else:
    return list(iterchain.from_iterable(list_of_lists))
  #  flat_list = [item for sublist in list_of_lists for item in sublist]


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
  encoded_string = urlencode(string, safe='-_.!~')
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

pattern = re.compile("|".join(GREEK2ENGLISH.keys()))
def greek2english(text: str) -> str:
  return pattern.sub(lambda m: GREEK2ENGLISH[re.escape(m.group(0))], text)

'''
def greek2english(text:str):
  for symbol, spelling in GREEK2ENGLISH.items():
    text = text.replace(symbol, spelling)
  return text
'''

def replace_non_unicode(text:str):
  normalized_text = ''.join(c for c in unicodedata.normalize('NFKD', text) if unicodedata.category(c) != 'Mn')
  return normalized_text


def normalize(s:str):
  text = greek2english(s)
  text = re.sub(r'[^a-zA-Z0-9\s]', ' ', text)  # Remove non-alphanumeric characters
  text = replace_non_unicode(text)
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
    return {"Authorization": "Bearer " + token}, time.time()

def print_error_info(x:Exception,thread_name =''):
  exc_type, exc_value, exc_traceback = sys.exc_info()
  traceback_list = traceback.extract_tb(exc_traceback)
  error_message = f'{thread_name} thread has finished with error "{x}"\n{x.__doc__}:'
  for tb_info in traceback_list:
    filename = tb_info.filename
    module_name = tb_info.name
    line_number = tb_info.lineno
    error_message += f"  - File: {filename}, Function: {module_name}, Line: {line_number}\n"
  print(error_message)


def run_tasks(tasks:list)->dict:
  '''
  Executes a list of tasks concurrently using ThreadPoolExecutor
  Args:
    tasks: A list of tuples, where each tuple contains a function and a tuple of its arguments. For example:
            [(func1, (arg1, arg2)), (func2, (arg1,))]
    if function has one argument convert it to tuple as (my_list,)
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
  

def multithread(big_list:list, func, **kwargs):
  '''
  input:
    big_list: list of items to be processed
    func: function to be applied to each item in the list
    max_workers: number of threads to use for processing
  output:
    list of results from applying func to each item in big_list
  '''
  max_workers = kwargs.pop('max_workers', 10)
  def chunk_func(chunk:list):
    return [func(item,**kwargs) for item in chunk]
  
  results = []
  with ThreadPoolExecutor(max_workers=max_workers) as ex:
    futures = [ex.submit(chunk_func, chunk) for i,chunk in list2chunks_generator(big_list,max_workers)]   
    [results.extend(f.result()) for f in as_completed(futures)]
  return results

'''
def list2chunks(input_list:list, num_chunks:int=0,chunk_size:int=0):
  list_length = len(input_list)
  if chunk_size:
    num_chunks = list_length // chunk_size
    if num_chunks*chunk_size < list_length:
      num_chunks += 1
  else:
    assert(num_chunks),"Either num_chunks or chunk_size must be specified"
    chunk_size = list_length // num_chunks # floor division

  remainder = list_length % num_chunks
  chunks = []
  start = 0
  for i in range(num_chunks):
    end = start + chunk_size
    if i < remainder:
      end += 1  # Distribute the remainder among the first 'remainder' chunks
    chunks.append(input_list[start:end])
    start = end
  return chunks
'''

def list2chunks_generator(input_list:list, num_chunks:int=0, chunk_size:int=0)->Generator[tuple[int,list],None,None]:
    '''
      hint: Generator[YieldType, SendType, ReturnType]
    '''
    list_length = len(input_list)
    if chunk_size:
      num_chunks = ceil(list_length / chunk_size)
    else:
      assert num_chunks > 0, "Either num_chunks or chunk_size must be specified and positive."

    if list_length == 0:
      return 0,[]# Yields nothing for empty list
    if num_chunks == 0:
      yield 0,input_list # Yield the whole list as one chunk

    base_chunk_size = list_length // num_chunks
    remainder = list_length % num_chunks
    start = 0
    for i in range(num_chunks):
      real_chunk_size = base_chunk_size + int(i < remainder) # distributing remainder among chunks
      if real_chunk_size == 0 and start < list_length:
        continue

      end = start + real_chunk_size
      yield i, input_list[start:end]
      start = end
      if start >= list_length:
        break


def bisect(data_list:list, criterion):
    """
    input:
      data_list - list of objects in ascending order
      criterion - function applied to the element in the list. Must return bool

    Returns:
      index of the first element object with true criterion.
    """
    low = 0
    high = len(data_list) - 1
    bisector_index = -1

    while low <= high:
      mid = (low + high) // 2
      if criterion(data_list[mid]):
        bisector_index = mid
        high = mid - 1 # This element's criterion is true, search lower
      else:
        low = mid + 1  # This element's criterion is not true, search higher

    return bisector_index


def atempt_request4(url:str,retries=10,sleeptime=5):
  req = urllib.request.Request(url=url)
  for attempt in range(retries):
    try:
      response = urllib.request.urlopen(req).read()
      return response
    except Exception as e:
      print(f'{e} on attempt {attempt} out of {retries} to obtain {url}')
      sleep(sleeptime)
      continue
  return ''



def most_frequent(data:list,if_data_empty_return=''):
    """
    output:
      most frequently occurring value in a list.
    """
    if not data:
        return None
    counts = Counter(data)
    most_common = counts.most_common(1)
    return most_common[0][0]


class Tee(object):
    def __init__(self, filename, mode="w"):
        self.file = open(filename, mode,encoding='utf-8')
        self.stdout = sys.stdout
        print(f'Runtime messages will be in {filename}')

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
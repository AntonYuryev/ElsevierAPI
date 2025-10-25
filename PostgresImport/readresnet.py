#!/bin/python3 -u

import logging, sys, re, os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from datetime import datetime
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import ElementTree
from lxml import etree as et
from timeit import default_timer as timer
from threading import Lock as lock

from src.cachingwriter import CachingWriter
from src.logging import start_logging
from src.myhash import myhash

# xml header
xmlheader = '<?xml version="1.0" encoding="UTF-8" standalone="no" ?>\n'
# tracks total database writes
linesread = 0
# approx total lines in current release
totalines = 1447002036
# reproducible seed
readversion = False

# fields in reference table
refcolumns = [
    "unique_id",
    "Authors",
    "BiomarkerType",
    "CellLineName",
    "CellObject",
    "CellType",
    "ChangeType",
    "Collaborator",
    "Company",
    "Condition",
    "DOI",
    "EMBASE",
    "c",
    "ESSN",
    "Experimental System",
    "Intervention",
    "ISSN",
    "Journal",
    "MedlineTA",
    "Mode of Action",
    "mref",
    #"msrc",
    "Sentence/msrc",
    "NCT ID",
    "Organ",
    "Organism",
    "Percent",
    "Phase",
    "Phenotype",
    "PII",
    "PMID",
    "PubVersion",
    "PubYear",
    "PUI",
    "pX",
    "QuantitativeType",
    "Source",
    "Start",
    "StudyType",
    "TextMods",
    "TextRef",
    "Tissue",
    "Title",
    "TrialStatus",
    "URL",
    "id",
]

versioncache = CachingWriter("version", False)
attrcache = CachingWriter("attr")
controlcache = CachingWriter("control")
refcache = CachingWriter("reference", False)  # no caching
nodecache = CachingWriter("node")
pathcache = CachingWriter("pathway")

start = timer()


# make refmapper
def makerefmap():
    result = {}
    count = 0
    for ref in refcolumns:
        names = ref.split('/')
        for n in names:
            result[n] = count
        count += 1

    return result


# make reference array, return an array
def makeref(refhash):
  result = [""] * len(refcolumns)
  result[0] = refhash
  return result

# refmap is map of reference items a map of item, index in array
refmap = makerefmap()


def refid(ref_row:list):
  '''
  output:
    reference identifier from TextRef
  '''
  ref = ref_row[refmap["TextRef"]].split("#")[0]
  return ref if re.search(r"\d", ref) else str(ref_row[0])


def parseVersion(xml):
    """write version records to version table"""
    with open('whatisgoingon.xml', 'w') as x:
        x.write(xml)
    resnet_elem = ElementTree(ET.fromstring(xml))
    with open("version.table", "w") as f:
        h = resnet_elem.findall(".//attr")
        for h in h:
            name = h.get("name")
            value = h.get("value")
            val = (name, value)
            logging.info("name:%s\tval:%s" % val)

            if name == "ReleaseDate":
                try:
                    date = datetime.strptime(value, "%B %d, %Y")
                except:
                    try:
                        date = datetime.strptime(value, "%b %d, %Y")
                    except:
                        date = None
                versioncache.write(("update", date), f)

            versioncache.write(val, f)


def indexAttribute(attr:ET.Element):
  '''
  output:
    hcode, name, value, index for attr
  '''
  # attr is an XML document
  name = attr.get("name")
  value = attr.get("value")
  index_str = attr.get("index")
  index = 0 if index_str is None else int(index_str)
  hcode = myhash(name + "|" + value)
  return hcode, name, value, index


def parseNode(node:ET.Element):
    attributes = []
    nodeRef = []
    nodeLocalId = {}
    urn = node.get("urn")
    local_id = node.get("local_id")
    nodehash = myhash(urn)
    nodeLocalId = {local_id: nodehash}
    nodeName = ""
    nodeType = ""

    for nodeattr in node.findall("./attr"):
      attr_tuple = indexAttribute(nodeattr)
      if attr_tuple[1] == "Name":
        nodeName = attr_tuple[2]
      elif attr_tuple[1] == "NodeType":
        nodeType = attr_tuple[2]
      else:
        hcode1, name1, value1,_ = attr_tuple
        attributes.append((hcode1, name1, value1))
        nodeRef.append(hcode1)

    thisnode = (nodehash, urn, nodeName, nodeType, nodeRef)
    return thisnode, attributes, nodeLocalId


def parseResnet(resnet_xml:str):
  # cut this big section out
  if "<attachments>" in resnet_xml:
    resnet_xml = re.sub("<attachments>.*</attachments>", "", resnet_xml, flags=re.M)

  resnet_elem = ElementTree(ET.fromstring(resnet_xml))
  resnet_root = resnet_elem.getroot()
  nodes = []
  attributes = []
  controls = []
  pathways = []
  isPathway = False
  refhash = -1
  invalid_rnef_attrs = set()
  references = []
  resnetHashes = []

  resnet_type = resnet_root.get("type")
  if resnet_type: # <resnet> is either Pathway or Group
    if resnet_type == "Pathway":
      isPathway = True
      for attr in resnet_elem.findall("./properties/attr"):
        attr_tuple = indexAttribute(attr)
        resnetHashes.append(attr_tuple[0])
        # keep those attributes connected to resnet properties
        hcode1, name1, value1,_ = attr_tuple
        attributes.append((hcode1, name1, value1)) # add to list to store in attr table
    else: # if not a pathway bail out early. resnet can be Group
      return [], [], [], [], [], {} #attributes, nodes, controls, [], pathways, invalid_attrs

  nodeLocalId = {}
  for node in resnet_elem.findall("./nodes/node"):  # node
    nodeRef = []
    urn = node.get("urn")
    local_id = node.get("local_id")
    nodehash = myhash(urn)
    nodeLocalId[local_id] = nodehash
    nodeName = ""
    nodeType = ""

    for nodeattr in node.findall("./attr"):
      attr_tuple = indexAttribute(nodeattr)
      if attr_tuple[1] == "Name":
        nodeName = attr_tuple[2]
      elif attr_tuple[1] == "NodeType":
        nodeType = attr_tuple[2]
      else:
        hcode1, name1, value1,_ = attr_tuple
        attributes.append((hcode1, name1, value1))
        nodeRef.append(hcode1)

      node_tuple = (nodehash, urn, nodeName, nodeType, nodeRef)
      nodes.append(node_tuple)

  controlHashes = []

  for control in resnet_elem.findall("./controls/control"):  # controls
    inref = []
    outref = []
    inoutref = []
    controlType = ""
    ontology = ""
    relationship = ""
    mechanism = ""
    effect = ""

    for controllink in control.findall("./link"):  # links
      ref = controllink.get("ref")
      ty = controllink.get("type")

      if ty == "in":
        inref.append(nodeLocalId[ref])
      elif ty == "out":
        outref.append(nodeLocalId[ref])
      elif ty == "in-out":
        inoutref.append(nodeLocalId[ref])
      else:
        logging.warning(f"Unknown link type found: {ty}")

    # refhash = myhash( str(inref, outref, inoutref) + str(random.random())
    #refhash = -1
    # tried using xml of the control for this snippet but generating the XML is very slow.
    # in some cases there may be more than one  for in and out
    # however these may be only for the lipidomics project.  If required
    # the data type could be arrays instead of integers
    localrefs = []
    index0based = False
    for controlattr in control.findall("./attr"):  # attributes of the control
      _, name, value, index = indexAttribute(controlattr)
      # these are specific to the control, and not to any indivitual references
      if name == "ControlType":
        controlType = value
      elif name == "Effect":
        effect = value
      elif name == "Mechanism":
        mechanism = value
      elif name == "Ontology":
        ontology = value
      elif name == "Relationship":
        relationship = value
      else:
        # index is an attribute of the control that says which reference# it belongs to
        # it can be 0-based or 1-based 
        # expand list to fit index
        if index == 0:
          index0based = True
      
        if index >= len(localrefs):
          extendby = index - len(localrefs) + 1
          new_refs = [makeref(refhash) for _ in range(extendby)] # list of lists
          localrefs.extend(new_refs) # add to list of references for this relationship
        xref = localrefs[index]
        try:
          xref[refmap[name]] = value
        except KeyError:
          #logging.error(e)
          #logging.error(f"{hcode}, {name}, {value}, {index}, {localrefs}")
          invalid_rnef_attrs.add(name)

    if not index0based:
      localrefs = localrefs[1:]
    # end of loop over attributes

    # use the absolute references for hash, not 'local' to enable combining unique controls
    chash = myhash(
        str(
            (
                inref,
                inoutref,
                outref,
                controlType,
                ontology,
                relationship,
                effect,
                mechanism,
            )
        )
    )
    
    control_tuple = (
        chash,
        inref,
        inoutref,
        outref,
        controlType,
        ontology,
        relationship,
        effect,
        mechanism,
        chash,
    )
    # assign non-unique hashes for s
    # now localrefs is an array of arrays, where we we need an array of tuples
    
    for i,refrow in enumerate(localrefs):
      refrow[-1] = chash # last data element links to control
      refrow[0] = ""
      refrow[0] = myhash(str(refrow)) # unique hash of this reference
      refrow.append(refid(refrow)) # extract ref identifier from TextRef to count unqiue references per relation
      localrefs[i] = tuple(refrow)

    [references.append(r) for r in localrefs]
    # end of this control in resnet_elem.findall("./controls/control")
    #
    # use the absolute references for hash, not 'local' to enable combining unique controls
    controlHashes.append(chash)
    controls.append(control_tuple)
  # end of loop over controls

  if isPathway:
    rurn = resnet_root.get("urn")
    rname = resnet_root.get("name")
    phash = myhash(str((refhash, rname, resnet_type, rurn, resnetHashes, controlHashes)))
    p = (phash, rname, resnet_type, rurn, resnetHashes, controlHashes)
    pathways.append(p)

  return attributes, nodes, controls, references, pathways, invalid_rnef_attrs


def precord(resnet:str,fname:str='',invalid_rnef_attrs:set={}):
  """parse record and create db records"""
  attributes, nodes, controls, references, pathways, invalid_attrs = parseResnet(resnet)
  invalid_rnef_attrs.update(invalid_attrs)

  with lock():
    with open(fname+"_control.table", "a",encoding='utf-8') as f:
      [controlcache.write(i, f) for i in controls]

    with open(fname+"_reference.table", "a",encoding='utf-8') as f:  # references is an array of arrays
      [refcache.write(i, f) for i in references]

    with open(fname+"_attr.table", "a",encoding='utf-8') as f:
      [attrcache.write(i, f) for i in attributes]

    with open(fname+"_node.table", "a",encoding='utf-8') as f:
      [nodecache.write(i, f) for i in nodes]

    with open(fname+"_pathway.table", "a",encoding='utf-8') as f:
      [pathcache.write(i, f) for i in pathways]

  return len(controls), len(pathways)


def read_resnet(path2rnef:str):
  logging.info(f"reading file {path2rnef}")
  fname = Path(path2rnef).stem
  start = datetime.now()
  max_workers = min(32,(os.cpu_count() or 1) + 4)
  with ThreadPoolExecutor(max_workers,thread_name_prefix='ParseRNEF') as executor:
    futures = []
    invalid_rnef_attrs = set()
    resnet_counter = 0
    context = et.iterparse(path2rnef, tag="resnet", recover=True) #huge_tree=True
    for action, elem in context:
      xml_str = et.tostring(elem).decode()
      future = executor.submit(precord,xml_str,fname,invalid_rnef_attrs)
      futures.append(future)
      #control_count, pathway_count = precord(xml_str,fname,invalid_rnef_attrs)
      resnet_counter += 1
      if resnet_counter > 0 and resnet_counter%5000 == 0:
        elapsed_time = datetime.now()-start
        print(f'Processed {resnet_counter} resnet sections in {elapsed_time}')
      elem.clear()
      while elem.getprevious() is not None:
        del elem.getparent()[0]
    del context

    control_counter = 0
    pathway_counter  = 0
    for f in as_completed(futures):
      control_count, pathway_count = f.result()
      control_counter += control_count
      pathway_counter += pathway_count
      if control_counter > 0 and control_counter%5000 == 0:
        elapsed_time = datetime.now()-start
        print(f'Processed {control_counter} controls and {pathway_counter} pathways in {elapsed_time}')

  print('Invalid attributes in RNEF:')
  print(invalid_rnef_attrs)
  logging.info("completed")
  attrcache.stats()
  controlcache.stats()
  refcache.stats()
  nodecache.stats()
  pathcache.stats()


def main(fname:str = ''):
  start_logging(folder=".", file=Path(fname).stem)
  global readversion
  if not fname:
    fname = sys.argv[1]
    if len(sys.argv) > 2 and len(sys.argv[2]) > 0:
      char = (sys.argv[2].lower())[0]
      readversion = char in ("t", "1")
      logging.info(f"readversion set to: {readversion}")

  read_resnet(fname)

#renin-angiot-system.rnef
#resnet18_mammal_07152025.rnef
main('C:/ResnetDump/PostgresTestData/resnet18_mammal_07152025.rnef')
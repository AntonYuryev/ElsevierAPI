import urllib.request, urllib.parse, os, time
from lxml import etree as ET
from urllib.error import HTTPError
from time import sleep
from ..utils import dir2flist, execution_time,pretty_xml,next_tag,dir2flist,replace_non_unicode,urn_encode
from ..ResnetAPI.NetworkxObjects import PSObject,PSRelation,AUTHORS,JOURNAL,PUBYEAR
from ..ResnetAPI.ResnetGraph import ResnetGraph,Reference,TITLE,SENTENCE,OBJECT_TYPE
from collections import defaultdict


BASE_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
CLINVAR_CACHE_DIR = os.path.join(os.getcwd(),'ENTELLECT_API/ElsevierAPI/NCBI/__clinvarcache__/')
CLINVAR_CACHE_BASE = 'clinvar'
CLINVAR_CLASSIFICATION = 'ClinVar classification'
CLINVAR_RNEF_PROPS = ['Name','Alias','Experimental System',CLINVAR_CLASSIFICATION]

class cvSNP(PSObject):
   def __init__(self,name:str,id:str,database='dbSNP'):
       super().__init__()
       self.set_property('Name',replace_non_unicode(name))
       self.Idendifires = {database:id}
       urn = 'urn:agi-gv-'+database.lower()+':'+id
       self.set_property('URN',urn)
       self.MAFs = dict() # 


def cache_path(basename=CLINVAR_CACHE_BASE):
    cache_files = dir2flist(CLINVAR_CACHE_DIR,fnames_has=basename,file_ext='.xml')
    return os.path.join(CLINVAR_CACHE_DIR,CLINVAR_CACHE_BASE+str(len(cache_files))+'.xml')


def gv2SNP(gv:PSObject,rsid2SNP:dict[str,cvSNP]):
    my_rsids = [gv.name()]+list(gv.get_props(['Alias']))
    for rsid in my_rsids:
      if rsid.startswith('rs'):
          return rsid2SNP[rsid]
    return cvSNP()


def rs2rcv(rsids:list):
  stepSize = 200
  rsids_len = len(rsids)
  rcv_ids = set()
  for i in range(0, rsids_len, stepSize):
    rsids_chunk = ','.join(s[2:] for s in rsids[i:i+stepSize])
    elink_params = {'id':rsids_chunk,'dbfrom':'snp','db':'clinvar'}
    req_l = urllib.request.Request(url=BASE_URL+'elink.fcgi?'+urllib.parse.urlencode(elink_params))
    # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?&db=clinvar&dbfrom=snp&id=
    response_l = urllib.request.urlopen(req_l).read()
    link2cv = ET.fromstring(response_l.strip())
    cvids = [e.text for e in link2cv.findall('LinkSet/LinkSetDb/Link/Id')]
    
  for i in range(0, len(cvids), stepSize):
    cvids_chunk = ','.join(s for s in cvids[i:i+stepSize])
    efecth_params = {'db':'clinvar','id':cvids_chunk}
    req = urllib.request.Request(url=BASE_URL+'esummary.fcgi?'+urllib.parse.urlencode(efecth_params))
    # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?&db=clinvar&id=
    for attempt in range(1,11):
      try:
          response = urllib.request.urlopen(req).read()
          break
      except HTTPError as e:
          print(f'Fetch EUtils attempt {attempt} failed with error {e}.  Will try again in 5 seconds')
          sleep(5)
          continue

    cvs = ET.fromstring(response.strip())
    rcvs = [e.text for e in cvs.findall('DocumentSummarySet/DocumentSummary/supporting_submissions/rcv/string')]
    rcv_ids.update(rcvs)
    print(f'Downloaded {i+stepSize} SNPs out of {len(cvids)}')
    sleep(1)
  return rcv_ids


def parse_measure(measure:ET._Element,only_rsids=set(),assertion_id='')->tuple[PSObject,PSObject,set]:
    '''
    output:
      gv, gene, molecular_consequences
    '''
    empty_return = tuple(PSObject(),PSObject(),set())
    gv = PSObject({OBJECT_TYPE:['GeneticVariant']})
    name_tag = measure.find('./Name')
    if name_tag is not None:
      gv_name = name_tag.text.strip()
      if gv_name:
        gv.update_with_value('Name',replace_non_unicode(gv_name))

      elemval_tag = name_tag.find('./ElementValue')
      if elemval_tag is not None:
        gvs_code = elemval_tag.text
        gv.update_with_value('Alias',replace_non_unicode(gvs_code))

    rs = ''
    gv_ids = dict()
    gv_xrefs = measure.findall('./XRef')
    for xref in gv_xrefs:
      db = str(xref.get('DB'))
      if db  == 'dbSNP':
        rs = 'rs'+xref.get('ID')
      else:
        identifier = xref.get('ID')
        db = db2attr(db)
        gv.update_with_value(db,identifier)
        gv_ids[db] = identifier
    
    if only_rsids and rs not in only_rsids:
      return empty_return

    if not gv.name():
      if rs:
        gv.update_with_value('Name',rs)
      elif gv_ids:
        db, identifier = next(iter(gv_ids.items()))
        gv.update_with_value('Name',db+':'+identifier)
      else:
        gv.update_with_value('Name',gv.get_prop("Alias"))

    if rs:
      gv.update_with_value('URN',urn_encode(rs,'agi-gv-dbsnp'))
    elif gv_ids:
        db,identifier = next(iter(gv_ids.items()))
        gv.update_with_value('URN',urn_encode(identifier,f'clinvar-{db}'))
    else:
        gvname = gv.name()
        if gvname:
          gv.update_with_value('URN',urn_encode(gv.name(),f'clinvar'))
        else:
          print(f'cannot create URN for SNP in  Clinvar record {assertion_id}')
          return empty_return


    gene_name = ''
    geneid = str()
    measured_relationship = measure.find('./MeasureRelationship')
    gene_ids = dict()
    if measured_relationship is None:
      elemval = measure.find('./Name/ElementValue')
      if elemval is not None:
        gene_name = elemval.text
      else:
        attr = measure.find('./AttributeSet/Attribute')
        if attr is not None:
          gene_name = attr.get('Accession')
    else:
      gene_name = measured_relationship.find('./Symbol/ElementValue').text
      assert(isinstance(measured_relationship,ET._Element))
      if measured_relationship.get('Type','') == 'within single gene':
        gene_xrefs = measured_relationship.findall('./XRef')
        for xref in gene_xrefs:
          assert(isinstance(xref,ET._Element))
          db = xref.get('DB','')
          if db == 'Gene':
            geneid = xref.get('ID')
          else:
            gene_ids[db] = xref.get('ID')

    if gene_name:
      gene = PSObject({'Name':[gene_name],OBJECT_TYPE:['Protein']})
      if geneid:
          gene.update_with_value('URN',urn_encode(geneid,'agi-llid'))
          gene.update_with_value('LocusLink ID',geneid)
      elif gene_ids:
          db,identifier = next(iter(gene_ids.items()))
          gene.update_with_value('URN',urn_encode(identifier,f'clinvar-{db}'))
      else:
          gene.update_with_value('URN',urn_encode(gene_name,'agi-prot'))
    else:
       gene = PSObject()

    molecular_consequences = {e.text for e in measure.findall('./AttributeSet/Attribute') if e.get('Type','') == 'MolecularConsequence'}
    
    return gv,gene,molecular_consequences


@staticmethod
def db2attr(db:str):
  normdb = str(db)
  comma_pos = normdb.find(',')
  if comma_pos != -1:
    normdb = normdb[comma_pos+1:].strip()
    normdb = normdb.replace('&amp;','')
    normdb = normdb.replace('&','')
  return replace_non_unicode(normdb[:47]+' ID') # max attribute name length = 50 


def parse_trait(assertion:ET._Element):
    trait = assertion.find('./TraitSet/Trait')
    assert(isinstance(trait,ET._Element))
    trait_type = trait.get('Type','')

    traid_ids = defaultdict(set)
    [traid_ids[db2attr(xr.get('DB'))].add(xr.get('ID')) for xr in trait.findall('./XRef') if xr.get('Type') != 'Allelic variant']
    [traid_ids[db2attr(xr.get('DB'))].add(xr.get('ID')) for xr in trait.findall('./Name/XRef') if xr.get('Type') != 'Allelic variant']

    aliases = [n.find('./ElementValue').text for n in trait.findall('./Name')]
    aliases = [a for a in aliases if a not in  ['not provided','not specified','AllHighlyPenetrant']]
    [aliases.append(n.find('./ElementValue').text) for n in trait.findall('./Symbol')]

    if aliases:
      trait_name_s = aliases[0]

      if trait_type in ['Disease','Finding','NamedProteinVariant']:
        objtype = 'Disease'
      elif trait_type in ['DrugResponse','BloodGroup']:
        objtype = 'CellProcess'
      elif trait_type == 'PhenotypeInstruction':
        objtype = ''
      else:
        print(f'Unknown type "{trait_type}" for trait "{trait_name_s}"')
        debug_info = pretty_xml(ET.tostring(assertion).decode())
        print(debug_info)

      if objtype:
        disease = PSObject({'Name':[trait_name_s],
                            'Alias':aliases,
                            'URN':[urn_encode(trait_name_s,'clinvar-disease')],
                            OBJECT_TYPE : [objtype]})
        disease.update({k:list(v) for k,v in traid_ids.items()})

        sentences = list()
        for elem in trait.findall('./AttributeSet/Attribute'):
          elem_type = elem.get('Type')
          if elem_type in ['public definition','disease mechanism']:
            s = elem.text
            if s != 'not provided':
              sentences.append(s)

        return disease, sentences # disease has trait IDs
      else:
       return PSObject(),[]
    else:
      return PSObject(),[]


def rcv2psrels(ClinVarSet:ET._Element,mapdic:dict[str,dict[str,dict[str,PSObject]]],only4rsids=set(),include_benign=False):
  '''
  input:
    only4rsids - filter by dbSNP rs identifiers
    mapdic - {OBJECT_TYPE:{propname:{provalue:PSObject}}}
  output:
    rel_fa:PSRelation = GV-FunctionalAssociation-Disease
    rel_gc:PSRelation = GV-GeneticChange-Gene
    where GV, Disease, Gene are mapped by mapdic
  '''
  gv = PSObject()
  disease = PSObject()
  gene = PSObject()
  rel_fa = PSRelation()
  rel_gc = PSRelation()

  title = ClinVarSet.find('./Title').text
  assertion = ClinVarSet.find('.ReferenceClinVarAssertion')
  assertion_id = assertion.get('ID')
  
  classifications = {x.text for x in assertion.findall('Classifications/GermlineClassification/Description')}
  #ClinVar contains benign and unknown significancs SNPs that do not have disease associations
  #Therefore gv2gene relations are added only if phenotypes exist for GV
  # and gv germline_classification does not contain "benign"
  if not include_benign:
    for c in classifications:
      if c in ['Benign','Likely benign','Uncertain significance']:
        return PSRelation(),PSRelation()

  measure_elem = assertion.find('./MeasureSet/Measure')
  if measure_elem is None:
    measure_elem = assertion.find('./GenotypeSet/MeasureSet/Measure')
  gv,gene,molecular_consequences = parse_measure(measure_elem,only4rsids,assertion_id)
  if not gv:
    return rel_fa,rel_gc
  
  if classifications:
    gv[CLINVAR_CLASSIFICATION] = list(classifications)

  method_elem = assertion.find('.ObservedIn/Method/MethodType')
  method = '' if method_elem is None else method_elem.text
  observed_data = assertion.find('.ObservedIn/ObservedData')
  sentence = ''
  if observed_data is None:
    gvdis_ref = Reference('Clinvar',assertion_id)
    gvgene_ref = Reference('Clinvar',assertion_id)
  else:
    citation = observed_data.find('./Citation')
    if citation is None:
      gvdis_ref = Reference('Clinvar',assertion_id)
      gvgene_ref = Reference('Clinvar',assertion_id)
    else:
      sentence = observed_data.find('./Attribute').text
      if sentence == 'not provided':
        sentence = ''
      
      pubmed_id = ''
      ref_ids = dict()

      for citation_id in citation.findall('./ID'):
        citation_source = citation_id.get('Source')
        if citation_source == 'PubMed':
          pubmed_id = citation_id.text
          break
        else:
          ref_ids[citation_source] = citation_id.text
      if pubmed_id:
        gvdis_ref = Reference('PMID',citation_id.text)
        gvdis_ref.Identifiers.update(ref_ids)
        gvgene_ref = Reference('PMID',citation_id.text)
        gvdis_ref[TITLE] = [title]
        gvgene_ref[TITLE] = [title]
      elif ref_ids:
        id_type,identifier = next(iter(ref_ids.items()))
        gvdis_ref = Reference(id_type,identifier)
        gvgene_ref = Reference(id_type,identifier)
        gvdis_ref[TITLE] = [title]
        gvgene_ref[TITLE] = [title]
      else:
        gvdis_ref = Reference('Clinvar',assertion_id)
        gvgene_ref = Reference('Clinvar',assertion_id)

  disease,sentences = parse_trait(assertion)
  if sentence:
    sentences.append(sentence)

  if disease:
    disease.remap(mapdic,['Name','Alias'])
      
    gvdis_ref.Identifiers['Clinvar'] = assertion_id
    gvdis_ref['Clinvar ID'] = [assertion_id]
    text_ref = gvdis_ref._make_textref()
    if TITLE not in gvdis_ref:
      gvdis_ref[TITLE] = ['ClinVar: public archive of interpretations of clinically relevant variants']
    gvdis_ref[AUTHORS] = ['Landrum, M.J. et al.']
    gvdis_ref[JOURNAL] = ['Nucleic Acids Res.']
    gvdis_ref[PUBYEAR] = [2016]
    gvd_snippet = dict()
    if sentence:
      gvd_snippet[SENTENCE] = set(sentences)
    else:
      gvd_snippet[SENTENCE] = {'Record from ClinVar (NCBI) includes this association between SNP (SNV) and phenotype'}
    
    if method:
      gvd_snippet['Experimental System'] = {method}
    if gvd_snippet:
      gvd_snippet['Source'] = {'ClinVar'}
      gvdis_ref.add_snippet(text_ref,gvd_snippet)
    refs = [gvdis_ref]

    dict4rel_fa = defaultdict(list)
    dict4rel_fa[OBJECT_TYPE].append('FunctionalAssociation')
    rel_fa = PSRelation.make_rel(gv,disease,dict(dict4rel_fa),refs,is_directional=False)

  if gene:
    gene.remap(mapdic,['LocusLink ID','Name','Alias'])
    mol_consequences = '. '.join(molecular_consequences)
    if mol_consequences:
      gvgene_ref.Identifiers['Clinvar'] = assertion_id
      text_ref = gvgene_ref._make_textref()
      gvg_snippet = dict()
      gvg_snippet[SENTENCE] = {mol_consequences}
      if gvg_snippet:
        gvg_snippet['Source'] = {'ClinVar'}
        gvgene_ref.add_snippet(text_ref,gvg_snippet)

    if TITLE not in gvdis_ref:
      gvgene_ref[TITLE] = ['ClinVar: public archive of interpretations of clinically relevant variants']
    gvgene_ref[AUTHORS] = ['Landrum, M.J. et al.']
    gvgene_ref[JOURNAL] = ['Nucleic Acids Res.']
    gvgene_ref[PUBYEAR] = [2016]

    dict4rel_gc = defaultdict(list)
    dict4rel_gc[OBJECT_TYPE].append('GeneticChange')
    rel_gc = PSRelation.make_rel(gv,gene,dict(dict4rel_gc),[gvgene_ref],is_directional=True)
  
  return rel_fa,rel_gc

   
def rcv2rn(ClinVarResultSet:ET._Element, mapdic:dict[str,dict[str,dict[str,PSObject]]], only4rsids:list,include_benign=False):
  all_rels = list()
  for ClinVarSet in ClinVarResultSet.findall('ClinVarSet'):
    rel_fa,rel_gc = rcv2psrels(ClinVarSet,mapdic,only4rsids,include_benign)
    if rel_fa:
      all_rels.append(rel_fa)
      if rel_gc:
        all_rels.append(rel_gc)

  clinvar_graph = ResnetGraph()
  clinvar_graph.add_psrels(all_rels)
  return clinvar_graph


def __rcvcache2rn(dump_dir:str,rsids:list,mapdic:dict[str, dict[str, dict[str, PSObject]]],include_benign=False):
  '''
    Reads ClinVar cache
  '''
  listing = dir2flist(dump_dir,file_ext='xml')
  clinvar_graph = ResnetGraph()
  for dump in listing:
    for record in next_tag(dump,'ClinVarResult-Set'):
      record_G = rcv2rn(record,mapdic,rsids,include_benign)
      clinvar_graph = clinvar_graph.compose(record_G)

  return clinvar_graph,mapdic


def downloadCV(_4rsids:list,mapdic:dict[str, dict[str, dict[str, PSObject]]],include_benign=False)->ResnetGraph:
    """
    input:
        list of rs identifiers
    output:
        ResnetGraph::gene2gv2dis 
    """
    start = time.time()
    gene2gv2dis,mapdic = __rcvcache2rn(CLINVAR_CACHE_DIR,_4rsids,mapdic,include_benign)
    assert(isinstance(gene2gv2dis,ResnetGraph))

    gvs = gene2gv2dis._psobjs_with('GeneticVariant',OBJECT_TYPE)
    rsids2download = list(set(_4rsids).difference(ResnetGraph.names(gvs)))
    print(f'{len(_4rsids)-len(rsids2download)} records were found in clinvar cache')
    rsids_len = len(rsids2download)
    if rsids_len:
      rcv_ids = list(rs2rcv(rsids2download))
      if rcv_ids:
        with open(cache_path(),'w',encoding='utf-8') as f:
          f.write('<batch>')
          params = {'db':'clinvar','rettype':'clinvarset'}
          stepSize = 200
          for i in range(0, len(rcv_ids), stepSize):
            ids = ','.join(rcv_ids[i:i+stepSize])
            params.update({'id':ids})
            req = urllib.request.Request(url=BASE_URL+'efetch.fcgi?'+urllib.parse.urlencode(params))
            for attempt in range(1,11):
              try:
                response = urllib.request.urlopen(req).read()
                break
              except HTTPError as e:
                print(f'Fetch EUtils attempt {attempt} failed with error {e}.  Will try again in 5 seconds')
                sleep(5)
                continue
              
            cvs = ET.fromstring(response.strip())
            f.write(pretty_xml(ET.tostring(cvs),True))
            chunk_G = rcv2rn(cvs,mapdic,rsids2download,include_benign)
            gene2gv2dis.add_graph(chunk_G)
          f.write('</batch>')
    print(f'Download was done in {execution_time(start)}')
    return gene2gv2dis


def cv_entprops4rnef(gene2gv2dis:ResnetGraph):
  diseases = gene2gv2dis._psobjs_with('Disease',OBJECT_TYPE)
  disease_dbis = {p for d in diseases for p in d.keys() if p.endswith(' ID')}
  return CLINVAR_RNEF_PROPS+list(disease_dbis)
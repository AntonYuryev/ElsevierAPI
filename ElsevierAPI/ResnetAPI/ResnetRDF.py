from argparse import Namespace
import stat
import rdflib as rdf
from .ResnetGraph import ResnetGraph
from .NetworkxObjects import PSObject,PSRelation,REFCOUNT
from ..ETM_API.references import Reference,PS_ID_TYPES, RELEVANCE
import json
from pyld import jsonld
from urllib.parse import quote
from itertools import combinations


EDM = rdf.Namespace("https://data.elsevier.com/schema/edm/")
HGRAPH = rdf.Namespace("https://data.elsevier.com/health/core/ex/")
SCHEMA = rdf.Namespace("https://schema.org/")
PUBMED = rdf.Namespace("https://pubmed.ncbi.nlm.nih.gov/")
DOI = rdf.Namespace("https://dx.doi.org/")
DBSNP = rdf.Namespace("https://www.ncbi.nlm.nih.gov/snp/")
LOINC = rdf.Namespace("https://loinc.org/")
PII = rdf.Namespace('https://www.sciencedirect.com/science/article/abs/pii/')
EMBASE = rdf.Namespace('https://www.embase.com/records?subaction=viewrecord&id=')
CT = rdf.Namespace('https://www.clinicaltrials.gov/ct2/show/')
GENE =  rdf.Namespace('https://www.ncbi.nlm.nih.gov/gene/')
GO =  rdf.Namespace('http://purl.obolibrary.org/obo/go_')
CAS =  rdf.Namespace('http://cas.org/')
PCCID =  rdf.Namespace('http://rdf.ncbi.nlm.nih.gov/pubchem/endpoint/cid')
MESH =  rdf.Namespace('hhttps://meshb.nlm.nih.gov/record/ui?name=mesh')
EC =  rdf.Namespace('https://enzyme.expasy.org/EC/ec')
UNIPROT = rdf.Namespace('https://www.uniprot.org/core')
RESNET = rdf.Namespace('https://mammal-profservices.pathwaystudio.com/')
ETM = rdf.Namespace('https://demo.elseviertextmining.com/')

context = {'edm':EDM,
                'schema':SCHEMA,
                'pubmed':PUBMED,
                'doi':DOI,
                'dbsnp':DBSNP,
                'loinc':LOINC,
                'hgraph':HGRAPH,
                'cas':CAS,
                'trial':CT,
                'embase':EMBASE,
                'go':GO,
                'uniprot':UNIPROT,
                'mesh':MESH,
                'ecnum':EC,
                'resnet':RESNET,
                'etm':ETM
                }
                
REFID2NAMESPACE = {'PMID':PUBMED, 'DOI':DOI, 'PII':PII, 'EMBASE':EMBASE, 'NCT ID':CT}

class ResnetRDF(rdf.Graph):
    pass
    def __init__(self):
        super().__init__()
        [self.bind(k,v) for k,v in context.items()]


    def add_reference(self,ref:Reference, to_reluri:str):
        was_added = False
        for i in PS_ID_TYPES:
            try:
                id_value = ref.Identifiers[i]
                namespace = REFID2NAMESPACE[i]
                ref_uri = namespace[id_value]
                ref_rel_uri = to_reluri+':'+i+':'+id_value
                was_added = True
                break
            except KeyError: continue
        
        if not was_added:
            try:
                title_uri = quote(ref['Title'][0])
                ref_uri = rdf.URIRef('https://scholar.google.com/'+title_uri)
                ref_rel_uri = to_reluri+':title:'+title_uri
            except KeyError:
                textref_uri = quote(ref.Identifiers['TextRef'])
                ref_uri = RESNET[textref_uri]
                ref_rel_uri = to_reluri+':textref:'+textref_uri

        self.add((to_reluri,SCHEMA['CreativeWork'],ref_rel_uri))
        self.add((ref_rel_uri,SCHEMA.identifier,ref_uri))
            #   [self.add(ref_uri, SCHEMA['author'],a) for a in ref[AUTHORS]]
            #   self.add((ref_uri, SCHEMA['title'],ref[TITLE][0]))
            #   self.add((ref_uri, SCHEMA['datePublished'],ref[PUBYEAR][0]))
            #   self.add((ref_uri, SCHEMA['periodical'],ref[JOURNAL][0]))
        try:
            relevance = ref[RELEVANCE][0]
            self.add((ref_rel_uri, ETM["relevance"], rdf.Literal(relevance, datatype=rdf.XSD.double)))
        except KeyError: pass
        try:
            pubtype = ref['PubTypes'][0]
            self.add((ref_rel_uri, ETM["pubtype"], rdf.Literal(pubtype)))
        except KeyError: pass


    def add_identifier(self,to_uri,id2namespace:dict):
        for id, nmspc in id2namespace.items():
            self.add((to_uri, SCHEMA['sameAs'], nmspc[id]))
    
    def add_identifiers_from_prop(self, from_property:str,for_obj:PSObject, in_namespc:str):
        obj_uri = self.__obj_uri(for_obj)
        try:
            for prop_val in for_obj[from_property]:
                identifier = quote(str(prop_val))
                self.add((obj_uri, SCHEMA['sameAs'], in_namespc[identifier]))
        except KeyError: return


    @staticmethod
    def __obj_uri(obj:PSObject): return RESNET[obj['URN'][0]]
    @staticmethod
    def __rel_uri(rel:PSRelation,resnet:ResnetGraph): return RESNET[resnet.rel_urn(rel)]

    def add_obj_prop(self,from_property:str,for_obj:PSObject,as_rdftype:str):
        obj_uri = self.__obj_uri(for_obj)
        try:
            for prop_val in for_obj[from_property]:
                annotation = quote(str(prop_val))
                self.add((obj_uri, as_rdftype, rdf.Literal(annotation)))
        except KeyError: return


    def add_rel_prop(self,from_property:str,for_rel:PSRelation,as_rdftype:str,resnet:ResnetGraph):
        obj_uri = self.__rel_uri(for_rel,resnet)
        for prop_val in for_rel[from_property]:
            annotation = quote(str(prop_val))
            self.add((obj_uri, as_rdftype, rdf.Literal(annotation)))


    def add_psobj(self,node:PSObject):
        self.add_obj_prop('ObjTypeName',node,rdf.RDF.type)
        self.add_obj_prop('Name',node,SCHEMA.name)
        self.add_obj_prop('Gene',node,RESNET['Gene'])
        self.add_identifiers_from_prop('LOING ID',node, LOINC)
        self.add_identifiers_from_prop('hGraph ID',node, ETM)
        return self.__obj_uri(node)


    def add_psrel(self,rel:PSRelation, resnet:ResnetGraph):
        self.add_rel_prop('ObjTypeName',rel,rdf.RDF.type,resnet)
        resnet.rel_name(rel)
        self.add_rel_prop('Name',rel,SCHEMA.name,resnet)
        self.add_rel_prop(REFCOUNT,rel,RESNET['ref_count'],resnet)
        [self.add_rel_prop(x,rel,RESNET[x],resnet) for x in rel.keys() if str(x).find('Type') >= 0] 
        # covers ChangeType, BiomarkerType, PubTypes
        return self.__rel_uri(rel,resnet)


    def add_triple(self, rel:PSRelation, resnet:ResnetGraph):
        regulators, targets = resnet.find_nodes(rel)
        rel_uri = self.add_psrel(rel, resnet)
        for reg in regulators:
            reg_uri = self.add_psobj(reg)
            for targ in targets:
                targ_uri = self.add_psobj(targ)
                self.add((reg_uri, rel_uri, targ_uri))
                
        if not targets:
            reg_pairs = combinations(regulators,2)
            for pair in reg_pairs:
                node1uri = self.__obj_uri(pair[0])
                node2uri = self.__obj_uri(pair[1])
                self.add((node1uri, rel_uri, node2uri))

        [self.add_reference(ref,rel_uri) for ref in rel._get_refs()]


    def load_resnet(self, resnet:ResnetGraph):
        for regulatorID,targetID,rel in resnet.edges.data('relation'):
            self.add_triple(rel,resnet)
            

    def to_json(self,fname,format="json-ld"):
        with open(fname, "w", encoding='utf-8') as f:
            byte_str = self.serialize(format=format,context=context,encoding='utf-8')
            json_dict = json.loads(byte_str.decode("utf-8"))
            json_compacted = jsonld.compact(json_dict, context)
            json.dump(json_compacted, f, indent=2, sort_keys=True,ensure_ascii=False)

    @classmethod
    def fromResnetGraph(cls,g:ResnetGraph):
        rdf = cls()
        rdf.load_resnet(g)
        return rdf
        

    def get_jsonld_id(self, node:PSObject):
        urn = node['URN'][0]
        end = 11
        if urn[8:end] == 'go:': return SCHEMA.identifier, GO[urn[end:]]

        end = 12
        if urn[8:end] == 'cas:': return SCHEMA.identifier, CAS[urn[end:]]
        elif urn[8:end] == 'enz:': return SCHEMA.identifier, EC[urn[end:]]

        end = 13
        if urn[8:end] == 'llid:': return SCHEMA.identifier, GENE[urn[end:]]
        elif urn[8:end] == 'smol:': return SCHEMA.substance, SCHEMA[self['Name'][0]]

        end = 14
        if urn[8:end]  == 'pccid:': return SCHEMA.identifier, PUBMED[urn[end:]]
        
        if urn[8:15]  == 'protfc:': return 'functional_class', RESNET[self['Name'][0]]

        if urn[8:16]  == 'meshdis:': return SCHEMA.identifier, MESH[urn[16:]]
        if urn[8:16]  == 'disease:': return rdf.RDF.MedicalCondition, RESNET[self['Name'][0]]
        
        if urn[8:17]  == 'gv-dbsnp:': return SCHEMA.identifier, GO[urn[17:]]
        if urn[8:18]  == 'gocomplex:': return SCHEMA.identifier, GO[urn[18:]]
        if urn[8:19]  == 'gocellproc:': return SCHEMA.identifier, GO[urn[19:]]
        
        return SCHEMA.name, SCHEMA[self['Name'][0]]

            

        

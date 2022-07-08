import rdflib as rdf
from .ResnetGraph import ResnetGraph
from .NetworkxObjects import PSObject,PSRelation,REFCOUNT
from ..ETM_API.references import Reference,PS_ID_TYPES,RELEVANCE,TITLE,PUBYEAR,JOURNAL,LOINCID,THRESHOLD, hGRAPHID,EDMID
import json
from pyld import jsonld
from urllib.parse import quote
from itertools import combinations
'''
The following namespaces are available by directly importing from rdflib:
BRICK,CSVW,DC,DCAT,DCMITYPE,DCTERMS,DCAM,DOAP,FOAF,ODRL2,ORG,OWL,PROF,PROV,
QB,RDF,RDFS,SDO,SH,SKOS,SOSA,SSN,TIME,VANN,VOID,XSD
'''
                
REFID2PREFIX = {'PMID':'pubmed', 'DOI':'doi', 'PII':'pii', 'EMBASE':'embase', 
                'NCT ID':'clintrial', 'PUI':'embase'}
REL_PROPS4RDF = ['BiomarkerType','ChangeType','Mechanism','PubTypes']
SAMEAS2PREFIX = {LOINCID:'loinc', hGRAPHID:'edm', EDMID:'edm'}


class ResnetRDF(rdf.Graph):
    pass
    def __init__(self):
        super().__init__()
        context_dict = dict(json.load(open('ElsevierAPI/ResnetAPI/rdf/ResnetRDFcontext.json')))
        self.context = {k:rdf.Namespace(v) for k,v in context_dict.items()}
        [self.bind(k,v) for k,v in self.context.items()]
        self.resnet = ResnetGraph()


    def get_context(self):
        return {t[0]:t[1] for t in self.namespaces()}


    def make_uri(self, prefix, term):
        try:
            namespc_url = self.context[prefix]
            return namespc_url[term]
        except KeyError:
            return rdf.URIRef(prefix + term)


    def __resnet_uri(self,term):
        return self.make_uri('resnet',term)

    def __edm_uri(self,term):
        return self.make_uri('edm',term)
    

    def urn2uri(self,node:PSObject):
        urn = node['URN'][0]
        end = 11
        if urn[8:end] == 'go:': return self.make_uri('go',urn[end:])

        end = 12
        if urn[8:end] == 'cas:': return self.make_uri('cas',urn[end:])
        elif urn[8:end] == 'enz:': return self.make_uri('ecnum',urn[end:])

        end = 13
        if urn[8:end] == 'llid:': return self.make_uri('entrezgeneid',urn[end:])
        #elif urn[8:end] == 'smol:': return RESNET[node['URN'][0]

        end = 14
        if urn[8:end]  == 'pccid:': return self.make_uri('pubmed',urn[end:])
        
        #if urn[8:15]  == 'protfc:': RESNET[node['URN'][0]

        if urn[8:16]  == 'meshdis:': 
            return self.make_uri('mesh',urn[16:])
        #if urn[8:16]  == 'disease:': return RESNET[node['Name'][0]
        
        if urn[8:17]  == 'gv-dbsnp:': return self.make_uri('dbsnp',urn[17:])
        if urn[8:18]  == 'gocomplex:': return self.make_uri('go',urn[18:])
        if urn[8:19]  == 'gocellproc:': return self.make_uri('go',urn[19:])
        
        return ''
        #return rdf.URIRef('http://www.google.com/'+quote(node['Name'][0]))


    def add_reference(self,ref:Reference, to_rel_uri:str):
        was_added = False
        for i in PS_ID_TYPES:
            try:
                id_value = ref.Identifiers[i]
                ref_uri =  self.make_uri(REFID2PREFIX[i], id_value)
                ref_rel_uri = to_rel_uri+':'+i.replace(' ','_')+':'+id_value
                # uri cannot have whitespaces
                was_added = True
                break
            except KeyError: continue
        
        if not was_added:
            try:
                title_uri = quote(ref['Title'][0])
                ref_uri = rdf.URIRef('https://scholar.google.com/'+title_uri)
                ref_rel_uri = to_rel_uri+':title:'+title_uri
            except KeyError:
                return
                textref_uri = quote(ref.Identifiers['TextRef'])
                ref_uri = self.__resnet_uri(textref_uri)
                ref_rel_uri = to_rel_uri+':textref:'+textref_uri
        
        self.add((to_rel_uri,self.__resnet_uri('HasEvidence'),ref_rel_uri))
        self.add((ref_rel_uri,rdf.RDF.type, self.__resnet_uri('Evidence')))
        self.add((ref_rel_uri, self.make_uri('schema','reference'), ref_uri))
        self.add((ref_uri, rdf.RDF.type, self.make_uri('schema','CreativeWork')))
        self.add((ref_uri, rdf.RDFS.seeAlso, rdf.Literal(ref_uri,datatype=rdf.XSD.anyURI)))
        try:
            self.add((ref_uri, self.make_uri('schema','title'),rdf.Literal(ref[TITLE][0])))
        except KeyError: pass
        try:
            self.add((ref_uri, self.make_uri('schema','datePublished'),rdf.Literal(ref[PUBYEAR][0])))
        except KeyError: pass
        try:
            self.add((ref_uri, self.make_uri('schema','periodical'),rdf.Literal(ref[JOURNAL][0])))
        except KeyError: pass
    #   [self.add(ref_uri, SCHEMA['author'],a) for a in ref[AUTHORS]] 
           
        try:
            relevance = ref[RELEVANCE][0]
            self.add((ref_rel_uri, self.make_uri('etm','relevance'), rdf.Literal(relevance, datatype=rdf.XSD.double)))
        except KeyError: pass
        try:
            pubtype = ref['PubTypes'][0]
            self.add((ref_rel_uri, self.make_uri('etm','pubtype'), rdf.Literal(pubtype)))
        except KeyError: pass


    def __obj_uri(self, obj:PSObject): 
        try:
            return self.__edm_uri(obj['EDM ID'][0])#rdf.Literal(self.__edm_uri(obj['EDM ID'][0]))
        except KeyError:
            return self.__resnet_uri(obj['URN'][0])


    def __rel_uri(self, rel:PSRelation,resnet:ResnetGraph): 
        rel_urn = resnet.rel_urn(rel)
        rel_uri = self.__resnet_uri(rel_urn)
        return rel_uri
 

    def add_obj_prop(self,from_property:str,for_obj:PSObject,as_rdftype,
                    using_namespace:rdf.Namespace=''):
        obj_uri = self.__obj_uri(for_obj)

        if isinstance(using_namespace, rdf.Namespace):
            try:
                for prop_val in for_obj[from_property]:
                    annotation = str(prop_val)
                    self.add((obj_uri, as_rdftype, using_namespace[quote(annotation)]))
                return
            except KeyError: return

        try:
            for prop_val in for_obj[from_property]:
                annotation = str(prop_val)
                self.add((obj_uri, as_rdftype, rdf.Literal(annotation)))
                return
        except KeyError: return


    def add_maf(self, gv:PSObject):
        obj_uri = self.__obj_uri(gv)
        gv_name = gv['Name'][0]
        try:
            for allele2freq in gv['Minor allele frequencies']:
                sep_pos = str(allele2freq).find('-',1)
                allele = allele2freq[:sep_pos]
                freq = float(allele2freq[sep_pos+1:])
                allele_uri = self.make_uri('dbsnp',gv_name+'-'+allele)
                self.add((obj_uri, self.make_uri('dbsnp','allele'), allele_uri))
                self.add((allele_uri, self.make_uri('dbsnp','MAF'), rdf.Literal(freq, datatype=rdf.XSD.float)))
                return
        except KeyError: return


    def obj_sameas(self,for_obj:PSObject):
        obj_uri = self.__obj_uri(for_obj)
        urnuri = self.urn2uri(for_obj)
        if urnuri:
            self.add((obj_uri, rdf.OWL.sameAs, urnuri))

        for prop_id, prefix in SAMEAS2PREFIX.items():
            try:
                for prop_val in for_obj[prop_id]:
                    sameas_uri = self.make_uri(prefix,quote(str(prop_val)))
                    self.add((obj_uri, rdf.OWL.sameAs, sameas_uri))
            except KeyError: continue
        return obj_uri


    def add_psobj(self,node:PSObject):
        obj_uri = self.obj_sameas(node)
        self.add_obj_prop('ObjTypeName',node,rdf.RDF.type,self.context['resnet'])
        self.add_obj_prop('Name',node,self.context['schema']['name'])
        self.add_obj_prop('Gene',node,self.__resnet_uri('gene'),self.context['resnet'])
        self.add_obj_prop(THRESHOLD,node,self.__resnet_uri('threshold'),self.context['resnet'])
        self.add_maf(node)
        return obj_uri
        

    def add_rel_prop(self,from_property:str,for_rel:PSRelation,as_rdftype:str,using_namespace:rdf.Namespace=''):
        if isinstance(using_namespace,rdf.Namespace):
            try:
                rel_prop = for_rel[from_property]
                obj_uri = self.__rel_uri(for_rel,self.resnet)
                for prop_val in rel_prop:
                    annotation = str(prop_val)
                    self.add((obj_uri, as_rdftype, rdf.Literal(using_namespace[quote(annotation)])))
                    return
            except KeyError: return
        
        try:
            rel_prop = for_rel[from_property]
            obj_uri = self.__rel_uri(for_rel,self.resnet)
            for prop_val in rel_prop:
                annotation = str(prop_val)
                self.add((obj_uri, as_rdftype, rdf.Literal(annotation)))
                return
        except KeyError: return


    def add_psrel(self,rel:PSRelation):
        rel_uri = self.__rel_uri(rel,self.resnet)
       # self.add_rel_prop('ObjTypeName',rel,rdf.RDF.type,resnet,self.context['resnet'])
        rel_obj_type = str(rel['ObjTypeName'][0])
        rel_obj_type = rel_obj_type+'Relation' if rel_obj_type[0].isupper() else rel_obj_type+'_relation'
        self.add((rel_uri, rdf.RDF.type, self.__resnet_uri(rel_obj_type)))

        self.resnet.rel_name(rel)
        self.add_rel_prop('Name',rel,self.context['schema']['name'])

        ref_count = rel[REFCOUNT][0]
        self.add((rel_uri, self.__resnet_uri('ref_count'), rdf.Literal(ref_count, datatype=rdf.XSD.integer)))

        [self.add_rel_prop(t,rel,self.__resnet_uri(t)) for t in REL_PROPS4RDF]
        return rel_uri


    def add_triple(self, rel:PSRelation):
        regulators, targets = self.resnet.find_nodes(rel)
        rel.load_references()
        
        rel_uri = self.add_psrel(rel)
        self.add((rel_uri, rdf.RDF.predicate, self.__resnet_uri(rel['ObjTypeName'][0])))
        self.add((rel_uri, self.make_uri('schema','description'), rdf.Literal(self.triple_description(regulators,rel,targets))))

        if targets:
            for reg in regulators:
                reg_uri = self.add_psobj(reg)
                self.add((rel_uri, rdf.RDF.subject, reg_uri))
                for targ in targets:
                    targ_uri = self.add_psobj(targ)
                    self.add((rel_uri, rdf.RDF.object, targ_uri))
        else:
            reg_pairs = list(combinations(regulators,2))
            for pair in reg_pairs:
                node1uri = self.add_psobj(pair[0])
                node2uri = self.add_psobj(pair[1])
                self.add((rel_uri, rdf.RDF.object, node1uri))
                self.add((rel_uri, rdf.RDF.object, node2uri))

        [self.add_reference(ref, rel_uri) for ref in rel._get_refs()]

    
    @staticmethod
    def triple_description(regulators:list,rel:PSRelation,targets=[]):
        rel_type = rel['ObjTypeName'][0]

        if rel_type == 'Biomarker': 
            disease_name = ','.join([x['Name'][0] for x in regulators])
            biomarker_name = ','.join([x['Name'][0] for x in targets])
            return biomarker_name + ' - biomarker of ' + disease_name
        if rel_type == 'GeneticChange': 
            disease_name = ','.join([x['Name'][0] for x in regulators])
            biomarker_name = ','.join([x['Name'][0] for x in targets])
            return biomarker_name + ' - genetic biomarker of ' + disease_name

        if rel_type == 'has_diagnostic_procedure':
            disease_name = ','.join([x['Name'][0] for x in regulators])
            medprocs_name = ','.join([x['Name'][0] for x in targets])
            return medprocs_name + ' is diagnostic procedure for ' + disease_name

        if rel_type == 'has_treatment_procedure':
            disease_name = ','.join([x['Name'][0] for x in regulators])
            medprocs_name = ','.join([x['Name'][0] for x in targets])
            return medprocs_name + ' is treatment for ' + disease_name
        

        if rel_type == 'FunctionalAssociation':
            if targets:
                disease_name = ','.join([x['Name'][0] for x in regulators])
                risk_factor = ','.join([x['Name'][0] for x in targets])
                return risk_factor+' is risk factor for '+disease_name

            if regulators[0]['ObjTypeName'][0] == 'Disease': 
                disease_name = regulators[0]['Name'][0]
                gv_name = regulators[1]['Name'][0]
            else:
                disease_name = regulators[1]['Name'][0]
                gv_name = regulators[0]['Name'][0]
            return gv_name + ' is linked to ' + disease_name

        return NotImplemented


    def to_json(self,fname,format="json-ld"):
        ctx = self.get_context()
        with open(fname, "w", encoding='utf-8') as f:
            byte_str = self.serialize(format=format,context=ctx,encoding='utf-8')
            json_dict = json.loads(byte_str.decode("utf-8"))
            json_compacted = jsonld.compact(json_dict, ctx)
            json.dump(json_compacted, f, indent=2, sort_keys=True,ensure_ascii=False)


    def to_jsons(self,format="json-ld"):
        ctx = self.get_context()
        byte_str = self.serialize(format=format,context=ctx,encoding='utf-8')
        json_dict = json.loads(byte_str.decode("utf-8"))
        json_compacted = jsonld.compact(json_dict, ctx)
        return json.dumps(json_compacted, indent=2, sort_keys=True,ensure_ascii=False)


    @classmethod
    def fromResnetGraph(cls,g:ResnetGraph):
        rdf = cls()
        rdf.resnet = g
        rdf.load_resnet()
        return rdf


    def load_resnet(self):
        for regulatorID,targetID,rel in self.resnet.edges.data('relation'):
            self.add_triple(rel)

            

        

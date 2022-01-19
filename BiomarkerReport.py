from ElsevierAPI import load_api_config
from ElsevierAPI.ResnetAPI.SemanticSearch import SemanticSearch
import ElsevierAPI.ResnetAPI.PathwayStudioGOQL as OQL
from ElsevierAPI.ResnetAPI.NetworkxObjects import PSObject
import pandas as pd
import numpy as np
import time
import networkx as nx

BIOMARKER_RELATIONS = ['Biomarker','StateChange','QuantitativeChange']
DISEASES = ['Disease','ClinicalParameter']
SPECIFICITY = 'Promiscuity'

class BiomarkersReport(SemanticSearch):
    pass
    
    def __init__(self, APIconfig):
        super().__init__(APIconfig)
        self.PageSize = 500
        self.add_rel_props(['Name','Sentence','PubYear','Title','RelationNumberOfReferences'])
        self.add_ent_props(['Name'])
        

    def find_diseases(self, disease_name:str):
        mapping_prop = 'Name'
        prop2obj = self.map_prop2entities([disease_name],mapping_prop,DISEASES,get_childs=True)
        if not prop2obj:
            mapping_prop = 'Alias'
            prop2obj = self.map_prop2entities([disease_name],mapping_prop,DISEASES,get_childs=True)
        
        # prop2obj = {propValue:{id:psObj}}
        if prop2obj:
            mapped_diseases = list(prop2obj[disease_name].values())
            self.Disease = PSObject(mapped_diseases[0])
        else:
            print('No diseases were found for %s' % disease_name)
            return

        disease_ids = self.Disease[self.__child_ids__] + self.Disease['Id']
        oql_query = 'SELECT Entity WHERE id = ({ids})'.format(ids=','.join(map(str,disease_ids)))
        request_name = 'Find children for {disease}'.format(disease=self.Disease['Name'][0])
        disease_graph = self.process_oql(oql_query, request_name)
        self.diseases = disease_graph.get_objects(DISEASES)

    def load_graph(self, disease_ids:list):
        oql_query = OQL.expand_entity_by_id(disease_ids,BIOMARKER_RELATIONS,['Protein','Complex','FunctionalClass','Cell'])
        self.process_oql(oql_query, 'Find protein and cell biomakers')

        select_metabolites = '(SELECT Entity WHERE objectType=SmallMol AND Class=\'Endogenous compound\')'
        oql_query = 'SELECT Relation WHERE NeighborOf (SELECT Entity WHERE id = ({ids})) AND NeighborOf {metabolites} AND objectType = ({biomarkers})'
        oql_query = oql_query.format(ids=','.join(map(str,disease_ids)),metabolites=select_metabolites, biomarkers=','.join(BIOMARKER_RELATIONS))
        self.process_oql(oql_query, 'Find metabolite biomakers')

        biomarker_ids = [y['Id'][0] for x,y in self.Graph.nodes(data=True) if y['ObjTypeName'][0] not in DISEASES]
        disease_ids =  [y['Id'][0] for x,y in self.Graph.nodes(data=True) if y['ObjTypeName'][0] in DISEASES]
        oql_query = OQL.connect_ids(biomarker_ids,disease_ids,['Regulation'])
        self.process_oql(oql_query, 'Find Regulation between biomarkers and diseases')


    def init_semantic_search (self, reset_pandas = False):
        BiomarkerNames =  [y['Name'][0] for x,y in self.Graph.nodes(data=True) if y['ObjTypeName'][0] not in DISEASES]
        BiomarkerScores = pd.DataFrame()
        BiomarkerScores['Name'] = np.array(BiomarkerNames)
        print('Will score %d biomarkers linked to %s' % (len(BiomarkerScores),self.Disease['Name'][0]))
        self.load_pandas(BiomarkerScores,prop_names_in_header=True) #maps entities in ['Name'] column

    def semantic_search(self):
        disease_name2ids = dict()
        for disease in self.diseases:
            dis_name = disease['Name'][0]
            disease_ids = list(self._get_obj_ids_by_props(disease['Id'],["Id"],only_obj_types=DISEASES)) #finds children for disease
            disease_name2ids[dis_name] = disease_ids

        for disease in self.diseases:
            dis_name = disease['Name'][0]
            self.set_how2connect(BIOMARKER_RELATIONS+['Regulation'],[],'')
            linked_entities_count = self.link2concept(dis_name,disease_name2ids[dis_name],relations=self.Graph)
            print('%d biomarkers are linked to %s' % (linked_entities_count,dis_name))


    def biomarker_specificity(self,biomarker_ids:tuple):
        oql_query = 'SELECT Entity WHERE Connected by (SELECT Relation WHERE objectType = ({rel_types})) to (SELECT Entity WHERE id = ({id}))'
        oql_query = oql_query.format(rel_types=','.join(BIOMARKER_RELATIONS), id=','.join(map(str,list(biomarker_ids))))
        disease_count = self.get_result_size(oql_query)
        return disease_count

    def add_specificity(self):
        biomarker_count = len(self.RefCountPandas.index)
        print('Finding biomarker specificity for %d biomarkers' % biomarker_count)
        start_specificity = time.time()
        self.RefCountPandas[SPECIFICITY] = self.RefCountPandas[self.__temp_id_col__].apply(lambda x: self.biomarker_specificity(x))
        print ('Specificity for %d biomarkers was found in %s' %(biomarker_count, self.execution_time(start_specificity)))

        #for idx in self.RefCountPandas.index:
         #   biomarker_ids = list(self.RefCountPandas.at[idx,self.__temp_id_col__])
          #  specificity = self.biomarker_specificity(biomarker_ids)
           # self.RefCountPandas.at[idx,'Specificity'] = specificity

 
    def weighted_counts(self, add_specificity=False):
        self.weight_prop = 'ObjTypeName'
        self.weight_dict = {'Biomarker':1.0,'StateChange':0.5,'QuantitativeChange':0.5,'Regulation':0.25}
        self.__print_refs__ = False

        self.semantic_search()
        self.RefCountPandas.drop(columns=[self.__mapped_by__,'Resnet name'],inplace=True)
        input_disease_col = self._col_name_prefix+self.Disease['Name'][0]
        if add_specificity:
            self.add_specificity()
            normalized_score_col = 'Normalized score'
            self.RefCountPandas[normalized_score_col] = self.RefCountPandas[input_disease_col]/self.RefCountPandas[SPECIFICITY]
            self.RefCountPandas = self.RefCountPandas.sort_values(by=[normalized_score_col],ascending=False)
            self.RefCountPandas.insert(2, normalized_score_col, self.RefCountPandas.pop(normalized_score_col))
            self.RefCountPandas.insert(3, input_disease_col, self.RefCountPandas.pop(input_disease_col))
            self.RefCountPandas.insert(4, SPECIFICITY, self.RefCountPandas.pop(SPECIFICITY))
        else:
            self.RefCountPandas = self.RefCountPandas.sort_values(by=[input_disease_col],ascending=False)
            self.RefCountPandas.insert(2, input_disease_col, self.RefCountPandas.pop(input_disease_col))
        

if __name__ == "__main__":
    start_time = time.time()
    APIconfig = load_api_config()
    bm = BiomarkersReport(APIconfig)
    bm.flush_dump()

    bm.find_diseases('Fibrosis')
    t_n = bm.Disease['Name'][0]
    disease_ids =  [y['Id'][0] for x,y in bm.Graph.nodes(data=True) if y['ObjTypeName'][0] in DISEASES]
    bm.load_graph(disease_ids)
    bm.init_semantic_search()

    bm.semantic_search()
    count_file = t_n+" biomarkers counts.tsv"
    ref_file = t_n+" biomarker references.tsv"
    bm.print_ref_count(count_file,referencesOut=ref_file,sep='\t')

    bm.drop_refcount_columns()
    bm.weighted_counts(add_specificity=True)
    count_file = t_n+" biomarkers weighted counts.tsv"
    bm.print_ref_count(count_file,sep='\t')

    print('Biomarkers for %s was found in %s' % (t_n, bm.execution_time(start_time)))
    print ('Semantic counts for each biomarker are in %s' % count_file)
    print('References supporting semantic counts are in %s' % ref_file)

    '''
    NormalizedCount = bm.normalize_counts()
    fout = t_n + ' biomarker normalized report.tsv'
    NormalizedCount.to_csv(fout, sep='\t', index=False,float_format='%g')
    print ('Ranked indications are in %s' % fout)
    '''
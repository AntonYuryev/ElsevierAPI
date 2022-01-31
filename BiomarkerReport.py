from os import sep
from ElsevierAPI import load_api_config
from ElsevierAPI.ResnetAPI.SemanticSearch import SemanticSearch
import ElsevierAPI.ResnetAPI.PathwayStudioGOQL as OQL
from ElsevierAPI.ResnetAPI.NetworkxObjects import PSObject
import pandas as pd
import numpy as np
import time

QUANTITATIVE_BIOMARKER_RELS = ['Biomarker','StateChange','QuantitativeChange','CellExpression']
GENETIC_BIOMARKER_RELS = ['GeneticChange']
DISEASES = ['Disease','ClinicalParameter']
QUANTITATIVE_BIOMARKERS = ['Protein','Complex','FunctionalClass','Cell']
GENETIC_BIOMARKERS = ['Protein']
SPECIFICITY = 'Promiscuity'
SOLUBLE_BIOMARKERS_TISSUES = ['blood', 'plasma', 'serum']

class BiomarkersReport(SemanticSearch):
    pass
    
    def __init__(self, APIconfig):
        super().__init__(APIconfig)
        self.PageSize = 500
        self.add_rel_props(['Name','Sentence','PubYear','Title','RelationNumberOfReferences'])
        self.add_ent_props(['Name'])
        self.score_genetic_biomarkers = False
        self.only_soluble = False

    def find_only_soluble(self):
        self.only_soluble = True
        self.score_genetic_biomarkers = False
        self.add_rel_props(['Tissue'])
        

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
        if self.score_genetic_biomarkers:
            biomarker_types = GENETIC_BIOMARKERS 
            biomarker_rel_types = GENETIC_BIOMARKER_RELS
        else:
            biomarker_types = QUANTITATIVE_BIOMARKERS
            biomarker_rel_types = QUANTITATIVE_BIOMARKER_RELS
        
        if self.only_soluble:
            oql_query = 'Select Relation WHERE NeighborOf (SELECT Entity WHERE id = (' + ','.join(map(str,disease_ids)) + '))'
            oql_query += " AND objectType = (" + ','.join(biomarker_rel_types) + ')'
            oql_query += ' AND NeighborOf (SELECT Entity WHERE objectType = (' + ','.join(biomarker_types) + ') AND "Cell Localization" = Secreted)'
            oql_query += ' AND Tissue = ({blood_tissues})'.format(blood_tissues=','.join(SOLUBLE_BIOMARKERS_TISSUES))
            query_name = 'Find soluble {types} biomarkers'.format(types=','.join(biomarker_types))
        else:
            oql_query = OQL.expand_entity_by_id(disease_ids,biomarker_rel_types,biomarker_types)
            query_name = 'Find {types} biomarkers'.format(types=','.join(biomarker_types))

        self.process_oql(oql_query, query_name)

        if not self.score_genetic_biomarkers:
            select_metabolites = '(SELECT Entity WHERE objectType=SmallMol AND Class=\'Endogenous compound\')'
            oql_query = 'SELECT Relation WHERE NeighborOf (SELECT Entity WHERE id = ({ids})) AND NeighborOf {metabolites} AND objectType = ({biomarkers})'
            query_name = 'Find metabolite biomakers'
            if self.only_soluble:
                oql_query = oql_query + ' AND Tissue = ({blood_tissues})'.format(blood_tissues=','.join(SOLUBLE_BIOMARKERS_TISSUES))
                query_name = 'Find metabolite biomakers in blood tissues'
            oql_query = oql_query.format(ids=','.join(map(str,disease_ids)),metabolites=select_metabolites, biomarkers=','.join(biomarker_rel_types))
            self.process_oql(oql_query, query_name)


        biomarker_ids = [y['Id'][0] for x,y in self.Graph.nodes(data=True) if y['ObjTypeName'][0] not in DISEASES]
        disease_ids =  [y['Id'][0] for x,y in self.Graph.nodes(data=True) if y['ObjTypeName'][0] in DISEASES]
        oql_query = OQL.connect_ids(biomarker_ids,disease_ids,['Regulation'])
        self.process_oql(oql_query, 'Find Regulation between biomarkers and diseases')

        if self.only_soluble:
            self.Graph.count_references()
            self.Graph.filter_references({'Tissue':SOLUBLE_BIOMARKERS_TISSUES},QUANTITATIVE_BIOMARKER_RELS)

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

        if self.score_genetic_biomarkers:
            biomarker_rel_types = QUANTITATIVE_BIOMARKER_RELS + GENETIC_BIOMARKER_RELS
        else:
            biomarker_rel_types = QUANTITATIVE_BIOMARKER_RELS

            biomarker_rel_types.append('Regulation')

        for disease in self.diseases:
            dis_name = disease['Name'][0]
            self.set_how2connect(biomarker_rel_types,[],'')
            linked_entities_count = self.link2concept(dis_name,disease_name2ids[dis_name],relations=self.Graph)
            print('%d biomarkers are linked to %s' % (linked_entities_count,dis_name))

        #adding biomarker type column
        self.name2objs, objid2names = bm.Graph.get_prop2obj_dic('Name', list(bm.RefCountPandas['Name']))
        bm.RefCountPandas['Type'] = bm.RefCountPandas['Name'].apply(lambda x: self.name2objs[x][0]['ObjTypeName'][0])

        biomarker_descr = ' genetic ' if bm.score_genetic_biomarkers else ' soluble '
        self.fout_prefix = bm.Disease['Name'][0]+biomarker_descr
        count_file = self.fout_prefix+" biomarkers counts.tsv"
        ref_file = self.fout_prefix+" biomarker references.tsv"
        self.print_ref_count(count_file,referencesOut=ref_file,sep='\t')
        print ('Semantic counts for each biomarker are in %s' % count_file)
        print('References supporting semantic counts are in %s' % ref_file)


    def biomarker_specificity(self,biomarker_ids:tuple):
        #biomarker specificity can be calculated only for quantitative biomarkers
        oql_query = 'SELECT Entity WHERE Connected by (SELECT Relation WHERE objectType = ({rel_types})) to (SELECT Entity WHERE id = ({id}))'
        oql_query = oql_query.format(rel_types=','.join(QUANTITATIVE_BIOMARKER_RELS), id=','.join(map(str,list(biomarker_ids))))
        disease_count = self.get_result_size(oql_query)
        return disease_count

    def add_specificity(self):
        biomarker_count = len(self.RefCountPandas.index)
        print('Finding biomarker specificity for %d biomarkers' % biomarker_count)
        start_specificity = time.time()
        self.RefCountPandas[SPECIFICITY] = self.RefCountPandas[self.__temp_id_col__].apply(lambda x: self.biomarker_specificity(x))
        print ('Specificity for %d biomarkers was found in %s' %(biomarker_count, self.execution_time(start_specificity)))

 
    def weighted_counts(self, add_specificity=False):
        self.weight_prop = 'ObjTypeName'
        self.weight_dict = {'Biomarker':1.0,'StateChange':0.5,'QuantitativeChange':0.5,'GeneticChange':0.5,'CellExpression':0.5,'Regulation':0.25}
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
        
        count_file = self.fout_prefix+" biomarkers normalized weighted scores.tsv" if calculate_specificity else self.fout_prefix+" biomarkers weighted scores.tsv"
        self.print_ref_count(count_file,sep='\t')
        print ('Weighted semantic counts for each biomarker are in %s' % count_file)

    def biomarker_disease_scores(self):
        biomarker2disease = pd.DataFrame()
        row_counter = 0
        ref_count_columns = [col for col in bm.RefCountPandas.columns if bm._col_name_prefix in col]
        for row in bm.RefCountPandas.index:
            biomarker = bm.RefCountPandas.iloc[row][0]
            for col_name in ref_count_columns:
                disease = col_name[len(bm._col_name_prefix):]
                refcount = bm.RefCountPandas.iloc[row][col_name]
                if refcount > 0.0000001:
                    biomarker2disease.at[row_counter,'Biomarker'] = biomarker
                    biomarker2disease.at[row_counter,'Disease'] = disease
                    biomarker2disease.at[row_counter,'Weighted Refcount'] = float(refcount)
                    row_counter += 1

        biomarker2disease.sort_values(by=['Weighted Refcount'],ascending=False, inplace=True)
        biomarker2disease['Type'] = biomarker2disease['Biomarker'].apply(lambda x: self.name2objs[x][0]['ObjTypeName'][0])
        count_file = self.fout_prefix+' biomarker-disease map.txt'
        biomarker2disease.to_csv(count_file,sep='\t', index=False)
        print ('Semantic counts for each biomarker-disease pair are in %s' % count_file)


if __name__ == "__main__":
    start_time = time.time()
    APIconfig = load_api_config()
    bm = BiomarkersReport(APIconfig)
    bm.flush_dump()
    bm.score_genetic_biomarkers = False
    bm.find_only_soluble()

    calculate_specificity = False if bm.score_genetic_biomarkers else False 
    # change here to True to calculate biomarker specificity
    # there is no scientific rational to calculate biomarker specificity for genetic biomarkers
    bm.find_diseases('Uterine neoplasm')
   
    disease_ids =  [y['Id'][0] for x,y in bm.Graph.nodes(data=True) if y['ObjTypeName'][0] in DISEASES]
    bm.load_graph(disease_ids)
    bm.init_semantic_search()
    bm.semantic_search()
    
    bm.drop_refcount_columns()
    bm.weighted_counts(add_specificity=calculate_specificity)

    bm.biomarker_disease_scores()
    print('Biomarkers for %s was found in %s' % (bm.Disease['Name'][0], bm.execution_time(start_time)))
    

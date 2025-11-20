from .DiseaseTargets import DiseaseTargets
from .DiseaseTargets import ANTAGONIST_TARGETS_WS,AGONIST_TARGETS_WS,REFERENCE_IDENTIFIERS,UNKNOWN_TARGETS_WS
from .ResnetGraph import ResnetGraph,PSObject,OBJECT_TYPE,EFFECT,PROTEIN_TYPES
from .SemanticSearch import RANK,time,PHENOTYPE_WORKSHEET
from .ResnetAPISession import APISession,OQL,NO_REL_PROPERTIES,BIBLIO_PROPERTIES,ALL_CHILDS
from ..pandas.panda_tricks import df,np
from ..FDA.fda_api import FDA
from ..ReaxysAPI.Reaxys_API import drugs2props
from numpy import nan_to_num
import networkx as nx
from .DrugTargetConfidence import DrugTargetConsistency
from ..utils import run_tasks,execution_time,os,DEFAULT_CONFIG_DIR



DRUG2TARGET_REGULATOR_SCORE = 'Regulator score'
PHARMAPENDIUM_ID = 'Marketed Drugs'
DRUG_CLASS = 'Drug class' # agonist/antagonist
# do not propagate througuh Binding relations for drug rank calculation:
DT_RELTYPES = ['DirectRegulation','Regulation','Expression','MolSynthesis','MolTransport']
DIETARY_SUPPLEMENT='Dietary Supplement'
REAXYS_FIELDS = {'IDE.CHA':'Charge',
                'IDE.MW':'Molecular Weight',
                'CALC.LOGP':'Lipophilicity',
                'CALC.TPSA':'TPSA', # Topological Polar Surface Area
                'DE.DE':'pKa'
                }

class Drugs4Targets(DiseaseTargets):
  pass
  def __init__(self,*args,**kwargs):
    """
    input:
      APIconfig - args[0]
      what2retrieve - default BIBLIO_PROPERTIES which is required to load self.set_target_disease_state()
      connect2server - default True
    """
    my_kwargs = {
            #'disease':[],
            'what2retrieve':BIBLIO_PROPERTIES,
            'strict_mode':True,
            'data_dir':'',
            'add_bibliography' : False,
            'target_types' : ['Protein'],
           # 'pathway_folders':[],
            #'pathways': [],
            'consistency_correction4target_rank':False,
            'DTfromDB' : False, # set to True to load drug-targets relations from database 
            # instead of using drug2targets cache. Use it if drug-target relations were edited in Pathway Studio
            #'add_inhibitors4':[], # list of tuples (reltype, [target_names])
            "max_childs" : 11,
            #'drug_groups':[],
            'use_in_children':False,
            'BBBP':False, # add Blood-Brain Barrier Permeability estimation to drugs
            'ontology_file' : os.path.join(os.getcwd(),DEFAULT_CONFIG_DIR,'ResnetAPI/ontology/Drugs MoAs.txt')
            }
    
    my_kwargs.update(kwargs)
    super().__init__(*args,**my_kwargs)
    self.target_uid2rank = dict() # {target_dbid:rank}
    self.direct_target2drugs = PSObject() # used for annotation of ANTAGONIST_TARGETS_WS,AGONIST_TARGETS_WS with drugs
    self.indirect_target2drugs = PSObject() # used for annotation of ANTAGONIST_TARGETS_WS,AGONIST_TARGETS_WS with drugs
    self.add_targets2drugs_ws = True
    
  def load_dt(self,**kwargs):
    self.dt_consist = DrugTargetConsistency(self.APIconfig,**kwargs)
    return self.dt_consist
############################ service functions: read/write, load from database, cache######################################            
  def report_name(self):
      rep_pred = 'suggested' if self.params['strict_mode'] else 'predicted'
      return f'{self._disease2str()} {rep_pred} targets,drugs'


  def report_path(self,extension=''):
      if extension:
          ext = extension if extension.find('.')>-1 else '.'+extension
      else:
          ext = ''
      return self.data_dir+self.report_name()+ext


  def __rank(self,targets:list[PSObject],_4drug:PSObject,with_effect:str,correct_by_consistency:bool)->dict[int,float]:
    '''
    output:
      target_uid2corrected_rank = {target_uid:rank},\n
      if "correct_by_consistency" rank from "self.target_uid2rank" corrected by "self.dt_consist" 
    '''
  # debug ranking for a drug here:
    # if _4drug.urn() in ['urn:agi-smol:jr-446','urn:agi-cas:1808916-28-0']:
    #    print('')
    targets_uids = ResnetGraph.uids(targets)
    target_uid2corrected_rank = {k:v for k,v in self.target_uid2rank.items() if k in targets_uids}
    meanpX = 6.7 # this value was calculated as average pX of all DirectRegulation from Reaxys
    for target in targets:
      target_uid = target.uid()
      dt_rels = self.drugs2targets._psrels4(_4drug.uid(),target_uid)
      if dt_rels: # dt_rels may be empty because not every target from df is linked to drug
        dt_rel = dt_rels[0]
        if dt_rel.isprimarytarget():
          pX = float(dt_rels[0].get_prop('Affinity',0,meanpX))
          pXcorrection = 1.5 + pX/12.0
          target_uid2corrected_rank[target_uid] *= pXcorrection
          
        if correct_by_consistency:
          consistency_correction = self.dt_consist.consistency_correction(_4drug,target,with_effect)
          target_uid2corrected_rank[target_uid] *= consistency_correction

    return target_uid2corrected_rank
  

  def add_rank(self,_2df:df,d2t_graph:ResnetGraph):
    """
Adds columns DRUG2TARGET_REGULATOR_SCORE to drug_df.  Used for SNEA
    """
    drug_objs = d2t_graph.psobjs_with(only_with_values=['SmallMol'])
    drug2rank = {d.name():float(d.get_prop(DRUG2TARGET_REGULATOR_SCORE,if_missing_return=0.0)) for d in drug_objs}

    df_copy = _2df.dfcopy()
    # intitalizing columns with default values if df was not visited:
    if DRUG2TARGET_REGULATOR_SCORE not in list(df_copy.columns):
      df_copy[DRUG2TARGET_REGULATOR_SCORE] = np.nan

    for idx in df_copy.index:
      drug_name = df_copy.at[idx,'Name']
      try:
        rank2add = drug2rank[drug_name]
        new_rank = nan_to_num(df_copy.at[idx,DRUG2TARGET_REGULATOR_SCORE]) + rank2add
        df_copy.at[idx,DRUG2TARGET_REGULATOR_SCORE] = new_rank
      except KeyError:
        continue

    df_copy._name_ = _2df._name_
    return df_copy
  
  
  def subtract_rank(self,my_df:df, antigraph:ResnetGraph):
      drug_objs = antigraph.psobjs_with(only_with_values=['SmallMol'])
      drug2rank = {d.name():float(d.get_prop(DRUG2TARGET_REGULATOR_SCORE,if_missing_return=0.0)) for d in drug_objs}
      df_copy = my_df.dfcopy()

      for idx in df_copy.index:
          drug_name = df_copy.at[idx,'Name']
          try:
              rank2subtract = drug2rank[drug_name]
              my_rank = nan_to_num(df_copy.at[idx,DRUG2TARGET_REGULATOR_SCORE])
              new_rank = my_rank-rank2subtract
              df_copy.at[idx,DRUG2TARGET_REGULATOR_SCORE] = new_rank
          except KeyError:
              continue

      return df_copy


  def __add_targets2(self,drug_df:df,drug2targets:ResnetGraph,prefix=''):
    """
Adds following columns to drug_df:\n
Directly inhibited targets',\n'Indirectly inhibited targets',\n'Directly activated targets',\n'Indirectly activated targets'
    """
    def __make_rows()->dict[str,list[str]]:
      columns2add2drug = dict()
      for  drug in drug2targets._get_nodes():
        if drug.objtype() == 'SmallMol': # node is dict here
          direct_targets, indirect_targets = drug2targets.direct_indirect_targets(drug.uid())

          direct_target_names = list()  #direct_targets are sorted by pX
          for target_id, px, consistency in direct_targets:
            direct_target = drug2targets._get_node(target_id)
            target_name = direct_target.name()
            self.direct_target2drugs.update_with_value(target_name,drug.name())
            target_name += ' ('
            if float(px) > 0.0:
                target_name += f'pX={px};'
            target_name += f'cohesion={consistency})'
            direct_target_names.append(target_name)
          direct_targets_str = ','.join(direct_target_names)

          indirect_target_names = list()
          for target_id, px, consistency in indirect_targets:
            indirect_target = drug2targets._get_node(target_id)
            target_name = indirect_target.name()
            self.indirect_target2drugs.update_with_value(target_name,drug.name())
            target_name += ' ('
            if float(px) > 0.0:
              target_name += f'pX={px};'
            target_name += f'cohesion={consistency})'
            indirect_target_names.append(target_name)
          indirect_targets_str = ','.join(indirect_target_names)

          if direct_targets_str or indirect_targets_str:
            drug_class = str(drug.get_prop(DRUG_CLASS))
            columns2add2drug[drug.name()] = [drug_class,direct_targets_str,indirect_targets_str]
      return columns2add2drug

    # making dataframe
    df_copy = drug_df.dfcopy()
    def __column_name(clas:str,is_direct:bool)->str:
      return f"{'Directly' if is_direct else 'Indirectly'} {'inhibited' if clas == 'antagonist' else 'activated'} targets"
    
    # intitalizing columns with default values if df was not visited:
    for drug_class in ['antagonist', 'agonist']:
      for is_direct in [True, False]:
        col_name = __column_name(drug_class, is_direct)
        if col_name not in df_copy.columns:
          df_copy[col_name] = '' 

    columns2add2drug = __make_rows()
    for idx in df_copy.index:
      drug_name = df_copy.at[idx,'Name']
      if drug_name in columns2add2drug:
        drug_class,direct_targets,indirect_targets = columns2add2drug[drug_name]
        if direct_targets:
          colname = __column_name(drug_class,True)
          if df_copy.at[idx,colname] == '':
            df_copy.loc[idx,colname] = (prefix+direct_targets).strip()
          else:
            df_copy.loc[idx,colname] = ','.join([df_copy.loc[idx,colname],prefix+direct_targets])
        
        if indirect_targets:
          colname = __column_name(drug_class,False)
          if df_copy.loc[idx,colname] == '':
            df_copy.loc[idx,colname] = (prefix+indirect_targets).strip()
          else:
            df_copy.loc[idx,colname] = ','.join([df_copy.loc[idx,colname],prefix+indirect_targets])

    df_copy._name_ = drug_df._name_
    return df_copy


  def __rank_drugs4(self,targets:list[PSObject],with_effect:str,correct_by_consistency:bool):
    """
    output:
      ResnetGraph with drugs annotated with DRUG2TARGET_REGULATOR_SCORE, DRUG_CLASS (agonist/antagonist)
    """
    my_dtG = self.drugs2targets.neighborhood(set(targets),[],DT_RELTYPES,[with_effect])
    my_drugs = my_dtG._psobjs_with('SmallMol','ObjTypeName')
    my_drugs = [d for d in my_drugs if d not in self._targets()] # to remove metabolite targets
    
    # initializing drug ranks
    drug2rank = {d.uid():[my_dtG.rank_regulator(d,self.__rank(targets,d,with_effect,correct_by_consistency))] for d in my_drugs}
    nx.set_node_attributes(my_dtG,drug2rank,DRUG2TARGET_REGULATOR_SCORE)

    drug_class = 'agonist' if with_effect == 'positive' else 'antagonist'
    druguid2class = {d.uid():[drug_class] for d in my_drugs}
    nx.set_node_attributes(my_dtG,druguid2class,DRUG_CLASS)
    
    return my_dtG


  def load_target_ranks(self,from_ranked_targets_df:df,for_antagonists=True):
    '''
    input:
      from_ranked_targets_df must have RANK column
    output:
      self.targets4agonists\n
      self.targets4antagonists\n
      self.target_uid2rank
    '''
    my_targets = self._targets()
    name2target = {n.name():n for n in my_targets}
    target_name_column = self.__resnet_name__ if self.__resnet_name__ in from_ranked_targets_df.columns else 'Name'
    
    df_targets = set()
    missing_targets_counter = 0
    for idx in from_ranked_targets_df.index:
      target_name = from_ranked_targets_df.at[idx,target_name_column]
      target_rank = from_ranked_targets_df.at[idx,RANK]
      if target_name in name2target:
        target = name2target[target_name]
        df_targets.add(target)
        target_uid = target.uid()
        self.target_uid2rank[target_uid] = float(target_rank)
        if target_uid in self.Graph: # self.Graph is empty in df for SNEA
          # boosting drug ranking by adding target children to ranking network with the rank of the parent
          children = self.Graph._get_node(target_uid).childs() # only self.Graph is guaranteed to have CHILDS annotation
          for child in children:
            if child.name() not in name2target:
              self.target_uid2rank[child.uid()] = float(target_rank) 
              df_targets.add(child)
      else:
        print(f'Node named "{target_name}" is not found')
        missing_targets_counter += 1

    if missing_targets_counter:
      print(f'{missing_targets_counter} targets out of {len(my_targets)} were not found in self.targets()')

    if for_antagonists:
        self.targets4antagonists = list(df_targets)
    else:
        self.targets4agonists = list(df_targets)
    return


  def init_drug_df(self, drugs:list[PSObject]):
      '''
      input:
        "drugs" used to load drugs from SNEA samples
      loads:
        self.RefCountPandas
      '''
      forbidden_drugs =['DMSO', 'glucose','nicotine']
      print('Loading drug ranking worksheet')
      drug2pharmapendium_id = dict()
      clean_drugs = list(drugs)
      # bad fix for PS ontology problem:
      for drug in drugs:
        drug_name = drug.name()
        if 'antigen' in drug_name:
          clean_drugs.remove(drug)
        if 'vaccine' in drug_name: 
          clean_drugs.remove(drug)
        if drug_name in forbidden_drugs: 
          clean_drugs.remove(drug)
        
        pp_id = str(drug.get_prop('PharmaPendium ID',if_missing_return='Dietary Supplement'))
        drug2pharmapendium_id[drug_name] = pp_id

      count = len(clean_drugs)
      clean_drugs = [d for d in clean_drugs if d not in self.disease_inducers]
      clean_drugs = [d for d in clean_drugs if "'" not in d.urn()] # drugs with "'" in URN are not searchable in PS

      new_session = self._clone_session(what2retrieve=NO_REL_PROPERTIES)
      PAINScompounds = set(new_session.process_oql(OQL.selectPAINs(),'Select PAINS compounds')._get_nodes())
      clean_drugs = [d for d in clean_drugs if d not in PAINScompounds]
      new_session.close_connection()
      
      removed_count = count-len(clean_drugs)
      print(f'{removed_count} drugs were removed because they are known to induce {self._disease2str()}')

      mapped_drugs,_ = self.load_dbids4(clean_drugs,with_props=['Reaxys ID'])
      drug_df = self.load_df(clean_drugs,max_childs=0,max_threads=25)
      self.RefCountPandas = drug_df.merge_dict(drug2pharmapendium_id,new_col=PHARMAPENDIUM_ID,map2column='Name')

      print('Initialized "Drugs" worksheet with %d drugs from database for ranking' % len(self.RefCountPandas))
      # self.RefCountPandas does not have self.__temp_id_col__ column at this point 
      # because ALL_CHILDS = 0. It will be added in link2disease_concepts()
      return mapped_drugs


  def regulatory_rank(self):
      '''
      output:
        df with new columns: DRUG2TARGET_REGULATOR_SCORE, 'Directly inhibited targets','Directly activated targets',
        'Indirectly inhibited targets', 'Indirectly activated targets'
      '''
      need_correction = self.params.get('consistency_correction4target_rank',True)
      antagonist_graph = self.__rank_drugs4(self.targets4antagonists,with_effect='negative',correct_by_consistency=need_correction)
      agonist_graph = self.__rank_drugs4(self.targets4agonists,with_effect='positive',correct_by_consistency=need_correction) 

      antagonist_antigraph = self.__rank_drugs4(self.targets4antagonists,with_effect='positive',correct_by_consistency=False)
      agonist_antigraph = self.__rank_drugs4(self.targets4agonists,with_effect='negative',correct_by_consistency=False)

      new_ranked_df = self.add_rank(self.RefCountPandas,antagonist_graph)
      new_ranked_df = self.add_rank(new_ranked_df,agonist_graph)
      new_ranked_df = self.subtract_rank(new_ranked_df,antagonist_antigraph)
      new_ranked_df = self.subtract_rank(new_ranked_df,agonist_antigraph)

      if self.add_targets2drugs_ws:
        new_ranked_df = self.__add_targets2(new_ranked_df,antagonist_graph)
        new_ranked_df = self.__add_targets2(new_ranked_df,agonist_graph)
        new_ranked_df = self.__add_targets2(new_ranked_df,antagonist_antigraph,prefix='\ndisease-promoting: ')
        new_ranked_df = self.__add_targets2(new_ranked_df,agonist_antigraph,prefix='\ndisease-promoting :')

      predicted_drugs = new_ranked_df.greater_than(0,DRUG2TARGET_REGULATOR_SCORE)['Name'].to_list()
      known_drugs = ResnetGraph.names(self.drugs_linked2disease)
      drugs2score = known_drugs + predicted_drugs
      new_ranked_df = new_ranked_df.filter_by(drugs2score,'Name')
      new_ranked_df[DRUG2TARGET_REGULATOR_SCORE] = new_ranked_df[DRUG2TARGET_REGULATOR_SCORE].fillna(0.0)

      disnames = ','.join(self.input_disease_names())
      new_ranked_df._name_ = f'Drugs for {disnames}'
      new_ranked_df = new_ranked_df.sortrows(by=[DRUG2TARGET_REGULATOR_SCORE])
      new_ranked_df.col2rank[DRUG2TARGET_REGULATOR_SCORE] = 1
      print('Found %d drugs for %d antagonist targets and %d agonist targets for "%s"' % 
          (len(new_ranked_df),len(self.targets4antagonists),len(self.targets4agonists),self.input_disease_names()),flush=True)
      return new_ranked_df
      

  def link2disease_concepts(self,in_drugdf:df):
    '''
    adds self.__temp_id_col__ column to in_drugdf and scores concepts for drugs
    '''
    if not self.input_diseases: return in_drugdf # for SNEA.make_drugs_df()
    if len(in_drugdf) < 2:
      print('Drugs worksheet is empty !!!!')
      return in_drugdf

    if self.__temp_id_col__ not in in_drugdf.columns: # case when drug_df was loaded from cache RNEF file:
    # max_childs_count must be zero (=ALL_CHILDS) to avoid removal of drugs with children 
    # and to merge with results of SemanticSearch.bibliography() future
      drug_df = self.add_temp_id(in_drugdf,max_childs=ALL_CHILDS,max_threads=25)
      before_mapping = len(drug_df)
      print(f'{before_mapping - len(drug_df)} rows were deleted because database identifier cannot be found for drug names')
    else:
      drug_df = in_drugdf.dfcopy()

    # add functions here for speed debuging.  drug_df has self.__temp_id_col__ column here:
    phenotypedf_rows = list()
    def concepts2rows(concepts:list[PSObject],concept_name:str):
      '''
      function creates info worksheet for concepts called "Phenotype"
      '''
      concepts_rows = []
      for concept in concepts:
        rank = concept.get_prop('rank')
        weight = concept.get_prop('regulator weight',if_missing_return=1.0)
        connectivity = concept.get_prop('Local connectivity')
        concepts_rows.append([concept_name,rank,concept.name(),len(concept.childs()),weight,connectivity])
      return concepts_rows
        
    kwargs = {'connect_by_rels':['Regulation'],
              'with_effects':['positive'],
              'boost_with_reltypes':['FunctionalAssociation','Regulation'],
              'clone2retrieve' : REFERENCE_IDENTIFIERS,
              'init_refstat' : False,
              'column_name': 'Cell processess to activate in '+ self._disease2str(),
              'add_relevance_concept_column':True}
    drug_df,concepts = self.score_concept('processes2activate',drug_df,**kwargs)[2:4]
    if concepts:
      phenotypedf_rows += concepts2rows(concepts,kwargs['column_name'])

    kwargs['with_effects'] = ['negative']
    kwargs['column_name'] = 'Cell processess to inhibit in '+ self._disease2str()
    kwargs['column_rank'] = drug_df.max_colrank()
    drug_df,concepts = self.score_concept('processes2inhibit',drug_df,**kwargs)[2:4]
    if concepts:
      phenotypedf_rows += concepts2rows(concepts,kwargs['column_name'])
    kwargs.pop('column_rank')
    
    # self.params['processes'] are formed by score_target_semantics()
    kwargs['with_effects'] = ['unknown']
    kwargs['column_name'] = 'Cell processess affected by '+ self._disease2str()
    drug_df,concepts = self.score_concept('processes',drug_df,**kwargs)[2:4]
    if concepts:
      phenotypedf_rows += concepts2rows(concepts,kwargs['column_name'])

    kwargs['with_effects'] = ['positive']
    kwargs['column_name'] = 'Cells to activate  in '+ self._disease2str()
    drug_df,concepts = self.score_concept('cells2activate',drug_df,**kwargs)[2:4]
    if concepts:
      phenotypedf_rows += concepts2rows(concepts,kwargs['column_name'])

    kwargs.pop('with_effects','')
    kwargs['column_name'] = 'Clinical parameters for '+ self._disease2str()
    drug_df,concepts = self.score_concept('clinical_parameters',drug_df,**kwargs)[2:4]
    if concepts:
      phenotypedf_rows += concepts2rows(concepts,kwargs['column_name'])

    if 'symptoms' in self.params:
      # calculating first symptoms aggravated by the drug. All columns except 'Aggravated symptoms" will be deleted
      # column  Refcount aggravate symptoms for ... will be subtaracted from symptoms columns 
      # and drugs known to exacerbate symptoms more than inhibit them will be deleted from drug_df
      kwargs['with_effects'] = ['positive']
      aggravate_symptoms_col = 'aggravate symptoms for '+ self._disease2str()
      kwargs['column_name'] = aggravate_symptoms_col
      kwargs['column_rank'] = drug_df.max_colrank()
      drug_df = self.score_concept('symptoms',drug_df,**kwargs)[2]
      drug_df = drug_df.dfcopy(rename2={'Relevant symptoms':'Aggravated symptoms'})

      symptoms_columns = []
      colname = 'inhibit symptoms for '+ self._disease2str()
      kwargs['column_name'] = colname
      kwargs['with_effects'] = ['negative']
      #kwargs['column_rank'] = 1
      drug_df,concepts = self.score_concept('symptoms',drug_df,**kwargs)[2:4]
      if concepts:
        phenotypedf_rows += concepts2rows(concepts,kwargs['column_name'])
        symptoms_columns.append(colname)

      kwargs['with_effects'] = ['unknown']
      colname = 'symptoms for '+ self._disease2str()
      kwargs['column_name'] = colname
      #kwargs['column_rank'] = 2
      drug_df,concepts = self.score_concept('symptoms',drug_df,**kwargs)[2:4]
      if concepts:
        phenotypedf_rows += concepts2rows(concepts,kwargs['column_name'])
        symptoms_columns.append(colname)

      aggravated_symptoms_refcount_col = self._refcount_colname(aggravate_symptoms_col)
      for col in symptoms_columns:
        refcount_symptom_col = self._refcount_colname(col)
        drug_df[refcount_symptom_col] -= drug_df[aggravated_symptoms_refcount_col]
        # removing drugs that aggravate symptoms more than inhibit them:
        drug_df = drug_df.greater_than(-1,refcount_symptom_col)
      
      #dropping temporary aggravated_symptoms columns:
      aggravated_symptoms_cols = self._refcount_columns(aggravate_symptoms_col)
      my_cols = drug_df.columns.drop(aggravated_symptoms_cols).to_list()
      drug_df = drug_df.dfcopy(my_cols)
      kwargs.pop('with_effects','')

    kwargs['column_name'] = 'regulation of diseases similar to '+ self._disease2str()
    drug_df,concepts = self.score_concept('similar_diseases',drug_df,**kwargs)[2:4]
    if concepts:
      phenotypedf_rows += concepts2rows(concepts,kwargs['column_name'])

    kwargs['add_relevance_concept_column'] = False
    kwargs['column_name'] = 'regulation of '+ self._disease2str()
    drug_df = self.score_concepts(self.input_diseases,drug_df,**kwargs)[2]

    kwargs['column_name'] = self._disease2str()+' clinical trials'
    kwargs['boost_with_reltypes'] = []
    drug_df = self.score_concepts(self.input_diseases,drug_df,**kwargs)[2]

    rank_col = 'Concept rank'
    children_col = '# children'
    concepts_df = df.from_rows(phenotypedf_rows,['Type',rank_col,'Name',children_col,'weight','# linked drugs (includes links to concept ontology children)'])
    concepts_df = concepts_df.sortrows(by=[rank_col,'Type','weight',children_col],ascending=[True,True,False,False])
    concepts_df._name_ = PHENOTYPE_WORKSHEET

    if self.params.get('BBBP',False):
      drug_df = self.__addBBBP(drug_df)

    self.add2report(concepts_df)
    return drug_df


  def make_raw_drug_df(self):
    self.load_target_ranks(self.report_pandas[AGONIST_TARGETS_WS],for_antagonists=False)
    self.load_target_ranks(self.report_pandas[ANTAGONIST_TARGETS_WS],for_antagonists=True) 
    ranked_df = self.regulatory_rank()
    # add here functions for speed debuging.  ranked_df does not have self.__temp_id_col__ column here
    if self.params.get('disease',''): 
    # case of drug repurposing using input disease network and targets from knowledge graph
      if self.params.get('debug',False):
        if self.params.get("use_in_children",False):
          d2cd = self.drug2childdose(ranked_df,multithread=False)
          ranked_df = ranked_df.merge_dict(d2cd,'Children dose',PHARMAPENDIUM_ID)

        if self.params.get('add_bibliography',False):
          # has to follow the same order as in "no debug" mode to have the same results for testing
          names2reflinks = self.RefStats.reflinks(ranked_df,'Name',self.names4tmsearch())
          full_drug_df = self.link2disease_concepts(ranked_df)
          refcount_col = self.tm_refcount_colname('Name',self.names4tmsearch())
          full_drug_df = self.RefStats.add_reflinks(names2reflinks,refcount_col,full_drug_df,'Name')
        else:
          full_drug_df = self.link2disease_concepts(ranked_df)

      elif self.params.get('add_bibliography',False):
        tasks = [(self.RefStats.reflinks,(ranked_df,'Name',self.names4tmsearch(),[],True))]
        tasks.append((self.link2disease_concepts,(ranked_df,)))
        if self.params.get("use_in_children",False):
          tasks.append((self.drug2childdose,(ranked_df,)))
        results = run_tasks(tasks)

        full_drug_df = results['link2disease_concepts']
        names2reflinks = results['reflinks']
        refcount_col = self.tm_refcount_colname('Name',self.names4tmsearch())
        full_drug_df = self.RefStats.add_reflinks(names2reflinks,refcount_col,full_drug_df,'Name')

        if self.params.get("use_in_children",False):
          d2cd = results['drug2childdose']
          full_drug_df = full_drug_df.merge_dict(d2cd,'Children dose',PHARMAPENDIUM_ID)
      else:
        tasks = [(self.link2disease_concepts,(ranked_df,))]
        if self.params.get("use_in_children",False):
          tasks.append((self.drug2childdose,(ranked_df,)))
        results = run_tasks(tasks)

        full_drug_df = results['link2disease_concepts']
        if self.params.get("use_in_children",False):
          d2cd = results['drug2childdose']
          full_drug_df = full_drug_df.merge_dict(d2cd,'Children dose',PHARMAPENDIUM_ID)
    else:
        full_drug_df = ranked_df # case of finding drugs from SNEA results

    raw_drug_df = self.make_count_df(full_drug_df,'rawDrugs')
    self.add2raw(raw_drug_df)
    return raw_drug_df


  def score_drugs(self,normalize=True):
    '''
    output:
      raw_drug_df ('rawDrugs'),ranked_drugs_df ('Drugs')
      adds  raw_drug_df to self.raw_data['rawDrugs'], adds ranked_drugs_df to self.report_pandas['Drugs]
      set SNEA "normalize" to False 
    '''
    start_time = time.time()
    raw_drug_df = self.make_raw_drug_df()
    # add here functions for debuging. raw_drug_df does not have self.__temp_id_col__ column here:

    if normalize:
      ranked_drugs_df,_ = self.normalize('rawDrugs','Drugs','Name')

      #reordering drug df columns:
      drugdf_columns = list(ranked_drugs_df.columns)
      drugdf_columns.remove(PHARMAPENDIUM_ID)
      drugdf_columns = [PHARMAPENDIUM_ID]+drugdf_columns
      ranked_drugs_df = ranked_drugs_df.dfcopy(drugdf_columns)
      
      ranked_drugs_df.tab_format['tab_color'] = 'green'
      ranked_drugs_df._name_ = 'Drugs'
      self.add2report(ranked_drugs_df)
      
      print("Drug ranking was done in %s" % execution_time(start_time), flush=True)
      print('Normalized worksheet named "Drugs" was added to report')
      return raw_drug_df,ranked_drugs_df
    else:
      print('Ranked drugs are in worksheet "rawDrugs" in raw data')
      return raw_drug_df
    

  def drugs4metabolites(self)->ResnetGraph:
    met_targets = self.params.get('metabolites2inhibit',[])
    if met_targets:
      oql = f'SELECT Relation WHERE Effect = negative AND NeighborOf upstream (SELECT Entity WHERE Name = ({met_targets}))\
          AND NeighborOf downstream ({OQL.select_drugs()})'
      drugs4metG = self.process_oql(oql,'Find drugs for metabolites')
      drugs = drugs4metG.regulators()
      print(f'Found {len(drugs)} drugs inhibiting {met_targets}')
      #self.Graph = self.Graph.compose(drugs4metG)
      return drugs4metG
    else:
      return ResnetGraph()
    

  def drugs4add_inhibitors4(self)->tuple[ResnetGraph,ResnetGraph]:
    reltype2targets = dict(self.params.get("add_inhibitors4",dict()))
    if not reltype2targets: return ResnetGraph(), ResnetGraph()
    targetnames4antagonists = self.report_pandas[ANTAGONIST_TARGETS_WS]['Name'].to_list()
    targetnames4agonists = self.report_pandas[AGONIST_TARGETS_WS]['Name'].to_list()
    bad_drugsG = ResnetGraph()
    good_drugsG = ResnetGraph()
    for reltype,target_names in reltype2targets.items():
      targets2expand = [n for n in self.__targets__ if n.name() in target_names]
      targets2expand_dbids = ResnetGraph.dbids(targets2expand)
      oql_query = f'SELECT Relation WHERE NeighborOf upstream (SELECT Entity WHERE id = ({targets2expand_dbids}))'
      oql_query += f' AND NeighborOf downstream ({OQL.select_drugs()}) AND objectType = {reltype}'
      expanded_targets_g = self.process_oql(oql_query,f'Find drugs modulating {target_names} by {reltype}')
      bad_drugsG =  bad_drugsG.compose(expanded_targets_g.subgraph_by_relprops(['positive'],[EFFECT]))
      good_rels = (expanded_targets_g.subgraph_by_relprops(['negative'],[EFFECT]))._psrels()
      # assigning synergistic effect to good drugs relations ensures use of target rank 
      # by self.__rank_drugs4 to score good drug.
      for rel in good_rels:
        t_name = rel.targets()[0].name()
        if t_name in targetnames4antagonists:
          rel[EFFECT] = ['negative']
        elif t_name in targetnames4agonists:
          rel[EFFECT] = ['positive']
        else: continue
        rel.urn(refresh=True)
      good_drugsG = good_drugsG.compose(ResnetGraph.from_rels(good_rels))

    return good_drugsG, bad_drugsG # bad drugs must de deleted from ranking 


  def drugs_induce_symptoms(self,minref=4):
    symptom_names = self.params.get('symptoms',[])
    if isinstance(symptom_names,dict):
      symptom_names = list(symptom_names.keys())

    processes2inhibit = self.params.get('processes2inhibit',[])
    if isinstance(processes2inhibit,dict):
      processes2inhibit = list(processes2inhibit.keys())

    symptom_names += processes2inhibit
    symptom_names += self.input_disease_names()

    bad_drugs = []
    if symptom_names:
      #select_symptoms = OQL.get_childs(symptom_names,['Name'],include_parents=True)
      select_symptoms = OQL.get_entities_by_props(symptom_names,['Name'])
      oql = f'SELECT Relation WHERE NeighborOf downstream ({OQL.select_drugs()}) \
        AND NeighborOf upstream ({select_symptoms}) AND Effect = positive'
      if minref:
        oql += f' AND RelationNumberOfReferences >= {minref}'
      bad_drugs += self.process_oql(oql)._psobjs_with('SmallMol',OBJECT_TYPE)

    processes2activate = self.params.get('processes2activate',[])
    if isinstance(processes2activate,dict):
      processes2activate = list(processes2activate.keys())

    if processes2activate:
      #select_processes = OQL.get_childs(processes2activate,['Name'],include_parents=True)
      select_processes = OQL.get_entities_by_props(processes2activate,['Name'])
      oql = f'SELECT Relation WHERE NeighborOf downstream ({OQL.select_drugs()}) \
        AND NeighborOf upstream ({select_processes}) AND Effect = negative'
      if minref:
        oql += f' AND RelationNumberOfReferences >= {minref}'
      bad_drugs += self.process_oql(oql)._psobjs_with('SmallMol',OBJECT_TYPE)

    print(f'Found {len(bad_drugs)} drugs inducing symptoms.  They will be deleted from the drug list')
    return set(bad_drugs)
    

  def select_drugs(self,limit2drugs:set[PSObject]=set())->list[PSObject]:
      '''
      input:
        limit2drugs - {PSObject}, used by SNEA make_drugs_df
      loads:
        self.drugs2targets - ResnetGraph with drugs linked to self._targets()
      output:
        list of drugs from self.drugs2targets
      '''
      DTfromDB = self.params.get('DTfromDB',False)
      self.drugs2targets = self.dt_consist.load_drug_graph(self._targets(),limit2drugs,DTfromDB).copy()
      self.drugs2targets = self.drugs2targets.compose(self.drugs4metabolites())
      # subtracting self._targets( to remove metabolite targets
      drugs_linked2targets = set(self.drugs2targets._psobjs_with('SmallMol','ObjTypeName')) - self._targets()
      print(f'Found {len(drugs_linked2targets)} drugs linked to {len(self._targets())} targets for ranking')

      drugs_withno_targets = []
      if hasattr(self,'drugs_linked2disease'):
        # adding drugs that are known to inhibit disease but have no links to disease targets
        drugs_withno_targets = list(set(self.drugs_linked2disease).difference(drugs_linked2targets))
        if drugs_withno_targets:
          my_api_session = APISession(self.APIconfig,NO_REL_PROPERTIES)
          my_api_session.entProps = ['Name','PharmaPendium ID']
          oql_query = 'SELECT Entity WHERE URN = ({props})'
          oql_query += " AND InOntology (SELECT Annotation WHERE Ontology='Pathway Studio Ontology' AND Relationship='is-a') under (SELECT OntologicalNode WHERE Name = drugs)"
          r_n = f'Find drugs linked to disease that have no targets linked to disease'
          drugs_withno_targets_urns = set(ResnetGraph.urns(drugs_withno_targets))
          drugs4disease = self.iterate_oql(oql_query,drugs_withno_targets_urns,request_name=r_n)
          drugs_withno_targets = drugs4disease._get_nodes()
          
          if drugs_withno_targets:
            print('Found additional %d drugs linked to "%s" but not linked to its targets' % 
                (len(drugs_withno_targets),self._disease2str()))
            
      good_drugsG, bad_drugsG = self.drugs4add_inhibitors4()
      if bad_drugsG:
        bad_drugs = set(bad_drugsG._psobjs_with('SmallMol','ObjTypeName'))
        self.drugs2targets.remove_nodes_from(ResnetGraph.uids(bad_drugs))
        print(f'Excluded {len(bad_drugs)} drugs that activate targets that must be inhibited')
      if good_drugsG:
        self.drugs2targets = self.drugs2targets.compose(good_drugsG)
        print(f'Add {len(good_drugsG)} drugs that inhibit targets that must be inhibited')

      selected_drugs = drugs_linked2targets - set(drugs_withno_targets)
     # selected_drugs = selected_drugs - self.drugs_induce_symptoms()
      return list(selected_drugs)
  

  def drug2childdose(self,rankedf:df,multithread = True):
    print('Adding use in children from FDA drug labels')
    marketed_drugs = rankedf[rankedf[PHARMAPENDIUM_ID].notna()][PHARMAPENDIUM_ID].to_list()
    marketed_drugs = [d for d in marketed_drugs if d and d != DIETARY_SUPPLEMENT]
    if marketed_drugs:
      fda_api = FDA()
      return fda_api.child_doses_mt(marketed_drugs) if multithread else fda_api.child_doses(marketed_drugs)
    else:
      return dict()
    

  def init_load_score(self):
    my_drugs = self.select_drugs()
    if self.params['debug']:
      DBdrugs = self.init_drug_df(my_drugs) # DBdrugs have "Molecular weight","Reaxys ID" properties
      if self.params.get('BBBP',False):
        f2d2prop = drugs2props(list(DBdrugs),REAXYS_FIELDS)
        for f,mapdict in f2d2prop.items():
          self.RefCountPandas[f] = self.RefCountPandas['Name'].map(mapdict)
      if self.params.get('consistency_correction4target_rank',True):
        load_fromdb = self.params.get('recalculate_dtconsistency',False)
        self.dt_consist.load_confidence_dict(self._targets(),{},load_fromdb)
        self.dt_consist.save_network()
      self.score_drugs()
      self.dt_consist.clear()
    elif self.params.get('consistency_correction4target_rank',True):
      load_fromdb = self.params.get('recalculate_dtconsistency',False)
      tasks = [(self.init_drug_df,(my_drugs,))] # need comma after my_drugs to make it a tuple
      tasks.append((self.dt_consist.load_confidence_dict,(self._targets(),{},load_fromdb)))
      DBdrugs = run_tasks(tasks)['init_drug_df']

      tasks = [(self.score_drugs,()),(self.dt_consist.save_network,())]
      if self.params.get('BBBP',False):
        tasks.append((drugs2props,(list(DBdrugs),REAXYS_FIELDS)))
        results = run_tasks(tasks)
        f2d2prop = results['drugs2props']
        for f,mapdict in f2d2prop.items():
          self.report_pandas['Drugs'][f] = self.report_pandas['Drugs']['Name'].map(mapdict)
      else:
        run_tasks(tasks)

      self.dt_consist.clear()
    else:
      self.init_drug_df(my_drugs)
      self.score_drugs()
    return


  def find_rank_targets(self):
    super().make_report()
    print('Found %d targets for antagonists, %d targets for agonists, %d targets with unknown disease state'%
    (len(self.report_pandas[ANTAGONIST_TARGETS_WS])-1, len(self.report_pandas[AGONIST_TARGETS_WS])-1,len(self.report_pandas[UNKNOWN_TARGETS_WS])-1))
    print('%d targets in total were ranked' % len(self._targets()))
    return self._targets()


  def __addBBBP(self,_2df:df):
    print('Adding "Brain-Plasma ratio" to report')
    bpG = ResnetGraph.fromRNEF('D:/Python/BBB/DrugsBrainPlasmaRatio.rnef')
    nodes_with_bp = bpG._get_nodes()
    print(f'Loaded {len(nodes_with_bp)} drugs with known brain-plasma ratio')
    urn2bp = dict()
    for n in nodes_with_bp:
      bp = n.get_prop('Average blood-plasma ratios')
      if bp:
        urn2bp[n.urn()] = float(bp)
        
    drug_df = _2df.dfcopy()
    drug_df['Brain-Plasma ratio (%)'] = drug_df['URN'].map(urn2bp)
    not_nan_count = drug_df['Brain-Plasma ratio (%)'].count()
    print(f"{not_nan_count} drugs were annotated with 'Brain-Plasma ratio (%)'")
    return drug_df


  def annotate_report(self):
    '''
    self.report_pandas['Drugs'] does not have self.__temp_id_col__ column at this point
    '''
    self.report_pandas[ANTAGONIST_TARGETS_WS] = self.report_pandas[ANTAGONIST_TARGETS_WS].merge_psobject(self.direct_target2drugs,'Directly Inhibited by','Name',values21cell=True)
    self.report_pandas[ANTAGONIST_TARGETS_WS] = self.report_pandas[ANTAGONIST_TARGETS_WS].merge_psobject(self.indirect_target2drugs,'Indirectly Inhibited by','Name',values21cell=True)
    
    self.report_pandas[AGONIST_TARGETS_WS] =  self.report_pandas[AGONIST_TARGETS_WS].merge_psobject(self.direct_target2drugs,'Directly Activated by','Name',values21cell=True)
    self.report_pandas[AGONIST_TARGETS_WS] = self.report_pandas[AGONIST_TARGETS_WS].merge_psobject(self.indirect_target2drugs,'Indirectly Activated by','Name',values21cell=True)

    drug_groups = self.params.get('drug_groups',[])
    self.add_groups(self.report_pandas['Drugs'],drug_groups)

    #best_drugs = self.report_pandas['Drugs'].smaller_than(0.05,RANK + ' expopvalue')
    self.add2report(self.report_pandas['Drugs'].column_stats('Groups',sep=','))
    return


  def refs2targets(self):
    if self.params['add_bibliography']:
      search_names = self.names4tmsearch()
      self.refs2report(ANTAGONIST_TARGETS_WS,search_names)
      self.refs2report(AGONIST_TARGETS_WS,search_names)
      self.refs2report(UNKNOWN_TARGETS_WS,search_names)


  def make_report(self):
    start_time = time.time()
    if self.params['debug']:
      self.find_rank_targets()
      self.refs2targets()
      self.load_dt()
      self.init_load_score()
      self.add_tm_bibliography_df()
    else:
      tasks = [(self.load_dt,())]
      tasks.append((self.find_rank_targets,()))
      run_tasks(tasks)

      if self.params['add_bibliography']:
        tasks = [(self.init_load_score,()),(self.refs2targets,())]
        run_tasks(tasks)
        self.add_tm_bibliography_df()
      else:
        self.init_load_score()
    
    self.annotate_report()
    self.add_params2df()
    print('Drug repurposing for %s was done in %s' % (self._disease2str(), execution_time(start_time)))

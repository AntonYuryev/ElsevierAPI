from ElsevierAPI import load_api_config, APISession
from ETM_API.references import PS_ID_TYPES
from  ElsevierAPI.ResnetAPI.ResnetGraph import PSObject,CHILDS
from  ElsevierAPI.pandas.panda_tricks import ExcelWriter,df

from TargetGeneticsReport import DATA_DIR
DATA_DIR = 'D:/Python/PMI/'

class TissueLocalization(APISession):
    pass
    def __init__(self, APIconfig):
        super().__init__(APIconfig['ResnetURL'], APIconfig['PSuserName'], APIconfig['PSpassword'])
        self.printed_ids = set()
        self.report = df(columns=['Top category','#References','Organs','Tissues','Cell types','Cell lines'])
        self.flush_dump_files()
        self.PageSize = 1000

    def get_orphan_ids(self):
        return set(self.Graph).difference(self.printed_ids)


    def get_columns(self, node:PSObject):
        sub_organs = set()
        sub_tissues = set()
        celltypes = set()
        cellines = set()


        for id in node[CHILDS]:
            if id in self.printed_ids: continue
            references = self.Graph.get_neighbors_refs([id])
            child = self.Graph._get_node(id)

            if 'Organ' in child['ObjTypeName']:
                sub_organs.add((child,len(references)))
            elif 'Tissue' in child['ObjTypeName']:
                sub_tissues.add((child,len(references)))
            elif 'SemanticConcept' in child['ObjTypeName']:
                cellines.add((child,len(references)))
            elif 'CellType' in child['ObjTypeName']:
                if str(child['URN'][0]).find('cellline') > 0:
                    cellines.add((child,len(references)))
                else:
                    celltypes.add((child,len(references)))
            else:
                print('%s has unknown anatomy type %s' % (child['Name'][0],child['ObjTypeName'][0]))

        sub_organs = list(sub_organs)
        sub_tissues = list(sub_tissues)
        celltypes = list(celltypes)
        cellines = list(cellines)

        sub_organs.sort(key=lambda x: x[1], reverse = True)
        sub_tissues.sort(key=lambda x: x[1], reverse = True)
        celltypes.sort(key=lambda x: x[1], reverse = True)
        cellines.sort(key=lambda x: x[1], reverse = True)

        sub_organ_str = ''
        if sub_organs:
            sub_organ_str = str(len(sub_organs))+' organs:'+','.join([s[0]['Name'][0]+'('+str(s[1])+')' for s in sub_organs])

        tissue_str=''
        if sub_tissues:
            tissue_str = str(len(sub_tissues))+' tissues:'+','.join([s[0]['Name'][0]+'('+str(s[1])+')' for s in sub_tissues])

        cell_str=''
        if celltypes:
            cell_str = str(len(celltypes))+' cell types:'+','.join([s[0]['Name'][0]+'('+str(s[1])+')' for s in celltypes])

        cell_line_str = ''
        if cellines:
            cell_line_str = str(len(cellines))+' cell lines:'+','.join([s[0]['Name'][0]+'('+str(s[1])+')' for s in cellines])
        
        
        return sub_organ_str, tissue_str, cell_str, cell_line_str


    def add2report(self,anatomy_obj:PSObject):
        if anatomy_obj['Id'][0] in self.printed_ids: return
        sub_organ_str, tissue_str, cell_str, cell_line_str = self.get_columns(anatomy_obj)

        obj_name = anatomy_obj['Name'][0]
        anatomy_branch_ids = anatomy_obj[CHILDS]+anatomy_obj['Id']

        relation_graph = self.Graph.get_neighbors_graph(set(anatomy_branch_ids))
        references = relation_graph.load_references()
        ref_count = len(references)
        self.report.loc[len(self.report)] = [obj_name,ref_count,sub_organ_str,tissue_str,cell_str,cell_line_str]
        self.printed_ids.update(anatomy_branch_ids)


    def __get_parents(self, for_child_ids:list, depth:int):
        parents = self.add_parents(for_child_ids=for_child_ids,depth=depth)
        if parents:
            parent_organs = [x for x in parents if x['ObjTypeName'][0] == 'Organ']
            if parent_organs: 
                return parent_organs
            else:
                parent_organs = [x for x in parents if x['ObjTypeName'][0] == 'Tissue']
                if parent_organs: 
                    return parent_organs
                else:
                    return [x for x in parents if x['ObjTypeName'][0] == 'SemanticConcept']

    def add_new_parents(self):
        orphan_ids = self.get_orphan_ids()
        new_parents = set()
        depth = 1
        parent_organs = self.__get_parents(orphan_ids,depth)
                
        while parent_organs:
            new_parents.update(parent_organs)
            new_parent_organs_ids = [x['Id'][0] for x in parent_organs]
            parent_organs_children = ps_api.Graph.get_children_ids(new_parent_organs_ids,at_depth=depth)
            orphan_ids = orphan_ids.difference(parent_organs_children)
            depth += 1
            parent_organs = self.__get_parents(orphan_ids,depth)

        return new_parents


    def print_report(self, to_file:str):
        self.report.sort_values(by=['#References'],ascending=False, inplace=True)
        self.report.to_csv(to_file, sep='\t',index=False)


    def load_graph(self, entity_name:str):
        get_entity = 'SELECT Entity WHERE (Name,Alias) = ({})'.format(entity_name)
        get_anatomy = 'SELECT Entity WHERE objectType = (Organ, Tissue, CellType)'
        oql_query = 'SELECT Relation WHERE objectType = (CellExpression,MolTransport) AND NeighborOf ({}) AND NeighborOf ({})'
        oql_query = oql_query.format(get_entity, get_anatomy)
        request_name = 'Find anatomical concepts for {}'.format(entity_name)
        self.add_rel_props(PS_ID_TYPES)

        self.process_oql(oql_query,request_name)
        self.load_children()


    def make_report(self):
        organs = ps_api.Graph.get_objects(['Organ'])
        organs.sort(key=lambda x: len(x[CHILDS]), reverse=True)
        [self.add2report(organ) for organ in organs]

        new_parents = list(ps_api.add_new_parents())
        [self.add2report(parent) for parent in new_parents]

        tissues = ps_api.Graph.get_objects(['Tissue'])
        [self.add2report(tissue) for tissue in tissues]

        cells = ps_api.Graph.get_objects(['CellType'])
        [self.add2report(cell) for cell in cells]


if __name__ == "__main__":
    entities = ['CNR1','CNR2','GPR55','GPR18','GPR119']
    xslx_file = DATA_DIR+'cannabinoid_receptors_tissue_localization.xlsx'
    writer = ExcelWriter(xslx_file, engine='xlsxwriter')

    for entity in entities:
        ps_api = TissueLocalization(load_api_config())
        ps_api.load_graph(entity)
        ps_api.make_report()
        ps_api.report.sort_values(by=['#References'],ascending=False, inplace=True)
        ps_api.report.df2excel(writer, sheet_name=entity+'_anatomy')

    writer.save()


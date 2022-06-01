import pandas as pd
from pandas import DataFrame as df
import rdfpandas
from ..ResnetAPI.ResnetGraph import PSObject


class df(pd.DataFrame):
    pass
    def __init__(self, p=pd.DataFrame(),**kwargs):
        if kwargs:
            p = pd.DataFrame(**kwargs)
        super(df,self).__init__(p)

    def format_vertical_headers(self):
        """Display a dataframe with vertical column headers"""
        styles = [dict(selector="th", props=[('width', '40px')]),
                dict(selector="th.col_heading",
                    props=[("writing-mode", "vertical-rl"),
                            ('transform', 'rotateZ(180deg)'), 
                            ('height', '290px'),
                            ('vertical-align', 'top')])]
        return (self.fillna('').style.set_table_styles(styles))


    def pandas2rdf(self, remap:dict):
        rdf_pandas = pd.DataFrame(self)
        for c in rdf_pandas.columns:
            rdf_pandas.rename(columns=remap[c])
            return rdfpandas.to_graph(rdf_pandas)


    def apply_and_concat(self, field, func, column_names):
            return pd.concat((self,self[field].apply(lambda cell: pd.Series(func(field,cell),index=column_names))),axis=1)

    #not yet tested
    def apply_and_concat2(self, concept_name, func, column_names):
            return pd.concat((self,self['Name'].apply(lambda cell: pd.Series(func(cell,concept_name),index=column_names))),axis=1)


    @classmethod
    def dict2pd(cls,dic:dict, key_colname:str, value_colname:str)->pd.DataFrame:
        """
        input - single value dic {str:str or float} 
        """
        new_df = cls(pd.DataFrame(columns=[key_colname,value_colname]))
        new_df = pd.DataFrame.from_dict({key_colname:list(dic.keys()), value_colname:list(dic.values())})
        return new_df

    def merge_dict(self, dict2add:dict, new_col:str, map2column:str, add_all=False, case_sensitive_match=False):
        pd2merge = df.dict2pd(dict2add,map2column,new_col)
        in2pd = df(self)
        if not case_sensitive_match:
            pd2merge[map2column].apply(lambda x: str(x).lower())
            in2pd[map2column].apply(lambda x: str(x).lower())

        how = 'outer' if add_all else 'left'
        merged_pd = in2pd.merge(pd2merge,how, on=map2column)
        return merged_pd

    @classmethod
    def psobj2pd(cls,obj:PSObject, key_colname:str, value_colname:str):
        key_col = list()
        value_col = list()
        for k,v_list in obj.items():
            key_col += [k]*len(v_list)
            value_col += v_list
            
        new_df = cls(pd.DataFrame(columns=[key_colname,value_colname]))
        new_df = pd.DataFrame.from_dict({key_colname:key_col,value_colname:value_col})
        return new_df


    def merge_psobject(self, obj:PSObject, new_col:str, map2column:str, 
            add_all=False, values21cell=False, sep=';', case_sensitive_match = False):

        in2pd = df(self)
        if not case_sensitive_match:
            in2pd[map2column].apply(lambda x: str(x).lower())

        if values21cell:
            if case_sensitive_match:
                merge_dict = {k:sep.join(v) for k,v in obj.items()}
            else:
                merge_dict = {str(k).lower():sep.join(v) for k,v in obj.items()}
                
            return in2pd.merge_dict(merge_dict,new_col,map2column,add_all)
        else:          
            obj_pd = df.psobj2pd(obj,map2column,new_col)
            how = 'outer' if add_all else 'left'
            return in2pd.merge(obj_pd,how, on=map2column)


    def get_rows(self, by_value1, in_column1, and_by_value2, in_column2):
        return self.loc[(self[in_column1] == by_value1) & (self[in_column2] == and_by_value2)]

    def get_cells(self,by_value1, in_column1, and_by_value2, in_column2, from_column3):
        return (self.loc[(self[in_column1] == by_value1) & (self[in_column2] == and_by_value2),from_column3]).iloc[0]

    def df2json(self, to_file:str, dir=''):
        dump_fname = to_file
        if  to_file[-5:] != '.json': dump_fname += '.json'
        dump_fname = dir+dump_fname
        with open(dump_fname, 'w') as jsonout:
            self.reset_index(inplace=True)
            jsonout.write(self.to_json(None,indent=2, orient='index'))


    def df2excel(self, writer:pd.ExcelWriter,sheet_name:str, index=False, vertical_header=False):
        self.to_excel(writer, sheet_name=sheet_name, startrow=1, header=False, index=index)
        workbook  = writer.book
        worksheet = writer.sheets[sheet_name]

        header_format = {
                    'bold': True,
                    'text_wrap': True,
                    'valign': 'vjustify',#'bottom',#'top', #vcenter',
                    'align': 'center' # 'left' 
                    #'fg_color': '#D7E4BC',
                    #'border': 1,
                    }

        header_height = 15
        if vertical_header:
            header_format.update({'rotation': '90'})
            header_height = 120

        format = workbook.add_format(header_format)
        worksheet.set_row(0,header_height,format)
        for col_num, value in enumerate(self.columns.values):
            worksheet.write(0, col_num, value,format)


        



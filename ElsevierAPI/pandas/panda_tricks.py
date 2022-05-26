import pandas as pd
from pandas import DataFrame as df
import rdfpandas
from ..ResnetAPI.ResnetGraph import PSObject

class df(pd.DataFrame):
    pass
    def __init__(self, p=pd.DataFrame()):
        super().__init__(p)

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
        return cls.from_dict({key_colname:list(dic.keys()), value_colname:list(dic.values())})

    def merge_dict(self, dict2add:dict, new_col:str, map2column:str, add_all=False):
        pd2merge = self.dict2pd(dict2add,map2column,new_col)
        pd2merge[map2column].apply(lambda x: str(x).lower())
        in2pd = df(self)
        in2pd[map2column].apply(lambda x: str(x).lower())

        how = 'outer' if add_all else 'left'
        merged_pd = in2pd.merge(pd2merge,how, on=map2column)
        return merged_pd

    def psobj2pd(self,obj:PSObject, key_colname:str, value_colname:str):
        self = df(columns=[key_colname,value_colname])
        _dict = dict()

        key_col = list()
        value_col = list()
        for k,v_list in obj.items():
            key_col += [k]*len(v_list)
            value_col += v_list
            _dict[k] = ';'.join(v_list)
            
        self = pd.DataFrame.from_dict({key_colname:key_col,value_colname:value_col})
        return _dict


    def merge_psobject(self, obj:PSObject, new_col:str, map2column:str, add_all=False, values21cell=False):
        _pd2, _dict = self.psobj2pd(obj,map2column,new_col)

        if values21cell:
            return self.merge_dict(self,_dict,new_col,map2column,add_all)
        else:
            _pd2[map2column].apply(lambda x: str(x).lower())
            in2pd = df(self)
            in2pd[map2column].apply(lambda x: str(x).lower())
            how = 'outer' if add_all else 'left'
            merged_pd = in2pd.merge(_pd2,how, on=map2column)
            return merged_pd


    def get_row(by_value1, in_column1, and_by_value2, in_column2, from_pd:df):
        return from_pd.loc[(from_pd[in_column1] == by_value1) & (from_pd[in_column2] == and_by_value2)]

    def get_cell(by_value1, in_column1, and_by_value2, in_column2, from_column3, in_pd:df):
        return (in_pd.loc[(in_pd[in_column1] == by_value1) & (in_pd[in_column2] == and_by_value2),from_column3]).iloc[0]

    def df2json(self, to_file:str, dir=''):
        dump_fname = to_file
        if  to_file[-5:] != '.json': dump_fname += '.json'
        dump_fname = dir+dump_fname
        with open(dump_fname, 'w') as jsonout:
            self.reset_index(inplace=True)
            jsonout.write(self.to_json(None,indent=2, orient='index'))


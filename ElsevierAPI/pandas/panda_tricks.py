import pandas as pd
from pandas.api.types import is_numeric_dtype
from pandas import ExcelWriter as ExcelWriter
import rdfpandas
from ..ResnetAPI.NetworkxObjects import PSObject
import string
from openpyxl import load_workbook

MIN_COLUMN_WIDTH = 0.65 # in inches
NUMBER_OF_REFERENCE = 'Number of references'

def make_hyperlink(identifier:str,url:str,display_str=''):
        # hyperlink in Excel does not work with long URLs, 
        display_str = display_str if display_str else identifier
        return '=HYPERLINK("'+url+identifier+'",\"{}\")'.format(display_str)

class df(pd.DataFrame):
    pass
    header_format = {
                    'bold': True,
                    'text_wrap': True,
                    'valign': 'vjustify',#'bottom',#'top', #vcenter',
                    'align': 'center', # 'left'
                    'height':15
                    #'fg_color': '#D7E4BC',
                    #'border': 1,
                    }

    column2format = dict() # {column_index:{'font_color':'blue'}}
    conditional_frmt = dict() # {area:{conditional_format}}
    

    def __init__(self, *args, **kwargs):
        df_name = kwargs.pop('name','')  
        pd.DataFrame.__init__(self, *args, **kwargs)
        self._name_ = df_name
        self.column2format = dict()
        self.conditional_frmt = dict()

        
    def copy_format(self, from_df:'df'):
        self.header_format = from_df.header_format
        self.column2format = from_df.column2format
        self.conditional_frmt = from_df.conditional_frmt


    def set_col_format(self,fmt:dict):
        self.column2format = fmt


    def add_column_format(self,column_name:str,fmt_property:str,fmt_value:str):
        try:
            column_format = self.column2format[column_name]
        except KeyError:
            column_format = dict()

        column_format.update({fmt_property:fmt_value})
        self.column2format[column_name] = column_format


    def set_hyperlink_color(self, column_names:list):
        [self.add_column_format(column_name,'font_color','blue') for column_name in column_names]


    def __copy_attrs(self,from_df:'df'):
        self._name_ = from_df._name_
        self.copy_format(from_df)


    def inch_width(self,col:str or int):
        my_column = col if isinstance(col,str) else str(self.columns[col])
        try:
            return self.column2format[my_column]['inch_width']
        except KeyError:
            return 0


    def set_inch_width(self,layout:dict):
        """
        layout = {column_index:inch_width}
        """
        for col_idx, width in layout.items():
            col_name = str(self.columns[col_idx])
            self.column2format[col_name] = {'inch_width':width}


    def inch_width_layout(self):
        layout = {i:self.inch_width(i) for i in range(0,len(self.columns))}
        return layout


    @classmethod
    def copy_df(cls, other:'df') ->'df':
        newdf = cls(other.copy())
        newdf.__copy_attrs(other)
        return newdf


    @classmethod
    def read(cls, *args, **kwargs):
        '''
        Input
        -----
        args[0] - input filname for reading Excel file add sheet_name=my_worksheet_name to kwargs
        kwargs = {sheet_name:str, read_formula:bool}
        '''
        df_name = kwargs.pop('name','')
        fname = str(args[0])
        extension = fname[fname.rfind('.')+1:]
        if extension == 'xlsx':
            read_formula = kwargs.pop('read_formula',False)
            if read_formula:
                wb = load_workbook(filename=fname)
                try:
                    sheet = wb[kwargs['sheet_name']]  
                except KeyError:
                    sheet_names = wb.get_sheet_names()
                    sheet = wb[sheet_names[0]]

                header_pos = kwargs.pop('header',0)
                skiprows = kwargs.pop('skiprows',0)
                _df = df(pd.DataFrame(sheet.values))
                _df.columns = _df.iloc[header_pos].to_list()
                _df = df(_df[header_pos+1+skiprows:])
                _df._name_ = df_name if df_name else sheet
                return _df
            try:
                _df = df(pd.read_excel(*args, **kwargs))
                _df._name_ = df_name
                return _df
            except FileNotFoundError: 
                raise FileNotFoundError
        elif extension in ('tsv','txt'):
            kwargs.pop('sheet_name','')
            kwargs.pop('read_formula',False)
            try:
                _df = df(pd.read_csv(fname,sep='\t', **kwargs))
                _df._name_ = df_name
                return _df
            except FileNotFoundError:
                raise FileNotFoundError
        elif extension == ('csv'):
            kwargs.pop('sheet_name','')
            kwargs.pop('read_formula',False)
            try:
                _df = df(pd.read_csv(fname,sep=',', **kwargs))
                _df._name_ = df_name
                return _df
            except FileNotFoundError:
                raise FileNotFoundError


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
    def from_dict(cls, dic:dict,*args, **kwargs):
        df_name = kwargs.pop('name','')
        new_pd = cls(pd.DataFrame.from_dict(dic,*args, **kwargs))
        new_pd._name_ = df_name
        return new_pd


    @classmethod 
    def from_rows(cls, rows:set or list, header:list):
        dic  = {i:list(row) for i,row in enumerate(rows)}
        new_pd = cls(pd.DataFrame.from_dict(dic,orient='index',columns=header))
        return new_pd


    @classmethod
    def dict2pd(cls,dic:dict, key_colname:str, value_colname:str)->'df':
        """
        input - single value dic {str:str or float} 
        """
        #new_df = cls(pd.DataFrame(columns=[key_colname,value_colname]))
        new_df = df(pd.DataFrame.from_dict({key_colname:list(dic.keys()), value_colname:list(dic.values())}))
        return new_df


    def merge_dict(self, dict2add:dict, new_col:str, map2column:str, add_all=False, case_sensitive_match=False):
        pd2merge = df.dict2pd(dict2add,map2column,new_col)
        in2pd = df.copy_df(self)
        if not case_sensitive_match:
            pd2merge[map2column].apply(lambda x: str(x).lower())
            in2pd[map2column].apply(lambda x: str(x).lower())

        how = 'outer' if add_all else 'left'
        merged_pd = df(in2pd.merge(pd2merge,how, on=map2column))
        merged_pd.__copy_attrs(self)
        return merged_pd


    def append_df(self, other:'df'):
        merged_df = df(pd.concat([self,other],ignore_index=True),name=self._name_)
        return merged_df

    def merge_df(self, *args, **kwargs):
        df_name = kwargs.pop('name',self._name_)
        merged_pd = df(self.merge(*args, **kwargs),name=df_name)
        merged_pd.__copy_attrs(self)
        return merged_pd

    @classmethod
    def psobj2pd(cls,obj:PSObject,key_colname:str,value_colname:str,df_name=str()):
        key_col = list()
        value_col = list()
        for k,v_list in obj.items():
            key_col += [k]*len(v_list)
            value_col += v_list
            
        new_df = cls(pd.DataFrame(columns=[key_colname,value_colname]))
        new_df = df(pd.DataFrame.from_dict({key_colname:key_col,value_colname:value_col}),name=df_name)
        return new_df


    def merge_psobject(self, obj:PSObject, new_col:str, map2column:str, 
            add_all=False, values21cell=False, sep=';', case_sensitive_match = False):

        in2pd = df.copy_df(self)
        if not case_sensitive_match:
            in2pd[map2column] = self[map2column].apply(lambda x: str(x).lower())

        if values21cell:
            if case_sensitive_match:
                merge_dict = {k:sep.join(v) for k,v in obj.items()}
            else:
                merge_dict = {str(k).lower():sep.join(v) for k,v in obj.items()}
                
            in2pd =  in2pd.merge_dict(merge_dict,new_col,map2column,add_all)
            in2pd[map2column] = self[map2column]
            return in2pd
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


    def __worksheet_area4(self, first_col:str,last_col:str,first_row=2,last_row=0):
        first_col_idx = self.columns.to_list().index(first_col)
        last_col_idx = self.columns.to_list().index(last_col)

        first_col_letter_idx = string.ascii_uppercase[first_col_idx]
        last_col_letter_idx = string.ascii_uppercase[last_col_idx]

        area = first_col_letter_idx+str(first_row)+':'+last_col_letter_idx
        if last_row:
            area += str(last_row)
        else:
            area += str(len(self))
        
        return area

    
    def set_conditional_frmt(self,conditional_format:dict, for_first_col:str,and_last_col:str,for_first_row=2,and_last_row=0):
        area = self.__worksheet_area4(for_first_col,and_last_col,for_first_row,and_last_row)
        self.conditional_frmt[area] = conditional_format


    def make_header_vertical(self,height=120):
        self.header_format['rotation'] = '90'
        self.header_format['height'] = height
        self.header_format['valign'] = 'vcenter'


    def make_header_horizontal(self):
        self.header_format.pop('rotation','not_found')
        self.header_format.pop('height','not_found')
        self.header_format.pop('valign','not_found')


    def df2excel(self, writer:ExcelWriter,sheet_name:str):
        self.to_excel(writer, sheet_name=sheet_name, startrow=1, header=False, index=False, float_format='%g')
        workbook  = writer.book
        worksheet = writer.sheets[sheet_name]

        # writing header
        header_height = self.header_format.pop('height',15)
        format = workbook.add_format(self.header_format)
        worksheet.set_row(0,header_height,format)
        for col_num, value in enumerate(self.columns.values):
            worksheet.write(0, col_num, value,format)

        # formating columns
        for idx, column_name in enumerate(self.columns.to_list()):
            try:
                col_fmt_dic = dict(self.column2format[column_name])
                width = col_fmt_dic.pop('width', 20)
                wrap = col_fmt_dic.pop('wrap_text', False)
                inch_width = col_fmt_dic.pop('inch_width', 1)
                column_format = workbook.add_format(col_fmt_dic)
                if wrap:
                    column_format.set_text_wrap()
                worksheet.set_column(idx, idx, width, column_format)
                # set_column(first_col, last_col, width, cell_format, options)
            except KeyError:
                continue

        for area, fmt in self.conditional_frmt.items():
            worksheet.conditional_format(area, fmt)


    def clean(self):
        clean_pd = self.dropna(how='all',subset=None)
        clean_pd = clean_pd.drop_duplicates()
        clean_pd = clean_pd.fillna('')
        clean_df = df(clean_pd, name=self._name_)
        clean_df.copy_format(self)
        return clean_df


    def filter_by(self,values:list, in_column:str):
        return df(self[self[in_column].isin(values)])


    def greater_than(self, value:float, in_column:str):
        return df(self[self[in_column] >= value])


    def drop_empty_columns(self, max_abs=0.00000000000000001, subset=list()):
        my_columns = subset if subset else self.columns
        my_numeric_columns = subset if subset else [x for x in my_columns if is_numeric_dtype(self[x])]
        no_value_cols = list()
        [no_value_cols.append(col) for col in my_numeric_columns if abs(self[col].max()) < max_abs]
        no_empty_cols = df(self.drop(columns = no_value_cols))
        no_empty_cols.copy_format(self)
        no_empty_cols._name_ = self._name_
        print('%s columns were dropped  because they have all values = 0' % no_value_cols)
        return no_empty_cols
        

    def table_layout(self, table_witdh=7.5):
        """
        Returns {col_idx:width} for default table_witdh=7.5
        """
        length_averages = dict()
        columns = list(self.columns)
        rownum = float(len(self.index))
        col_idx = 0
        for col in columns:
            col_values = list(self[col])+[col]
            col_values = list(map(str,col_values))
            split_col_values = list()
            for v in col_values:
                sentences = v.split('\n')
                split_col_values += sentences
            col_lengths = list(map(len, split_col_values))
            length_averages[col_idx] = float(sum(col_lengths))/rownum
            col_idx +=1
        
        scale_factor = table_witdh/sum(length_averages.values())
        normalized_layout = {k:v*scale_factor for k,v in length_averages.items()}

        # columns cannot be too narrow
        #min_width = 0.65
        [normalized_layout.update({i:MIN_COLUMN_WIDTH}) for i,w in normalized_layout.items() if w < MIN_COLUMN_WIDTH]

        # trying to fit the table into the page by narrowing very wide columns
        new_table_width = sum(list(normalized_layout.values()))
        width_excess = new_table_width-table_witdh
        for col_idx,col_width in normalized_layout.items():
            if col_width > 3:
                new_width = col_width - width_excess
                if new_width > 2:
                    normalized_layout[col_idx]= new_width
                break

        new_table_width = sum(list(normalized_layout.values()))
        width_excess = new_table_width-table_witdh
        if width_excess > 0:
            #distributing width excess across all columns wider than min_col_width
            wide_col_idxes = list()
            sum_width2trim = 0.0
            for idx,width in normalized_layout.items():
                if width > MIN_COLUMN_WIDTH:
                    wide_col_idxes.append(idx)
                    sum_width2trim =  sum_width2trim + width

            scale_factor = sum_width2trim/(sum_width2trim+width_excess)
        
            for col_idx,col_width in normalized_layout.items():
                if col_idx in wide_col_idxes:
                    normalized_layout[col_idx] = col_width*scale_factor
        
        #table_width = sum(list(normalized_layout.values()))
        # assert table_width <= 6.8
        return normalized_layout


    def sort_columns_by_list(self,only_columns:list):
        my_columns = [c for c in only_columns if c in self.columns]
        new_df = df(self[my_columns])
        new_df.copy_format(self)
        return new_df


    def clean4doc(self,max_row:int=None,only_columns=[],ref_limit=dict(), as_str=True):    
        clean_df = self.sort_columns_by_list(only_columns) if only_columns else df.copy_df(self)
        clean_df.dropna(how='all',subset=None,inplace=True)
        clean_df.drop_duplicates(inplace=True)
        clean_df.fillna('', inplace=True)

        if max_row: clean_df = df(clean_df.head(max_row))

        if ref_limit:
            for col_name, cutoff in ref_limit.items():
                clean_df = df(clean_df.loc[clean_df[col_name] >= cutoff])

        if as_str:
            for col in clean_df.columns:
                if clean_df[col].dtype in ['float','float64']:
                    clean_df[col] = clean_df[col].map(lambda x: '%2.2f' % x)
                clean_df[col] = clean_df[col].astype(str)
        return clean_df


    @staticmethod
    def read_clean(*args,**kwargs):
        """
        Input
        -----
        file name i sin args[0]
        'only_columns' = [col_names]
        'ref_limit' = [col_name:reflimit]
        'as_str' = True - formats floats to string as %2.2f'
        max_row = 0 - read only first max_row rows
        other *args,**kwargs as in pandas.read_excel, pandas.read_excel
        """
        only_columns = kwargs.pop('only_columns',[])
        ref_limit = kwargs.pop('ref_limit',dict())
        as_str = kwargs.pop('as_str',True)
        max_row = kwargs.pop('max_row',0)
        try:
            table_df = df.read(*args,**kwargs)
            if table_df.empty:
                print('worksheet %s in %s is empty' % (table_df._name_,args[0]))
                return df()
            
            if only_columns:
                if isinstance(only_columns[0],int):
                    select_columns = [v for i,v in enumerate(table_df.columns.to_list()) if i in only_columns]
                else:
                    select_columns = only_columns
            else:
                select_columns = []

            return table_df.clean4doc(max_row,select_columns,ref_limit,as_str) 
        except FileNotFoundError: 
            raise FileNotFoundError


    @staticmethod
    def formula2hyperlink(formula:str):
        '''
        extracts hyperlink from HYPERLINK formula in Excel
        '''
        values = formula[formula.find('(')+1:-2]
        url, txt = values.split(',')
        return url.strip(' "'), txt.strip(' "')


    def remove_rows_by(self, values:list, in_column:str):
        clean_df = df(self[~self[in_column].isin(values)])
        clean_df.copy_format(self)
        return clean_df
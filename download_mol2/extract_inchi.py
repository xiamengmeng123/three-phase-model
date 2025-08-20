###step 2 从docs.csv中抽取inchi信息并进行排序及字符串处理，得到标准的分子id,inchi的inchi.csv文件
###修改文件夹名csv_name和开始序号df.index
import os 
import pandas as pd
# import csv


cif_path = "D:/vspractice/cif_file/"
csv_name = "3-excludeMgLi" 
df = pd.read_csv(f'{cif_path}docs-{csv_name}.csv')

def extract_inchi(row):
    # start_index = row.find('InChI=') + len('InChI=')
    start_index = row.find('InChI=')
    end_index = row.find("'",start_index)
    inchi_key = row[start_index:end_index]
    return inchi_key 
    

#提取第29列数据，并对每行进行字符串处理
column_29 = df.iloc[:,29].apply(extract_inchi)
column_29.to_csv(f'{cif_path}inchi-{csv_name}.csv',index=False,header=False)   #csv里有逗号需要加双引号""以防止解析成分隔符，仍看做整体
# column_29.to_csv('inchi-5-excludeMgLi.csv',index=False,quoting=csv.QUOTE_NONE,escapechar='\\')   #这里加转义符可以把双引号去掉，但没必要

#删除重复行
df = pd.read_csv(f'{cif_path}inchi-{csv_name}.csv',header=None)
df.drop_duplicates(inplace=True)
df.reset_index(drop=True,inplace=True)
df.index += 128167
df.to_csv(f'{cif_path}inchi-{csv_name}.csv',index=True,header=False)






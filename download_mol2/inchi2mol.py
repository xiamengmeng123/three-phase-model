###step 3 从inchi.csv文件读取分子id和inchi信息，并转换成mol格式文件
###修改变量ele和df的输入文件
import  os
import pandas as pd
from rdkit import Chem
import rdkit.Chem

cif_path = "D:/vspractice/cif_file/"
ele = "F:/MP/3456/"
# ele = os.path.join(cif_path,"three_ele")

df = pd.read_csv(f'{cif_path}inchi-3-excludeMgLi.csv',header=None)
#遍历DataFrame中的每一行
invalid_inchi_list = []
for index,row in df.iterrows():
    file_name = f'{row[0]}.mol'
    INCHI = row[1]
    try:
        mol = Chem.MolFromInchi(INCHI)
        if mol is None:
            invalid_inchi_list.append(INCHI)
        else:
            molblock = Chem.MolToMolBlock(mol)
            with open(f'{ele}/{file_name}','w') as file:
                file.write(molblock)
    except Exception as e:
        print(f"ERROR processing InChI: {INCHI}")
        invalid_inchi_list.append(INCHI)

if invalid_inchi_list:
    print("Invalid InChI:")
    invalid_df = pd.DataFrame(invalid_inchi_list,columns=['Invalid InChI'])
    invalid_df.to_csv(f'{cif_path}invalid_inchi.csv',index=False)

#单个inchi转mol，用于查看invalid_inchi.csv里的无效inchi
# INCHI = "InChI=1S/C3H3FN4O/c1-8-6-2(4)3(9)5-7-8/h1H3"
# mol = Chem.MolFromInchi(INCHI)
# molblock = Chem.MolToMolBlock(mol)
# print(molblock)
# with open(f'{ele}/xx.mol','w') as file:
#     file.write(molblock)

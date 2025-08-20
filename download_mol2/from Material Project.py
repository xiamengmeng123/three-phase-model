###step1 读取Material Project的搜索分子的docs信息,生成docs.csv文件,该过程较慢
import os
from mp_api.client import MPRester
import rdkit
from rdkit import Chem
import rdkit.Chem
import matplotlib.pyplot as plt
import pandas as pd
from requests.exceptions import ReadTimeout
from requests.adapters import HTTPAdapter
import requests

# 在绘图之前关闭交互模式
plt.ioff()

save_directory = "D:/vspractice/cif_file/"

#搜索分子的docs信息
with MPRester(api_key='wLy5pd6Ks27Xy76HygnKMJqBgtB4DWHn') as mpr:
    try:
        # docs = mpr.molecules.summary.search(molecule_ids=["08dfdd83319a24bf2a0c31db2f83e909-C6H7N1O2-1-2","1b9596a8277f19bb553b58d2e67a1f2f-C6H7N1O2-1-2"])
        # docs = mpr.molecules.summary.search(molecule_ids=["14e9a3c9c48c9cc3b6b9c7b92401335e-C2F3N1S2-m2-1","f9979ab58248dc4966cbb56481d647d2-C1Mg1N2O1S1-m1-2","4b78910e8f83b3b7d4b60b6f4d479fac-C4F1H2Li1O4-0-1"])#,fields=["molecule_ids","inchi"])
        
        #常用元素CHONFMgClPSLiB
        #docs = mpr.molecules.summary.search(nelements=4,exclude_elements=["Mg","S","P","Li"],fields=["molecule_id","inchi"])
        docs = mpr.molecules.summary.search(nelements=3,exclude_elements=["Mg","Li"],fields=["molecule_id","inchi"])
        df = pd.DataFrame(docs)
        df.to_csv(f'{save_directory}docs.csv',index=False)
        # id_list = []
        # icc = 0
        # ele = os.path.join(save_directory,"four_ele")
        # #将inchi格式转换成mol文件
        # for doc in docs[0,10]:#获取材料索引号        
        #     id_list.append(doc.inchi)
        #     INCHI = doc.inchi
        #     if INCHI is not None:
        #         print('=============='+INCHI+'==========')
        #         mol=rdkit.Chem.inchi.MolFromInchi(INCHI)
        #         if mol is not None:
        #             icc += 1
        #             molblock = Chem.MolToMolBlock(mol)
        #             # print(molblock)
        #             with open(f'{ele}/{icc}.mol','w') as file:
        #                 file.write(molblock)
        #             # print(molblock,file=open('{icc}.mol','w+'))
        #             print(icc,'  ',doc.molecule_id,'  ',doc.inchi,file=open(f'{save_directory}inchi_list.txt','a'))
        #         else:
        #             None
        #     else:
        #         print('INCHI is None,skipping thie entry')
    except ReadTimeout as e:
        print(f"Read timeout occurred: {e}")
        print("Retrying...")






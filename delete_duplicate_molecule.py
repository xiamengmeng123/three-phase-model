import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import numpy as np
from tqdm import tqdm  # 用于显示进度条

# 读取CSV文件（无表头）
df = pd.read_csv('D:/vspractice/csv/mol_unimol_fep.csv', header=None, names=['molecule_id', 'SMILES', 'target'])
print(f"原始分子数量: {len(df)}")

# 函数：生成分子指纹
def get_fingerprint(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)

# 生成分子指纹
print("正在生成分子指纹...")
df['fingerprint'] = df['SMILES'].apply(get_fingerprint)

# 移除无效的SMILES
df = df.dropna(subset=['fingerprint']).reset_index(drop=True)
print(f"有效分子数量: {len(df)}")

# 设置相似度阈值 (0.8 = 80%相似度)
SIMILARITY_THRESHOLD = 0.8

# 初始化保留列表
keep_mask = np.ones(len(df), dtype=bool)  # 初始全部保留

# 使用进度条进行逐对比较
print("正在进行分子相似性比较...")
for i in tqdm(range(len(df))):
    if not keep_mask[i]:  # 如果分子已被标记为删除，则跳过
        continue
        
    fp_i = df.iloc[i]['fingerprint']
    
    # 比较当前分子与所有后续分子
    for j in range(i + 1, len(df)):
        if not keep_mask[j]:  # 跳过已标记删除的分子
            continue
            
        fp_j = df.iloc[j]['fingerprint']
        similarity = DataStructs.TanimotoSimilarity(fp_i, fp_j)
        
        # 如果相似度超过阈值，标记为删除
        if similarity > SIMILARITY_THRESHOLD:
            keep_mask[j] = False

# 创建去重后的DataFrame
deduplicated_df = df[keep_mask].copy()
deduplicated_df = deduplicated_df.drop(columns=['fingerprint'])

# 保存结果
print(f"去重后分子数量: {len(deduplicated_df)}")
deduplicated_df.to_csv('D:/vspractice/csv/mol_free_energy.csv', index=False, header=False)
print("结果已保存")
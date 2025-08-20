import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

# === 路径设置 ===
base_dir = "/mmx/unimol"
cluster_file = f"{base_dir}/csv/train_fep.csv" 
experiment_file = f"{base_dir}/csv/test_fep.csv"         
output_file = "external_train_vs_test_similarity.csv"   # 输出结果csv
histogram_file = "similarity_distribution.png"         # 相似度分布图

# === 1. 读取代表分子（无表头） ===
df_cluster = pd.read_csv(cluster_file, header=None)
cluster_smiles = df_cluster[0].tolist()
# cluster_smiles = df_cluster.iloc[:, 1].dropna().tolist()

# 过滤代表分子中的无效SMILES，计算指纹
valid_cluster_smiles = []
cluster_fps = []
for smi in cluster_smiles:
    mol = Chem.MolFromSmiles(smi)
    if mol is not None:
        valid_cluster_smiles.append(smi)
        cluster_fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048))
    else:
        print(f"Warning: Invalid cluster SMILES filtered out: {smi}")

print(f"有效代表分子数量: {len(valid_cluster_smiles)}")

# === 2. 读取分子（有表头，第二列为SMILES） ===
df_exp = pd.read_csv(experiment_file,header=None,dtype=str)
# exp_smiles_list = df_exp.iloc[:, 1].tolist()
exp_smiles_list = (
    df_exp.iloc[:, 0]
        .dropna()
        .str.strip()
        .tolist()
)

# 过滤分子中的无效SMILES
valid_exp_smiles = []
for smi in exp_smiles_list:
    mol = Chem.MolFromSmiles(smi)
    if mol is not None:
        valid_exp_smiles.append(smi)
    else:
        print(f"Warning: Invalid experiment SMILES filtered out: {smi}")

print(f"有效分子数量: {len(valid_exp_smiles)}")

# === 3. 计算每个分子与所有代表分子的最大 Tanimoto 相似度 ===
results = []

for exp_smi in tqdm(valid_exp_smiles, desc="Comparing similarities"):
    exp_mol = Chem.MolFromSmiles(exp_smi)
    exp_fp = AllChem.GetMorganFingerprintAsBitVect(exp_mol, 2, 2048)
    sims = [DataStructs.TanimotoSimilarity(exp_fp, cluster_fp) for cluster_fp in cluster_fps]
    max_sim = max(sims) if sims else None
    results.append((exp_smi, max_sim))

# === 4. 保存相似度结果 ===
df_result = pd.DataFrame(results, columns=["experiment_smiles", "tanimoto_similarity"])
df_result.to_csv(output_file, index=False)
print(f"相似度结果已保存到 {output_file}")

# === 5. 绘制相似度分布直方图 ===
# df_result = pd.read_csv(output_file)
df_valid = df_result[df_result["tanimoto_similarity"].notna()]
df_valid = df_valid[df_valid["tanimoto_similarity"] <= 0.8].copy()


plt.figure(figsize=(8, 5))  
sns.histplot(df_valid["tanimoto_similarity"], bins=20, kde=True, color="steelblue")
plt.xlabel("Max Tanimoto Similarity")
plt.ylabel("Number of  Molecules")
# plt.title("Similarity Distribution: Experiment vs. Training set")
plt.axvline(0.85, color='red', linestyle='--', label="Typical Similarity Threshold")
plt.legend()
plt.tight_layout()
plt.savefig(histogram_file, dpi=300)
print(f"相似度分布图已保存到 {histogram_file}")

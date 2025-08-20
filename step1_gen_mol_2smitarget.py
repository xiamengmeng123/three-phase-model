########
#若 xx.csv 存在：输出 molecule_id, SMILES, TARGET 三列。

#若 xx.csv 不存在：输出仅 molecule_id, SMILES 两列

import os
import pandas as pd
from rdkit import Chem
import numpy as np

# === 路径参数 ===
base_dir = ''
energy_path = f"{base_dir}/csv/xx.csv"
mol_data_dir = "mol_data"
out_csv = "csv/mol_interaction_energy.csv"

def mol_to_smiles(mol_file: str) -> str | None:
    """将 .mol 文件转换为 SMILES 字符串"""
    try:
        mol = Chem.MolFromMolFile(mol_file, sanitize=True)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol)
    except Exception as e:
        print(f"[WARN] 读取失败 {mol_file}: {e}")
        return None

def build_from_existing_csv(path: str) -> pd.DataFrame:
    """根据已有的 molecule_id-TARGET 表补全 SMILES 列"""
    df = pd.read_csv(path, names=["molecule_id", "TARGET"], header=None)
    df.dropna(how="all", inplace=True)
    df.dropna(subset=["molecule_id", "TARGET"], how="any", inplace=True)
    df.drop_duplicates(subset="molecule_id", inplace=True)
    df.sort_values(by="TARGET", inplace=True)
    df.reset_index(drop=True, inplace=True)

    smiles_list = []
    for mol_id in df["molecule_id"]:
        mol_file = os.path.join(mol_data_dir, f"{int(mol_id)}.mol")
        smiles = mol_to_smiles(mol_file)
        smiles_list.append(smiles)

    df.insert(1, "SMILES", smiles_list)
    return df

def build_from_mol_folder(folder: str) -> pd.DataFrame:
    """如果 xx.csv 不存在，从 mol 文件夹生成 molecule_id‑SMILES 表"""
    data = []
    for fname in os.listdir(folder):
        if fname.endswith(".mol"):
            mol_id = os.path.splitext(fname)[0]
            mol_file = os.path.join(folder, fname)
            smiles = mol_to_smiles(mol_file)
            if smiles is not None:
                data.append([mol_id, smiles])
    df = pd.DataFrame(data, columns=["molecule_id", "SMILES"])
    df.sort_values(by="molecule_id", inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df

def main():
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)

    if os.path.isfile(energy_path):
        print(f"[INFO] 检测到 {energy_path}，使用带 TARGET 的处理逻辑。")
        df_out = build_from_existing_csv(energy_path)
        df_out.to_csv(out_csv, header=None, index=False)
    else:
        print(f"[INFO] 未检测到 {energy_path}，从 mol_data 文件夹生成 molecule_id 与 SMILES。")
        df_out = build_from_mol_folder(mol_data_dir)
        df_out.to_csv(out_csv, header=None, index=False)

    print(f"[SUCCESS] 结果保存至 {out_csv}，共 {len(df_out)} 条分子。")
    print(df_out.head())

if __name__ == "__main__":
    main()

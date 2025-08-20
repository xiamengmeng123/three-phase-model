##yhrun -N 1 -p v100 --cpus-per-gpu=4 --gpus-per-node=1  python step2_regression.py
import pandas as pd
from rdkit import Chem
from unimol_tools import MolPredict

# 路径
validation_path = '/mmx/unimol/csv/experiment.csv'
output_csv = '/mmx/unimol/csv/experiment_cleaned.csv'
save_dir = '/mmx/unimol/script/save_predict/'

# ---------- 1. 读取并重命名列 ----------
# data = pd.read_csv(validation_path, usecols=['molecule_id', 'SMILES'])
# data.columns = ["molecule_id", "SMILES"]
data = pd.read_csv(validation_path, header=None, names=["molecule_id", "SMILES"])


# ---------- 2. SMILES 合法性筛选 ----------
def is_valid_smiles(smi):
    try:
        return Chem.MolFromSmiles(smi) is not None
    except:
        return False

data = data[data['SMILES'].apply(is_valid_smiles)].reset_index(drop=True)
data.to_csv(output_csv, index=False)  # 保存筛选后的数据

# ---------- 3. 加载模型并预测 ----------
clf = MolPredict(load_model='./clf_interaction')
validation_pred = clf.predict(output_csv)

# ---------- 4. 结果保存 ----------
data['pred'] = validation_pred
clf.save_predict(data=data, dir=save_dir, prefix='step2_prediction_experiment')

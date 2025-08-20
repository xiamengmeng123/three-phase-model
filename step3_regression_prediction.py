#yhrun -N 1 -p v100 --cpus-per-gpu=4 --gpus-per-node=1  python step3_regression.py
import pandas as pd
from unimol_tools import MolPredict

# 1. 读取初始预测结果
validation_path = '/mmx/unimol/script/save_predict/step2_prediction_experiment.predict.1.csv'
# data = pd.read_csv(validation_path, usecols=['molecule_id', 'SMILES', 'pred'])
data = pd.read_csv(validation_path, usecols=['molecule_id', 'SMILES','pred'])

# 2. 筛选 pred < -100 的样本 可视化时紧密结合
filtered_data = data[data['pred'] < -100].reset_index(drop=True)

# 3. 保存为临时 CSV 用于下一阶段模型
temp_input_path = '/mmx/unimol/script/save_predict/temp_filtered_input.csv'
filtered_data[['molecule_id', 'SMILES']].to_csv(temp_input_path, index=False)

# 4. 加载第二阶段模型
clf = MolPredict(load_model='./clf_fep')

# 5. 用第二阶段模型进行 FEP 预测
fep_pred = clf.predict(temp_input_path)
filtered_data['fep'] = fep_pred

# 6. 按 fep 排序
filtered_data = filtered_data.sort_values(by='fep').reset_index(drop=True)

# 7. 保存最终预测结果
clf.save_predict(
    data=filtered_data,
    dir='/mmx/unimol/script/save_predict/',
    prefix='step3_prediction_experiment'
)

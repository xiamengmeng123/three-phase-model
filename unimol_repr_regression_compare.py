#########  对相互作用能数据进行线性回归和Uni-Mol训练，并对比残差图
#输入：包含相互作用能的csv文件
#输出：残差图
#平台：bohrium使用


import  os 
import pandas as pd 
import numpy as np
data = pd.read_csv('csv/mol_interaction_energy.csv',sep=',').iloc[:,[1,2]]#.sample(100)
print("--------------------Original data---------------------")
print(data.head())
# data.columns = ["molecule_id","SMILES","TARGET"]  
data.columns = ["SMILES","TARGET"]

#将数据集的80%设置为训练数据集，20%设置为测试数据集
train_fraction = 0.8
train_data = data.sample(frac=train_fraction,random_state=1)
train_data.to_csv("csv/train.csv",index=False)
test_data = data.drop(train_data.index)
test_data.to_csv("csv/test.csv",index=False)

#设定训练/测试目标
train_y = np.array(train_data["TARGET"].values.tolist())
test_y = np.array(test_data["TARGET"].values.tolist())


#创建一个用来存后面结果的results
results = {}

#可视化结果
import matplotlib.pyplot as plt 
import seaborn as sns
font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 15}
plt.hist(train_data["TARGET"], bins=20, label="Train Data")
plt.hist(test_data["TARGET"], bins=20, label="Test Data")
plt.ylabel("Count", fontdict=font)
plt.xlabel("energy", fontdict=font)
plt.legend()
# plt.show()


#2D指纹表征和unimol分子repr表征向量
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from unimol_tools import UniMolRepr

def calculate_2dqsar_repr(smiles):
    mol = Chem.MolFromSmiles(smiles)
    #计算morgan指纹（半径为3，长度为1024位）
    fp = AllChem.GetMorganFingerprintAsBitVect(mol,3,nBits=1024)
    return np.array(fp)


#单个数据分子表示读取
# def calculate_unimol_qsar_repr(smiles):
#     clf = UniMolRepr(data_type='molecule', remove_hs=False)
#     smiles_list = [smiles]
#     unimol_repr = clf.get_repr(smiles_list, return_atomic_reprs=True)
#     fp = np.array(unimol_repr['cls_repr'])
#     return fp

# train_data["unimol_qsar_mr"] = train_data["SMILES"].apply(calculate_unimol_qsar_repr)
# test_data["unimol_qsar_mr"] = test_data["SMILES"].apply(calculate_unimol_qsar_repr)

#多个数据连续分子表示读取
def calculate_unimol_qsar_repr(data):
    clf = UniMolRepr(data_type='molecule', remove_hs=False)
    smiles_list = data["SMILES"].tolist()  #100个数据
    repr_dict = clf.get_repr(smiles_list,return_atomic_reprs=True)  # 这里的 clf 需要提前定义
    unimol_repr_list = np.array(repr_dict['cls_repr'])
    unimol_repr = repr_dict['cls_repr']
    return unimol_repr

# #2D表示
# train_data["2dqsar_mr"] = train_data["SMILES"].apply(calculate_2dqsar_repr) 
# test_data["2dqsar_mr"] = test_data["SMILES"].apply(calculate_2dqsar_repr)

# unimol表示
train_data["unimol_qsar_mr"] = calculate_unimol_qsar_repr(train_data)  #(100,512)数组转成pandas的列
test_data["unimol_qsar_mr"] = calculate_unimol_qsar_repr(test_data)


print(train_data["unimol_qsar_mr"].iloc[:1].values.tolist())

#分类
import numpy as np 
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.neural_network import MLPRegressor
from xgboost import XGBRegressor
from sklearn.linear_model import Ridge, Lasso, ElasticNet
from sklearn.svm import SVR
from sklearn.ensemble import GradientBoostingRegressor
from lightgbm import LGBMRegressor
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score,mean_absolute_error
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV

#将训练和测试数据转换成numpy数组
train_x = np.array(train_data["unimol_qsar_mr"].values.tolist())
train_y = np.array(train_data["TARGET"].values.tolist())
test_x = np.array(test_data["unimol_qsar_mr"].values.tolist())
test_y = np.array(test_data["TARGET"].values.tolist())
print(train_x.shape)

# # 重塑数据维度
train_x = np.reshape(train_x, (train_x.shape[0], -1))
test_x = np.reshape(test_x, (test_x.shape[0], -1))

# #数据归一化
# scaler = StandardScaler()
# train_x_normalized = scaler.fit_transform(train_x)
# test_x_normalized = scaler.transform(test_x)


#加入超参数
regressors = [
    ("Linear", LinearRegression(),{}), # 线性回归模型  Linear Regression
    ("RR", Ridge(random_state=42),{'alpha': [0.1, 1, 10]}), # 岭回归模型 Ridge Regression
    ("Lasso", Lasso(random_state=42),{'alpha': [0.1, 1, 10]}), # Lasso回归模型 Lasso Regression
    ("ER", ElasticNet(random_state=42),{'alpha': [0.1, 1, 10],'l1_ratio': [0.1,0.5,0.9]}), # ElasticNet回归模型 ElasticNet Regression
    ("SV", SVR(),{'C': [0.1,1,10]}),  # 支持向量回归模型 Support Vector
    ("K-NN", KNeighborsRegressor(),{'n_neighbors': [3,5,7]}),  # K-最近邻回归模型 K-Nearest Neighbors
    ("DT", DecisionTreeRegressor(random_state=42),{'max_depth': [3,5,7]}),  # 决策树回归模型 Decision Tree
    ("RF", RandomForestRegressor(random_state=42),{'n_estimators': [50,100,200]}), # 随机森林回归模型 Random Forest
    ("GB", GradientBoostingRegressor(random_state=42),{'n_estimators': [50,100,200]}), # 梯度提升回归模型 Gradient Boosting
    ("XGB", XGBRegressor(random_state=42),{'n_estimators': [50,100,200]}), # XGBoost回归模型  XGBoost
    ("LGBM", LGBMRegressor(random_state=42),{'num_leaves': [31,63,127]}), # LightGBM回归模型 LightGBM
    ("MLP", MLPRegressor( # 多层感知器（神经网络）回归模型 Multi-layer Perceptron
        hidden_layer_sizes=(128,64,32),
        learning_rate_init=0.0001,
        activation='relu', solver='adam', 
        max_iter=10000, random_state=42),{'hidden_layer_sizes': [(64,),(128,64,32)]}),
]

with open('csv/performance_regression.txt','a') as f:
    f.write(f"===============================================\n")

# 对每个回归器进行训练和预测，并计算各项性能指标
for name, regressor, param_grid in regressors:
    # 对每个回归器进行交叉验证网格搜索
    grid_search = GridSearchCV(regressor,param_grid,cv=5) #使用5折交叉验证
    grid_search.fit(train_x,train_y)
    best_params = grid_search.best_params_
    best_model = grid_search.best_estimator_
    # 训练回归器
    # regressor.fit(train_x, train_y)
    best_model.fit(train_x,train_y)
    # 预测训练数据和测试数据
    pred_train_y = best_model.predict(train_x)
    pred_test_y = best_model.predict(test_x)
    # 将预测结果添加到训练数据和测试数据中
    train_data[f"{name}_pred"] = pred_train_y
    test_data[f"{name}_pred"] = pred_test_y
    # 计算测试数据的性能指标
    mse = mean_squared_error(test_y, pred_test_y)
    mae = mean_absolute_error(test_y,pred_test_y)
    r2 = r2_score(test_y,pred_test_y)
    se = abs(test_y - pred_test_y)
    results[f"{name}"] = {"MSE": mse, "MAE": mae, "R^2": r2, "error": se}
    print(f"[unimol-QSAR][{name}]\tMSE:{mse:.4f}\tMAE:{mae:.4f}\tR^2:{r2:.4f}")

       
    with open('csv/performance_regression.txt','a') as f:
        f.write(f"[unimol-QSAR][{name}]\tMSE:{mse:.4f}\tMAE:{mae:.4f}\tR^2:{r2:.4f}\n")

# 绘制残差图
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
import seaborn as sns

# 绘制残差图
residuals_data = []
for name, result in results.items():
    # if name.startswith("unimol-QSAR"):
    model_residuals = pd.DataFrame({"Model": name, "Error": result["error"]})
    residuals_data.append(model_residuals)

residuals_df = pd.concat(residuals_data, ignore_index=True)
residuals_df.sort_values(by="Error", ascending=True, inplace=True)
model_order = residuals_df.groupby("Model")["Error"].median().sort_values(ascending=True).index

# 使用seaborn绘制violinplot
plt.figure(figsize=(10, 7), dpi=300)
font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 15}
sns.boxplot(y="Model", x="Error", data=residuals_df, order=model_order)
plt.yticks(rotation=45,fontsize=10)
plt.xlabel("Abs Error", fontdict=font)
plt.ylabel("Models", fontdict=font)
plt.savefig("figure/cancha_plot.png",dpi=300)

# 针对每个回归器绘制实际和预测值散点图
for name, regressor, param_grid in regressors:
    model_name = name
    pred_values = test_data[f"{model_name}_pred"]
    actual_values = test_y

    plt.figure(figsize=(8,6))
    plt.scatter(actual_values,pred_values,color='blue',alpha=0.5)

    #添加拟合曲线
    z = np.polyfit(actual_values,pred_values,1)
    p = np.poly1d(z)
    plt.plot(actual_values,p(actual_values),color='red')

    plt.title(f"classification for {model_name}")
    plt.xlabel("Actual values")
    plt.ylabel("Predicted values")
    plt.savefig(f"figure/{model_name}_regression.png",dpi=300)

######uni-mol模型训练
from unimol_tools import MolTrain,MolPredict
from sklearn.metrics import mean_squared_error

# clf = MolTrain(task='regression',
#                 data_type='molecule',
#                 epochs=100,
#                 learning_rate=0.0001,
#                 batch_size=32,
#                 early_stopping=10,
#                 metrics='mse',
#                 split='random',
#                 save_path='./clf',
#                 )
# clf.fit('csv/train.csv')

predm = MolPredict(load_model='./clf_interaction')
pred_train_y = predm.predict('csv/train.csv').reshape(-1)
pred_test_y = predm.predict('csv/test.csv').reshape(-1)

mse = mean_squared_error(test_y,pred_test_y)
mae = mean_absolute_error(test_y,pred_test_y)
r2 = r2_score(test_y,pred_test_y)
se = abs(test_y - pred_test_y)
results[f'Uni-Mol'] = {"MSE": mse, "MAE": mae, "R^2": r2, "error": se}
print(f"[Uni-Mol]\tMSE:{mse:.4f}\tMAE:{mae:.4f}\tR^2:{r2:.4f}")

with open('csv/performance_regression.txt','a') as f:
        f.write(f"[Uni-Mol]\tMSE:{mse:.4f}\tMAE:{mae:.4f}\tR^2:{r2:.4f}\n")

#绘制残差图
residuals_data = []
for name, result in results.items():
    if name.startswith("Uni-Mol"):
        model_residuals = pd.DataFrame({"Model": name, "Error": result["error"]})
        residuals_data.append(model_residuals)
residuals_df = pd.concat(residuals_data, ignore_index=True)
residuals_df.sort_values(by="Error",ascending=True,inplace=True)
model_order = residuals_df.groupby("Model")["Error"].median().sort_values(ascending=True).index

# 使用searborn 绘制Uni-Mol模型的violinplot
plt.figure(figsize=(10,7),dpi=300)
font = {'family': 'serif',
        'color': 'black',
        'weight': 'normal',
        'size': 8}

sns.boxplot(y="Model",x="Error",data=residuals_df, order=model_order)
plt.yticks(rotation=45,fontsize=10)
plt.xlabel("Abs Error", fontdict=font)
plt.ylabel("Models",fontdict=font)
plt.savefig(f"figure/Uni-Mol_regression.png",dpi=300)


#横向比较残差图
residuals_data = []
for name,result in results.items():
    # if result["MSE"] > 10:
    #     continue
    model_residuals = pd.DataFrame({"Model": name,"Error": result["error"]})
    residuals_data.append(model_residuals)

residuals_df = pd.concat(residuals_data, ignore_index=True)
residuals_df.sort_values(by="Error",ascending=True,inplace=True)
model_order = residuals_df.groupby("Model")["Error"].median().sort_values(ascending=True).index

plt.figure(figsize=(10,7),dpi=300)
font = {'family': 'serif',
        'color': 'black',
        'weight': 'normal',
        'size': 15}
sns.boxplot(y="Model",x="Error",data=residuals_df,order=model_order)
plt.yticks(rotation=45,fontsize=10)
plt.xlabel("Abs Error", fontdict=font)
plt.ylabel("Models",fontdict=font)
plt.savefig(f"figure/cancha_vs_regression.png",dpi=300)
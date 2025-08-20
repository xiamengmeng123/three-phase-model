############## 用Uni-Mol训练框架，并画出测试图结果

import pandas as pd
from unimol_tools import MolTrain,MolPredict
from sklearn.metrics import mean_squared_error,r2_score
import matplotlib.pyplot as plt 
import numpy as np 
import joblib


####数据处理
data = pd.read_csv('mol_free_energy.csv',sep=',').iloc[:,[1,2]]#.sample(100)#.iloc[0:48236,[1,2]]#
print(data.head())
data.columns = ["SMILES","TARGET"]

#将数据集的80%设置为训练数据集，20%设置为测试数据集
train_fraction = 0.8
train_data = data.sample(frac=train_fraction,random_state=1)
train_data.to_csv("train_fep_2.csv",index=False)
test_data = data.drop(train_data.index)
test_data.to_csv("test_fep_2.csv",index=False)

##训练模型初始化
clf = MolTrain(task='regression',
                data_type='molecule',
                epochs=100,
                learning_rate=0.0001,
                batch_size=32,
                early_stopping=10,
                metrics='mse',
                split='random',
                save_path='./clf_fep',  #!!!
                )

####模型训练
clf.fit('train_fep_2.csv')

####基于SMILES文件输入模式的预测
clf = MolPredict(load_model='./clf_fep')  #模型保存  !!!
test_path = 'test_fep_2.csv'
test_pred = clf.predict(test_path)

####统计结果并画图
df = pd.read_csv(test_path,header='infer')
test_target = df['TARGET'].values

#########异常点开始
#####找出测试集中异常点的索引
residuals = np.abs(test_target - test_pred.flatten())
threshold = np.percentile(residuals,98)   #使用第98百分位数作为阈值
print(threshold)
noise_indices = np.where(residuals > threshold)[0]

######删除测试集中异常点
cleaned_test_data = test_data.drop(test_data.index[noise_indices])
cleaned_test_data.to_csv("test_fep_2_cleaned.csv",index=False)

#####基于删除异常点的数据进行预测
test_pred = clf.predict('test_fep_2_cleaned.csv')
test_clean_path = 'test_fep_2_cleaned.csv'

####统计删除异常点结果并画图
df = pd.read_csv(test_clean_path,header='infer')
test_target = df['TARGET'].values
##################################异常点结束

rmse_test = np.sqrt(mean_squared_error(test_target,test_pred.flatten()))
R2_test = r2_score(test_target,test_pred.flatten())


fig, ax = plt.subplots(figsize=(5,5),dpi=300)
xmin = min(test_pred.flatten().min(),test_target.min())
xmax = max(test_pred.flatten().max(),test_target.max())
ymin = xmin
ymax = xmax

#####绘制测试图
ax.scatter(test_target,test_pred.flatten(),alpha=0.2,s=10,c='red',label='Test')
# ax.scatter(test_target[noise_indices],test_pred.flatten()[noise_indices],c='blue',label='Noise Ponints')


ax.text(0.6,0.11,"RMSE (Test) = " + "%.3f"%(rmse_test),fontsize=10,transform=ax.transAxes)
ax.text(0.6,0.07,"R$^{2}$ (Test) = " + "%.3f"%(R2_test),fontsize=10,transform=ax.transAxes)

plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)

ax.set_xlabel('target')
ax.set_ylabel('predict',labelpad=0,fontsize=10)
# ax.set_title('free energy of target-predict',fontsize=14)   #!!!

_ = ax.plot([xmin,xmax],[ymin,ymax],c='k',ls='--')
ax.legend(loc='upper left')

plt.savefig('figure/clf_predict_free_energy.png',dpi=300) 




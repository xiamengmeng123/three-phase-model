import pandas as pd
from unimol_tools import MolTrain,MolPredict
from sklearn.metrics import mean_squared_error,r2_score
import matplotlib.pyplot as plt 
import numpy as np 


####数据处理
data = pd.read_csv('csv/mol_interaction_energy.csv',sep=',').iloc[:,[1,2]]#.sample(100)
print(data.head())
data.columns = ["SMILES","TARGET"]

#将数据集的80%设置为训练数据集，20%设置为测试数据集
train_fraction = 0.8
train_data = data.sample(frac=train_fraction,random_state=1)
train_data.to_csv("csv/train.csv",index=False,header=['SMILES','TARGET'])
test_data = data.drop(train_data.index)
test_data.to_csv("csv/test.csv",index=False,header=['SMILES','TARGET'])

###训练模型初始化
clf = MolTrain(task='regression',
                data_type='molecule',
                epochs=100,
                learning_rate=0.0001,
                batch_size=32,
                early_stopping=10,
                metrics='mse',
                split='random',
                save_path='./clf_interaction',
                )

####模型训练
clf.fit('csv/train.csv')

####基于SMILES文件输入模式的预测
clf = MolPredict(load_model='./clf_interaction')
test_path = 'csv/test.csv'
test_pred = clf.predict(test_path)


####统计结果并画图
df = pd.read_csv(test_path,header='infer')
test_target = df['TARGET'].values


rmse_test = np.sqrt(mean_squared_error(test_target,test_pred.flatten()))
R2_test = r2_score(test_target,test_pred.flatten())

fig, ax = plt.subplots(figsize=(5,5),dpi=300)
xmin = min(test_pred.flatten().min(),test_target.min())
xmax = max(test_pred.flatten().max(),test_target.max())
ymin = xmin
ymax = xmax

ax.scatter(test_target,test_pred.flatten(),alpha=0.2,s=10,c='red',label='Test')
# ax.scatter(test_target[noise_indices],test_pred.flatten()[noise_indices],c='blue',label='Noise Ponints')

ax.text(0.6,0.11,"RMSE (Test) = " + "%.3f"%(rmse_test),fontsize=10,transform=ax.transAxes)
ax.text(0.6,0.07,"R$^{2}$ (Test) = " + "%.3f"%(R2_test),fontsize=10,transform=ax.transAxes)

plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)

ax.set_xlabel('target')
ax.set_ylabel('predict',labelpad=0)

_ = ax.plot([xmin,xmax],[ymin,ymax],c='k',ls='--')
ax.legend(loc='upper left')

plt.savefig('figure/clf_predict_interaction.png',dpi=300)





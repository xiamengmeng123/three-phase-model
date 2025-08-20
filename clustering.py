# #######
# #输入：energy.csv,包含两列，["molecule_id","SMILES","target"]
# #输出：分子聚类图和提取簇的典型官能团
import os 
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Draw
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from rdkit.Chem import rdFMCS
# from rdkit.Chem.Draw import rdMolDraw2D
import plotly.express as px
import seaborn as sns
import random
import matplotlib.pyplot as plt 
from unimol_tools import UniMolRepr


# 单个数据依次产生
# def calculate_unimol_qsar_repr(smiles):  
#     clf = UniMolRepr(data_type='molecule', remove_hs=False)
#     smiles_list = [smiles]
#     unimol_repr = clf.get_repr(smiles_list, return_atomic_reprs=True)
#     fp = np.array(unimol_repr['cls_repr'])
#     return fp

#多个数据连续产生
def calculate_unimol_qsar_repr(data):
    clf = UniMolRepr(data_type='molecule', remove_hs=False)
    smiles_list = data["SMILES"].tolist()  #500个数据
    repr_dict = clf.get_repr(smiles_list,return_atomic_reprs=True)  # 这里的 clf 需要提前定义
    unimol_repr_list = np.array(repr_dict['cls_repr'])
    unimol_repr = repr_dict['cls_repr']
    return unimol_repr

def calculate_unimol_qsar_repr_atomic(data):   #没找到降维方法
    clf = UniMolRepr(data_type='molecule', remove_hs=False)
    smiles_list = data["SMILES"].tolist()  #500个数据
    repr_dict = clf.get_repr(smiles_list,return_atomic_reprs=True)  # 这里的 clf 需要提前定义
    unimol_repr_list = np.array(repr_dict['atomic_reprs'])
    unimol_repr = repr_dict['atomic_reprs']
    # a = unimol_repr_list[4204].shape   #形状不一样，统计固定为(120，512)
    # print(f"{a}")
    
    #填充数组到固定shape
    # target_length = 120*512
    # unimol_repr = np.zeros((len(repr_dict['atomic_reprs']),target_length))
    # print(unimol_repr.shape)
    # for i, arr in enumerate(repr_dict['atomic_reprs']):
    #     arr_length = len(arr)
    #     unimol_repr[i, :arr_length] = arr
    #unimol_repr = unimol_repr_list.reshape(unimol_repr_list.shape[0],-1)   #将三维变成二维，并第一个维度大小不变，第二个维度自动换算，总维度不变 
    return unimol_repr


def calculate_2dqsar_repr(smiles):
    mol = Chem.MolFromSmiles(smiles)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 3, nBits=1024)
    return np.array(fp)

def tsne_dimension_reduction(data):
    # X = np.array(data["unimol_qsar_mr"].values.tolist())   ##特征向量改变！
    X = np.array(data["unimol_qsar_mr"].apply(lambda x: np.ravel(x)).tolist())  
    # X = np.array(data["unimol_qsar_mr_atomic"].apply(lambda x: np.ravel(x)).tolist())  #维度不同，不能用数组,atomic无法使用
    tsne = TSNE(n_components=2, random_state=42)
    X_tsne = tsne.fit_transform(X)
    return X_tsne

def kmeans_clustering(data, n_clusters=30):
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    data["Cluster"] = kmeans.fit_predict(data[["Dimension 1", "Dimension 2"]])
    return data

def save_cluster_images(tsne_data, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for cluster_id in tsne_data['Cluster'].unique():
        cluster_data = tsne_data[tsne_data["Cluster"] == cluster_id]
        cluster_dir = os.path.join(output_dir, f"Cluster_{cluster_id}")
        os.makedirs(cluster_dir, exist_ok=True)
        for index, row in cluster_data.iterrows():
            mol = Chem.MolFromSmiles(row["SMILES"])
            Draw.MolToFile(mol, os.path.join(cluster_dir, f"{row['molecule_id']}.png"))

def add_cluster_image_path(tsne_data, output_dir):
    image_paths = []
    for index, row in tsne_data.iterrows():
        cluster_id = row["Cluster"]
        image_path = os.path.join(output_dir, f"Cluster_{cluster_id}", f"{row['molecule_id']}.png")
        image_paths.append(image_path)
    tsne_data["Image_Path"] = image_paths
    return tsne_data

def extract_common_substructure(tsne_data, n_class):
    for i in range(n_class):
        cluster_data = tsne_data[tsne_data["Cluster"] == i]
        smiles_list = cluster_data["SMILES"].tolist()
        mol_list = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]
        
        similarity_matrix = np.zeros((len(mol_list), len(mol_list)))
        common_substructure_count = np.zeros((len(mol_list), len(mol_list)))
        
        for j in range(len(mol_list)):
            for k in range(j+1, len(mol_list)):
                similarity = DataStructs.TanimotoSimilarity(AllChem.GetMorganFingerprintAsBitVect(mol_list[j], 3, nBits=1024),
                                                            AllChem.GetMorganFingerprintAsBitVect(mol_list[k], 3, nBits=1024))
                similarity_matrix[j, k] = similarity
                similarity_matrix[k, j] = similarity
                
                mcs = rdFMCS.FindMCS([mol_list[j], mol_list[k]])
                common_substructure_count[j, k] = len(mol_list[j].GetSubstructMatches(Chem.MolFromSmarts(mcs.smartsString)))
                common_substructure_count[k, j] = common_substructure_count[j, k]
        
        max_common_substructure_index = np.unravel_index(np.argmax(common_substructure_count), common_substructure_count.shape)
        max_common_substructure = Chem.MolFromSmarts(rdFMCS.FindMCS([mol_list[max_common_substructure_index[0]], mol_list[max_common_substructure_index[1]]]).smartsString)
        
        img = Draw.MolToImage(max_common_substructure, size=(300, 300))
        img.save(f"figure/cluster_unimol_cls/cluster_{i}_common_substructure.png")



# 读取原始数据
data = pd.read_csv("csv/mol_free_energy.csv", sep=',', header=None).head(500)
data.columns = ["molecule_id", "SMILES", "TARGET"]

# 提取特征指纹
# data["unimol_qsar_mr"] = data["SMILES"].apply(calculate_unimol_qsar_repr)  #单个数据提取方式
data["unimol_qsar_mr"] = calculate_unimol_qsar_repr(data)   #多个数据提取方式
data["unimol_qsar_mr_atomic"] = calculate_unimol_qsar_repr_atomic(data)
# data.to_csv("csv/unimol_cls_feature_500.csv", index=False)
# data.to_csv('csv/gfkd_mmx_data.csv',index=False,header=True)

# t-SNE降维
X_tsne = tsne_dimension_reduction(data)
tsne_data = pd.DataFrame(data=X_tsne, columns=["Dimension 1", "Dimension 2"])
tsne_data = pd.concat([data[["molecule_id", "SMILES", "TARGET"]], tsne_data], axis=1)
# tsne_data = pd.concat([data[["SMILES"]], tsne_data], axis=1)
tsne_data.replace([np.inf, -np.inf], np.nan, inplace=True)
tsne_data.dropna(inplace=True)


# K均值聚类分析
tsne_data = kmeans_clustering(tsne_data, n_clusters=30)
tsne_data.to_csv("figure/cluster_unimol_cls/clustered_data.csv", index=False)

# 保存簇的文件图
output_dir = "figure/cluster_unimol_cls/cluster_images/"
save_cluster_images(tsne_data, output_dir)

# 添加簇的文件图路径到tsne_data列
tsne_data = add_cluster_image_path(tsne_data, output_dir)

# 提取簇的典型官能团
extract_common_substructure(tsne_data, n_class=30)

#t-sne可视化配置
marker_list = ['o', 's', '^', 'D','v', 'p', '*', 'X','+','h','>','<','H','D','s','d']
classes = tsne_data["Cluster"].unique()
class_list = np.unique(tsne_data['Cluster'])
n_class = len(class_list)
palette = sns.color_palette("hls",n_class)
sns.palplot(palette)

random.seed(1234)
random.shuffle(marker_list)
random.shuffle(palette)

#可视化展示
#不同的符号表示不同的标注
plt.figure(figsize=(14,14))
for idx,cluster_id in enumerate(class_list):#遍历每个簇
    color = palette[idx]
    marker = marker_list[idx % len(marker_list)]

    #找到所有标注簇为当前簇的图像索引
    indices = tsne_data['Cluster']==cluster_id
    plt.scatter(tsne_data.loc[indices, 'Dimension 1'], tsne_data.loc[indices, 'Dimension 2'], color=color, marker=marker, label=f'{cluster_id}')
    
plt.title('KMeans Clustering on t-SNE Data')
plt.legend(fontsize=16, markerscale=1, bbox_to_anchor=(1, 1))
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small')  # 调整图例位置和字体大小
plt.xticks([])
plt.yticks([])
plt.savefig("figure/cluster_unimol_cls/cluster_plot.png",dpi=300)  # 保存图片为PNG格式


# 使用Plotly创建交互式t-SNE图
fig = px.scatter(tsne_data, x='Dimension 1', y='Dimension 2', color='Cluster',symbol='Cluster',labels='Cluster', 
                 hover_data=['molecule_id', 'TARGET', 'Image_Path'],opacity=0.8,width=1000,height=600)
fig.update_layout(margin=dict(l=0,r=0,b=0,t=0))
# 保存交互式t-SNE图html格式
fig.write_html("figure/cluster_unimol_cls/tsne_plot.html")



    

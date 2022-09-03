
from sklearn.neighbors import KDTree
from metrics import *
from keras.models import Model
from keras.layers import Dense, Input
import scipy.sparse
from scipy import sparse as sp
import networkx as nx
import glob
import itertools
import os
import numpy as np
import random
import pandas as pd
from sklearn.model_selection import train_test_split
import pickle as pkl
import scipy.sparse


#' @param num.cc Number of canonical vectors to calculate
#' @param seed.use Random seed to set.
#' @importFrom SVD
# num_cc 是最终图中每个节点的特征向量长度，
def embed(data1, data2, num_cc=20):
    random.seed(123)
    object1 = Scale(data1)
    object2 = Scale(data2)
    mat3 = np.matmul(np.matrix(object1).transpose(), np.matrix(object2))
    a = SVD(mat=mat3, num_cc=int(num_cc))
    embeds_data = np.concatenate((a[0], a[1]))
    ind = np.where(
        [embeds_data[:, col][0] < 0 for col in range(embeds_data.shape[1])])[0]
    embeds_data[:, ind] = embeds_data[:, ind] * (-1)

    embeds_data = pd.DataFrame(embeds_data)
    embeds_data.index = np.concatenate(
        (np.array(data1.columns), np.array(data2.columns)))
    embeds_data.columns = ['D_' + str(i) for i in range(num_cc)]
    d = a[2]
    #' d = np.around(a[2], 3)  #.astype('int')
    return embeds_data, d


def Embed(data_use1, data_use2, features, count_names, num_cc):

    ###check if variance of each gene are 0, if so, delete it.
    #features = checkFeature(data_use1, features)
    #features = checkFeature(data_use2, features)

    data1 = data_use1.loc[features, ]
    data2 = data_use2.loc[features, ]
    embed_results = embed(data1=data1, data2=data2, num_cc=num_cc)
    cell_embeddings = np.matrix(embed_results[0])
    combined_data = data1.merge(data2,
                                left_index=True,
                                right_index=True,
                                how='inner')
    new_data1 = combined_data.loc[count_names, ].dropna()
    # loadings=loadingDim(new.data1,cell.embeddings)
    loadings = pd.DataFrame(np.matmul(np.matrix(new_data1), cell_embeddings))
    loadings.index = new_data1.index
    return embed_results, loadings


def checkFeature(data_use, features):
    data1 = data_use.loc[features, ]
    feature_var = data1.var(1)
    Var_features = features[np.where(feature_var != 0)[0]]
    return Var_features

#dist输出n个点和所有n个点的top-k的最小距离，ind输出输出每一个n个点和所有n个点的top-k的最小距离的点的index
def kNN(data, k, query=None):
    tree = KDTree(data)
    if query is None:
        query = data
    dist, ind = tree.query(query, k)
    return dist, ind


#输出cell1和cell2中自己与自己的kNN结果、和对方kNN的结果，k表示取top-k个距离最近的点
#' @param cell_embedding : pandas data frame
def KNN(cell_embedding, cells1, cells2, k):
    embedding_cells1 = cell_embedding.loc[cells1, ]
    embedding_cells2 = cell_embedding.loc[cells2, ]
    nnaa = kNN(embedding_cells1, k=k + 1)
    nnbb = kNN(embedding_cells2, k=k + 1)
    nnab = kNN(data=embedding_cells2, k=k, query=embedding_cells1)
    nnba = kNN(data=embedding_cells1, k=k, query=embedding_cells2)
    return nnaa, nnab, nnba, nnbb, cells1, cells2

#mutual nearest neighbours，输出了两列点，每一行为两个点的index，相互连线
def MNN(neighbors, colnames, num):
    max_nn = np.array([neighbors[1][1].shape[1], neighbors[2][1].shape[1]])
    if ((num > max_nn).any()):
        num = np.min(max_nn)
        # convert cell name to neighbor index
    cells1 = colnames
    cells2 = colnames
    nn_cells1 = neighbors[4]
    nn_cells2 = neighbors[5]
    cell1_index = [
        list(nn_cells1).index(i) for i in cells1 if (nn_cells1 == i).any()
    ]
    cell2_index = [
        list(nn_cells2).index(i) for i in cells2 if (nn_cells2 == i).any()
    ]
    ncell = range(neighbors[1][1].shape[0])
    ncell = np.array(ncell)[np.in1d(ncell, cell1_index)]
    # initialize a list
    mnn_cell1 = [None] * (len(ncell) * 5)
    mnn_cell2 = [None] * (len(ncell) * 5)
    idx = -1
    for cell in ncell:
        neighbors_ab = neighbors[1][1][cell, 0:5]
        mutual_neighbors = np.where(
            neighbors[2][1][neighbors_ab, 0:5] == cell)[0]
        for i in neighbors_ab[mutual_neighbors]:
            idx = idx + 1
            mnn_cell1[idx] = cell
            mnn_cell2[idx] = i
    mnn_cell1 = mnn_cell1[0:(idx + 1)]
    mnn_cell2 = mnn_cell2[0:(idx + 1)]
    import pandas as pd
    mnns = pd.DataFrame(np.column_stack((mnn_cell1, mnn_cell2)))
    mnns.columns = ['cell1', 'cell2']
    return mnns

#intra = TRUE时，输入都为real-ST，所以输出不需要合并data1和data2，只保留一个即可
def AE_dim_reduction(data1,data2,intra = False):
    encoding_dim = 50
    feature_num = int(np.shape(data1)[0])
    diff = feature_num - encoding_dim
    input_img = Input(shape=(feature_num,))

    # encoder layers
    encoded = Dense(round(diff/2)+encoding_dim, activation='relu')(input_img)
    encoded = Dense(round(diff/4)+encoding_dim, activation='relu')(encoded)
    encoded = Dense(round(diff/8)+encoding_dim, activation='relu')(encoded)
    encoder_output = Dense(encoding_dim)(encoded)

    # decoder layers
    decoded = Dense(round(diff/8)+encoding_dim, activation='relu')(encoder_output)
    decoded = Dense(round(diff/4)+encoding_dim, activation='relu')(decoded)
    decoded = Dense(round(diff/2)+encoding_dim, activation='relu')(decoded)
    decoded = Dense(feature_num, activation='tanh')(decoded)

    # construct the autoencoder model
    autoencoder = Model(input=input_img, output=decoded)

    # construct the encoder model for plotting
    encoder = Model(input=input_img, output=encoder_output)

    # compile autoencoder
    autoencoder.compile(optimizer='adam', loss='mse')
    if intra:
        train = data1.transpose()
    else:
        train = pd.concat([data1, data2], axis=1).transpose()
    # training
    autoencoder.fit(train, train,
                    epochs=5,
                    batch_size=200,
                    shuffle=True)

    # output encoded results
    results = pd.DataFrame(encoder.predict(train).transpose())
    if intra:
        results.columns = np.array(data1.columns)
    else:
        results.columns = np.concatenate((np.array(data1.columns), np.array(data2.columns)))

    results.index = ['D_' + str(i) for i in range(np.shape(results)[0])]
    results = results.transpose()
    return results

#用embed的方式和直接使用所有feature的方式，分别得到的所有edges，做交集，即为最终filter之后的结果
def filterEdge(edges, neighbors, mats, features, k_filter):
    nn_cells1 = neighbors[4]
    nn_cells2 = neighbors[5]
    mat1 = mats.loc[features, nn_cells1].transpose()
    mat2 = mats.loc[features, nn_cells2].transpose()
    cn_data1 = l2norm(mat1)
    cn_data2 = l2norm(mat2)
    nn = kNN(data=cn_data2.loc[nn_cells2, ],
             query=cn_data1.loc[nn_cells1, ],
             k=k_filter)
    # position = [
    #     np.where(
    #         edges.loc[:, "cell2"][x] == nn[1][edges.loc[:, 'cell1'][x], ])[0]
    #     for x in range(edges.shape[0])
    # ]
    # nps = np.concatenate(position, axis=0)
    # fedge = edges.iloc[nps, ]
    index = []
    for x in range(edges.shape[0]):
        if len(np.where(edges.loc[:, "cell2"][x] == nn[1][edges.loc[:, 'cell1'][x], ])[0]) > 0:
            index.append(x)
    index = np.unique(index)
    fedge = edges.iloc[index,]
    print("\t Finally identified ", fedge.shape[0], " MNN edges")
    return (fedge)

#final_edges前两列为两组细胞，每一行是两个细胞间的一个连线
def Link_graph(count_list,
               norm_list,
               scale_list,
               features,
               combine,
               intra,
               k_filter=50 ):
    all_edges = []
    for row in combine:
        i = row[0]
        j = row[1]
        counts1 = count_list[i]
        counts2 = count_list[j]
        norm_data1 = norm_list[i]
        norm_data2 = norm_list[j]
        scale_data1 = scale_list[i]
        scale_data2 = scale_list[j]
        rowname = counts1.index.to_list()

        encoded = AE_dim_reduction(scale_data1,scale_data2,intra)
        cell_embedding, loading = Embed(data_use1=scale_data1,
                                        data_use2=scale_data2,
                                        features=features,
                                        count_names=rowname,
                                        num_cc=3)
        norm_embedding = l2norm(mat=cell_embedding[0])

        #用autoencoder降维方法可代替Embed和l2norm两个function

        cells1 = counts1.columns
        cells2 = counts2.columns
        neighbor = KNN(cell_embedding=encoded,cells1=cells1,cells2=cells2,k=50)
        #neighbor = KNN(cell_embedding=cell_embedding[0],cells1=cells1,cells2=cells2,k=5)
        mnn_edges = MNN(neighbors=neighbor,
                        colnames=cell_embedding[0].index,
                        num= 50)

        Mat = pd.concat([norm_data1, norm_data2], axis=1)
        final_edges = filterEdge(edges=mnn_edges,
                                 neighbors=neighbor,
                                 mats=Mat,
                                 features=features,
                                 k_filter=k_filter)

        all_edges.append(final_edges)
    return all_edges

def Link_Graph(outputdir,mode):
    path0 = os.path.join(os.getcwd(), outputdir)

    #' import processed data
    files1 = glob.glob(path0 + "/ST_count/*.csv")
    files1.sort()
    count_list = []
    for df in files1:
        print(df)
        count_list.append(pd.read_csv(df, index_col=0))

    files2 = glob.glob(path0 + "/ST_norm/*.csv")
    files2.sort()
    norm_list = []
    for df in files2:
        print(df)
        norm_list.append(pd.read_csv(df, index_col=0))

    files3 = glob.glob(path0 + "/ST_scale/*.csv")
    files3.sort()
    scale_list = []
    for df in files3:
        print(df)
        scale_list.append(pd.read_csv(df, index_col=0))

    files4 = glob.glob(path0 + "/ST_label/*.csv")
    files4.sort()
    label_list = []
    for df in files4:
        print(df)
        label_list.append(pd.read_csv(df, index_col=0))

    fpath = os.path.join(path0, 'Variable_features.csv')
    features = pd.read_csv(fpath, index_col=False).values.flatten()

    N = len(count_list)
    if (N == 1):
        combine = pd.Series([(0, 0)])
    else:
        combin = list(itertools.product(list(range(N)), list(range(N))))
        index = [i for i, x in enumerate([i[0] < i[1] for i in combin]) if x]
        combine = pd.Series(combin)[index]

    link1 = Link_graph(count_list=count_list,
                       norm_list=norm_list,
                       scale_list=scale_list,
                       features=features,
                       intra= False,
                       combine=combine)

    graph_pse_real = link1[0].reset_index()
    graph_pse_real.to_csv('./Datadir/Linked_graph1.csv')

    return label_list

#' data preperation
def input_data(DataDir,mode):
    label_list = Link_Graph(outputdir='Infor_Data',mode = mode)
    DataPath1 = '{}/Pseudo_ST1.csv'.format(DataDir)
    DataPath2 = '{}/Real_ST2.csv'.format(DataDir)

    #' read the data
    pse_st = pd.read_csv(DataPath1, index_col=0, sep=',')
    real_st_data = pd.read_csv(DataPath2, index_col=0, sep=',')

    pse_st_label = label_list[0]
    real_st_label = label_list[1]
    celltypes = pse_st_label.columns

    random.seed(123)

    ### if real-wrold data sets, real_st_label is all 0
    pse_train_data, pse_val_data, pse_train_label, pse_val_label = train_test_split(
        pse_st, pse_st_label, test_size=0.1, random_state=1)
    pse_test_data = pse_val_data
    pse_test_label = pse_val_label

    #' save objects

    PIK = "{}/datasets.dat".format(DataDir)
    res = [
        pse_train_data, pse_test_data, pse_val_data, pse_train_label, pse_test_label,
        pse_val_label, real_st_data, real_st_label,celltypes
    ]

    with open(PIK, "wb") as f:
        pkl.dump(res, f)

    print('load data succesfully....')


def load_data(datadir,mode):
    input_data(datadir,mode = mode)
    PIK = "{}/datasets.dat".format(datadir)
    with open(PIK, "rb") as f:
        objects = pkl.load(f)
    pse_train_data, pse_test_data, pse_val_data, pse_train_label, pse_test_label,pse_val_label, real_st_data, real_st_label,celltypes = tuple(objects)

    pl_train_data = pd.concat([pse_train_data, real_st_data])
    pl_train_label = pd.concat([pse_train_label, real_st_label],)

    pl_train_data = np.array(pl_train_data)
    pse_test_data = np.array(pse_test_data)
    pse_val_data = np.array(pse_val_data)
    pl_train_label = np.array(pl_train_label)
    pse_test_label = np.array(pse_test_label)
    pse_val_label = np.array(pse_val_label)

    #' convert pandas data frame to csr_matrix format
    pl_train_data = scipy.sparse.csr_matrix(pl_train_data.astype('Float64'))
    pse_val_data = scipy.sparse.csr_matrix(pse_val_data.astype('Float64'))
    pse_test_data = scipy.sparse.csr_matrix(pse_test_data.astype('Float64'))

    #' @param M; the number of labeled pseduoST samples in training set
    pse_train_data_len = len(pse_train_data)

    #' 4) get the feature object by combining training, test, valiation sets

    features = sp.vstack((sp.vstack((pl_train_data, pse_val_data)), pse_test_data)).tolil()
    features = preprocess_features(features)

    #' 5) Given cell type, generate three sets of labels with the same dimension


    all_labels = np.concatenate(
        [np.concatenate([pl_train_label, pse_val_label]), pse_test_label])
    all_labels = pd.DataFrame(all_labels)

    #' new label with binary values

    idx_train = range(pse_train_data_len)
    idx_pred = range(pse_train_data_len, len(pl_train_label))
    idx_val = range(len(pl_train_label), len(pl_train_label) + len(pse_val_label))
    idx_test = range(
        len(pl_train_label) + len(pse_val_label),
        len(pl_train_label) + len(pse_val_label) + len(pse_test_label))

    pse_train_mask = sample_mask(idx_train, all_labels.shape[0])
    real_train_mask = sample_mask(idx_pred, all_labels.shape[0])
    val_mask = sample_mask(idx_val, all_labels.shape[0])
    test_mask = sample_mask(idx_test, all_labels.shape[0])

    labels_binary_train = np.zeros(all_labels.shape)
    labels_binary_val = np.zeros(all_labels.shape)
    labels_binary_test = np.zeros(all_labels.shape)

    labels_binary_train[pse_train_mask, :] = all_labels.iloc[pse_train_mask, :]
    labels_binary_val[val_mask, :] = all_labels.iloc[val_mask, :]
    labels_binary_test[test_mask, :] = all_labels.iloc[test_mask, :]

    #' ----- construct adjacent matrix ---------
    #id_graph1一共三列，第一列是index，第二列是pseudo-ST的cell index，第三列是real-ST的cell index
    id_graph1 = pd.read_csv('{}/Linked_graph1.csv'.format(datadir),
                            index_col=0,
                            sep=',')

    ###for ablation(connection)
    if mode == 'pseudo':
        id_graph2 = pd.read_csv('{}/Linked_graph2.csv'.format(datadir),
                            sep=',',
                            index_col=0)
    else:
        print('load coor files!!!!')
        id_graph2 = pd.read_csv('Infor_Data/mindis_cell_indices.csv',
                                sep=',',index_col=None)

    #' --- map index ----
    pse_val_data = pd.DataFrame(pse_val_data)
    pse_test_data = pd.DataFrame(pse_test_data)

    fake1 = np.array([-1] * len(real_st_data.index))
    index_no_real = np.concatenate((pse_train_data.index, fake1, pse_val_data.index,
                             pse_test_data.index)).flatten()

    fake2 = np.array([-1] * len(pse_train_data))
    fake3 = np.array([-1] * (len(pse_val_data) + len(pse_test_data)))
    index_no_pse = np.concatenate((fake2, np.array(real_st_data.index), fake3)).flatten()

    #' ---------------------------------------------
    #'  intra-graph(id_grp1和id_grp2每组点为矩阵中的对角线对称点）
    #' ---------------------------------------------
    #cells的两列全都是real-ST data

    ###for ablation(connection)
    cells = id_graph2[['cell1','cell2']] + len(pse_train_data)
    id_grp1 = np.array([0,0])
    id_grp2 = np.array([0,0])

    for i in range(len(id_graph2)):
        id_grp1 = np.row_stack((id_grp1,[cells.iloc[i,0], cells.iloc[i,1]]))
        id_grp2 = np.row_stack((id_grp2,[cells.iloc[i,1], cells.iloc[i,0]]))


    #' ---------------------------------------------
    #'  inter-graph(id_gp1和id_gp2的每组点为矩阵中的对角线对称点）
    #' ---------------------------------------------
    #cells第一列为pseudo，第二列为real
    cells = id_graph1[['cell1','cell2']]
    cells.iloc[:,1] = cells.iloc[:,1] + len(pse_train_data)
    for i in range(len(id_graph1)):
        if id_graph1.iloc[i, 0] < len(pse_train_data):
            id_graph1.iloc[i, 0] = id_graph1.iloc[i, 0] + len(pse_train_data)
        elif (id_graph1.iloc[i, 0] >= len(pse_train_data)) and (id_graph1.iloc[i, 0] < (len(pse_train_data) + len(pse_val_data))):
            id_graph1.iloc[i, 0] = id_graph1.iloc[i, 0] + len(pse_train_data) + len(pse_val_data)
        elif id_graph1.iloc[i, 0] >= (len(pse_train_data) + len(pse_val_data) + len(pse_test_data)):
            id_graph1.iloc[i, 0] = id_graph1.iloc[i, 0] + len(pse_train_data) + len(pse_val_data) + len(pse_test_data)

    id_gp1 = np.array([0,0])
    id_gp2 = np.array([0,0])

    for i in range(len(id_graph1)):
        id_gp1 = np.row_stack((id_gp1,[cells.iloc[i,0], cells.iloc[i,1]]))
        id_gp2 = np.row_stack((id_gp2,[cells.iloc[i,1], cells.iloc[i,0]]))


    ##matrix为邻接矩阵
    matrix = np.identity(len(all_labels))


    ###for ablation(connection)
    for i in range(len(id_grp1)):
       matrix[id_grp1[i][0],id_grp1[i][1]] = 1
       matrix[id_grp2[i][0],id_grp2[i][1]] = 1

    for i in range(len(id_gp1)):
        matrix[id_gp1[i][0],id_gp1[i][1]] = 1
        matrix[id_gp2[i][0],id_gp2[i][1]] = 1
    # for i in range(len(all_labels)):
    #    matrix[i][i] = 0

    adj = graph(matrix)
    adj = nx.adjacency_matrix(nx.from_dict_of_lists(adj))

    print("assign input coordinatly....")
    return adj, features, labels_binary_train, labels_binary_val, \
           labels_binary_test, pse_train_mask, real_train_mask, val_mask, test_mask, all_labels, all_labels,celltypes

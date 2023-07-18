import sys
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
import seaborn as sns
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from scipy.cluster.hierarchy import ward, fcluster, cophenet
import scipy.cluster.hierarchy as shc
from sklearn.cluster import AgglomerativeClustering
import numpy as np
import os
import errno
import pacmap
from itertools import zip_longest
import constants
from pathlib import Path

def parse_command_line():
    n_clust = 4
    try:
        matrix_name = sys.argv[1] # matrix file
        snapshots_name = sys.argv[2] # snapshots file
        output_name = sys.argv[3] # output folder name
        
        if len(sys.argv) == 5:
            n_clust = int(sys.argv[4])

        os.mkdir(output_name)
    except IndexError:
        print("Must have at least three arguments: <distance_matrix_file> <snapshots_file> <output_folder_name> [n_cluster]")
        exit(1)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise

    data = pd.read_csv(matrix_name)
    data.drop(data.columns[[-1,]], axis=1, inplace=True) # drop evolution column
    n_trics = len(data.columns) + 1
    patients = pd.read_csv(snapshots_name)

    # remove patietns with less than 3 appointments
    patients.dropna(subset = 'ALSFRS-R', inplace = True)
    counts = patients['REF'].value_counts()
    mask = counts >= 3
    filtered_patients = patients[patients['REF'].isin(counts[mask].index)]
    filtered_patients = filtered_patients.groupby('REF').first().reset_index()

    return data, patients,filtered_patients, output_name, n_clust

def parse_config_file():
    constants.get_config(sys.argv[1])
    try:
        matrix_name = constants.MATRIX_FILE # matrix file
        snapshots_name = constants.SNAPSHOTS_FILE # snapshots file
        output_name = constants.OUTPUT_FOLDER# output folder name
        n_clust = constants.N_CLUST 
        features =  list(constants.TEMPORAL_FEATURES.keys())
        ref_feature = constants.REF_FEATURE
        os.mkdir(output_name)
    except OSError as exc:
        try:
            Path(output_name + 'trajectories/').mkdir(parents=True, exist_ok=True) 
            Path(output_name + 'visualizations/').mkdir(parents=True, exist_ok=True) 
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise

    data = pd.read_csv(matrix_name)
    data.drop(data.columns[[-1,]], axis=1, inplace=True) # drop evolution column
    n_trics = len(data.columns) + 1
    patients = pd.read_csv(snapshots_name)

    if constants.DATA_PREPROCESSING:
        #patients.dropna(subset = list(constants.MAIN_FEATURE.keys()), inplace = True)

        # remove patietns with less than MIN_APP appointments
        counts = patients[ref_feature].value_counts()
        mask = counts >= constants.MIN_APP
        filtered_patients = patients[patients[ref_feature].isin(counts[mask].index)]
        filtered_patients = filtered_patients.groupby(ref_feature).first().reset_index()
    
    return data, patients,filtered_patients, output_name, n_clust, features, ref_feature

def get_color_list(n_clust):
    colormap = plt.cm.get_cmap('rainbow')
    colors = [colormap(i/n_clust) for i in range(n_clust)]
    return colors

def tsne(tempData, labels, output_name, n_clust):
    """
    Computes the tsne dimensionality reduction 

    Parameters
    ----------
    tempData: data to plot
    labels: labels of each patient in data
    output_name: name of the output folder
    n_clust: number of clusters in the data
    """
    colors = get_color_list(n_clust)
    tsne = TSNE(n_components=2, random_state=0)
    X_2d = tsne.fit_transform(tempData)

    new = tempData.copy()
    new['tsne-2d-one'] = X_2d[:,0]
    new['tsne-2d-two'] = X_2d[:,1]
    fig = plt.figure(figsize=(16,10))
    sns.scatterplot(
        x = "tsne-2d-one", y = "tsne-2d-two",
        hue = labels,
        palette = colors,
        data = new,
        legend = "full"
)
    fig.savefig(output_name + 'visualizations/tsne.pdf')

def pacmap_func(tempData, labels,output_name, n_clust):
    """
    Computes the pacmap dimensionality reduction 

    Parameters
    ----------
    tempData: data to plot
    labels: labels of each patient in data
    output_name: name of the output folder
    n_clust: number of clusters in the data
    """
    colors = get_color_list(n_clust)
    print(colors)
    embedding = pacmap.PaCMAP(n_components=2, random_state=0)
    X_2d = embedding.fit_transform(tempData.values)

    new = tempData.copy()
    new['pacmap-2d-one'] = X_2d[:,0]
    new['pacmap-2d-two'] = X_2d[:,1]

    plt.figure(figsize=(16,10))
    sns.scatterplot(
        x = "pacmap-2d-one", y = "pacmap-2d-two",
        hue = labels,
        palette = colors,
        data = new,
        legend = "full"
    )

    plt.savefig(output_name + 'visualizations/pacmap.pdf')

def hierarchical_clustering(data,output_name, n_clust):
    """
    Computes the agglomerative clustering 

    Parameters
    ----------
    data: data to cluster
    output_name: name of the output folder
    n_clust: number of clusters in the data
    """
    plt.ylabel('distance')
    clusters = shc.linkage(data, method="ward", metric="euclidean")

    shc.dendrogram(clusters, p = 20, truncate_mode = 'lastp', # show only the last p merged clusters
                show_leaf_counts = False) 
    plt.gcf()
    plt.savefig(output_name + '/dendrogram.pdf')

    Ward_model = AgglomerativeClustering(n_clusters= n_clust, metric='euclidean', linkage='ward')
    Ward_model.fit(data)

    print('Silhouette Score: ', silhouette_score(data, Ward_model.labels_, metric='euclidean'))
    print('Calinski Harabasz Score: ', calinski_harabasz_score(data, Ward_model.labels_))
    print('Davies Bouldin Score: ', davies_bouldin_score(data, Ward_model.labels_))

    return Ward_model.labels_


def slope(x1, y1, x2, y2):
    m = (y2-y1)/(x2-x1)
    return m

def format_mogp_axs(ax, max_x=8, x_step=1.0, y_label=[0,24,48], y_minmax=(-3, 53)):
    ax.set_xlim([0, max_x])
    ax.set_xticks(np.arange(0, max_x + 1, x_step))
    ax.set_yticks(y_label)
    ax.set_ylim(y_minmax)
    return ax


def simple_trajectories(clusters, n_clust, features, ref_feature,output_name):
    """
    Computes the trajectories of each clustering in the temporal features

    Parameters
    ----------
    clusters: list of n_clust lists with each list comprising the snapshots of the patients in the corresponding cluster
    output_name: name of the output folder
    n_clust: number of clusters in the data
    """
    colors = get_color_list(n_clust)

    for feature in features:
        fig, ax = plt.subplots()

        max_val =0
        for j in range(n_clust):
            prg = clusters[j].groupby(ref_feature)[feature]           

            lst = []
            for _, group in prg:
                lst.append(group.values)
            
            transposed = [list(filter(None,i)) for i in zip_longest(*lst)]

            n_samples = []
  
            means = []
            ci = []
            max_val = 0
            for values in transposed:
                if values != []:
                    if max_val < np.nanmax(values):
                        max_val = np.nanmax(values)
                    n_samples.append(len(values))
                    ci = np.append(ci, 1.96 * np.nanstd(values)/np.sqrt(len(values)))
                    means = np.append(means,    np.nanmean(values))
            app = np.arange(0,len(means))
            #num_pat = 'n = {}'.format(len(clusters[j].groupby('REF')))
            
            ax.plot(app, means, marker = '.', color = colors[j], label = 'Cluster ' + str(j+1))
            ax.fill_between(app, (means-ci), (means+ci), color=colors[j], alpha=.1)

        
            slope_value=[]
            for i in range(1,6):
                v=slope(i-1,means[i-1], i, means[i])
                slope_value.append(v)
                plt.text( i-0.5 , (means[i] + means[i-1])/2, str(round(v,2)), fontsize=5, color = colors[j])
                
                plt.text(i, means[i] + 0.001, str(n_samples[i-1]), fontsize=8, color = colors[j], fontweight= 'bold')

        format_mogp_axs(ax, 5, 1, y_label=[0,max_val/2,max_val], y_minmax = (0, max_val+1))

        plt.xlabel("Appointments")
        plt.ylabel(str(feature))
        plt.legend()
        fig.savefig( output_name + 'trajectories/'+ str(feature) + '.pdf')
from utils import parse_command_line, hierarchical_clustering, tsne, pacmap_func,simple_trajectories,parse_config_file
import pandas as pd 

if __name__ == "__main__":

    data, patients, filtered_patients, output_name, n_clust, features, ref_feature = parse_config_file()

    #find the clusters in data
    labels = hierarchical_clustering(data, output_name, n_clust)


    #visualize the representations
    tsne(data, labels, output_name, n_clust)
    pacmap_func(data,labels, output_name, n_clust)

    df = pd.DataFrame() 
    df['Patient_ID'] = filtered_patients[ref_feature].values
    df['Labels'] = labels
    df.to_csv(output_name + '/labels.csv') 

    #plot the trajectories of each cluster
    patients = pd.merge(patients, df[['Patient_ID', 'Labels']].rename(columns={'Patient_ID':ref_feature}), on = ref_feature)
    patients = patients[patients['Labels'].notna()]
    clusters = []
    for i in range(n_clust):
        clusters.append(patients.loc[patients['Labels'] == i])
    simple_trajectories(clusters, n_clust, features,ref_feature, output_name) 
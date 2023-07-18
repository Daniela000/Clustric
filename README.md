# Clustric

Identifying disease progression patterns and groups of similar progressors is becoming relevant to the understanding of disease behaviours and to improving personalized treatments. Amyotrophic Lateral Sclerosis (ALS) is a neurodegenerative disease with patients manifesting heterogeneous temporal progressions. In this study we propose a novel approach, ClusTric, to learn comprehensive patterns from triclustering and use these patterns as features to obtain groups of patients with similar progressions. This is achieved from an agglomerative clustering process. We performed experiments using the Lisbon ALS clinic dataset containing data from patients' follow-ups. Our method showed to have relevant results regarding the flexibility of the analysis against the state-of-the-art.This method can be applyed to the 3W data of any disease, obtaining similar results.

### Triclustering

Obtain the triclustering representations in:

```
https://git.lasige.di.fc.ul.pt/dfsoares/bictric_project
```

### Agglomerative Clustering

Agglomerative clustering with Ward's linkage is applyed to the triclustering matrices produced by the previous triclustering algorithm.

Run Clustric:
```
python3 .\src\clustric.py <config_file>
```

#### Config File

```
MATRIX_FILE: "<path_to_triclustering_similarity_matrix>
SNAPSHOTS_FILE: <path_to_snapshots_file>
OUTPUT_FOLDER: <output_folder_name>
N_CLUST: <number_of_clusters>  
MIN_APP: <minimal_number_of_appointments>             
DATA_PREPROCESSING: True #remove patients with less than MIN_APP
REF_FEATURE: <feature_to_identify_each_patient>

TEMPORAL_FEATURES: <list_temporal_features>
```

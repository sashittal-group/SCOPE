import numpy as np
import pandas as pd
from itertools import combinations
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


def get_kmeans_cluster_labels(cell_fraction_F, k, num_restarts=10, verbose=False, store=True):

    X = cell_fraction_F.T.to_numpy()

    best_score = np.inf
    best_labels = None
    
    if verbose: print("Running k-Means with k=", k)
    for i in range(num_restarts):
        kmeans = KMeans(n_clusters=k, random_state=i * 100)
        labels = kmeans.fit_predict(X)
        score = silhouette_score(X, labels)
        if score < best_score:
            best_score = score
            best_labels = labels
        
    cluster_labelling = pd.DataFrame({
        'mutation': cell_fraction_F.columns.to_list(),
        'mutation_group': best_labels
    })

    if store:
        cluster_labelling.to_csv(f"../outputs/scope/laks/mutation_clusters/kmeans_k_{k}.csv")
    
    return cluster_labelling, best_score
    


def get_cluster_labels_and_scores_for_range(cell_fraction_F, k_range, num_restarts=10, verbose=False, store=True):
    scores = []
    cluster_labels = []
    
    for k in k_range:
        cluster_labelling, score = get_kmeans_cluster_labels(cell_fraction_F, k, num_restarts, verbose, store)
        scores.append(score)
        cluster_labels.append(cluster_labelling)
    
    return cluster_labels, scores


def get_cell_fraction_boundaries(cell_fraction_F, mutation_group_labels, thres=0.75):

    mutation_groups = sorted(mutation_group_labels["mutation_group"].unique())
    
    medians = np.zeros((cell_fraction_F.shape[0], len(mutation_groups)))
    upper_percentiles = np.zeros((cell_fraction_F.shape[0], len(mutation_groups)))
    lower_percentiles = np.zeros((cell_fraction_F.shape[0], len(mutation_groups)))

    for i in range(len(mutation_groups)):
        mut_grp = mutation_groups[i]
        mutations_in_group = mutation_group_labels[mutation_group_labels["mutation_group"] == mut_grp]["mutation"]

        df = cell_fraction_F[mutations_in_group]
        F_clone = df.to_numpy()

        med = np.median(F_clone, axis=1)
        upp = np.percentile(F_clone, axis=1, q=int(100*thres))
        low = np.percentile(F_clone, axis=1, q=int(100*(1-thres)))

        medians[:, i] = med
        upper_percentiles[:, i] = upp
        lower_percentiles[:, i] = low
    
    F_bar = pd.DataFrame(medians, index=cell_fraction_F.index, columns=mutation_groups)
    F_hi  = pd.DataFrame(upper_percentiles, index=cell_fraction_F.index, columns=mutation_groups)
    F_lo  = pd.DataFrame(lower_percentiles, index=cell_fraction_F.index, columns=mutation_groups)

    return F_bar, F_hi, F_lo


def filter_mutations_with_multiple_subclonality(
    F_plus,
    F_minus,
    cell_fractions_F,
    mutation_group_labels,
    threshold,
):
    clusters_with_multiple_subclonal_mutation_groups = {}

    for group in F_plus.columns:
        mask = (F_plus[group] < 1) & (F_minus[group] > 0)
        subclonal_clusters = F_plus.index[mask].tolist()

        if len(subclonal_clusters) > 1:
            clusters_with_multiple_subclonal_mutation_groups[group] = list(combinations(subclonal_clusters, 2))

    mutation_group_labels = mutation_group_labels.copy()

    new_clone_id = F_plus.shape[1]
    for group, cluster_pairs in clusters_with_multiple_subclonal_mutation_groups.items():

        mutations_in_group = mutation_group_labels.loc[
            mutation_group_labels["mutation_group"] == group, "mutation"
        ]

        if len(mutations_in_group) == 0:
            continue

        for cn_a, cn_b in cluster_pairs:
            F_sub = cell_fractions_F.loc[[cn_a, cn_b], mutations_in_group]
            vals = F_sub.to_numpy().T
            filter_mask = (
                (vals[:, 0] > threshold) & (vals[:, 0] < 1 - threshold) &
                (vals[:, 1] > threshold) & (vals[:, 1] < 1 - threshold)
            )

            mutations_to_update = mutations_in_group[filter_mask]
            mutation_group_labels.loc[
                mutation_group_labels["mutation"].isin(mutations_to_update),
                "mutation_group"
            ] = new_clone_id
            new_clone_id += 1
    
    return mutation_group_labels

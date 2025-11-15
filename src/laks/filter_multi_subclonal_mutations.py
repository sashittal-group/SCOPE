import pandas as pd
import numpy as np
import argparse
import sys

def get_mutation_group_cn_clusters_to_filter(F_plus, F_minus):
    cn_clusters = F_plus.index
    mutations_groups = F_plus.columns
    filtered_mutation_groups = {}

    for mutation_group in mutations_groups:
        subclonal = []
        for cn_cluster in cn_clusters:
            if F_plus.loc[cn_cluster, mutation_group] < 1 and F_minus.loc[cn_cluster, mutation_group] > 0:
                subclonal.append(cn_cluster)
        
        if len(subclonal) > 1:
            mutation_group_pairs = []
            for i in range(len(subclonal)):
                for j in range(i+1, len(subclonal)):
                    mutation_group_pairs.append((subclonal[i], subclonal[j]))
            
            filtered_mutation_groups[mutation_group] = mutation_group_pairs
        
    return filtered_mutation_groups



def filter(F_prime, mutation_group_to_recluster, mutation_group_labels, cn_a, cn_b,thres, recluster_id):

    mutations_of_int = mutation_group_labels[mutation_group_labels['clone'] == mutation_group_to_recluster]['mutation']

    F_int = F_prime.loc[[cn_a, cn_b], mutations_of_int]
    X_int = F_int.to_numpy().T
    mask = (X_int[:, 0] > thres) & (X_int[:, 0] < 1 - thres) & (X_int[:, 1] > thres) & (X_int[:, 1] < 1 - thres)
    mask_df = pd.DataFrame(mask, index=F_int.T.index)
    mask_df.index.name = 'mutation'
    mask_df.columns = ['seg']
    mutation_group_labels = pd.merge(mutation_group_labels, mask_df, on='mutation', how='left')
    df = mutation_group_labels
    print("UPDATED", len(df[df["seg"] == True]))
    df["clone"] = np.where(df["seg"] == True, recluster_id, df["clone"])

    return mutation_group_labels[['mutation', 'clone']]

# def filter():
#     F_bar = pd.read_csv(f"../data/laks/scope/F_bar.csv", index_col=0)
#     kmeans_clones = pd.read_csv(f"../data/laks/scope/kmeans_labels.csv", index_col=0)
#     F = pd.read_csv("../data/laks/scope/F.csv", index_col=0)

# mutation_groups_to_recluster = [0, 3, 5, 7, 8, 10]

if __name__ == "__main__":
    # parser = argparse.ArgumentParser()
    # parser.add_argument('--F_prime_filepath', type=str, help='F_prime filepath')
    # parser.add_argument('--clustering_labels_filepath', type=str, help='clustering labels filepath')
    # parser.add_argument('-o', type=str, help='output filepath')
    # args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    
    # main(args)

    for k in range(7, 26):
        for threshold in (0.60, 0.65, 0.70, 0.72, 0.75, 0.80):
            for filter_threshold in [0.1, 0.2, 0.25, 0.3, 0.4]:
    # for k in range(12, 13):
    #     for threshold in [0.75]:
            #  for filter_threshold in [0.2]:
                F_prime_filepath = f"data/laks/scope/F.csv"
                F_plus_filepath = f"data/laks/scope/cell_fractions/k_{k}_t_{threshold}/F_plus.csv"
                F_minus_filepath = f"data/laks/scope/cell_fractions/k_{k}_t_{threshold}/F_minus.csv"
                kmeans_labels_filepath = f"data/laks/scope/kmeans_labels/k_{k}.csv"
                outpath = f"data/laks/scope/kmeans_labels/filtered/k_{k}_t_{threshold}_f_{filter_threshold}.csv"
                
                F_prime = pd.read_csv(F_prime_filepath, index_col=0)
                F_plus = pd.read_csv(F_plus_filepath, index_col=0)
                F_plus.columns = F_plus.columns.astype(int)
                F_minus = pd.read_csv(F_minus_filepath, index_col=0)
                F_minus.columns = F_minus.columns.astype(int)
                kmeans_labels = pd.read_csv(kmeans_labels_filepath, index_col=0)
                
                filtered_mutation_groups = get_mutation_group_cn_clusters_to_filter(F_plus, F_minus)
                
                print(filtered_mutation_groups)

                updated_labels = kmeans_labels
                recluster_id = F_plus.shape[1]

                for mutation_group, cn_cluster_pairs in filtered_mutation_groups.items():
                    for cn_a, cn_b in cn_cluster_pairs:
                        updated_labels = filter(F_prime=F_prime, mutation_group_to_recluster=mutation_group,
                                                mutation_group_labels=updated_labels, cn_a=cn_a, cn_b=cn_b,
                                                thres=filter_threshold, recluster_id=recluster_id)
                    
                    recluster_id += 1

                updated_labels.to_csv(outpath)

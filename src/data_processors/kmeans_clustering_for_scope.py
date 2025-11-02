import argparse
import sys
import pandas as pd
import numpy as np
import os
from sklearn.cluster import KMeans
import numpy as np
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt


def get_F_boundaries(F, mutation_group_mapping, thres=0.75):

    mutation_groups = sorted(mutation_group_mapping["mutation_group"].unique())
    
    medians = np.zeros((F.shape[0], len(mutation_groups)))
    upper_percentiles = np.zeros((F.shape[0], len(mutation_groups)))
    lower_percentiles = np.zeros((F.shape[0], len(mutation_groups)))

    for i in range(len(mutation_groups)):
        mut_grp = mutation_groups[i]
        mutations_in_group = mutation_group_mapping[mutation_group_mapping["mutation_group"] == mut_grp]["mutation"]

        df = F.loc[:, mutations_in_group]
        F_clone = df.to_numpy()

        med = np.median(F_clone, axis=1)
        upp = np.percentile(F_clone, axis=1, q=int(100*thres))
        low = np.percentile(F_clone, axis=1, q=int(100*(1-thres)))

        medians[:, i] = med
        upper_percentiles[:, i] = upp
        lower_percentiles[:, i] = low
    
    F_bar = pd.DataFrame(medians, index=F.index, columns=mutation_groups)
    F_hi  = pd.DataFrame(upper_percentiles, index=F.index, columns=mutation_groups)
    F_lo  = pd.DataFrame(lower_percentiles, index=F.index, columns=mutation_groups)

    return F_bar, F_hi, F_lo


def main(args):

    input_prefix = args.i
    output_root_prefix = args.o
    start_k = args.start_k
    end_k = args.end_k
    runs_per_k = args.runs_per_k
    
    input_folder = input_prefix.split('/')[-2]
    output_prefix = output_root_prefix + input_folder
    os.makedirs(output_prefix, exist_ok=True)

    df_total = pd.read_parquet(f"{input_prefix}_read_count.parquet")
    df_variant = pd.read_parquet(f"{input_prefix}_variant_count.parquet")
    df = pd.read_parquet(f"{input_prefix}_character_matrix_without_noise.parquet")
    df_copy_number = pd.read_parquet(f"{input_prefix}_copy_numbers.parquet")
    
    df_cluster = df[['cluster_id']]

    df_total = df_total.join(df_cluster, how='left')
    df_variant = df_variant.join(df_cluster, how='left')
    df_copy_number = df_copy_number.join(df_cluster, how='left')

    total_table = df_total.groupby('cluster_id').sum(numeric_only=True)
    alt_table = df_variant.groupby('cluster_id').sum(numeric_only=True)
    df_cn = df_copy_number.groupby('cluster_id').mean().astype(int)

    vaf = alt_table.div(total_table).replace(np.nan, 0)
    F = (vaf * df_cn).clip(upper=1)

    ranges = range(start_k, end_k + 1, 1)
    
    X = F.T.to_numpy()

    silhouette_scores = []
    all_labels = []

    for k in ranges:
        # print(k)
        for i in range(runs_per_k):
            kmeans = KMeans(n_clusters=k, random_state=i * 100)
            labels = kmeans.fit_predict(X)
            score = silhouette_score(X, labels)
            silhouette_scores.append(score)
            all_labels.append(labels)

    sillhouettes_scores_per_k = np.array(silhouette_scores).reshape((len(ranges), runs_per_k))
    best_silhouettes_for_k = np.min(sillhouettes_scores_per_k, axis=1)

    plt.figure(figsize=(4, 3))
    plt.plot(ranges, best_silhouettes_for_k, 'bo-')
    plt.xlabel('Number of clusters (k)')
    plt.ylabel('Silhouette Score')
    plt.ylim((0.15, 0.65))
    plt.savefig(f"{output_prefix}/silhouette_scores.svg")
    plt.close()

    ranges = np.array(ranges)
    best_k = ranges[np.argmin(best_silhouettes_for_k)]

    index_of_best_k = int(np.where(ranges == best_k)[0][0])
    silhouette_scores_of_best_k = sillhouettes_scores_per_k[index_of_best_k, :]
    best_label_index_of_best_best = np.argmin(silhouette_scores_of_best_k)
    labels_per_k = np.array(all_labels).reshape((len(ranges), runs_per_k, -1))
    best_labels = labels_per_k[index_of_best_k, best_label_index_of_best_best]

    kmeans_labels = pd.DataFrame({
        'mutation': F.columns.to_list(),
        'mutation_group': best_labels
    })
    
    F_bar, F_hi, F_lo = get_F_boundaries(F, kmeans_labels, 0.75)

    mutation_group_sizes = kmeans_labels.groupby('mutation_group').size()

    F_hi.to_csv(f"{output_prefix}/F_plus.csv")
    F_lo.to_csv(f"{output_prefix}/F_minus.csv")
    F_bar.to_csv(f"{output_prefix}/F_bar.csv")
    mutation_group_sizes.to_csv(f"{output_prefix}/clone_sizes.csv")


    clones = sorted(kmeans_labels["mutation_group"].unique())
    n_clones = len(clones)

    # Set up subplot grid: 2 columns
    n_cols = 4
    n_rows = (n_clones + 5) // n_cols

    fig, axes = plt.subplots(
        nrows=n_rows, ncols=n_cols,
        figsize=(12, 4 * n_rows),
        sharex=True
    )

    axes = axes.flatten()

    # Define explicit CN cluster order
    cn_order = sorted(F.index.unique())[::-1]  # or use a fixed list if known: ['A','B','C','D','E','F','G','H','I']

    for i, clone in enumerate(clones):
        ax = axes[i]

        mutations_in_clone = kmeans_labels[kmeans_labels["mutation_group"] == clone]["mutation"].to_list()

        df = F.loc[:, mutations_in_clone]
        df = df.loc[cn_order]

        ax.boxplot(df.T.values, vert=False, labels=df.index)

        ax.set_ylabel("CN Clusters")
        ax.set_xlabel("CF")
        ax.set_title(f"{df.shape[1]} mutations from\nkmeans clustering clone {clone}")

    # Hide any unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.savefig(f"{output_prefix}/clone_cluster_cell_fractions.svg")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='input prefix', default='sample')
    parser.add_argument('-o', type=str, help='output prefix', default = 'scope_input')
    parser.add_argument('--start_k', type=int, help='start k', default = 5)
    parser.add_argument('--end_k', type=int, help='end k', default = 20)
    parser.add_argument('--runs_per_k', type=int, help='runs_per k', default = 5)
    parser.add_argument('-v', action='store_true', default=False)
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    
    main(args)

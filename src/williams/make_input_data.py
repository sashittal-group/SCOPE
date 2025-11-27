import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import KMeans
import os
import sys
import argparse
from sklearn.metrics import silhouette_score

from scope.process_data import get_cell_fraction_boundaries


def make_input_for_sample(patient_id):
    print("PATIENT_ID", patient_id)

    df_snv_counts = pd.read_csv(f"data/williams/scratch/{patient_id}/snv_counts.csv")
    df_cn_clones = pd.read_csv(f"data/williams/scratch/{patient_id}/hscn_clones.csv")

    df_clones = pd.read_csv(f"data/williams/signatures_dataset/clones_trees/{patient_id}_clones.tsv", sep='\t')

    df = pd.merge(df_snv_counts, df_clones, how='left', on='cell_id')

    df['mutation'] = df['chr'].astype(str) + ":" + df["start"].astype(str) + ":" + df["ref"] + ":" + df["alt"]

    total_table = df.pivot_table(
            index='clone_id',
            columns='mutation',
            values='total_counts',
            aggfunc='sum',
            fill_value=0
        )
    
    alt_table = df.pivot_table(
            index='clone_id',
            columns='mutation',
            values='alt_counts',
            aggfunc='sum',
            fill_value=0
        )
    
    cn_table = df.pivot_table(
            index='clone_id',
            columns='mutation',
            values='state',
            aggfunc='mean',
            fill_value=0
        ).astype(int)
    
    vaf_table = alt_table.div(total_table).fillna(0)
    F = (vaf_table * cn_table).clip(upper=1)

    os.makedirs(f"outputs/scope/williams/{patient_id}", exist_ok=True)

    ranges = range(10, 21, 1)
    runs_per_k = 5

    X = F.T.to_numpy()

    silhouette_scores = []
    all_labels = []


    for k in ranges:
        print(k)
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
    plt.title('Silhouette Method')
    plt.savefig(f"outputs/scope/williams/{patient_id}/silhouette_scores.svg")
    plt.close()

    ranges = np.array(ranges)

    for best_k in ranges:

    # best_k = ranges[np.argmin(best_silhouettes_for_k)]

        index_of_best_k = int(np.where(ranges == best_k)[0][0])
        silhouette_scores_of_best_k = sillhouettes_scores_per_k[index_of_best_k, :]
        best_label_index_of_best_best = np.argmin(silhouette_scores_of_best_k)
        labels_per_k = np.array(all_labels).reshape((len(ranges), runs_per_k, -1))
        best_labels = labels_per_k[index_of_best_k, best_label_index_of_best_best]

        kmeans_labels = pd.DataFrame({
            'mutation': F.columns.to_list(),
            'mutation_group': best_labels
        })
        os.makedirs(f"outputs/scope/williams/{patient_id}/mutation_clusters/k_{best_k}", exist_ok=True)
        kmeans_labels.to_csv(f"outputs/scope/williams/{patient_id}/mutation_clusters/k_{best_k}/kmeans_labels.csv", index=False)

        clusters = sorted(kmeans_labels["mutation_group"].unique())
        n_clusters = len(clusters)

        n_cols = 4
        n_rows = (n_clusters + 5) // n_cols

        fig, axes = plt.subplots(
            nrows=n_rows, ncols=n_cols,
            figsize=(12, 4 * n_rows),
            sharex=True
        )

        axes = axes.flatten()

        cn_order = sorted(F.index.unique())[::-1]  # or use a fixed list if known: ['A','B','C','D','E','F','G','H','I']
        print(cn_order)

        for i, cluster in enumerate(clusters):
            ax = axes[i]

            mutations_in_cluster = kmeans_labels[
                kmeans_labels["mutation_group"] == cluster
            ]["mutation"]

            df = F[mutations_in_cluster]

            # Ensure consistent y-axis order (A→I)
            df = df.loc[cn_order]
            
            # Horizontal boxplot per clone, sorted A→I
            ax.boxplot(df.T.values, vert=False, tick_labels=df.index)

            ax.set_ylabel("CN Clusters")
            ax.set_xlabel("CF")
            ax.set_title(f"{df.shape[1]} mutations from\nk-means clone {cluster}")

        # Hide any unused subplots
        for j in range(i + 1, len(axes)):
            fig.delaxes(axes[j])

        plt.tight_layout()
        plt.savefig(f"outputs/scope/williams/{patient_id}/mutation_clusters/k_{best_k}//cell_fractions.svg")
        plt.close()
        F_bar, F_hi, F_lo = get_cell_fraction_boundaries(F, kmeans_labels, 0.75)

        F_hi.to_csv(f"outputs/scope/williams/{patient_id}/mutation_clusters/k_{best_k}/F_plus.csv")
        F_lo.to_csv(f"outputs/scope/williams/{patient_id}/mutation_clusters/k_{best_k}/F_minus.csv")

    df_cn_clones_patient = pd.read_csv(f"data/williams/scratch/{patient_id}/hscn_clones.csv")

    df = df_cn_clones_patient

    merged = df.merge(
        df,
        on=['chr', 'start', 'end'],
        suffixes=('_1', '_2')
    )
    merged = merged[merged['clone_id_1'] != merged['clone_id_2']]

    cn_clusters = sorted(df[~df['clone_id'].isna()]['clone_id'].unique())

    allele = 'Maj'

    num_cn_clusters = len(cn_clusters)

    conficts = np.zeros((num_cn_clusters, num_cn_clusters), dtype=int)

    for i in range(num_cn_clusters):
        for j in range(num_cn_clusters):
            clone_i = cn_clusters[i]
            clone_j = cn_clusters[j]
            num_conflicts = len(merged[(merged['clone_id_1'] == clone_i) & (merged['clone_id_2'] == clone_j) & (merged[f'{allele}_1'] > 0) & (merged[f'{allele}_2'] == 0)])
            
            conficts[i][j] = num_conflicts

    major_conficts_df = pd.DataFrame(conficts, index=cn_clusters, columns=cn_clusters)

    allele = 'Min'

    num_cn_clusters = len(cn_clusters)

    conficts = np.zeros((num_cn_clusters, num_cn_clusters), dtype=int)

    for i in range(num_cn_clusters):
        for j in range(num_cn_clusters):
            clone_i = cn_clusters[i]
            clone_j = cn_clusters[j]
            num_conflicts = len(merged[(merged['clone_id_1'] == clone_i) & (merged['clone_id_2'] == clone_j) & (merged[f'{allele}_1'] > 0) & (merged[f'{allele}_2'] == 0)])
            
            conficts[i][j] = num_conflicts

    minor_conficts_df = pd.DataFrame(conficts, index=cn_clusters, columns=cn_clusters)
    conflicts_df = major_conficts_df + minor_conficts_df

    conflicts_df.to_csv(f"outputs/scope/williams/{patient_id}/loh_conflicts.csv")
    print(patient_id)
    print(conflicts_df)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', type=str, help='sample', default='sample')
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    sample = args.sample
    
    make_input_for_sample(sample)

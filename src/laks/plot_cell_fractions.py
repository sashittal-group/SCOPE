import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys
import pandas as pd

def plot_cell_fractions(F_prime, clustering_labels, output_filepath, clone_column_name="clone", method_name="k-means"):

    clusters = sorted(clustering_labels[clone_column_name].unique())
    n_clusters = len(clusters)

    n_cols = 4
    n_rows = (n_clusters + 5) // n_cols

    fig, axes = plt.subplots(
        nrows=n_rows, ncols=n_cols,
        figsize=(12, 4 * n_rows),
        sharex=True
    )

    axes = axes.flatten()

    cn_order = sorted(F_prime.index.unique())[::-1]  # or use a fixed list if known: ['A','B','C','D','E','F','G','H','I']
    print(cn_order)

    for i, cluster in enumerate(clusters):
        ax = axes[i]
        mutations_in_cluster = clustering_labels[clustering_labels[clone_column_name] == cluster]["mutation"]
        df = F_prime[mutations_in_cluster]
        df = df.loc[cn_order]
        F_clust = df.to_numpy()

        thres = 0.75
        median = np.median(F_clust, axis=1)
        upper_percentiles = np.percentile(F_clust, axis=1, q=int(100 * thres))
        lower_percentiles = np.percentile(F_clust, axis=1, q=int(100 * (1 - thres)))

        # Horizontal boxplot per clone, sorted A→I
        ax.boxplot(df.T.values, vert=False, tick_labels=df.index)

        ax.set_ylabel("CN Clusters")
        ax.set_xlabel("CF")
        ax.set_title(f"{df.shape[1]} mutations from\n{method_name} clone {cluster}")

    # Hide any unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.savefig(output_filepath)


def main(args):
    F_prime_filepath = args.F_prime_filepath
    clustering_labels_filepath = args.clustering_labels_filepath
    output_filepath = args.o

    F_prime = pd.read_csv(F_prime_filepath, index_col=0)
    clustering_labels = pd.read_csv(clustering_labels_filepath, index_col=0)

    plot_cell_fractions(F_prime, clustering_labels, output_filepath)


def main1():

    for k in range(7, 26):
        F_prime_filepath = f"data/laks/scope/F.csv"
        clustering_labels_filepath = f"data/laks/scope/kmeans_labels/k_{k}.csv"
        output_filepath = f"data/laks/scope/cell_fraction_figures/k_{k}.pdf"

        F_prime = pd.read_csv(F_prime_filepath, index_col=0)
        clustering_labels = pd.read_csv(clustering_labels_filepath, index_col=0)

        plot_cell_fractions(F_prime, clustering_labels, output_filepath)


def main2():
    for k in range(7, 26):
        for threshold in (0.60, 0.65, 0.70, 0.72, 0.75, 0.80):
            for filter_threshold in [0.1, 0.2, 0.25, 0.3, 0.4]:
                print(k, threshold, filter_threshold)
                F_prime_filepath = f"data/laks/scope/F.csv"
                clustering_labels_filepath = f"data/laks/scope/kmeans_labels/filtered/k_{k}_t_{threshold}_f_{filter_threshold}.csv"
                output_filepath = f"data/laks/scope/cell_fraction_figures/filtered/k_{k}_t_{threshold}_f_{filter_threshold}.pdf"

                F_prime = pd.read_csv(F_prime_filepath, index_col=0)
                clustering_labels = pd.read_csv(clustering_labels_filepath, index_col=0)

                plot_cell_fractions(F_prime, clustering_labels, output_filepath)




if __name__ == "__main__":
    # parser = argparse.ArgumentParser()
    # parser.add_argument('--F_prime_filepath', type=str, help='F_prime filepath')
    # parser.add_argument('--clustering_labels_filepath', type=str, help='clustering labels filepath')
    # parser.add_argument('-o', type=str, help='output filepath')
    # args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    
    # main(args)

    main2()



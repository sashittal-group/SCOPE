
import numpy as np
import pandas as pd
import argparse
import sys
import os


def get_F_boundaries(F_prime, mutation_group_mapping, thres=0.75, clone_column_name="clone"):

    mutation_groups = sorted(mutation_group_mapping[clone_column_name].unique())
    
    medians = np.zeros((F_prime.shape[0], len(mutation_groups)))
    upper_percentiles = np.zeros((F_prime.shape[0], len(mutation_groups)))
    lower_percentiles = np.zeros((F_prime.shape[0], len(mutation_groups)))

    for i in range(len(mutation_groups)):
        mut_grp = mutation_groups[i]
        mutations_in_group = mutation_group_mapping[mutation_group_mapping[clone_column_name] == mut_grp]["mutation"]

        df = F_prime[mutations_in_group]
        F_clone = df.to_numpy()

        med = np.median(F_clone, axis=1)
        upp = np.percentile(F_clone, axis=1, q=int(100*thres))
        low = np.percentile(F_clone, axis=1, q=int(100*(1-thres)))

        medians[:, i] = med
        upper_percentiles[:, i] = upp
        lower_percentiles[:, i] = low
    
    F_bar = pd.DataFrame(medians, index=F_prime.index, columns=mutation_groups)
    F_plus  = pd.DataFrame(upper_percentiles, index=F_prime.index, columns=mutation_groups)
    F_minus  = pd.DataFrame(lower_percentiles, index=F_prime.index, columns=mutation_groups)

    return F_bar, F_plus, F_minus


def main(args):
    F_prime_filepath = args.F_prime_filepath
    clustering_labels_filepath = args.clustering_labels_filepath
    threshold = args.t
    output_dir = args.output_dir

    write_Fs_to_file(F_prime_filepath, clustering_labels_filepath, threshold, output_dir)

    

def write_Fs_to_file(F_prime_filepath, clustering_labels_filepath, threshold, output_dir):
    F_prime = pd.read_csv(F_prime_filepath, index_col=0)
    clustering_labels = pd.read_csv(clustering_labels_filepath, index_col=0)

    F_bar, F_plus, F_minus = get_F_boundaries(F_prime, clustering_labels, thres=threshold)

    os.makedirs(output_dir, exist_ok=True)

    F_bar.to_csv(f"{output_dir}/F_bar.csv")
    F_plus.to_csv(f"{output_dir}/F_plus.csv")
    F_minus.to_csv(f"{output_dir}/F_minus.csv")


def main2():
    for k in range(7, 26):
        for threshold in (0.60, 0.65, 0.70, 0.72, 0.75, 0.80):
            F_prime_filepath = f"data/laks/scope/F.csv"
            clustering_labels_filepath = f"data/laks/scope/kmeans_labels/k_{k}.csv"
            output_dir = f"data/laks/scope/cell_fractions/k_{k}_t_{threshold}"

            write_Fs_to_file(F_prime_filepath, clustering_labels_filepath, threshold, output_dir)


def main3():
    for k in range(7, 26):
        for threshold in (0.60, 0.65, 0.70, 0.72, 0.75, 0.80):
            for filter_threshold in [0.1, 0.2, 0.25, 0.3, 0.4]:
                print(k, threshold, filter_threshold)
                F_prime_filepath = f"data/laks/scope/F.csv"
                clustering_labels_filepath = f"data/laks/scope/kmeans_labels/filtered/k_{k}_t_{threshold}_f_{filter_threshold}.csv"
                output_dir = f"data/laks/scope/cell_fractions/filtered/k_{k}_t_{threshold}_f_{filter_threshold}"

                write_Fs_to_file(F_prime_filepath, clustering_labels_filepath, threshold, output_dir)


if __name__ == "__main__":
    # parser = argparse.ArgumentParser()
    # parser.add_argument('--F_prime_filepath', type=str, help='F_prime filepath')
    # parser.add_argument('--clustering_labels_filepath', type=str, help='clustering labels filepath')
    # parser.add_argument('-o', type=str, help='output directory')
    # parser.add_argument('-t', type=float, help='threshold')
    # parser.add_argument('-output-postfix', type=str, help='output-postfix')
    # args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    
    # main(args)

    main3()

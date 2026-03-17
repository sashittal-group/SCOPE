import argparse
import sys
import pandas as pd
import numpy as np
import os
from sklearn.cluster import KMeans
import numpy as np
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
import os


# TODO: Use arguments to get input and ouput files instead of using a new script
def main(args):

    STR = args.s

    OUT_DIR = args.o
    INP_DIR = args.i
    POST_DIR = args.p
    GROUND_TRUTH_DIR = args.g

    F_bar = pd.read_csv(f"data/simulation/{INP_DIR}/{STR}/F_bar.csv", index_col=0)
    kmeans_clones = pd.read_csv(f"data/simulation/{INP_DIR}/{STR}/kmeans_clones.csv", index_col=0)

    X = pd.read_csv(f"data/simulation/{OUT_DIR}/{STR}/solution_0/X.csv", index_col=0)
    mutation_groups_selected = X[X['0'] > 0.5].index.tolist()
    mutation_groups_not_selected = X[X['0'] < 0.5].index.tolist()

    if len(mutation_groups_not_selected) > 0:

        input_prefix = f"data/simulation/{GROUND_TRUTH_DIR}/{STR}/sim"

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

        kmeans_clones_X = pd.merge(kmeans_clones, X, left_on="mutation_group", right_index=True, how="left")

        mutations_not_taken = kmeans_clones_X[kmeans_clones_X['0'] < 0.5]

        F_not = F.loc[:, mutations_not_taken['mutation']]

        F_bar_t = F_bar.iloc[:, mutation_groups_selected]

        df1 = F_not
        df2 = F_bar_t

        df1_values = df1.values.T  # shape: (num_cols_df1, num_rows)
        df2_values = df2.values.T  # shape: (num_cols_df2, num_rows)

        # compute pairwise euclidean distances between all columns
        dist_matrix = cdist(df1_values, df2_values, metric='euclidean')

        # make it a DataFrame
        dist_df = pd.DataFrame(
            dist_matrix,
            index=df1.columns,   # columns of first df
            columns=df2.columns  # columns of second df
        )

        min_col_per_row = dist_df.idxmin(axis=1)
        min_col_df = min_col_per_row.rename('min_column').reset_index()

        min_col_df['mutation'] = min_col_df['index'].str[1:].astype(int)

        mutation_group_update = pd.merge(kmeans_clones_X, min_col_df[['index', 'min_column']], left_on='mutation', right_on='index', how='left')
        mutation_group_update['mutation_group'] = mutation_group_update['min_column'].combine_first(
            mutation_group_update['mutation_group']
        )

        mutation_group_update['mutation_group'] = mutation_group_update['mutation_group'].astype(str).str.strip().astype(int)

        mutation_group_update = mutation_group_update[['mutation', 'mutation_group']]
        os.makedirs(f"data/simulation/{POST_DIR}/{STR}", exist_ok=True)
        mutation_group_update.to_csv(f"data/simulation/{POST_DIR}/{STR}/kmeans_cleaned_clones.csv")
    
    else:
        kmeans_clones.to_csv(f"data/simulation/{POST_DIR}/{STR}/kmeans_cleaned_clones.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', type=str, help='STR', default='sample')
    parser.add_argument('-i', type=str, help='INP_DIR', default='scope_input_kmeans')
    parser.add_argument('-o', type=str, help='OUT_DIR', default='scope_output_kmeans')
    parser.add_argument('-p', type=str, help='POST_DIR', default='scope_post_kmeans')
    parser.add_argument('-g', type=str, help='GROUND_TRUTH', default='ground_truth')
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    
    main(args)
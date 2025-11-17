import pandas as pd
from scipy.spatial.distance import cdist
import os

def recluster(mutation_groups_X, F_bar, mutations_not_selected, mutation_groups_selected, F_prime):
    
    F_not = F_prime.loc[:, mutations_not_selected]
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

    print(dist_df)

    min_col_per_row = dist_df.idxmin(axis=1)
    min_col_df = min_col_per_row.rename('min_column').reset_index()
    min_col_df['mutation'] = min_col_df['index']

    print(min_col_df)
    
    mutation_group_update = pd.merge(mutation_groups_X, min_col_df[['index', 'min_column']], left_on='mutation', right_on='index', how='left')
    mutation_group_update['clone'] = mutation_group_update['min_column'].combine_first(
        mutation_group_update['clone']
    )

    mutation_group_update['clone'] = mutation_group_update['clone'].astype(str).str.strip().astype(int)
    mutation_group_update = mutation_group_update[['mutation', 'clone']]

    return mutation_group_update




def recluster_io(k, t, f):
    X = pd.read_csv(f"data/laks/scope/outputs/filtered/k_{k}_t_{t}_f_{f}/X/solution_0.csv", index_col=0)
    mutation_groups_selected = X[X['0'] > 0.5].index.tolist()

    print(mutation_groups_selected)

    mutation_groups = pd.read_csv(f"data/laks/scope/kmeans_labels/filtered/k_{k}_t_{t}_f_{f}.csv", index_col=0)

    mutation_groups_X = pd.merge(mutation_groups, X, left_on="clone", right_index=True, how="left")
    mutations_not_selected = mutation_groups_X[mutation_groups_X['0'] < 0.5]['mutation'].to_list()

    print(len(mutations_not_selected))

    F_prime = pd.read_csv("data/laks/scope/F.csv", index_col=0)
    F_bar = pd.read_csv(f"data/laks/scope/cell_fractions/filtered/k_{k}_t_{t}_f_{f}/F_bar.csv", index_col=0)

    mutation_group_updated = recluster(mutation_groups_X, F_bar, mutations_not_selected, mutation_groups_selected, F_prime)

    mutation_group_updated.to_csv(f"data/laks/scope/reclustered_labels/filtered/k_{k}_t_{t}_f_{f}.csv")


if __name__ == "__main__":
    # for k in range(7, 26):
    #     for threshold in (0.60, 0.65, 0.70, 0.72, 0.75):
    #         for filter_threshold in [0.1, 0.2, 0.25, 0.3, 0.4]:
    #             try:
    #                 print(k, threshold, filter_threshold)
    #                 recluster_io(k, threshold, filter_threshold)
    #             except Exception as e:
    #                 print(e)

    k = 13
    t = 0.7
    f = 0.2
    recluster_io(k, t, f)

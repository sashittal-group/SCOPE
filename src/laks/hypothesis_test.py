import numpy as np
import pandas as pd
import os
from scipy.stats import betabinom, combine_pvalues


min_copy_number = 0
max_copy_number = 11

ROOT_DIR = "scratch/data"
def get_ground_truth():
    df_cell_cn = pd.read_csv(f"{ROOT_DIR}/ov2295_cell_cn.csv.gz")
    binsize = int(df_cell_cn['end'][0])
    df_cell_snv = pd.read_csv(f"{ROOT_DIR}/ov2295_snv_counts.csv.gz", low_memory=False)
    df = df_cell_snv
    df["bin_start"] = df["coord"] // binsize * binsize + 1
    df = pd.merge(df, df_cell_cn, left_on=['cell_id', 'chrom', 'bin_start'], right_on=['cell_id', 'chr', 'start'], how='left')
    df_cell_clusters = pd.read_csv(f"{ROOT_DIR}/ov2295_clone_clusters.csv.gz")
    df["mutation"] = df["chrom"].astype(str) + ":" + df["coord"].astype(str) + ":" + df["ref"] + ":" + df["alt"]
    df = pd.merge(df, df_cell_clusters, on='cell_id', how='left')

    return df

df = get_ground_truth()

def get_alt_and_total(df_comb, clone_labels):

    clusters = np.array(list("ABCDEFGHI"))
    clones = np.sort(clone_labels['clone'].unique())

    n_clusters = len(clusters)
    n_clones = clones.max() + 1

    # Map cluster letters → indices
    cluster_to_idx = {c: i for i, c in enumerate(clusters)}

    # Clone numbers are already integers, so identity mapping is enough
    clone_to_idx = {c: c for c in clones}

    # Filter once
    valid = df_comb['state'].between(1, max_copy_number)

    df = df_comb.loc[valid, ['state','clone_id','clone','alt_counts','total_counts']].copy()

    # Convert to integer indices
    df['ci'] = df['clone_id'].map(cluster_to_idx)
    df['cj'] = df['clone'].map(clone_to_idx)

    df = df.dropna(subset=['ci','cj'])
    df['ci'] = df['ci'].astype(np.int32)
    df['cj'] = df['cj'].astype(np.int32)
    df['state'] = df['state'].astype(np.int32)

    # Allocate arrays
    alts = np.zeros((max_copy_number + 1, n_clusters, n_clones), dtype=np.int32)
    totals = np.zeros_like(alts)

    # Accumulate using vectorized indexed addition
    np.add.at(alts,  (df['state'], df['ci'], df['cj']), df['alt_counts'])
    np.add.at(totals,(df['state'], df['ci'], df['cj']), df['total_counts'])

    return alts, totals


# def get_alt_and_total(df_comb, clone_labels):

#     clusters = list("ABCDEFGHI")
#     clones = sorted(clone_labels['clone'].unique().tolist())
    

#     n_clusters = len(clusters)
#     n_clones = max(clones) + 1

#     cluster_to_idx = {c: i for i, c in enumerate(clusters)}
#     clone_to_idx = {c: c for i, c in enumerate(clones)}

#     alts = np.zeros((max_copy_number + 1, n_clusters, n_clones), dtype=np.int32)
#     totals = np.zeros_like(alts)

#     for _, row in df_comb.iterrows():
#         state = int(row['state'])
#         if 1 <= state <= max_copy_number:
#             try:
#                 ci = cluster_to_idx.get(row['clone_id'])
#                 cj = clone_to_idx.get(row['clone'])
#                 if ci is not None and cj is not None:
#                     alts[state, ci, cj] += row['alt_counts']
#                     totals[state, ci, cj] += row['total_counts']
#             except Exception as e:
#                 print(e)
#                 print(row)
    
#     return alts, totals



def hypothesis_test(alts, totals, dispersion_param=200, error_rate=0.01):

    pvals = np.ones_like(alts, dtype=np.float64)  # fill with 1s by default
    
    for copy_number in range(1, max_copy_number + 1):
        alt = alts[copy_number]
        total = totals[copy_number]

        # Effective success probability after accounting for error
        true_p = 1.0 / copy_number
        p_eff = true_p * (1 - error_rate)

        # Beta-binomial parameters
        ado_alpha = p_eff * dispersion_param
        ado_beta = dispersion_param * (1 - p_eff)

        # Compute p-values (vectorized)
        pvals_slice = betabinom.cdf(alt, total, ado_alpha, ado_beta)
        pvals[copy_number] = pvals_slice

    return pvals


def combine_weighted_pvalues(pvals, totals, eps=1e-15):
    """
    Combine p-values across copy numbers for each cluster-clone pair using weighted Stouffer.
    
    Parameters
    ----------
    pvals : np.ndarray
        3D array [copy_number, n_clusters, n_clones] of p-values
    totals : np.ndarray
        3D array [copy_number, n_clusters, n_clones] of weights (e.g., total_counts)
    eps : float
        Small value to clip p-values away from 0 or 1 to avoid infinities
    
    Returns
    -------
    combined : np.ndarray
        2D array [n_clusters, n_clones] of combined p-values
    """
    
    n_states, n_clusters, n_clones = pvals.shape
    combined = np.ones((n_clusters, n_clones))
    
    start_copy_number = 2  # skip 0 if you’re indexing by copy_number
    
    for i in range(n_clusters):
        for j in range(n_clones):
            pvals_ij = pvals[start_copy_number:, i, j]
            weights_ij = totals[start_copy_number:, i, j]

            # Mask valid entries with finite p-values and positive weights
            mask = np.isfinite(pvals_ij) & (weights_ij > 0)
            if np.any(mask):
                # Clip p-values to avoid -inf/+inf
                pvals_safe = np.clip(pvals_ij[mask], eps, 1 - eps)
                
                # Weighted Stouffer combination
                _, p_comb = combine_pvalues(
                    pvals_safe,
                    method='stouffer',
                    weights=weights_ij[mask]
                )
                combined[i, j] = p_comb
            else:
                combined[i, j] = np.nan
    return combined



def hypothesis_test_io(k, t, f):
    mutation_groups = pd.read_csv(f"data/laks/scope/reclustered_labels/filtered/k_{k}_t_{t}_f_{f}.csv", index_col=0)
    df_copy = df.copy()
    df_scope = pd.merge(df_copy, mutation_groups, on='mutation', how='left')
    alt_scope, total_scope = get_alt_and_total(df_scope, mutation_groups)
    pvals = hypothesis_test(alt_scope, total_scope)
    pvals_scope_combined = combine_weighted_pvalues(pvals, total_scope)
    df_out = pd.DataFrame(pvals_scope_combined, index=list("ABCDEFGHI"))
    os.makedirs(f"data/laks/scope/hypothesis_tests/filtered_reclustered/k_{k}_t_{t}_f_{f}", exist_ok=True)
    df_out.to_csv(f"data/laks/scope/hypothesis_tests/filtered_reclustered/k_{k}_t_{t}_f_{f}/combined_pvals.csv")



if __name__ == "__main__":
    for k in range(7, 26):
        for threshold in (0.60, 0.65, 0.70, 0.72, 0.75):
            for filter_threshold in [0.1, 0.2, 0.25, 0.3, 0.4]:
                try:
                    print(k, threshold, filter_threshold)
                    hypothesis_test_io(k, threshold, filter_threshold)
                except Exception as e:
                    print(e)
            
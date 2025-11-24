import numpy as np
import pandas as pd
from scipy.stats import betabinom, combine_pvalues


MIN_COPY_NUMBER = 0
MAX_COPY_NUMBER = 11


def get_alt_and_total_counts_for_cn_cluster_mutation_group(df_cell_mutation_cn_state, mutation_group_labels):

    df_comb = pd.merge(df_cell_mutation_cn_state, mutation_group_labels, on='mutation', how='left')

    clusters = np.array(list("ABCDEFGHI"))
    mutation_groups = np.sort(mutation_group_labels['mutation_group'].unique())

    n_clusters = len(clusters)
    n_mutation_groups = mutation_groups.max() + 1

    cluster_to_idx = {c: i for i, c in enumerate(clusters)}

    valid = df_comb['state'].between(1, MAX_COPY_NUMBER)

    df = df_comb.loc[valid, ['state','clone_id','mutation_group','alt_counts','total_counts']].copy()

    df['ci'] = df['clone_id'].map(cluster_to_idx)
    df['cj'] = df['mutation_group']

    df = df.dropna(subset=['ci','cj'])
    df['ci'] = df['ci'].astype(int)
    df['cj'] = df['cj'].astype(int)
    df['state'] = df['state'].astype(int)

    alt_counts = np.zeros((MAX_COPY_NUMBER + 1, n_clusters, n_mutation_groups), dtype=int)
    total_counts = np.zeros_like(alt_counts)

    np.add.at(alt_counts,  (df['state'], df['ci'], df['cj']), df['alt_counts'])
    np.add.at(total_counts,(df['state'], df['ci'], df['cj']), df['total_counts'])

    return alt_counts, total_counts


def hypothesis_test(alt_counts, total_counts, dispersion_param=300, error_rate=0.01):

    pvals = np.ones_like(alt_counts, dtype=float)
    
    for copy_number in range(1, MAX_COPY_NUMBER + 1):
        alt = alt_counts[copy_number]
        total = total_counts[copy_number]

        true_p = 1.0 / copy_number
        p_eff = true_p * (1 - error_rate)

        ado_alpha = p_eff * dispersion_param
        ado_beta = dispersion_param * (1 - p_eff)

        pvals_slice = betabinom.cdf(alt, total, ado_alpha, ado_beta)
        pvals[copy_number] = pvals_slice

    return pvals


def combine_weighted_pvalues(pvals, totals, eps=1e-15):
    
    n_states, n_clusters, n_clones = pvals.shape
    combined = np.ones((n_clusters, n_clones))
    
    start_copy_number = 2 
    
    for i in range(n_clusters):
        for j in range(n_clones):
            pvals_ij = pvals[start_copy_number:, i, j]
            weights_ij = totals[start_copy_number:, i, j]

            mask = np.isfinite(pvals_ij) & (weights_ij > 0)
            if np.any(mask):
                pvals_safe = np.clip(pvals_ij[mask], eps, 1 - eps)
                
                _, p_comb = combine_pvalues(
                    pvals_safe,
                    method='stouffer',
                    weights=weights_ij[mask]
                )
                combined[i, j] = p_comb
            else:
                combined[i, j] = np.nan
    return combined

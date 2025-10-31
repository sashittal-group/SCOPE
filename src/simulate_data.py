import pandas as pd
import sys
import argparse
import math
from scipy.stats import betabinom
import networkx as nx
import matplotlib.pyplot as plt 
import numpy as np

from simulation_helpers import draw_tree, draw_clone_tree, writeDOT, tree_to_newick, merge_cell_leaves
from phylogeny_utils import plot_spectral_clustering 


def main(args):

    np.random.seed(args.s)

    n_mutation_groups = args.m
    nclusters = args.p
    max_losses = args.k
    max_cn = args.maxcn
    mutation_rate = args.l
    mutation_grp_mean_size = args.size
    nnodes = n_mutation_groups + nclusters
    nbins = 1000
    
    event_order = list(np.random.permutation(['c'] * n_mutation_groups + ['d'] * (nclusters - 1)))

    cn_event_count = 1
    snv_event_count = 0
    for i in range(len(event_order)):
        if event_order[i] == 'c':
            event_order[i] = f"c{snv_event_count}"
            snv_event_count += 1
        else:
            event_order[i] = f"d{cn_event_count}"
            cn_event_count += 1

    # build tree
    T = nx.DiGraph() # mutation tree
    Tc = nx.DiGraph() # copy number tree

    # add root nodes
    T.add_node('root')
    Tc.add_node(0)

    B = np.zeros((n_mutation_groups + nclusters, n_mutation_groups + 1), dtype=int)
    num_children = np.zeros((nnodes, 1), dtype=int)

    for node_index, event in enumerate(event_order):
        nprev_mutations = sum([1 for x in event_order[:node_index] if x.startswith('c')])
        node_index += 1

        prev_num_children = num_children[:node_index].flatten()
        weights_a = (np.arange(start=1, stop=node_index+1, step=1, dtype=float) ** 3)
        weights_b = 1 / (1 + prev_num_children * prev_num_children)
        weights = weights_a * weights_b
        weights /= weights.sum()
        
        T_nodes = list(T.nodes)
        if event.startswith('d') and nprev_mutations > 0:
            choice_indices = []
            for node_index1, event1 in enumerate(T_nodes):
                if not (event1.startswith('d') or event1 == 'root'):
                    choice_indices.append(node_index1)
            
            weights_norm = weights[choice_indices]
            weights_norm = weights_norm / weights_norm.sum()
            parent_node_index = np.random.choice(choice_indices, p=weights_norm)
        else:
            parent_node_index = np.random.choice(range(node_index), p=weights)

        parent_node = T_nodes[parent_node_index]
        T.add_edge(parent_node, event)
        num_children[parent_node_index] += 1
        B[node_index, :] = B[parent_node_index, :]
    
        if event.startswith('d'):
            cluster_id = int(event.lstrip('d'))
            B[node_index, -1] = cluster_id
            parent_cluster_id = B[parent_node_index, -1]
            Tc.add_edge(parent_cluster_id, cluster_id)

        elif event.startswith('c'):
            mutation = int(event.lstrip('c'))
            B[node_index, mutation] = 1

        if args.v:
            print(parent_node_index, parent_node, node_index, event)
    
    mutation_group_sizes = np.random.poisson(lam=mutation_grp_mean_size, size=n_mutation_groups).tolist()
    n_mutations = sum(mutation_group_sizes)

    mutation_list = np.arange(n_mutations)
    mutation_list = list(np.random.permutation(mutation_list))
    
    mutation_group = np.repeat(np.arange(len(mutation_group_sizes)), mutation_group_sizes)

    df_mutation_group = pd.DataFrame({
        'mutation': mutation_list,
        'mutation_group': mutation_group
    }).sort_values('mutation').reset_index(drop=True)

    Rb = np.zeros((nclusters, nbins), dtype=int)
    Rb[0, :] = np.random.randint(max_cn - max_losses - 1, size = nbins) + max_losses + 1
    for cluster_id in range(nclusters):
        for bin in range(nbins):
            Rb[cluster_id, bin] = np.random.randint(max_cn - 1) + 1
    
    mut_to_bin = { mut: np.random.randint(0, nbins) for mut in range(n_mutations) }
    mut_bins = np.array([mut_to_bin[mut] for mut in range(n_mutations)])
    R = Rb[:, mut_bins]

    # check that all copy number states are non-zero positive
    assert(len(np.where(R == 0)[0]) == 0)
    
    # check all SNV losses are supported by CNVs
    
    if args.v:
        print('-'*50)
        print('copy number states')
        print('-'*50)
        print(R)    

    # assign cells and generate character-state matrix
    leaf_indices = []
    for idx, node in enumerate(T.nodes):
        if len(T[node]) == 0:
            leaf_indices.append(idx)    
    nleaves = len(leaf_indices)
    
    ncells = args.n
    assert(ncells > nclusters)
    print('n_mutation_groups:', n_mutation_groups, 'n_cells:', ncells, 'n_leaves:', nleaves)
    assert(ncells > nleaves)
    cell_assignment = np.random.randint(n_mutation_groups + nclusters, size=ncells-nleaves)
    complete_cell_assignment = list(cell_assignment) + leaf_indices
    
    Bcell = B[complete_cell_assignment, :]

    # observed matrix
    A = B.copy()
    for mutation in range(n_mutation_groups):
        A[A[:,mutation] > 1, mutation] = 0
    Acell = A[complete_cell_assignment, :]

    mut_groups = df_mutation_group["mutation_group"].to_numpy()
    Acell_muts = Acell[:, mut_groups]
    Acell_muts = np.hstack((Acell_muts, Acell[:, [-1]]))
    
    # cell tree
    celltree = T.copy()
    for cell_id, assigned_node_index in enumerate(complete_cell_assignment):
        celltree.add_edge(list(T.nodes)[assigned_node_index], f's{cell_id}')
        
    # generate read counts
    mean_coverage = args.cov
    fp_rate = args.a
    fn_rate = args.b
    ado_precision = args.ado

    cluster_ids = Acell_muts[:, -1]
    nvariant = Acell_muts[:, :-1]

    ntotal = R[cluster_ids[:, None], np.arange(n_mutations)[None, :]]

    latent_vaf = nvariant / ntotal

    Rtotal = np.random.poisson(mean_coverage, size=(ncells, n_mutations))
    post_error_vaf = fp_rate + (1 - fp_rate - fn_rate) * latent_vaf
    ado_alpha = post_error_vaf * ado_precision
    ado_beta = ado_precision * (1 - post_error_vaf)

    Vcount = betabinom.rvs(Rtotal, ado_alpha, ado_beta)
    CNs = ntotal
    
    # generate the binarized mutation matrix
    vaf_threshold = args.vafthreshold
    variant_read_threshold = args.readthreshold
    VAF_mat = Vcount / Rtotal
    mutation_mat = ((VAF_mat >= vaf_threshold) & (Vcount >= variant_read_threshold)).astype(int)
    mutation_mat = np.hstack((mutation_mat, Acell_muts[:,-1][:,np.newaxis]))
    
    # introduce missing entries
    Acell_missing = Acell_muts.copy()
    Rtotal_missing = Rtotal.copy()
    Vcount_missing = Vcount.copy()
    Acell_noisy = mutation_mat.copy()

    missing_rate = args.d
    n_entries = ncells * n_mutations
    nmissing = math.floor(missing_rate * n_entries)
    selected_cell_indices = np.random.randint(ncells, size=nmissing)
    selected_character_indices = np.random.randint(n_mutations, size=nmissing)
    Acell_missing[selected_cell_indices, selected_character_indices] = -1
    Rtotal_missing[selected_cell_indices, selected_character_indices] = 0
    Vcount_missing[selected_cell_indices, selected_character_indices] = 0
    Acell_noisy[selected_cell_indices, selected_character_indices] = -1
    
    # write ground truth files
    prefix = args.o
    with open(f'{prefix}_tree_edgelist.csv', 'w') as out:
        for edge in celltree.edges:
            out.write(f'{edge[0]},{edge[1]}\n')

    with open(f'{prefix}_tree.newick', 'w') as out:
        out.write(tree_to_newick(celltree) + ';')
    
    writeDOT(celltree, f'{prefix}_tree.dot')

    draw_tree(T, f"{prefix}_T.svg")
    draw_tree(merge_cell_leaves(celltree), f'{prefix}_merged_tree.svg')
    draw_clone_tree(Tc, f'{prefix}_cn_tree.svg')

    mutation_group_cn_cluster_columns = [f'c{idx}' for idx in range(n_mutation_groups)] + ['cluster_id']
    mutation_columns = [f'm{idx}' for idx in range(n_mutations)]
    mutation_cn_cluster_columns = mutation_columns  + ['cluster_id']
    cell_indices = [f's{idx}' for idx in range(ncells)]

    df_B = pd.DataFrame(B, index=list(T.nodes),
                        columns = mutation_group_cn_cluster_columns, dtype=int)            
    df_Bcell = pd.DataFrame(Bcell, index=cell_indices, columns = mutation_group_cn_cluster_columns, dtype=int)            
    df_Acell = pd.DataFrame(Acell_muts, index=cell_indices, columns = mutation_cn_cluster_columns, dtype=int)    
    df_Acell_noisy = pd.DataFrame(Acell_noisy, index=cell_indices, columns = mutation_cn_cluster_columns, dtype=int)
    df_Rtotal = pd.DataFrame(Rtotal, index=cell_indices, columns = mutation_columns, dtype=int)
    df_Vcount = pd.DataFrame(Vcount, index=cell_indices, columns = mutation_columns, dtype=int)    
    df_Rtotal_missing = pd.DataFrame(Rtotal_missing, index=cell_indices, columns = mutation_columns, dtype=int)
    df_Vcount_missing = pd.DataFrame(Vcount_missing, index=cell_indices, columns = mutation_columns, dtype=int)    
    df_CNs = pd.DataFrame(CNs, index=cell_indices, columns = mutation_columns, dtype=int)

    df_mut_to_bin = pd.DataFrame.from_dict(mut_to_bin, orient='index', columns=['bin'])

    df_mut_to_bin.to_parquet(f'{prefix}_mutation_to_bin_mapping.parquet', index=True)
    # df_B.to_parquet(f'{prefix}_multi_state_tree_node_character_matrix.parquet', index=True)
    df_Bcell.to_parquet(f'{prefix}_multi_state_character_matrix.parquet', index=True)
    df_Acell.to_parquet(f'{prefix}_character_matrix_without_noise.parquet', index=True)
    # df_Acell_noisy.to_parquet(f'{prefix}_character_matrix.parquet', index=True)

    # df_Rtotal.to_csv(f'{prefix}_read_count_without_missing.csv')
    # df_Vcount.to_csv(f'{prefix}_variant_count_without_missing.csv')

    df_Rtotal_missing.to_parquet(f'{prefix}_read_count.parquet', index=True)
    df_Vcount_missing.to_parquet(f'{prefix}_variant_count.parquet', index=True)
    df_CNs.to_parquet(f'{prefix}_copy_numbers.parquet', index=True)
    df_mutation_group.to_parquet(f'{prefix}_mutation_group.parquet', index=False)


    # Generate cell fraction
    df_total = df_Rtotal_missing
    df_variant = df_Vcount_missing
    df = df_Acell
    df_cluster = df[['cluster_id']]
    df_total = df_total.join(df_cluster, how='left')
    total_table = df_total.groupby('cluster_id').sum(numeric_only=True)
    df_variant = df_variant.join(df_cluster, how='left')
    alt_table = df_variant.groupby('cluster_id').sum(numeric_only=True)
    vaf = alt_table.div(total_table).replace(np.nan, 0)
    df_copy_number = df_CNs
    df_copy_number = df_copy_number.join(df_cluster, how='left')
    df_cn = df_copy_number.groupby('cluster_id').mean()
    F = (vaf * df_cn).clip(upper=1)

    # print(F)
    
    plot_spectral_clustering(F, n_clusters=F.shape[0], filepath=f"{prefix}_spectral_clustering.svg")

    clones = sorted(df_mutation_group["mutation_group"].unique())
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

        mutations_in_clone = df_mutation_group[df_mutation_group["mutation_group"] == clone]["mutation"].to_list()

        df = F.iloc[:, mutations_in_clone]
        df = df.loc[cn_order]

        ax.boxplot(df.T.values, vert=False, labels=df.index)

        ax.set_ylabel("CN Clusters")
        ax.set_xlabel("CF")
        ax.set_title(f"{df.shape[1]} mutations from\nground truth clone {clone}")

    # Hide any unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.savefig(f"{prefix}_clone_cluster_cell_fractions.svg")
    plt.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', type=int, help='number of samples [5]', default = 5)
    parser.add_argument('-m', type=int, help='number of SNV mutations [5]', default = 5)
    parser.add_argument('-p', type=int, help='number of clusters [1]', default = 1)
    parser.add_argument('-k', type=int, help='number of SNV losses per character [0]', default = 0)
    parser.add_argument('-o', type=str, help='output prefix', default='sample')
    parser.add_argument('-s', type=int, help='seed [0]', default = 0)
    parser.add_argument('--size', type=int, help='mutation group size [100]', default = 100)
    parser.add_argument('-d', type=float, help='missing data rate [0.0]', default=0)
    parser.add_argument('--cov', type=float, help='coverage of read count [50]', default = 50)
    parser.add_argument('-a', type=float, help='false positive error rate', default = 0)
    parser.add_argument('-b', type=float, help='false negative error rate', default = 0)
    parser.add_argument('--ado', type=float, help='precision parameter for ado [15]', default = 15)
    parser.add_argument('--maxcn', type=float, help='maximum allowed copy number [8]', default = 8)
    parser.add_argument('--readthreshold', type=int, help='variant read count threshold for generating the mutation matrix [5]', default=5)
    parser.add_argument('--vafthreshold', type=float, help='VAF threshold for generating the mutation matrix [0.1]', default = 0.1)
    parser.add_argument('-l', type=float, help='rate of mutation loss [0.8]', default = 0.8)
    parser.add_argument('-v', action='store_true', default=False)
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    
    main(args)

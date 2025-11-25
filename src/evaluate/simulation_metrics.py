import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

from scope.phylogeny_utils import generate_perfect_phylogeny, add_clusters_to_mutation_T, read_phertilizer_tree

def calculate_pairwise_ancestry_accuracy_for_scope(
        SIMULATION_STR, 
        type="GIVEN_MUTATION_CLUSTERING", 
        remove_unassigned_groups=False
    ):

    ## TODO: Provide files as arguments instead of filepaths
    GROUND_TRUTH_DIR = f"../data/simulation/ground_truth/{SIMULATION_STR}"

    ground_truth_mutation_group = pd.read_parquet(f"{GROUND_TRUTH_DIR}/sim_mutation_group.parquet")
    ground_truth_mutation_group.set_index('mutation', inplace=True)

    if type == "GIVEN_MUTATION_CLUSTERING":
        INPUT_FOLDER = "scope_input"
        OUTPUT_FOLDER = "scope_output"
        
        methods_mutation_group = ground_truth_mutation_group.copy()
    else:
        if type == "KMEANS_GIVEN_K":
            INPUT_FOLDER = "scope_input_kmeans_known_k"
            OUTPUT_FOLDER = "scope_output_kmeans_known_k"
        elif type == "KMEANS":
            INPUT_FOLDER = "scope_input_kmeans"
            OUTPUT_FOLDER = "scope_output_kmeans"
        else:
            raise ValueError("Unknown type", type)

        SCOPE_INPUT_DIR = f"../data/simulation/{INPUT_FOLDER}/{SIMULATION_STR}"
        
        methods_mutation_group = pd.read_csv(f"{SCOPE_INPUT_DIR}/kmeans_clones.csv", index_col=0)
        methods_mutation_group["mutation"] = methods_mutation_group["mutation"].str[1:].astype(int)
        methods_mutation_group.set_index("mutation", inplace=True)

    merged_mutation_group = pd.merge(ground_truth_mutation_group, methods_mutation_group, how='left', on='mutation', suffixes=['_ground_truth', '_method'])

    mutation_group_pair_counts = merged_mutation_group.value_counts(['mutation_group_ground_truth', 'mutation_group_method']).reset_index(name='mutation_counts')

    total_ground_truth_mutation_groups = int(ground_truth_mutation_group['mutation_group'].max()) + 1
    total_methods_mutation_groups = int(methods_mutation_group['mutation_group'].max()) + 1
    total_mutations = int(ground_truth_mutation_group.index.max()) + 1

    Tree = pd.read_csv(f"{GROUND_TRUTH_DIR}/sim_tree_edgelist.csv", header=None)
    Tree_mut = Tree[~Tree[1].str.startswith('s')]
    T_truth = nx.DiGraph()

    for index, row in Tree_mut.iterrows():
        parent, child = row[0], row[1]
        T_truth.add_edge(parent, child)

    B_mat = pd.read_csv(f"../data/simulation/{OUTPUT_FOLDER}/{SIMULATION_STR}/solution_0/B.csv", index_col=0)
    solT, _ = generate_perfect_phylogeny(B_mat)

    us = []; vs = []; rels = []

    for i in range(total_ground_truth_mutation_groups):
        for j in range(total_ground_truth_mutation_groups):
            us.append(i)
            vs.append(j)
            if i == j: rels.append(1)
            elif nx.has_path(T_truth, f'c{i}', f'c{j}'): rels.append(3)
            elif nx.has_path(T_truth, f'c{j}', f'c{i}'): rels.append(4)
            else: rels.append(2)

    ground_truth_ancestry = pd.DataFrame({ 'u1': us, 'v1': vs, 'rel1': rels })

    us = []; vs = []; rels = []

    for i in range(total_methods_mutation_groups):
        for j in range(total_methods_mutation_groups):
            us.append(i)
            vs.append(j)
            if i == j: rels.append(1)
            elif not solT.has_node(str(i)) or not solT.has_node(str(j)): rels.append(0) 
            elif nx.has_path(solT, str(i), str(j)): rels.append(3)
            elif nx.has_path(solT, str(j), str(i)): rels.append(4)
            else: rels.append(2)

    methods_ancestry = pd.DataFrame({ 'u2': us, 'v2': vs, 'rel2': rels })

    merged_ancestry = pd.merge(ground_truth_ancestry, methods_ancestry, how='cross')

    calc = pd.merge(merged_ancestry, mutation_group_pair_counts, left_on=['u1', 'u2'], right_on=['mutation_group_ground_truth', 'mutation_group_method'], how='left')
    calc.rename(columns={'mutation_counts': 'u1u2'}, inplace=True)
    calc.drop(columns=['mutation_group_ground_truth', 'mutation_group_method'], inplace=True)
    calc = pd.merge(calc, mutation_group_pair_counts, left_on=['v1', 'v2'], right_on=['mutation_group_ground_truth', 'mutation_group_method'], how='left')
    calc.rename(columns={'mutation_counts': 'v1v2'}, inplace=True)
    calc.drop(columns=['mutation_group_ground_truth', 'mutation_group_method'], inplace=True)
    calc['u1u2'] = calc['u1u2'].fillna(0)
    calc['v1v2'] = calc['v1v2'].fillna(0)
    calc['prod'] = calc['u1u2'] * calc['v1v2']

    if remove_unassigned_groups: calc = calc[calc['rel2'] != 0]

    calc['match'] = (calc['rel1'] == calc['rel2']).astype(int)

    total = (total_mutations * ( total_mutations - 1) ) // 2
    matched_prod = (calc['match'] * calc['prod']).sum()
    correct = (matched_prod - total_mutations) // 2
    
    accuracy = correct / total

    return total, correct, accuracy


def calculate_pairwise_ancestry_accuracy_for_scope_post(
        SIMULATION_STR, 
        type="KMEANS", 
        remove_unassigned_groups=False
    ):

    GROUND_TRUTH_DIR = f"../data/simulation/ground_truth/{SIMULATION_STR}"

    ground_truth_mutation_group = pd.read_parquet(f"{GROUND_TRUTH_DIR}/sim_mutation_group.parquet")
    ground_truth_mutation_group.set_index('mutation', inplace=True)

    
    if type == "KMEANS_GIVEN_K":
        INPUT_FOLDER = "scope_input_kmeans_known_k"
        OUTPUT_FOLDER = "scope_output_kmeans_known_k"
        POST_KMEANS_DIR = "scope_post_kmeans_known_k"
    elif type == "KMEANS":
        INPUT_FOLDER = "scope_input_kmeans"
        OUTPUT_FOLDER = "scope_output_kmeans"
        POST_KMEANS_DIR = "scope_post_kmeans"
    else:
        raise ValueError("Unknown type", type)

    SCOPE_INPUT_DIR = f"../data/simulation/{INPUT_FOLDER}/{SIMULATION_STR}"
    
    methods_mutation_group = pd.read_csv(f"../data/simulation/{POST_KMEANS_DIR}/{SIMULATION_STR}/kmeans_cleaned_clones.csv", index_col=0)
    methods_mutation_group["mutation"] = methods_mutation_group["mutation"].str[1:].astype(int)
    methods_mutation_group.set_index("mutation", inplace=True)

    merged_mutation_group = pd.merge(ground_truth_mutation_group, methods_mutation_group, how='left', on='mutation', suffixes=['_ground_truth', '_method'])

    mutation_group_pair_counts = merged_mutation_group.value_counts(['mutation_group_ground_truth', 'mutation_group_method']).reset_index(name='mutation_counts')

    total_ground_truth_mutation_groups = int(ground_truth_mutation_group['mutation_group'].max()) + 1
    total_methods_mutation_groups = int(methods_mutation_group['mutation_group'].max()) + 1
    total_mutations = int(ground_truth_mutation_group.index.max()) + 1

    Tree = pd.read_csv(f"{GROUND_TRUTH_DIR}/sim_tree_edgelist.csv", header=None)
    Tree_mut = Tree[~Tree[1].str.startswith('s')]
    T_truth = nx.DiGraph()

    for index, row in Tree_mut.iterrows():
        parent, child = row[0], row[1]
        T_truth.add_edge(parent, child)

    B_mat = pd.read_csv(f"../data/simulation/{OUTPUT_FOLDER}/{SIMULATION_STR}/solution_0/B.csv", index_col=0)
    solT, _ = generate_perfect_phylogeny(B_mat)

    us = []; vs = []; rels = []

    for i in range(total_ground_truth_mutation_groups):
        for j in range(total_ground_truth_mutation_groups):
            us.append(i)
            vs.append(j)
            if i == j: rels.append(1)
            elif nx.has_path(T_truth, f'c{i}', f'c{j}'): rels.append(3)
            elif nx.has_path(T_truth, f'c{j}', f'c{i}'): rels.append(4)
            else: rels.append(2)

    ground_truth_ancestry = pd.DataFrame({ 'u1': us, 'v1': vs, 'rel1': rels })

    us = []; vs = []; rels = []

    for i in range(total_methods_mutation_groups):
        for j in range(total_methods_mutation_groups):
            us.append(i)
            vs.append(j)
            if i == j: rels.append(1)
            elif not solT.has_node(str(i)) or not solT.has_node(str(j)): rels.append(0) 
            elif nx.has_path(solT, str(i), str(j)): rels.append(3)
            elif nx.has_path(solT, str(j), str(i)): rels.append(4)
            else: rels.append(2)

    methods_ancestry = pd.DataFrame({ 'u2': us, 'v2': vs, 'rel2': rels })

    merged_ancestry = pd.merge(ground_truth_ancestry, methods_ancestry, how='cross')

    calc = pd.merge(merged_ancestry, mutation_group_pair_counts, left_on=['u1', 'u2'], right_on=['mutation_group_ground_truth', 'mutation_group_method'], how='left')
    calc.rename(columns={'mutation_counts': 'u1u2'}, inplace=True)
    calc.drop(columns=['mutation_group_ground_truth', 'mutation_group_method'], inplace=True)
    calc = pd.merge(calc, mutation_group_pair_counts, left_on=['v1', 'v2'], right_on=['mutation_group_ground_truth', 'mutation_group_method'], how='left')
    calc.rename(columns={'mutation_counts': 'v1v2'}, inplace=True)
    calc.drop(columns=['mutation_group_ground_truth', 'mutation_group_method'], inplace=True)
    calc['u1u2'] = calc['u1u2'].fillna(0)
    calc['v1v2'] = calc['v1v2'].fillna(0)
    calc['prod'] = calc['u1u2'] * calc['v1v2']

    if remove_unassigned_groups: calc = calc[calc['rel2'] != 0]

    calc['match'] = (calc['rel1'] == calc['rel2']).astype(int)

    total = (total_mutations * ( total_mutations - 1) ) // 2
    matched_prod = (calc['match'] * calc['prod']).sum()
    correct = (matched_prod - total_mutations) // 2
    
    accuracy = correct / total

    return total, correct, accuracy


## TODO: Refactor code to reuse common parts of code
def calculate_pairwise_ancestry_accuracy_for_phertilizer(SIMULATION_STR, collapse=False):

    GROUND_TRUTH_DIR = f"../data/simulation/ground_truth/{SIMULATION_STR}"
    PHERTILIZER_OUTPUT_DIR = f"../data/simulation/phertilizer_output/{SIMULATION_STR}"

    mutation_group = pd.read_parquet(f"{GROUND_TRUTH_DIR}/sim_mutation_group.parquet")
    mutation_group.set_index('mutation', inplace=True)

    df_snv_clusters = pd.read_csv(f"{PHERTILIZER_OUTPUT_DIR}/snv_clusters.csv")
    df_snv_clusters['mutation'] = df_snv_clusters['mutation'].str[3:].astype(int)
    df_snv_clusters.set_index('mutation', inplace=True)
    df_snv_clusters.sort_index(inplace=True)

    df_merged = pd.merge(mutation_group, df_snv_clusters, left_index=True, right_index=True, how='left')
    df_merged['cluster'] = df_merged['cluster'].fillna(-1).astype(int)

    pair_counts = df_merged.value_counts(['mutation_group', 'cluster'])

    total_mutation_groups = int(df_merged['mutation_group'].max()) + 1
    total_mutation_groups_pred = int(df_merged['cluster'].max()) + 1
    total_mutations = int(df_merged.index.max()) + 1
    total_mutations_pred = (df_merged['cluster'] != -1).sum()

    # print(total_mutations_pred)

    Tree = pd.read_csv(f"{GROUND_TRUTH_DIR}/sim_tree_edgelist.csv", header=None)
    Tree_mut = Tree[~Tree[1].str.startswith('s')]
    T_truth = nx.DiGraph()

    for index, row in Tree_mut.iterrows():
        parent, child = row[0], row[1]
        T_truth.add_edge(parent, child)

    T_phert = nx.DiGraph()

    edges = []
    leaves = []

    with open(f"{PHERTILIZER_OUTPUT_DIR}/tree.txt", "r") as f:
        lines = f.readlines()
        
        n_edges = int(lines[0].split()[0])
        
        for line in lines[1:n_edges+1]:
            u, v = map(int, line.split())
            edges.append((u, v))
            T_phert.add_edge(u, v)
        

    mutation_group_ancestral_relation = np.zeros((total_mutation_groups, total_mutation_groups), dtype=int)

    for i in range(total_mutation_groups):
        for j in range(total_mutation_groups):
            if i == j: continue
            if nx.has_path(T_truth, f'c{i}', f'c{j}'):
                mutation_group_ancestral_relation[i][j] = 1

    us = []
    vs = []
    rel = []

    for i in range(total_mutation_groups):
        for j in range(total_mutation_groups):
            us.append(i)
            vs.append(j)
            if i == j: rel.append(1)
            elif nx.has_path(T_truth, f'c{i}', f'c{j}'): rel.append(3)
            elif nx.has_path(T_truth, f'c{j}', f'c{i}'): rel.append(4)
            else: rel.append(2)

    true_ancestry_df = pd.DataFrame({
        'u1': us,
        'v1': vs,
        'rel1': rel 
    })

    us = []
    vs = []
    rel = []

    for i in range(total_mutation_groups_pred):
        for j in range(total_mutation_groups_pred):
            us.append(i)
            vs.append(j)
            if i == j: rel.append(1)
            elif not T_phert.has_node(i) or not T_phert.has_node(j): rel.append(0)
            elif nx.has_path(T_phert, i, j): rel.append(3)
            elif nx.has_path(T_phert, j, i): rel.append(4)
            else: rel.append(2)

    pred_ancestry_df = pd.DataFrame({
        'u2': us,
        'v2': vs,
        'rel2': rel 
    })

    merged_ancestry = pd.merge(true_ancestry_df, pred_ancestry_df, how='cross')

    merged_ancestry_ = pd.merge(merged_ancestry, pair_counts, left_on=['u1', 'u2'], right_on=['mutation_group', 'cluster'], how='left')
    merged_ancestry_.rename(columns={'count': 'u1u2'}, inplace=True)
    merged_ancestry_ = pd.merge(merged_ancestry_, pair_counts, left_on=['v1', 'v2'], right_on=['mutation_group', 'cluster'], how='left')
    merged_ancestry_.rename(columns={'count': 'v1v2'}, inplace=True)
    merged_ancestry_ = merged_ancestry_.fillna(0)

    merged_ancestry_["prod"] = merged_ancestry_["u1u2"] * merged_ancestry_["v1v2"]
    merged_ancestry_["match"] = (merged_ancestry_["rel1"] == merged_ancestry_["rel2"]).astype(int)

    merged_ancestry_ = merged_ancestry_[merged_ancestry_["rel2"] != 0]

    # TODO: Fix the difference addition for scope    
    total_pairs = (total_mutations * (total_mutations - 1)) // 2

    match_prod = (merged_ancestry_["prod"] * merged_ancestry_["match"]).sum()
    correct_pairs = (match_prod - total_mutations_pred) // 2

    phert_acc = correct_pairs / total_pairs

    return total_pairs, correct_pairs, phert_acc


def calculate_pairwise_ancestry_accuracy_for_sbmclone(SIMULATION_STR):
    GROUND_TRUTH_DIR = f"../data/simulation/ground_truth/{SIMULATION_STR}"

    ground_truth_mutation_group = pd.read_parquet(f"{GROUND_TRUTH_DIR}/sim_mutation_group.parquet")
    ground_truth_mutation_group.set_index('mutation', inplace=True)

    Tree = pd.read_csv(f"{GROUND_TRUTH_DIR}/sim_tree_edgelist.csv", header=None)
    Tree_mut = Tree[~Tree[1].str.startswith('s')]
    T_truth = nx.DiGraph()

    for index, row in Tree_mut.iterrows():
        parent, child = row[0], row[1]
        T_truth.add_edge(parent, child)
    
    SBMCLONE_OUT_DIR = f"../data/simulation/sbmclone_output/{SIMULATION_STR}"

    with open(f"{SBMCLONE_OUT_DIR}/cluster-assignments.txt", "r") as f:
        lines = f.readlines()

    row_assignments = [int(x) for x in lines[0].strip().split(",")]
    col_assignments = [int(x) for x in lines[1].strip().split(",")]

    methods_mutation_group = pd.DataFrame({
        "mutation": np.arange(len(col_assignments)),
        "mutation_group": col_assignments
    })
    methods_mutation_group.set_index("mutation", inplace=True)
    methods_mutation_group['mutation_group'] = methods_mutation_group['mutation_group'] - 1

    merged_mutation_group = pd.merge(ground_truth_mutation_group, methods_mutation_group, how='left', on='mutation', suffixes=['_ground_truth', '_method'])

    mutation_group_pair_counts = merged_mutation_group.value_counts(['mutation_group_ground_truth', 'mutation_group_method']).reset_index(name='mutation_counts')

    total_ground_truth_mutation_groups = int(ground_truth_mutation_group['mutation_group'].max()) + 1
    total_methods_mutation_groups = int(methods_mutation_group['mutation_group'].max()) + 1
    total_mutations = int(ground_truth_mutation_group.index.max()) + 1

    us = []; vs = []; rels = []

    for i in range(total_ground_truth_mutation_groups):
        for j in range(total_ground_truth_mutation_groups):
            us.append(i)
            vs.append(j)
            if i == j: rels.append(1)
            elif nx.has_path(T_truth, f'c{i}', f'c{j}'): rels.append(3)
            elif nx.has_path(T_truth, f'c{j}', f'c{i}'): rels.append(4)
            else: rels.append(2)

    ground_truth_ancestry = pd.DataFrame({ 'u1': us, 'v1': vs, 'rel1': rels })
    
    blockmatrix = pd.read_csv(f"{SBMCLONE_OUT_DIR}/blockmatrix.csv", header=None)
    B_mat = (blockmatrix > 0.01).astype(int)
    solT, _ = generate_perfect_phylogeny(B_mat)

    us = []; vs = []; rels = []

    for i in range(total_methods_mutation_groups):
        for j in range(total_methods_mutation_groups):
            us.append(i)
            vs.append(j)
            if i == j: rels.append(1)
            elif not solT.has_node(str(i)) or not solT.has_node(str(j)): rels.append(0) 
            elif nx.has_path(solT, str(i), str(j)): rels.append(3)
            elif nx.has_path(solT, str(j), str(i)): rels.append(4)
            else: rels.append(2)

    methods_ancestry = pd.DataFrame({ 'u2': us, 'v2': vs, 'rel2': rels })

    merged_ancestry = pd.merge(ground_truth_ancestry, methods_ancestry, how='cross')

    calc = pd.merge(merged_ancestry, mutation_group_pair_counts, left_on=['u1', 'u2'], right_on=['mutation_group_ground_truth', 'mutation_group_method'], how='left')
    calc.rename(columns={'mutation_counts': 'u1u2'}, inplace=True)
    calc.drop(columns=['mutation_group_ground_truth', 'mutation_group_method'], inplace=True)
    calc = pd.merge(calc, mutation_group_pair_counts, left_on=['v1', 'v2'], right_on=['mutation_group_ground_truth', 'mutation_group_method'], how='left')
    calc.rename(columns={'mutation_counts': 'v1v2'}, inplace=True)
    calc.drop(columns=['mutation_group_ground_truth', 'mutation_group_method'], inplace=True)
    calc['u1u2'] = calc['u1u2'].fillna(0)
    calc['v1v2'] = calc['v1v2'].fillna(0)
    calc['prod'] = calc['u1u2'] * calc['v1v2']

    # if remove_unassigned_groups: calc = calc[calc['rel2'] != 0]

    calc['match'] = (calc['rel1'] == calc['rel2']).astype(int)

    total = (total_mutations * ( total_mutations - 1) ) // 2
    matched_prod = (calc['match'] * calc['prod']).sum()
    ## TODO: Need to do - total_pred_mutations instead of total_mutations
    correct = (matched_prod - total_mutations) // 2

    accuracy = correct / total

    # return total, correct, accuracy, calc, solT
    return total, correct, accuracy


def make_results_df(fn):
    # parameters
    ncells = [1000, 5000, 10000]
    n_mutation_groups = [5, 10, 15]
    mutation_group_sizes = [100, 500, 1000]
    nclusters = [5, 10, 15]
    coverages = [0.02, 0.05, 0.1]
    seeds = np.arange(5)

    ncells_list = []
    n_mutation_groups_list = []
    mutation_group_sizes_list = []
    n_clusters_list = []
    coverages_list = []
    seeds_list = []
    accuracies_list = []
    failures_list = []

    for ncell in ncells:
        for n_mutation_group in n_mutation_groups:
            for mutation_group_size in mutation_group_sizes:
                for coverage in coverages:
                    for ncluster in nclusters:
                        for seed in seeds:
                            ncells_list.append(ncell)
                            n_mutation_groups_list.append(n_mutation_group)
                            mutation_group_sizes_list.append(mutation_group_size)
                            coverages_list.append(coverage)
                            n_clusters_list.append(ncluster)
                            seeds_list.append(seed)
                            try:
                                SIMULATION_STR = f'n{ncell}_m{n_mutation_group}_size{mutation_group_size}_cov{str(coverage)}_p{ncluster}_s{seed}'
                                _, _, accuracy = fn(SIMULATION_STR)
                                accuracies_list.append(accuracy)
                                failures_list.append(None)
                            except Exception as e:
                                accuracies_list.append(None)
                                failures_list.append(e)

    results_df = pd.DataFrame({
        'ncells': ncells_list,
        'n_mutation_groups': n_mutation_groups_list,
        'mutation_group_size': mutation_group_sizes_list,
        'coverage': coverages_list,
        'n_clusters': n_clusters_list,
        'seed': seeds_list,
        'accuracy': accuracies_list,
        'error': failures_list,
    })

    return results_df


def nearest_d_or_root_for_c_nodes(G, root):
    result = {}

    for cnode in [n for n in G.nodes if str(n).startswith('c')]:

        current = cnode
        ancestor = None

        while current != root:
            preds = list(G.predecessors(current)) if G.is_directed() else [
                n for n in G.neighbors(current) if nx.shortest_path_length(G, root, n) < nx.shortest_path_length(G, root, current)
            ]
            if not preds:
                break
            parent = preds[0]

            if str(parent).startswith('d'):
                ancestor = parent
                break

            current = parent

        if ancestor is None:
            ancestor = root

        result[cnode] = ancestor

    return result


def nearest_CN_or_root_for_nodes(G, root):
    result = {}

    for cnode in [n for n in G.nodes if isinstance(n, (int, np.integer))]:

        current = cnode
        ancestor = None

        while current != root:
            preds = list(G.predecessors(current)) if G.is_directed() else [
                n for n in G.neighbors(current) if nx.shortest_path_length(G, root, n) < nx.shortest_path_length(G, root, current)
            ]
            if not preds:
                break
            parent = preds[0]

            if str(parent).startswith('CN_'):
                ancestor = parent
                break

            current = parent

        if ancestor is None:
            ancestor = root

        result[cnode] = ancestor

    return result


def cluster_placement_accuracy(STR, type="GIVEN_CLUSTERING"):

    if type == "GIVEN_CLUSTERING":
        GROUND_TRUTH_DIR = f"../data/simulation/ground_truth/{STR}"
        SCOPE_OUT_DIR = f"../data/simulation/scope_output/{STR}"
        SCOPE_IN_DIR = f"../data/simulation/scope_input/{STR}"
    elif type == "KMEANS_GIVEN_K":
        GROUND_TRUTH_DIR = f"../data/simulation/ground_truth/{STR}"
        SCOPE_OUT_DIR = f"../data/simulation/scope_output_kmeans_known_k/{STR}"
        SCOPE_IN_DIR = f"../data/simulation/scope_input_kmeans_known_k/{STR}"
    elif type == "KMEANS":
        GROUND_TRUTH_DIR = f"../data/simulation/ground_truth/{STR}"
        SCOPE_OUT_DIR = f"../data/simulation/scope_output_kmeans/{STR}"
        SCOPE_IN_DIR = f"../data/simulation/scope_input_kmeans/{STR}"
    elif type == "KMEANS_GIVEN_K_POST":
        GROUND_TRUTH_DIR = f"../data/simulation/ground_truth/{STR}"
        SCOPE_OUT_DIR = f"../data/simulation/scope_output_kmeans_known_k/{STR}"
        SCOPE_IN_DIR = f"../data/simulation/scope_post_kmeans_known_k/{STR}"
    elif type == "KMEANS_POST":
        GROUND_TRUTH_DIR = f"../data/simulation/ground_truth/{STR}"
        SCOPE_OUT_DIR = f"../data/simulation/scope_output_kmeans/{STR}"
        SCOPE_IN_DIR = f"../data/simulation/scope_post_kmeans/{STR}"
    else:
        raise ValueError(f"type unknown {type}")

    df_mut_mut_group = pd.read_parquet(f"{GROUND_TRUTH_DIR}/sim_mutation_group.parquet")

    Tree = pd.read_csv(f"{GROUND_TRUTH_DIR}/sim_tree_edgelist.csv", header=None)
    Tree_mut = Tree[~Tree[1].str.startswith('s')]
    T_truth = nx.DiGraph()

    for index, row in Tree_mut.iterrows():
        parent, child = row[0], row[1]
        T_truth.add_edge(parent, child)

    cluster_map = nearest_d_or_root_for_c_nodes(T_truth, 'root')

    cluster_df = pd.DataFrame(list(cluster_map.items()), columns=['mutation_group', 'cn_cluster'])
    cluster_df['cn_cluster'] = cluster_df['cn_cluster'].replace('root', 'd0')
    cluster_df['mutation_group'] = cluster_df['mutation_group'].str[1:].astype(int)
    cluster_df['cn_cluster'] = cluster_df['cn_cluster'].str[1:].astype(int)

    df_mut_cn_cluster = pd.merge(df_mut_mut_group, cluster_df, on='mutation_group', how='left')

    B = pd.read_csv(f"{SCOPE_OUT_DIR}/solution_0/B.csv", index_col=0)
    X = pd.read_csv(f"{SCOPE_OUT_DIR}/solution_0/X.csv", index_col=0)
    G = pd.read_csv(f"{SCOPE_OUT_DIR}/solution_0/G.csv", index_col=0)
    G.columns = G.columns.astype(int)
    B.columns = B.columns.astype(int)
    X.index = X.index.astype(int)

    solT_mut, _ = generate_perfect_phylogeny(B)
    fixed_T = add_clusters_to_mutation_T(solT_mut, X, G, B)

    cluster_map_pred = nearest_CN_or_root_for_nodes(fixed_T, 'root')

    cluster_df_pred = pd.DataFrame(list(cluster_map_pred.items()), columns=['mutation_group', 'cn_cluster'])
    cluster_df_pred['cn_cluster'] = cluster_df_pred['cn_cluster'].str[3:].map(lambda x: ord(x) - ord('A'))

    if type == "GIVEN_CLUSTERING":
        df_mut_mut_group_pred = pd.read_parquet(f"{GROUND_TRUTH_DIR}/sim_mutation_group.parquet")
    else:
        FILE_NAME = "kmeans_clones"
        if type[-4:] == "POST": FILE_NAME = "kmeans_cleaned_clones" 
        df_mut_mut_group_pred = pd.read_csv(f"{SCOPE_IN_DIR}/{FILE_NAME}.csv", index_col=0)
        df_mut_mut_group_pred['mutation'] = df_mut_mut_group_pred['mutation'].str[1:].astype(int)
    
    df_mut_cn_cluster_pred = pd.merge(df_mut_mut_group_pred, cluster_df_pred, on='mutation_group', how='left')

    df_merged = pd.merge(df_mut_cn_cluster, df_mut_cn_cluster_pred, on='mutation', how='left', suffixes=['_ground_truth', '_method'])
    df_merged['cn_cluster_method'].fillna(-1)
    df_merged['match'] = (df_merged['cn_cluster_ground_truth'] == df_merged['cn_cluster_method']).astype(int)

    total = int(len(df_merged))
    correct = int(df_merged['match'].sum())
    accuracy = correct / total

    return total, correct, accuracy


def find_lca(T, nodes, depths):
    if not nodes:
        return None
    ancestors_list = [set(nx.ancestors(T, n)).union({n}) for n in nodes]
    common_ancestors = set.intersection(*ancestors_list)
    lca = max(common_ancestors, key=lambda n: depths[n])
    return lca

def compute_distances_to_lcas(T, nodes, lca_nodes):
    import numpy as np
    
    dist_matrix = pd.DataFrame(np.inf, index=nodes, columns=lca_nodes, dtype=float)

    for lca in lca_nodes:
        lengths = nx.single_source_shortest_path_length(T, lca)
        for node, d in lengths.items():
            if node in nodes:
                dist_matrix.at[node, lca] = d

    return dist_matrix


def cluster_placement_accuracy_phertilizer(STR):
    GROUND_TRUTH_DIR = f"../data/simulation/ground_truth/{STR}"
    PHERTILIZER_OUT_DIR = f"../data/simulation/phertilizer_output/{STR}"
    PHERTILIZER_IN_DIR = f"../data/simulation/phertilizer_input/{STR}"

    df_mut_mut_group = pd.read_parquet(f"{GROUND_TRUTH_DIR}/sim_mutation_group.parquet")

    Tree = pd.read_csv(f"{GROUND_TRUTH_DIR}/sim_tree_edgelist.csv", header=None)
    Tree_mut = Tree[~Tree[1].str.startswith('s')]
    T_truth = nx.DiGraph()

    for index, row in Tree_mut.iterrows():
        parent, child = row[0], row[1]
        T_truth.add_edge(parent, child)
    
    cluster_map = nearest_d_or_root_for_c_nodes(T_truth, 'root')

    cluster_df = pd.DataFrame(list(cluster_map.items()), columns=['mutation_group', 'cn_cluster'])
    cluster_df['cn_cluster'] = cluster_df['cn_cluster'].replace('root', 'd0')
    cluster_df['mutation_group'] = cluster_df['mutation_group'].str[1:].astype(int)
    cluster_df['cn_cluster'] = cluster_df['cn_cluster'].str[1:].astype(int)

    df_mut_cn_cluster = pd.merge(df_mut_mut_group, cluster_df, on='mutation_group', how='left')

    T_phert = read_phertilizer_tree(f"{PHERTILIZER_OUT_DIR}/tree.txt")

    phertilizer_cell_cluster = pd.read_csv(f"{PHERTILIZER_OUT_DIR}/cell_clusters.csv")

    phertilizer_snv_cluster = pd.read_csv(f"{PHERTILIZER_OUT_DIR}/snv_clusters.csv")
    phertilizer_snv_cluster['mutation'] = phertilizer_snv_cluster['mutation'].str[3:].astype(int)

    ground_truth_cell_cluster = pd.read_parquet(f"{GROUND_TRUTH_DIR}/sim_character_matrix_without_noise.parquet")
    ground_truth_cell_cluster = ground_truth_cell_cluster.loc[:, 'cluster_id']

    clone_cell_cluster = pd.merge(phertilizer_cell_cluster, ground_truth_cell_cluster, right_index=True, left_on='cell', how='left')
    clone_cell_cluster.rename(columns={'cluster': 'clone'}, inplace=True)

    clone_cluster_counts = (
        clone_cell_cluster.groupby(["clone", "cluster_id"])
        .size()
        .reset_index(name="cell_count")
    )

    clone_totals = clone_cluster_counts.groupby("clone")["cell_count"].transform("sum")
    clone_cluster_counts["proportion"] = clone_cluster_counts["cell_count"] / clone_totals
 
    clone_cluster_filt = clone_cluster_counts[clone_cluster_counts['cell_count'] >= 3]

    n_clusters = int(ground_truth_cell_cluster.max()) + 1
    n_clones_pred = max(int(phertilizer_snv_cluster['cluster'].max()), int(phertilizer_cell_cluster['cluster'].max())) + 1

    start_cluster = np.zeros(shape=(n_clones_pred, n_clusters), dtype=int)

    for _, row in clone_cluster_filt.iterrows():
        start_cluster[int(row["clone"]), int(row["cluster_id"])] = 1

    start_cluster_df = pd.DataFrame(
        start_cluster,
        index=np.arange(n_clones_pred),
        columns=np.arange(n_clusters)
    )

    clones_per_cluster = {
        cluster: list(start_cluster_df.index[start_cluster_df[cluster] > 0])
        for cluster in start_cluster_df.columns
    }

    # --- Step 1: Detect root ---
    root_nodes = [n for n, deg in T_phert.in_degree() if deg == 0]
    if len(root_nodes) != 1:
        raise ValueError(f"Expected a single root, found {len(root_nodes)} nodes with in-degree 0")
    root_node = root_nodes[0]

    # --- Step 2: Compute depths from root ---
    depths = nx.shortest_path_length(T_phert, source=root_node)

    # --- Step 3: Prepare clone -> cluster mapping ---
    clones_per_cluster = {
        cluster: list(start_cluster_df.index[start_cluster_df[cluster] > 0])
        for cluster in start_cluster_df.columns
    }

    # --- Step 4: Define LCA helper ---
    def find_lca(T, nodes, depths):
        if not nodes:
            return None
        ancestors_list = [set(nx.ancestors(T, n)).union({n}) for n in nodes]
        common_ancestors = set.intersection(*ancestors_list)
        lca = max(common_ancestors, key=lambda n: depths[n])
        return lca

    # --- Step 5: Compute LCA per cluster ---
    cluster_lca = {}
    for cluster, clones in clones_per_cluster.items():
        cluster_lca[cluster] = find_lca(T_phert, clones, depths)

    # --- Step 6: Return as DataFrame ---
    lca_df = pd.DataFrame.from_dict(cluster_lca, orient='index', columns=['LCA_clone'])

    dist_df = compute_distances_to_lcas(T_phert, T_phert.nodes(), lca_df['LCA_clone'].unique())

    nearest_lcas = []

    for node in dist_df.index:
        if node in dist_df.columns:
            dist_df.at[node, node] = np.inf

    for idx, row in dist_df.iterrows():
        min_dist = row.min()
        nearest = row.index[row == min_dist][0]
        nearest_lcas.append({'target_node': idx, 'nearest_lcas': nearest})

    nearest_lca_df = pd.DataFrame(nearest_lcas)

    clone_cluster_mat = np.zeros((n_clones_pred, n_clusters), dtype=int)

    for _, row in nearest_lca_df.iterrows():
        node = row['target_node']
        nearest_lca = row['nearest_lcas']

        for i, _ in lca_df[lca_df['LCA_clone'] == nearest_lca].iterrows():
            parent = i
            clone_cluster_mat[node][parent] = 1

    clone_cluster_mat = clone_cluster_mat | start_cluster

    subclonal_df = pd.DataFrame(
        clone_cluster_mat,
        index=np.arange(n_clones_pred),
        columns=np.arange(n_clusters)
    )

    phertilizer_snv_cluster_subclonal = pd.merge(phertilizer_snv_cluster, subclonal_df, left_on='cluster', right_index=True, how='left')
    phertilizer_snv_cluster_subclonal.set_index('mutation', inplace=True)
    phertilizer_snv_cluster_subclonal = phertilizer_snv_cluster_subclonal.sort_index()

    df_merged = pd.merge(df_mut_cn_cluster, phertilizer_snv_cluster_subclonal, left_on='mutation', right_index=True, how='left')
    df_merged = df_merged.fillna(0)

    cols = list(range(n_clusters))
    df_pred = df_merged[cols]
    subclonal_pred = df_pred.to_numpy() 

    gt = np.zeros((len(df_merged), n_clusters), dtype=int)
    for cluster in range(n_clusters):
        gt[df_merged['cn_cluster'] == cluster, cluster] = 1
    
    correct_pred = gt * subclonal_pred 

    corrects = correct_pred.sum(axis=1)
    corrects

    pred_true = subclonal_pred.sum(axis=1)
    pred_true

    accuracies = corrects / pred_true

    accuracies = np.nan_to_num(accuracies, nan=0.0, posinf=0.0, neginf=0.0)

    accuracy = np.mean(accuracies)

    return 0, 0, float(accuracy)


def calculate_runtime(STR, type):

    if type == "SCOPE":
        file_sub = "scope_output_kmeans"
    elif type == "KMEANS":
        file_sub = "scope_input_kmeans"
    elif type == "PHERTILIZER":
        file_sub = "phertilizer_output"
    elif type == "SBMCLONE":
        file_sub = "sbmclone_output"
    
    filename = f"../data/simulation/{file_sub}/{STR}/benchmark"

    df = pd.read_csv(filename, sep='\t')
    if df.empty:
        print("File is empty.")
        return

    if pd.isna(df.at[0, 'cpu_time']):
        df.at[0, 'cpu_time'] = df.at[0, 's']
    
    row = df.iloc[0]


    return row['cpu_time'], row['max_vms'], row['s']
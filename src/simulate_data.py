#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 2022

@author: Palash Sashittal
"""

import pandas as pd
import sys
import argparse
import itertools
import math
import numpy as np
from scipy.stats import betabinom
from collections import defaultdict

import networkx as nx

def connected_node_groups(G, nodes):
    subG = G.subgraph(nodes).to_undirected()
    groups = list(nx.connected_components(subG))

    ordered_groups = []
    for group in groups:
        subG = G.subgraph(group)
        order = list(nx.topological_sort(subG))
        ordered_groups.append(order)
    return ordered_groups
    

def merge_linear_paths(G):
    G = G.copy()
    nodes = []
    for n in G.nodes:
        if str(n).startswith("c"):
            event_children = 0
            for child in list(G.successors(n)):
                if not str(child).startswith('s'):
                    event_children += 1
            if event_children <= 1:
                nodes.append(n)

    connected_nodes = connected_node_groups(G, nodes)

    for group in connected_nodes:
        main_node = group[0]
        for node in group[1:]:
            G = nx.contracted_nodes(G, main_node, node, self_loops=False)
        
        new_name = f"{main_node} (+{len(group)})"
        nx.relabel_nodes(G, {main_node: new_name}, copy=False)
    
    return G


def merge_cell_leaves(G: nx.DiGraph):
    G = G.copy()
    parent_mut = defaultdict(list)
    for n in G.nodes:
        if str(n).startswith("s"):
            p = list(G.predecessors(n))[0]
            parent_mut[p].append(str(n))

    for mut, cells in parent_mut.items():
        for cell in cells:
            G.remove_node(cell)
        G.add_edge(mut, f"{cells[0]} (+{len(cells)})")

    return G


def writeDOT(T, dot_file):
    with open(dot_file, 'w') as output:

        output.write(f'digraph N {{\n')
        output.write(f"\toverlap=\"false\"\n")
        output.write(f"\trankdir=\"TB\"\n")

        idx_dict = {}
        idx = 0
        for node in T.nodes:
            idx_dict[node] = idx
            output.write(f'\t{idx} [label=\"{node}\", style=\"bold\"];\n')
            idx += 1
        
        for edge in T.edges:
            output.write(f"\t{idx_dict[edge[0]]} -> {idx_dict[edge[1]]} [style=\"bold\"];\n")
        
        output.write(f'}}')    

def tree_to_newick(T, root=None):
    if root is None:
        roots = list(filter(lambda p: p[1] == 0, T.in_degree()))
        assert 1 == len(roots)
        root = roots[0][0]
    subgs = []
    while len(T[root]) == 1:
        root = list(T[root])[0]
    for child in T[root]:
        while len(T[child]) == 1:
            child = list(T[child])[0]
        if len(T[child]) > 0:
            child_newick = tree_to_newick(T, root=child)
            if child_newick != '()':
                subgs.append(child_newick)
        else:
            if child.startswith('s'):
                subgs.append(child)
    # return "(" + ','.join(map(str, subgs)) + ")"
    if len(subgs) == 1:
        return str(subgs[0])
    else:
        return "(" + ','.join(map(str, subgs)) + ")"

def draw_clone_tree(T, filename):
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import networkx as nx

    # Assign a color to each unique first letter
    letters = sorted({str(node)[0] for node in T.nodes})
    cmap = plt.cm.get_cmap('tab10', len(letters))
    letter_to_color = {letter: mcolors.to_hex(cmap(i)) for i, letter in enumerate(letters)}

    # Map each node to its corresponding color
    node_colors = [letter_to_color[str(node)[0]] for node in T.nodes]

    # Compute layout and draw
    pos = nx.nx_agraph.graphviz_layout(T, prog='dot')
    plt.figure(figsize=(6, 6))
    nx.draw(
        T, pos, with_labels=True, arrows=True,
        node_size=1000, node_color=node_colors,
        font_size=10, font_weight='bold', arrowsize=10,
        alpha=0.5
    )

    # Optional: add a legend for first-letter colors
    for letter, color in letter_to_color.items():
        plt.scatter([], [], color=color, label=letter)
    plt.legend(title="First Letter", loc="upper left")

    plt.savefig(filename)
    plt.close()



def draw_tree(T, filename):
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import networkx as nx

    # Assign a color to each unique first letter
    letters = sorted({str(node)[0] for node in T.nodes})
    cmap = plt.cm.get_cmap('tab10', len(letters))
    letter_to_color = {letter: mcolors.to_hex(cmap(i)) for i, letter in enumerate(letters)}

    # Map each node to its corresponding color
    node_colors = [letter_to_color[str(node)[0]] for node in T.nodes]

    # Compute layout and draw
    pos = nx.nx_agraph.graphviz_layout(T, prog='dot')
    plt.figure(figsize=(12, 24))
    nx.draw(
        T, pos, with_labels=True, arrows=True,
        node_size=800, node_color=node_colors,
        font_size=10, font_weight='bold', arrowsize=10,
        alpha=0.5
    )

    # Optional: add a legend for first-letter colors
    for letter, color in letter_to_color.items():
        plt.scatter([], [], color=color, label=letter)
    plt.legend(title="First Letter", loc="upper left")

    plt.savefig(filename)
    plt.close()

  

def main(args):

    np.random.seed(args.s)

    T = nx.DiGraph() # mutation tree
    Tc = nx.DiGraph() # copy number tree

    # add root nodes
    T.add_node('root')
    Tc.add_node(0)

    # build tree
    ncharacters = args.m
    nmutations = ncharacters
    nclusters = args.p
    max_losses = args.k
    max_cn = args.maxcn
    mutation_rate = args.l
    nnodes = ncharacters + nclusters
    
    character_list = [f'c{character_index}' for character_index in range(ncharacters)]
    cluster_list = [f'd{cluster_index}' for cluster_index in range(1, nclusters)]
    # TODO: Move CN events to front to avoid no mutations from a CN
    first_n_event = int(nnodes * 1)
    permutated_characters = np.random.permutation(character_list).astype(str).tolist()
    event_order = np.random.permutation(permutated_characters[:first_n_event] + cluster_list).astype(str).tolist() + \
        np.random.permutation(permutated_characters[first_n_event:]).astype(str).tolist()

    cn_event = 1
    for i in range(len(event_order)):
        if event_order[i].startswith("d"):
            event_order[i] = f"d{cn_event}"
            cn_event += 1

    loss_counter = np.zeros((ncharacters, 1))
    loss_dictionary = {f'd{cluster_index}': [] for cluster_index in range(1, nclusters)}
    
    B = np.zeros((ncharacters + nclusters, ncharacters + 1), dtype=int)
    R = np.zeros((nclusters, ncharacters), dtype=int)    
    # TODO: Map mutations to bins and assign CN for bins
    R[0, :] = np.random.randint(max_cn - max_losses - 1, size = ncharacters) + max_losses + 1

    Tc_parent = {
        "d2": "d1",
        "d3": "d2",
        "d6": "d5",
        "d7": "d5",
        "d8": "d5",
    }
    
    num_children = np.zeros((nnodes, 1), dtype=int)
    for node_index, event in enumerate(event_order):
        nprev_mutations = sum([1 for x in event_order[:node_index] if x.startswith('c')])
        node_index += 1

        prev_num_children = num_children[:node_index].flatten()
        weights_a = (np.arange(start=1, stop=node_index+1, step=1, dtype=float) ** 10)
        weights_b = 1 / (0.01 + 10 * prev_num_children * prev_num_children)  # prefer fewer children
        weights = weights_a * weights_b
        weights /= weights.sum()
        
        T_nodes = list(T.nodes)
        if event.startswith('d') and nprev_mutations > 0:
            choice_indices = []
            for node_index1, event1 in enumerate(T_nodes):
                if not (event1.startswith('d') or event1 == 'root'):
                    choice_indices.append(node_index1)
            
            if event in Tc_parent:
                parent_cluster = Tc_parent[event]
                filtered_choice_indices = []
                for i in choice_indices:
                    node = T_nodes[i]
                    if nx.has_path(T, parent_cluster, node):
                        filtered_choice_indices.append(i)
                choice_indices = filtered_choice_indices

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
            R[cluster_id, :] = R[parent_cluster_id, :]

            for mutation in range(ncharacters):
                if B[parent_node_index, mutation] == 1 and loss_counter[mutation] < max_losses:
                    if np.random.rand() < mutation_rate:
                        B[node_index, mutation] = loss_counter[mutation] + 2
                        loss_counter[mutation] += 1
                        loss_dictionary[event].append(mutation)
                        R[cluster_id, mutation] -= 1

        elif event.startswith('c'):
            mutation = int(event.lstrip('c'))
            B[node_index, mutation] = 1

        if args.v:
            print(parent_node_index, parent_node, node_index, event)

    # randomize the copy number states for mutations that have never been lost
    # TODO: Place into clusters.
    for mutation in range(nmutations):
        if loss_counter[mutation] == 0:
            for cluster_id in range(nclusters):
                R[cluster_id, mutation] = np.random.randint(max_cn - 1) + 1

    # check that all copy number states are non-zero positive
    assert(len(np.where(R == 0)[0]) == 0)
    
    # check all SNV losses are supported by CNVs
    for cn_edge in Tc.edges:
        for mutation in loss_dictionary[f'd{cn_edge[1]}']:
            assert(R[cn_edge[0], mutation] > R[cn_edge[1], mutation])    

    if args.v:
        print('-'*50)
        print('loss counter')
        print('-'*50)
        print(loss_counter)
        print('-'*50)
        print('loss dictionary')
        print('-'*50)
        print(loss_dictionary)
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
    print(ncharacters, ncells, nleaves)
    if ncells > nleaves:
        cell_assignment = np.random.randint(ncharacters, size=ncells-nleaves)
        complete_cell_assignment = list(cell_assignment) + leaf_indices
    else:
        cell_assignment = np.random.randint(ncharacters, size=ncells-1)
        complete_cell_assignment = list(cell_assignment) + leaf_indices[-1:]
    
    Bcell = B[complete_cell_assignment, :]
        
    # observed matrix
    A = B.copy()
    for mutation in range(ncharacters):
        A[A[:,mutation] > 1, mutation] = 0
    Acell = A[complete_cell_assignment, :]
    
    # cell tree
    celltree = T.copy()
    for cell_id, assigned_node_index in enumerate(complete_cell_assignment):
        celltree.add_edge(list(T.nodes)[assigned_node_index], f's{cell_id}')

        
    # generate read counts
    mean_coverage = args.cov
    fp_rate = args.a
    fn_rate = args.b
    ado_precision = args.ado

    Rtotal = np.zeros((ncells, nmutations), dtype=int)
    Vcount = np.zeros((ncells, nmutations), dtype=int)
    CNs = np.zeros((ncells, nmutations), dtype=int)
    for cell in range(ncells):
        for mutation in range(nmutations):
            cluster_id = Acell[cell, -1]
            nvariant = Acell[cell, mutation]
            ntotal = R[cluster_id, mutation]

            latent_vaf = nvariant / ntotal

            nreads = np.random.poisson(mean_coverage)
            Rtotal[cell, mutation] = int(nreads)

            post_error_vaf = fp_rate + (1 - fp_rate - fn_rate) * latent_vaf
            ado_alpha = post_error_vaf * ado_precision
            ado_beta = ado_precision * (1 - post_error_vaf)
            nvariant_reads = betabinom.rvs(nreads, ado_alpha, ado_beta)

            Vcount[cell, mutation] = int(nvariant_reads)
            CNs[cell, mutation] = ntotal

    # generate the binarized mutation matrix
    vaf_threshold = args.vafthreshold
    variant_read_threshold = args.readthreshold
    VAF_mat = Vcount / Rtotal
    mutation_mat = ((VAF_mat >= vaf_threshold) & (Vcount >= variant_read_threshold)).astype(int)
    mutation_mat = np.hstack((mutation_mat, Acell[:,-1][:,np.newaxis]))
    
    # introduce missing entries
    Acell_missing = Acell.copy()
    Rtotal_missing = Rtotal.copy()
    Vcount_missing = Vcount.copy()
    Acell_noisy = mutation_mat.copy()

    missing_rate = args.d
    n_entries = ncells * ncharacters
    nmissing = math.floor(missing_rate * n_entries)
    selected_cell_indices = np.random.randint(ncells, size=nmissing)
    selected_character_indices = np.random.randint(ncharacters, size=nmissing)
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

    draw_tree(celltree, f'{prefix}_tree.png')

    draw_tree(merge_cell_leaves(merge_linear_paths(celltree)), f'{prefix}_merged_tree.png')

    draw_clone_tree(Tc, f'{prefix}_copy_number_tree.png')

    df_B = pd.DataFrame(B, index=list(T.nodes),
                        columns = [f'c{idx}' for idx in range(ncharacters)] + ['cluster_id'], dtype=int)            
    df_Bcell = pd.DataFrame(Bcell, index=[f's{idx}' for idx in range(ncells)],
                            columns = [f'c{idx}' for idx in range(ncharacters)] + ['cluster_id'], dtype=int)            
    df_Acell = pd.DataFrame(Acell, index=[f's{idx}' for idx in range(ncells)],
                            columns = [f'c{idx}' for idx in range(ncharacters)] + ['cluster_id'], dtype=int)    
    df_Acell_noisy = pd.DataFrame(Acell_noisy, index=[f's{idx}' for idx in range(ncells)],
                                  columns = [f'c{idx}' for idx in range(ncharacters)] + ['cluster_id'], dtype=int)
    
    df_Rtotal = pd.DataFrame(Rtotal, index=[f's{idx}' for idx in range(ncells)],
                        columns = [f'c{idx}' for idx in range(ncharacters)], dtype=int)
    df_Vcount = pd.DataFrame(Vcount, index=[f's{idx}' for idx in range(ncells)],
                        columns = [f'c{idx}' for idx in range(ncharacters)], dtype=int)    
    df_Rtotal_missing = pd.DataFrame(Rtotal_missing, index=[f's{idx}' for idx in range(ncells)],
                        columns = [f'c{idx}' for idx in range(ncharacters)], dtype=int)
    df_Vcount_missing = pd.DataFrame(Vcount_missing, index=[f's{idx}' for idx in range(ncells)],
                        columns = [f'c{idx}' for idx in range(ncharacters)], dtype=int)    
    
    df_CNs = pd.DataFrame(CNs, index=[f's{idx}' for idx in range(ncells)],
                        columns = [f'c{idx}' for idx in range(ncharacters)], dtype=int)    

    df_B.to_csv(f'{prefix}_multi_state_tree_node_character_matrix.csv')
    df_Bcell.to_csv(f'{prefix}_multi_state_character_matrix.csv')
    df_Acell.to_csv(f'{prefix}_character_matrix_without_noise.csv')
    df_Acell_noisy.to_csv(f'{prefix}_character_matrix.csv')
    
    df_Rtotal.to_csv(f'{prefix}_read_count_without_missing.csv')
    df_Vcount.to_csv(f'{prefix}_variant_count_without_missing.csv')
    df_Rtotal_missing.to_csv(f'{prefix}_read_count.csv')
    df_Vcount_missing.to_csv(f'{prefix}_variant_count.csv')
    df_CNs.to_csv(f'{prefix}_copy_numbers.csv')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', type=int, help='number of samples [5]', default = 5)
    parser.add_argument('-m', type=int, help='number of SNV mutations [5]', default = 5)
    parser.add_argument('-p', type=int, help='number of clusters [1]', default = 1)
    parser.add_argument('-k', type=int, help='number of SNV losses per character [0]', default = 0)
    parser.add_argument('-o', type=str, help='output prefix', default='sample')
    parser.add_argument('-s', type=int, help='seed [0]', default = 0)
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
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np


def generate_perfect_phylogeny(df_binary):

    solT_mut = nx.DiGraph()
    solT_mut.add_node('root')

    solT_cell = nx.DiGraph()
    solT_cell.add_node('root')

    df_binary = df_binary[df_binary.sum().sort_values(ascending=False).index]    

    for cell_id, row in df_binary.iterrows():
        if cell_id == 'root':
            continue

        curr_node = 'root'
        for column in df_binary.columns[row.values == 1]:
            if column in solT_mut[curr_node]:
                curr_node = column
            else:
                if column in solT_mut.nodes:
                    raise NameError(f'{column} is being repeated')
                solT_mut.add_edge(curr_node, column)
                solT_cell.add_edge(curr_node, column)
                curr_node = column

        solT_cell.add_edge(curr_node, cell_id)   

    return solT_mut, solT_cell



def build_phylogeny(B):
    G = nx.DiGraph()
    root = "root"
    G.add_node(root, label="root")

    taxa_sets = []
    for cell_id, row in B.iterrows():
        if not row.any():
            continue
        chars = frozenset([j for j, val in enumerate(row) if val == 1])
        taxa_sets.append((cell_id, chars))

    sorted_cells_with_taxa_sets = sorted(taxa_sets, key=lambda x: len(x[1]))
    node_map = {frozenset(): root}

    for cell_id, taxa_set in sorted_cells_with_taxa_sets:
        parent = max([p for p in node_map if p.issubset(taxa_set)], key=len)
        parent_cell = node_map[parent]

        diff = taxa_set - parent  # mutations gained at this node
        mut_diffs = [B.columns[j] for j in sorted(diff)]
        edge_name = "\n".join(str(mut_diffs)) if mut_diffs else "empty"
        node_name = cell_id

        G.add_node(node_name, label=node_name)
        G.add_edge(parent_cell, node_name, label=edge_name)

        node_map[taxa_set] = node_name

    return G

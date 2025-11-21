import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import SpectralBiclustering
import numpy as np
import seaborn as sns
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt



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
        for column in df_binary.columns[row.values > 0.5]:
            if solT_mut.has_edge(curr_node, column):
                curr_node = column
            else:
                if column in solT_mut.nodes:
                    raise NameError(f'{column} is being repeated')
                solT_mut.add_edge(curr_node, column)
                solT_cell.add_edge(curr_node, column)
                curr_node = column

        if not str(cell_id).isdigit():
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



def draw_clone_tree(T, edge_labels=None, figsize=(6, 4), filepath: str = None):

    # Create labels with concise formatting
    labels = {node: str(node) for node in T.nodes()}

    # Use graphviz 'dot' layout for hierarchical structure
    pos = nx.nx_agraph.graphviz_layout(T, prog="dot", args="-Gnodesep=.5 -Granksep=50")

    cna_nodes = {n for n in T.nodes if str(n).startswith("CN")}

    # --- Step 2: assign color group by closest ancestor letter ---
    color_by_ancestor = {}

    # Perform BFS from each letter node to assign descendants
    for cna_node in cna_nodes:
        for node in nx.descendants(T, cna_node):
            # If node not yet assigned or this letter is closer (shallower)
            if node not in color_by_ancestor:
                color_by_ancestor[node] = cna_node
            else:
                # choose the closer ancestor (shorter path)
                old_letter = color_by_ancestor[node]
                if nx.shortest_path_length(T, cna_node, node) < nx.shortest_path_length(T, old_letter, node):
                    color_by_ancestor[node] = cna_node

    # Letter nodes are their own color
    for cna_node in cna_nodes:
        color_by_ancestor[cna_node] = cna_node

    # --- Step 3: map letters to colors ---
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors

    unique_letters = sorted(cna_nodes)
    cmap = cm.get_cmap('tab10', len(unique_letters))
    color_map = {letter: mcolors.to_hex(cmap(i)) for i, letter in enumerate(unique_letters)}

    node_colors = [color_map.get(color_by_ancestor.get(n, 'root'), '#cccccc') for n in T.nodes]

    # Create figure with adjusted size
    plt.figure(figsize=figsize)

    # Draw the tree with node colors
    nx.draw(
        T,
        pos,
        labels=labels,
        font_size=7,
        font_weight="bold",
        node_size=1000,
        edge_color="#555555",
        width=1.5,
        arrows=True,
        node_color=node_colors
    )

    if edge_labels: nx.draw_networkx_edge_labels(T, pos, edge_labels=edge_labels, font_size=6)

    plt.tight_layout()
    if filepath: plt.savefig(filepath)
    else: plt.show()
    plt.close()


def add_clusters_to_mutation_T(T: nx.DiGraph, X: pd.DataFrame, G: pd.DataFrame, B: pd.DataFrame,  ):
    selected_muts = X[ X > 0.5 ].dropna().index.tolist()
    depths = dict(nx.single_source_shortest_path_length(T, 'root'))

    def get_cluster_name(cluster):
        return f"CN_{cluster}"

    G_sub = G.loc[:, selected_muts]
    for cluster, muts in G_sub.iterrows():
        T.add_node(get_cluster_name(cluster))
        gained_muts = muts.index[muts > 0.5].tolist()

        if gained_muts:    
            root_muts = []
            for m in gained_muts:
                has_parent = False
                for mo in gained_muts:
                    if m == mo: continue

                    if np.all(B[mo] + 0.5 >= B[m]): # TODO: Rethink
                        has_parent = True
                        break

                if not has_parent:
                    root_muts.append(m)
            
            root_muts = np.array(root_muts)

            for m in root_muts:
                old_parent = next(T.predecessors(m))
                T.remove_edge(old_parent, m)
                T.add_edge(get_cluster_name(cluster), m)
                if not T.has_edge(old_parent, cluster):
                    T.add_edge(old_parent, get_cluster_name(cluster))
        
        else:
            B_cluster = B.loc[cluster]
            muts = B_cluster[B_cluster > 0.5].index.tolist()

            if muts: deepest_mut = max(muts, key=lambda n: depths.get(n, -1))
            else: deepest_mut = 'root'

            T.add_edge(deepest_mut, get_cluster_name(cluster))
    
    return T


def fix_T(B: pd.DataFrame, G: pd.DataFrame, T):

    nodes = T.nodes()
   
    for cluster, muts in G.iterrows():
        gained_muts = muts.index[muts == 1].tolist()

        gained_muts = [m for m in gained_muts if m in nodes]

        root_muts = []

        for m in gained_muts:
            has_parent = False
            for mo in gained_muts:
                if m == mo: continue

                if np.all(B[mo] >= B[m]):
                    has_parent = True
                    break

            if not has_parent:
                root_muts.append(m)
        
        root_muts = np.array(root_muts)

        for m in root_muts:
            old_parent = next(T.predecessors(m))
            T.remove_edge(old_parent, m)
            T.add_edge(cluster, m)

    return T

def canonical_form(G):
    if not nx.is_tree(G):
        raise ValueError("Graph must be a tree")
    if "root" not in G:
        raise ValueError("Tree must contain a node named 'root'")
    return _encode_tree(G, "root", None)

def _encode_tree(G, node, parent):
    children = [n for n in G.neighbors(node) if n != parent]
    encoded_children = sorted(_encode_tree(G, c, node) for c in children)
    # include the current node label in the encoding
    return f"{node}({''.join(encoded_children)})"


def read_phertilizer_tree(filename):

    with open(filename) as f:
        lines = f.read().strip().splitlines()

    num_edges = int(lines[0].split()[0])
    edge_lines = lines[1:1 + num_edges]
    edges = [tuple(map(int, line.split())) for line in edge_lines]

    num_leaves = int(lines[1 + num_edges].split()[0])
    leaf_lines = lines[2 + num_edges:2 + num_edges + num_leaves]
    leaves = [int(line.strip()) for line in leaf_lines]

    G = nx.DiGraph()
    G.add_edges_from(edges)

    return G

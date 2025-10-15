import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import SpectralBiclustering
import numpy as np
import seaborn as sns
import pandas as pd



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



def plot_spectral_clustering(F):
    model = SpectralBiclustering(n_clusters=7, method='log', random_state=0)
    model.fit(F.values)

    frac_biclust = F.iloc[np.argsort(model.row_labels_)]
    frac_biclust = frac_biclust.iloc[:, np.argsort(model.column_labels_)]

    sns.heatmap(frac_biclust, cmap='coolwarm')

    # Get where the cluster boundaries are
    row_order = np.argsort(model.row_labels_)
    col_order = np.argsort(model.column_labels_)

    row_clusters, row_counts = np.unique(model.row_labels_[row_order], return_counts=True)
    col_clusters, col_counts = np.unique(model.column_labels_[col_order], return_counts=True)

    row_lines = np.cumsum(row_counts)[:-1]
    col_lines = np.cumsum(col_counts)[:-1]

    # Plot with boundaries
    plt.figure(figsize=(6, 5))
    ax = sns.heatmap(frac_biclust, cmap='viridis', cbar=True)

    # Draw horizontal lines
    for r in row_lines:
        ax.axhline(r, color='white', lw=2)

    # Draw vertical lines
    for c in col_lines:
        ax.axvline(c, color='white', lw=2)

    plt.title("Spectral Biclustering of CF")
    plt.xlabel("Mutations")
    plt.ylabel("CN Cluster IDs")
    plt.show()


import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches



def draw_clone_tree(T):

    # Create labels with concise formatting
    labels = {node: str(node) for node in T.nodes()}

    # Use graphviz 'dot' layout for hierarchical structure
    pos = nx.nx_agraph.graphviz_layout(T, prog="dot", args="-Gnodesep=.5 -Granksep=50")

    letter_nodes = {n for n in T.nodes if str(n).isalpha() and len(str(n)) == 1}

    # --- Step 2: assign color group by closest ancestor letter ---
    color_by_ancestor = {}

    # Perform BFS from each letter node to assign descendants
    for letter in letter_nodes:
        for node in nx.descendants(T, letter):
            # If node not yet assigned or this letter is closer (shallower)
            if node not in color_by_ancestor:
                color_by_ancestor[node] = letter
            else:
                # choose the closer ancestor (shorter path)
                old_letter = color_by_ancestor[node]
                if nx.shortest_path_length(T, letter, node) < nx.shortest_path_length(T, old_letter, node):
                    color_by_ancestor[node] = letter

    # Letter nodes are their own color
    for letter in letter_nodes:
        color_by_ancestor[letter] = letter

    # --- Step 3: map letters to colors ---
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors

    unique_letters = sorted(letter_nodes)
    cmap = cm.get_cmap('tab10', len(unique_letters))
    color_map = {letter: mcolors.to_hex(cmap(i)) for i, letter in enumerate(unique_letters)}

    node_colors = [color_map.get(color_by_ancestor.get(n, 'root'), '#cccccc') for n in T.nodes]

    # Create figure with adjusted size
    plt.figure(figsize=(6, 4))

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


    plt.tight_layout()
    plt.show()


def add_clusters_to_clonal_T(T: nx.DiGraph, X: pd.DataFrame, G: pd.DataFrame, B: pd.DataFrame,  ):
    selected_muts = X[X>0.5].dropna().index.tolist()
    depths = dict(nx.single_source_shortest_path_length(T, 'root'))

    G_sub = G.loc[:, selected_muts]
    for cluster, muts in G_sub.iterrows():
        T.add_node(cluster)
        gained_muts = muts.index[muts == 1].tolist()
        if gained_muts:
            
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
                if not T.has_edge(old_parent, cluster):
                    T.add_edge(old_parent, cluster)
        
        else:
            B_cluster = B.loc[cluster]
            muts = B_cluster[B_cluster == 1].index.tolist()

            deepest_mut = max(muts, key=lambda n: depths.get(n, -1))

            T.add_edge(deepest_mut, cluster)
            # children = list(T.successors(deepest_mut))
            # for child in children:
            #     print("reassigning child", cluster, deepest_mut, child)
            #     if T.has_edge(deepest_mut, child):
            #         T.remove_edge(deepest_mut, child)
                # T.add_edge(cluster, child)
    
    return T

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

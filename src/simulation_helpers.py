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




def draw_clone_tree(T: nx.DiGraph, filename: str):
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
        node_size=2000, node_color=node_colors,
        font_size=15, arrowsize=20,
    )

    plt.savefig(filename)
    plt.close()

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import SpectralBiclustering
import numpy as np
import seaborn as sns


def plot_spectral_clustering(F, filepath=None):
    n_clusters = F.shape[0]
    model = SpectralBiclustering(n_clusters=n_clusters, method='log', random_state=0)
    model.fit(F.values)

    frac_biclust = F.iloc[np.argsort(model.row_labels_)]
    frac_biclust = frac_biclust.iloc[:, np.argsort(model.column_labels_)]

    # Get where the cluster boundaries are
    row_order = np.argsort(model.row_labels_)
    col_order = np.argsort(model.column_labels_)

    row_clusters, row_counts = np.unique(model.row_labels_[row_order], return_counts=True)
    col_clusters, col_counts = np.unique(model.column_labels_[col_order], return_counts=True)

    row_lines = np.cumsum(row_counts)[:-1]
    col_lines = np.cumsum(col_counts)[:-1]

    plt.figure(figsize=(6, 5))
    ax = sns.heatmap(frac_biclust, cmap='viridis', cbar=True)

    for r in row_lines:
        ax.axhline(r, color='white', lw=2)

    for c in col_lines:
        ax.axvline(c, color='white', lw=2)

    plt.title("Spectral Biclustering of Cell Fractions")
    plt.xlabel("Mutations")
    plt.ylabel("CN Cluster IDs")
    
    if filepath is None: plt.show()
    else: plt.savefig(filepath)
    plt.close()

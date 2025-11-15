import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

def run_kmeans_on_laks(F_prime, k):

    runs_per_k = 10

    X = F_prime.T.to_numpy()

    best_score = np.inf
    best_labels = None

    print("Running K-Means with k=", k)
    for i in range(runs_per_k):
        kmeans = KMeans(n_clusters=k, random_state=i * 100)
        labels = kmeans.fit_predict(X)
        score = silhouette_score(X, labels)
        if score < best_score:
            best_score = score
            best_labels = labels
    
    kmeans_labels = pd.DataFrame({
        'mutation': F_prime.columns.to_list(),
        'clone': best_labels
    })
    kmeans_labels

    return kmeans_labels


if __name__ == "__main__":
    F = pd.read_csv("data/laks/scope/F.csv", index_col=0)
    for k in range(7, 26):
        kmeans_labels = run_kmeans_on_laks(F, k)
        kmeans_labels.to_csv(f"data/laks/scope/kmeans_labels/k_{k}.csv")


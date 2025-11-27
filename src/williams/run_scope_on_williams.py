import os
import argparse
import sys
import pandas as pd

from scope.solve_ilp import solve_cppme
from scope.phylogeny_utils import generate_perfect_phylogeny, draw_clone_tree, add_clusters_to_mutation_T, canonical_form



def run_for_sample(SAMPLE_ID, k):

    PATH = f"outputs/scope/williams/{SAMPLE_ID}/mutation_clusters/k_{k}"

    loh_count_threshold = 50

    kmeans_labels = pd.read_csv(f"{PATH}/kmeans_labels.csv", index_col=0)
    cluster_weights = kmeans_labels.groupby("mutation_group").size().tolist()

    cn_cell_df = pd.read_csv(f"data/williams/signatures_dataset/clones_trees/{SAMPLE_ID}_clones.tsv", sep='\t')

    loh_conflicts = pd.read_csv(f"outputs/scope/williams/{SAMPLE_ID}/loh_conflicts.csv", index_col=0)
    loh_conflicts = (loh_conflicts > loh_count_threshold).astype(int)

    counts = cn_cell_df['clone_id'].value_counts()

    print(counts)

    proportions = counts / counts.sum()

    print(proportions)

    n_clusters = len(cn_cell_df['clone_id'].unique())

    threshold = (1 / 3) * (1 / n_clusters)
    filtered = proportions[proportions > threshold]
    if len(filtered) <= 3:
        threshold = (1 / 5) * (1 / n_clusters)
        filtered = proportions[proportions > threshold]

    print(threshold)
    
    # filtered = counts[counts >= 50]
    filtered = filtered.sort_index().index.to_list()
    print(filtered)

    F_hi = pd.read_csv(f"{PATH}/F_plus.csv", index_col=0)
    F_lo = pd.read_csv(f"{PATH}/F_minus.csv", index_col=0)

    filtered_existing = [f for f in filtered if f in F_hi.index and f in F_lo.index]

    F_plus = F_hi.loc[filtered_existing, :]
    F_minus = F_lo.loc[filtered_existing, :]
    loh_conflicts = (
        loh_conflicts
        .reindex(index=filtered_existing, columns=filtered_existing, fill_value=0)
    )
    print(loh_conflicts)

    forbidden_pairs = []

    if loh_conflicts is not None:
        for i in range(loh_conflicts.shape[0]):
            for j in range(loh_conflicts.shape[1]):
                if loh_conflicts.iloc[i, j] > 0.5: forbidden_pairs.append((j, i))
        forbidden_pairs = sorted(forbidden_pairs)
    

    print(F_plus)

    cluster_not_at_root = False
    only_mutation_tree_variant = True

    solutions, best_val = solve_cppme(F_plus=F_plus, F_minus=F_minus, cluster_weights=cluster_weights, 
                only_mutation_tree_variant=only_mutation_tree_variant, cluster_not_at_root=cluster_not_at_root,
                forbidden_pairs=forbidden_pairs)


    out_path = f"outputs/scope/williams/{SAMPLE_ID}/solutions_mut/k_{k}"
    os.makedirs(out_path, exist_ok=True)

    for i, solution in enumerate(solutions):

        if solution is None:
            break
        
        print("SOLUTION", i)

        
        solution_path = f"{out_path}/solution_{i}"
        os.makedirs(solution_path, exist_ok=True)

        X, B, U, F, G = solution
        
        X.to_csv(f"{solution_path}/X.csv")
        B.to_csv(f"{solution_path}/B.csv")
        U.to_csv(f"{solution_path}/U.csv")
        F.to_csv(f"{solution_path}/F.csv")
        G.to_csv(f"{solution_path}/G.csv")

        S, _ = generate_perfect_phylogeny(B)
        T = add_clusters_to_mutation_T(S, X, G, B)

        draw_clone_tree(T, filepath=f"{solution_path}/T.pdf")

    with open(f"{out_path}/summary.txt", 'w') as f:
        print(f"#UNIQUE SOLUTIONS: {len(solutions)} WITH VALUE {best_val}", file=f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', type=str, help='sample', default='sample')
    parser.add_argument('-k', type=int, help='k', default=15)
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    sample = args.sample
    k = args.k
    
    run_for_sample(sample, k)

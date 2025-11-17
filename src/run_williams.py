import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import KMeans
import numpy as np
from sklearn.metrics import silhouette_score

from src.phylogeny_utils import plot_spectral_clustering

from src.solve_ilp_2 import solve_cncff
import gurobipy as gp
from src.phylogeny_utils import generate_perfect_phylogeny, draw_clone_tree, add_clusters_to_clonal_T, canonical_form
import os

import argparse
import sys



def get_F_boundaries(F, mutation_group_mapping, thres=0.75):

    mutation_groups = sorted(mutation_group_mapping["clone"].unique())
    
    medians = np.zeros((F.shape[0], len(mutation_groups)))
    upper_percentiles = np.zeros((F.shape[0], len(mutation_groups)))
    lower_percentiles = np.zeros((F.shape[0], len(mutation_groups)))

    for i in range(len(mutation_groups)):
        mut_grp = mutation_groups[i]
        mutations_in_group = mutation_group_mapping[mutation_group_mapping["clone"] == mut_grp]["mutation"]

        df = F[mutations_in_group]
        F_clone = df.to_numpy()

        med = np.median(F_clone, axis=1)
        upp = np.percentile(F_clone, axis=1, q=int(100*thres))
        low = np.percentile(F_clone, axis=1, q=int(100*(1-thres)))

        medians[:, i] = med
        upper_percentiles[:, i] = upp
        lower_percentiles[:, i] = low
    
    F_bar = pd.DataFrame(medians, index=F.index, columns=mutation_groups)
    F_hi  = pd.DataFrame(upper_percentiles, index=F.index, columns=mutation_groups)
    F_lo  = pd.DataFrame(lower_percentiles, index=F.index, columns=mutation_groups)

    return F_bar, F_hi, F_lo


def run_for_sample(SAMPLE_ID):

    PATH = f"data/williams/scratch/{SAMPLE_ID}"

    kmeans_labels = pd.read_csv(f"{PATH}/kmeans_labels.csv", index_col=0)
    cluster_weights = kmeans_labels.groupby("clone").size().tolist()

    cn_cell_df = pd.read_csv(f"data/williams/signatures_dataset/clones_trees/{SAMPLE_ID}_clones.tsv", sep='\t')

    counts = cn_cell_df['clone_id'].value_counts()
    proportions = counts / counts.sum()

    n_clusters = len(cn_cell_df['clone_id'].unique())

    threshold = (1 / 3) * (1 / n_clusters)

    filtered = proportions[proportions > threshold]
    filtered = filtered.sort_index().index.to_list()
    print(filtered)

    F_hi = pd.read_csv(f"{PATH}/F_hi.csv", index_col=0)
    F_lo = pd.read_csv(f"{PATH}/F_lo.csv", index_col=0)

    filtered_existing = [f for f in filtered if f in F_hi.index and f in F_lo.index]

    F_plus = F_hi.loc[filtered_existing, :]
    F_minus = F_lo.loc[filtered_existing, :]

    n_clones = F_hi.shape[1] + 1

    cluster_not_at_root = False
    only_mutation_tree_variant = True

    while n_clones >= 1:
        print(f"TRYING {n_clones} CLONES")
        solution, model = solve_cncff(F_plus, F_minus, n_clones=n_clones,
                                        cluster_weights=cluster_weights, time_limit=60 * 60 * 24, 
                                        cluster_not_at_root=cluster_not_at_root,
                                        only_mutation_tree_variant=only_mutation_tree_variant)
        if solution is None: n_clones -= 1
        else: break
    
    X, _, _, _, _ = solution
    X = X.to_numpy().T[0]

    print(f"THERE IS A SOLUTION WITH {n_clones} CLONES")

    found_Bs = []
    solutions = []
    best_val = 0
    for i in range(1000):    
        solution, model = solve_cncff(F_plus, F_minus, n_clones=n_clones,
                                                cluster_weights=cluster_weights, time_limit=60 * 60 * 24, 
                                                cluster_not_at_root=cluster_not_at_root,
                                                found_Bs=found_Bs, only_mutation_tree_variant=only_mutation_tree_variant, Xs=X)
        if solution is None:
            break
        
        print("SOLUTION", i)
        found_B = solution[1].astype(int).to_numpy()
        found_Bs.append(found_B)
        solutions.append(solution)
        best_val = max(best_val, model.ObjVal)

    unique_solutions = []
    solution_strs = {}

    out_path = f"data/williams/scratch/{SAMPLE_ID}/scope_mut_2"
    os.makedirs(out_path, exist_ok=True)

    for i, solution in enumerate(solutions):
        try:
            solution_path = f"{out_path}/solution_{i}"
            os.makedirs(solution_path, exist_ok=True)

            X, B, U, F, G = solution
            
            X.to_csv(f"{solution_path}/X.csv")
            B.to_csv(f"{solution_path}/B.csv")
            U.to_csv(f"{solution_path}/U.csv")
            F.to_csv(f"{solution_path}/F.csv")
            G.to_csv(f"{solution_path}/G.csv")

            solT_mut, _ = generate_perfect_phylogeny(B)
            fixed_T = add_clusters_to_clonal_T(solT_mut, X, G, B)
            T_code = canonical_form(fixed_T)

            draw_clone_tree(fixed_T, f"{solution_path}/T.pdf")

            if T_code not in solution_strs:
                print(i)
                solution_strs[T_code] = i
                unique_solutions.append(solution)
            else:
                print(i, 'same as', solution_strs[T_code])
            
        except Exception as e:
            print(e)

    with open(f"{out_path}/summary.txt", 'w') as f:
        print(f"#UNIQUE SOLUTIONS: {len(unique_solutions)} WITH VALUE {best_val}", file=f)




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', type=str, help='sample', default='sample')
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    sample = args.sample
    
    run_for_sample(sample)

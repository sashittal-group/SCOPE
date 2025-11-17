from src.solve_ilp_2 import solve_cncff
import pandas as pd
import os
from src.phylogeny_utils import *
import argparse
import os


def run_scope(F_plus, F_minus, mutation_cluster_labels, loh_conflicts=None, cluster_not_at_root=False):
    cluster_weights = mutation_cluster_labels.groupby("clone").size().tolist()

    restriction_pairs = []

    if loh_conflicts is not None:
        for i in range(loh_conflicts.shape[0]):
            for j in range(loh_conflicts.shape[1]):
                if loh_conflicts.iloc[i, j] > 0.5: restriction_pairs.append((j, i))
        restriction_pairs = sorted(restriction_pairs)

    total_clones = F_plus.shape[1] + 1
    n_clones = total_clones

    while n_clones >= 1:
        print(f"TRYING {n_clones} CLONES")
        solution, model = solve_cncff(F_plus, F_minus, n_clones=n_clones,
                                        cluster_weights=cluster_weights, time_limit=60 * 60 * 24, 
                                        cluster_parent_restriction_pairs=restriction_pairs, 
                                        cluster_not_at_root=cluster_not_at_root)
        if solution is None: n_clones -= 1
        else: break

    print(f"THERE IS A SOLUTION WITH {n_clones} CLONES")

    X, _, _, _, _ = solution
    X = X.to_numpy().T[0]

    print(f"THERE IS A SOLUTION WITH {n_clones} CLONES")

    found_Bs = []
    solutions = []
    best_val = 0
    for i in range(1000):    
        solution, model = solve_cncff(F_plus, F_minus, n_clones=n_clones,
                                                cluster_weights=cluster_weights, time_limit=24* 60 * 10, 
                                                cluster_parent_restriction_pairs=restriction_pairs, 
                                                cluster_not_at_root=cluster_not_at_root,
                                                found_Bs=found_Bs, Xs=X)
        if solution is None:
            break
        
        print("SOLUTION", i)
        found_B = solution[1].astype(int).to_numpy()
        found_Bs.append(found_B)
        solutions.append(solution)
        best_val = max(best_val, model.ObjVal)
        
    return solutions, best_val


def output_scope(F_plus_filepath, F_minus_filepath, mutation_cluster_filepath, cluster_not_at_root, loh_conflicts_filepath, output_directory):
    ## OUTPUT THE SOLUTIONS
    F_hi = pd.read_csv(F_plus_filepath, index_col=0)
    F_lo = pd.read_csv(F_minus_filepath, index_col=0)
    mutation_cluster_labels = pd.read_csv(mutation_cluster_filepath, index_col=0)
    loh_conflicts = pd.read_csv(loh_conflicts_filepath, index_col=0)
    loh_conflicts = (loh_conflicts > 50).astype(int)

    solutions, best_val = run_scope(F_hi, F_lo, mutation_cluster_labels, loh_conflicts, cluster_not_at_root=cluster_not_at_root)

    unique_solutions = []
    solution_strs = {}

    for i, solution in enumerate(solutions):
        try:
            X_path = f"{output_directory}/X"
            B_path = f"{output_directory}/B"
            U_path = f"{output_directory}/U"
            F_path = f"{output_directory}/F"
            G_path = f"{output_directory}/G"
            T_fig_path = f"{output_directory}/T_figure"

            os.makedirs(X_path, exist_ok=True)
            os.makedirs(B_path, exist_ok=True)
            os.makedirs(U_path, exist_ok=True)
            os.makedirs(F_path, exist_ok=True)
            os.makedirs(G_path, exist_ok=True)
            os.makedirs(T_fig_path, exist_ok=True)

            X, B, U, F, G = solution
            
            X.to_csv(f"{X_path}/solution_{i}.csv")
            B.to_csv(f"{B_path}/solution_{i}.csv")
            U.to_csv(f"{U_path}/solution_{i}.csv")
            F.to_csv(f"{F_path}/solution_{i}.csv")
            G.to_csv(f"{G_path}/solution_{i}.csv")

            solT_mut, _ = generate_perfect_phylogeny(B)
            fixed_T = add_clusters_to_clonal_T(solT_mut, X, G, B)
            T_code = canonical_form(fixed_T)

            draw_clone_tree(fixed_T, f"{T_fig_path}/solution_{i}.pdf")
            
            if T_code not in solution_strs:
                print(i)
                solution_strs[T_code] = i
                unique_solutions.append(solution)
            else:
                print(i, 'same as', solution_strs[T_code])
            
        except Exception as e:
            print(e)

    with open(f"{output_directory}/summary.txt", 'w') as f:
        print(f"#UNIQUE SOLUTIONS: {len(unique_solutions)} WITH VALUE {best_val}", file=f)
        
# if __name__ == "__main__":

#     loh_conflicts_filepath = "data/laks/scope/loh-counts.csv"

#     for k in range(7, 26):
#         for threshold in (0.60, 0.65, 0.70, 0.72, 0.75, 0.80):
#             F_plus_filepath = f"data/laks/scope/cell_fractions/k_{k}_t_{threshold}/F_plus.csv"
#             F_minus_filepath = f"data/laks/scope/cell_fractions/k_{k}_t_{threshold}/F_minus.csv"
#             mutation_labels_filepath = f"data/laks/scope/kmeans_labels/k_{k}.csv"
#             output_directory = f"data/laks/scope/outputs/k_{k}_t_{threshold}"

#             os.makedirs(output_directory, exist_ok=True)

#             output_scope(F_plus_filepath=F_plus_filepath, F_minus_filepath=F_minus_filepath, mutation_cluster_filepath=mutation_labels_filepath,
#                          cluster_not_at_root=False, loh_conflicts_filepath=loh_conflicts_filepath, output_directory=output_directory) 



def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--loh_conflicts_filepath",
        type=str,
        required=True
    )
    parser.add_argument(
        "--F_plus_filepath",
        type=str,
        required=True
    )
    parser.add_argument(
        "--F_minus_filepath",
        type=str,
        required=True
    )
    parser.add_argument(
        "--mutation_labels_filepath",
        type=str,
        required=True
    )
    parser.add_argument(
        "--output_directory",
        type=str,
        required=True
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    os.makedirs(args.output_directory, exist_ok=True)

    output_scope(
        F_plus_filepath=args.F_plus_filepath,
        F_minus_filepath=args.F_minus_filepath,
        mutation_cluster_filepath=args.mutation_labels_filepath,
        cluster_not_at_root=False,
        loh_conflicts_filepath=args.loh_conflicts_filepath,
        output_directory=args.output_directory,
    )

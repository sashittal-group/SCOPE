import argparse
import sys
import pandas as pd
import os
import string
import gurobipy as gp

from src.solve_ilp_2 import solve_cncff
from src.phylogeny_utils import *


def main(args):

    input_prefix = args.i
    output_root_prefix = args.o
    input_folder = input_prefix.split('/')[-1]
    output_prefix = f"{output_root_prefix}/{input_folder}"

    F_plus = pd.read_csv(f"{input_prefix}/F_plus.csv", index_col=0)
    F_minus = pd.read_csv(f"{input_prefix}/F_minus.csv", index_col=0)
    clone_sizes = pd.read_csv(f"{input_prefix}/clone_sizes.csv", index_col=0)

    clone_sizes_list = clone_sizes['0'].to_list()

    F_plus.index = list(string.ascii_uppercase[:len(F_plus)])
    F_minus.index = list(string.ascii_uppercase[:len(F_minus)])

    print(F_plus)
    print(F_minus)
    print(clone_sizes)
    print(clone_sizes_list)

    total_clones = len(clone_sizes_list) + 1 # +1 for root node

    n_clones = total_clones 
    while n_clones >= 1:
        solutions, best_objective, status = solve_cncff(F_plus, F_minus, n_clones=n_clones, n_solutions=1,
                                cluster_weights=clone_sizes_list, time_limit=2*60)
        if status == gp.GRB.TIME_LIMIT:
            solutions, best_objective, status = solve_cncff(F_plus, F_minus, n_clones=n_clones, n_solutions=1,
                                cluster_weights=clone_sizes_list, time_limit=10*60)
        # if status == gp.GRB.TIME_LIMIT:
        #     solutions, best_objective, status = solve_cncff(F_plus, F_minus, n_clones=total_clones, n_solutions=1,
        #                         cluster_weights=clone_sizes_list, time_limit=30*60)
        print("CLONES:", n_clones, ", SOLUTION COUNT:", len(solutions))
        if len(solutions) > 0:
            break
        n_clones -= 1
    
    unique_solutions = []
    solution_strs = {}

    for i, solution in enumerate(solutions):
        try:
            solution_path = f"{output_prefix}/solution_{i}"
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

            draw_clone_tree(fixed_T, f"{solution_path}/T.svg")

            if T_code not in solution_strs:
                print(i)
                solution_strs[T_code] = i
                unique_solutions.append(solution)
            else:
                print(i, 'same as', solution_strs[T_code])
            
        except Exception as e:
            print(e)
    
    print("UNIQUE SOLUTIONS:", len(unique_solutions))
    
    with open(f"{output_prefix}/summary.txt", 'w') as f:
        print(f"#UNIQUE SOLUTIONS: {len(unique_solutions)} WITH VALUE {best_objective}", file=f)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='input prefix', default='sample')
    parser.add_argument('-o', type=str, help='output prefix', default = 'scope_input')
    parser.add_argument('-v', action='store_true', default=False)
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    
    main(args)

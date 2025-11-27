import argparse
import sys
import pandas as pd
import os
import string
import gurobipy as gp

from src.scope.solve_ilp import solve_cppme
from src.scope.phylogeny_utils import *


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


    solutions, best_value = solve_cppme(F_plus, F_minus, cluster_weights=clone_sizes_list)
    
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

            S, _ = generate_perfect_phylogeny(B)
            T = add_clusters_to_mutation_T(S, X, G, B)
            T_code = canonical_form(T)

            draw_clone_tree(T, filepath=f"{solution_path}/T.svg")

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
        print(f"#UNIQUE SOLUTIONS: {len(unique_solutions)} WITH VALUE {best_value}", file=f)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='input prefix', default='sample')
    parser.add_argument('-o', type=str, help='output prefix', default = 'scope_input')
    parser.add_argument('-v', action='store_true', default=False)
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    
    main(args)
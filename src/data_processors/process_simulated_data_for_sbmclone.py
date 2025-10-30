import argparse
import sys
import pandas as pd
import os
import numpy as np

def main(args):

    input_prefix = args.i
    output_root_prefix = args.o
    mutation_threshold = args.threshold

    input_folder = input_prefix.split('/')[-2]
    output_prefix = output_root_prefix + input_folder
    os.makedirs(output_prefix, exist_ok=True)

    df_read_counts = pd.read_csv(f"{input_prefix}_read_count.csv", index_col=0)
    df_variant_counts = pd.read_csv(f"{input_prefix}_variant_count.csv", index_col=0)

    vaf = df_variant_counts.div(df_read_counts).replace(np.nan, 0)

    is_mutated = (vaf >= mutation_threshold).astype(int)

    rows, cols = np.where(is_mutated.values > 0.5)

    sbm_input = pd.DataFrame({
        0: rows,
        1: cols
    })

    sbm_input.to_csv(f"{output_prefix}/matrix.csv", index=False, header=False)




    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='input prefix', default='sample')
    parser.add_argument('-o', type=str, help='output prefix', default = 'scope_input')
    parser.add_argument('--threshold', type=float, help='VAF threshold for mutation', default = 0.1)
    parser.add_argument('-v', action='store_true', default=False)
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    
    main(args)

import argparse
import sys
import pandas as pd
import numpy as np
import os

from scope.process_data import get_cell_fraction_boundaries

def main(args):

    input_prefix = args.i
    clone_table_file = args.clone_table
    output_root_prefix = args.o

    df_total = pd.read_parquet(f"{input_prefix}_read_count.parquet")
    df_variant = pd.read_parquet(f"{input_prefix}_variant_count.parquet")
    df = pd.read_parquet(f"{input_prefix}_character_matrix_without_noise.parquet")
    df_copy_number = pd.read_parquet(f"{input_prefix}_copy_numbers.parquet")
    df_mutation_group = pd.read_parquet(f"{input_prefix}_mutation_group.parquet")

    df_cluster = df[['cluster_id']]

    df_total = df_total.join(df_cluster, how='left')
    df_variant = df_variant.join(df_cluster, how='left')
    df_copy_number = df_copy_number.join(df_cluster, how='left')

    total_table = df_total.groupby('cluster_id').sum(numeric_only=True)
    alt_table = df_variant.groupby('cluster_id').sum(numeric_only=True)
    df_cn = df_copy_number.groupby('cluster_id').mean().astype(int)

    vaf = alt_table.div(total_table).replace(np.nan, 0)

    F = (vaf * df_cn).clip(upper=1)

    if clone_table_file is None: return
    
    df_clone = pd.read_parquet(f"{input_prefix}{clone_table_file}")
    F_bar, F_hi, F_lo = get_cell_fraction_boundaries(F, df_clone, 0.75)

    input_folder = input_prefix.split('/')[-2]

    output_prefix = output_root_prefix + input_folder

    os.makedirs(output_prefix, exist_ok=True)

    mutation_group_sizes = df_mutation_group.groupby('mutation_group').size()

    F_hi.to_csv(f"{output_prefix}/F_plus.csv")
    F_lo.to_csv(f"{output_prefix}/F_minus.csv")
    F_bar.to_csv(f"{output_prefix}/F_bar.csv")
    mutation_group_sizes.to_csv(f"{output_prefix}/clone_sizes.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='input prefix', default='sample')
    parser.add_argument('-o', type=str, help='output prefix', default = 'scope_input')
    parser.add_argument('--clone-table', type=str, help='clone table file', default = None)
    parser.add_argument('-v', action='store_true', default=False)
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    
    main(args)
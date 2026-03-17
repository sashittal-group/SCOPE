import argparse
import sys
import pandas as pd
import numpy as np
import os

def main(args):

    input_prefix = args.i
    output_root_prefix = args.o

    input_folder = input_prefix.split('/')[-2]
    output_prefix = output_root_prefix + input_folder
    os.makedirs(output_prefix, exist_ok=True)

    df_read_counts = pd.read_parquet(f"{input_prefix}_read_count.parquet")
    df_variant_counts = pd.read_parquet(f"{input_prefix}_variant_count.parquet")
    df_mut_bin = pd.read_parquet(f"{input_prefix}_mutation_to_bin_mapping.parquet")

    df_read_counts = df_read_counts.loc[:, (df_read_counts != 0).any()]

    df_total_long = df_read_counts.reset_index().melt(id_vars='index', var_name='mutation', value_name='total')
    df_total_long = df_total_long.rename(columns={'index': 'cell_id'})
    
    df_var_long = df_variant_counts.reset_index().melt(id_vars='index', var_name='mutation', value_name='variant')
    df_var_long = df_var_long.rename(columns={'index': 'cell_id'})

    df_long = pd.merge(df_total_long, df_var_long, on=['cell_id', 'mutation'])

    df_long = df_long[df_long['total'] != 0]
    df_long['cell_id'] = df_long['cell_id'].str[1:].astype(int)
    df_long['mutation'] = df_long['mutation'].str[1:].astype(int)

    df_mut_bin.index.name = "mutation"
    df_mut_bin_unq = df_mut_bin.drop_duplicates(subset=["bin"], keep="first")
    df_mut_bin_unq['bin_id'] = np.arange(len(df_mut_bin_unq))

    df_mut_bin = pd.merge(df_mut_bin, df_mut_bin_unq, on='bin', how='left')
    df_mut_bin.index.name = "mutation"

    df_long = pd.merge(df_long, df_mut_bin, on="mutation", how='left')

    os.makedirs(output_prefix, exist_ok=True)       

    pharming_read_counts = pd.DataFrame({
        'segment': df_long['bin_id'],
        'mutation': df_long['mutation'],
        'cell': df_long['cell_id'],
        'varReads': df_long['variant'],
        'totReads': df_long['total'],
    })
    pharming_read_counts.to_csv(f"{output_prefix}/read_counts.tsv", index=False, sep='\t')

    df_copy_numbers = pd.read_parquet(f"{input_prefix}_copy_numbers.parquet")
    df_copy_numbers.index = df_copy_numbers.index.str[1:].astype(int)
    df_copy_numbers.columns = df_copy_numbers.columns.str[1:].astype(int)
    df_copy_numbers = df_copy_numbers.loc[:, df_mut_bin_unq.index]

    df_copy_numbers_long = df_copy_numbers.reset_index().melt(id_vars='index', var_name='mutation', value_name='copy')
    df_copy_numbers_long.rename(columns={'index': 'cell'}, inplace=True)

    df_copy_numbers_long_cn = pd.merge(df_copy_numbers_long, df_mut_bin_unq, on='mutation', how='left')

    pharming_copy_numbers = pd.DataFrame({
        'segment': df_copy_numbers_long_cn['bin_id'],
        'cell': df_copy_numbers_long_cn['cell'],
        'copiesX': df_copy_numbers_long_cn['copy'],
        'copiesY': [0] * len(df_copy_numbers_long_cn),
    })
    pharming_copy_numbers.to_csv(f"{output_prefix}/copy_numbers.csv", index=False)

    unique_segments = df_mut_bin_unq['bin_id'].unique()
    unique_segments = sorted(unique_segments)

    with open(f"{output_prefix}/segments.txt", "w") as f:
        for seg in unique_segments:
            f.write(f"{seg}\n")



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='input prefix', default='sample')
    parser.add_argument('-o', type=str, help='output prefix', default = 'scope_input')
    parser.add_argument('-v', action='store_true', default=False)
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    
    main(args)
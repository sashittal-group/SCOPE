import argparse
import sys
import pandas as pd
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

    df_read_counts.columns.name = "mutation"
    df_variant_counts.columns.name = "mutation"

    df_total_long = df_read_counts.stack().rename("total").reset_index()
    df_var_long = df_variant_counts.stack().rename("variant").reset_index()

    df_long = pd.merge(df_total_long, df_var_long, on='mutation', how='left')
    df_long = df_long.rename(columns={'index': 'cell_id'})
    df_long = df_long[df_long['total'] != 0]

    df_mut_bin["mutation"] = "m" + df_mut_bin.index.astype(str)

    phertilizer_snv_counts = pd.DataFrame({
        'chr': ['X'] * df_long.shape[0],
        'mutation': df_long['mutation'],
        'cell_id': df_long['cell_id'],
        'alternate': ['X'] * df_long.shape[0],
        'var_count': df_long['variant'],
        'total_count': df_long['total']
    })
    phertilizer_snv_counts.to_csv(f"{output_prefix}/snv_counts.tsv", index=False, header=False, sep='\t')

    df_long = pd.merge(df_long, df_mut_bin, on="mutation", how='left')

    df_binned_counts = (
        df_long.rename(columns={"cell_id": "cell"})
        .pivot(index="cell", columns="bin", values="total")
        .fillna(0)
        .astype(int)
    )
    df_binned_counts.to_csv(f"{output_prefix}/binned_read_counts.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='input prefix', default='sample')
    parser.add_argument('-o', type=str, help='output prefix', default = 'scope_input')
    parser.add_argument('-v', action='store_true', default=False)
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    
    main(args)

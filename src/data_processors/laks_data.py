import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class Laks():

    def __init__(self, DATA_DIR):
        self.df_clone_snvs = pd.read_csv(f"{DATA_DIR}/ov2295_clone_snvs.csv.gz", low_memory=False)
        self.df_clone_cn = pd.read_csv(f"{DATA_DIR}/ov2295_clone_cn.csv.gz")

        bin_sizes = (self.df_clone_cn["end"] - self.df_clone_cn["start"] + 1).unique()
        if len(bin_sizes) > 1:
            raise ValueError(f"There are multiple ({bin_sizes}) bin sizes. There should be only 1 bin size.")

        self.bin_size = bin_sizes

        self.df_clone_snvs["bin_start"] = self.df_clone_snvs["coord"] // self.bin_size * self.bin_size + 1
        self.df_clone_snvs_cn = pd.merge(
            self.df_clone_snvs,
            self.df_clone_cn,
            left_on=["clone_id", "chrom", "bin_start"],
            right_on=["clone_id", "chr", "start"],
            how='left',
        )
        
        df = self.df_clone_snvs_cn
        df["mutation"] = df["chrom"].astype(str) + ":" + df["coord"].astype(str) + ":" + df["ref"].astype(str) + ":" + df["alt"].astype(str)
        
        self.total_table = self.df_clone_snvs_cn.pivot_table(
            index='clone_id',
            columns='mutation',
            values='total_counts',
            aggfunc='sum',
            fill_value=0
        )
        self.alt_table = self.df_clone_snvs_cn.pivot_table(
            index='clone_id',
            columns='mutation',
            values='alt_counts',
            aggfunc='sum',
            fill_value=0
        )
        self.vaf_table = self.alt_table.div(self.total_table).replace(np.nan, 0)

        self.copy_number_table = self.df_clone_snvs_cn.pivot_table(
            index='clone_id',
            columns='mutation',
            values='total_cn',
            aggfunc='mean',
            fill_value=0
        ).astype(int)

        self.F = (self.vaf_table * self.copy_number_table).clip(upper=1)


    def get_cell_fractions(self):
        return self.F


import os 
import pandas as pd


def main(peak_cov_meta_file, store_dir):
    # read file
    peak_cov_meta_df = pd.read_csv(peak_cov_meta_file, index_col=[0,1,2])
    peak_cov_meta_df_length = pd.Series(peak_cov_meta_df.index.get_level_values(2) - peak_cov_meta_df.index.get_level_values(1))
    peak_cov_meta_df_length.index = peak_cov_meta_df.index
    # normalize by enhancer length
    peak_cov_meta_df_norm = peak_cov_meta_df.div(peak_cov_meta_df_length, axis=0)*1000
    # save file
    store_file = os.path.join(store_dir, "meta_cov_norm.csv")
    peak_cov_meta_df_norm.to_csv(store_file, index=True)
    return 


if __name__ == "__main__":
    peak_cov_meta_file = "/data5/deepro/starrseq/main_library/4_quality_control_peaks/data/peak_cov/meta_cov.csv"
    store_dir = "/data5/deepro/starrseq/main_library/9_rnaseq_comparison/data"
    main(peak_cov_meta_file, store_dir)
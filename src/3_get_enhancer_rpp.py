import os 
import pandas as pd


def main(peak_cov_meta_file, store_dir):
    # read file
    peak_cov_meta_df = pd.read_csv(peak_cov_meta_file, index_col=[0,1,2])
    # drop duplicates based on index - fixes an issue about missing regions
    peak_cov_meta_df = peak_cov_meta_df.loc[~peak_cov_meta_df.index.duplicated(keep="first"), :]
    # normalize by input plasmid count
    for i in range(1, 4):
        ri_cols = [c for c in peak_cov_meta_df.columns if ((c[-2:]==f"R{i}") and (c[:2] != "IN"))]
        peak_cov_meta_df.loc[:, ri_cols] = peak_cov_meta_df.loc[:, ri_cols].div(peak_cov_meta_df[f'IN_R{i}'], axis=0)
    # save file
    store_file = os.path.join(store_dir, "meta_rpp.csv")
    peak_cov_meta_df.iloc[:, 3:].to_csv(store_file, index=True)
    return 


if __name__ == "__main__":
    peak_normcov_meta_file = "/data5/deepro/starrseq/main_library/9_rnaseq_comparison/data/meta_cov_norm.csv"
    store_dir = "/data5/deepro/starrseq/main_library/9_rnaseq_comparison/data"
    main(peak_normcov_meta_file, store_dir)

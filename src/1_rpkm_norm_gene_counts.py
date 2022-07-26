import os 
import pandas as pd


def main(gene_counts_file, gene_annot_file, store_dir):
    # read file
    gene_counts_meta_df = pd.read_csv(gene_counts_file)
    gene_annot_df = pd.read_csv(gene_annot_file)
    # check that their ensemble ids are the same and in order
    assert (gene_annot_df.iloc[:, 0] == gene_counts_meta_df.iloc[:, 0]).all()
    # calculate gene length
    gene_annot_df["length"] = gene_annot_df.end - gene_annot_df.start
    gene_counts_meta_df = gene_counts_meta_df.set_index(list(gene_counts_meta_df.columns[:3]))
    gene_annot_df.index = gene_counts_meta_df.index
    # normalize by gene length
    gene_counts_meta_df_norm = gene_counts_meta_df.div(gene_annot_df.length, axis=0)*1000
    # normalize by library size
    gene_counts_meta_df_norm = gene_counts_meta_df_norm.div(gene_counts_meta_df.sum(), axis=1)*1e6
    # save file
    store_file = os.path.join(store_dir, "meta_counts_norm.csv")
    gene_counts_meta_df_norm.to_csv(store_file, index=True)
    return 


if __name__ == "__main__":
    gene_counts_meta_file = "/data5/deepro/starrseq/main_library/7_rnaseq/data/de/gene_counts.csv"
    gene_annot_file = "/data5/deepro/starrseq/main_library/9_rnaseq_comparison/data/gene_annotations.csv"
    store_dir = "/data5/deepro/starrseq/main_library/9_rnaseq_comparison/data"
    main(gene_counts_meta_file, gene_annot_file, store_dir)

import os
import pandas as pd


def main(gtf_file, store_dir):
    # read gene annotation file
    gene_annot_df = pd.read_csv(gtf_file, comment="#", sep="\t", header=None, low_memory=False)
    # only keep genes
    gene_annot_df = gene_annot_df.loc[gene_annot_df.iloc[:,2]=="gene"]
    # take out ensemble id and gene name from the information column
    gene_annot_df["chrm"] = gene_annot_df.iloc[:,0]
    gene_annot_df["start"] = gene_annot_df.iloc[:,3]
    gene_annot_df["end"] = gene_annot_df.iloc[:, 4]
    gene_annot_df["ensemble_id"] = [x[0].split()[1].replace('"', "") for x in gene_annot_df.iloc[:, 8].str.split(";") ]
    gene_annot_df["gene_name"] = [x[2].split()[1].replace('"', "") for x in gene_annot_df.iloc[:, 8].str.split(";") ]
    # save file
    store_file = os.path.join(store_dir, "gene_annotations.csv")
    gene_annot_df.loc[:, ["ensemble_id", "gene_name", "chrm", "start", "end"]].to_csv(store_file, index=False)
    return


if __name__ == "__main__":
    gtf_file = "/data5/deepro/genomes/hg38/gencode.v38.annotation.gtf" # from http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
    store_dir = "/data5/deepro/starrseq/main_library/9_rnaseq_comparison/data"
    main(gtf_file, store_dir)
import os
import pandas as pd
from scipy import stats


###################
# filepath parser #
###################

def get_egmap_file(egmap_dir, lib_short, method, filebase):
    return os.path.join(egmap_dir, lib_short, method, f"{filebase}.tsv")

def get_peak_file(peak_dir, lib_short, filebase):
    return os.path.join(peak_dir, lib_short, f"{filebase}.bed")

def get_rnaseq_file(rnaseq_dir, lib_short):
    return os.path.join(rnaseq_dir, f"{lib_short}_vs_CC.tsv")

def get_egmap_df(egmap_file):
    egmap_df = pd.read_csv(egmap_file, sep="\t", skiprows=3, skipfooter=10, engine="python")
    egmap_df = egmap_df.loc[egmap_df.iloc[:, 0] == "Ensembl Genes"]
    return egmap_df

def get_peak_df(peak_file):
    peak_df = pd.read_csv(peak_file, sep="\t", header=None)
    return peak_df

def get_rnaseq_df(rnaseq_de_file):
    rnaseq_de_df = pd.read_csv(rnaseq_de_file, sep="\t")
    return rnaseq_de_df

def get_peak_cov_files(peak_cov_dir, lib_short):
    in_peakcov_file = os.path.join(peak_cov_dir, lib_short, "input.csv")
    out_peakcov_file = os.path.join(peak_cov_dir, lib_short, "output.csv")
    return in_peakcov_file, out_peakcov_file

################
# t-test stats #
################

def get_t_stats_and_percent_de_genes(rnaseq_df, egmap_df):
    # find enhancer mapped genes that were present in rnaseq data
    rnaseq_selected = rnaseq_df.loc[rnaseq_df.gene_symbol.str.lower().isin(egmap_df.Genes.str.lower())]
    rnaseq_others = rnaseq_df.loc[~rnaseq_df.gene_symbol.str.lower().isin(egmap_df.Genes.str.lower())]
    # calculate percent de genes
    percent_de = len(rnaseq_selected.loc[rnaseq_selected.FDR<0.05])*100/len(rnaseq_selected)
    # calculate t-test of logfc
    test_pval = stats.ttest_ind(rnaseq_others.logFC, rnaseq_selected.logFC).pvalue
    return percent_de, test_pval

import datatable
import matplotlib
import numpy
import os
import pandas
import scipy
import seaborn 
import heatmap_grammar

# critical parameter setting!
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams['figure.figsize'] = 15, 15
pandas.options.mode.chained_assignment = None  # default='warn'

# setup current direcotry
data_file_directory = "~/Documents/Bioinformatics_analysis/eIF4G-analysis/eIF4G_data"
output_directory = "~/Documents/Bioinformatics_analysis/eIF4G-analysis/eIF4G_output"

## acquire the data
TCGA_GTEX_RNAseq = datatable.fread(os.path.join(os.path.expanduser(data_file_directory),
"TcgaTargetGtex_RSEM_Hugo_norm_count"), header=True).to_pandas()
TCGA_GTEX_RNAseq = TCGA_GTEX_RNAseq.set_index(['sample']).T


TCGA_CNV = datatable.fread(os.path.join(os.path.expanduser(output_directory),"ProcessedData","TCGA_CNV.csv"), header=True).to_pandas()
TCGA_CNV = TCGA_CNV.set_index([TCGA_CNV.columns[0]])
## too few homdel samples for correlation analysis, combine homdel and del samples for correlation
TCGA_CNV1 = TCGA_CNV.replace(['HOMDEL'], 'DEL')


def get_EIF_CNV_RNAseq_combined_data(gene_name):
    TCGA_CNV_gene = TCGA_CNV1[gene_name]
    col_name = str(gene_name)+str('_CNV')
    # merge RNAseq and CNV data based on sample names
    df = pandas.merge(TCGA_GTEX_RNAseq, 
                      TCGA_CNV_gene, 
                      how='left',
                      left_index=True, 
                      right_index=True,
                      suffixes=("", "_CNV"))
    # annotate the CNV status from GTEX data as "NAT"
    df.loc[df.index.str.contains('GTEX'), col_name] = "NAT"
    # remove samples without CNV status, such as target samples 
    df = df[df[col_name].notna()]
    return (df)

#TCGA_GTEX_RNAseq_CNV_EIF4G1 = get_EIF_CNV_RNAseq_combined_data("EIF4G1")

def rnaseq_cnv_corr_sum(df, gene_name):
    y = pandas.DataFrame()
    col_name = str(gene_name)+str('_CNV')
    CNV_status = pandas.Series(["AMP", "DUP", "DIPLOID", "DEL", 'NAT'])
    for CNV in CNV_status:
        print(CNV)
        # select rows based on the CNV status
        df1 = df.loc[df[col_name] == CNV]
        # remove the CNV status column to facilitate the correlation analysis
        df1 = df1.drop(col_name, axis=1)
        x = df1.corrwith(df1[gene_name])
        x.name = CNV
        y = pandas.concat([y, x], axis=1)
    y.dropna(inplace=True)
    y.index.names = ["Gene_Symbol"]
    # pandas.DataFrame.ge: get Greater than or equal to of dataframe and other
    # pandas.DataFrame.le: get Less than or equal to of dataframe and other
    z = y[y.ge(0.4).any(axis=1)|y.le(-0.4).any(axis=1)]
    return (y, z)

EIF4G1_CNV_COR_RNAseq, EIF4G1_CNV_COR_RNAseq_sig = rnaseq_cnv_corr_sum(df = get_EIF_CNV_RNAseq_combined_data("EIF4G1"), gene_name = "EIF4G1")
EIF4G1_CNV_COR_RNAseq_sig.to_csv(os.path.join(os.path.expanduser(output_directory),"ProcessedData","EIF4G1_CNV_COR_RNAseq_sig.csv"))

EIF4A2_CNV_COR_RNAseq, EIF4A2_CNV_COR_RNAseq_sig = rnaseq_cnv_corr_sum(df = get_EIF_CNV_RNAseq_combined_data("EIF4A2"), gene_name = "EIF4A2")
EIF4A2_CNV_COR_RNAseq_sig.to_csv(os.path.join(os.path.expanduser(output_directory),"ProcessedData","EIF4A2_CNV_COR_RNAseq_sig.csv"))

EIF3E_CNV_COR_RNAseq, EIF3E_CNV_COR_RNAseq_sig = rnaseq_cnv_corr_sum(df = get_EIF_CNV_RNAseq_combined_data("EIF3E"), gene_name = "EIF3E")
EIF3E_CNV_COR_RNAseq_sig.to_csv(os.path.join(os.path.expanduser(output_directory),"ProcessedData","EIF3E_CNV_COR_RNAseq_sig.csv"))

EIF3H_CNV_COR_RNAseq, EIF3H_CNV_COR_RNAseq_sig = rnaseq_cnv_corr_sum(df = get_EIF_CNV_RNAseq_combined_data("EIF3H"), gene_name = "EIF3H")
EIF3H_CNV_COR_RNAseq_sig.to_csv(os.path.join(os.path.expanduser(output_directory),"ProcessedData","EIF3H_CNV_COR_RNAseq_sig.csv"))



h = seaborn.clustermap(
    EIF4G1_CNV_COR_RNAseq_sig,
    figsize=(10, 10),
    #method="centroid",
    #metric="euclidean",
    tree_kws=dict(linewidths=0.5, colors=(0.2, 0.2, 0.4)),
    cmap="coolwarm",
    #cbar_pos=(.2, .2, .03, .4),
    yticklabels=False
)
matplotlib.pyplot.show()



## this function calculate pvalue, but not useful for cor analysis

def rnaseq_cnv_corr(df, gene_name):
    x = pandas.DataFrame()
    col_name = str(gene_name)+str('_CNV')
    CNV_status = pandas.Series(["AMP", "DUP", "DIPLOID", "DEL", 'NAT'])
    # list all gene names from the column names
    gene_list = list(df.columns.values)
    gene_list.remove('EIF4G1_CNV')
    for CNV in CNV_status:
        df1 = df.loc[df[col_name] == CNV]
        y = pandas.DataFrame()
        for gene in gene_list:
            # select rows based on the CNV status
            r, p = scipy.stats.pearsonr(df1[gene_name], df1[gene])
            a = [[r,p]]
            t = pandas.DataFrame(a, columns =[str(CNV), str(CNV)+"pvalue"], index=[gene])
            y = pandas.concat([y,t],axis=0)
        #y.dropna(inplace=True)
        #z = y.loc[y['pvalue'] < 0.05
        x = pandas.concat([y, x], axis=1)
    z = x[x.ge(0.4).any(axis=1)|x.le(-0.4).any(axis=1)]
    return (x)

rnaseq_cor =  rnaseq_cnv_corr(df = TCGA_GTEX_RNAseq_CNV_EIF4G1, 
                              gene_name = "EIF4G1")


rnaseq_cor_new = rnaseq_cor[(rnaseq_cor['AMP'] > 0.4) | (rnaseq_cor['DUP'] > 0.4) | 
                            (rnaseq_cor['DIPLOID'] > 0.4) | (rnaseq_cor['DEL'] > 0.4) | 
                            (rnaseq_cor['NAT'] > 0.4) |
                            (rnaseq_cor['AMP'] < -0.4) | (rnaseq_cor['DUP'] < -0.4) | 
                            (rnaseq_cor['DIPLOID'] < -0.4) | (rnaseq_cor['DEL'] < -0.4) | 
                            (rnaseq_cor['NAT'] < -0.4)]
rnaseq_cor_new1 = rnaseq_cor_new.dropna()

h = seaborn.clustermap(
    rnaseq_cor_new1[['AMP', 'DUP', 'DIPLOID', 'DEL', 'NAT']],
    figsize=(10, 10),
    #method="centroid",
    #metric="euclidean",
    tree_kws=dict(linewidths=0.5, colors=(0.2, 0.2, 0.4)),
    cmap="coolwarm",
    #cbar_pos=(.2, .2, .03, .4),
    yticklabels=False
)
matplotlib.pyplot.show()




















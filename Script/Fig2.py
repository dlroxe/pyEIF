import datatable
import matplotlib
import numpy
import os
import pandas
import scipy
import seaborn 
import sklearn
import umap
import statannot
import umap.plot
import sklearn.preprocessing

#import heatmap_grammar

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
                                                "TcgaTargetGtex_RSEM_Hugo_norm_count"), 
                                   header=True).to_pandas()
TCGA_GTEX_RNAseq = TCGA_GTEX_RNAseq.set_index(['sample']).T


TCGA_CNV = datatable.fread(os.path.join(os.path.expanduser(output_directory),
                                        "ProcessedData","TCGA_CNV.csv"), 
                           header=True).to_pandas()
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
    # remove all row when gene value is zero
    df = df[df[gene_name] != 0]
    # annotate the CNV status from GTEX data as "NAT"
    df.loc[df.index.str.contains('GTEX'), col_name] = "NAT"
    # remove samples without CNV status, such as target samples 
    df = df[df[col_name].notna()]
    return (df)

#TCGA_GTEX_RNAseq_CNV_EIF4G1 = get_EIF_CNV_RNAseq_combined_data(gene_name = "EIF4G1")
EIF4G1_CNV_RNA = get_EIF_CNV_RNAseq_combined_data("EIF4G1")

def plot_RNAseq_CNV_status(df, gene_name):
    #seaborn.set(rc={'figure.figsize':(12,10)})
    #seaborn.set(font_scale=2)
    dft = df.groupby(['EIF4G1_CNV'])['EIF4G1_CNV'].count()
    order = ["NAT", "AMP", "DUP", "DIPLOID", "DEL"]
    matplotlib.pyplot.clf()
    fig, ax = matplotlib.pyplot.subplots(figsize=(10, 10))
    seaborn.set_style("ticks")
    ax = seaborn.violinplot(x='EIF4G1_CNV', y=gene_name, data=df, order=order)
    statannot.add_stat_annotation(ax, data=df, 
                                  x='EIF4G1_CNV', y=gene_name, order=order,
                                  box_pairs=[("NAT", "AMP"), ("NAT", "DUP"), ("NAT", "DIPLOID"), 
                                             ("AMP", "DUP"), ("AMP", "DIPLOID"), ("DUP", "DIPLOID")],
                                  test='Mann-Whitney', text_format='star', loc='inside', verbose=2)
    ax.set_xticklabels([("Normal"+ '\nn = ' + str(dft["NAT"])),
                        ("Amp"+ '\nn = ' + str(dft["AMP"])),
                        ("Dup"+ '\nn = ' + str(dft["DUP"])),
                        ("Diploid"+ '\nn = ' + str(dft["DIPLOID"])),
                        ("Del"+ '\nn = ' + str(dft["DEL"]))])
    ax.set_xlabel("EIF4G1 CNV status")
    ax.set_ylabel(str(gene_name)+" RNA expression")
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(os.path.join(os.path.expanduser(output_directory), 
                                           "Fig2", (str(gene_name)+'expression_CNV.pdf')),
                              dpi=300,
                              edgecolor='w',
                              #transparent=True,
                              bbox_inches='tight')   
    matplotlib.pyplot.show()


plot_RNAseq_CNV_status(df = EIF4G1_CNV_RNA, gene_name = "BRCA1")
plot_RNAseq_CNV_status(df = EIF4G1_CNV_RNA, gene_name = "RANGAP1")
plot_RNAseq_CNV_status(df = EIF4G1_CNV_RNA, gene_name = "ORC1")
plot_RNAseq_CNV_status(df = EIF4G1_CNV_RNA, gene_name = "CDC20")
plot_RNAseq_CNV_status(df = EIF4G1_CNV_RNA, gene_name = "SREBF2")
plot_RNAseq_CNV_status(df = EIF4G1_CNV_RNA, gene_name = "HMGCR")
plot_RNAseq_CNV_status(df = EIF4G1_CNV_RNA, gene_name = "HMGCS1")

plot_RNAseq_CNV_status(df = EIF4G1_CNV_RNA, gene_name = "EIF4G1")
plot_RNAseq_CNV_status(df = EIF4G1_CNV_RNA, gene_name = "EIF4E")
plot_RNAseq_CNV_status(df = EIF4G1_CNV_RNA, gene_name = "EIF4A2")
plot_RNAseq_CNV_status(df = EIF4G1_CNV_RNA, gene_name = "EIF3E")
plot_RNAseq_CNV_status(df = EIF4G1_CNV_RNA, gene_name = "CENPI")
plot_RNAseq_CNV_status(df = EIF4G1_CNV_RNA, gene_name = "LAMA3")
plot_RNAseq_CNV_status(df = EIF4G1_CNV_RNA, gene_name = "ITGA6")

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

EIF4G1_CNV_COR_RNAseq, EIF4G1_CNV_COR_RNAseq_sig = rnaseq_cnv_corr_sum(df = get_EIF_CNV_RNAseq_combined_data("EIF4G1"), 
                                                                       gene_name = "EIF4G1")
EIF4G1_CNV_COR_RNAseq_sig.to_csv(os.path.join(os.path.expanduser(output_directory),
                                              "ProcessedData",
                                              "EIF4G1_CNV_COR_RNAseq_sig.csv"))

EIF4A2_CNV_COR_RNAseq, EIF4A2_CNV_COR_RNAseq_sig = rnaseq_cnv_corr_sum(df = get_EIF_CNV_RNAseq_combined_data("EIF4A2"), 
                                                                       gene_name = "EIF4A2")
EIF4A2_CNV_COR_RNAseq_sig.to_csv(os.path.join(os.path.expanduser(output_directory),
                                              "ProcessedData",
                                              "EIF4A2_CNV_COR_RNAseq_sig.csv"))

EIF3E_CNV_COR_RNAseq, EIF3E_CNV_COR_RNAseq_sig = rnaseq_cnv_corr_sum(df = get_EIF_CNV_RNAseq_combined_data("EIF3E"), 
                                                                     gene_name = "EIF3E")
EIF3E_CNV_COR_RNAseq_sig.to_csv(os.path.join(os.path.expanduser(output_directory),
                                             "ProcessedData",
                                             "EIF3E_CNV_COR_RNAseq_sig.csv"))

EIF3H_CNV_COR_RNAseq, EIF3H_CNV_COR_RNAseq_sig = rnaseq_cnv_corr_sum(df = get_EIF_CNV_RNAseq_combined_data("EIF3H"), 
                                                                     gene_name = "EIF3H")
EIF3H_CNV_COR_RNAseq_sig.to_csv(os.path.join(os.path.expanduser(output_directory),
                                             "ProcessedData",
                                             "EIF3H_CNV_COR_RNAseq_sig.csv"))


h = seaborn.clustermap(
    EIF4G1_CNV_COR_RNAseq_sig,
    #EIF4G1_CNV_COR_RNAseq_sig.drop(['NAT'], axis=1),
    figsize=(10, 10),
    #method="centroid",
    #metric="euclidean",
    tree_kws=dict(linewidths=0.5, colors=(0.2, 0.2, 0.4)),
    cmap="coolwarm",
    #cbar_pos=(.2, .2, .03, .4),
    yticklabels=False
)
matplotlib.pyplot.show()


### umap analysis on cluster genes from heatmap
EIF4G1_CNV_RNA_COR_cluster1 = pandas.read_csv(os.path.join(os.path.expanduser(output_directory),
                                                           "ProcessedData", 
                                                           "CNV_RNAseq_COR_cluster1.csv")) 
EIF4G1_CNV_RNA_COR_cluster2 = pandas.read_csv(os.path.join(os.path.expanduser(output_directory),
                                                           "ProcessedData", 
                                                           "CNV_RNAseq_COR_cluster2.csv"))
EIF4G1_CNV_RNA_COR_cluster3 = pandas.read_csv(os.path.join(os.path.expanduser(output_directory),
                                                           "ProcessedData", 
                                                           "CNV_RNAseq_COR_cluster3.csv")) 
EIF4G1_CNV_RNA_COR_cluster4 = pandas.read_csv(os.path.join(os.path.expanduser(output_directory),
                                                           "ProcessedData", 
                                                           "CNV_RNAseq_COR_cluster4.csv")) 
EIF4G1_CNV_RNA_COR_cluster5 = pandas.read_csv(os.path.join(os.path.expanduser(output_directory),
                                                           "ProcessedData", 
                                                           "CNV_RNAseq_COR_cluster5.csv")) 

def get_RNAse_CNV_cluster(df, gene_name):
    cluster_RNAseq = TCGA_GTEX_RNAseq[df["gene"]]
    #gene_name = "EIF4G1"
    TCGA_CNV_gene = TCGA_CNV1[gene_name]
    col_name = str(gene_name)+str('_CNV')
    TCGA_CNV_gene.rename((str(gene_name)+'_CNV'), inplace = True)

    # merge RNAseq and CNV data based on sample names
    df = pandas.merge(cluster_RNAseq, 
                      TCGA_CNV_gene, 
                      how='left',
                      left_index=True, 
                      right_index=True,
                      suffixes=("", "_CNV"))
    # annotate the CNV status from GTEX data as "NAT"
    df.loc[df.index.str.contains('GTEX'), col_name] = "NAT"
    #df["EIF4G1_CNV"]
    # remove samples without CNV status, such as TARGET samples 
    df = df[df[col_name].notna()]
    df1 = df.drop([col_name], axis=1).values
    return(df, df1)

df, df1 = get_RNAse_CNV_cluster(EIF4G1_CNV_RNA_COR_cluster5, "EIF4G1")

### UMAP =======================================================
def plot_umap(df, df1, cluster):
    scaled_data = sklearn.preprocessing.StandardScaler().fit_transform(df1)
    embedding = umap.UMAP(n_neighbors=10, min_dist=0.3,).fit_transform(scaled_data)
    embedding.shape
    fig, ax = matplotlib.pyplot.subplots(figsize=(10,10))
    seaborn.set_style("ticks")
    ax = seaborn.scatterplot(x=embedding[:, 0],
                             y=embedding[:, 1], 
                             hue=df.EIF4G1_CNV,
                             hue_order = ['AMP', 'DUP', 'DIPLOID','DEL','NAT'],
                             legend='full', 
                             palette= dict({'AMP':'red',
                                               'DUP':'orange',
                                               'DIPLOID': 'grey',
                                               'DEL': 'blue',
                                               'NAT': 'green'}), 
                             alpha = 1/5,
                             edgecolor="none",
                             s=10)
    ax.legend(loc='lower left')
    ax.set_title('UMAP with '+ str(cluster)+ 'gene')
    ax.set_xlabel("UMAP_1")
    ax.set_ylabel("UMAP_2")
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(os.path.join(os.path.expanduser(output_directory), 
                                           "Fig2", (str(cluster)+' UMAP.pdf')),
                              dpi=300, alpha = 1/5,
                              edgecolor='none',
                              #transparent=True,
                              bbox_inches='tight')   
plot_umap(df, df1,"cluster5")

### UMAP =======================================================

mapper = umap.UMAP().fit(df1)
umap.plot.points(mapper)
umap.plot.points(mapper, labels=df.EIF4G1_CNV)

## PCA  ======================================
from sklearn.decomposition import PCA
penguins_data=df.select_dtypes(numpy.number)
penguins_data.head()
penguins_info=df.select_dtypes(exclude='float')
penguins_info.head()
cnv=penguins_info.EIF4G1_CNV.tolist()
pca = sklearn.decomposition.PCA(n_components=4)
penguins_pca= pca.fit_transform(penguins_data)
pc_df = pandas.DataFrame(data = penguins_pca , 
        columns = ['PC1', 'PC2','PC3', 'PC4'])
pc_df.head()

pc_df['CNV']=cnv
pc_df_amp = pc_df.loc[pc_df['CNV'] == "AMP"]
pc_df_dup = pc_df.loc[pc_df['CNV'] == "DUP"]
pc_df_dip = pc_df.loc[pc_df['CNV'] == "DIPLOID"]
pc_df_del = pc_df.loc[pc_df['CNV'] == "DEL"]
pc_df_nat = pc_df.loc[pc_df['CNV'] == "NAT"]

pc_df.head()
pca.explained_variance_ratio_
seaborn.set_theme(style='white')
color_dict = dict({'DIPLOID':'grey',
                   'NAT':'green',
                   'DUP': 'orange',
                   'AMP': 'red',
                   'DEL': 'blue'})


matplotlib.pyplot.figure(figsize=(12,10))
with seaborn.plotting_context("notebook",font_scale=1.25):
    seaborn.scatterplot(x="PC1", y="PC2",
                    data=pc_df, 
                    hue="CNV",
                    palette=color_dict,
                    alpha = 1/5,
                    s=10)
# control x and y limits
matplotlib.pyplot.ylim(-100, 100)
matplotlib.pyplot.xlim(-100, 80)
matplotlib.pyplot.show()


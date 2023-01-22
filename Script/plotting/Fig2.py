import datatable
import matplotlib
import numpy
import os
import pandas
import polars
import seaborn
import sklearn
import umap
import umap.plot
import statannot
import sklearn.preprocessing

from absl import app
from absl import flags
from absl import logging
from typing import Tuple

# critical parameter setting!
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams['figure.figsize'] = 15, 15
pandas.options.mode.chained_assignment = 'raise'  # default='warn'

# TODO(dlroxe): Try to define FLAGs in a single place, and reference them from
#               everywhere that needs them.
FLAGS = flags.FLAGS
flags.DEFINE_string('data_directory',
                    # '~/Documents/Bioinformatics_analysis/eIF4G-analysis/eIF4G_data',
                    os.path.join('~', 'Desktop', 'pyeif_data'),
                    'parent dir for data files')
flags.DEFINE_string('output_directory',
                    # '~/Documents/Bioinformatics_analysis/eIF4G-analysis/eIF4G_output',
                    os.path.join('~', 'Desktop', 'pyeif_output'),
                    'parent dir for output')

_SCATTERPLOT_COLOR_DICT = {
  'AMP': 'red',
  'DEL': 'blue',
  'DIPLOID': 'grey',
  'DUP': 'orange',
  'NAT': 'green',
}


def get_EIF_CNV_RNAseq_combined_data(TCGA_CNV1, TCGA_GTEX_RNAseq, gene_name):
  TCGA_CNV_gene = TCGA_CNV1[gene_name]
  col_name = f'{gene_name}_CNV'
  # merge RNAseq and CNV data based on sample names
  df = pandas.merge(
    TCGA_GTEX_RNAseq,
    TCGA_CNV_gene,
    how='left',
    left_index=True,
    right_index=True,
    suffixes=("", "_CNV"),
  )
  # remove all row when gene value is zero
  df = df[df[gene_name] != 0]
  # annotate the CNV status from GTEX data as "NAT"
  df.loc[df.index.str.contains('GTEX'), col_name] = "NAT"
  # remove samples without CNV status, such as target samples
  df = df[df[col_name].notna()]
  return df


def plot_RNAseq_CNV_status(df, gene_name):
  # seaborn.set(rc={'figure.figsize':(12,10)})
  # seaborn.set(font_scale=2)
  dft = df.groupby(['EIF4G1_CNV'])['EIF4G1_CNV'].count()
  order = ["NAT", "AMP", "DUP", "DIPLOID", "DEL"]
  matplotlib.pyplot.clf()
  fig, ax = matplotlib.pyplot.subplots(figsize=(10, 10))
  seaborn.set_style("ticks")
  ax = seaborn.violinplot(x='EIF4G1_CNV', y=gene_name, data=df, order=order)
  statannot.add_stat_annotation(
    ax,
    box_pairs=[
      ("NAT", "AMP"),
      ("NAT", "DUP"),
      ("NAT", "DIPLOID"),
      ("AMP", "DUP"),
      ("AMP", "DIPLOID"),
      ("DUP", "DIPLOID"),
    ],
    data=df,
    loc='inside',
    order=order,
    test='Mann-Whitney',
    text_format='star',
    verbose=2,
    x='EIF4G1_CNV',
    y=gene_name,
  )
  ax.set_xticklabels([
    f'Normal\nn = {dft["NAT"]}',
    f'Amp\nn = {dft["AMP"]}',
    f'Dup\nn = {dft["DUP"]}',
    f'Diploid\nn = {dft["DIPLOID"]}',
    f'Del\nn = {dft["DEL"]}',
  ])
  ax.set_xlabel("EIF4G1 CNV status")
  ax.set_ylabel(f'{gene_name} RNA expression')
  matplotlib.pyplot.tight_layout()
  matplotlib.pyplot.savefig(
    _abspath(
      FLAGS.output_directory, 'Fig2', f'{gene_name}-expression-CNV.pdf'),
    bbox_inches='tight',
    dpi=300,
    edgecolor='w',
    # transparent=True,
  )
  matplotlib.pyplot.show()


def rnaseq_cnv_corr_sum(df, gene_name):
  y = pandas.DataFrame()
  col_name = f'{gene_name}_CNV'
  for CNV in pandas.Series(["AMP", "DUP", "DIPLOID", "DEL", 'NAT']):
    logging.info('rnaseq_cnv_corr_sum() considering %s:%s', gene_name, CNV)
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
  z = y[y.ge(0.4).any(axis=1) | y.le(-0.4).any(axis=1)]
  return y, z


def get_RNAse_CNV_cluster(TCGA_CNV1, TCGA_GTEX_RNAseq, df, gene_name):
  cluster_RNAseq = TCGA_GTEX_RNAseq[df["gene"]]
  # gene_name = "EIF4G1"
  TCGA_CNV_gene = TCGA_CNV1[gene_name]
  col_name = f'{gene_name}_CNV'
  TCGA_CNV_gene.rename(col_name, inplace=True)

  # merge RNAseq and CNV data based on sample names
  df = pandas.merge(
    cluster_RNAseq,
    TCGA_CNV_gene,
    how='left',
    left_index=True,
    right_index=True,
    suffixes=("", "_CNV"),
  )
  # annotate the CNV status from GTEX data as "NAT"
  df.loc[df.index.str.contains('GTEX'), col_name] = "NAT"
  # df["EIF4G1_CNV"]
  # remove samples without CNV status, such as TARGET samples
  df = df[df[col_name].notna()]
  df1 = df.drop([col_name], axis=1).values
  return df, df1


### UMAP =======================================================
def plot_umap(df, df1, cluster):
  scaled_data = sklearn.preprocessing.StandardScaler().fit_transform(df1)
  # Per https://umap-learn.readthedocs.io/en/latest/reproducibility.html,
  # UMAP output can be made stable/reproducible by specifying a value for
  # random_state.  This reduces performance because it constrains the
  # amount of parallel execution (threading) that the UMAP library may
  # employ.
  #
  # For those unfamiliar, an approximate primer:  'Threads' permit a single
  # program to divide work among multiple CPU cores.  When this is done, there
  # are no general guarantees that the cores will always execute the same
  # instructions, in the same order relative to each other, in the same
  # amount of time, from one program execution to the next.  That is, threads
  # offer a trade-off between performance and reproducibility.  From a technical
  # standpoint, there are various ways of managing thread execution to trade
  # back some performance in exchange for reproducibility of specific outcomes.
  #
  # At any rate, in the case of UMAP, 1) it uses threads "under the covers"
  # and 2) it offers the random_state "knob" to trade back some performance
  # for reproducible plots.  For the particular task at hand here, no
  # *significant* loss has been observed.
  #
  # However, a "suggested waiver" on that web page points out that such
  # reproducibility comes with the risk of producing a result that is unduly
  # influenced by the random_state value.
  # TODO(dlroxe): It might be useful to produce plots both with and without
  #               a random seed.  That way, others can exactly reproduce the
  #               reproducible plot; but they can also generate and
  #               qualitatively compare a non-reproducible plot without such
  #               artifact risk.
  embedding = umap.UMAP(
    n_neighbors=10, min_dist=0.3, random_state=42).fit_transform(scaled_data)
  logging.info('UMAP embedding shape: %s', embedding.shape)
  fig, ax = matplotlib.pyplot.subplots(figsize=(10, 10))
  seaborn.set_style("ticks")
  ax = seaborn.scatterplot(
    # alpha=1 / 5,
    edgecolor="none",
    hue=df.EIF4G1_CNV,
    hue_order=['AMP', 'DUP', 'DIPLOID', 'DEL', 'NAT'],
    legend='full',
    palette=_SCATTERPLOT_COLOR_DICT,
    s=10,
    x=embedding[:, 0],
    y=embedding[:, 1],
  )
  ax.legend(loc='lower left')
  ax.set_title(f'UMAP with {cluster} gene')
  ax.set_xlabel("UMAP_1")
  ax.set_ylabel("UMAP_2")
  matplotlib.pyplot.tight_layout()
  # N.B. The '-fixed' suffix below is meant to signify a reproducible plot,
  #      i.e. one made with a specified random_state value.
  matplotlib.pyplot.savefig(
    _abspath(FLAGS.output_directory, "Fig2", f'{cluster}-UMAP-fixed.pdf'),
    # alpha=1 / 5,
    bbox_inches='tight',
    dpi=300,
    edgecolor='none',
    # transparent=True,
  )


def _read_cluster_csv(filename: str) -> pandas.DataFrame:
  return datatable.fread(
    _abspath(FLAGS.output_directory, "ProcessedData", filename)).to_pandas()


# TODO(dlroxe): This function is cropping up in several places.  Find a
#               One True Home for it.  Note, this version is a little nicer
#               than the others, because it takes varargs and internalizes
#               an os.path.join() employed by callers elsewhere.
def _abspath(*args):
  # What an absurd incantation to resolve "~"; but OK. Thanks, StackOverflow.
  return os.path.abspath(
    os.path.expanduser(os.path.expandvars(os.path.join(*args))))


def main(unused_argv):
  logging.info('commencing work for Fig2')
  ## acquire the data
  TCGA_GTEX_RNAseq = datatable.fread(
    _abspath(FLAGS.data_directory, "TcgaTargetGtex_RSEM_Hugo_norm_count"),
    header=True).to_pandas()
  TCGA_GTEX_RNAseq = TCGA_GTEX_RNAseq.set_index(['sample']).T

  tcga_cnv_path = _abspath(
    FLAGS.output_directory, 'ProcessedData', 'TCGA-CNV-thresholds.csv')
  logging.info(f'reading: {tcga_cnv_path}')
  TCGA_CNV = datatable.fread(tcga_cnv_path, header=True).to_pandas()

  TCGA_CNV.set_index([TCGA_CNV.columns[0]], inplace=True)
  # too few homdel samples for correlation analysis;
  # combine homdel and del samples for correlation
  TCGA_CNV1 = TCGA_CNV.replace(['HOMDEL'], 'DEL')

  # TCGA_GTEX_RNAseq_CNV_EIF4G1 = get_EIF_CNV_RNAseq_combined_data(gene_name = "EIF4G1")
  EIF4G1_CNV_RNA = get_EIF_CNV_RNAseq_combined_data(
    TCGA_CNV1, TCGA_GTEX_RNAseq, "EIF4G1")

  # TODO(dlroxe): Make the list of genes flag-configurable.
  logging.info('commencing CNV-status violin plots')
  for gene in (
      "BRCA1", "RANGAP1", "ORC1", "CDC20", "SREBF2", "HMGCR", "HMGCS1",
      "EIF4G1", "EIF4E", "EIF4A2", "EIF3E", "CENPI", "LAMA3", "ITGA6"):
    plot_RNAseq_CNV_status(df=EIF4G1_CNV_RNA, gene_name=gene)

  # TODO(dlroxe): The next several lines seem to collect data and write it to
  #               CSV files.  Consider moving them to init_data.py, and adding
  #               tests.
  def get_combined_data(
      gene_name: str) -> Tuple[pandas.DataFrame, pandas.DataFrame]:
    complete, significant = rnaseq_cnv_corr_sum(
      df=get_EIF_CNV_RNAseq_combined_data(
        TCGA_CNV1, TCGA_GTEX_RNAseq, gene_name),
      gene_name=gene_name)

    # The write_csv() in polars is much faster than the to_csv() in pandas.
    polars.from_pandas(significant).write_csv(
      _abspath(FLAGS.output_directory,
               'ProcessedData',
               f'{gene_name}_CNV_COR_RNAseq_sig.csv'))

    return complete, significant

  logging.info('collecting combined RNASeq data')
  # TODO(dlroxe): Make gene names flag-configurable.
  EIF3E_CNV_COR_RNAseq, EIF3E_CNV_COR_RNAseq_sig = get_combined_data("EIF3E")
  EIF3H_CNV_COR_RNAseq, EIF3H_CNV_COR_RNAseq_sig = get_combined_data("EIF3H")
  EIF4A2_CNV_COR_RNAseq, EIF4A2_CNV_COR_RNAseq_sig = get_combined_data("EIF4A2")
  EIF4G1_CNV_COR_RNAseq, EIF4G1_CNV_COR_RNAseq_sig = get_combined_data("EIF4G1")

  # TODO(dlroxe): This plots just the EIF4G1 'significant' dataframe.
  #               If plotting all the other genes similarly is desirable, then
  #               this plotting code should be moved inside get_combined_data().
  #               N.B. In that case, the 4 calls to get_combined_data() above
  #               should be collapsed into a loop.
  h = seaborn.clustermap(
    EIF4G1_CNV_COR_RNAseq_sig,
    # EIF4G1_CNV_COR_RNAseq_sig.drop(['NAT'], axis=1),
    # cbar_pos=(.2, .2, .03, .4),
    cmap="coolwarm",
    figsize=(10, 10),
    # method="centroid",
    # metric="euclidean",
    tree_kws=dict(linewidths=0.5, colors=(0.2, 0.2, 0.4)),
    yticklabels=False
  )
  matplotlib.pyplot.show()

  ### umap analysis on cluster genes from heatmap
  # The CSV files used here are generated by .get_cluster_genes() in Fig2.R.
  # TODO(dlroxe): See about generating the .csv files in init_data.py.
  logging.info('reading cluster CSVs')
  EIF4G1_CNV_RNA_COR_cluster1 = _read_cluster_csv('CNV_RNAseq_COR_cluster1.csv')
  EIF4G1_CNV_RNA_COR_cluster2 = _read_cluster_csv('CNV_RNAseq_COR_cluster2.csv')
  EIF4G1_CNV_RNA_COR_cluster3 = _read_cluster_csv('CNV_RNAseq_COR_cluster3.csv')
  EIF4G1_CNV_RNA_COR_cluster4 = _read_cluster_csv('CNV_RNAseq_COR_cluster4.csv')
  EIF4G1_CNV_RNA_COR_cluster5 = _read_cluster_csv('CNV_RNAseq_COR_cluster5.csv')

  df, df1 = get_RNAse_CNV_cluster(
    TCGA_CNV1, TCGA_GTEX_RNAseq, EIF4G1_CNV_RNA_COR_cluster5, "EIF4G1")

  logging.info('plotting umap')
  plot_umap(df, df1, "cluster5")

  ### UMAP =======================================================

  mapper = umap.UMAP().fit(df1)
  umap.plot.points(mapper)
  umap.plot.points(mapper, labels=df.EIF4G1_CNV)

  ## PCA  ======================================
  penguins_data = df.select_dtypes(numpy.number)
  penguins_data.head()
  penguins_info = df.select_dtypes(exclude='float')
  penguins_info.head()
  cnv = penguins_info.EIF4G1_CNV.tolist()
  pca = sklearn.decomposition.PCA(n_components=4)
  penguins_pca = pca.fit_transform(penguins_data)
  pc_df = pandas.DataFrame(
    data=penguins_pca, columns=['PC1', 'PC2', 'PC3', 'PC4'])
  pc_df.head()

  pc_df['CNV'] = cnv
  pc_df_amp = pc_df.loc[pc_df['CNV'] == "AMP"]
  pc_df_dup = pc_df.loc[pc_df['CNV'] == "DUP"]
  pc_df_dip = pc_df.loc[pc_df['CNV'] == "DIPLOID"]
  pc_df_del = pc_df.loc[pc_df['CNV'] == "DEL"]
  pc_df_nat = pc_df.loc[pc_df['CNV'] == "NAT"]

  logging.info('pc_df.head():\n%s', pc_df.head())
  logging.info(
    'pca.explained_variance_ratio_: %s', pca.explained_variance_ratio_)
  seaborn.set_theme(style='white')

  logging.info('making PCA scatterplot')
  matplotlib.pyplot.figure(figsize=(12, 10))
  with seaborn.plotting_context("notebook", font_scale=1.25):
    seaborn.scatterplot(x="PC1", y="PC2",
                        data=pc_df,
                        hue="CNV",
                        palette=_SCATTERPLOT_COLOR_DICT,
                        alpha=1 / 5,
                        s=10)
  # control x and y limits
  matplotlib.pyplot.xlim(-100, 80)
  matplotlib.pyplot.ylim(-100, 100)
  matplotlib.pyplot.show()
  logging.info('Fig2 work complete')


if __name__ == "__main__":
  app.run(main)

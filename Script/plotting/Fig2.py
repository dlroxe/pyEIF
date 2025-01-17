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
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
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

# TODO(dlroxe): Instead of defining generate_clustermap_plots as a plain
#               boolean flag, perhaps define a stringlist flag that specifies
#               which genes should yield clustermap plots.
flags.DEFINE_boolean('generate_clustermap_plots', False,
                     'When true, causes PCA plots to be generated alongside '
                     'UMAP plots.')
flags.DEFINE_integer('umap_random_state', None,
                     'Force the use of a specific random seed for UMAP plots. '
                     'This forces the creation of UMAP plots that are '
                     'perfectly reproducible for a given state value, '
                     'but risks the creation of a UMAP plot that shows good '
                     'separation as a consequence of the state value rather '
                     'than as a consequence of the data being analyzed.')

_SCATTERPLOT_COLOR_DICT = {
  'AMP': 'red',
  'DEL': 'blue',
  'DIPLOID': 'grey',
  'DUP': 'orange',
  'NAT': 'green',
}


# tcga_cnv_merged_del_thresholds: a dataframe with CNV thresholds,
# in which HOMDEL has been replaced by DEL, because there are not
# enough HOMDEL values for useful analysis, but HOMDEL values ought
# to be counted in a DEL+HOMDEL analysis.
#
# tcga_gtex_rnaseq: contents of TcgaTargetGtex_RSEM_Hugo_norm_count
# def get_eif_cnv_rnaseq_combined_data
def get_eif_cnv_rnaseq_combined_data(
    tcga_cnv_merged_del_thresholds: pandas.DataFrame,
    tcga_gtex_rnaseq: pandas.DataFrame,
    gene_name: str):
  # merge RNAseq and CNV data based on sample names
  # TODO(dlroxe): in get_rnaseq_cnv_cluster() below,
  #               tcga_cnv_merged_del_thresholds[gene_name] was renamed to
  #               col_name... was that intended here, too?
  #               (Perhaps that function and this one should be combined.)
  df = pandas.merge(
    tcga_gtex_rnaseq,
    tcga_cnv_merged_del_thresholds[gene_name],
    how='left',
    left_index=True,
    right_index=True,
    suffixes=("", "_CNV"),
  )
  # remove all row when gene value is zero
  df = df[df[gene_name] != 0]
  # annotate the CNV status from GTEX data as "NAT"
  col_name = f'{gene_name}_CNV'
  df.loc[df.index.str.contains('GTEX'), col_name] = "NAT"
  # remove samples without CNV status, such as target samples
  df = df[df[col_name].notna()]
  return df


def plot_rnaseq_cnv_status(df, gene_name):
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
  for cnv in pandas.Series(["AMP", "DUP", "DIPLOID", "DEL", 'NAT']):
    logging.info('rnaseq_cnv_corr_sum() considering %s:%s', gene_name, cnv)
    # select rows based on the CNV status
    df1 = df.loc[df[col_name] == cnv]
    # remove the CNV status column to facilitate the correlation analysis
    df1 = df1.drop(col_name, axis=1)
    x = df1.corrwith(df1[gene_name])
    x.name = cnv
    y = pandas.concat([y, x], axis=1)
  y.dropna(inplace=True)
  y.index.names = ["Gene_Symbol"]
  # pandas.DataFrame.ge: get Greater than or equal to of dataframe and other
  # pandas.DataFrame.le: get Less than or equal to of dataframe and other
  z = y[y.ge(0.4).any(axis=1) | y.le(-0.4).any(axis=1)]
  return y, z


# tcga_cnv_merged_del_thresholds: a dataframe with CNV thresholds,
# in which HOMDEL has been replaced by DEL, because there are not
# enough HOMDEL values for useful analysis, but HOMDEL values ought
# to be counted in a DEL+HOMDEL analysis.
#
# tcga_gtex_rnaseq: contents of TcgaTargetGtex_RSEM_Hugo_norm_count
# def get_eif_cnv_rnaseq_combined_data
def get_rnaseq_cnv_cluster(
    tcga_cnv_merged_del_thresholds: pandas.DataFrame,
    tcga_gtex_rnaseq: pandas.DataFrame,
    df,
    gene_name):
  cluster_rnaseq = tcga_gtex_rnaseq[df["gene"]]
  # gene_name = "EIF4G1"
  tcga_cnv_gene = tcga_cnv_merged_del_thresholds[gene_name]
  col_name = f'{gene_name}_CNV'
  tcga_cnv_gene.rename(col_name, inplace=True)

  # merge RNAseq and CNV data based on sample names
  df = pandas.merge(
    cluster_rnaseq,
    tcga_cnv_gene,
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
  #
  # The behavior of this function is FLAG-configurable.  When the script is
  # executed with --umap_random_state=<int value>, then <int value> is supplied
  # to the random_state argument below, which will force a reproducible (but
  # possibly misleading) result.  By default (i.e. when unspecified), 'None'
  # is passed, which yields a non-reproducible result that is not at
  # statistical risk.
  embedding = umap.UMAP(
    n_neighbors=10, min_dist=0.3, random_state=FLAGS.umap_random_state
  ).fit_transform(scaled_data)
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
  filename_suffix = 'UMAP'
  if FLAGS.umap_random_state:
    filename_suffix = f'UMAP-fixed-{FLAGS.umap_random_state}'
  matplotlib.pyplot.savefig(
    _abspath(
      FLAGS.output_directory, "Fig2", f'{cluster}-{filename_suffix}.pdf'),
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
  tcga_gtex_rnaseq = datatable.fread(
    _abspath(FLAGS.data_directory, "TcgaTargetGtex_RSEM_Hugo_norm_count"),
    header=True).to_pandas()
  tcga_gtex_rnaseq = tcga_gtex_rnaseq.set_index(['sample']).T

  tcga_cnv_path = _abspath(
    FLAGS.output_directory, 'ProcessedData', 'TCGA-CNV-thresholds.csv')
  logging.info(f'reading: {tcga_cnv_path}')
  tcga_cnv = datatable.fread(tcga_cnv_path, header=True).to_pandas()

  tcga_cnv.set_index([tcga_cnv.columns[0]], inplace=True)
  # too few homdel samples for correlation analysis;
  # combine homdel and del samples for correlation
  tcga_cnv1 = tcga_cnv.replace(['HOMDEL'], 'DEL')

  # TCGA_GTEX_RNAseq_CNV_EIF4G1 = get_EIF_CNV_RNAseq_combined_data(gene_name = "EIF4G1")
  eif4g1_cnv_rna = get_eif_cnv_rnaseq_combined_data(
    tcga_cnv1, tcga_gtex_rnaseq, "EIF4G1")

  # TODO(dlroxe): Make the list of genes flag-configurable.
  logging.info('commencing CNV-status violin plots')
  for gene in (
      "BRCA1", "RANGAP1", "ORC1", "CDC20", "SREBF2", "HMGCR", "HMGCS1",
      "EIF4G1", "EIF4E", "EIF4A2", "EIF3E", "CENPI", "LAMA3", "ITGA6"):
    plot_rnaseq_cnv_status(df=eif4g1_cnv_rna, gene_name=gene)

  # TODO(dlroxe): The next several lines seem to collect data and write it to
  #               CSV files.  Consider moving them to init_data.py, and adding
  #               tests.
  def get_write_and_plot_combined_data(
      local_gene_name: str) -> Tuple[pandas.DataFrame, pandas.DataFrame]:
    logging.info('fetching and plotting cluster data for %s', local_gene_name)
    complete, significant = rnaseq_cnv_corr_sum(
      df=get_eif_cnv_rnaseq_combined_data(
        tcga_cnv1, tcga_gtex_rnaseq, local_gene_name),
      gene_name=local_gene_name)

    # The write_csv() in polars is much faster than the to_csv() in pandas.
    polars.from_pandas(significant).write_csv(
      _abspath(FLAGS.output_directory,
               'ProcessedData',
               f'{local_gene_name}_CNV_COR_RNAseq_sig.csv'))

    if FLAGS.generate_clustermap_plots:
      h = seaborn.clustermap(
        significant,
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

    return complete, significant

  logging.info('collecting combined RNASeq data')
  # TODO(dlroxe): Make gene names flag-configurable.
  for gene_name in ('EIF3E', 'EIF3H', 'EIF4A2', 'EIF4G1'):
    get_write_and_plot_combined_data(gene_name)

  ### umap analysis on cluster genes from heatmap
  # The CSV files used here are generated by .get_cluster_genes() in Fig2.R.
  # TODO(dlroxe): See about generating the .csv files in init_data.py.
  # TODO(dlroxe): Just generate the UMAP for all .csv files below, in a loop?
  logging.info('reading cluster CSVs')
  eif4g1_cnv_rna_cor_cluster1 = _read_cluster_csv('CNV_RNAseq_COR_cluster1.csv')
  eif4g1_cnv_rna_cor_cluster2 = _read_cluster_csv('CNV_RNAseq_COR_cluster2.csv')
  eif4g1_cnv_rna_cor_cluster3 = _read_cluster_csv('CNV_RNAseq_COR_cluster3.csv')
  eif4g1_cnv_rna_cor_cluster4 = _read_cluster_csv('CNV_RNAseq_COR_cluster4.csv')
  eif4g1_cnv_rna_cor_cluster5 = _read_cluster_csv('CNV_RNAseq_COR_cluster5.csv')

  df, df1 = get_rnaseq_cnv_cluster(
    tcga_cnv1, tcga_gtex_rnaseq, eif4g1_cnv_rna_cor_cluster5, "EIF4G1")

  logging.info('plotting umap')
  plot_umap(df, df1, "cluster5")

  # UMAP =======================================================

  mapper = umap.UMAP().fit(df1)
  umap.plot.points(mapper)
  umap.plot.points(mapper, labels=df.EIF4G1_CNV)

  # PCA  ======================================
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
  logging.info('pc_df.head():\n%s', pc_df.head())
  for cnv_var in ('AMP', 'DUP', 'DIPLOID', 'DEL', 'NAT'):
    pc_df_var = pc_df.loc[pc_df['CNV'] == cnv_var]
    logging.info('pc_df_%s.head():\n%s', cnv_var, pc_df_var.head())

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

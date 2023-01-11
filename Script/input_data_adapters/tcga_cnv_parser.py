"""
This file provides abstractions for TCGA CNV data.

The file presumes that certain TCGA data have been downloaded and are locally
available (in unzipped form):

https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz
https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz
https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz

Locations of these files must be specified when the TcgaCnvParser class is
instantiated.

"""

from absl import logging
from scipy import stats
from typing import List, Optional

import datatable
import os
import pandas
import sys

sys.path += ['input_data_adapters']

import org_hs_eg_db_lookup


class TcgaCnvParser:
  """
  Methods and configuration settings for parsing TCGA CNV data.
  """

  # This variable is left public because it encodes the implicit meaning for the
  # values in Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.
  cnv_code_mappings = {
    2.0: 'AMP',
    1.0: 'DUP',
    0.0: 'DIPLOID',
    -1.0: 'DEL',
    -2.0: 'HOMDEL',
  }

  def __init__(
      self,
      data_directory: str,
      output_directory: str,
      cnv_data_by_gene_thresholds: str,
      cnv_data_by_gene_values: str,
      cnv_data_phenotypes: str,
  ):
    self._data_directory = data_directory
    self._output_directory = output_directory
    self._cnv_data_by_gene_thresholds = cnv_data_by_gene_thresholds
    self._cnv_data_by_gene_values = cnv_data_by_gene_values
    self._cnv_data_phenotypes = cnv_data_phenotypes

  def get_tcga_cnv_value(self, raw_data_file: str = None) -> pandas.DataFrame:
    """
    Reads raw_data_file and returns a related dataframe.

    The input file contains raw data in the following form:

                                 TCGA-A5-A0GI-01  TCGA-S9-A7J2-01  TCGA-06-0150-01  ...   TCGA-DD-A115-01
    Sample                                                                          ...
    ACAP3                                  0.0             -1.0              0.0    ...             0.0
    ACTRT2                                 0.0             -1.0              0.0    ...             0.0
    AGRN                                   0.0             -1.0              0.0    ...             0.0
    ANKRD65                                0.0             -1.0              0.0    ...             0.0
    ATAD3A                                 0.0             -1.0              0.0    ...             0.0

    The rows are genes, and the columns are samples.  This function transposes
    the data and selects certain genes.  For example, for certain EIF genes, it
    returns a dataframe of this form:

    Sample          EIF4G1 EIF3E EIF3H
    TCGA-A5-A0GI-01    0.0   0.0   0.0
    TCGA-S9-A7J2-01    0.0   0.0   0.0
    TCGA-06-0150-01    0.0   0.0   0.0
    ...                ...   ...   ...
    TCGA-DD-A115-01    0.0  -1.0  -1.0

    :param raw_data_file: the name of a file (relative to the configured data
      directory) containing raw data
    :return: a data frame with samples as rows.
    """
    input_file = os.path.join(self._data_directory, raw_data_file)
    logging.info('reading from %s', input_file)

    # Unfortunately, pandas documentation suggests that chaining is prone
    # to failure:
    # http://pandas.pydata.org/pandas-docs/dev/user_guide/indexing.html#returning-a-view-versus-a-copy
    #
    # So, 'df' is referenced repeatedly.  OTOH, this approach is probably
    # more memory efficient.
    df = datatable.fread(file=input_file).to_pandas()
    df.sort_values(by=['Sample'], inplace=True)
    df.set_index('Sample', inplace=True)
    return df.transpose()

  def get_tcga_cnv(
      self,
      values_data_frame: Optional[pandas.DataFrame] = None) -> pandas.DataFrame:
    """
    Returns get_tcga_cnv_value(), with numeric cell values replaced by labels.

    Sample output for a selection of EIF genes:

    Sample            EIF4G1    EIF3E    EIF3H
    TCGA-A5-A0GI-01  DIPLOID  DIPLOID  DIPLOID
    TCGA-S9-A7J2-01  DIPLOID  DIPLOID  DIPLOID
    TCGA-06-0150-01  DIPLOID  DIPLOID  DIPLOID
    ...                  ...      ...      ...
    TCGA-DD-A115-01  DIPLOID      DEL      DEL

    :param values_data_frame: if None, then the function uses the value
           returned by
           get_tcga_cnv_value(
           'Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes')
    :return: a data frame with samples as rows, selected genes as columns, and
     string labels as cell values.
    """
    # Note the .replace() call, which just applies the dict, and is very quick.
    if values_data_frame is None:
      values_data_frame = self.get_tcga_cnv_value(
        raw_data_file=self._cnv_data_by_gene_thresholds)
    values_data_frame.replace(self.cnv_code_mappings, inplace=True)
    return values_data_frame

  # TODO(dlroxe): Probably it's worth documenting the join() semantics more
  #               carefully, particularly regarding the indices, in
  #               merge_cnv_phenotypes().
  def merge_cnv_phenotypes(
      self,
      cnv_data: Optional[pandas.DataFrame] = None,
      phenotype_data: Optional[pandas.DataFrame] = None) -> pandas.DataFrame:
    """
    Merges TCGA 'sample type' and 'primary disease' phenotypes with CNV data.

    For example, CNV data might include this row:

    Sample             gene1    gene2    gene3  ...   gene4    gene5    gene6
    TCGA-A5-A0GI-01  DIPLOID  DIPLOID  DIPLOID  ... DIPLOID  DIPLOID  DIPLOID

    Phenotype data might include this row:

                       sample.type                        primary_disease
    Sample
    TCGA-A5-A0GI-01  Primary Tumor  uterine corpus endometrioid carcinoma

    The merged data would look like this:

    Sample             gene1    gene2    gene3  ...   gene6    sample.type                        primary_disease
    TCGA-A5-A0GI-01  DIPLOID  DIPLOID  DIPLOID  ... DIPLOID  Primary Tumor  uterine corpus endometrioid carcinoma


    :param cnv_data: a dataframe obtained from
      get_tcga_cnv() or get_tcga_value()
    :param phenotype_data: a dataframe based derived from data referenced by
      FLAGS.cnv_data_phenotypes
    :return: a merged dataframe that combines CNV value/threshold data with CNV
     phenotype data.
    """
    cnv = self.get_tcga_cnv() if cnv_data is None else cnv_data

    if phenotype_data is None:
      phenotype_data = datatable.fread(
        file=os.path.join(self._data_directory, self._cnv_data_phenotypes)
      ).to_pandas()[['sample', 'sample_type', '_primary_disease']]
      phenotype_data.rename(
        columns={
          'sample': 'Sample',
          'sample_type': 'sample.type',
          '_primary_disease': 'primary_disease',
        }, inplace=True)
      phenotype_data.sort_values(by=['Sample'], inplace=True)
      phenotype_data.set_index('Sample', inplace=True)

    return cnv.join(phenotype_data, how='inner')

  @classmethod
  def get_top_genes(
      cls,
      df: pandas.DataFrame,
      labels: List[str],
      percent: int,
      genedb_handle: org_hs_eg_db_lookup.OrgHsEgDbLookup,
  ) -> pandas.DataFrame:
    sample_count = len(df.index)  # save number of samples before alterations

    # make a copy; then modify the copy in-place
    df = pandas.DataFrame(df, copy=True)
    df.index.name = 'rowname'
    df.reset_index(inplace=True)

    # melt() makes another copy, which is then modified in-place
    df = df.melt(id_vars=['rowname'], var_name='Gene', value_name='Value',
                 ignore_index=True)
    df.set_index('rowname', inplace=True)
    df = pandas.DataFrame(df.loc[lambda x: x['Value'].isin(labels)])
    df = df.groupby(by='Gene').count()
    df = df.apply(lambda x: 100 * x / sample_count)
    df = pandas.DataFrame(df.loc[lambda x: x['Value'] > percent])
    df.reset_index(inplace=True)
    df = df.assign(
      entrez=lambda x: x['Gene'].apply(
        genedb_handle.translate_gene_symbol_to_entrez_id))
    return df

  # TODO(dlroxe): Fix up the function docstring below.
  def cooccurance_analysis(self, df: pandas.DataFrame, gene01: str, gene02: str,
                           cnv_spec: List[str]) -> None:
    """
    For example, 'sheet 1' should have something like this:

                          EIF3H_AMP_DUP    EIF3H_NO_AMP_DUP
        EIF4G1_AMP_DUP    2220             1233
        EIF4G1_NO_AMP_DUP 2528             4864

    That is: 2220 samples are AMP|DUP for G1, AND are either AMP|DUP for 3H.
    """
    df = df[[gene01, gene02]]

    # TODO(dlroxe): Consider using 'crosstab' for this.
    gene01y_gene02y = len(
      df[df[gene01].isin(cnv_spec) & df[gene02].isin(cnv_spec)])
    gene01y_gene02n = len(
      df[df[gene01].isin(cnv_spec) & ~df[gene02].isin(cnv_spec)])
    gene01n_gene02y = len(
      df[~df[gene01].isin(cnv_spec) & df[gene02].isin(cnv_spec)])
    gene01n_gene02n = len(
      df[~df[gene01].isin(cnv_spec) & ~df[gene02].isin(cnv_spec)])

    def row_or_col_name(gene: str, matches_cnv_spec: bool) -> str:
      """Returns e.g. EIF4G1_AMP_DUP, EIF4G1_NO_AMP_DUP, EIF4G1_AMP, etc."""
      components = [
        x for x in
        ([gene, None if matches_cnv_spec else 'NO'] + cnv_spec) if x]
      return ' '.join(components)

    eif = pandas.DataFrame(
      data={
        row_or_col_name(gene=gene02, matches_cnv_spec=True):
          [gene01y_gene02y, gene01n_gene02y],
        row_or_col_name(gene=gene02, matches_cnv_spec=False):
          [gene01y_gene02n, gene01n_gene02n],
      },
      index=[
        row_or_col_name(gene=gene01, matches_cnv_spec=True),
        row_or_col_name(gene=gene01, matches_cnv_spec=False),
      ])

    logging.info('got adjusted counts:\n%s', eif)

    odds_ratio, p_value = stats.fisher_exact(eif, alternative='greater')
    fisher = pandas.DataFrame(
      data={'Odds Ratio': [odds_ratio], 'P Value': [p_value]})
    logging.info('got fisher test:\n%s', fisher)

    chi_sq, p_value = stats.chisquare(eif)
    chi_test = pandas.DataFrame(
      data={'Chi-Squared': [chi_sq], 'P Value': [p_value]})
    logging.info('got chi-sq test:\n%s', chi_test)

    excel_output_file = os.path.join(
      self._output_directory, "Fig1",
      '_'.join([gene01, gene02] + cnv_spec) + '.xlsx')
    with pandas.ExcelWriter(path=excel_output_file) as writer:
      # TODO(dlroxe): rownames=True, then False, False
      eif.to_excel(writer, sheet_name='1')
      fisher.to_excel(writer, sheet_name='Fisheroneside')
      chi_test.to_excel(writer, sheet_name='chi_test')

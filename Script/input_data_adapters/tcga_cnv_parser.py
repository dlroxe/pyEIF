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
import org_hs_eg_db_lookup
import pandas


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
      data_directory: Optional[str],
      output_directory: Optional[str],
      cnv_data_by_gene_thresholds: Optional[str],
      cnv_data_by_gene_values: Optional[str],
      cnv_data_phenotypes: Optional[str],
  ):
    # member variables copied from constructor args
    self._data_directory: Optional[str] = data_directory
    self._output_directory: Optional[str] = output_directory
    self._cnv_data_by_gene_thresholds: Optional[
      str] = cnv_data_by_gene_thresholds
    self._cnv_data_by_gene_values: Optional[str] = cnv_data_by_gene_values
    self._cnv_data_phenotypes: Optional[str] = cnv_data_phenotypes

    # computed values
    self._raw_values_data: Optional[
      pandas.DataFrame] = self._init_raw_values_data()
    self._threshold_raw_data: Optional[
      pandas.DataFrame] = self._init_threshold_raw_data()  # -2, -1, etc.
    self._threshold_data: Optional[
      pandas.DataFrame] = self._init_threshold_data()  # HOMDEL, DEL, etc.
    self._melted_threshold_data: Optional[
      pandas.DataFrame] = self._init_melted_threshold_data()
    self._phenotype_data: Optional[
      pandas.DataFrame] = self._init_phenotype_data()

  # As usual for Python, the leading '_' for this and the other _init...()
  # functions that follow is a signal to programmers that they ought not to
  # call this function from outside the class.  In this case, the reason is
  # that reading the CNV data from disk is computationally expensive, and the
  # idea is to arrange the program in such a way as to guarantee that the
  # expensive operation is performed exactly once (specifically, during
  # the __init__() call when an instance of TcgaCnvParser is constructed.
  def _init_raw_values_data(self) -> Optional[pandas.DataFrame]:
    if self._cnv_data_by_gene_values is None:
      logging.warning('no raw CNV data file specified')
      return None

    raw_data_file = self._abspath(
      os.path.join(self._data_directory, self._cnv_data_by_gene_values))

    if not os.path.exists(raw_data_file):
      logging.warning('raw CNV data file does not exist: %s', raw_data_file)
      return None

    logging.info('reading from %s', raw_data_file)

    # Unfortunately, pandas documentation suggests that chaining is prone
    # to failure:
    #
    # http://pandas.pydata.org/pandas-docs/dev/user_guide/indexing.html#returning-a-view-versus-a-copy
    #
    # So, 'df' is referenced repeatedly.  OTOH, this approach is probably
    # more memory efficient.
    #
    # Values are (well) within the range [-5.0, 5.0], so it is tempting to
    # save memory by converting from 'float64' to 'float32'.  However, they are
    # stored with 3 decimal places, and use of 'float32' results in loss of
    # precision (e.g. '0.010' -> '0.010002').  So, the default 'float64' width
    # is retained.
    df = datatable.fread(file=raw_data_file).to_pandas()
    df.sort_values(by=['Sample'], inplace=True)
    df.set_index('Sample', inplace=True)
    return df.transpose()

  def _init_threshold_raw_data(self) -> Optional[pandas.DataFrame]:
    if self._cnv_data_by_gene_thresholds is None:
      logging.warning('no CNV threshold data file specified')
      return None

    raw_data_file = self._abspath(
      os.path.join(self._data_directory, self._cnv_data_by_gene_thresholds))

    return self._read_tcga_cnv_values(raw_data_file)

  def _init_threshold_data(self) -> Optional[pandas.DataFrame]:
    # Note the .replace() call, which just applies the dict, and is very quick.
    #
    # Because there are only a handful of allowed cell values, setting the
    # 'category' data type yields a substantial memory savings, from 678,819
    # bytes per column for the full data set, to 11,324 bytes per column (that's
    # a factor of 60).
    if self._threshold_raw_data is None:
      logging.warning('no raw threshold data is available')
      return None
    values_data_frame = pandas.DataFrame(self._threshold_raw_data, copy=True)
    logging.info('initializing thresholds from:\n%s', self._threshold_raw_data)
    values_data_frame.replace(self.cnv_code_mappings, inplace=True)
    return values_data_frame.astype('category')

  def _init_melted_threshold_data(self):
    if self._threshold_data is None:
      logging.warning('no threshold data is available')
      return None

    df = pandas.DataFrame(self._threshold_data, copy=True)
    df.index.name = 'rowname'
    df.reset_index(inplace=True)
    df = df.melt(
      id_vars=['rowname'], var_name='Gene', value_name='Value',
      ignore_index=True)
    df.set_index('rowname', inplace=True)
    return df

  def _init_phenotype_data(self) -> Optional[pandas.DataFrame]:
    if not self._cnv_data_phenotypes:
      logging.warning('no phenotype data is available')
      return None

    phenotype_file = self._abspath(os.path.join(
      self._data_directory, self._cnv_data_phenotypes))

    phenotype_data = (
      self._read_csv(
        filename=phenotype_file)[['sample', 'sample_type', '_primary_disease']]
      .astype('string')  # accurate, and mem-efficient vs. default 'object'
    )
    phenotype_data.rename(
      columns={
        'sample': 'Sample',
        'sample_type': 'sample.type',
        '_primary_disease': 'primary_disease',
      }, inplace=True)
    phenotype_data.sort_values(by=['Sample'], inplace=True)
    phenotype_data.set_index('Sample', inplace=True)
    logging.info('initialized phenotype data as:\n%s', phenotype_data)
    return phenotype_data

  # TODO(dlroxe):  This function is cropping up in a few places; try to find
  #                a One True Home for it.
  @staticmethod
  def _abspath(path):
    # What an absurd incantation to resolve "~"; but OK. Thanks, StackOverflow.
    return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))

  @classmethod
  def _read_tcga_cnv_values(cls, raw_data_file: str) -> pandas.DataFrame:
    """
    Reads raw_data_file and returns a transposed and indexed dataframe.

    :param raw_data_file: the name of a file (relative to the configured data
      directory) containing raw data
    :return: a data frame with samples as rows.
    """

    # Unfortunately, pandas documentation suggests that chaining is prone
    # to failure:
    #
    # http://pandas.pydata.org/pandas-docs/dev/user_guide/indexing.html#returning-a-view-versus-a-copy
    #
    # So, 'df' is referenced repeatedly.  OTOH, this approach is probably
    # more memory efficient.
    #
    # Values are (well) within the range [-5.0, 5.0], so it is tempting to
    # save memory by converting from 'float64' to 'float32'.  However, they are
    # stored with 3 decimal places, and use of 'float32' results in loss of
    # precision (e.g. '0.010' -> '0.010002').  So, the default 'float64' width
    # is retained.  (One supposes that all the values could be multiplied by
    # 1,000 and retained in 16-bit integers, but conversions seem like more
    # trouble than they're worth until and unless memory usage becomes a far
    # more critical concern.)
    df = cls._read_csv(raw_data_file)
    df.sort_values(by=['Sample'], inplace=True)
    df.set_index('Sample', inplace=True)
    return df.transpose()

  @staticmethod
  def _read_csv(filename: str) -> Optional[pandas.DataFrame]:
    if os.path.exists(filename):
      logging.info('reading data from: %s', filename)
      return datatable.fread(file=filename).to_pandas()

    logging.error('file does not exist: %s', filename)
    return None

  def get_tcga_cnv_value(self) -> pandas.DataFrame:
    """
    Returns TCGA CNV data organized with genes as columns and samples as rows.

    For example, TCGA provides threshold data in this format:

              TCGA-A5-A0GI-01  TCGA-S9-A7J2-01  TCGA-06-0150-01  ...   TCGA-DD-A115-01
    Sample                                                       ...
    ACAP3               0.0             -1.0              0.0    ...             0.0
    ACTRT2              0.0             -1.0              0.0    ...             0.0
    AGRN                0.0             -1.0              0.0    ...             0.0
    ANKRD65             0.0             -1.0              0.0    ...             0.0
    ATAD3A              0.0             -1.0              0.0    ...             0.0

    The rows are genes, and the columns are samples.  This function returns a
    transposed copy of the data:

    Sample            ACAP3  ACTRT2   AGRN   ANKRD65  ATAD3A
    TCGA-A5-A0GI-01    0.0     0.0    0.0       0.0     0.0
    TCGA-S9-A7J2-01   -1.0    -1.0   -1.0      -1.0    -1.0
    TCGA-06-0150-01    0.0     0.0    0.0       0.0     0.0
    ...                ...     ...    ...       ...     ...
    TCGA-DD-A115-01    0.0     0.0    0.0       0.0     0.0

    :return: a data frame with samples as rows.  The returned data is always
             a freshly-copied dataframe, which can be manipulated without fear
             of altering the underlying data.
    """
    return pandas.DataFrame(self._raw_values_data, copy=True)

  def get_tcga_cnv_threshold_categories(self) -> pandas.DataFrame:
    """
    Returns CNV threshold data , with numeric cell values replaced by labels.

    Sample output for a selection of EIF genes:

    Sample            EIF4G1    EIF3E    EIF3H
    TCGA-A5-A0GI-01  DIPLOID  DIPLOID  DIPLOID
    TCGA-S9-A7J2-01  DIPLOID  DIPLOID  DIPLOID
    TCGA-06-0150-01  DIPLOID  DIPLOID  DIPLOID
    ...                  ...      ...      ...
    TCGA-DD-A115-01  DIPLOID      DEL      DEL

    :return: a data frame with samples as rows, selected genes as columns, and
     string labels as cell values.  Every call to this function returns a fresh
     copy of the data, which may be manipulated without concern for the
     integrity of the underlying data.
    """
    return pandas.DataFrame(self._threshold_data, copy=True)

  # TODO(dlroxe): Probably it's worth documenting the join() semantics more
  #               carefully, particularly regarding the indices, in
  #               merge_cnv_phenotypes().
  def merge_cnv_thresholds_and_phenotypes(self) -> pandas.DataFrame:
    """
    Merges TCGA 'sample type' and 'primary disease' phenotypes with CNV
    threshold data.

    For example, CNV threshold data might include this row:

    Sample             gene1    gene2    gene3  ...   gene4    gene5    gene6
    TCGA-A5-A0GI-01  DIPLOID  DIPLOID  DIPLOID  ... DIPLOID  DIPLOID  DIPLOID

    Phenotype data might include this row:

                       sample.type                        primary_disease
    Sample
    TCGA-A5-A0GI-01  Primary Tumor  uterine corpus endometrioid carcinoma

    The merged data would look like this:

    Sample             gene1    gene2    gene3  ...   gene6    sample.type                        primary_disease
    TCGA-A5-A0GI-01  DIPLOID  DIPLOID  DIPLOID  ... DIPLOID  Primary Tumor  uterine corpus endometrioid carcinoma

    :return: a merged dataframe that combines CNV value/threshold data with CNV
     phenotype data.
    """
    cnv = self.get_tcga_cnv_threshold_categories()  # this is a fresh copy
    return cnv.join(self._phenotype_data, how='inner')

  # TODO(dlroxe): Document this method.
  def get_top_genes(
      self,
      labels: List[str],
      percent: int,
      genedb_handle: org_hs_eg_db_lookup.OrgHsEgDbLookup,
  ) -> pandas.DataFrame:
    sample_count = len(self._threshold_data.index)
    df = self._melted_threshold_data
    df = df.loc[lambda x: x['Value'].isin(labels)]
    df = df.groupby(by='Gene').count()
    df = df.apply(lambda x: 100 * x / sample_count)
    df = df.loc[lambda x: x['Value'] > percent]
    df.reset_index(inplace=True)

    # Synthesize a new column named 'entrez', which contains the Entrez ID
    # corresponding to 'Gene', or pandas.NA if no Entrez ID can be determined.
    #
    # The new column is specified with a nullable int64 type (that is, an int64
    # that permits <NA> values).  The capital-I 'Int64' is the documented Pandas
    # type alias for this purpose.  Values in this column will either be 64-bit
    # integers or 'pandas.NA', which appears in printed form as '<NA>'.
    #
    # See also "Nullable integer data type" in the Pandas docs:
    #
    # https://pandas.pydata.org/docs/user_guide/integer_na.html
    #
    # See also this article regarding Entrez IDs and 64-bit ints:
    #
    # https://ncbiinsights.ncbi.nlm.nih.gov/2021/09/02/64-bit-gis/
    df = df.assign(
      entrez=lambda x: x['Gene'].apply(
        genedb_handle.translate_gene_symbol_to_entrez_id))
    return df.astype({'entrez': 'Int64'})

  # TODO(dlroxe): Fix up the function docstring below.
  def co_occurrence_analysis(self, df: pandas.DataFrame, gene01: str,
                             gene02: str,
                             cnv_spec: List[str]) -> None:
    """
    For example, 'sheet 1' should have something like this:

                          EIF3H_AMP_DUP    EIF3H_NO_AMP_DUP
        EIF4G1_AMP_DUP    2220             1233
        EIF4G1_NO_AMP_DUP 2528             4864

    That is: 2220 samples are AMP|DUP for G1, AND are either AMP|DUP for 3H.
    """
    df = df[[gene01, gene02]]

    # TODO(dlroxe): Consider using 'crosstab' for this (see also
    #               'contingency table'.
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

#!/usr/bin/python3
"""
Read and parse raw data related to the pyEIF project.

This file may be imported as a module or run as a stand-alone executable.

The program requires data, which should be separately downloaded into a
data directory.  The data location may be specified by a command-line flag.

./init_data.py --data_directory=~/Desktop/pyeif_data

For a description of available flags, execute with the --help option:

./init_data.py --help

This will show additional options, such as a way to specify names of
individual data files.  Typically, it should not be necessary to specify
them, because by default the program references files using the same
names they are given in the repositories where they are officially
maintained.

"""
from absl import app
from absl import flags
from typing import List, Optional
import datatable
import os
import pandas
from scipy import stats

FLAGS = flags.FLAGS
flags.DEFINE_string('data_directory', os.path.join('~', 'Desktop', 'pyeif_data'), 'parent dir for data files')
flags.DEFINE_string('output_directory', os.path.join('~', 'Desktop', 'pyeif_output'), 'parent dir for output')
flags.DEFINE_string('cnv_data_by_gene_values', 'Gistic2_CopyNumber_Gistic2_all_data_by_genes',
                    'the path, relative to data_directory, where raw data values can be '
                    'found for tissue samples with gene copy numbers.')
flags.DEFINE_string('cnv_data_by_gene_thresholds', 'Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes',
                    'the path, relative to data_directory, where threshold data values can be '
                    'found for tissue samples with gene copy numbers.  These values are constrained to be integers '
                    'in the range [-2, 2].')
flags.DEFINE_string('cnv_data_phenotypes', 'TCGA_phenotype_denseDataOnlyDownload.tsv',
                    'the path, relative to data_directory, where phenotype data can be '
                    'found for tissue samples named in cnv_data_by_gene_values and '
                    'cnv_data_by_gene_thresholds.')


class TcgaCnvParser:
  """
  Methods and configuration settings for parsing TCGA CNV data.
  """

  # This variable is left public because it encodes the implicit meaning for the values
  # in Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.
  cnv_code_mappings = {
    2.0: 'AMP',
    1.0: 'DUP',
    0.0: 'DIPLOID',
    -1.0: 'DEL',
    -2.0: 'HOMDEL',
  }

  def __init__(self):
    pass

  @classmethod
  def get_tcga_cnv_value(cls, raw_data_file: Optional[str] = None) -> pandas.DataFrame:
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

    The rows are genes, and the columns are samples.  This function transposes the data and selects certain genes.
    For example, for certain EIF genes, it returns a dataframe of this form:

    Sample          EIF4G1 EIF3E EIF3H
    TCGA-A5-A0GI-01    0.0   0.0   0.0
    TCGA-S9-A7J2-01    0.0   0.0   0.0
    TCGA-06-0150-01    0.0   0.0   0.0
    ...                ...   ...   ...
    TCGA-DD-A115-01    0.0  -1.0  -1.0

    :param raw_data_file: the name of a file (relative to the configured data directory) containing raw data
    :return: a data frame with samples as rows.
    """
    return datatable.fread(file=os.path.join(FLAGS.data_directory, raw_data_file)
                           ).to_pandas().sort_values(by=['Sample']).set_index('Sample').transpose()

  @classmethod
  def get_tcga_cnv(cls, values_data_frame: Optional[pandas.DataFrame] = None) -> pandas.DataFrame:
    """
    Returns the output of get_tcga_cnv_value(), but with numeric cell values replaced by labels.

    Sample output for a selection of EIF genes:

    Sample            EIF4G1    EIF3E    EIF3H
    TCGA-A5-A0GI-01  DIPLOID  DIPLOID  DIPLOID
    TCGA-S9-A7J2-01  DIPLOID  DIPLOID  DIPLOID
    TCGA-06-0150-01  DIPLOID  DIPLOID  DIPLOID
    ...                  ...      ...      ...
    TCGA-DD-A115-01  DIPLOID      DEL      DEL

    :param values_data_frame: if None, then the function uses the value
           returned by get_tcga_cnv_value('Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes')
    :return: a data frame with samples as rows, selected genes as columns, and string labels as cell values.
    """
    # Note the .replace() call.  It just applies the dictionary, and is very quick.
    if values_data_frame is None:
      values_data_frame = cls.get_tcga_cnv_value(raw_data_file=FLAGS.cnv_data_by_gene_thresholds)
    return values_data_frame.replace(cls.cnv_code_mappings)

  # TODO(dlroxe): Probably it's worth documenting the join() semantics more carefully, particularly regarding the
  # indices, in merge_cnv_phenotypes().
  @classmethod
  def merge_cnv_phenotypes(cls, cnv_data: Optional[pandas.DataFrame] = None,
                           phenotype_data: Optional[pandas.DataFrame] = None) -> pandas.DataFrame:
    """
    Merges TCGA phenotype data, specifically 'sample type' and 'primary disease', with CNV data.

    For example, CNV data might include this row:

    Sample          7SK|ENSG00000232512.2 7SK|ENSG00000249352.3 7SK|ENSG00000254144.2  ... snoZ6|ENSG00000264452.1 snoZ6|ENSG00000266692.1 snosnR66
    TCGA-A5-A0GI-01               DIPLOID               DIPLOID               DIPLOID  ...                 DIPLOID                 DIPLOID  DIPLOID

    Phenotype data might include this row:

                       sample.type                        primary_disease
    Sample
    TCGA-A5-A0GI-01  Primary Tumor  uterine corpus endometrioid carcinoma

    The merged data would look like this:

                    7SK|ENSG00000232512.2 7SK|ENSG00000249352.3 7SK|ENSG00000254144.2  ... snosnR66    sample.type                        primary_disease
    TCGA-A5-A0GI-01               DIPLOID               DIPLOID               DIPLOID  ...  DIPLOID  Primary Tumor  uterine corpus endometrioid carcinoma


    :param cnv_data: a dataframe obtained from get_tcga_cnv() or get_tcga_value()
    :return: a merged dataframe that combines CNV value/threshold data with CNV phenotype data.
    """
    cnv = TcgaCnvParser.get_tcga_cnv() if cnv_data is None else cnv_data

    if phenotype_data is None:
      phenotype_data = datatable.fread(file=os.path.join(FLAGS.data_directory, FLAGS.cnv_data_phenotypes)
                                       ).to_pandas()[['sample', 'sample_type', '_primary_disease']].rename(
        columns={
          'sample': 'Sample',
          'sample_type': 'sample.type',
          '_primary_disease': 'primary_disease',
        }).sort_values(by=['Sample']).set_index('Sample')

    return cnv.join(phenotype_data, how='inner')

  @classmethod
  def get_top_genes(cls, df: pandas.DataFrame, labels: List[str], percent: int) -> pandas.DataFrame:
    sample_number = len(df.index)
    df.index.name = 'rowname'
    return (  # Extra outer parens here permit easy formatting that starts each chained function call on its own line.
      df.reset_index()
      .melt(id_vars=['rowname'], var_name='Gene', value_name='Value', ignore_index=True)
      .set_index('rowname')
      .loc[lambda x: x['Value'].isin(labels)]
      .groupby(by='Gene')
      .count()
      .apply(lambda x: 100 * x / sample_number)
      .loc[lambda x: x['Value'] > percent]
    )
    # TODO(dlroxe): equivalent of:
    '''
    dplyr::mutate(entrez = AnnotationDbi::mapIds(org.Hs.eg.db,
                                                     keys = as.character(.data$Gene),
                                                     column = "ENTREZID",
                                                     keytype = "SYMBOL",
                                                     multiVals = "first"
        ))
    '''

  @classmethod
  def coocurrance_analysis(cls, df: pandas.DataFrame, gene01: str, gene02: str, cnv_1: str, cnv_2: str) -> None:
    def _rename_gene(gene: str) -> str:
      """Returns a name based on 'gene' and the CNV parameters."""
      modifier = '' if gene in (cnv_1, cnv_2) else 'NO'
      return gene + modifier + cnv_1 + cnv_2

    new_gene01 = _rename_gene(gene01)
    new_gene02 = _rename_gene(gene02)
    eif = (
      df[[gene01, gene02]]  # TODO(dlroxe): mimic "all_of()", i.e. error unless both genes are in df
      .rename(columns={
        gene01: new_gene01,
        gene02: new_gene02,
      })
    )
    print(f'got adjusted eif\n{eif}')

    odds_ratio, p_value = stats.fisher_exact(eif, alternative='greater')
    fisher = pandas.DataFrame(data={'Odds Ratio': [odds_ratio], 'P Value': [p_value]})

    chi_sq, p_value = stats.chisquare(eif[new_gene01], eif[new_gene02])
    chi_test = pandas.DataFrame(data={'Chi-Squared': [chi_sq], 'P Value': [p_value]})

    with pandas.ExcelWriter(
        path=os.path.join(FLAGS.output_directory, "Fig1", gene01 + gene02 + cnv_1 + cnv_2 + '.xlsx')) as writer:
      eif.to_excel(writer, sheet_name='1')  # TODO(dlroxe): rownames=True
      fisher.to_excel(writer, sheet_name='Fisheroneside')  # TODO(dlroxe): rownames=False
      chi_test.to_excel(writer, sheet_name='chi_test')  # TODO(dlroxe): rownames=False


def main(argv):
  eif_genes = ["EIF4G1", "EIF3E", "EIF3H", ]

  all_data = TcgaCnvParser.get_tcga_cnv_value(raw_data_file=FLAGS.cnv_data_by_gene_values)
  print(f'all data\n{all_data}')

  raw_threshold_data = TcgaCnvParser.get_tcga_cnv_value(raw_data_file=FLAGS.cnv_data_by_gene_thresholds)
  all_threshold_data = TcgaCnvParser.get_tcga_cnv(values_data_frame=raw_threshold_data)
  eif_threshold_data = all_threshold_data[eif_genes]
  print(f'eif threshold data\n{eif_threshold_data}')

  merged_phenotyped_data = TcgaCnvParser.merge_cnv_phenotypes(all_threshold_data)
  print(f'all threshold data, merged with phenotypes:\n{merged_phenotyped_data}')

  # TOP_AMP_PATH, TOP_GAIN_PATH, TOP_HOMDEL_PATH are omitted for the time being,
  # because pathway analysis is harder in Python than in R.
  top_amp_genes = TcgaCnvParser.get_top_genes(df=all_threshold_data, labels=["AMP"], percent=5)
  print(f'top amp genes:\n{top_amp_genes}')

  top_gain_genes = TcgaCnvParser.get_top_genes(df=all_threshold_data, labels=["DUP", "AMP"], percent=30)
  print(f'top gain genes:\n{top_gain_genes}')

  top_homdel_genes = TcgaCnvParser.get_top_genes(df=all_threshold_data, labels=["HOMDEL"], percent=5)
  print(f'top homdel genes:\n{top_homdel_genes}')

  top_genes_base_path = os.path.join(FLAGS.output_directory, "Fig1")
  with pandas.ExcelWriter(path=os.path.join(top_genes_base_path, 'TOP_AMP_genes.xlsx')) as writer:
    top_amp_genes.to_excel(writer, sheet_name='1')  # TODO(dlroxe): rownames=True

  with pandas.ExcelWriter(path=os.path.join(top_genes_base_path, 'TOP_GAIN_genes.xlsx')) as writer:
    top_gain_genes.to_excel(writer, sheet_name='1')  # TODO(dlroxe): rownames=True

  with pandas.ExcelWriter(path=os.path.join(top_genes_base_path, 'TOP_HOMDEL_genes.xlsx')) as writer:
    top_homdel_genes.to_excel(writer, sheet_name='1')  # TODO(dlroxe): rownames=True

  print(f'"top genes" analyses are written to {top_genes_base_path}.')

  # TODO(dlroxe): The following function calls all crash.  I'm can't quite discern the intent of the original R code.
  print('attempting coocurrance analysis 1')
  TcgaCnvParser.coocurrance_analysis(df=all_threshold_data, gene01="EIF4G1", gene02="EIF3E", cnv_1="AMP", cnv_2="AMP")

  print('attempting coocurrance analysis 2')
  TcgaCnvParser.coocurrance_analysis(df=all_threshold_data, gene01="EIF4G1", gene02="EIF3E", cnv_1="AMP", cnv_2="DUP")

  print('attempting coocurrance analysis 3')
  TcgaCnvParser.coocurrance_analysis(df=all_threshold_data, gene01="EIF4G1", gene02="EIF3H", cnv_1="AMP", cnv_2="AMP")

  print('attempting coocurrance analysis 4')
  TcgaCnvParser.coocurrance_analysis(df=all_threshold_data, gene01="EIF4G1", gene02="EIF3H", cnv_1="AMP", cnv_2="DUP")
  print('processing complete')


if __name__ == "__main__":
  app.run(main)

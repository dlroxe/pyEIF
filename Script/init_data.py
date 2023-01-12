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
from absl import logging

import datatable
import os
import pandas
import sys

sys.path += ['input_data_adapters']

import org_hs_eg_db_lookup
import tcga_cnv_parser

FLAGS = flags.FLAGS
flags.DEFINE_string('data_directory',
                    os.path.join('~', 'Desktop', 'pyeif_data'),
                    'parent dir for data files')
flags.DEFINE_string('output_directory',
                    os.path.join('~', 'Desktop', 'pyeif_output'),
                    'parent dir for output')
flags.DEFINE_string('cnv_data_by_gene_values',
                    'Gistic2_CopyNumber_Gistic2_all_data_by_genes',
                    'the path, relative to data_directory, where raw data '
                    'values can be found for tissue samples with gene copy '
                    'numbers.')
flags.DEFINE_string('cnv_data_by_gene_thresholds',
                    'Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes',
                    'the path, relative to data_directory, where threshold '
                    'data values can be found for tissue samples with gene '
                    'copy numbers.  These values are constrained to be '
                    'integers in the range [-2, 2].')
flags.DEFINE_string('cnv_data_phenotypes',
                    'TCGA_phenotype_denseDataOnlyDownload.tsv',
                    'the path, relative to data_directory, where phenotype '
                    'data can be found for tissue samples named in '
                    'cnv_data_by_gene_values and cnv_data_by_gene_thresholds.')
flags.DEFINE_string('org_hs_eg_sqlite',
                    'org.Hs.eg.db/inst/extdata/org.Hs.eg.sqlite',
                    'the path, relative to data_directory, where the database '
                    'underlying the "org.Hs.eg.db" package may be found. '
                    'This package may be obtained via the "Source Package" '
                    'link at '
                    'https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html '
                    'and deployed to the data_directory using "tar -xzvf".')

# Note, these flags designate files in *output_directory*.
flags.DEFINE_string('top_amp_path',
                    'ProcessedData/top-amp-path.csv',
                    'the path, relative to output_directory, of the '
                    'csv-formatted data frame that describes gene pathways as '
                    'evaluated by Fig1.R.')
flags.DEFINE_string('top_gain_path',
                    'ProcessedData/top-gain-path.csv',
                    'the path, relative to output_directory, of the '
                    'csv-formatted data frame that describes gene pathways as '
                    'evaluated by Fig1.R.')
flags.DEFINE_string('top_homdel_path',
                    'ProcessedData/top-homdel-path.csv',
                    'the path, relative to output_directory, of the '
                    'csv-formatted data frame that describes gene pathways as '
                    'evaluated by Fig1.R.')

# force exceptions for bad chains of pandas operations
pandas.options.mode.chained_assignment = 'raise'


def _abspath(path):
  # What an absurd incantation to resolve "~"; but OK. Thanks, StackOverflow.
  return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))


def main(unused_argv):
  # TODO(dlroxe): Decide whether logging to stdout is desired in the long term.
  logging.get_absl_handler().python_handler.stream = sys.stdout

  eif_genes = ["EIF4G1", "EIF3E", "EIF3H", ]

  parser = tcga_cnv_parser.TcgaCnvParser(
    data_directory=FLAGS.data_directory,
    output_directory=FLAGS.output_directory,
    cnv_data_by_gene_thresholds=FLAGS.cnv_data_by_gene_thresholds,
    cnv_data_by_gene_values=FLAGS.cnv_data_by_gene_values,
    cnv_data_phenotypes=FLAGS.cnv_data_phenotypes,
  )

  all_data = parser.get_tcga_cnv_value(
    raw_data_file=FLAGS.cnv_data_by_gene_values)
  logging.info('all data\n%s', all_data)

  raw_threshold_data = parser.get_tcga_cnv_value(
    raw_data_file=FLAGS.cnv_data_by_gene_thresholds)
  all_threshold_data = parser.get_tcga_cnv(
    values_data_frame=raw_threshold_data)
  eif_threshold_data = all_threshold_data[eif_genes]
  logging.info('eif threshold data\n%s', eif_threshold_data)

  merged_phenotyped_data = parser.merge_cnv_phenotypes(
    all_threshold_data)
  logging.info('all threshold data, merged with phenotypes:\n%s',
               merged_phenotyped_data)

  # TOP_AMP_PATH, TOP_GAIN_PATH, TOP_HOMDEL_PATH are obtained from R-generated
  # CSV files for the time being, because pathway analysis is harder in Python
  # than in R.
  #
  # The relevant R code invokes ReactomePA::enrichPathway(),
  # which is implemented using 'enricher_internal()', which is defined here:
  #
  # https://github.com/YuLab-SMU/DOSE/blob/37572b5a462843dd2478ecf4bcf583bbedd1a357/R/enricher_internal.R
  #
  # It may be possible to build similar functionality in Python on top of the
  # same SQLite database used by the Reactome package:
  #
  # reactome_db_handle = reactome_lookup.ReactomeLookup(
  #  _abspath(os.path.join(FLAGS.data_directory, FLAGS.reactome_sqlite)))

  org_hs_eg_db_handle = org_hs_eg_db_lookup.OrgHsEgDbLookup(
    _abspath(os.path.join(FLAGS.data_directory, FLAGS.org_hs_eg_sqlite)))

  top_amp_path = datatable.fread(
    file=os.path.join(FLAGS.output_directory, FLAGS.top_amp_path)).to_pandas()

  sample_count = len(all_threshold_data.index)
  melted_threshold_data = parser.melt_threshold_data(all_threshold_data)

  top_amp_genes = parser.get_top_genes(
    sample_count=sample_count,
    df=melted_threshold_data, labels=["AMP"], percent=5,
    genedb_handle=org_hs_eg_db_handle)
  logging.info('top amp genes:\n%s', top_amp_genes)

  top_gain_path = datatable.fread(
    file=os.path.join(FLAGS.output_directory, FLAGS.top_gain_path)).to_pandas()
  top_gain_genes = parser.get_top_genes(
    sample_count=sample_count,
    df=melted_threshold_data,
    labels=["DUP", "AMP"],
    percent=30,
    genedb_handle=org_hs_eg_db_handle)
  logging.info('top gain genes:\n%s', top_gain_genes)

  top_homdel_path = datatable.fread(
    file=os.path.join(FLAGS.output_directory,
                      FLAGS.top_homdel_path)).to_pandas()
  top_homdel_genes = parser.get_top_genes(
    sample_count=sample_count,
    df=melted_threshold_data,
    labels=["HOMDEL"], percent=5,
    genedb_handle=org_hs_eg_db_handle)
  logging.info('top homdel genes:\n%s', top_homdel_genes)

  top_genes_base_path = os.path.join(FLAGS.output_directory, "Fig1")
  with pandas.ExcelWriter(
      path=os.path.join(top_genes_base_path, 'TOP_AMP_genes.xlsx')) as writer:
    # TODO(dlroxe): rownames=True
    top_amp_genes.to_excel(writer, sheet_name='1')
    top_amp_path.to_excel(writer, sheet_name='2')

  with pandas.ExcelWriter(
      path=os.path.join(top_genes_base_path, 'TOP_GAIN_genes.xlsx')) as writer:
    # TODO(dlroxe): rownames=True
    top_gain_genes.to_excel(writer, sheet_name='1')
    top_gain_path.to_excel(writer, sheet_name='2')

  with pandas.ExcelWriter(path=os.path.join(top_genes_base_path,
                                            'TOP_HOMDEL_genes.xlsx')) as writer:
    # TODO(dlroxe): rownames=True
    top_homdel_genes.to_excel(writer, sheet_name='1')
    top_homdel_path.to_excel(writer, sheet_name='2')

  logging.info('"top genes" analyses have been written under %s.',
               top_genes_base_path)

  logging.info('attempting cooccurance analysis 1')
  parser.cooccurance_analysis(df=all_threshold_data, gene01="EIF4G1",
                              gene02="EIF3E", cnv_spec=["AMP", ])

  logging.info('attempting coocurrance analysis 2')
  parser.cooccurance_analysis(df=all_threshold_data, gene01="EIF4G1",
                              gene02="EIF3E", cnv_spec=["AMP", "DUP"])

  logging.info('attempting coocurrance analysis 3')
  parser.cooccurance_analysis(df=all_threshold_data, gene01="EIF4G1",
                              gene02="EIF3H", cnv_spec=["AMP", ])

  logging.info('attempting coocurrance analysis 4')
  parser.cooccurance_analysis(df=all_threshold_data, gene01="EIF4G1",
                              gene02="EIF3H", cnv_spec=["AMP", "DUP"])

  logging.info('processing complete')


if __name__ == "__main__":
  app.run(main)

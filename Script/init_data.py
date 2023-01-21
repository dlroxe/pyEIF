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
from typing import List

# TODO(dlroxe): This package structure is an improvement on the
#               prior use of "sys.path += ['./input_data_adapters']",
#               but probably still could use some thought.  It would be
#               nice, for example, not to have to resort to 'from' syntax.
from input_data_adapters import org_hs_eg_db_lookup
from input_data_adapters import tcga_cnv_parser

import datatable
import os
import pandas
import sys

FLAGS = flags.FLAGS
flags.DEFINE_string('data_directory',
                    # '~/Documents/Bioinformatics_analysis/eIF4G-analysis/eIF4G_data',
                    os.path.join('~', 'Desktop', 'pyeif_data'),
                    'parent dir for data files')
flags.DEFINE_string('output_directory',
                    # '~/Documents/Bioinformatics_analysis/eIF4G-analysis/eIF4G_output',
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


def _get_precomputed_paths():
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
  def read_top_path(file):
    loc = _abspath(os.path.join(FLAGS.output_directory, file))
    if os.path.exists(loc):
      return datatable.fread(loc).to_pandas()
    else:
      logging.warning('could not read: %s', loc)
      return None

  return (
    read_top_path(FLAGS.top_amp_path),
    read_top_path(FLAGS.top_gain_path),
    read_top_path(FLAGS.top_homdel_path),
  )


def _abspath(path):
  # What an absurd incantation to resolve "~"; but OK. Thanks, StackOverflow.
  return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))


def main(unused_argv):
  # TODO(dlroxe): The next line causes the logging module to write to stdout
  #               rather than to its traditional log file.
  #               Decide whether logging to stdout is desired in the long term.
  logging.get_absl_handler().python_handler.stream = sys.stdout

  eif_genes = ["EIF4G1", "EIF3E", "EIF3H", ]

  parser = tcga_cnv_parser.TcgaCnvParser(
    data_directory=FLAGS.data_directory,
    output_directory=FLAGS.output_directory,
    cnv_data_by_gene_thresholds=FLAGS.cnv_data_by_gene_thresholds,
    cnv_data_by_gene_values=FLAGS.cnv_data_by_gene_values,
    cnv_data_phenotypes=FLAGS.cnv_data_phenotypes,
  )
  # All the expensive datatable.fread() work, and parsing, for CNV data, is now
  # done.  The processed data can now be accessed repeatedly without significant
  # cost.

  all_data = parser.get_tcga_cnv_value()
  logging.info('all data\n%s', all_data)

  all_threshold_data = parser.get_tcga_cnv_threshold_categories()
  eif_threshold_data = all_threshold_data[eif_genes]
  logging.info('eif threshold data\n%s', eif_threshold_data)

  merged_phenotyped_data = parser.merge_cnv_thresholds_and_phenotypes()
  logging.info('all threshold data, merged with phenotypes:\n%s',
               merged_phenotyped_data)

  org_hs_eg_db_handle = org_hs_eg_db_lookup.OrgHsEgDbLookup(
    _abspath(os.path.join(FLAGS.data_directory, FLAGS.org_hs_eg_sqlite)))

  top_amp_paths, top_gain_paths, top_homdel_paths = _get_precomputed_paths()

  def top_genes(cnv_spec: List[str], percent: int) -> pandas.DataFrame:
    return parser.get_top_genes(cnv_spec, percent, org_hs_eg_db_handle)

  top_amp_genes = top_genes(["AMP"], 5)
  logging.info('top amp genes:\n%s', top_amp_genes)

  top_gain_genes = top_genes(["DUP", "AMP"], 30)
  logging.info('top gain genes:\n%s', top_gain_genes)

  top_homdel_genes = top_genes(["HOMDEL"], 5)
  logging.info('top homdel genes:\n%s', top_homdel_genes)

  def write_excel(name, genes, pathways):
    top_genes_base_path = os.path.join(FLAGS.output_directory, "Fig1")
    with pandas.ExcelWriter(
        path=os.path.join(top_genes_base_path, name)) as writer:
      # TODO(dlroxe): rownames=True
      genes.to_excel(writer, sheet_name='1')
      if pathways is not None:
        pathways.to_excel(writer, sheet_name='2')

  write_excel('TOP_AMP_genes.xlsx', top_amp_genes, top_amp_paths)
  write_excel('TOP_GAIN_genes.xlsx', top_gain_genes, top_gain_paths)
  write_excel('TOP_HOMDEL_genes.xlsx', top_homdel_genes, top_homdel_paths)

  logging.info('co-occurrence analysis: 4G1,3E + AMP')
  parser.co_occurrence_analysis("EIF4G1", "EIF3E", ["AMP", ])

  logging.info('co-occurrence analysis: 4G1,3E + AMP,DUP')
  parser.co_occurrence_analysis("EIF4G1", "EIF3E", ["AMP", "DUP"])

  logging.info('co-occurrence analysis 4G1,3H + AMP')
  parser.co_occurrence_analysis("EIF4G1", "EIF3H", ["AMP", ])

  logging.info('attempting co-occurrence analysis 4G1,3H + AMP,DUP')
  parser.co_occurrence_analysis("EIF4G1", "EIF3H", ["AMP", "DUP"])

  logging.info('processing complete')


if __name__ == "__main__":
  app.run(main)

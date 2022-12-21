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

FLAGS = flags.FLAGS
flags.DEFINE_string('data_directory', '~/Desktop/pyeif_data', 'parent dir for data files')
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
  def get_tcga_cnv(cls, values_data_frame: Optional[datatable.Frame] = None) -> pandas.DataFrame:
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
  def merge_cnv_phenotypes(cls, cnv_data: Optional[pandas.DataFrame] = None) -> pandas.DataFrame:
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

    phenotypes = datatable.fread(file=os.path.join(FLAGS.data_directory, FLAGS.cnv_data_phenotypes)
                                 ).to_pandas()[['sample', 'sample_type', '_primary_disease']].rename(
      columns={
        'sample': 'Sample',
        'sample_type': 'sample.type',
        '_primary_disease': 'primary_disease',
      }).sort_values(by=['Sample']).set_index('Sample')

    return cnv.join(phenotypes, how='inner')


# TODO(dlroxe): most code below this point is commented-out with docstring-style quotes, until it can be
# translated from R to Python.  After that is done, it probably should/will be reorganized into a class
# structure, as well.

def get_top_genes(df, label, percent):
  pass
  """
  sample_number <- nrow(df)

  TOP_AMP <- df %>%
    tibble::rownames_to_column(var = "rowname") %>%
    reshape2::melt(id.vars = "rowname", variable.name = "Gene") %>%
    dplyr::filter(value %in% label) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(Percent = n() / sample_number * 100) %>%
    dplyr::filter(Percent > percent) %>%
    droplevels() %>%
    dplyr::mutate(entrez = AnnotationDbi::mapIds(org.Hs.eg.db,
                                                 keys = as.character(.data$Gene),
                                                 column = "ENTREZID",
                                                 keytype = "SYMBOL",
                                                 multiVals = "first"
    ))
  return(TOP_AMP)
}
"""


def coocurrance_analysis(df, gene01, gene02, cnv_1, cnv_2):
  pass


"""
coocurrance_analysis <- function(df, gene01, gene02, cnv_1, cnv_2) {
  EIF <- df %>%
    dplyr::select(dplyr::all_of(c(gene01, gene02))) %>%
    mutate(
      !!gene01 := ifelse((!!as.name(gene01) == cnv_1 | !!as.name(gene01) == cnv_2),
        !!paste(gene01, cnv_1, cnv_2),
        !!paste(gene01, "NO", cnv_1, cnv_2)
      ),
      !!gene02 := ifelse((!!as.name(gene02) == cnv_1 | !!as.name(gene02) == cnv_2),
        !!paste(gene02, cnv_1, cnv_2),
        !!paste(gene02, "NO", cnv_1, cnv_2)
      )
    )
  file_name <- paste(file.path(output_directory, "Fig1"),
    "/",
    gene01,
    gene02,
    cnv_1,
    cnv_2,
    ".xlsx",
    sep = ""
  )
  xlsx::write.xlsx(as.data.frame.matrix(table(EIF[, gene01], EIF[, gene02])),
    file = file_name,
    sheetName = "1", row.names = TRUE
  )
  xlsx::write.xlsx(fisher.test(EIF[, gene01], EIF[, gene02], alternative = "greater") %>%
    broom::tidy(),
  file = file_name,
  sheetName = "Fisheroneside", append = TRUE, row.names = FALSE
  )
  xlsx::write.xlsx(chisq.test(EIF[, gene01], EIF[, gene02]) %>%
    broom::tidy(),
  file = file_name,
  sheetName = "chitest", append = TRUE, row.names = FALSE
  )
}

## function calling ============================================================
#
xlsx::write.xlsx2(get_top_genes(df = TCGA_CNV, label = "AMP", 5),
                  file = paste(file.path(output_directory, "Fig1"),
                               "/TOP_AMP_genes.xlsx",
                               sep = ""
                  ),
                  sheetName = "1", row.names = TRUE
)

xlsx::write.xlsx2(as.data.frame(TOP_AMP_PATH@result),
                  file = paste(file.path(output_directory, "Fig1"),
                               "/TOP_AMP_genes.xlsx",
                               sep = ""
                  ),
                  sheetName = "2", append = TRUE, row.names = FALSE
)

#
xlsx::write.xlsx2(get_top_genes(df = TCGA_CNV, label = c("DUP","AMP"), 30),
                  file = paste(file.path(output_directory, "Fig1"),
                               "/TOP_GAIN_genes.xlsx",
                               sep = ""
                  ),
                  sheetName = "1", row.names = TRUE
)

xlsx::write.xlsx2(as.data.frame(TOP_GAIN_PATH@result),
                  file = paste(file.path(output_directory, "Fig1"),
                               "/TOP_GAIN_genes.xlsx",
                               sep = ""
                  ),
                  sheetName = "2", append = TRUE, row.names = FALSE
)

#
xlsx::write.xlsx2(get_top_genes(df = TCGA_CNV, label = "HOMDEL", 5),
                  file = paste(file.path(output_directory, "Fig1"),
                               "/TOP_HOMDEL_genes.xlsx",
                               sep = ""
                  ),
                  sheetName = "1", row.names = TRUE
)

xlsx::write.xlsx2(as.data.frame(TOP_HOMDEL_PATH@result),
                  file = paste(file.path(output_directory, "Fig1"),
                               "/TOP_HOMDEL_genes.xlsx",
                               sep = ""
                  ),
                  sheetName = "2", append = TRUE, row.names = FALSE
)

##
coocurrance_analysis(df = TCGA_CNV,
                     gene01 = "EIF4G1",
                     gene02 = "EIF3E",
                     cnv_1 = "AMP",
                     cnv_2 = "AMP")  

coocurrance_analysis(df = TCGA_CNV,
                     gene01 = "EIF4G1",
                     gene02 = "EIF3E",
                     cnv_1 = "AMP",
                     cnv_2 = "DUP")  

coocurrance_analysis(df = TCGA_CNV,
                     gene01 = "EIF4G1",
                     gene02 = "EIF3H",
                     cnv_1 = "AMP",
                     cnv_2 = "AMP")  

coocurrance_analysis(df = TCGA_CNV,
                     gene01 = "EIF4G1",
                     gene02 = "EIF3H",
                     cnv_1 = "AMP",
                     cnv_2 = "DUP")  
"""


def main(argv):
  eif_genes = ["EIF4G1", "EIF3E", "EIF3H", ]

  all_data = TcgaCnvParser.get_tcga_cnv_value(raw_data_file=FLAGS.all_data_by_genes)
  print(f'all data\n{all_data}')

  raw_threshold_data = TcgaCnvParser.get_tcga_cnv_value(raw_data_file=FLAGS.cnv_data_by_gene_thresholds)
  all_threshold_data = TcgaCnvParser.get_tcga_cnv(values_data_frame=raw_threshold_data)
  eif_threshold_data = all_threshold_data[eif_genes]
  print(f'eif threshold data\n{eif_threshold_data}')

  merged_phenotyped_data = TcgaCnvParser.merge_cnv_phenotypes(all_threshold_data)
  print(f'all threshold data, merged with phenotypes:\n{merged_phenotyped_data}')


if __name__ == "__main__":
  app.run(main)

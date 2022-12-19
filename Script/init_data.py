#!/usr/bin/python3
import datatable
import os

# TODO(dlroxe): Make a command-line argument with this or some other value as a default.
data_file_directory = '~/Desktop/pyeif_data'

# TODO(dlroxe): Make a class, and store this as a class variable.

cnv_code_mappings = {
  2: 'AMP',
  1: 'DUP',
  0: 'DIPLOID',
  -1: 'DEL',
  -2: 'HOMDEL',
}


def get_tcga_cnv_value(genes_of_interest):
  """
  Reads 'Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes' and returns a related dataframe.

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

  :param genes_of_interest: a list of genes; in this example ['EIF4G1', 'EIF3E', 'EIF3H']
  :return: a data frame with samples as rows and selected genes as columns.
  """
  df = datatable.fread(
    file=os.path.join(data_file_directory,
                      'Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes')
  ).to_pandas().set_index('Sample').transpose()[genes_of_interest]
  return df


def get_tcga_cnv(genes_of_interest, values_data_frame=None):
  """
  Returns the output of get_tcga_cnv_value(), but with numeric cell values replaced by labels.

  Sample output for a selection of EIF genes:

  Sample            EIF4G1    EIF3E    EIF3H
  TCGA-A5-A0GI-01  DIPLOID  DIPLOID  DIPLOID
  TCGA-S9-A7J2-01  DIPLOID  DIPLOID  DIPLOID
  TCGA-06-0150-01  DIPLOID  DIPLOID  DIPLOID
  ...                  ...      ...      ...
  TCGA-DD-A115-01  DIPLOID      DEL      DEL

  :param genes_of_interest: a list of genes; in this example ['EIF4G1', 'EIF3E', 'EIF3H']
  :param values_data_frame: if None, then the function uses the value returned by get_tcga_cnv_value(genes_of_interest)
  :return: a data frame with samples as rows, selected genes as columns, and string labels as cell values.
  """
  df = values_data_frame or get_tcga_cnv_value(genes_of_interest)
  for col in genes_of_interest:
    for (code, label) in cnv_code_mappings.items():
      df[col].mask(df[col] == code, label, inplace=True)
  return df


def initialize_cnv_data():
  # rlang::env_binding_unlock(parent.env(environment()), nms = NULL)
  pass
  """
  CNV_path <- file.path(output_directory, "ProcessedData","TCGA_CNV.csv")
  CNV_value_path <- file.path(output_directory, "ProcessedData","TCGA_CNV_value.csv")

  # if (!file.exists(CNV_path)):
  # ...call get_tcga_cnv() and store the result
  # ...write the result of get_tcga_cnv() to the file
  # else:
  # ...read the data from the file instead of calling
  #    get_tcga_cnv()

  # if (!file.exists(CNV_value_path)):
  # ...call get_tcga_cnv_value() and store the result
  # ...write the result of get_tcga_cnv_value() to the file
  # else:
  # ...read the data from the file instead of calling
  #    get_tcga_cnv_value()

  TCGA_sampletype <- readr::read_tsv(file.path(
    data_file_directory,
    "TCGA_phenotype_denseDataOnlyDownload.tsv"
  )) %>%
    as.data.frame() %>%
    dplyr::distinct(.data$sample, .keep_all = TRUE) %>%
    stats::na.omit() %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "sample") %>%
    dplyr::select("sample_type", "_primary_disease") %>%
    dplyr::rename(
      "sample.type" = "sample_type",
      "primary.disease" = "_primary_disease"
    )

  # if (!file.exists(file.path(output_directory,"ProcessedData","TCGA_CNV_sampletype.csv"))):
  #
  #   assign("TCGA_CNV_sampletype",
  #          merge(TCGA_CNV,
  #                TCGA_sampletype,
  #                by    = "row.names",
  #                all.x = TRUE
  #          ) %>%
  #            dplyr::filter(.data$sample.type != "Solid Tissue Normal") %>%
  #            tibble::remove_rownames() %>%
  #            tibble::column_to_rownames(var = "Row.names"),
  #          envir = parent.env(environment())
  #   )
  # ...then write the file

  # else:
  # ...read the file and use it to assing TCGA_CNV_sampletype

  """


# initialize_cnv_data()

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


def main():
  eif_genes = ["EIF4G1", "EIF3E", "EIF3H", ]
  df = get_tcga_cnv(genes_of_interest=eif_genes)
  print(df)


if __name__ == "__main__":
  main()

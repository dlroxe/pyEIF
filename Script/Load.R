library(AnnotationDbi)
library(corrr)
library(cowplot)
library(data.table)
library(DescTools)
library(dplyr)
library(ggnewscale)
library(enrichplot)
library(forcats)
library(ggplot2)
library(likert)
library(magrittr) # %>%
library(org.Hs.eg.db)
library(ReactomePA)
library(reshape2)
library(stats)
library(stringr)
library(tibble)
library(tidyverse)
library(xlsx)

# "AnnotationDbi", "corrr", "cowplot", "data.table", "DescTools", "dplyr", "ggnewscale", "enrichplot", "forcats", "ggplot2", "likert", "magrittr", "org.Hs.eg.db", "ReactomePA", "reshape2", "stats", "stringr", "tibble", "tidyverse", "xlsx"


data_file_directory <- "~/Documents/Bioinformatics_analysis/eIF4G-analysis/eIF4G_data"
output_directory <- "~/Documents/Bioinformatics_analysis/eIF4G-analysis/eIF4G_output"

initialize_dir <- function() {
  dir.create(file.path(output_directory, "Fig1"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "Fig2"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "Fig3"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "Fig4"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "ProcessedData"),
             showWarnings = FALSE, recursive = TRUE)
}

initialize_dir()  

color <- function() {
  qual_col_pals <- RColorBrewer::brewer.pal.info[
    RColorBrewer::brewer.pal.info$category == "qual", ]
  
  return(unlist(mapply(
    RColorBrewer::brewer.pal,
    qual_col_pals$maxcolors,
    rownames(qual_col_pals)
  )))
}

initialize_format <- function() {
  
  assign("black_bold_tahoma_7",
         element_text(color = "black",
                      face = "bold",
                      size = 7),
         envir = parent.env(environment()))
  
  
  assign("black_bold_12",
         element_text(color = "black",
                      face = "bold",
                      size = 12),
         envir = parent.env(environment()))
  
  assign("black_bold_12_45",
         element_text(
           color = "black",
           face = "bold",
           size = 12,
           angle = 45,
           hjust = 1
         ),
         envir = parent.env(environment()))
  
  assign("black_bold_14",
         element_text(color = "black",
                      face = "bold",
                      size = 14),
         envir = parent.env(environment()))
  
  assign("black_bold_16",
         element_text(color = "black",
                      face = "bold",
                      size = 16),
         envir = parent.env(environment()))
  
  assign("black_bold_16_right",
         element_text(color = "black",
                      face = "bold",
                      size = 16,
                      angle = 90),
         envir = parent.env(environment()))
  
  assign("black_bold_16_45",
         element_text(
           color = "black",
           face = "bold",
           size = 16,
           angle = 45,
           hjust = 1
         ),
         envir = parent.env(environment()))
  
  assign("black_bold_16_90",
         element_text(
           color = "black",
           face = "bold",
           size = 16,
           angle = 90,
           hjust = 1,
           vjust = 0.5
         ),
         envir = parent.env(environment()))
  
  assign("black_bold_18",
         element_text(
           color = "black",
           face = "bold",
           size = 18
         ),
         envir = parent.env(environment()))
  
  assign("col_vector",
         color(),
         envir = parent.env(environment()))
}

initialize_format()

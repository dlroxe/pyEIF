library(AnnotationDbi)
library(corrr)
library(cowplot)
library(data.table)
library(DescTools)
library(dplyr)
library(forcats)
library(likert)
library(ggplot2)
library(magrittr) # %>%
library(org.Hs.eg.db)
library(stats)
library(tibble)
library(reshape2)
library(enrichplot)
library(ggnewscale)
library(stringr)
library(ReactomePA)
library(xlsx)

#data_file_directory <- "~/Downloads/NRF2_data"
data_file_directory <- "~/Documents/Bioinformatics_analysis/eIF4G-analysis/eIF4G_data"

output_directory <- "~/Documents/Bioinformatics_analysis/eIF4G-analysis/eIF4G_output"

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

## CNV analysis of all TCGA tumors =============================================
.get_TCGA_CNV <- function() {
  .TCGA_pancancer <- data.table::fread(
    file.path(
      data_file_directory,
      "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes"
    ),
    data.table = FALSE
  ) %>%
    tibble::as_tibble() %>%
    # as.data.frame(.) %>%
    dplyr::distinct(.data$Sample, .keep_all = TRUE) %>%
    na.omit() %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "Sample")
  
  # transpose function from the data.table keeps numeric values as numeric.
  .TCGA_pancancer_transpose <- data.table::transpose(.TCGA_pancancer)
  # get row and column names in order
  rownames(.TCGA_pancancer_transpose) <- colnames(.TCGA_pancancer)
  colnames(.TCGA_pancancer_transpose) <- rownames(.TCGA_pancancer)
  
  .TCGA_pancancer_transpose[.TCGA_pancancer_transpose == 2] <- "AMP"
  .TCGA_pancancer_transpose[.TCGA_pancancer_transpose == 1] <- "DUP"
  .TCGA_pancancer_transpose[.TCGA_pancancer_transpose == 0] <- "DIPLOID"
  .TCGA_pancancer_transpose[.TCGA_pancancer_transpose == -1] <- "DEL"
  .TCGA_pancancer_transpose[.TCGA_pancancer_transpose == -2] <- "HOMDEL"
  
  return(.TCGA_pancancer_transpose)
}

.get_TCGA_CNV_value <- function() {
  .TCGA_pancancer <- fread(
    file.path(
      data_file_directory,
      "Gistic2_CopyNumber_Gistic2_all_data_by_genes"
    ),
    data.table = FALSE
  ) %>%
    tibble::as_tibble() %>%
    # as.data.frame(.) %>%
    dplyr::distinct(.data$Sample, .keep_all = TRUE) %>%
    stats::na.omit() %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "Sample")
  
  # transpose function keeps numeric values as numeric.
  .TCGA_pancancer_transpose <- data.table::transpose(.TCGA_pancancer)
  # get row and colnames in order
  rownames(.TCGA_pancancer_transpose) <- colnames(.TCGA_pancancer)
  colnames(.TCGA_pancancer_transpose) <- rownames(.TCGA_pancancer)
  
  return(.TCGA_pancancer_transpose)
}

initialize_cnv_data <- function() {
  #rlang::env_binding_unlock(parent.env(environment()), nms = NULL)
  
  CNV_path <- file.path(output_directory, "ProcessedData","TCGA_CNV.csv")
  CNV_value_path <- file.path(output_directory, "ProcessedData","TCGA_CNV_value.csv")
  
  if (!file.exists(CNV_path)) {
    assign("TCGA_CNV",
           .get_TCGA_CNV(),
           envir = parent.env(environment())
    )
    data.table::fwrite(TCGA_CNV, CNV_path ,row.names = TRUE)
  } else {
    assign("TCGA_CNV",
           data.table::fread(CNV_path, data.table = FALSE, header = TRUE) %>%
             tibble::as_tibble() %>%
             #tibble::remove_rownames() %>%
             tibble::column_to_rownames(var = "V1"),
           envir = parent.env(environment())
    )
  }
  
  if (!file.exists(CNV_value_path)) {
    assign("TCGA_CNV_value",
           .get_TCGA_CNV(),
           envir = parent.env(environment())
    )
    data.table::fwrite(.get_TCGA_CNV_value(), 
                       CNV_value_path, 
                       row.names = TRUE)
  }  else {
  assign("TCGA_CNV_value",
         data.table::fread(CNV_value_path,
                           data.table = FALSE) %>%
           tibble::as_tibble() %>%
           #tibble::remove_rownames() %>%
           tibble::column_to_rownames(var = "V1"),
         envir = parent.env(environment())
  )
  }
  
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
  
  if (!file.exists(file.path(
    output_directory,
    "ProcessedData",
    "TCGA_CNV_sampletype.csv"
  ))) {
    assign("TCGA_CNV_sampletype",
           merge(TCGA_CNV,
                 TCGA_sampletype,
                 by    = "row.names",
                 all.x = TRUE
           ) %>%
             dplyr::filter(.data$sample.type != "Solid Tissue Normal") %>%
             tibble::remove_rownames() %>%
             tibble::column_to_rownames(var = "Row.names"),
           envir = parent.env(environment())
    )
    
    data.table::fwrite(TCGA_CNV_sampletype, file.path(
      output_directory,
      "ProcessedData",
      "TCGA_CNV_sampletype.csv"
    ),row.names = TRUE)
    
  } else {
    
    assign("TCGA_CNV_sampletype",
           data.table::fread(file.path(
             output_directory,
             "ProcessedData",
             "TCGA_CNV_sampletype.csv"),data.table = FALSE, header = TRUE
           ) %>%
             tibble::as_tibble() %>%
             #tibble::remove_rownames() %>%
             tibble::column_to_rownames(var = "V1"),
           envir = parent.env(environment())
    )
  }
  
  #rlang::env_binding_lock(parent.env(environment()), nms = NULL)
  
  return(NULL)
}

initialize_cnv_data()

get_top_genes <- function(df, label, percent) {
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

## pathway enrichment plot to analyze top CNV genes ============================
plot_enriched_pathway <- function(top_genes, label) {
  
  TOP_AMP_PATH <- unname(dplyr::pull(top_genes, entrez)) %>%
    ReactomePA::enrichPathway(pvalueCutoff = 0.15, 
                              readable = TRUE) 
  
  p1 <- ggplot(
    TOP_AMP_PATH,
    showCategory = 10,
    aes(GeneRatio, forcats::fct_reorder(Description, GeneRatio))
  ) +
    geom_segment(aes(xend = 0, yend = Description)) +
    geom_point(aes(color = p.adjust, size = Count)) +
    #scale_colour_gradient2() +
    scale_color_viridis_c(guide = guide_colorbar(reverse = TRUE)) +
    scale_size_continuous(range = c(2, 10)) +
    xlab("Gene Ratio") +
    ylab(NULL) +
    ggtitle(paste("The Most Enriched Pathways in Top", label, "Genes")) +
    # Modify labels of ggplot2 barplot
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 25)) +
    theme_bw() +
    theme(
      plot.title = black_bold_16,
      axis.title = black_bold_16,
      axis.text.x = black_bold_16,
      axis.text.y = black_bold_16,
      panel.grid = element_blank(),
      legend.title = black_bold_16,
      legend.text = black_bold_16,
      strip.text = black_bold_16
    )
  
  ggplot2::ggsave(
    path = file.path(output_directory, "Fig1"),
    filename = paste("Top", label, "dotplot.pdf"),
    plot = p1,
    width = 9,
    height = 9,
    useDingbats = FALSE)
  
  p2 <- clusterProfiler::cnetplot(TOP_AMP_PATH,
                                  node_label = "all",
                                  categorySize = "pvalue",
                                  colorEdge = TRUE,
                                  foldChange = dplyr::pull(top_genes, 
                                                           var = Percent, 
                                                           name = entrez),
                                  showCategory = 8,
                                  cex_label_category = 1,
                                  cex_gene	= 1.5,
                                  cex_label_gene = 0.8
  ) +
    guides(edge_color = "none")
  
  ggplot2::ggsave(
    path = file.path(output_directory, "Fig1"),
    filename = paste("Top", label, "cnetplot.pdf"),
    plot = p2,
    width = 9,
    height = 9,
    useDingbats = FALSE)
  
  plot2by2 <- cowplot::plot_grid(p1, p2, nrow = 2, rel_heights = c(1,1.2), labels = LETTERS[1:2])
  print(plot2by2)
  ggplot2::ggsave(
    path = file.path(output_directory, "Fig1"),
    filename = paste("Top", label, "in all TCGA.pdf"),
    plot = plot2by2,
    width = 10,
    height = 18,
    useDingbats = FALSE
  )
  
  return(TOP_AMP_PATH)
}

TOP_AMP_PATH <- plot_enriched_pathway(
  top_genes = get_top_genes(TCGA_CNV, c("AMP"), 5),
  label = "AMP"
)

TOP_GAIN_PATH <- plot_enriched_pathway(
  top_genes = get_top_genes(TCGA_CNV, c("DUP", "AMP"), percent = 30),
  label = "Gain"
)

TOP_HOMDEL_PATH <- plot_enriched_pathway(
  get_top_genes(TCGA_CNV, "HOMDEL", 5),
  "HOMDEL"
)

## matrix plot to compare CNV correlation ======================================
TOP_AMP_genes <- TOP_AMP_PATH@result %>% 
  dplyr::slice_head(n = 5) %>% 
  as.data.frame() %>%
  dplyr::select("Description","geneID") %>%
  dplyr::mutate(geneID = strsplit(geneID,"/"))

TOP_GAIN_genes <- TOP_GAIN_PATH@result %>% 
  dplyr::slice_head(n = 10) %>% 
  as.data.frame() %>%
  dplyr::select("Description","geneID") %>%
  dplyr::mutate(geneID = strsplit(geneID,"/"))

TOP_HOMDEL_genes <- TOP_HOMDEL_PATH@result %>% 
  dplyr::slice_head(n = 5) %>% 
  as.data.frame() %>%
  dplyr::select("Description","geneID") %>%
  dplyr::mutate(geneID = strsplit(geneID,"/"))

matrix_plot_CNV <- function(df, label) {
  # M <- stats::cor(df)
  testRes <- corrplot::cor.mtest(df, conf.level = 0.95)
  
  p1 <- corrplot::corrplot(stats::cor(df),
                           method = "color",
                           outline = T,
                           addgrid.col = "darkgray",
                           order = "hclust",
                           # addrect = 5,
                           rect.col = "black",
                           rect.lwd = 5, cl.pos = "b",
                           tl.col = "black",
                           tl.cex = 1,
                           cl.cex = 1,
                           addCoef.col = "white",
                           number.digits = 2,
                           number.cex = 0.5,
                           col = colorRampPalette(c("darkred", "white", "midnightblue"))(100)
  )
  print(p1) # print correlation matrix on the screen# save correlation plot as a pdf file
  pdf(
    file.path(output_directory, "Fig1", 
              paste("top", label, "corrmatrix.pdf")),
    width = 9,
    height = 9,
    useDingbats = FALSE
  )
  corrplot::corrplot(stats::cor(df),
                     method = "color",
                     outline = T,
                     addgrid.col = "darkgray",
                     order = "hclust",
                     # addrect = 5,
                     rect.col = "black",
                     rect.lwd = 5, cl.pos = "b",
                     tl.col = "black",
                     tl.cex = 1,
                     cl.cex = 1,
                     addCoef.col = "white",
                     number.digits = 2,
                     number.cex = 0.5,
                     col = colorRampPalette(c("darkred", "white", "midnightblue"))(100)
  )
  dev.off()
  
  return(NULL)
}

## Likert plot to compare CNV across cancer types ==============================
plot_CNV_Likert <- function(df, GeneID) {
  test1 <- df %>%
    dplyr::select(dplyr::all_of(GeneID), "primary.disease") %>%
    dplyr::mutate_if(is.character, as.factor) %>%
    dplyr::mutate_at(
      dplyr::vars(tidyselect::matches(GeneID)),
      forcats::fct_relevel,
      c("AMP", "DUP", "DIPLOID", "DEL", "HOMDEL")
    ) %>%
    tibble::rownames_to_column(var = "rowname") %>%
    reshape2::dcast( # primary.disease ~ !!as.name(GeneID),
      as.formula(paste("primary.disease ~", GeneID)),
      value.var = "primary.disease",
      fun.aggregate = length
    ) 
  
  p1 <- HH::likert(primary.disease ~ .,
                   data = test1,
                   positive.order = TRUE,
                   ylab = NULL,
                   main = paste(GeneID, "CNV status"),
                   # auto.key = list(columns = 1, reverse.rows = T),
                   as.percent = T,
                   # xlim=c(-40,-20,0,20,40,60,80,100),
                   # borders = list(),# <- This draws borders around the bars in the plot
                   strip = FALSE,
                   par.strip.text = list(cex = .7))
  
  print(p1)
  
  pdf(file = file.path(paste(file.path(output_directory, "Fig1")), 
                       paste0(GeneID, "_CNV_pancancer.pdf")),
      width = 9,
      height = 9,
      useDingbats = FALSE
  )
  print(p1) ### important to include print, otherwise the resulting plot.pdf is unreable
  
  dev.off()
  
  return(NULL)
  
}

plot_CNV_Likert(df = TCGA_CNV_sampletype, GeneID = "EIF4A2")

lapply(c("EIF4G1", "PABPC1", "EIF2B5", "EIF3E", "EIF3H", "EIF4A1","EIF4E"), 
       df = TCGA_CNV_sampletype, plot_CNV_Likert)


## bar plot to compare CNV percentage ==========================================
plot_Top_Genes_CNV <- function(df, label, string, color) {
  
  sample_number <- nrow(df)
  
  TOP_Genes <- df %>%
    tibble::rownames_to_column(var = "rowname") %>%
    melt(id.vars = "rowname", , variable.name = "Gene") %>%
    dplyr::filter(value == label) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(Percent = n() / sample_number * 100) %>%
    dplyr::filter(Percent > 5) %>%
    droplevels() %>%
    dplyr::filter(stringr::str_detect(Gene, string)) %>%
    dplyr::arrange(dplyr::desc(Percent)) 
  
  p1 <- ggplot(TOP_Genes, aes(x = reorder(Gene, -Percent), y = Percent)) +
    geom_bar(stat = "identity", width = 0.75, fill = color) +
    geom_text(aes(label = round(Percent, digits = 2)), vjust = -0.3) +
    labs(
      # x = "Translation initiation genes",
      y = paste(label, "frequency in all TCGA cancer %")
    ) +
    theme_bw() +
    theme(
      plot.title = black_bold_12,
      axis.title.x = element_blank(),
      axis.title.y = black_bold_12,
      axis.text.x = black_bold_12_45,
      axis.text.y = black_bold_12,
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.text = element_blank()
    )
  print(p1)
  ggplot2::ggsave(
    path = file.path(output_directory, "Fig1"),
    filename = paste("top", string, label, "in pancancer.pdf"),
    plot = p1,
    width = 7,
    height = 7,
    useDingbats = FALSE
  )
}

## co-occurrence  ======================================
CNV_vennplot <- function(df, gene01, gene02, gene03, gene04, cnv_1, cnv_2) {
  EIF <- df[, c(gene01, gene02, gene03, gene04)] %>%
  dplyr::mutate(
    !!gene01 := (!!as.name(gene01) == cnv_1 | !!as.name(gene01) == cnv_2),
     !!gene02 := (!!as.name(gene02)  == cnv_1 | !!as.name(gene02)  == cnv_2),
     !!gene03 := (!!as.name(gene03)  == cnv_1 | !!as.name(gene03)  == cnv_2),
     !!gene04 := (!!as.name(gene04)  == cnv_1 | !!as.name(gene04)  == cnv_2))
    # ":", "\n" split the final ratio name in two lines
    # for the plotting purpose.
  b <- limma::vennCounts(EIF)
  limma::vennDiagram(b)
  ## eulerr generates area-proportional Euler diagrams that display set
  ## relationships (intersections, unions, and disjoints) with circles.
  pos.Venn2 <- eulerr::euler(
    c(
      EIF4G1 = b[9, "Counts"], # EIF4E
      EIF4A2 = b[5, "Counts"], # EIF4G1
      EIF3E = b[3, "Counts"], # EIF4A1
      EIF3H = b[2, "Counts"], # EIF4EBP1
      "EIF4G1&EIF4A2" = b[13, "Counts"],
      "EIF4G1&EIF3E" = b[11, "Counts"],
      "EIF4G1&EIF3H" = b[10, "Counts"],
      "EIF4A2&EIF3E" = b[7, "Counts"],
      "EIF4A2&EIF3H" = b[6, "Counts"],
      "EIF3E&EIF3H" = b[4, "Counts"],
      "EIF4G1&EIF4A2&EIF3E" = b[15, "Counts"],
      "EIF4G1&EIF4A2&EIF3H" = b[14, "Counts"],
      "EIF4G1&EIF3E&EIF3H" = b[12, "Counts"],
      "EIF4A2&EIF3E&EIF3H" = b[8, "Counts"],
      "EIF4G1&EIF4A2&EIF3E&EIF3H" = b[16, "Counts"]
    ),
    # shape = "ellipse"
  )
  p2 <- plot(pos.Venn2,
    # key = TRUE,
    # main = paste(tissue_type, sample_type, CORs_type),
    lwd = 0,
    fill = c(
      "#999999", "#009E73",
      "#56B4E9", "#E69F00"
    ),
    # quantities = list(type = c("percent", "counts")),
    quantities = list(cex = 1.25),
    labels = list(
      labels = c(gene01, gene02, gene03, gene04),
      cex = 1.25
    )
  )
  print(p2)
  ggplot2::ggsave(
    path = file.path(output_directory, "Fig1"),
    filename = paste0(cnv_1, cnv_2, "Venn.pdf"),
    plot = p2,
    width = 6,
    height = 6,
    useDingbats = FALSE
  )
}

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
matrix_plot_CNV(df = TCGA_CNV_value %>%
                  dplyr::select(dplyr::all_of(
                    c("MYC",
                      unique(unlist(TOP_AMP_genes$geneID))))), 
                label = "AMP")

matrix_plot_CNV(df = TCGA_CNV_value %>%
                  dplyr::select(dplyr::all_of(
                    c("CDKN2A", "CDKN2B", 
                      unique(unlist(TOP_HOMDEL_genes$geneID))))), 
                label = "HOMDEL")

plot_Top_Genes_CNV(
  df = TCGA_CNV,
  label = "AMP",
  string = "RPL|EIF|PABPC1|MYC",
  color = "red"
)

plot_Top_Genes_CNV(
  df = TCGA_CNV,
  label = "HOMDEL",
  string = "IFNA|CDKN2A|CDKN2B",
  color = "blue"
)

##
CNV_vennplot(df = TCGA_CNV, 
             gene01 = "EIF4G1", 
             gene02 = "EIF4A2", 
             gene03 = "EIF3E", 
             gene04 = "EIF3H", 
             cnv_1  = "AMP", 
             cnv_2  = "DUP")

CNV_vennplot(df = TCGA_CNV, 
             gene01 = "EIF4G1", 
             gene02 = "EIF4A2", 
             gene03 = "EIF3E", 
             gene04 = "EIF3H", 
             cnv_1  = "AMP", 
             cnv_2  = "AMP")

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

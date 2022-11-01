library(shiny)
library(shinyjs)
library(tidyverse)
library(readxl)
library(shinythemes)
library(svglite)
library(pheatmap)
library(Seurat)
library(SingleCellExperiment)
library(scater)
# library(gridExtra)


# import sample list , raw gene count and TMM normalized expression matrix
genes.counts <- read.delim("./merged_counts/genes.counts.matrix", row.names = 1, check.names = FALSE )
genes.TMM.EXPR <-
  read.delim(
    "./merged_counts/genes.TMM.EXPR.matrix",
    row.names = 1,
    check.names = FALSE
  )

sample_info <-
  read.table("./merged_counts/sample.info",
    row.names = "sample",
    header = TRUE
  )
sample_info_exp <- read_csv("./merged_counts/sample_info_exp.csv")

WT_FL_0_2day_TMM_sample_exp <- read_csv("./merged_counts/WT_FL_0_2day_TMM_sample_exp.csv")

iso_exp_tpm <- read_delim(
  "data/isoforms(HUA-ccNET).fpkm_table",
  "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)
AD_pro_19920 <- read_excel("data/AD-pro-19920.xlsx")
blastp_AD1_HAU_v1_0_vs_arabidopsis_1 <-
  read_delim(
    "./bla_go/blastp_AD1_HAU_v1.0_vs_arabidopsis.1.txt",
    "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )
genes2Go <- read_excel("./bla_go/genes2Go.xlsx",
  col_names = T,
  skip = 1
)

# careate an integrated exp and info table
# tran_exp <- as_tibble(t(genes.TMM.EXPR), rownames = 'sample')
### arrange by time
sample_info <-
  sample_info %>%
  rownames_to_column(var = "sample") %>%
  group_by(strain) %>%
  arrange(time, .by_group = TRUE)
#### trans time to labels
trans_time <- function(str, pattr) {
  mat_obj <- str_match(str, pattern = pattr)
  times <- if_else(mat_obj[[3]] == "minus", "-", "")
  day <- if_else(mat_obj[[5]] == "12", "N", "L")
  out_str <- str_c(times, mat_obj[[4]], "DPA", "-", day)
  return(out_str)
}
sample_info <-
  sample_info %>%
  ungroup() %>%
  mutate(labs = map_chr(
    .$sample,
    trans_time, "(\\w+)-(\\w*)(\\d)-DPA(\\d+)h-(\\d)"
  ))

genes.TMM.EXPR <- genes.TMM.EXPR %>% select(sample_info$sample)

## single cell
# load("data/mergeSCE.RData") 转移到另一个硬盘
load("/data1/wangy/data/R_single_cell/mergeSCE.RData")
load("data/all_strain_sce.Rdata")

## ID convert
WHU <- read.delim(file = "./data/Gh_whu_to_hau.bed", header = F)
colnames(WHU) <- c("Chrom_In", "Start_In", "End_In", "ID_In", "Identity", "Strand_In",
                   "Chrom_Out", "Start_Out", "End_Out", "ID_Out", "Score", "Strand_Out")
ZJU <- read.delim(file = "./data/Gh_zju_to_hau.bed", header = F)
colnames(ZJU) <- c("Chrom_In", "Start_In", "End_In", "ID_In", "Identity", "Strand_In",
                   "Chrom_Out", "Start_Out", "End_Out", "ID_Out", "Score", "Strand_Out")
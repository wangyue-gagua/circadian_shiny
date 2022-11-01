library(tidyverse)
library(circacompare)
library(ggplot2)
library(data.table)
library(MetaCycle)
setwd("WT_FL_0_2dpa")
# salmon data import
# salmon_WT_FL_0_2day_TMM_sample_exp <- read_csv("salmon_WT_FL_0_2day_TMM_sample_exp.csv")
salmon_WT_FL_0_2day_TMM_sample_exp <- fread("salmon_WT_FL_0_2day_TMM_sample_exp.csv")
salmon_WT_FL_0_2day_TMM_sample_exp <- salmon_WT_FL_0_2day_TMM_sample_exp %>% as_tibble()

salmon_WT_FL_0_2day_genes_TMM_EXPR <- read_delim("salmon_WT_FL_0_2day_genes.TMM.EXPR.matrix",
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
) %>% column_to_rownames(var = "transcriptId")

salmon_WT_0_2day_genes_TMM_EXPR <- salmon_WT_FL_0_2day_genes_TMM_EXPR %>%
  select(contains("WT"))

salmon_FL_0_2day_genes_TMM_EXPR <- salmon_WT_FL_0_2day_genes_TMM_EXPR %>%
  select(contains("FL"))

salmon_WT_0_2day_genes_TMM_EXPR_mergeRep <- salmon_WT_0_2day_genes_TMM_EXPR
salmon_FL_0_2day_genes_TMM_EXPR_mergeRep <- salmon_FL_0_2day_genes_TMM_EXPR


sample_list <- salmon_WT_FL_0_2day_TMM_sample_exp %>%
  group_by(strain) %>%
  arrange(time) %>%
  select(sample)


for (i in 1:length(sample_list$strain)) {
  # "FL_0DPA_1am1.quant"
  name <- str_remove(sample_list$sample[[i]], "\\d\\.quant$")
  if (sample_list$strain[[i]] == "WT") {
    salmon_WT_0_2day_genes_TMM_EXPR_mergeRep <- salmon_WT_0_2day_genes_TMM_EXPR_mergeRep %>%
      mutate(!!name := rowMeans(select(salmon_WT_0_2day_genes_TMM_EXPR, contains(name))))
  }
  if (sample_list$strain[[i]] == "FL") {
    salmon_FL_0_2day_genes_TMM_EXPR_mergeRep <- salmon_FL_0_2day_genes_TMM_EXPR_mergeRep %>%
      mutate(!!name := rowMeans(select(salmon_FL_0_2day_genes_TMM_EXPR, contains(name))))
  }
}

salmon_WT_0_2day_genes_TMM_EXPR_mergeRep_selected <- salmon_WT_0_2day_genes_TMM_EXPR_mergeRep[37:length(salmon_WT_0_2day_genes_TMM_EXPR_mergeRep)] %>%
  rownames_to_column(var = "transcriptId")
salmon_FL_0_2day_genes_TMM_EXPR_mergeRep_selected <- salmon_FL_0_2day_genes_TMM_EXPR_mergeRep[37:length(salmon_FL_0_2day_genes_TMM_EXPR_mergeRep)] %>%
  rownames_to_column(var = "transcriptId")

library(MetaCycle)
write.csv(salmon_WT_0_2day_genes_TMM_EXPR_mergeRep_selected, file = "salmon_WT_0_2day_genes_TMM_EXPR_mergeRep_selected.csv", row.names = F)
write.csv(salmon_FL_0_2day_genes_TMM_EXPR_mergeRep_selected, file = "salmon_FL_0_2day_genes_TMM_EXPR_mergeRep_selected.csv", row.names = F)
meta2d(
  infile = "salmon_WT_0_2day_genes_TMM_EXPR_mergeRep_selected.csv",
  filestyle = "csv",
  timepoints = seq(1, 69, 4),
  outdir = "salmon_WT_meta2d",
  parallelize = TRUE,
  nCores = 10,
  minper = 24,
  maxper = 24,
  ARSdefaultPer = 24,
  cycMethod = "JTK",
)

meta2d(
  infile = "salmon_FL_0_2day_genes_TMM_EXPR_mergeRep_selected.csv",
  filestyle = "csv",
  timepoints = seq(1, 69, 4),
  outdir = "salmon_FL_meta2d",
  parallelize = TRUE,
  nCores = 10,
  minper = 24,
  maxper = 24,
  ARSdefaultPer = 24,
  cycMethod = "JTK",
)

salmon_WT_meta2d <- read_csv("salmon_WT_meta2d/JTKresult_salmon_WT_0_2day_genes_TMM_EXPR_mergeRep_selected.csv")
salmon_FL_meta2d <- read_csv("salmon_FL_meta2d/JTKresult_salmon_FL_0_2day_genes_TMM_EXPR_mergeRep_selected.csv")


## 将时间序列分割为两个部分分别为 seq(1, 45, 4)和 seq(25, 69, 4)
salmon_WT_0_2day_genes_TMM_EXPR_mergeRep_selected <- read_csv("salmon_WT_0_2day_genes_TMM_EXPR_mergeRep_selected.csv")
salmon_FL_0_2day_genes_TMM_EXPR_mergeRep_selected <- read_csv("salmon_FL_0_2day_genes_TMM_EXPR_mergeRep_selected.csv")

salmon_WT_0_2day_genes_TMM_EXPR_mergeRep_selected_section1 <- salmon_WT_0_2day_genes_TMM_EXPR_mergeRep_selected %>% 
select(c(1, 2:13))
salmon_WT_0_2day_genes_TMM_EXPR_mergeRep_selected_section2 <- salmon_WT_0_2day_genes_TMM_EXPR_mergeRep_selected %>%
select(c(1, 8:19))

salmon_FL_0_2day_genes_TMM_EXPR_mergeRep_selected_section1 <- salmon_FL_0_2day_genes_TMM_EXPR_mergeRep_selected %>%
select(c(1, 2:13))
salmon_FL_0_2day_genes_TMM_EXPR_mergeRep_selected_section2 <- salmon_FL_0_2day_genes_TMM_EXPR_mergeRep_selected %>%
select(c(1, 8:19))

write.csv(salmon_WT_0_2day_genes_TMM_EXPR_mergeRep_selected_section1, file = "salmon_WT_0_2day_genes_TMM_EXPR_mergeRep_selected_section1.csv", row.names = F)
write.csv(salmon_WT_0_2day_genes_TMM_EXPR_mergeRep_selected_section2, file = "salmon_WT_0_2day_genes_TMM_EXPR_mergeRep_selected_section2.csv", row.names = F)
write.csv(salmon_FL_0_2day_genes_TMM_EXPR_mergeRep_selected_section1, file = "salmon_FL_0_2day_genes_TMM_EXPR_mergeRep_selected_section1.csv", row.names = F)
write.csv(salmon_FL_0_2day_genes_TMM_EXPR_mergeRep_selected_section2, file = "salmon_FL_0_2day_genes_TMM_EXPR_mergeRep_selected_section2.csv", row.names = F)

meta2d(
  infile = "salmon_WT_0_2day_genes_TMM_EXPR_mergeRep_selected_section1.csv",
  filestyle = "csv",
  timepoints = seq(1, 45, 4),
  outdir = "salmon_WT_meta2d_section1",
  parallelize = TRUE,
  nCores = 10,
  minper = 24,
  maxper = 24,
  ARSdefaultPer = 24,
  cycMethod = "JTK",
)

meta2d(
  infile = "salmon_WT_0_2day_genes_TMM_EXPR_mergeRep_selected_section2.csv",
  filestyle = "csv",
  timepoints = seq(25, 69, 4),
  outdir = "salmon_WT_meta2d_section2",
  parallelize = TRUE,
  nCores = 10,
  minper = 24,
  maxper = 24,
  ARSdefaultPer = 24,
  cycMethod = "JTK",
)

meta2d(
  infile = "salmon_FL_0_2day_genes_TMM_EXPR_mergeRep_selected_section1.csv",
  filestyle = "csv",
  timepoints = seq(1, 45, 4),
  outdir = "salmon_FL_meta2d_section1",
  parallelize = TRUE,
  nCores = 10,
  minper = 24,
  maxper = 24,
  ARSdefaultPer = 24,
  cycMethod = "JTK",
)

meta2d(
  infile = "salmon_FL_0_2day_genes_TMM_EXPR_mergeRep_selected_section2.csv",
  filestyle = "csv",
  timepoints = seq(25, 69, 4),
  outdir = "salmon_FL_meta2d_section2",
  parallelize = TRUE,
  nCores = 10,
  minper = 24,
  maxper = 24,
  ARSdefaultPer = 24,
  cycMethod = "JTK",
)

salmon_WT_meta2d_section1 <- read_csv("salmon_WT_meta2d_section1/JTKresult_salmon_WT_0_2day_genes_TMM_EXPR_mergeRep_selected_section1.csv")
salmon_WT_meta2d_section2 <- read_csv("salmon_WT_meta2d_section2/JTKresult_salmon_WT_0_2day_genes_TMM_EXPR_mergeRep_selected_section2.csv")
salmon_FL_meta2d_section1 <- read_csv("salmon_FL_meta2d_section1/JTKresult_salmon_FL_0_2day_genes_TMM_EXPR_mergeRep_selected_section1.csv")
salmon_FL_meta2d_section2 <- read_csv("salmon_FL_meta2d_section2/JTKresult_salmon_FL_0_2day_genes_TMM_EXPR_mergeRep_selected_section2.csv")

salmon_WT_meta2d_section12 <- full_join(salmon_WT_meta2d_section1, salmon_WT_meta2d_section2, by = "CycID", suffix = c("_section1", "_section2"))
Ghir_D01G016160
salmon_WT_meta2d_section12 %>% filter(str_detect(CycID, "Ghir_D01G016160"))
salmonPlotRepCirca("Ghir_D01G016160.1")
# cosinor2 Comparison of cosinor parameters of two populations

library(cosinor2)
library(ggplot2) 
library(circacompare)
set.seed(42)

metaInfo_salmon_WT_FL_0_2day_TMM_sample_exp <- salmon_WT_FL_0_2day_TMM_sample_exp %>% select(1:6)

salmonPlotRepCirca <- function(str, alia_name = "") {
    df <- cbind(metaInfo_salmon_WT_FL_0_2day_TMM_sample_exp, measure = salmon_WT_FL_0_2day_TMM_sample_exp[str])
    colnames(df)[7] <- "measure"


rects <-
    data.frame(
        xstart = c(14.5, 38.5, 62.5),
        xend = c(23, 47, 71)
    )
ggplot(data = df, aes(time, measure)) +
    geom_point(aes(col = strain)) +
    geom_smooth(aes(group = interaction(as.factor(replicate), strain), color = strain), span = 0.3) +
    facet_wrap(~replicate, nrow = 2) +
    scale_x_continuous(breaks = seq(1, 69, 4), labels = df$labs[1:(length(df$labs) /
        4)]) +
    geom_rect(
        data = rects,
        aes(
            xmin = xstart,
            xmax = xend,
            ymin = 0,
            ymax = Inf
        ),
        inherit.aes = FALSE,
        alpha = 0.2
    ) +
    labs(title = str, subtitle = alia_name) +
    ylab("expression level (TMM)") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
}
library(grid)
salmon_grid_plot_WT_FL_0_2dpa <- function(a_chr_list, alias=NULL){
  ## input a id list and a gene name, output a grid plot
  gobs <- map(a_chr_list, salmonPlotRepCirca, alias)
  new_gobs <- map(gobs[seq_along(gobs)], `+`, theme(axis.title.x = element_blank(), axis.text.x = element_blank()))
  new_gobs[-1] <- gobs[-1]
  
  gridExtra::grid.arrange(grobs=new_gobs)
}


salmon_grid_plot_WT_FL_0_2dpa(filter(salmon_FL_WT_meta2d, geneId == "Ghir_A01G000210")$CycID, alias = "Ghir_A01G000210")
salmon_grid_plot_WT_FL_0_2dpa(filter(salmon_FL_WT_meta2d, geneId == "Ghir_D04G000830")$CycID)

# 
salmon_FL_WT_meta2d_hasIsoform_all_has_rhythm_duplicate_gene_tbl <- read_tsv("salmon_FL_WT_meta2d_hasIsoform_all_has_rhythm_duplicate_gene_tbl.tsv")

# import merged_counts/txiCountMatrix.tsv
WT_FL_0_2day_genes_TMM_EXPR <- read_ts
salmon_WT_FL_0_2day_count_sample_exp <- read_tsv("../merged_counts/txiCountMatrix.tsv")
head(salmon_WT_FL_0_2day_count_sample_exp)
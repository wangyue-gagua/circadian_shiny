plotRepCirca("Ghir_D12G009360")
my_cir_plot_WT_FL_0_2dpa("Ghir_D12G009360")

# WT_FL_0_2day_TMM_sample_exp_detrended <-  WT_FL_0_2day_TMM_sample_exp  %>% group_by(strain, replicate, period) %>%
#     mutate(across(where(is.double), ~ . / (mean(.) + .Machine$double.eps)))

write_csv(WT_FL_0_2day_TMM_sample_exp_detrended, "./mediumDataSave/WT_FL_0_2day_TMM_sample_exp_detrended.csv")
WT_FL_0_2day_TMM_sample_exp_detrended <- fread("./mediumDataSave/WT_FL_0_2day_TMM_sample_exp_detrended.csv") %>%
 as_tibble()

# randomSixCluster1Genens <- c("Ghir_D12G020190", "Ghir_A12G005280", "Ghir_A13G013150", "Ghir_D07G006230", "Ghir_D02G019940", "Ghir_D08G002570")

p1 <- plotRepCircaDetrended("Ghir_D12G020190")
p2 <- plotRepCircaDetrended("Ghir_A12G005280")
p3 <- plotRepCircaDetrended("Ghir_A13G013150")
p4 <- plotRepCircaDetrended("Ghir_D07G006230")
p5 <- plotRepCircaDetrended("Ghir_D02G019940")
p6 <- plotRepCircaDetrended("Ghir_D08G002570")

# patchwork collect legend
p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(guides = "collect")
ggsave("figure/geneRepCirca/WT_FL_0_2day_detrended_cluster1_Genes_6gene_detrended.pdf")

WT_0_2day_genes_TMM_EXPR_mergeRep_selected <- fread("WT_0_2day_genes_TMM_EXPR_mergeRep_selected.csv") %>% as_tibble()
head(WT_0_2day_genes_TMM_EXPR_mergeRep_selected)
head(WT_FL_0_2day_TMM_sample_exp_detrended)
# WT_FL_0_2day_genes_TMM_EXPR_detrended <- WT_FL_0_2day_TMM_sample_exp_detrended %>%
#     arrange(time) %>%
#     select(1, 7:ncol(WT_FL_0_2day_TMM_sample_exp_detrended)) %>%
#     column_to_rownames("sample") %>%
#     t() %>%
#     as.data.frame() %>%
#     rownames_to_column("Geneid")
# write_csv(WT_FL_0_2day_genes_TMM_EXPR_detrended, "./mediumDataSave/WT_FL_0_2day_genes_TMM_EXPR_detrended.csv")
WT_FL_0_2day_genes_TMM_EXPR_detrended <- fread("./mediumDataSave/WT_FL_0_2day_genes_TMM_EXPR_detrended.csv") %>%
    as_tibble() %>%
    column_to_rownames("Geneid")

sample_list <- WT_FL_0_2day_TMM_sample_exp %>%
    group_by(strain) %>%
    arrange(time) %>%
    select(sample)

WT_0_2day_genes_TMM_EXPR_detrended <- WT_FL_0_2day_genes_TMM_EXPR_detrended %>% select(contains("WT"))
FL_0_2day_genes_TMM_EXPR_detrended <- WT_FL_0_2day_genes_TMM_EXPR_detrended %>% select(contains("FL"))

WT_0_2day_genes_TMM_EXPR_mergeRep_detrended <- WT_0_2day_genes_TMM_EXPR_detrended
FL_0_2day_genes_TMM_EXPR_mergeRep_detrended <- FL_0_2day_genes_TMM_EXPR_detrended

for (i in 1:length(sample_list$strain)) {
    name <- str_remove(sample_list$sample[[i]], "\\d$")
    if (sample_list$strain[[i]] == "WT") {
        WT_0_2day_genes_TMM_EXPR_mergeRep_detrended <- WT_0_2day_genes_TMM_EXPR_mergeRep_detrended %>%
            mutate(!!name := rowMeans(select(WT_0_2day_genes_TMM_EXPR_detrended, contains(name))))
    }
    if (sample_list$strain[[i]] == "FL") {
        FL_0_2day_genes_TMM_EXPR_mergeRep_detrended <- FL_0_2day_genes_TMM_EXPR_mergeRep_detrended %>%
            mutate(!!name := rowMeans(select(FL_0_2day_genes_TMM_EXPR_detrended, contains(name))))
    }
}

WT_0_2day_genes_TMM_EXPR_mergeRep_selected_detrended <- WT_0_2day_genes_TMM_EXPR_mergeRep_detrended[37:length(WT_0_2day_genes_TMM_EXPR_mergeRep_detrended)] %>%
    rownames_to_column(var = "Geneid")
FL_0_2day_genes_TMM_EXPR_mergeRep_selected_detrended <- FL_0_2day_genes_TMM_EXPR_mergeRep_detrended[37:length(FL_0_2day_genes_TMM_EXPR_mergeRep_detrended)] %>%
    rownames_to_column(var = "Geneid")

# write_csv(WT_0_2day_genes_TMM_EXPR_mergeRep_selected_detrended, "./mediumDataSave/WT_0_2day_genes_TMM_EXPR_mergeRep_selected_detrended.csv")
# write_csv(FL_0_2day_genes_TMM_EXPR_mergeRep_selected_detrended, "./mediumDataSave/FL_0_2day_genes_TMM_EXPR_mergeRep_selected_detrended.csv")
WT_0_2day_genes_TMM_EXPR_mergeRep_selected_detrended <- read_csv("./mediumDataSave/WT_0_2day_genes_TMM_EXPR_mergeRep_selected_detrended.csv")
library(MetaCycle)
meta2d(
  infile = "./mediumDataSave/WT_0_2day_genes_TMM_EXPR_mergeRep_selected_detrended.csv",
  filestyle = "csv",
  timepoints = seq(1, 69, 4),
  outdir = "WT_meta2d_detrended",
  parallelize = TRUE,
  nCores = 10,
  minper = 24,
  maxper = 24,
  ARSdefaultPer = 24,
  cycMethod = "JTK",
)

meta2d(
  infile = "./mediumDataSave/FL_0_2day_genes_TMM_EXPR_mergeRep_selected_detrended.csv",
  filestyle = "csv",
  timepoints = seq(1, 69, 4),
  outdir = "FL_meta2d_detrended",
  parallelize = TRUE,
  nCores = 10,
  minper = 24,
  maxper = 24,
  ARSdefaultPer = 24,
  cycMethod = "JTK",
)

WT_meta2d_detrended <- read_csv("WT_meta2d_detrended/JTKresult_WT_0_2day_genes_TMM_EXPR_mergeRep_selected_detrended.csv")
FL_meta2d_detrended <- read_csv("FL_meta2d_detrended/JTKresult_FL_0_2day_genes_TMM_EXPR_mergeRep_selected_detrended.csv")


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

plotRepCirca("Ghir_D04G020730", "GhMYB25")
ggsave("figure/geneRepCirca/Ghir_D04G020730_GhMYB25.pdf", width = 10, height = 10)
plotRepCirca("Ghir_D12G017660", "GhMYB25-like")
ggsave("figure/geneRepCirca/Ghir_D12G017660_GhMYB25-like.pdf", width = 10, height = 10)
my_cir_plot_WT_FL_0_2dpa("Ghir_D12G017660", "GhMYB25-like")
plotRepCircaDetrended("Ghir_D12G017660", "GhMYB25-like")
ggsave("figure/geneRepCirca/Ghir_D12G017660_GhMYB25-like_detrended.pdf", width = 10, height = 10)

## MYB25-like MYB25相位差
WT_FL_0_2day_TMM_sample_exp %>%
    select(1:6, "Ghir_D04G020730", "Ghir_D12G017660") %>%
    pivot_longer(cols = c("Ghir_D04G020730", "Ghir_D12G017660"), names_to = "genes", values_to = "expression") %>%
    ggplot(aes(time, expression)) +
    geom_point(aes(col = strain)) +
    geom_smooth(aes(group = interaction(as.factor(replicate), strain, genes), color = strain, linetype = genes), span = 0.3) +
    facet_wrap(~replicate, nrow = 2) +
    scale_x_continuous(breaks = seq(1, 69, 4), ) +
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
    geom_vline(xintercept = c(6, 18), linetype = "dotted") +
    labs(title = "Phase diff between MYB25 and MYB25-like", subtitle = NULL) +
    ylab("relative expression/(TMM)") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave("figure/geneRepCirca/Ghir_D04G020730_Ghir_D12G017660_phaseDiff.pdf", width = 10, height = 10)

## 查找Detrend后metacycle检测结果的变化
WT_dat_detrended <- WT_meta2d_detrended %>% filter(BH.Q < 0.01)
FL_dat_detrended <- FL_meta2d_detrended %>% filter(BH.Q < 0.01)
WT_circ_genes_detrended <- WT_dat_detrended$CycID
FL_circ_genes_detrended <- FL_dat_detrended$CycID

WT_circ_genes_specific_detrended <- setdiff(WT_circ_genes_detrended, FL_circ_genes_detrended)
FL_circ_genes_specific_detrended <- setdiff(FL_circ_genes_detrended, WT_circ_genes_detrended)
WT_FL_circ_genes_common_detrended <- intersect(WT_circ_genes_detrended, FL_circ_genes_detrended)
writeLines(WT_circ_genes_specific_detrended, "mediumDataSave/WT_circ_genes_specific_detrended.txt") # 3441
writeLines(FL_circ_genes_specific_detrended, "mediumDataSave/FL_circ_genes_specific_detrended.txt") # 4965
writeLines(WT_FL_circ_genes_common_detrended, "mediumDataSave/WT_FL_circ_genes_common_detrended.txt") # 3370
library(VennDiagram)
venn.diagram(
    x = list(
        WT_detrended = WT_circ_genes_detrended,
        FL_detrended = FL_circ_genes_detrended,
        WT = WT_circ_genes,
        FL = FL_circ_genes
    ),
    col = "transparent",
    fill = c("blue", "green", "yellow", "grey50"),
    alpha = 0.50,
    cex = 1.2, # 每个区域label名称的大小
    cat.col = c("darkblue", "darkgreen", "orange", "grey50"), # 分类颜色
    cat.cex = 1.2, # 每个分类名称大小
    cat.dist = 0.07,
    cat.pos = 0, #
    cat.fontfamily = "serif", # 分类字体
    rotation.degree = 270, # 旋转角度
    margin = 0.2,
    filename = "figure/samplesCor/circ_genes.png",
    imagetype = "png",
)

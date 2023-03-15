# 对第一次测序数据进行metaCycle分析
setwd("WT_FL_0_2dpa")
# data import
library(tidyverse)
library(data.table)

WT_FL_minus2_plus2_day_TMM_sample_exp <- fread("../merged_counts/sample_info_exp.csv") %>% as_tibble()

WT_FL_minus2_plus2_day_genes_TMM_EXPR <- read.table("../merged_counts/genes.TMM.EXPR.matrix")
head(colnames(WT_FL_minus2_plus2_day_genes_TMM_EXPR))
head(WT_FL_minus2_plus2_day_TMM_sample_exp$sample)

# 修正列名 FL-minus2-DPA24h-1 WT.0.DPA12h.1 为 FL_minus2_DPA24h_1
str_replace_all("FL-minus2-DPA24h-1", "-", "_")
WT_FL_minus2_plus2_day_TMM_sample_exp <- WT_FL_minus2_plus2_day_TMM_sample_exp %>%
    mutate(sample = str_replace_all(sample, "-", "_"))

original_colnames <- colnames(WT_FL_minus2_plus2_day_genes_TMM_EXPR)
mutated_colnames <- str_replace_all(original_colnames, "\\.", "_")
colnames(WT_FL_minus2_plus2_day_genes_TMM_EXPR) <- mutated_colnames

WT_minus2_plus2_day_genes_TMM_EXPR <- WT_FL_minus2_plus2_day_genes_TMM_EXPR %>% select(contains("WT"))
head(WT_minus2_plus2_day_genes_TMM_EXPR)
FL_minus2_plus2_day_genes_TMM_EXPR <- WT_FL_minus2_plus2_day_genes_TMM_EXPR %>% select(contains("FL"))

sample_list <- WT_FL_minus2_plus2_day_TMM_sample_exp %>%
    group_by(strain) %>%
    arrange(time) %>%
    select(sample)
head(sample_list)

WT_minus2_plus2_genes_TMM_EXPR_mergeRep <- WT_minus2_plus2_day_genes_TMM_EXPR
FL_minus2_plus2_genes_TMM_EXPR_mergeRep <- FL_minus2_plus2_day_genes_TMM_EXPR

for (i in 1:length(sample_list$strain)) {
    name <- str_remove(sample_list$sample[[i]], "_\\d$")
    if (sample_list$strain[[i]] == "WT") {
        WT_minus2_plus2_genes_TMM_EXPR_mergeRep <- WT_minus2_plus2_genes_TMM_EXPR_mergeRep %>%
            mutate(!!name := rowMeans(select(WT_minus2_plus2_day_genes_TMM_EXPR, contains(name))))
    }
    if (sample_list$strain[[i]] == "FL") {
        FL_minus2_plus2_genes_TMM_EXPR_mergeRep <- FL_minus2_plus2_genes_TMM_EXPR_mergeRep %>%
            mutate(!!name := rowMeans(select(FL_minus2_plus2_day_genes_TMM_EXPR, contains(name))))
    }
}


WT_minus2_plus2_genes_TMM_EXPR_mergeRep_selected <- WT_minus2_plus2_genes_TMM_EXPR_mergeRep[19:length(WT_minus2_plus2_genes_TMM_EXPR_mergeRep)] %>% rownames_to_column(var = "Geneid")
FL_minus2_plus2_genes_TMM_EXPR_mergeRep_selected <- FL_minus2_plus2_genes_TMM_EXPR_mergeRep[19:length(FL_minus2_plus2_genes_TMM_EXPR_mergeRep)] %>% rownames_to_column(var = "Geneid")
# 拆分为3个2周期 每个周期内5个时间点
WT_genes_TMM_EXPR_mergeRep_selected_P1 <- WT_minus2_plus2_genes_TMM_EXPR_mergeRep_selected %>% select(1:6)
WT_genes_TMM_EXPR_mergeRep_selected_P2 <- WT_minus2_plus2_genes_TMM_EXPR_mergeRep_selected %>% select(1, 4:8)
WT_genes_TMM_EXPR_mergeRep_selected_P3 <- WT_minus2_plus2_genes_TMM_EXPR_mergeRep_selected %>% select(1, 6:10)
FL_genes_TMM_EXPR_mergeRep_selected_P1 <- FL_minus2_plus2_genes_TMM_EXPR_mergeRep_selected %>% select(1:6)
FL_genes_TMM_EXPR_mergeRep_selected_P2 <- FL_minus2_plus2_genes_TMM_EXPR_mergeRep_selected %>% select(1, 4:8)
FL_genes_TMM_EXPR_mergeRep_selected_P3 <- FL_minus2_plus2_genes_TMM_EXPR_mergeRep_selected %>% select(1, 6:10)

# 保存至csv文件
write_csv(WT_genes_TMM_EXPR_mergeRep_selected_P1, "mediumDataSave/firstSeqExp/WT_genes_TMM_EXPR_mergeRep_selected_P1.csv")
write_csv(WT_genes_TMM_EXPR_mergeRep_selected_P2, "mediumDataSave/firstSeqExp/WT_genes_TMM_EXPR_mergeRep_selected_P2.csv")
write_csv(WT_genes_TMM_EXPR_mergeRep_selected_P3, "mediumDataSave/firstSeqExp/WT_genes_TMM_EXPR_mergeRep_selected_P3.csv")
write_csv(FL_genes_TMM_EXPR_mergeRep_selected_P1, "mediumDataSave/firstSeqExp/FL_genes_TMM_EXPR_mergeRep_selected_P1.csv")
write_csv(FL_genes_TMM_EXPR_mergeRep_selected_P2, "mediumDataSave/firstSeqExp/FL_genes_TMM_EXPR_mergeRep_selected_P2.csv")
write_csv(FL_genes_TMM_EXPR_mergeRep_selected_P3, "mediumDataSave/firstSeqExp/FL_genes_TMM_EXPR_mergeRep_selected_P3.csv")

write_csv(WT_minus2_plus2_genes_TMM_EXPR_mergeRep_selected, "mediumDataSave/firstSeqExp/WT_minus2_plus2_genes_TMM_EXPR_mergeRep_selected.csv")
write_csv(FL_minus2_plus2_genes_TMM_EXPR_mergeRep_selected, "mediumDataSave/firstSeqExp/FL_minus2_plus2_genes_TMM_EXPR_mergeRep_selected.csv")

library(MetaCycle)
meta2d(
    infile = "mediumDataSave/firstSeqExp/WT_genes_TMM_EXPR_mergeRep_selected_P1.csv",
    filestyle = "csv",
    timepoints = seq(1, 49, 12),
    outdir = "FirstSeqMeta2d",
    parallelize = TRUE,
    nCores = 10,
    minper = 24,
    maxper = 24,
    ARSdefaultPer = 24,
    cycMethod = "JTK",
)

meta2d(
    infile = "mediumDataSave/firstSeqExp/WT_genes_TMM_EXPR_mergeRep_selected_P2.csv",
    filestyle = "csv",
    timepoints = seq(1, 49, 12),
    outdir = "FirstSeqMeta2d",
    parallelize = TRUE,
    nCores = 10,
    minper = 24,
    maxper = 24,
    ARSdefaultPer = 24,
    cycMethod = "JTK",
)

meta2d(
    infile = "mediumDataSave/firstSeqExp/WT_genes_TMM_EXPR_mergeRep_selected_P3.csv",
    filestyle = "csv",
    timepoints = seq(1, 49, 12),
    outdir = "FirstSeqMeta2d",
    parallelize = TRUE,
    nCores = 10,
    minper = 24,
    maxper = 24,
    ARSdefaultPer = 24,
    cycMethod = "JTK",
)

meta2d(
    infile = "mediumDataSave/firstSeqExp/FL_genes_TMM_EXPR_mergeRep_selected_P1.csv",
    filestyle = "csv",
    timepoints = seq(1, 49, 12),
    outdir = "FirstSeqMeta2d",
    parallelize = TRUE,
    nCores = 10,
    minper = 24,
    maxper = 24,
    ARSdefaultPer = 24,
    cycMethod = "JTK",
)

meta2d(
    infile = "mediumDataSave/firstSeqExp/FL_genes_TMM_EXPR_mergeRep_selected_P2.csv",
    filestyle = "csv",
    timepoints = seq(1, 49, 12),
    outdir = "FirstSeqMeta2d",
    parallelize = TRUE,
    nCores = 10,
    minper = 24,
    maxper = 24,
    ARSdefaultPer = 24,
    cycMethod = "JTK",
)

meta2d(
    infile = "mediumDataSave/firstSeqExp/FL_genes_TMM_EXPR_mergeRep_selected_P3.csv",
    filestyle = "csv",
    timepoints = seq(1, 49, 12),
    outdir = "FirstSeqMeta2d",
    parallelize = TRUE,
    nCores = 10,
    minper = 24,
    maxper = 24,
    ARSdefaultPer = 24,
    cycMethod = "JTK",
)

# 全周期的metacycle结果
meta2d(
    infile = "mediumDataSave/firstSeqExp/WT_minus2_plus2_genes_TMM_EXPR_mergeRep_selected.csv",
    filestyle = "csv",
    timepoints = seq(1, 97, 12),
    outdir = "FirstSeqMeta2d",
    parallelize = TRUE,
    nCores = 10,
    minper = 24,
    maxper = 24,
    ARSdefaultPer = 24,
    cycMethod = "JTK",
)

meta2d(
    infile = "mediumDataSave/firstSeqExp/FL_minus2_plus2_genes_TMM_EXPR_mergeRep_selected.csv",
    filestyle = "csv",
    timepoints = seq(1, 97, 12),
    outdir = "FirstSeqMeta2d",
    parallelize = TRUE,
    nCores = 10,
    minper = 24,
    maxper = 24,
    ARSdefaultPer = 24,
    cycMethod = "JTK",
)

## 分周期读取，选取每个周期均ADJ.P > 0.5的基因作为无节律基因
WT_genes_TMM_EXPR_mergeRep_selected_P1_meta_2d <- read_csv("FirstSeqMeta2d/JTKresult_WT_genes_TMM_EXPR_mergeRep_selected_P1.csv")
WT_genes_TMM_EXPR_mergeRep_selected_P1_unCircle_genes <- WT_genes_TMM_EXPR_mergeRep_selected_P1_meta_2d %>%
    filter(ADJ.P > 0.5) %>%
    pull(CycID)

WT_genes_TMM_EXPR_mergeRep_selected_P2_meta_2d <- read_csv("FirstSeqMeta2d/JTKresult_WT_genes_TMM_EXPR_mergeRep_selected_P2.csv")
WT_genes_TMM_EXPR_mergeRep_selected_P2_unCircle_genes <- WT_genes_TMM_EXPR_mergeRep_selected_P2_meta_2d %>%
    filter(ADJ.P > 0.5) %>%
    pull(CycID)

WT_genes_TMM_EXPR_mergeRep_selected_P3_meta_2d <- read_csv("FirstSeqMeta2d/JTKresult_WT_genes_TMM_EXPR_mergeRep_selected_P3.csv")
WT_genes_TMM_EXPR_mergeRep_selected_P3_unCircle_genes <- WT_genes_TMM_EXPR_mergeRep_selected_P3_meta_2d %>%
    filter(ADJ.P > 0.5) %>%
    pull(CycID)

FL_genes_TMM_EXPR_mergeRep_selected_P1_meta_2d <- read_csv("FirstSeqMeta2d/JTKresult_FL_genes_TMM_EXPR_mergeRep_selected_P1.csv")
FL_genes_TMM_EXPR_mergeRep_selected_P1_unCircle_genes <- FL_genes_TMM_EXPR_mergeRep_selected_P1_meta_2d %>%
    filter(ADJ.P > 0.5) %>%
    pull(CycID)

FL_genes_TMM_EXPR_mergeRep_selected_P2_meta_2d <- read_csv("FirstSeqMeta2d/JTKresult_FL_genes_TMM_EXPR_mergeRep_selected_P2.csv")
FL_genes_TMM_EXPR_mergeRep_selected_P2_unCircle_genes <- FL_genes_TMM_EXPR_mergeRep_selected_P2_meta_2d %>%
    filter(ADJ.P > 0.5) %>%
    pull(CycID)

FL_genes_TMM_EXPR_mergeRep_selected_P3_meta_2d <- read_csv("FirstSeqMeta2d/JTKresult_FL_genes_TMM_EXPR_mergeRep_selected_P3.csv")
FL_genes_TMM_EXPR_mergeRep_selected_P3_unCircle_genes <- FL_genes_TMM_EXPR_mergeRep_selected_P3_meta_2d %>%
    filter(ADJ.P > 0.5) %>%
    pull(CycID)

anyUncircaGenes <- reduce(
    list(
        WT_genes_TMM_EXPR_mergeRep_selected_P1_unCircle_genes,
        WT_genes_TMM_EXPR_mergeRep_selected_P2_unCircle_genes,
        WT_genes_TMM_EXPR_mergeRep_selected_P3_unCircle_genes,
        FL_genes_TMM_EXPR_mergeRep_selected_P1_unCircle_genes,
        FL_genes_TMM_EXPR_mergeRep_selected_P2_unCircle_genes,
        FL_genes_TMM_EXPR_mergeRep_selected_P3_unCircle_genes
    ),
    intersect
)
# 将基因列表存入中间文件方便读取
# writeLines(anyUncircaGenes, "mediumDataSave/firstSeqExp/anyUncircaGenes.txt")
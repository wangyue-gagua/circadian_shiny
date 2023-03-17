# 第一次测序结果进行circacompare节律差异比较
library(tidyverse)
library(data.table)
# WT_FL_0_2day_TMM_sample_exp <- read_csv("../merged_counts/WT_FL_0_2day_TMM_sample_exp.csv")
WT_FL_minus2_plus2_day_TMM_sample_exp <- fread("../merged_counts/sample_info_exp.csv") %>% as_tibble()

WT_FL_minus2_plus2_day_genes_TMM_EXPR <- read.table("../merged_counts/genes.TMM.EXPR.matrix")

library(circacompare)
set.seed(42)
# 读取强节律基因列表

WT_FL_minus2_plus2_circa_genes <- readLines("mediumDataSave/firstSeqExp/WT_FL_minus2_plus2_circa_genes.txt")

metaInfo_WT_FL_minus2_plus2_day_TMM_sample_exp <- WT_FL_minus2_plus2_day_TMM_sample_exp %>% select(1:6)

WT_FL_minus2_plus2_circa_genes_compare <- tibble(
    geneId = WT_FL_minus2_plus2_circa_genes,
    messorDiff = numeric(length(WT_FL_minus2_plus2_circa_genes)),
    amplitudeDiff = numeric(length(WT_FL_minus2_plus2_circa_genes)),
    phaseDiff = numeric(length(WT_FL_minus2_plus2_circa_genes)),
    messorDiffPvalue = numeric(length(WT_FL_minus2_plus2_circa_genes)),
    amplitudeDiffPvalue = numeric(length(WT_FL_minus2_plus2_circa_genes)),
    phaseDiffPvalue = numeric(length(WT_FL_minus2_plus2_circa_genes))
)

for (i in seq(1, length(WT_FL_minus2_plus2_circa_genes))) {
    geneId <- WT_FL_minus2_plus2_circa_genes[i]
    df <- cbind(metaInfo_WT_FL_minus2_plus2_day_TMM_sample_exp, WT_FL_minus2_plus2_day_TMM_sample_exp[geneId])
    colnames(df)[7] <- "measure"
    outDf <- circacompare_mixed(
        x = df,
        col_time = "time",
        col_group = "strain",
        col_outcome = "measure",
        col_id = "replicate",
        control = list(grouped_params = c("k", "alpha", "phi"), random_params = c("phi1")),
        period = 24
    )
    if (is.null(outDf)) {
        next
    }
    messorDiff <- outDf$summary %>%
        filter(parameter == "Mesor difference estimate") %>%
        pull(value)
    amplitudeDiff <- outDf$summary %>%
        filter(parameter == "Amplitude difference estimate") %>%
        pull(value)
    phaseDiff <- outDf$summary %>%
        filter(parameter == "Phase difference estimate") %>%
        pull(value)
    messorDiffPvalue <- outDf$summary %>%
        filter(parameter == "P-value for mesor difference") %>%
        pull(value)
    amplitudeDiffPvalue <- outDf$summary %>%
        filter(parameter == "P-value for amplitude difference") %>%
        pull(value)
    phaseDiffPvalue <- outDf$summary %>%
        filter(parameter == "P-value for difference in phase") %>%
        pull(value)
    # 每次输出一行并拼接
    WT_FL_minus2_plus2_circa_genes_compare$messorDiff[[i]] <- messorDiff
    WT_FL_minus2_plus2_circa_genes_compare$amplitudeDiff[[i]] <- amplitudeDiff
    WT_FL_minus2_plus2_circa_genes_compare$phaseDiff[[i]] <- phaseDiff
    WT_FL_minus2_plus2_circa_genes_compare$messorDiffPvalue[[i]] <- messorDiffPvalue
    WT_FL_minus2_plus2_circa_genes_compare$amplitudeDiffPvalue[[i]] <- amplitudeDiffPvalue
    WT_FL_minus2_plus2_circa_genes_compare$phaseDiffPvalue[[i]] <- phaseDiffPvalue
    print(i)
}

# 输出结果
write_tsv(WT_FL_minus2_plus2_circa_genes_compare, "mediumDataSave/firstSeqExp/WT_FL_minus2_plus2_circa_genes_compare.tsv")
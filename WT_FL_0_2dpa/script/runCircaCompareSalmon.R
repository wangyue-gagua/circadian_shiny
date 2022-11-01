library(tidyverse)
salmon_WT_FL_0_2day_TMM_sample_exp <- read_csv("salmon_WT_FL_0_2day_TMM_sample_exp.csv")

salmon_WT_FL_0_2day_genes_TMM_EXPR <- read_delim("salmon_WT_FL_0_2day_genes.TMM.EXPR.matrix",
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
) %>% column_to_rownames(var = "transcriptId")
# cosinor2 Comparison of cosinor parameters of two populations

library(cosinor2)
library(ggplot2)
library(circacompare)
set.seed(42)
# 合并FL和WT的meta2d结果，并通过circacompare添加mesor和amplitude以及phase的差异和p值
salmon_WT_meta2d <- read_csv("salmon_WT_meta2d/JTKresult_salmon_WT_0_2day_genes_TMM_EXPR_mergeRep_selected.csv")
salmon_FL_meta2d <- read_csv("salmon_FL_meta2d/JTKresult_salmon_FL_0_2day_genes_TMM_EXPR_mergeRep_selected.csv")

salmon_FL_WT_meta2d <- left_join(salmon_FL_meta2d, salmon_WT_meta2d, by = "CycID", suffix = c("_FL", "_WT")) %>% 
mutate(geneId = map(CycID, ~str_split(., "\\.")[[1]][1]), hasIsoformDup = duplicated(geneId))

salmon_FL_WT_meta2d_hasIsoform_geneId <- salmon_FL_WT_meta2d %>%
  filter(hasIsoformDup == TRUE) %>% pull(geneId)

salmon_FL_WT_meta2d <- salmon_FL_WT_meta2d %>% mutate(hasIsoform = ifelse(geneId %in% salmon_FL_WT_meta2d_hasIsoform_geneId, TRUE, FALSE))
salmon_FL_WT_meta2d_hasIsoform <- salmon_FL_WT_meta2d %>%
filter(hasIsoform == TRUE)
salmon_FL_WT_meta2d_hasIsoform_all_has_rhythm <- salmon_FL_WT_meta2d_hasIsoform %>% 
filter(BH.Q_WT < 0.01 & BH.Q_FL < 0.01)

metaInfo_salmon_WT_FL_0_2day_TMM_sample_exp <- salmon_WT_FL_0_2day_TMM_sample_exp %>% select(1:6)

# 首先筛选出节律基因，避免模型不收敛
salmon_FL_WT_meta2d_hasIsoform_all_has_rhythm <- salmon_FL_WT_meta2d_hasIsoform_all_has_rhythm %>%
  mutate(messorDiff = c(0), amplitudeDiff = c(0), phaseDiff = c(0), messorDiffPvalue = c(0), amplitudeDiffPvalue = c(0), phaseDiffPvalue = c(0))

salmon_FL_WT_meta2d_hasIsoform_all_has_rhythm_duplicate_gene <- salmon_FL_WT_meta2d_hasIsoform_all_has_rhythm$geneId[duplicated(salmon_FL_WT_meta2d_hasIsoform_all_has_rhythm$geneId)]
salmon_FL_WT_meta2d_hasIsoform_all_has_rhythm_duplicate_gene_tbl <- salmon_FL_WT_meta2d_hasIsoform_all_has_rhythm %>%
  filter(geneId %in% salmon_FL_WT_meta2d_hasIsoform_all_has_rhythm_duplicate_gene)



for (i in seq(1, length(salmon_FL_WT_meta2d_hasIsoform_all_has_rhythm_duplicate_gene_tbl$CycID))) {
  geneId <- salmon_FL_WT_meta2d_hasIsoform_all_has_rhythm_duplicate_gene_tbl$CycID[[i]]
  df <- cbind(metaInfo_salmon_WT_FL_0_2day_TMM_sample_exp, salmon_WT_FL_0_2day_TMM_sample_exp[geneId])
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
  if(is.null(outDf)) {
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
  salmon_FL_WT_meta2d_hasIsoform_all_has_rhythm_duplicate_gene_tbl$messorDiff[[i]] <- messorDiff
  salmon_FL_WT_meta2d_hasIsoform_all_has_rhythm_duplicate_gene_tbl$amplitudeDiff[[i]] <- amplitudeDiff
  salmon_FL_WT_meta2d_hasIsoform_all_has_rhythm_duplicate_gene_tbl$phaseDiff[[i]] <- phaseDiff
  salmon_FL_WT_meta2d_hasIsoform_all_has_rhythm_duplicate_gene_tbl$messorDiffPvalue[[i]] <- messorDiffPvalue
  salmon_FL_WT_meta2d_hasIsoform_all_has_rhythm_duplicate_gene_tbl$amplitudeDiffPvalue[[i]] <- amplitudeDiffPvalue
  salmon_FL_WT_meta2d_hasIsoform_all_has_rhythm_duplicate_gene_tbl$phaseDiffPvalue[[i]] <- phaseDiffPvalue
  print(i)
}

write_csv(salmon_FL_WT_meta2d_hasIsoform_all_has_rhythm_duplicate_gene_tbl, "salmon_FL_WT_meta2d_hasIsoform_all_has_rhythm_duplicate_gene_tbl.csv")
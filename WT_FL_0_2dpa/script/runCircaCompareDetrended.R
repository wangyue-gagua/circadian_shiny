# 对detrend后数据进行circacompare分析

# data import
library(tidyverse)

WT_FL_0_2day_TMM_sample_exp_detrended <- data.table::fread("./mediumDataSave/WT_FL_0_2day_TMM_sample_exp_detrended.csv") %>%
    as_tibble()
head(WT_FL_0_2day_TMM_sample_exp_detrended)
library(circacompare)
set.seed(42)
# 合并FL和WT的meta2d结果，并通过circacompare添加mesor和amplitude以及phase的差异和p值
WT_meta2d <- read_csv("WT_meta2d_detrended/JTKresult_WT_0_2day_genes_TMM_EXPR_mergeRep_selected_detrended.csv")
FL_meta2d <- read_csv("FL_meta2d_detrended/JTKresult_FL_0_2day_genes_TMM_EXPR_mergeRep_selected_detrended.csv")

FL_WT_meta2d <- left_join(FL_meta2d, WT_meta2d, by = "CycID", suffix = c("_FL", "_WT"))

metaInfo_WT_FL_0_2day_TMM_sample_exp_detrended <- WT_FL_0_2day_TMM_sample_exp_detrended %>% select(1:6)

# 首先筛选出节律基因，避免模型不收敛
FL_WT_meta2d_circa <- FL_WT_meta2d %>%
  filter(BH.Q_FL < 0.01 & BH.Q_WT < 0.01) %>%
  mutate(messorDiff = c(0), amplitudeDiff = c(0), phaseDiff = c(0), messorDiffPvalue = c(0), amplitudeDiffPvalue = c(0), phaseDiffPvalue = c(0))

for (i in seq(1, length(FL_WT_meta2d_circa$CycID))) {
  geneId <- FL_WT_meta2d_circa$CycID[[i]]
  df <- cbind(metaInfo_WT_FL_0_2day_TMM_sample_exp_detrended, WT_FL_0_2day_TMM_sample_exp_detrended[geneId])
  colnames(df)[7] <- "measure"
  outDf <- circacompare_mixed(
    x = df,
    col_time = "time",
    col_group = "strain",
    col_outcome = "measure",
    col_id = "replicate",
    control = list(grouped_params = c("alpha", "phi"), random_params = c("phi1")), # k即mesor，因为detrend后mesor恒为1，故会导致模型不收敛，故不考虑k
    period = 24
  )
  if(is.null(outDf)) {
    next
  }
  amplitudeDiff <- outDf$summary %>%
    filter(parameter == "Amplitude difference estimate") %>%
    pull(value)
  phaseDiff <- outDf$summary %>%
    filter(parameter == "Phase difference estimate") %>%
    pull(value)
  amplitudeDiffPvalue <- outDf$summary %>%
    filter(parameter == "P-value for amplitude difference") %>%
    pull(value)
  phaseDiffPvalue <- outDf$summary %>%
    filter(parameter == "P-value for difference in phase") %>%
    pull(value)
  FL_WT_meta2d_circa$amplitudeDiff[[i]] <- amplitudeDiff
  FL_WT_meta2d_circa$phaseDiff[[i]] <- phaseDiff
  FL_WT_meta2d_circa$amplitudeDiffPvalue[[i]] <- amplitudeDiffPvalue
  FL_WT_meta2d_circa$phaseDiffPvalue[[i]] <- phaseDiffPvalue
  print(i)
}

write_tsv(FL_WT_meta2d_circa, "mediumDataSave/detrended/WT_FL_0_2day_TMM_sample_exp_detrended_circacompare.tsv")
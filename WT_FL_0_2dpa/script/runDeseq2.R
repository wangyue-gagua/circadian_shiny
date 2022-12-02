## DEseq2差异表达分析
library(DESeq2)
library(tidyverse)

WT_FL_0_2day_count_matrix <- read_tsv("../merged_counts/WT_FL_0_2day_genes.counts.matrix") %>% column_to_rownames(var = "...1")
WT_FL_0_2day_TMM_sample_exp <- data.table::fread("../merged_counts/WT_FL_0_2day_TMM_sample_exp.csv") %>%
    as_tibble()
sampleInfo <- WT_FL_0_2day_TMM_sample_exp %>%
    select(c(1:6)) %>%
    column_to_rownames(var = "sample") %>%
    mutate(phase = time %% 24) %>%
    mutate_all(factor)

# dds <- DESeqDataSetFromMatrix(
# )
# deseqDDS <- DESeq(dds)
# res <- results(deseqDDS, contrast = c("time
#     countData = WT_FL_0_2day_count_matrix,
#     colData = sampleInfo,
#     design = ~ time + strain", "25", "69"))
# head(res)
# resultsNames(deseqDDS)
# WT_FL_0_2day_WTvsFL_res <- results(deseqDDS, contrast = c("strain", "WT", "FL"))
# WT_FL_0_2day_WTvsFL_res_tible <- WT_FL_0_2day_WTvsFL_res %>%
#     as.data.frame() %>%
#     rownames_to_column(var = "Geneid") %>%
#     drop_na()

# head(WT_FL_0_2day_WTvsFL_res_tible)
# write_tsv(WT_FL_0_2day_WTvsFL_res_tible, "mediumDataSave/deseq2/WT_FL_0_2day_WTvsFL.tsv")


## 对每个时间点分别进行差异表达分析WT与FL差异

timePointWTFLDeseq <- function(timePoint) {
    timePointSampleInfo <- sampleInfo %>% filter(time == timePoint)
    timePointCountMatrix <- WT_FL_0_2day_count_matrix %>%
        select(rownames(timePointSampleInfo)) %>%
        filter(rowSums(.) > 0)
    timePointDeseq <- DESeqDataSetFromMatrix(
        countData = timePointCountMatrix,
        colData = timePointSampleInfo,
        design = ~strain
    )
    timePointDeseq <- DESeq(timePointDeseq)
    print("DE　完了")
    timePointRes <- results(timePointDeseq, contrast = c("strain", "WT", "FL"))
    timePointResTable <- timePointRes %>%
        as.data.frame() %>%
        rownames_to_column(var = "Geneid") %>%
        drop_na()
    write_tsv(timePointResTable, paste0("mediumDataSave/deseq2/WT_FL_0_2day_", timePoint, "h_WTvsFL.tsv"))
}

for (i in unique(sort(sampleInfo$time))) {
    timePointWTFLDeseq(i)
}

## 对每个phase分别进行差异表达分析WT与FL差异
phaseWTFLDeseq <- function(phasePoint) {
    phaseSampleInfo <- sampleInfo %>% filter(phase == phasePoint)
    phaseCountMatrix <- WT_FL_0_2day_count_matrix %>%
        select(rownames(phaseSampleInfo)) %>%
        filter(rowSums(.) > 0)
    phaseDeseq <- DESeqDataSetFromMatrix(
        countData = phaseCountMatrix,
        colData = phaseSampleInfo,
        design = ~strain
    )
    phaseDeseq <- DESeq(phaseDeseq)
    print("DE　完了")
    phaseRes <- results(phaseDeseq, contrast = c("strain", "WT", "FL"))
    phaseResTable <- phaseRes %>%
        as.data.frame() %>%
        rownames_to_column(var = "Geneid") %>%
        drop_na()
    write_tsv(phaseResTable, paste0("mediumDataSave/deseq2/WT_FL_0_2day_", phasePoint, "Phase_WTvsFL.tsv"))
}

for (i in unique(sort(sampleInfo$phase))) {
    phaseWTFLDeseq(i)
}

## 对每个period分别进行差异表达分析WT与FL差异
periodPointWTFLDeseq <- function(periodPoint) {
    periodSampleInfo <- sampleInfo %>% filter(period == periodPoint)
    periodCountMatrix <- WT_FL_0_2day_count_matrix %>%
        select(rownames(periodSampleInfo)) %>%
        filter(rowSums(.) > 0)
    periodDeseq <- DESeqDataSetFromMatrix(
        countData = periodCountMatrix,
        colData = periodSampleInfo,
        design = ~strain
    )
    periodDeseq <- DESeq(periodDeseq)
    print("DE　完了")
    periodRes <- results(periodDeseq, contrast = c("strain", "WT", "FL"))
    periodResTable <- periodRes %>%
        as.data.frame() %>%
        rownames_to_column(var = "Geneid") %>%
        drop_na()
    write_tsv(periodResTable, paste0("mediumDataSave/deseq2/WT_FL_0_2day_", periodPoint, "Period_WTvsFL.tsv"))
}

for (i in unique(sort(sampleInfo$period))) {
    periodPointWTFLDeseq(i)
}

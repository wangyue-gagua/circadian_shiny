# 第一次测序与第二次测序检测基因数对比
first_sample_info_exp <- data.table::fread("../merged_counts/sample_info_exp.csv") %>%
    as_tibble()
first_sample_list <- first_sample_info_exp %>%
    group_by(strain) %>%
    arrange(time, .by_group = T) %>%
    pull(sample)
first_count_matrix <- read_tsv("../merged_counts/genes.counts.matrix") %>%
    rename(Geneid = `...1`) %>%
    select(Geneid, !!first_sample_list)

second_sample_info_exp <- data.table::fread("../merged_counts/WT_FL_0_2day_TMM_sample_exp.csv") %>%
    as_tibble()
second_sample_list <- second_sample_info_exp %>%
    group_by(strain) %>%
    arrange(time, .by_group = T) %>%
    pull(sample)


second_count_matrix <- read_tsv("../merged_counts/WT_FL_0_2day_genes.counts.matrix") %>%
    rename(Geneid = `...1`) %>%
    select(Geneid, !!second_sample_list)


first_count_matrix_longer <- first_count_matrix %>%
    pivot_longer(cols = -Geneid, names_to = "sample", values_to = "count") %>%
    mutate(sample = factor(sample, levels = first_sample_list))

second_count_matrix_longer <- second_count_matrix %>%
    pivot_longer(cols = -Geneid, names_to = "sample", values_to = "count") %>%
    mutate(sample = factor(sample, levels = second_sample_list))

stat_box_data <- function(y, upper_limit = log10(max(first_count_matrix_longer$count)) * 1.15) {
    return(
        data.frame(
            y = 0.95 * upper_limit,
            label = length(y[10^y > 3])
        )
    )
}
first_count_matrix_longer_pdf <- first_count_matrix_longer %>% ggplot(aes(x = sample, y = log10(count))) +
    geom_boxplot() +
    stat_summary(
        fun.data = stat_box_data,
        geom = "text",
        size = 2,
        color = "red",
        hjust = 0.5,
        vjust = 0.9
    ) +
    labs(title = "First sequencing Gene Counts") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("figure/geneExpDiff/first_count_matrix_longer.pdf", first_count_matrix_longer_pdf, width = 10, height = 4)

second_stat_box_data <- function(y, upper_limit = log10(max(second_count_matrix_longer$count)) * 1.15) {
    return(
        data.frame(
            y = 0.95 * upper_limit,
            label = length(y[10^y > 3])
        )
    )
}
second_count_matrix_longer_pdf <- second_count_matrix_longer %>% ggplot(aes(x = sample, y = log10(count))) +
    geom_boxplot() +
    stat_summary(
        fun.data = second_stat_box_data,
        geom = "text",
        size = 2,
        color = "red",
        hjust = 0.5,
        vjust = 0.9
    ) +
    labs(title = "Second sequencing Gene Counts") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("figure/geneExpDiff/second_count_matrix_longer.pdf", second_count_matrix_longer_pdf, width = 20, height = 4)

library(ggplot2)
library(DESeq2)
library(ggrepel)
library(patchwork)
# 差异分析火山图
padjCutOff <- 0.0001
## period 0，1，2
period0_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_0Period_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Period0Plot <- deResVolcano(period0_de, "Period 0 WT vs FL")

period1_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_1Period_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Period1Plot <- deResVolcano(period1_de, "Period 1 WT vs FL")

period2_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_2Period_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Period2Plot <- deResVolcano(period2_de, "Period 2 WT vs FL")

ggsave("figure/deseq2/Period0_1_2_WTvsFL.pdf", Period0Plot + Period1Plot + Period2Plot + plot_layout(guides = "collect"), width = 20, height = 6)

## phase 1，5，9，13，17，21
phase1_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_1Phase_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Phase1Plot <- deResVolcano(phase1_de, "Phase 1 WT vs FL")

phase5_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_5Phase_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Phase5Plot <- deResVolcano(phase5_de, "Phase 5 WT vs FL")

phase9_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_9Phase_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Phase9Plot <- deResVolcano(phase9_de, "Phase 9 WT vs FL")

phase13_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_13Phase_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Phase13Plot <- deResVolcano(phase13_de, "Phase 13 WT vs FL")

phase17_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_17Phase_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Phase17Plot <- deResVolcano(phase17_de, "Phase 17 WT vs FL")

phase21_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_21Phase_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Phase21Plot <- deResVolcano(phase21_de, "Phase 21 WT vs FL")

ggsave("figure/deseq2/Phase1_5_9_13_17_21_WTvsFL.pdf",
    Phase1Plot + Phase5Plot + Phase9Plot + Phase13Plot + Phase17Plot + Phase21Plot + plot_layout(ncol = 3, guides = "collect"),
    width = 20, height = 12
)

## timePoint seq(1, 69, 4)
time1_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_1h_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Time1Plot <- deResVolcano(time1_de, "Time 1 WT vs FL")

time5_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_5h_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Time5Plot <- deResVolcano(time5_de, "Time 5 WT vs FL")

time9_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_9h_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Time9Plot <- deResVolcano(time9_de, "Time 9 WT vs FL")

time13_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_13h_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Time13Plot <- deResVolcano(time13_de, "Time 13 WT vs FL")

time17_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_17h_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Time17Plot <- deResVolcano(time17_de, "Time 17 WT vs FL")

time21_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_21h_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Time21Plot <- deResVolcano(time21_de, "Time 21 WT vs FL")

ggsave("figure/deseq2/Time1_5_9_13_17_21_WTvsFL.pdf",
    Time1Plot + Time5Plot + Time9Plot + Time13Plot + Time17Plot + Time21Plot + plot_layout(ncol = 3, guides = "collect"),
    width = 20, height = 12
)

time25_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_25h_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Time25Plot <- deResVolcano(time25_de, "Time 25 WT vs FL")

time29_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_29h_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Time29Plot <- deResVolcano(time29_de, "Time 29 WT vs FL")

time33_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_33h_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Time33Plot <- deResVolcano(time33_de, "Time 33 WT vs FL")

time37_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_37h_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Time37Plot <- deResVolcano(time37_de, "Time 37 WT vs FL")

time41_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_41h_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Time41Plot <- deResVolcano(time41_de, "Time 41 WT vs FL")

time45_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_45h_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Time45Plot <- deResVolcano(time45_de, "Time 45 WT vs FL")

ggsave("figure/deseq2/Time25_29_33_37_41_45_WTvsFL.pdf",
    Time25Plot + Time29Plot + Time33Plot + Time37Plot + Time41Plot + Time45Plot + plot_layout(ncol = 3, guides = "collect"),
    width = 20, height = 12
)

time49_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_49h_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Time49Plot <- deResVolcano(time49_de, "Time 49 WT vs FL")

time53_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_53h_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Time53Plot <- deResVolcano(time53_de, "Time 53 WT vs FL")

time57_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_57h_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Time57Plot <- deResVolcano(time57_de, "Time 57 WT vs FL")

time61_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_61h_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Time61Plot <- deResVolcano(time61_de, "Time 61 WT vs FL")

time65_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_65h_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Time65Plot <- deResVolcano(time65_de, "Time 65 WT vs FL")

time69_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_69h_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
Time69Plot <- deResVolcano(time69_de, "Time 69 WT vs FL")

ggsave("figure/deseq2/Time49_53_57_61_65_69_WTvsFL.pdf",
    Time49Plot + Time53Plot + Time57Plot + Time61Plot + Time65Plot + Time69Plot + plot_layout(ncol = 3, guides = "collect"),
    width = 20, height = 12
)

## all
all_de <- read_tsv("mediumDataSave/deseq2/WT_FL_0_2day_WTvsFL.tsv") %>%
    mutate(change = ifelse(padj > padjCutOff | abs(log2FoldChange) < 1, "ns",
        ifelse(log2FoldChange > 1, "up", "down")
    ))
AllPlot <- deResVolcano(all_de, "All WT vs FL")
ggsave("figure/deseq2/All_WTvsFL.pdf", AllPlot, width = 10, height = 10)

all_de_WT_high_geneList <- all_de %>%
    filter(change == "up") %>%
    pull(Geneid)
writeLines(all_de_WT_high_geneList, "mediumDataSave/deseq2/all_de_WT_high_geneList.txt")

all_de_WT_low_geneList <- all_de %>%
    filter(change == "down") %>%
    pull(Geneid)
writeLines(all_de_WT_low_geneList, "mediumDataSave/deseq2/all_de_WT_low_geneList.txt")

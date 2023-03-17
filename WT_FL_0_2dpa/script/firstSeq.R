library(tidyverse)
library(data.table)
# WT_FL_0_2day_TMM_sample_exp <- read_csv("../merged_counts/WT_FL_0_2day_TMM_sample_exp.csv")
WT_FL_minus2_plus2_day_TMM_sample_exp <- fread("../merged_counts/sample_info_exp.csv") %>% as_tibble()

WT_FL_minus2_plus2_day_genes_TMM_EXPR <- read.table("../merged_counts/genes.TMM.EXPR.matrix")
head(WT_FL_minus2_plus2_day_genes_TMM_EXPR)
WT_FL_minus2_plus2_day_TMM_sample_exp$time


WT_FL_minus2_plus2_day_auxinRelatedGeneScaledMtx <- WT_FL_minus2_plus2_day_TMM_sample_exp %>%
    dplyr::select(c(1:6, {{ auxinRelatedGeneList }})) %>%
    mutate_if(is.double, scale) %>%
    pivot_longer(cols = starts_with("Ghir"), names_to = "genes", values_to = "expression") %>%
    left_join(auxinRelatedTibble, by = c("genes" = "geneId"))

WT_FL_minus2_plus2_day_auxinRelatedGeneScaledPlot <- WT_FL_minus2_plus2_day_auxinRelatedGeneScaledMtx %>%
    ggplot(aes(time, expression)) +
    geom_point(aes(col = genes)) +
    geom_smooth(aes(group = interaction(as.factor(replicate), strain, genes), color = genes, linetype = strain), span = 0.3) +
    scale_x_continuous(breaks = seq(-48, 48, 12)) +
    # facet_wrap(~replicate, nrow = 2) +
    # geom_rect(
    #     data = data.frame(
    #                 xstart = c(14.5, 38.5, 62.5),
    #                 xend = c(23, 47, 71)
    #             ),
    #     aes(
    #         xmin = xstart,
    #         xmax = xend,
    #         ymin = -Inf,
    #         ymax = Inf
    #     ),
    #     inherit.aes = FALSE,
    #     alpha = 0.2
    # ) +
    # geom_vline(xintercept = c(6, 18), linetype = "dotted") +
    labs(title = "Phase diff between Auxin related cascading regulation", subtitle = NULL) +
    ylab("Z-score of relative expression/(TMM)") +
    facet_grid(
        Families ~ replicate,
        scales = "free_y"
    ) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave("figure/geneRepCirca/WT_FL_minus2_plus2_day_auxinRelatedGeneScaledPlot.pdf", WT_FL_minus2_plus2_day_auxinRelatedGeneScaledPlot, width = 10, height = 10)

blastPResult <- read_tsv("../bla_go/blastp_AD1_HAU_v1.0_vs_arabidopsis.1.txt")
head(blastPResult)
# Query  Ghir_A13G001130.1_HAU-AD1_v1.0
cytokinRelatedGenes <- blastPResult %>%
    filter(Description %>% str_detect("cytokinin")) %>%
    mutate(Geneid = str_extract(Query, "Ghir_[A-Z0-9]+")) %>%
    select(Geneid, Description)
cytokinRelatedGenesList <- cytokinRelatedGenes$Geneid
WT_FL_minus2_plus2_day_cytokinRelatedGeneScaledMtx <- WT_FL_minus2_plus2_day_TMM_sample_exp %>%
    dplyr::select(c(1:6, {{ cytokinRelatedGenesList }})) %>%
    mutate_if(is.double, scale) %>%
    pivot_longer(cols = starts_with("Ghir"), names_to = "genes", values_to = "expression") %>%
    left_join(cytokinRelatedGenes, by = c("genes" = "Geneid"))
head(WT_FL_minus2_plus2_day_cytokinRelatedGeneScaledMtx$labs)

WT_FL_minus2_plus2_day_cytokinRelatedGeneScaledPlot <- WT_FL_minus2_plus2_day_cytokinRelatedGeneScaledMtx %>%
    ggplot(aes(time, expression)) +
    geom_point(aes(col = genes)) +
    geom_smooth(aes(group = interaction(as.factor(replicate), strain, genes), color = genes, linetype = strain), span = 0.3) +
    scale_x_continuous(breaks = seq(-48, 48, 12)) +
    # facet_wrap(~replicate, nrow = 2) +
    geom_rect(
        data = data.frame(
            xstart = seq(-36, 36, 24),
            xend = seq(-24, 48, 24)
        ),
        aes(
            xmin = xstart,
            xmax = xend,
            ymin = -Inf,
            ymax = Inf
        ),
        inherit.aes = FALSE,
        alpha = 0.2
    ) +
    # geom_vline(xintercept = c(6, 18), linetype = "dotted") +
    labs(title = "Phase diff between Cytokinin related cascading regulation", subtitle = NULL) +
    ylab("Z-score of relative expression/(TMM)") +
    facet_grid(
        Description ~ replicate,
        scales = "free_y"
    ) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave("figure/geneRepCirca/WT_FL_minus2_plus2_day_cytokinRelatedGeneScaledPlot.pdf", WT_FL_minus2_plus2_day_cytokinRelatedGeneScaledPlot, width = 10, height = 10)

## ABA
abaRelatedGenes <- blastPResult %>%
    filter(Description %>% str_detect("ABA")) %>%
    mutate(Geneid = str_extract(Query, "Ghir_[A|D][A-Z0-9]+")) %>%
    select(Geneid, Description)
abaRelatedGenesList <- abaRelatedGenes$Geneid %>% na.omit()
WT_FL_minus2_plus2_day_abaRelatedGeneScaledMtx <- WT_FL_minus2_plus2_day_TMM_sample_exp %>%
    dplyr::select(c(1:6, {{ abaRelatedGenesList }})) %>%
    mutate_if(is.double, scale) %>%
    pivot_longer(cols = starts_with("Ghir"), names_to = "genes", values_to = "expression") %>%
    left_join(abaRelatedGenes, by = c("genes" = "Geneid"))
WT_FL_minus2_plus2_day_abaRelatedGeneScaledPlot <- WT_FL_minus2_plus2_day_abaRelatedGeneScaledMtx %>%
    ggplot(aes(time, expression)) +
    geom_point(aes(col = genes)) +
    geom_smooth(aes(group = interaction(as.factor(replicate), strain, genes), color = genes, linetype = strain), span = 0.3) +
    scale_x_continuous(breaks = seq(-48, 48, 12)) +
    # facet_wrap(~replicate, nrow = 2) +
    geom_rect(
        data = data.frame(
            xstart = seq(-36, 36, 24),
            xend = seq(-24, 48, 24)
        ),
        aes(
            xmin = xstart,
            xmax = xend,
            ymin = -Inf,
            ymax = Inf
        ),
        inherit.aes = FALSE,
        alpha = 0.2
    ) +
    # geom_vline(xintercept = c(6, 18), linetype = "dotted") +
    labs(title = "Phase diff between ABA related cascading regulation", subtitle = NULL) +
    ylab("Z-score of relative expression/(TMM)") +
    facet_grid(
        Description ~ replicate,
        scales = "free_y"
    ) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave("figure/geneRepCirca/WT_FL_minus2_plus2_day_abaRelatedGeneScaledPlot.pdf", WT_FL_minus2_plus2_day_abaRelatedGeneScaledPlot, width = 10, height = 10)

## jasmonate
jasmonateRelatedGenes <- blastPResult %>%
    filter(Description %>% str_detect("jasmonate")) %>%
    mutate(Geneid = str_extract(Query, "Ghir_[A|D][A-Z0-9]+")) %>%
    select(Geneid, Description)
jasmonateRelatedGenesList <- jasmonateRelatedGenes$Geneid %>% na.omit()
WT_FL_minus2_plus2_day_jasmonateRelatedGeneScaledMtx <- WT_FL_minus2_plus2_day_TMM_sample_exp %>%
    dplyr::select(c(1:6, {{ jasmonateRelatedGenesList }})) %>%
    mutate_if(is.double, scale) %>%
    pivot_longer(cols = starts_with("Ghir"), names_to = "genes", values_to = "expression") %>%
    left_join(jasmonateRelatedGenes, by = c("genes" = "Geneid"))
WT_FL_minus2_plus2_day_jasmonateRelatedGeneScaledPlot <- WT_FL_minus2_plus2_day_jasmonateRelatedGeneScaledMtx %>%
    ggplot(aes(time, expression)) +
    geom_point(aes(col = genes)) +
    geom_smooth(aes(group = interaction(as.factor(replicate), strain, genes), color = genes, linetype = strain), span = 0.3) +
    scale_x_continuous(breaks = seq(-48, 48, 12)) +
    # facet_wrap(~replicate, nrow = 2) +
    geom_rect(
        data = data.frame(
            xstart = seq(-36, 36, 24),
            xend = seq(-24, 48, 24)
        ),
        aes(
            xmin = xstart,
            xmax = xend,
            ymin = -Inf,
            ymax = Inf
        ),
        inherit.aes = FALSE,
        alpha = 0.2
    ) +
    # geom_vline(xintercept = c(6, 18), linetype = "dotted") +
    labs(title = "Phase diff between jasmonate related cascading regulation", subtitle = NULL) +
    ylab("Z-score of relative expression/(TMM)") +
    facet_grid(
        Description ~ replicate,
        scales = "free_y"
    ) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave("figure/geneRepCirca/WT_FL_minus2_plus2_day_jasmonateRelatedGeneScaledPlot.pdf", WT_FL_minus2_plus2_day_jasmonateRelatedGeneScaledPlot, width = 10, height = 10)

## gibberellin
gibberellinRelatedGenes <- blastPResult %>%
    filter(Description %>% str_detect("gibberellin")) %>%
    mutate(Geneid = str_extract(Query, "Ghir_[A|D][A-Z0-9]+")) %>%
    select(Geneid, Description)
gibberellinRelatedGenesList <- gibberellinRelatedGenes$Geneid %>% na.omit()
WT_FL_minus2_plus2_day_gibberellinRelatedGeneScaledMtx <- WT_FL_minus2_plus2_day_TMM_sample_exp %>%
    dplyr::select(c(1:6, {{ gibberellinRelatedGenesList }})) %>%
    mutate_if(is.double, scale) %>%
    pivot_longer(cols = starts_with("Ghir"), names_to = "genes", values_to = "expression") %>%
    left_join(gibberellinRelatedGenes, by = c("genes" = "Geneid"))
WT_FL_minus2_plus2_day_gibberellinRelatedGeneScaledPlot <- WT_FL_minus2_plus2_day_gibberellinRelatedGeneScaledMtx %>%
    ggplot(aes(time, expression)) +
    geom_point(aes(col = genes)) +
    geom_smooth(aes(group = interaction(as.factor(replicate), strain, genes), color = genes, linetype = strain), span = 0.3) +
    scale_x_continuous(breaks = seq(-48, 48, 12)) +
    # facet_wrap(~replicate, nrow = 2) +
    geom_rect(
        data = data.frame(
            xstart = seq(-36, 36, 24),
            xend = seq(-24, 48, 24)
        ),
        aes(
            xmin = xstart,
            xmax = xend,
            ymin = -Inf,
            ymax = Inf
        ),
        inherit.aes = FALSE,
        alpha = 0.2
    ) +
    # geom_vline(xintercept = c(6, 18), linetype = "dotted") +
    labs(title = "Phase diff between gibberellin related cascading regulation", subtitle = NULL) +
    ylab("Z-score of relative expression/(TMM)") +
    facet_grid(
        Description ~ replicate,
        scales = "free_y"
    ) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave("figure/geneRepCirca/WT_FL_minus2_plus2_day_gibberellinRelatedGeneScaledPlot.pdf", WT_FL_minus2_plus2_day_gibberellinRelatedGeneScaledPlot, width = 10, height = 10)

# 用计算WT FL spearman相关系数的方法寻找差异基因
# 还是仅采取强节律基因
WT_FL_minus2_plus2_circa_genes <- readLines("mediumDataSave/firstSeqExp/WT_FL_minus2_plus2_circa_genes.txt")

WT_FL_minus2_plus2_circa_genes_TMM_sample_exp <- WT_FL_minus2_plus2_day_TMM_sample_exp %>%
    dplyr::select(c(1:6, {{ WT_FL_minus2_plus2_circa_genes }}))

WT_FL_minus2_plus2_circa_genes_TMM_EXPR <- WT_FL_minus2_plus2_day_genes_TMM_EXPR[WT_FL_minus2_plus2_circa_genes, ]

WT_FL_minus2_plus2_circa_genes_TMM_EXPR_cor <- WT_FL_minus2_plus2_circa_genes_TMM_EXPR %>% 
    rownames_to_column("geneId") %>% 
    rowwise() %>% 
    mutate(spearman_cor = cor.test(
        c_across(contains("WT")),
        c_across(contains("FL")),
        method = "spearman"
    )$estimate)

## 保存相关系数临时文件
# write_csv(WT_FL_minus2_plus2_circa_genes_TMM_EXPR_cor, "mediumDataSave/firstSeqExp/WT_FL_minus2_plus2_circa_genes_TMM_EXPR_cor.csv")
## 提取WT_FL_minus2_plus2_circa_genes_TMM_EXPR_cor$spearman_cor的四分位
WT_FL_minus2_plus2_circa_genes_TMM_EXPR_cor_quantile <- quantile(WT_FL_minus2_plus2_circa_genes_TMM_EXPR_cor$spearman_cor, probs = c(0.25, 0.75))
WT_FL_minus2_plus2_circa_genes_TMM_EXPR_cor_quantile

WT_FL_minus2_plus2_circa_genes_TMM_EXPR_cor %>% filter(
    spearman_cor == max(WT_FL_minus2_plus2_circa_genes_TMM_EXPR_cor$spearman_cor)
) %>% pull(geneId)

WT_FL_minus2_plus2_circa_genes_TMM_EXPR_cor %>% filter(
    spearman_cor == min(WT_FL_minus2_plus2_circa_genes_TMM_EXPR_cor$spearman_cor)
) %>% pull(geneId)

WT_FL_minus2_plus2_circa_genes_TMM_EXPR_rev_genes <- WT_FL_minus2_plus2_circa_genes_TMM_EXPR_cor %>% 
    filter(spearman_cor < 0) %>% 
    pull(geneId)

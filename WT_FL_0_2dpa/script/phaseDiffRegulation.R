# 探究相位差异与基因表达调控之间的关系
## 首先， 探究已知级联调控基因的相位差异
### Ghir_A08G011910 WD40
### Ghir_A08G017600 Ghir_D08G018460 GhMYB2
### Ghir_D05G013110 Ghir_A05G013360 GhGL2
### Ghir_D01G012940 Ghir_D12G017660 MYB25-like
library(readxl)
library(ggplot2)

myb25CascadeGenes <- read_excel("mediumDataSave/cascadeGenes/Myb25Cascade.xlsx", col_names = F)
colnames(myb25CascadeGenes) <- c("gene", "description")
myb25CascadeGenes_df <- separate(myb25CascadeGenes, description, c("GhirAlias", "AtAlias"), sep = "/") %>%
    separate_rows(gene, sep = "/")
GhirAliasLevelTbl <- tibble(
    GhirAlias = c("WD40", "GhMYB2", "GhGL2", "GhMYB25-like", "GhMYB25", "GhHD-1", "GhHOX3-RT-F", "GhEXPASIN_A1", "GhLGO"),
    cascadeLevel = c("L1", "L1", "L1", "L2", "L3", "L3", "L3", "L4", "L4")
)
myb25CascadeGenes_df <- myb25CascadeGenes_df %>%
    left_join(GhirAliasLevelTbl, by = "GhirAlias")
head(myb25CascadeGenes_df)
myb25CascadeGeneList <- myb25CascadeGenes_df$gene
myb25CascadeGenesMtx <- WT_FL_0_2day_TMM_sample_exp %>%
    dplyr::select(c(1:6, {{ myb25CascadeGeneList }})) %>%
    pivot_longer(cols = starts_with("Ghir"), names_to = "genes", values_to = "expression") %>%
    left_join(myb25CascadeGenes_df, by = c("genes" = "gene"))
head(myb25CascadeGenesMtx)
rects <-
    data.frame(
        xstart = c(14.5, 38.5, 62.5),
        xend = c(23, 47, 71)
    )
WD40_ExpA1_cascading_regulation <- myb25CascadeGenesMtx %>%
    ggplot(aes(time, expression)) +
    geom_point(aes(col = genes)) +
    geom_smooth(aes(group = interaction(as.factor(replicate), strain, genes), color = genes, linetype = strain), span = 0.3) +
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
    # geom_vline(xintercept = c(6, 18), linetype = "dotted") +
    labs(title = "Phase diff between WD40-EXPA1 cascading regulation", subtitle = NULL) +
    ylab("relative expression/(TMM)") +
    facet_grid(
        cascadeLevel + factor(GhirAlias,
            levels = c("WD40", "GhMYB2", "GhGL2", "GhMYB25-like", "GhMYB25", "GhHD-1", "GhHOX3-RT-F", "GhEXPASIN_A1", "GhLGO")
        ) ~ replicate,
        scales = "free_y"
    ) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave("figure/geneRepCirca/WD40_ExpA1_cascading_regulation.pdf", WD40_ExpA1_cascading_regulation, width = 10, height = 10)

## 生长素相关基因
auxinRelatedTibble <- tibble(
    geneId = c(
        "Ghir_A10G002560", "Ghir_D10G003340",
        "Ghir_A05G002590", "Ghir_A06G020380", "Ghir_A09G024180", "Ghir_D05G002730", "Ghir_D05G004880", "Ghir_D09G023320",
        "Ghir_A01G008830", "Ghir_A05G040190",
        "Ghir_A03G022510", "Ghir_D13G010120"
    ),
    Families = c(
        rep("ARF", 2),
        rep("Aux/IAA", 6),
        rep("GH3", 2),
        rep("SAUR", 2)
    )
)
auxinRelatedGeneList <- auxinRelatedTibble$geneId
auxinRelatedGeneMtx <- WT_FL_0_2day_TMM_sample_exp %>%
    dplyr::select(c(1:6, {{ auxinRelatedGeneList }})) %>%
    pivot_longer(cols = starts_with("Ghir"), names_to = "genes", values_to = "expression") %>%
    left_join(auxinRelatedTibble, by = c("genes" = "geneId"))

auxinRelatedPlot <- auxinRelatedGeneMtx %>%
    ggplot(aes(time, expression)) +
    geom_point(aes(col = genes)) +
    geom_smooth(aes(group = interaction(as.factor(replicate), strain, genes), color = genes, linetype = strain), span = 0.3) +
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
    # geom_vline(xintercept = c(6, 18), linetype = "dotted") +
    labs(title = "Phase diff between Auxin related cascading regulation", subtitle = NULL) +
    ylab("relative expression/(TMM)") +
    facet_grid(
        Families + factor(genes,
            levels = c(
                "Ghir_A10G002560", "Ghir_D10G003340",
                "Ghir_A05G002590", "Ghir_A06G020380", "Ghir_A09G024180", "Ghir_D05G002730", "Ghir_D05G004880", "Ghir_D09G023320",
                "Ghir_A01G008830", "Ghir_A05G040190",
                "Ghir_A03G022510", "Ghir_D13G010120"
            )
        ) ~ replicate,
        scales = "free_y"
    ) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave("figure/geneRepCirca/auxinRelatedPlot.pdf", auxinRelatedPlot, width = 10, height = 10)
head(auxinRelatedGeneMtx)

## normalize
auxinRelatedGeneScaledMtx <-  WT_FL_0_2day_TMM_sample_exp %>%
    dplyr::select(c(1:6, {{ auxinRelatedGeneList }}))  %>% 
    mutate_if(is.double, scale) %>% 
    pivot_longer(cols = starts_with("Ghir"), names_to = "genes", values_to = "expression") %>%
    left_join(auxinRelatedTibble, by = c("genes" = "geneId"))

auxinRelatedGeneScaledPlot <- auxinRelatedGeneScaledMtx %>% 
    ggplot(aes(time, expression)) +
    geom_point(aes(col = genes)) +
    geom_smooth(aes(group = interaction(as.factor(replicate), strain, genes), color = genes, linetype = strain), span = 0.3) +
    facet_wrap(~replicate, nrow = 2) +
    scale_x_continuous(breaks = seq(1, 69, 4), ) +
    geom_rect(
        data = rects,
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
    labs(title = "Phase diff between Auxin related cascading regulation", subtitle = NULL) +
    ylab("Z-score of relative expression/(TMM)") +
    facet_grid(
        Families  ~ replicate,
        scales = "free_y"
    ) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave("figure/geneRepCirca/auxinRelatedGeneScaledPlot.pdf", auxinRelatedGeneScaledPlot, width = 10, height = 10)


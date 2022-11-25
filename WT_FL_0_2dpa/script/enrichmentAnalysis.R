library(clusterProfiler)
# library(AnnotationHub)
# hub <- AnnotationHub()
# q <- query(hub, "Gossypium")
# GossypiumDB <- q[["AH85399"]]
# GossypiumDB
# AnnotationDbi::columns(GossypiumDB)
# AnnotationDbi::keytypes(GossypiumDB)
# AnnotationDbi::keys(GossypiumDB, "SYMBOL")[1:10]

# 构建对照表
library(tidyverse)
library(readxl)
genes2Go <- read_excel("../bla_go/genes2Go.xlsx",
    col_names = T, skip = 1
) %>%
    mutate(
        geneid = str_extract(Query, "Ghir_\\w{10}"),
        domain = str_extract(Description, "(.*)(?=:)"), Specific_Description = str_extract(Description, "(?<=:)(.*)$")
    ) %>%
    dplyr::select("Match", "geneid", "Description", domain, Specific_Description)
head(genes2Go)
BP_genes2Go <- genes2Go %>%
    filter(domain == "Biological Process") %>%
    dplyr::select("Match", "geneid", Specific_Description)

CC_genes2Go <- genes2Go %>%
    filter(domain == "Cellular Component") %>%
    dplyr::select("Match", "geneid", Specific_Description)

MF_genes2Go <- genes2Go %>%
    filter(domain == "Molecular Function") %>%
    dplyr::select("Match", "geneid", Specific_Description)

BP_id2desc <- BP_genes2Go %>% dplyr::select("Match", Specific_Description)
CC_id2desc <- CC_genes2Go %>% dplyr::select("Match", Specific_Description)
MF_id2desc <- MF_genes2Go %>% dplyr::select("Match", Specific_Description)

## cluster1基因的富集分析
cluster1Genes  <- read_lines("figure/circ_WT_FL_TMM_FL_specific_mtx_kmeans_seed1_cluster1_Genes.txt")

enresult_cluster1Genes_BP <- enricher(cluster1Genes,
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.1,
    TERM2GENE = BP_genes2Go,
    TERM2NAME = BP_id2desc
)
barplot(enresult_cluster1Genes_BP, showCategory = 20 , xlab = "Count") + ggtitle("cluster1Genes_BP")
ggsave("figure/enresult_cluster1Genes_BP.pdf")
dotplot(enresult_cluster1Genes_BP, showCategory = 20) + ggtitle("cluster1Genes_BP")
ggsave("figure/enresult_cluster1Genes_BP_dotplot.pdf")

enresult_cluster1Genes_CC <- enricher(cluster1Genes,
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.1,
    TERM2GENE = CC_genes2Go,
    TERM2NAME = CC_id2desc
)
barplot(enresult_cluster1Genes_CC, showCategory = 20 , xlab = "Count") + ggtitle("cluster1Genes_CC")
ggsave("figure/enresult_cluster1Genes_CC.pdf")
dotplot(enresult_cluster1Genes_CC, showCategory = 20) + ggtitle("cluster1Genes_CC")
ggsave("figure/enresult_cluster1Genes_CC_dotplot.pdf")

enresult_cluster1Genes_MF <- enricher(cluster1Genes,
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.1,
    TERM2GENE = MF_genes2Go,
    TERM2NAME = MF_id2desc
)
barplot(enresult_cluster1Genes_MF, showCategory = 20 , xlab = "Count") + ggtitle("cluster1Genes_MF")
ggsave("figure/enresult_cluster1Genes_MF.pdf")
dotplot(enresult_cluster1Genes_MF, showCategory = 20) + ggtitle("cluster1Genes_MF")
ggsave("figure/enresult_cluster1Genes_MF_dotplot.pdf")

## WT circa gene
WT_specific_circGenes <- read_lines("mediumDataSave/WT_circ_genes_specific.txt")

enresult_WT_specific_circGenes_BP <- enricher(WT_specific_circGenes,
    pvalueCutoff = 0.5,
    qvalueCutoff = 0.5,
    TERM2GENE = BP_genes2Go,
    TERM2NAME = BP_id2desc
)
DNA_replication <- enresult_WT_specific_circGenes_BP %>% as_tibble() %>% filter(Description == "DNA replication")
write_tsv(DNA_replication, "mediumDataSave/WT_specific_GOBP_DNA_replication.txt")
dotplot(enresult_WT_specific_circGenes_BP, showCategory = 20 ) + ggtitle("WT_specific_circGenes_BP")
ggsave("figure/enresult_WT_specific_circGenes_BP.pdf")
enresult_WT_specific_circGenes_MF <- enricher(WT_specific_circGenes,
    pvalueCutoff = 0.5,
    qvalueCutoff = 0.5,
    TERM2GENE = MF_genes2Go,
    TERM2NAME = MF_id2desc
)
dotplot(enresult_WT_specific_circGenes_MF, showCategory = 20) + ggtitle("WT_specific_circGenes_MF")

# FL circa gene
FL_specific_circGenes <- read_lines("mediumDataSave/FL_circ_genes_specific.txt")

enresult_FL_specific_circGenes_BP <- enricher(FL_specific_circGenes,
    pvalueCutoff = 0.5,
    qvalueCutoff = 0.5,
    TERM2GENE = BP_genes2Go,
    TERM2NAME = BP_id2desc
)
dotplot(enresult_FL_specific_circGenes_BP, showCategory = 20) + ggtitle("FL_specific_circGenes_BP")
ggsave("figure/enresult_FL_specific_circGenes_BP.pdf")

# WT FL circa gene
WT_FL_specific_circGenes <- read_lines("mediumDataSave/WT_FL_circ_genes.txt")

enresult_WT_FL_specific_circGenes_BP <- enricher(WT_FL_specific_circGenes,
    pvalueCutoff = 0.5,
    qvalueCutoff = 0.5,
    TERM2GENE = BP_genes2Go,
    TERM2NAME = BP_id2desc
)
dotplot(enresult_WT_FL_specific_circGenes_BP, showCategory = 20) + ggtitle("WT_FL_specific_circGenes_BP")
ggsave("figure/enresult_WT_FL_specific_circGenes_BP.pdf")

enrichplot::upsetplot(enresult_WT_FL_specific_circGenes_BP, categorySize="pvalue")

# testTerm <- enresult_WT_FL_specific_circGenes_BP$Description[1:5]
# enrichplot::pmcplot(enresult_WT_FL_specific_circGenes_BP$Description[1:5], 2010:2022)


## 剪接事件富集分析
### A3SS WT FL
A3SS_WT_genes <- read_lines("mediumDataSave/cycleGenes/A3SS_WT_cycle_genes.txt") %>% map_chr(~str_split(., "_\\d+$", n = 2)[[1]][1])
A3SS_FL_genes <- read_lines("mediumDataSave/cycleGenes/A3SS_FL_cycle_genes.txt") %>% map_chr(~str_split(., "_\\d+$", n = 2)[[1]][1])
enresult_A3SS_WT_genes_BP <- enricher(A3SS_WT_genes,
    pvalueCutoff = 0.5,
    qvalueCutoff = 0.5,
    TERM2GENE = BP_genes2Go,
    TERM2NAME = BP_id2desc
)
dotplot(enresult_A3SS_WT_genes_BP, showCategory = 20) + ggtitle("A3SS_WT_genes_BP")
ggsave("figure/enrichRes/enresult_A3SS_WT_genes_BP.pdf")

enresult_A3SS_FL_genes_BP <- enricher(A3SS_FL_genes,
    pvalueCutoff = 0.5,
    qvalueCutoff = 0.5,
    TERM2GENE = BP_genes2Go,
    TERM2NAME = BP_id2desc
)
dotplot(enresult_A3SS_FL_genes_BP, showCategory = 20) + ggtitle("A3SS_FL_genes_BP")
ggsave("figure/enrichRes/enresult_A3SS_FL_genes_BP.pdf")
### A5SS WT FL
A5SS_WT_genes <- read_lines("mediumDataSave/cycleGenes/A5SS_WT_cycle_genes.txt") %>% map_chr(~str_split(., "_\\d+$", n = 2)[[1]][1])
A5SS_FL_genes <- read_lines("mediumDataSave/cycleGenes/A5SS_FL_cycle_genes.txt") %>% map_chr(~str_split(., "_\\d+$", n = 2)[[1]][1])

enresult_A5SS_WT_genes_BP <- enricher(A5SS_WT_genes,
    pvalueCutoff = 0.5,
    qvalueCutoff = 0.5,
    TERM2GENE = BP_genes2Go,
    TERM2NAME = BP_id2desc
)
dotplot(enresult_A5SS_WT_genes_BP, showCategory = 20) + ggtitle("A5SS_WT_genes_BP")
ggsave("figure/enrichRes/enresult_A5SS_WT_genes_BP.pdf")

enresult_A5SS_FL_genes_BP <- enricher(A5SS_FL_genes,
    pvalueCutoff = 0.5,
    qvalueCutoff = 0.5,
    TERM2GENE = BP_genes2Go,
    TERM2NAME = BP_id2desc
)
dotplot(enresult_A5SS_FL_genes_BP, showCategory = 20) + ggtitle("A5SS_FL_genes_BP")
ggsave("figure/enrichRes/enresult_A5SS_FL_genes_BP.pdf")
### MXE WT FL
MXE_WT_genes <- read_lines("mediumDataSave/cycleGenes/MXE_WT_cycle_genes.txt") %>% map_chr(~str_split(., "_\\d+$", n = 2)[[1]][1])
MXE_FL_genes <- read_lines("mediumDataSave/cycleGenes/MXE_FL_cycle_genes.txt") %>% map_chr(~str_split(., "_\\d+$", n = 2)[[1]][1])
### RI WT FL
RI_WT_genes <- read_lines("mediumDataSave/cycleGenes/RI_WT_cycle_genes.txt") %>% map_chr(~str_split(., "_\\d+$", n = 2)[[1]][1])
RI_FL_genes <- read_lines("mediumDataSave/cycleGenes/RI_FL_cycle_genes.txt") %>% map_chr(~str_split(., "_\\d+$", n = 2)[[1]][1])

enresult_RI_WT_genes_BP <- enricher(RI_WT_genes,
    pvalueCutoff = 0.5,
    qvalueCutoff = 0.5,
    TERM2GENE = BP_genes2Go,
    TERM2NAME = BP_id2desc
)
dotplot(enresult_RI_WT_genes_BP, showCategory = 20) + ggtitle("RI_WT_genes_BP")
ggsave("figure/enrichRes/enresult_RI_WT_genes_BP.pdf")

enresult_RI_FL_genes_BP <- enricher(RI_FL_genes,
    pvalueCutoff = 0.5,
    qvalueCutoff = 0.5,
    TERM2GENE = BP_genes2Go,
    TERM2NAME = BP_id2desc
)
dotplot(enresult_RI_FL_genes_BP, showCategory = 20) + ggtitle("RI_FL_genes_BP")
ggsave("figure/enrichRes/enresult_RI_FL_genes_BP.pdf")
### SE WT FL
SE_WT_genes <- read_lines("mediumDataSave/cycleGenes/SE_WT_cycle_genes.txt") %>% map_chr(~str_split(., "_\\d+$", n = 2)[[1]][1])
SE_FL_genes <- read_lines("mediumDataSave/cycleGenes/SE_FL_cycle_genes.txt") %>% map_chr(~str_split(., "_\\d+$", n = 2)[[1]][1])

enresult_SE_WT_genes_BP <- enricher(SE_WT_genes,
    pvalueCutoff = 0.5,
    qvalueCutoff = 0.5,
    TERM2GENE = BP_genes2Go,
    TERM2NAME = BP_id2desc
)
dotplot(enresult_SE_WT_genes_BP, showCategory = 20) + ggtitle("SE_WT_genes_BP")
ggsave("figure/enrichRes/enresult_SE_WT_genes_BP.pdf")

enresult_SE_FL_genes_BP <- enricher(SE_FL_genes,
    pvalueCutoff = 0.5,
    qvalueCutoff = 0.5,
    TERM2GENE = BP_genes2Go,
    TERM2NAME = BP_id2desc
)
dotplot(enresult_SE_FL_genes_BP, showCategory = 20) + ggtitle("SE_FL_genes_BP")
ggsave("figure/enrichRes/enresult_SE_FL_genes_BP.pdf")

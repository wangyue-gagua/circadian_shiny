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
library(dplyr)
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

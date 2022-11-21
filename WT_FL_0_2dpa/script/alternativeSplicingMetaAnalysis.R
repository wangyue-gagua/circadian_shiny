# 分析meta2d结果
library(pheatmap)
library(ggplot2)
## A3SS WT
rMATs_merge_A3SS_WT_META2d <- read_csv("METACYCLE/rMATs_merge_A3SS_WT_META2d/meta2d_rMATs_merge_A3SS_WT_exp.csv")
A3SS_WT_cycle_genes <- rMATs_merge_A3SS_WT_META2d %>% filter(meta2d_BH.Q < 0.01) %>% pull(CycID)
rMATs_merge_A3SS_WT_exp_NaRemoved <- read_csv("./mediumDataSave/removeNaMeta/rMATs_merge_A3SS_WT_exp.csv")
A3SS_WT_cycle_exp <- rMATs_merge_A3SS_WT_exp_NaRemoved %>% filter(GeneID %in% A3SS_WT_cycle_genes)

## A3SS FL
rMATs_merge_A3SS_FL_META2d <- read_csv("METACYCLE/rMATs_merge_A3SS_FL_META2d/meta2d_rMATs_merge_A3SS_FL_exp.csv")
A3SS_FL_cycle_genes <- rMATs_merge_A3SS_FL_META2d %>% filter(meta2d_BH.Q < 0.01) %>% pull(CycID)
rMATs_merge_A3SS_FL_exp_NaRemoved <- read_csv("./mediumDataSave/removeNaMeta/rMATs_merge_A3SS_FL_exp.csv")
A3SS_FL_cycle_exp <- rMATs_merge_A3SS_FL_exp_NaRemoved %>% filter(GeneID %in% A3SS_FL_cycle_genes)

### A3SS WT vs FL
A3SS_WT_FL_cycle_exp <- full_join(A3SS_WT_cycle_exp, A3SS_FL_cycle_exp, by = "GeneID", suffix = c("_WT", "_FL")) %>% 
column_to_rownames("GeneID")
A3SS_WT_FL_cycle_exp[is.na(A3SS_WT_FL_cycle_exp)] <- 0
pheatmap(
    A3SS_WT_FL_cycle_exp,
    cluster_cols = FALSE,
    # cluster_rows = FALSE,
    show_rownames = FALSE,
    show_colnames = TRUE,
    scale = "row",
    filename = "figure/alter/A3SS_WT_FL_cycle_exp.pdf"
)

## A5SS WT
rMATs_merge_A5SS_WT_META2d <- read_csv("METACYCLE/rMATs_merge_A5SS_WT_META2d/meta2d_rMATs_merge_A5SS_WT_exp.csv")
A5SS_WT_cycle_genes <- rMATs_merge_A5SS_WT_META2d %>% filter(meta2d_BH.Q < 0.01) %>% pull(CycID)
rMATs_merge_A5SS_WT_exp_NaRemoved <- read_csv("./mediumDataSave/removeNaMeta/rMATs_merge_A5SS_WT_exp.csv")
A5SS_WT_cycle_exp <- rMATs_merge_A5SS_WT_exp_NaRemoved %>% filter(GeneID %in% A5SS_WT_cycle_genes)

## A5SS FL
rMATs_merge_A5SS_FL_META2d <- read_csv("METACYCLE/rMATs_merge_A5SS_FL_META2d/meta2d_rMATs_merge_A5SS_FL_exp.csv")
A5SS_FL_cycle_genes <- rMATs_merge_A5SS_FL_META2d %>% filter(meta2d_BH.Q < 0.01) %>% pull(CycID)
rMATs_merge_A5SS_FL_exp_NaRemoved <- read_csv("./mediumDataSave/removeNaMeta/rMATs_merge_A5SS_FL_exp.csv")
A5SS_FL_cycle_exp <- rMATs_merge_A5SS_FL_exp_NaRemoved %>% filter(GeneID %in% A5SS_FL_cycle_genes)

### A5SS WT vs FL
A5SS_WT_FL_cycle_exp <- full_join(A5SS_WT_cycle_exp, A5SS_FL_cycle_exp, by = "GeneID", suffix = c("_WT", "_FL")) %>%
column_to_rownames("GeneID")
A5SS_WT_FL_cycle_exp[is.na(A5SS_WT_FL_cycle_exp)] <- 0
pheatmap(
    A5SS_WT_FL_cycle_exp,
    cluster_cols = FALSE,
    # cluster_rows = FALSE,
    show_rownames = FALSE,
    show_colnames = TRUE,
    scale = "row",
    filename = "figure/alter/A5SS_WT_FL_cycle_exp.pdf"
)

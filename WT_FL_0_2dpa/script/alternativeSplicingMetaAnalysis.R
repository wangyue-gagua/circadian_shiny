# 分析meta2d结果
library(pheatmap)
library(ggplot2)
## A3SS WT
rMATs_merge_A3SS_WT_META2d <- read_csv("METACYCLE/rMATs_merge_A3SS_WT_META2d/meta2d_rMATs_merge_A3SS_WT_exp.csv")
A3SS_WT_cycle_genes <- rMATs_merge_A3SS_WT_META2d %>%
    filter(meta2d_BH.Q < 0.01) %>%
    pull(CycID)
rMATs_merge_A3SS_WT_exp_NaRemoved <- read_csv("./mediumDataSave/removeNaMeta/rMATs_merge_A3SS_WT_exp.csv")
A3SS_WT_cycle_exp <- rMATs_merge_A3SS_WT_exp_NaRemoved %>% filter(GeneID %in% A3SS_WT_cycle_genes)

## A3SS FL
rMATs_merge_A3SS_FL_META2d <- read_csv("METACYCLE/rMATs_merge_A3SS_FL_META2d/meta2d_rMATs_merge_A3SS_FL_exp.csv")
A3SS_FL_cycle_genes <- rMATs_merge_A3SS_FL_META2d %>%
    filter(meta2d_BH.Q < 0.01) %>%
    pull(CycID)
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
    main = "A3SS WT vs FL, N = 491",
    filename = "figure/alter/A3SS_WT_FL_cycle_exp.pdf"
)

## A5SS WT
rMATs_merge_A5SS_WT_META2d <- read_csv("METACYCLE/rMATs_merge_A5SS_WT_META2d/meta2d_rMATs_merge_A5SS_WT_exp.csv")
A5SS_WT_cycle_genes <- rMATs_merge_A5SS_WT_META2d %>%
    filter(meta2d_BH.Q < 0.01) %>%
    pull(CycID)
rMATs_merge_A5SS_WT_exp_NaRemoved <- read_csv("./mediumDataSave/removeNaMeta/rMATs_merge_A5SS_WT_exp.csv")
A5SS_WT_cycle_exp <- rMATs_merge_A5SS_WT_exp_NaRemoved %>% filter(GeneID %in% A5SS_WT_cycle_genes)

## A5SS FL
rMATs_merge_A5SS_FL_META2d <- read_csv("METACYCLE/rMATs_merge_A5SS_FL_META2d/meta2d_rMATs_merge_A5SS_FL_exp.csv")
A5SS_FL_cycle_genes <- rMATs_merge_A5SS_FL_META2d %>%
    filter(meta2d_BH.Q < 0.01) %>%
    pull(CycID)
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
    main = "A5SS WT vs FL, N = 248",
    filename = "figure/alter/A5SS_WT_FL_cycle_exp.pdf"
)

## MXE WT
rMATs_merge_MXE_WT_META2d <- read_csv("METACYCLE/rMATs_merge_MXE_WT_META2d/meta2d_rMATs_merge_MXE_WT_exp.csv")
MXE_WT_cycle_genes <- rMATs_merge_MXE_WT_META2d %>%
    filter(meta2d_BH.Q < 0.01) %>%
    pull(CycID)
rMATs_merge_MXE_WT_exp_NaRemoved <- read_csv("./mediumDataSave/removeNaMeta/rMATs_merge_MXE_WT_exp.csv")
MXE_WT_cycle_exp <- rMATs_merge_MXE_WT_exp_NaRemoved %>% filter(GeneID %in% MXE_WT_cycle_genes)

## MXE FL
rMATs_merge_MXE_FL_META2d <- read_csv("METACYCLE/rMATs_merge_MXE_FL_META2d/meta2d_rMATs_merge_MXE_FL_exp.csv")
MXE_FL_cycle_genes <- rMATs_merge_MXE_FL_META2d %>%
    filter(meta2d_BH.Q < 0.01) %>%
    pull(CycID)
rMATs_merge_MXE_FL_exp_NaRemoved <- read_csv("./mediumDataSave/removeNaMeta/rMATs_merge_MXE_FL_exp.csv")
MXE_FL_cycle_exp <- rMATs_merge_MXE_FL_exp_NaRemoved %>% filter(GeneID %in% MXE_FL_cycle_genes)

### MXE WT vs FL
MXE_WT_FL_cycle_exp <- full_join(MXE_WT_cycle_exp, MXE_FL_cycle_exp, by = "GeneID", suffix = c("_WT", "_FL")) %>%
    column_to_rownames("GeneID")
MXE_WT_FL_cycle_exp[is.na(MXE_WT_FL_cycle_exp)] <- 0
pheatmap(
    MXE_WT_FL_cycle_exp,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    show_rownames = FALSE,
    show_colnames = TRUE,
    scale = "row",
    main = "MXE WT vs FL, N = 1",
    filename = "figure/alter/MXE_WT_FL_cycle_exp.pdf"
)

## RI WT
rMATs_merge_RI_WT_META2d <- read_csv("METACYCLE/rMATs_merge_RI_WT_META2d/meta2d_rMATs_merge_RI_WT_exp.csv")
RI_WT_cycle_genes <- rMATs_merge_RI_WT_META2d %>%
    filter(meta2d_BH.Q < 0.01) %>%
    pull(CycID)
rMATs_merge_RI_WT_exp_NaRemoved <- read_csv("./mediumDataSave/removeNaMeta/rMATs_merge_RI_WT_exp.csv")
RI_WT_cycle_exp <- rMATs_merge_RI_WT_exp_NaRemoved %>% filter(GeneID %in% RI_WT_cycle_genes)

## RI FL
rMATs_merge_RI_FL_META2d <- read_csv("METACYCLE/rMATs_merge_RI_FL_META2d/meta2d_rMATs_merge_RI_FL_exp.csv")
RI_FL_cycle_genes <- rMATs_merge_RI_FL_META2d %>%
    filter(meta2d_BH.Q < 0.01) %>%
    pull(CycID)
rMATs_merge_RI_FL_exp_NaRemoved <- read_csv("./mediumDataSave/removeNaMeta/rMATs_merge_RI_FL_exp.csv")
RI_FL_cycle_exp <- rMATs_merge_RI_FL_exp_NaRemoved %>% filter(GeneID %in% RI_FL_cycle_genes)

### RI WT vs FL
RI_WT_FL_cycle_exp <- full_join(RI_WT_cycle_exp, RI_FL_cycle_exp, by = "GeneID", suffix = c("_WT", "_FL")) %>%
    column_to_rownames("GeneID")
RI_WT_FL_cycle_exp[is.na(RI_WT_FL_cycle_exp)] <- 0
pheatmap(
    RI_WT_FL_cycle_exp,
    cluster_cols = FALSE,
    # cluster_rows = FALSE,
    show_rownames = FALSE,
    show_colnames = TRUE,
    scale = "row",
    main = "RI WT vs FL, N = 775",
    filename = "figure/alter/RI_WT_FL_cycle_exp.pdf"
)

## SE WT
rMATs_merge_SE_WT_META2d <- read_csv("METACYCLE/rMATs_merge_SE_WT_META2d/meta2d_rMATs_merge_SE_WT_exp.csv")
SE_WT_cycle_genes <- rMATs_merge_SE_WT_META2d %>%
    filter(meta2d_BH.Q < 0.01) %>%
    pull(CycID)
rMATs_merge_SE_WT_exp_NaRemoved <- read_csv("./mediumDataSave/removeNaMeta/rMATs_merge_SE_WT_exp.csv")
SE_WT_cycle_exp <- rMATs_merge_SE_WT_exp_NaRemoved %>% filter(GeneID %in% SE_WT_cycle_genes)

## SE FL
rMATs_merge_SE_FL_META2d <- read_csv("METACYCLE/rMATs_merge_SE_FL_META2d/meta2d_rMATs_merge_SE_FL_exp.csv")
SE_FL_cycle_genes <- rMATs_merge_SE_FL_META2d %>%
    filter(meta2d_BH.Q < 0.01) %>%
    pull(CycID)
rMATs_merge_SE_FL_exp_NaRemoved <- read_csv("./mediumDataSave/removeNaMeta/rMATs_merge_SE_FL_exp.csv")
SE_FL_cycle_exp <- rMATs_merge_SE_FL_exp_NaRemoved %>% filter(GeneID %in% SE_FL_cycle_genes)

### SE WT vs FL
SE_WT_FL_cycle_exp <- full_join(SE_WT_cycle_exp, SE_FL_cycle_exp, by = "GeneID", suffix = c("_WT", "_FL")) %>%
    column_to_rownames("GeneID")
SE_WT_FL_cycle_exp[is.na(SE_WT_FL_cycle_exp)] <- 0
pheatmap(
    SE_WT_FL_cycle_exp,
    cluster_cols = FALSE,
    # cluster_rows = FALSE,
    show_rownames = FALSE,
    show_colnames = TRUE,
    scale = "row",
    main = "SE WT vs FL, N = 461",
    filename = "figure/alter/SE_WT_FL_cycle_exp.pdf"
)

## 统计cycle事件比例

alterEventCount <- tibble(
    event = c("A3SS", "A5SS", "MXE", "RI", "SE"),
    WT = c(length(A3SS_WT_cycle_genes), length(A5SS_WT_cycle_genes), length(MXE_WT_cycle_genes), length(RI_WT_cycle_genes), length(SE_WT_cycle_genes)),
    FL = c(length(A3SS_FL_cycle_genes), length(A5SS_FL_cycle_genes), length(MXE_FL_cycle_genes), length(RI_FL_cycle_genes), length(SE_FL_cycle_genes))
) %>%
    pivot_longer(cols = c(WT, FL), names_to = "strain", values_to = "count")

alterEventCount %>% ggplot(aes(x = event, y = count, fill = strain)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw() +
    xlab("Event") +
    ylab("Number of events") +
    ggtitle("Number of events with circadian pattern in WT and FL")

ggsave("figure/alter/alterEventCount.pdf", width = 6, height = 4)

## 合并事件
merge_WT_exp <- rbind(A3SS_WT_cycle_exp, A5SS_WT_cycle_exp, RI_WT_cycle_exp, SE_WT_cycle_exp) %>%
    column_to_rownames("GeneID")
pheatmap(
    merge_WT_exp,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    show_rownames = FALSE,
    gaps_row = c(length(A3SS_WT_cycle_genes), length(A5SS_WT_cycle_genes), length(RI_WT_cycle_genes), length(SE_WT_cycle_genes)),
    scale = "row",
    main = "WT, N = 1,165",
    filename = "figure/alter/merge_WT_cycle_exp.pdf"
)

merge_FL_exp <- rbind(A3SS_FL_cycle_exp, A5SS_FL_cycle_exp, RI_FL_cycle_exp, SE_FL_cycle_exp) %>%
    column_to_rownames("GeneID")

library(ComplexHeatmap)
dend <- cluster_within_group(t(as.matrix(merge_WT_exp)), rep(
    c("A3SS", "A5SS", "RI", "SE"),
    c(length(A3SS_WT_cycle_genes), length(A5SS_WT_cycle_genes), length(RI_WT_cycle_genes), length(SE_WT_cycle_genes))
))

# library(circlize)

col_anno <- HeatmapAnnotation(ZT = rep(rep(c("light", "night"), c(4, 2)), 3), col = list(ZT = c("light" = "grey", "night" = "black")))

pdf("figure/alter/complex_merge_WT_cycle_exp.pdf")
ht <- Heatmap(t(scale(t(as.matrix(merge_WT_exp)))),
    cluster_columns = FALSE,
    row_split = c(rep("A3SS", length(A3SS_WT_cycle_genes)), rep("A5SS", length(A5SS_WT_cycle_genes)), rep("RI", length(RI_WT_cycle_genes)), rep("SE", length(SE_WT_cycle_genes))),
    show_row_names = FALSE,
    top_annotation = col_anno,
    gap = unit(5, "mm"),
    name = "WT Proportion"
    # cluster_rows = dend,
    # row_split = 4
)
top_annotation <- draw(ht)
dev.off()

pdf("figure/alter/complex_merge_FL_cycle_exp.pdf")
ht <- Heatmap(t(scale(t(as.matrix(merge_FL_exp)))),
    cluster_columns = FALSE,
    row_split = c(rep("A3SS", length(A3SS_FL_cycle_genes)), rep("A5SS", length(A5SS_FL_cycle_genes)), rep("RI", length(RI_FL_cycle_genes)), rep("SE", length(SE_FL_cycle_genes))),
    show_row_names = FALSE,
    top_annotation = col_anno,
    gap = unit(5, "mm"),
    name = "FL Proportion"
    # cluster_rows = dend,
    # row_split = 4
)
top_annotation <- draw(ht)
dev.off()


# select row by row order
# sort(row_order(ht)$A3SS)

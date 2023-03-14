# 此脚本用于得到metacycle分析结果后的后续分析

## 寻找不受节律调控的基因
### 1. 读取detrend后的表达量数据
WT_0_2day_genes_TMM_EXPR_mergeRep_selected_detrended <- read_csv("./mediumDataSave/WT_0_2day_genes_TMM_EXPR_mergeRep_selected_detrended.csv")
FL_0_2day_genes_TMM_EXPR_mergeRep_selected_detrended <- read_csv("./mediumDataSave/FL_0_2day_genes_TMM_EXPR_mergeRep_selected_detrended.csv")
#### 7w个基因不能直接画热图，看来还是别画总体热图了 Error in hclust(d, method = method) : size cannot be NA nor exceed 65536 
### pheatmap绘制detrend后WT FL表达量热图
library(pheatmap)
WT_FL_0_2day_genes_TMM_EXPR_mergeRep_selected_detrended <- left_join(WT_0_2day_genes_TMM_EXPR_mergeRep_selected_detrended,
    FL_0_2day_genes_TMM_EXPR_mergeRep_selected_detrended,
    by = "Geneid"
) %>% column_to_rownames("Geneid") %>% 
    mutate(means = rowMeans(., na.rm = T), 
           sd = apply(., 1, sd, na.rm = T))

### 2. 读取detrended后metacycle结果
WT_meta2d_detrended <- read_csv("WT_meta2d_detrended/JTKresult_WT_0_2day_genes_TMM_EXPR_mergeRep_selected_detrended.csv")
FL_meta2d_detrended <- read_csv("FL_meta2d_detrended/JTKresult_FL_0_2day_genes_TMM_EXPR_mergeRep_selected_detrended.csv")

WT_uncirca_genes <- WT_meta2d_detrended %>% filter(ADJ.P > 0.5) %>% pull(CycID)
FL_uncirca_genes <- FL_meta2d_detrended %>% filter(ADJ.P > 0.5) %>% pull(CycID)

WT_FL_both_uncirca_genes <- intersect(WT_uncirca_genes, FL_uncirca_genes)

### 3. 读取原始表达量数据，detrend会使得基因表达量的均值变为1
#### geneExpDiff.R 引入second_count_matrix metacyle.R 引入metaInfo_WT_FL_0_2day_TMM_sample_exp

second_count_matrix_withMeanAndSd <- second_count_matrix %>%
    column_to_rownames("Geneid") %>%
    mutate(means = rowMeans(., na.rm = T), 
           sd = apply(., 1, sd, na.rm = T)) %>% 
    rownames_to_column("Geneid")
#### 存入文件
write_csv(second_count_matrix_withMeanAndSd, "./mediumDataSave/second_count_matrix_withMeanAndSd.csv")
second_count_matrix_withMeanAndSd <- read_csv("./mediumDataSave/second_count_matrix_withMeanAndSd.csv")
second_detected_geneCnt <- second_count_matrix_withMeanAndSd %>% 
    filter(means > 10) %>% 
    pull(Geneid)

WT_FL_both_uncirca_detected_genes <- intersect(WT_FL_both_uncirca_genes, second_detected_geneCnt)

second_count_matrix_WT_FL_both_uncirca_detected_genes <- second_count_matrix %>% 
    filter(Geneid %in% WT_FL_both_uncirca_detected_genes) %>% 
    column_to_rownames("Geneid")

pheatmap(
    second_count_matrix_WT_FL_both_uncirca_detected_genes,
    cluster_rows = T,
    cluster_cols = F,
    show_rownames = F,
    scale = "row",
    filename = "figure/detrended/WT_FL_both_uncirca_detected_genes_pheatmap.pdf",
)

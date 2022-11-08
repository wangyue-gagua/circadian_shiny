library(PCAtools)
sampleTimeTibble <- WT_FL_0_2day_TMM_sample_exp %>%
  select(time, sample, strain, replicate) %>%
  group_by(strain) %>%
  arrange(time, .by_group = T)

WT_FL_0_2day_genes_TMM_EXPR_TimeSorted <- WT_FL_0_2day_genes_TMM_EXPR %>% select(sampleTimeTibble$sample)
pcaMeta <- sampleTimeTibble %>% column_to_rownames("sample")

p <- pca(WT_FL_0_2day_genes_TMM_EXPR_TimeSorted, meta = pcaMeta)
screeplot(p)

biplot(p, x = 'PC1', y = 'PC2', shape = 'strain', colby = "time",
    legendPosition = 'top', lab = NULL, showLoadings = TRUE,
    title = 'PCA of WT FL 0-2dpa', subtitle = 'TMM normalized expression')
ggsave("figure/samplesCor/PCA_WT_FL_0_2day_TMM_EXPR_TimeSorted.pdf", width = 10, height = 10)
# attributes(p)
# p$metadata

plotloadings(p, labSize = 3)
plotRepCirca("Ghir_A10G015270")
ggsave("figure/geneRepCirca/Ghir_A10G015270_repCirca.pdf", width = 10, height = 10)

head(WT_FL_0_2day_genes_TMM_EXPR_TimeSorted)
hc <- hclust(dist(WT_FL_0_2day_genes_TMM_EXPR_TimeSorted), method = "average")
plot(hc, main = "Hierarchical Clustering Dendrogram", xlab = "Genes", sub = "", cex = 0.5)
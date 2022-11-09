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
library(fastcluster)
hc <- fastcluster::hclust(dist(t(WT_FL_0_2day_genes_TMM_EXPR_TimeSorted)))
# BiocManager::install("YuLab-SMU/treedataverse")
install.packages("https://cran.r-project.org/src/contrib/Archive/rvcheck/rvcheck_0.1.8.tar.gz", repos = NULL) # 降级，不然报错
BiocManager::install("ggtree")
library(ggtree)
ggtree(hc)
clus <- cutree(hc, 4)
g <- split(names(clus), clus)

p <- ggtree(hc, linetype='dashed')
clades <- sapply(g, function(n) MRCA(p, n))
p <- groupClade(p, clades, group_name='subtree') + aes(color=subtree)

d <- data.frame(label = names(clus), 
                  time = sampleTimeTibble$time,
                  phase = factor(sampleTimeTibble$time %% 24),
                  strain = sampleTimeTibble$strain)

p %<+% d + 
  layout_dendrogram() + 
  geom_tippoint(aes(fill=phase, x=x+.5),
                size=5, color = "black", shape = 21) + 
  geom_tiplab(aes(label=time), size=3, vjust = 0, hjust = .5, color='black') +
  geom_tiplab(angle=90, hjust=1, offset=-10, show.legend=FALSE) + 
  scale_color_brewer(palette='Set1', breaks=1:4) +
  theme_dendrogram(plot.margin=margin(6,6,80,6)) +
  theme(legend.position=c(.9, .6))


ggsave("figure/samplesCor/hclust_WT_FL_0_2day_TMM_EXPR_TimeSorted_noStrainShape.pdf", width = 10, height = 10)

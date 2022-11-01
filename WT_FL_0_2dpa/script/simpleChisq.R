library(tidyverse)
TE <- read_tsv("/data1/wangy/data/ribo/result.txt")

TE_diff <- TE %>% filter(log2FoldChange > 1 & pvalue < 0.05)

AD <- read_tsv("/data1/wangy/data/ribo/WT_all_genes_table.txt") %>% mutate(
  gene_id = map(X5, function(str){str_extract(str, "Ga\\d{2}g\\d{5}")})
)

res <- intersect(TE_diff$geneid, AD$gene_id)


chisq.test(matrix(c(6, 29, 2288, 34503), nrow = 2))
fisher.test(matrix(c(6, 29, 2288, 34503), nrow = 2))

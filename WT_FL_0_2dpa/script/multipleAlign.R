# visualize multiple align result
library(ggmsa)
library(patchwork)
Ghir_D04G000830_protein_seq <- "mediumDataSave/misc/Ghir_D04G000830_isoform_protein_align.fa"

## 总长度1107
ggmsa(Ghir_D04G000830_protein_seq, start = 1, end = 1107, char_width = 0.5, seq_name = T)
p1 <- ggmsa(Ghir_D04G000830_protein_seq, start = 1, end = 50, char_width = 0.5, seq_name = T)
p2 <- ggmsa(Ghir_D04G000830_protein_seq, start = 250, end = 300, char_width = 0.5, seq_name = T)
p3 <- ggmsa(Ghir_D04G000830_protein_seq, start = 1057, end = 1107, char_width = 0.5, seq_name = T)
p1 / p2 / p3
ggsave("figure/alter/genesAlign/Ghir_D04G000830_protein.pdf", width = 10, height = 10)

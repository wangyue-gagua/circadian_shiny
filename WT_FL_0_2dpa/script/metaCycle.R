setwd("WT_FL_0_2dpa")

# data import
library(tidyverse)
library(data.table)
# WT_FL_0_2day_TMM_sample_exp <- read_csv("../merged_counts/WT_FL_0_2day_TMM_sample_exp.csv")
WT_FL_0_2day_TMM_sample_exp <- fread("../merged_counts/WT_FL_0_2day_TMM_sample_exp.csv")
WT_FL_0_2day_TMM_sample_exp  <- WT_FL_0_2day_TMM_sample_exp %>% as_tibble()

WT_FL_0_2day_genes_TMM_EXPR <- read_delim("../merged_counts/WT_FL_0_2day_genes.TMM.EXPR.matrix",
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
) %>% column_to_rownames(var = "Geneid")

WT_0_2day_genes_TMM_EXPR <- WT_FL_0_2day_genes_TMM_EXPR %>%
  select(contains("WT"))

FL_0_2day_genes_TMM_EXPR <- WT_FL_0_2day_genes_TMM_EXPR %>%
  select(contains("FL"))

WT_0_2day_genes_TMM_EXPR_mergeRep <- WT_0_2day_genes_TMM_EXPR
FL_0_2day_genes_TMM_EXPR_mergeRep <- FL_0_2day_genes_TMM_EXPR


sample_list <- WT_FL_0_2day_TMM_sample_exp %>%
  group_by(strain) %>%
  arrange(time) %>%
  select(sample)

for (i in 1:length(sample_list$strain)) {
  name <- str_remove(sample_list$sample[[i]], "\\d$")
  if (sample_list$strain[[i]] == "WT") {
    WT_0_2day_genes_TMM_EXPR_mergeRep <- WT_0_2day_genes_TMM_EXPR_mergeRep %>%
      mutate(!!name := rowMeans(select(WT_0_2day_genes_TMM_EXPR, contains(name))))
  }
  if (sample_list$strain[[i]] == "FL") {
    FL_0_2day_genes_TMM_EXPR_mergeRep <- FL_0_2day_genes_TMM_EXPR_mergeRep %>%
      mutate(!!name := rowMeans(select(FL_0_2day_genes_TMM_EXPR, contains(name))))
  }
}

WT_0_2day_genes_TMM_EXPR_mergeRep_selected <- WT_0_2day_genes_TMM_EXPR_mergeRep[37:length(WT_0_2day_genes_TMM_EXPR_mergeRep)] %>%
  rownames_to_column(var = "Geneid")
FL_0_2day_genes_TMM_EXPR_mergeRep_selected <- FL_0_2day_genes_TMM_EXPR_mergeRep[37:length(FL_0_2day_genes_TMM_EXPR_mergeRep)] %>%
  rownames_to_column(var = "Geneid")

library(MetaCycle)
write.csv(WT_0_2day_genes_TMM_EXPR_mergeRep_selected, file = "WT_0_2day_genes_TMM_EXPR_mergeRep_selected.csv", row.names = F)
write.csv(FL_0_2day_genes_TMM_EXPR_mergeRep_selected, file = "FL_0_2day_genes_TMM_EXPR_mergeRep_selected.csv", row.names = F)
meta2d(
  infile = "WT_0_2day_genes_TMM_EXPR_mergeRep_selected.csv",
  filestyle = "csv",
  timepoints = seq(1, 69, 4),
  outdir = "WT_meta2d",
  parallelize = TRUE,
  nCores = 10,
  minper = 24,
  maxper = 24,
  ARSdefaultPer = 24,
  cycMethod = "JTK",
)

meta2d(
  infile = "FL_0_2day_genes_TMM_EXPR_mergeRep_selected.csv",
  filestyle = "csv",
  timepoints = seq(1, 69, 4),
  outdir = "FL_meta2d",
  parallelize = TRUE,
  nCores = 10,
  minper = 24,
  maxper = 24,
  ARSdefaultPer = 24,
  cycMethod = "JTK",
)
# downstream analysis data loading
library(tidyverse)
WT_meta2d <- read_csv("WT_meta2d/JTKresult_WT_0_2day_genes_TMM_EXPR_mergeRep_selected.csv")
FL_meta2d <- read_csv("FL_meta2d/JTKresult_FL_0_2day_genes_TMM_EXPR_mergeRep_selected.csv")

WT_dat <- WT_meta2d %>% filter(BH.Q < 0.01)
FL_dat <- FL_meta2d %>% filter(BH.Q < 0.01)

WT_0_2day_genes_TMM_EXPR_mergeRep_selected <- read_csv("WT_0_2day_genes_TMM_EXPR_mergeRep_selected.csv", )
FL_0_2day_genes_TMM_EXPR_mergeRep_selected <- read_csv("FL_0_2day_genes_TMM_EXPR_mergeRep_selected.csv", )
# tempWT <- WT_0_2day_genes_TMM_EXPR_mergeRep_selected %>% mutate(expSum = rowSums(select(., -Geneid))) # WT 中检测到50000w个基因表达
# sum(tempWT$expSum) / nrow(WT_0_2day_genes_TMM_EXPR_mergeRep_selected) # 267.1212
# tempFL <- FL_0_2day_genes_TMM_EXPR_mergeRep_selected %>% mutate(expSum = rowSums(select(., -Geneid))) # FL 中检测到50000w个基因表达
# sum(tempFL$expSum) / nrow(FL_0_2day_genes_TMM_EXPR_mergeRep_selected) # 247.5492
# FL_0_2day_genes_TMM_EXPR_mergeRep_selected %>% mutate(expSum = rowSums(select(., -Geneid))) %>% filter(expSum > 3) %>% nrow() # FL 中检测到49954w个基因表达
circ_WT_TMM_mtx <- WT_0_2day_genes_TMM_EXPR_mergeRep_selected %>%
  filter(Geneid %in% WT_dat$CycID) %>%
  column_to_rownames(var = "Geneid")
circ_FL_TMM_mtx <- FL_0_2day_genes_TMM_EXPR_mergeRep_selected %>%
  filter(Geneid %in% FL_dat$CycID) %>%
  column_to_rownames(var = "Geneid")

# plot
library(pheatmap)
pheatmap(circ_WT_TMM_mtx,
  show_rownames = F,
  cluster_cols = F,
  scale = "row",
  filename = "circ_WT_TMM_mtx.pdf",
)

pheatmap(circ_FL_TMM_mtx,
  show_rownames = F,
  cluster_cols = F,
  scale = "row",
  filename = "circ_FL_TMM_mtx.pdf",
)

## specific gene list
WT_circ_genes <- WT_dat$CycID
FL_circ_genes <- FL_dat$CycID

WT_circ_genes_specific <- setdiff(WT_circ_genes, FL_circ_genes)
FL_circ_genes_specific <- setdiff(FL_circ_genes, WT_circ_genes)
WT_FL_circ_genes <- intersect(WT_circ_genes, FL_circ_genes)
write_lines(WT_circ_genes_specific, "mediumDataSave/WT_circ_genes_specific.txt")
write_lines(FL_circ_genes_specific, "mediumDataSave/FL_circ_genes_specific.txt")
write_lines(WT_FL_circ_genes, "mediumDataSave/WT_FL_circ_genes.txt")
## pheatmap FL specifc合并WT FL
circ_WT_TMM_FL_specific_mtx <- WT_0_2day_genes_TMM_EXPR_mergeRep_selected %>%
  filter(Geneid %in% FL_circ_genes_specific)
circ_FL_TMM_FL_specific_mtx <- FL_0_2day_genes_TMM_EXPR_mergeRep_selected %>%
  filter(Geneid %in% FL_circ_genes_specific)

circ_WT_FL_TMM_FL_specific_mtx <- full_join(circ_WT_TMM_FL_specific_mtx, circ_FL_TMM_FL_specific_mtx, by = "Geneid") %>% 
  column_to_rownames(var = "Geneid")

pheatmap(circ_WT_FL_TMM_FL_specific_mtx,
  show_rownames = F,
  cluster_cols = F,
  scale = "row",
  filename = "figure/circ_WT_FL_TMM_FL_specific_mtx.pdf",
)
### kmeans cluster
set.seed(1)
circ_WT_FL_TMM_FL_specific_mtx_kmeans_seed1 <- pheatmap::pheatmap(circ_WT_FL_TMM_FL_specific_mtx,
  show_rownames = T,
  cluster_cols = F,
  scale = "row",
  kmeans_k = 8,
  # filename = "figure/circ_WT_FL_TMM_FL_specific_mtx_kmeans_seed1.pdf",
)
circ_WT_FL_TMM_FL_specific_mtx_kmeans_seed1_cluster <- circ_WT_FL_TMM_FL_specific_mtx_kmeans_seed1$kmeans$cluster
circ_WT_FL_TMM_FL_specific_mtx_kmeans_seed1_cluster1_Genes <- names(circ_WT_FL_TMM_FL_specific_mtx_kmeans_seed1_cluster[circ_WT_FL_TMM_FL_specific_mtx_kmeans_seed1_cluster==1])
# write_lines(circ_WT_FL_TMM_FL_specific_mtx_kmeans_seed1_cluster1_Genes, "figure/circ_WT_FL_TMM_FL_specific_mtx_kmeans_seed1_cluster1_Genes.txt")

circ_WT_TMM_FL_specific_cluster1_mtx <- WT_0_2day_genes_TMM_EXPR_mergeRep_selected %>%
  filter(Geneid %in% circ_WT_FL_TMM_FL_specific_mtx_kmeans_seed1_cluster1_Genes)
circ_FL_TMM_FL_specific_cluster1_mtx <- FL_0_2day_genes_TMM_EXPR_mergeRep_selected %>%
  filter(Geneid %in% circ_WT_FL_TMM_FL_specific_mtx_kmeans_seed1_cluster1_Genes)

circ_WT_FL_TMM_FL_specific_cluster1_mtx <- full_join(circ_WT_TMM_FL_specific_cluster1_mtx, circ_FL_TMM_FL_specific_cluster1_mtx, by = "Geneid") %>%
  column_to_rownames(var = "Geneid")

pheatmap(circ_WT_FL_TMM_FL_specific_cluster1_mtx,
  show_rownames = F,
  cluster_cols = F,
  scale = "row",
  filename = "figure/circ_WT_FL_TMM_FL_specific_cluster1_mtx.pdf",
)

## pheatmap WT specifc合并WT FL
circ_WT_TMM_WT_specific_mtx <- WT_0_2day_genes_TMM_EXPR_mergeRep_selected %>%
  filter(Geneid %in% WT_circ_genes_specific)
circ_FL_TMM_WT_specific_mtx <- FL_0_2day_genes_TMM_EXPR_mergeRep_selected %>%
  filter(Geneid %in% WT_circ_genes_specific)

circ_WT_FL_TMM_WT_specific_mtx <- full_join(circ_WT_TMM_WT_specific_mtx, circ_FL_TMM_WT_specific_mtx, by = "Geneid") %>%
  column_to_rownames(var = "Geneid")

pheatmap(circ_WT_FL_TMM_WT_specific_mtx,
  show_rownames = F,
  cluster_cols = F,
  scale = "row",
  filename = "figure/circ_WT_FL_TMM_WT_specific_mtx.pdf",
)


# 对meta结果中的节律基因做GO富集分析
library(readxl)
genes2Go <- read_excel("../bla_go/genes2Go.xlsx",
  col_names = T, skip = 1
) %>%
  mutate(geneid = str_extract(Query, "Ghir_\\w{10}")) %>%
  select("Match", "geneid", "Description")

id2desc <- genes2Go %>% select("Match", "Description")



library(VennDiagram)
venn.diagram(
  x = list(
    WT = WT_circ_genes,
    FL = FL_circ_genes
  ),
  filename = "circ_genes.png",
  imagetype = "png"
)

library(clusterProfiler)
enresult_WT_specific <- enricher(WT_circ_genes_specific,
  pvalueCutoff = 0.5,
  qvalueCutoff = 0.5,
  TERM2GENE = genes2Go,
  TERM2NAME = id2desc
)
pdf("WT_circ_genes_specific_enrich.pdf")
dotplot(enresult_WT_specific, title = "WT_circ_genes_specific_enrich")
dev.off()

enresult_FL_specific <- enricher(FL_circ_genes_specific,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  TERM2GENE = genes2Go,
  TERM2NAME = id2desc
)
pdf(file = "FL_circ_genes_specific_enrich.pdf")
dotplot(enresult_FL_specific, title = "FL_circ_genes_specific_enrich")
dev.off()

enresult_WT_FL <- enricher(WT_FL_circ_genes,
  pvalueCutoff = 0.5,
  qvalueCutoff = 0.5,
  TERM2GENE = genes2Go,
  TERM2NAME = id2desc
)
pdf("WT_FL_circ_genes_enrich.pdf")
dotplot(enresult_WT_FL, title = "WT_FL_circ_genes_enrich")
dev.off()


# 寻找看家基因，通过计算方差与均值
temp <- WT_0_2day_genes_TMM_EXPR_mergeRep_selected %>%
  rowwise() %>%
  mutate(mean = mean(c_across(where(is.numeric))), sd = sd(c_across(where(is.numeric))))

conserve <- temp %>% filter(sd < 0.5 && mean >= 10)


# C3 specific genes
c3_specific_genes_tbl <- read_tsv("/data1/wangy/data/ref/Gh/C3_specific_genes.txt")

c3_specific_circ_genes <- intersect(FL_circ_genes, c3_specific_genes_tbl$Gid)

#
WT_meta2d %>% filter(CycID == "Ghir_D11G029140")
FL_meta2d %>% filter(CycID == "Ghir_D11G029140")


# (Ax+b)sinωx+k

df <- tibble(x = seq(1, 10, by = 0.01), y = (sin(x) + 2) * (x))

ggplot(df, aes(x, y)) +
  geom_line() +
  geom_point()


Ghir_D01G016160_exp <- WT_0_2day_genes_TMM_EXPR_mergeRep_selected %>%
  filter(Geneid == "Ghir_D01G016160") %>%
  select(-Geneid)
plot(x = seq(1, 69, by = 4), y = Ghir_D01G016160_exp)
plot(x = seq(1, 69, by = 4), y = Ghir_D01G016160_exp / seq(1, 69, by = 4))

#
exp <- WT_0_2day_genes_TMM_EXPR_mergeRep_selected %>%
  filter(Geneid == "Ghir_D01G016160") %>%
  select(-Geneid) %>%
  as_vector()
df <- tibble(x = seq(1, 69, by = 4), y = exp)
plot(df$x, df$y)
nlc <- nls.control(maxiter = 1000)
exp_model <- nls(y ~ x * (sin(a * x) + b), data = df, start = list(a = 1, b = 1), control = nlc)
exp_model_sum <- summary(exp_model)

ggplot(df, aes(x, y)) +
  geom_point() +
  geom_line(aes(x, fitted(exp_model)))

mean(resid(exp_model))

attributes(exp_model_sum)
exp_model_sum[[1]]
sd(resid(exp_model))

lhy_exp <- WT_0_2day_genes_TMM_EXPR_mergeRep_selected %>%
  filter(Geneid == "Ghir_D12G012900") %>%
  select(-Geneid) %>%
  as_vector()
lhy_df <- tibble(x = seq(1, 69, by = 4), y = lhy_exp)
plot(lhy_df$x, lhy_df$y)

lhy_exp_model <- nls(y ~ x * (sin(a * x) + b), data = lhy_df, start = list(a = 1, b = 1), control = nlc)
summary(lhy_exp_model)
sd(resid(lhy_exp_model))

ggplot(lhy_df, aes(x, y)) +
  geom_point() +
  geom_line(aes(x, fitted(lhy_exp_model)))

# 按前后段表达量平均值进行筛选
WT_2dpa_high_exp_table <- WT_0_2day_genes_TMM_EXPR_mergeRep_selected %>%
  filter((WT_2DPA_9pm + WT_2DPA_1am + WT_2DPA_5am) / (WT_0DPA_9am + WT_0DPA_1pm + WT_0DPA_5pm) > 10 & (WT_0DPA_9am + WT_0DPA_1pm + WT_0DPA_5pm) > 1)


WT_2dpa_high_exp_TMM_mtx <- WT_0_2day_genes_TMM_EXPR_mergeRep_selected %>%
  filter(Geneid %in% WT_2dpa_high_exp_table$Geneid) %>%
  column_to_rownames(var = "Geneid")

# plot
library(pheatmap)
pheatmap(WT_2dpa_high_exp_TMM_mtx,
  show_rownames = F,
  cluster_cols = F,
  # scale = "row",
  # filename = "WT_2dpa_high_exp_TMM_mtx.pdf",
)

FL_2dpa_high_exp_TMM_mtx <- FL_0_2day_genes_TMM_EXPR_mergeRep_selected %>%
  filter(Geneid %in% WT_2dpa_high_exp_table$Geneid) %>%
  column_to_rownames(var = "Geneid")

pheatmap(FL_2dpa_high_exp_TMM_mtx,
  show_rownames = F,
  cluster_cols = F,
  # cluster_rows = F,
  # scale = "row",
  filename = "FL_2dpa_high_exp_TMM_mtx.pdf",
)
WT_2dpa_high_exp_C3_specific_genes <- intersect(WT_2dpa_high_exp_table$Geneid, c3_specific_genes_tbl$Gid)
WT_2dpa_high_exp_circa_genes <- intersect(WT_2dpa_high_exp_table$Geneid, FL_circ_genes)

WT_2dpa_high_exp_C3_circa_genes <- intersect(WT_2dpa_high_exp_circa_genes, c3_specific_circ_genes)

chisq.test(matrix(c(39, 236, 300, 69624), nrow = 2))
chisq.test(matrix(c(17, 258, 5319, 64605), nrow = 2))


## 样本间相关性
library(ggcorrplot)
cor_df <- cor(WT_FL_0_2day_genes_TMM_EXPR)
cor_plot <- ggcorrplot(cor_df, tl.cex = 4, lab = T, lab_size = 1)
cor_plot_scale <- ggcorrplot(cor_df, tl.cex = 4, lab = T, lab_size = 1,
 digits = 3, title = "样本间相关性", outline.color = "transparent") +
scale_fill_gradient(limit = c(0.25, 1), low = "white", high = "red")

ggsave("figure/cor_plot.pdf", cor_plot, device = "pdf")
ggsave("figure/cor_plot_scale.pdf", cor_plot_scale, device = "pdf")

WT_2dpa_high_exp_TMM_mtx_genes <- rownames(WT_2dpa_high_exp_TMM_mtx)

WT_2dpa_high_exp_TMM_mtx_genes_chrom <- str_extract(WT_2dpa_high_exp_TMM_mtx_genes, "[AD]\\d{2}")


WT_2dpa_high_exp_TMM_mtx_genes_chrom_mtx <- as_tibble(table(WT_2dpa_high_exp_TMM_mtx_genes_chrom)) %>%
  mutate(chrom = WT_2dpa_high_exp_TMM_mtx_genes_chrom, WT_ratio = n / sum(n))
WT_2dpa_high_exp_TMM_mtx_genes_chrom_mtx %>% ggplot() +
  geom_col(aes(chrom, WT_ratio)) +
  # decline axis text
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


all_genes_id <- FL_meta2d$CycID
all_genes_chrom <- str_extract(all_genes_id, "[AD]\\d{2}")
all_genes_chrom_mtx <- as_tibble(table(all_genes_chrom)) %>%
  mutate(chrom = all_genes_chrom, all_ratio = n / sum(n))

all_merge_chrom_mtx <- inner_join(all_genes_chrom_mtx, WT_2dpa_high_exp_TMM_mtx_genes_chrom_mtx, by = "chrom") %>%
  mutate(WT_all_ratio = WT_ratio / all_ratio)

all_merge_chrom_mtx %>% ggplot() +
  geom_col(aes(chrom, WT_all_ratio)) +
  geom_hline(yintercept = 1, color = "red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

WT_0_2day_genes_TMM_EXPR_mergeRep_selected_ratio <- tibble(
  Geneid = WT_0_2day_genes_TMM_EXPR_mergeRep_selected$Geneid,
  WT_0DPA_9am_WT_0DPA_1pm = (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_0DPA_1pm + .Machine$double.eps) / (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_0DPA_9am + .Machine$double.eps),
  WT_0DPA_1pm_WT_0DPA_5pm = (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_0DPA_5pm + .Machine$double.eps) / (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_0DPA_1pm + .Machine$double.eps),
  WT_0DPA_5pm_WT_0DPA_9pm = (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_0DPA_9pm + .Machine$double.eps) / (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_0DPA_5pm + .Machine$double.eps),
  WT_0DPA_9pm_WT_0DPA_1am = (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_0DPA_1am + .Machine$double.eps) / (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_0DPA_9pm + .Machine$double.eps),
  WT_0DPA_1am_WT_0DPA_5am = (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_0DPA_5am + .Machine$double.eps) / (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_0DPA_1am + .Machine$double.eps),
  WT_0DPA_5am_WT_1DPA_9am = (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_1DPA_9am + .Machine$double.eps) / (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_0DPA_5am + .Machine$double.eps),
  WT_1DPA_9am_WT_1DPA_1pm = (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_1DPA_1pm + .Machine$double.eps) / (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_1DPA_9am + .Machine$double.eps),
  WT_1DPA_1pm_WT_1DPA_5pm = (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_1DPA_5pm + .Machine$double.eps) / (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_1DPA_1pm + .Machine$double.eps),
  WT_1DPA_5pm_WT_1DPA_9pm = (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_1DPA_9pm + .Machine$double.eps) / (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_1DPA_5pm + .Machine$double.eps),
  WT_1DPA_9pm_WT_1DPA_1am = (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_1DPA_1am + .Machine$double.eps) / (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_1DPA_9pm + .Machine$double.eps),
  WT_1DPA_1am_WT_1DPA_5am = (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_1DPA_5am + .Machine$double.eps) / (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_1DPA_1am + .Machine$double.eps),
  WT_1DPA_5am_WT_2DPA_9am = (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_2DPA_9am + .Machine$double.eps) / (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_1DPA_5am + .Machine$double.eps),
  WT_2DPA_9am_WT_2DPA_1pm = (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_2DPA_1pm + .Machine$double.eps) / (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_2DPA_9am + .Machine$double.eps),
  WT_2DPA_1pm_WT_2DPA_5pm = (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_2DPA_5pm + .Machine$double.eps) / (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_2DPA_1pm + .Machine$double.eps),
  WT_2DPA_5pm_WT_2DPA_9pm = (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_2DPA_9pm + .Machine$double.eps) / (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_2DPA_5pm + .Machine$double.eps),
  WT_2DPA_9pm_WT_2DPA_1am = (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_2DPA_1am + .Machine$double.eps) / (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_2DPA_9pm + .Machine$double.eps),
  WT_2DPA_1am_WT_2DPA_5am = (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_2DPA_5am + .Machine$double.eps) / (WT_0_2day_genes_TMM_EXPR_mergeRep_selected$WT_2DPA_1am + .Machine$double.eps),
)
WT_0_2day_genes_TMM_EXPR_mergeRep_selected_ratio_mtx <- WT_0_2day_genes_TMM_EXPR_mergeRep_selected_ratio %>%
  rowwise() %>%
  mutate(mean = mean(c_across(where(is.numeric))), sd = sd(c_across(where(is.numeric)))) %>%
  filter(sd > 0)
WT_0_2day_genes_TMM_EXPR_mergeRep_selected_ratio_mtx_sd_quant <- quantile(WT_0_2day_genes_TMM_EXPR_mergeRep_selected_ratio_mtx$sd, c(0.25, 0.5, 0.75))
WT_0_2day_genes_TMM_EXPR_mergeRep_selected_ratio_mtx %>% filter(sd < WT_0_2day_genes_TMM_EXPR_mergeRep_selected_ratio_mtx_sd_quant[[1]])


WT_0_2day_genes_TMM_EXPR_mergeRep_selected_ratio_2dayhighexp <- mtx <- WT_0_2day_genes_TMM_EXPR_mergeRep_selected_ratio_mtx %>%
  filter(Geneid %in% rownames(WT_2dpa_high_exp_TMM_mtx)) %>%
  arrange(sd) %>%
  head(n = 30)

WT_0_2day_genes_TMM_EXPR_mergeRep_selected_ratio_2dayhighexp_chrom <- str_extract(WT_0_2day_genes_TMM_EXPR_mergeRep_selected_ratio_2dayhighexp$Geneid, "[AD]\\d{2}")
WT_0_2day_genes_TMM_EXPR_mergeRep_selected_ratio_2dayhighexp_chrom_mtx <- as_tibble(table(WT_0_2day_genes_TMM_EXPR_mergeRep_selected_ratio_2dayhighexp_chrom)) %>%
  mutate(chrom = WT_0_2day_genes_TMM_EXPR_mergeRep_selected_ratio_2dayhighexp_chrom, WT_ratio = n / sum(n))
WT_0_2day_genes_TMM_EXPR_mergeRep_selected_ratio_2dayhighexp_chrom_mtx %>% ggplot() +
  geom_col(aes(chrom, WT_ratio)) +
  # decline axis text
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

all_merge_chrom_mtx <- left_join(all_genes_chrom_mtx, WT_0_2day_genes_TMM_EXPR_mergeRep_selected_ratio_2dayhighexp_chrom_mtx, by = "chrom") %>%
  mutate(WT_all_ratio = WT_ratio / all_ratio)

all_merge_chrom_mtx %>% ggplot() +
  geom_col(aes(chrom, WT_all_ratio)) +
  geom_hline(yintercept = 1, color = "red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

WT_0_2day_genes_TMM_EXPR_mergeRep_selected_ratio_2dayhighexp_matrix <- WT_0_2day_genes_TMM_EXPR_mergeRep_selected %>%
  filter(Geneid %in% WT_0_2day_genes_TMM_EXPR_mergeRep_selected_ratio_2dayhighexp$Geneid) %>%
  column_to_rownames(var = "Geneid")
pheatmap(WT_0_2day_genes_TMM_EXPR_mergeRep_selected_ratio_2dayhighexp_matrix,
  show_rownames = T,
  cluster_cols = F,
  scale = "row",
  # filename = "WT_2dpa_high_exp_TMM_mtx.pdf",
)



# cosinor2 Comparison of cosinor parameters of two populations

library(cosinor2)
library(ggplot2)
library(circacompare)
set.seed(42)
# 合并FL和WT的meta2d结果，并通过circacompare添加mesor和amplitude以及phase的差异和p值
FL_WT_meta2d <- left_join(FL_meta2d, WT_meta2d, by = "CycID", suffix = c("_FL", "_WT"))

metaInfo_WT_FL_0_2day_TMM_sample_exp <- WT_FL_0_2day_TMM_sample_exp %>% select(1:6)

# 首先筛选出节律基因，避免模型不收敛
FL_WT_meta2d_circa <- FL_WT_meta2d %>%
  filter(BH.Q_FL < 0.01 & BH.Q_WT < 0.01) %>%
  mutate(messorDiff = c(0), amplitudeDiff = c(0), phaseDiff = c(0), messorDiffPvalue = c(0), amplitudeDiffPvalue = c(0), phaseDiffPvalue = c(0))

for (i in seq(1, length(FL_WT_meta2d_circa$CycID))) {
  geneId <- FL_WT_meta2d_circa$CycID[[i]]
  df <- cbind(metaInfo_WT_FL_0_2day_TMM_sample_exp, WT_FL_0_2day_TMM_sample_exp[geneId])
  colnames(df)[7] <- "measure"
  outDf <- circacompare_mixed(
    x = df,
    col_time = "time",
    col_group = "strain",
    col_outcome = "measure",
    col_id = "replicate",
    control = list(grouped_params = c("k", "alpha", "phi"), random_params = c("phi1")),
    period = 24
  )
  if(is.null(outDf)) {
    next
  }
  messorDiff <- outDf$summary %>%
    filter(parameter == "Mesor difference estimate") %>%
    pull(value)
  amplitudeDiff <- outDf$summary %>%
    filter(parameter == "Amplitude difference estimate") %>%
    pull(value)
  phaseDiff <- outDf$summary %>%
    filter(parameter == "Phase difference estimate") %>%
    pull(value)
  messorDiffPvalue <- outDf$summary %>%
    filter(parameter == "P-value for mesor difference") %>%
    pull(value)
  amplitudeDiffPvalue <- outDf$summary %>%
    filter(parameter == "P-value for amplitude difference") %>%
    pull(value)
  phaseDiffPvalue <- outDf$summary %>%
    filter(parameter == "P-value for difference in phase") %>%
    pull(value)
  FL_WT_meta2d_circa$messorDiff[[i]] <- messorDiff
  FL_WT_meta2d_circa$amplitudeDiff[[i]] <- amplitudeDiff
  FL_WT_meta2d_circa$phaseDiff[[i]] <- phaseDiff
  FL_WT_meta2d_circa$messorDiffPvalue[[i]] <- messorDiffPvalue
  FL_WT_meta2d_circa$amplitudeDiffPvalue[[i]] <- amplitudeDiffPvalue
  FL_WT_meta2d_circa$phaseDiffPvalue[[i]] <- phaseDiffPvalue
  print(i)
}

write_tsv(FL_WT_meta2d_circa, "FL_WT_meta2d_circa.tsv")
FL_WT_meta2d_circa <- read_tsv("FL_WT_meta2d_circa.tsv")

FL_WT_meta2d_circa %>% filter(abs(phaseDiff) > 8) # 相位差最大的两个基因 Ghir_D11G030190 Ghir_A10G017040
# 统计相位差异分布
phaseDiffTable <- tibble(
  PhaseDiff = c("<1", "1 ~ 4", "4 ~ 8", "8 ~ 12"),
  nums = c(sum(abs(FL_WT_meta2d_circa$phaseDiff) < 1),
           sum(abs(FL_WT_meta2d_circa$phaseDiff) >= 1 & abs(FL_WT_meta2d_circa$phaseDiff) < 4),
           sum(abs(FL_WT_meta2d_circa$phaseDiff) >= 4 & abs(FL_WT_meta2d_circa$phaseDiff) < 8),
           sum(abs(FL_WT_meta2d_circa$phaseDiff) >= 8 & abs(FL_WT_meta2d_circa$phaseDiff) < 12)),
)
label_value <- paste(' h(', round(phaseDiffTable$nums/sum(phaseDiffTable$nums) * 100, 1), '%)', sep = '')
label <- paste(phaseDiffTable$PhaseDiff, label_value, sep = '')
ggplot(phaseDiffTable, aes(x = "Content", y = nums, fill = PhaseDiff)) +
  geom_bar(stat = "identity") +
  coord_polar(theta = "y") +
  scale_fill_discrete(labels = label) +
  labs(x = '', y = '', title = "相位差异分布 N = 1149") +
  theme(axis.text = element_blank()) +
  theme(panel.grid=element_blank())


## plot
source("script/plotFunctionUtil.R")
# plot "Ghir_D12G020190" "Ghir_A12G005280" "Ghir_A13G013150" "Ghir_D07G006230" "Ghir_D02G019940" "Ghir_D08G002570"
p1 <- plotRepCirca("Ghir_D12G020190")
p2 <- plotRepCirca("Ghir_A12G005280", "AT5G26830.1 / The mRNA is cell-to-cell mobile")
p3 <- plotRepCirca("Ghir_A13G013150", "AT5G61170")
p4 <- plotRepCirca("Ghir_D07G006230", "AT4G11630")
p5 <- plotRepCirca("Ghir_D02G019940", "AT2G43090 / The mRNA is cell-to-cell mobile")
p6 <- plotRepCirca("Ghir_D08G002570", "AT5G57330")

p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(guides = "collect")
ggsave("figure/circ_WT_FL_TMM_FL_specific_mtx_kmeans_seed1_cluster1_Genes_6gene_test.pdf")

## core circadian genes
my_cir_plot_WT_FL_0_2dpa("Ghir_A04G012270", "LNK1")
ggsave("figure/my_cir_plot_WT_FL_0_2dpa/Ghir_A04G012270_LNK1.pdf")
my_cir_plot_WT_FL_0_2dpa("Ghir_A05G028660", "LNK2")
ggsave("figure/my_cir_plot_WT_FL_0_2dpa/Ghir_A05G028660_LNK2.pdf")
my_cir_plot_WT_FL_0_2dpa("Ghir_D05G009940", "PRR3")
ggsave("figure/my_cir_plot_WT_FL_0_2dpa/Ghir_D05G009940_PRR3.pdf")
PRR5 <- my_grid_plot_WT_FL_0_2dpa(
    str_split("Ghir_A11G001650/Ghir_D11G001640/Ghir_A12G025940/Ghir_D12G025960/Ghir_A05G042880", "/")[[1]],
    "PRR5"
)
ggsave("figure/my_cir_plot_WT_FL_0_2dpa/Ghir_A11G001650_Ghir_D11G001640_Ghir_A12G025940_Ghir_D12G025960_Ghir_A05G042880_PRR5.pdf", PRR5)
PRR7 <- my_grid_plot_WT_FL_0_2dpa(
    str_split("Ghir_A09G016930/Ghir_D09G016400/Ghir_A11G035130/Ghir_D11G036010/Ghir_D04G008730/Ghir_A05G035540", "/")[[1]],
    "PRR7"
)
ggsave("figure/my_cir_plot_WT_FL_0_2dpa/Ghir_A09G016930_Ghir_D09G016400_Ghir_A11G035130_Ghir_D11G036010_Ghir_D04G008730_Ghir_A05G035540_PRR7.pdf", PRR7)
PRR9 <- my_grid_plot_WT_FL_0_2dpa(
    str_split("Ghir_A11G010900/Ghir_D11G010820", "/")[[1]],
    "PRR9"
)
ggsave("figure/my_cir_plot_WT_FL_0_2dpa/Ghir_A11G010900_Ghir_D11G010820_PRR9.pdf", PRR9)

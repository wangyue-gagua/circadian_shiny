plotRepCirca("Ghir_D12G009360")
my_cir_plot_WT_FL_0_2dpa("Ghir_D12G009360")

WT_FL_0_2day_TMM_sample_exp_detrended <-  WT_FL_0_2day_TMM_sample_exp  %>% group_by(strain, replicate, period) %>%
    mutate(across(where(is.double), ~ . / (mean(.) + .Machine$double.eps)))

write_csv(WT_FL_0_2day_TMM_sample_exp_detrended, "./mediumDataSave/WT_FL_0_2day_TMM_sample_exp_detrended.csv")

# randomSixCluster1Genens <- c("Ghir_D12G020190", "Ghir_A12G005280", "Ghir_A13G013150", "Ghir_D07G006230", "Ghir_D02G019940", "Ghir_D08G002570")

p1 <- plotRepCircaDetrended("Ghir_D12G020190")
p2 <- plotRepCircaDetrended("Ghir_A12G005280")
p3 <- plotRepCircaDetrended("Ghir_A13G013150")
p4 <- plotRepCircaDetrended("Ghir_D07G006230")
p5 <- plotRepCircaDetrended("Ghir_D02G019940")
p6 <- plotRepCircaDetrended("Ghir_D08G002570")

# patchwork collect legend
p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(guides = "collect")
ggsave("figure/geneRepCirca/WT_FL_0_2day_detrended_cluster1_Genes_6gene_detrended.pdf")

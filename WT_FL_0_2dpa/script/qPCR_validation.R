library(readxl)
pRR5_qPCR_path <- "qPCR_rawData/PRR5_2_5_qPCR.xlsx"

pRR5_qPCR_tbl <- read_xlsx(pRR5_qPCR_path)

pRR5_qPCR_tbl_cq <- pRR5_qPCR_tbl %>%
    mutate(Sample = toupper(Sample)) %>%
    select(Sample, Target, Cq) %>%
    mutate(
        strain = toupper(str_split_i(Sample, "_", 1)),
        day = str_split_i(Sample, "_", 2),
        phase = str_split_i(Sample, "_", 3),
        rep = rep(c(1, 2, 3), length(pRR5_qPCR_tbl$Sample) / 3)
    ) %>%
    pivot_wider(names_from = Target, values_from = Cq) %>%
    mutate(deltaCq = PRR5 - UBQ)

controlList <- rep(pRR5_qPCR_tbl_cq$deltaCq[1:3], length(pRR5_qPCR_tbl_cq$deltaCq) / 3)
pRR5_qPCR_tbl_quant <- pRR5_qPCR_tbl_cq %>%
    mutate(
        deDeltaCq = deltaCq - controlList,
        quant = 2^(-deDeltaCq)
    ) %>%
    group_by(Sample, strain, day, phase) %>%
    summarise(means = mean(quant), sd = sd(quant)) %>%
    mutate(timeString = str_c(day, "Dpa", phase))
# head(pRR5_qPCR_tbl_quant)


pRR5_qPCR_tbl_quant %>% ggplot(aes(timeString, means, color = strain)) +
    geom_point() +
    geom_errorbar(aes(ymin = means - sd, ymax = means + sd), width = 0.2) +
    geom_line(aes(group = strain)) +
    labs(x = "Time", y = "Relative Expression", title = "PRR5 qPCR validation") +
    theme_bw()

ggsave("figure/qPCR_validation/PRR5_2_5.pdf")

# 去除离群值
pRR5_qPCR_tbl_cq_normal <- pRR5_qPCR_tbl_cq %>%
    mutate(
        deDeltaCq = deltaCq - controlList,
        quant = 2^(-deDeltaCq)
    ) %>%
    group_by(Sample, strain, day, phase) %>%
    mutate(
        upQuant = quantile(deltaCq, 0.75) + 1.5 * IQR(deltaCq),
        downQuant = quantile(deltaCq, 0.25) - 1.5 * IQR(deltaCq),
    ) %>%
    filter(deltaCq < upQuant & deltaCq > downQuant)


pRR5_qPCR_tbl_cq_filtered <- pRR5_qPCR_tbl_cq_normal %>%
    ungroup() %>%
    slice(-c(2, 6, 9, 11, 15, 18, 19, 22, 27, 28, 33, 34)) %>%
    mutate(deltaCq = PRR5 - UBQ)

slice(pRR5_qPCR_tbl_cq_normal, c(1, 28))

filterCntList <- rep(pRR5_qPCR_tbl_cq_filtered$deltaCq[1:2], length(pRR5_qPCR_tbl_cq_filtered$deltaCq) / 2)

pRR5_qPCR_tbl_cq_filtered_plot <- pRR5_qPCR_tbl_cq_filtered %>%
    ungroup() %>%
    mutate(
        deDeltaCq = deltaCq - filterCntList,
        quant = 2^(-deDeltaCq)
    ) %>%
    group_by(Sample, strain, day, phase) %>%
    summarise(means = mean(quant), sd = sd(quant)) %>%
    mutate(timeString = str_c(day, "Dpa", phase))

pRR5_qPCR_tbl_cq_filtered_plot %>% ggplot(aes(timeString, means, color = strain)) +
    geom_point() +
    geom_errorbar(aes(ymin = means - sd, ymax = means + sd), width = 0.2) +
    geom_line(aes(group = strain)) +
    labs(x = "Time", y = "Relative Expression", title = "PRR5 qPCR validation") +
    theme_bw()

ggsave("figure/qPCR_validation/PRR5_2_5_filtered.pdf")


## PCR second experiment validation PRR5

PRR5_second_1Period_path <- "qPCR_rawData/pRR5_second_1Period.xls"

PRR5_second_1Period_tbl <- read_xls(PRR5_second_1Period_path)

PRR5_second_1Period_tbl_cq <- PRR5_second_1Period_tbl %>%
    mutate(Sample = toupper(Sample)) %>%
    select(Sample, Target, Cq) %>%
    mutate(
        strain = toupper(str_split_i(Sample, "_", 1)),
        day = str_split_i(Sample, "_", 2),
        phase = str_split_i(Sample, "_", 3),
        rep = rep(c(1, 2, 3), length(PRR5_second_1Period_tbl$Sample) / 3)
    ) %>%
    pivot_wider(names_from = Target, values_from = Cq) %>%
    mutate(deltaCq = PRR5 - UBQ) %>%
    arrange(phase)


controlList <- rep(PRR5_second_1Period_tbl_cq$deltaCq[4:6], length(PRR5_second_1Period_tbl_cq$deltaCq) / 3) # FL_0_1 似乎更稳定

PRR5_second_1Period_tbl_cq_normal <- PRR5_second_1Period_tbl_cq %>%
    mutate(
        deDeltaCq = deltaCq - controlList,
        quant = 2^(-deDeltaCq)
    ) %>%
    group_by(Sample, strain, day, phase) %>%
    mutate(
        upQuant = quantile(deltaCq, 0.75) + 1.5 * IQR(deltaCq),
        downQuant = quantile(deltaCq, 0.25) - 1.5 * IQR(deltaCq),
    ) %>%
    filter(deltaCq < upQuant & deltaCq > downQuant)

# 绘图
PRR5_second_1Period_tbl_cq_normal_plot <- PRR5_second_1Period_tbl_cq %>%
    ungroup() %>%
    mutate(
        deDeltaCq = deltaCq - controlList,
        quant = 2^(-deDeltaCq)
    ) %>%
    group_by(Sample, strain, day, phase) %>%
    summarise(means = mean(quant), sd = sd(quant)) %>%
    mutate(timeString = str_c(day, "Dpa", phase))

PRR5_second_1Period_tbl_cq_normal_plot %>% ggplot(aes(timeString, means, color = strain)) +
    geom_point() +
    geom_errorbar(aes(ymin = means - sd, ymax = means + sd), width = 0.2) +
    geom_line(aes(group = strain)) +
    labs(x = "Time", y = "Relative Expression", title = "PRR5 qPCR Validation 1 Period") +
    theme_bw()

ggsave("figure/qPCR_validation/PRR5_second_1Period.pdf")

## 怀疑PRR5 PCR定量与RNA seq的峰值差异与Ubquitin的表达量呈现周期性差异有关

### 运行metacycle.R脚本获取表达量
### PRR5: Ghir_A05G042880  ubiquitin4: Ghir_D13G015430
WT_RNA_seq_PRR5 <- WT_0_2day_genes_TMM_EXPR_mergeRep_selected %>%
    filter(Geneid == "Ghir_A05G042880") %>%
    select(c(1:7)) %>%
    pivot_longer(cols = -Geneid, names_to = "Sample", values_to = "quant")
WT_RNA_seq_UBQ4 <- WT_0_2day_genes_TMM_EXPR_mergeRep_selected %>%
    filter(Geneid == "Ghir_D13G015430") %>%
    select(c(1:7)) %>%
    pivot_longer(cols = -Geneid, names_to = "Sample", values_to = "quant")

FL_RNA_seq_PRR5 <- FL_0_2day_genes_TMM_EXPR_mergeRep_selected %>%
    filter(Geneid == "Ghir_A05G042880") %>%
    select(c(1:7)) %>%
    pivot_longer(cols = -Geneid, names_to = "Sample", values_to = "quant")
FL_RNA_seq_UBQ4 <- FL_0_2day_genes_TMM_EXPR_mergeRep_selected %>%
    filter(Geneid == "Ghir_D13G015430") %>%
    select(c(1:7)) %>%
    pivot_longer(cols = -Geneid, names_to = "Sample", values_to = "quant")

WT_RNA_seq_PRR5_UBQ_relative <- tibble(
    Sample = factor(WT_RNA_seq_PRR5$Sample, levels = WT_RNA_seq_PRR5$Sample),
    relativeQuant = WT_RNA_seq_PRR5$quant / WT_RNA_seq_UBQ4$quant,
    scaledRelativeExpression = scale(WT_RNA_seq_PRR5$quant / WT_RNA_seq_UBQ4$quant)[, 1],
    PRR5_original = WT_RNA_seq_PRR5$quant,
    UBQ_original = WT_RNA_seq_UBQ4$quant
) %>%
    mutate(
        strain = str_split_i(Sample, "_", 1),
    )

WT_RNA_seq_PRR5_UBQ_relative %>%
    ggplot(aes(x = Sample, y = relativeQuant)) +
    geom_point() +
    geom_line(aes(group = strain)) +
    labs(x = "Sample", y = "Relative Expression", title = "PRR5 Relative Expression PRR5/UBQ") +
    theme_bw()

ggsave("figure/qPCR_validation/PRR5_RNA_seq_relative.pdf")

WT_RNA_seq_PRR5_UBQ_relative %>%
    ggplot(aes(x = Sample, y = PRR5_original)) +
    geom_point() +
    geom_line(aes(group = strain)) +
    labs(x = "Sample", y = "Relative Expression", title = "PRR5 Relative Expression PRR5 Original") +
    theme_bw()

ggsave("figure/qPCR_validation/PRR5_RNA_seq_original.pdf")

WT_RNA_seq_PRR5_UBQ_relative %>%
    ggplot(aes(x = Sample, y = UBQ_original)) +
    geom_point() +
    geom_line(aes(group = strain)) +
    labs(x = "Sample", y = "Relative Expression", title = "PRR5 Relative Expression UBQ Original") +
    theme_bw()

ggsave("figure/qPCR_validation/UBQ_RNA_seq_original.pdf")


WT_RNA_seq_PRR5_UBQ_relative %>%
    ggplot(aes(x = Sample, y = scaledRelativeExpression)) +
    geom_point() +
    geom_line(aes(group = strain)) +
    labs(x = "Sample", y = "Relative Expression", title = "PRR5 Relative Expression PRR5/UBQ") +
    theme_bw()

#### FL
FL_RNA_seq_PRR5_UBQ_allType <- tibble(
    Sample = factor(FL_RNA_seq_PRR5$Sample, levels = FL_RNA_seq_PRR5$Sample),
    relativeQuant = FL_RNA_seq_PRR5$quant / FL_RNA_seq_UBQ4$quant,
    scaledRelativeExpression = scale(FL_RNA_seq_PRR5$quant / FL_RNA_seq_UBQ4$quant)[, 1],
    PRR5_original = FL_RNA_seq_PRR5$quant,
    UBQ_original = FL_RNA_seq_UBQ4$quant
) %>%
    mutate(
        strain = str_split_i(Sample, "_", 1),
    )

FL_RNA_seq_PRR5_UBQ_allType %>%
    ggplot(aes(x = Sample, y = relativeQuant)) +
    geom_point() +
    geom_line(aes(group = strain)) +
    labs(x = "Sample", y = "Relative Expression", title = "PRR5 Relative Expression PRR5/UBQ") +
    theme_bw()

ggsave("figure/qPCR_validation/FL_PRR5_RNA_seq_relative.pdf")

FL_RNA_seq_PRR5_UBQ_allType %>%
    ggplot(aes(x = Sample, y = PRR5_original)) +
    geom_point() +
    geom_line(aes(group = strain)) +
    labs(x = "Sample", y = "Relative Expression", title = "PRR5 Relative Expression PRR5 Original") +
    theme_bw()

ggsave("figure/qPCR_validation/FL_PRR5_RNA_seq_original.pdf")

FL_RNA_seq_PRR5_UBQ_allType %>%
    ggplot(aes(x = Sample, y = UBQ_original)) +
    geom_point() +
    geom_line(aes(group = strain)) +
    labs(x = "Sample", y = "Relative Expression", title = "PRR5 Relative Expression UBQ Original") +
    theme_bw()

ggsave("figure/qPCR_validation/FL_UBQ_RNA_seq_original.pdf")


FL_PRR5_UBQ_RNASEQ_PCR <- tibble(
    Sample = FL_RNA_seq_PRR5_UBQ_allType$Sample,
    PRR5_RNAseq = FL_RNA_seq_PRR5_UBQ_allType$scaledRelativeExpression,
    PRR5_qPCR = scale(filter(PRR5_second_1Period_tbl_cq_normal_plot, strain == "FL")$means)[,1]
) %>% 
    pivot_longer(cols = -Sample, names_to = "RNAseq_qPCR", values_to = "quant")


FL_PRR5_UBQ_RNASEQ_PCR %>% 
    ggplot(aes(x = Sample, y = quant, color = RNAseq_qPCR)) +
    geom_point() +
    geom_line(aes(group = RNAseq_qPCR)) +
    labs(x = "RNAseq_qPCR", y = "Relative Expression", title = "PRR5 Relative Expression PRR5/UBQ") +
    theme_bw()

ggsave("figure/qPCR_validation/FL_PRR5_RNA_seq_qPCR_contrast.pdf")


# PCR validation PHD1 with SCF as ref gene
PHD1_SCF_path <- "qPCR_rawData/PHD1_SCFUBQ.xlsx"

PHD1_SCF_tbl <- read_excel(PHD1_SCF_path)

PHD1_SCF_tbl_cq <- PHD1_SCF_tbl %>%
    mutate(Sample = toupper(Sample)) %>%
    select(Sample, Target, Cq) %>%
    mutate(
        strain = toupper(str_split_i(Sample, "_", 1)),
        day = str_split_i(Sample, "_", 2),
        phase = str_split_i(Sample, "_", 3),
        rep = rep(c(1, 2, 3), length(PHD1_SCF_tbl$Sample) / 3)
    ) %>%
    pivot_wider(names_from = Target, values_from = Cq) %>%
    mutate(deltaCq = PHD1 - SCF_UBQ)

PHD1_SCF_controlList <- rep(PHD1_SCF_tbl_cq$deltaCq[1:3], length(PHD1_SCF_tbl_cq$deltaCq) / 3)

PHD1_SCF_qPCR_tbl_quant <- PHD1_SCF_tbl_cq %>%
    mutate(
        deDeltaCq = deltaCq - PHD1_SCF_controlList,
        quant = 2^(-deDeltaCq)
    ) %>%
    group_by(Sample, strain, day, phase) %>%
    summarise(means = mean(quant, na.rm = TRUE), sd = sd(quant, na.rm = TRUE)) %>%
    mutate(timeString = str_c(day, "Dpa", phase))

PHD1_SCF_qPCR_tbl_quant %>%
    ggplot(aes(x = timeString, y = means, color = strain)) +
    geom_point() +
    geom_errorbar(aes(ymin = means - sd, ymax = means + sd), width = 0.2) +
    geom_line(aes(group = strain)) +
    labs(x = "Time", y = "Relative Expression", title = "PHD1 Relative Expression PHD1/SCF") +
    theme_bw()

ggsave("figure/qPCR_validation/PHD1_SCF_qPCR.pdf")

# PCR validation PHD1 with SCF as ref gene / repeat experiment
PHD1_SCF_Rep_path <- "qPCR_rawData/PHD1_SCFUBQ_REP2.xlsx"
PHD1_SCF_Rep_tbl <- readxl::read_excel(PHD1_SCF_Rep_path)

PHD1_SCF_Rep_tbl_cq <- PHD1_SCF_Rep_tbl %>%
    mutate(Sample = toupper(Sample)) %>%
    select(Sample, Target, Cq) %>%
    mutate(
        strain = toupper(str_split_i(Sample, "_", 1)),
        day = str_split_i(Sample, "_", 2),
        phase = str_split_i(Sample, "_", 3),
        rep = rep(c(1, 2, 3), length(PHD1_SCF_Rep_tbl$Sample) / 3)
    ) %>%
    pivot_wider(names_from = Target, values_from = Cq) %>%
    mutate(deltaCq = PHD1 - SCF_UBQ)

PHD1_SCF_Rep_controlList <- rep(PHD1_SCF_Rep_tbl_cq$deltaCq[1:3], length(PHD1_SCF_Rep_tbl_cq$deltaCq) / 3)

PHD1_SCF_Rep_qPCR_tbl_quant <- PHD1_SCF_Rep_tbl_cq %>%
    mutate(
        deDeltaCq = deltaCq - PHD1_SCF_Rep_controlList,
        quant = 2^(-deDeltaCq)
    ) %>%
    group_by(Sample, strain, day, phase) %>%
    summarise(means = mean(quant, na.rm = TRUE), sd = sd(quant, na.rm = TRUE)) %>%
    mutate(timeString = str_c(day, "Dpa", phase))

PHD1_SCF_Rep_qPCR_tbl_quant %>%
    ggplot(aes(x = timeString, y = means, color = strain)) +
    geom_point() +
    geom_errorbar(aes(ymin = means - sd, ymax = means + sd), width = 0.2) +
    geom_line(aes(group = strain)) +
    labs(x = "Time", y = "Relative Expression", title = "PHD1 Relative Expression PHD1/SCF") +
    theme_bw()

ggsave("figure/qPCR_validation/PHD1_SCF_qPCR_Rep.pdf")

# PCR validation PHD1 with SCF as ref gene / repeat 2 experiment
PHD1_SCFUBQ_REP3_with_PRR5_path <- "qPCR_rawData/PHD1_SCFUBQ_REP3_with_PRR5.xlsx"
PHD1_SCFUBQ_REP3_with_PRR5_tbl <- readxl::read_excel(PHD1_SCFUBQ_REP3_with_PRR5_path) %>% 
    mutate(Sample = toupper(Sample)) %>%
    mutate(
        strain = toupper(str_split_i(Sample, "_", 1)),
        day = str_split_i(Sample, "_", 2),
        phase = str_split_i(Sample, "_", 3),
    )

# 拆分PHD1与PRR5
PHD1_SCFUBQ_REP3 <- PHD1_SCFUBQ_REP3_with_PRR5_tbl %>%
        filter(day == 1) %>% 
        select(Sample, Target, Cq, strain, day, phase)
with_PRR5 <- PHD1_SCFUBQ_REP3_with_PRR5_tbl %>%
        filter(day == 0) %>% 
        select(Sample, Target, Cq, strain, day, phase)

PHD1_SCFUBQ_REP3_cq <- PHD1_SCFUBQ_REP3 %>%
    mutate(
        rep = rep(c(1, 2, 3), length(PHD1_SCFUBQ_REP3$Sample) / 3)
    ) %>%
    pivot_wider(names_from = Target, values_from = Cq) %>%
    mutate(deltaCq = PHD1 - SCF_UBQ)

PHD1_SCFUBQ_REP3_controlList <- rep(PHD1_SCFUBQ_REP3_cq$deltaCq[1:3], length(PHD1_SCFUBQ_REP3_cq$deltaCq) / 3)
PHD1_SCFUBQ_REP3_quant <- PHD1_SCFUBQ_REP3_cq %>%
    mutate(
        deDeltaCq = deltaCq - PHD1_SCFUBQ_REP3_controlList,
        quant = 2^(-deDeltaCq)
    ) %>%
    group_by(Sample, strain, day, phase) %>%
    summarise(means = mean(quant, na.rm = TRUE), sd = sd(quant, na.rm = TRUE)) %>%
    mutate(timeString = str_c(day, "Dpa", phase))

PHD1_SCFUBQ_REP3_quant %>%
    ggplot(aes(x = timeString, y = means, color = strain)) +
    geom_point() +
    geom_errorbar(aes(ymin = means - sd, ymax = means + sd), width = 0.2) +
    geom_line(aes(group = strain)) +
    labs(x = "Time", y = "Relative Expression", title = "PHD1 Relative Expression PHD1/SCF") +
    theme_bw()

ggsave("figure/qPCR_validation/PHD1_SCF_qPCR_Rep3.pdf")

with_PRR5_cq <- with_PRR5 %>% 
    mutate(
        rep = rep(c(1, 2), length(with_PRR5$Sample) / 2)
    ) %>%
    pivot_wider(names_from = Target, values_from = Cq) %>%
    mutate(deltaCq = PRR5 - SCF_UBQ)

with_PRR5_controlList <- rep(with_PRR5_cq$deltaCq[1:2], length(with_PRR5_cq$deltaCq) / 2)
with_PRR5_quant <- with_PRR5_cq %>%
    mutate(
        deDeltaCq = deltaCq - with_PRR5_controlList,
        quant = 2^(-deDeltaCq)
    ) %>%
    group_by(Sample, strain, day, phase) %>%
    summarise(means = mean(quant, na.rm = TRUE), sd = sd(quant, na.rm = TRUE)) %>%
    mutate(timeString = str_c(day, "Dpa", phase))
with_PRR5_quant %>% 
    ggplot(aes(x = timeString, y = means, color = strain)) +
    geom_point() +
    geom_errorbar(aes(ymin = means - sd, ymax = means + sd), width = 0.2) +
    geom_line(aes(group = strain)) +
    labs(x = "Time", y = "Relative Expression", title = "PRR5 Relative Expression PRR5/SCF") +
    theme_bw()
ggsave("figure/qPCR_validation/PRR5_SCF_qPCR_Rep3.pdf")

# 新一批qPCR 1个周期，每个时间点2次重复  
# FL_1_N 引物：PHD1 PHD2 SFCU
# FL_0_N 引物：CCA1 PRR9 SFCU

# 读取数据
FL_1_0_N_path <- "qPCR_rawData/PHD1_PHD2_CCA1_PRR9_SCFU_FL0AND1PERIOD.xlsx"
FL_1_0_N_tbl <- readxl::read_excel(FL_1_0_N_path) %>% 
    mutate(Sample = toupper(Sample)) %>%
    mutate(
        strain = toupper(str_split_i(Sample, "_", 1)),
        day = str_split_i(Sample, "_", 2),
        phase = str_split_i(Sample, "_", 3),
    )

# 拆分两个周期
FL_1_N <- FL_1_0_N_tbl %>%
        filter(day == 1) %>% 
        select(Sample, Target, Cq, strain, day, phase)

FL_0_N <- FL_1_0_N_tbl %>%
        filter(day == 0) %>% 
        select(Sample, Target, Cq, strain, day, phase)

FL_0_N_cq <- FL_0_N %>% 
    mutate(
        rep = rep(c(1, 2), length(FL_0_N$Sample) / 2)
    ) %>%
    pivot_wider(names_from = Target, values_from = Cq) %>%
    mutate(CCA1_deltaCq = CCA1 - SCF_UBQ, PRR9_deltaCq = PRR9 - SCF_UBQ)

FL_0_N_controlList <- rep(FL_0_N_cq$CCA1_deltaCq[1:2], length(FL_0_N_cq$CCA1_deltaCq) / 2)
FL_0_N_quant <- FL_0_N_cq %>%
    mutate(
        CCA1_deDeltaCq = CCA1_deltaCq - FL_0_N_controlList,
        CCA1_quant = 2^(-CCA1_deDeltaCq),
        PRR9_deDeltaCq = PRR9_deltaCq - FL_0_N_controlList,
        PRR9_quant = 2^(-PRR9_deDeltaCq)
    ) %>%
    group_by(Sample, strain, day, phase) %>%
    summarise(CCA1_mean = mean(CCA1_quant, na.rm = TRUE),
     CCA1_sd = sd(CCA1_quant, na.rm = TRUE),
     PRR9_mean = mean(PRR9_quant, na.rm = TRUE),
     PRR9_sd = sd(PRR9_quant, na.rm = TRUE)) %>%
    mutate(timeString = str_c(day, "Dpa", phase))

FL_0_N_quant %>%
    ggplot(aes(x = timeString, y = CCA1_mean, color = strain)) +
    geom_point() +
    geom_errorbar(aes(ymin = CCA1_mean - CCA1_sd, ymax = CCA1_mean + CCA1_sd), width = 0.2) +
    geom_line(aes(group = strain)) +
    labs(x = "Time", y = "Relative Expression", title = "CCA1 Relative Expression CCA1/SCF") +
    theme_bw()

ggsave("figure/qPCR_validation/CCA1_SCF_qPCR_Period0_1point_NA.pdf")

FL_0_N_quant %>% 
    ggplot(aes(x = timeString, y = PRR9_mean, color = strain)) +
    geom_point() +
    geom_errorbar(aes(ymin = PRR9_mean - PRR9_sd, ymax = PRR9_mean + PRR9_sd), width = 0.2) +
    geom_line(aes(group = strain)) +
    labs(x = "Time", y = "Relative Expression", title = "PRR9 Relative Expression PRR9/SCF") +
    theme_bw()

ggsave("figure/qPCR_validation/PRR9_SCF_qPCR_Period0_1point_NA.pdf")

## FL_1_N
FL_1_N_cq <- FL_1_N %>% 
    mutate(
        rep = rep(c(1, 2), length(FL_1_N$Sample) / 2)
    ) %>%
    pivot_wider(names_from = Target, values_from = Cq) %>%
    mutate(PHD1_deltaCq = PHD1 - SCF_UBQ, PHD2_deltaCq = PHD2 - SCF_UBQ)

FL_1_N_controlList <- rep(FL_1_N_cq$PHD1_deltaCq[3:4], length(FL_1_N_cq$PHD1_deltaCq) / 2) # 因为第一个时间点数据缺失
FL_1_N_quant <- FL_1_N_cq %>%
    mutate(
        PHD1_deDeltaCq = PHD1_deltaCq - FL_1_N_controlList,
        PHD1_quant = 2^(-PHD1_deDeltaCq),
        PHD2_deDeltaCq = PHD2_deltaCq - FL_1_N_controlList,
        PHD2_quant = 2^(-PHD2_deDeltaCq)
    ) %>%
    group_by(Sample, strain, day, phase) %>%
    summarise(PHD1_mean = mean(PHD1_quant, na.rm = TRUE),
     PHD1_sd = sd(PHD1_quant, na.rm = TRUE),
     PHD2_mean = mean(PHD2_quant, na.rm = TRUE),
     PHD2_sd = sd(PHD2_quant, na.rm = TRUE)) %>%
    mutate(timeString = str_c(day, "Dpa", phase))

FL_1_N_quant %>%
    ggplot(aes(x = timeString, y = PHD1_mean, color = strain)) +
    geom_point() +
    geom_errorbar(aes(ymin = PHD1_mean - PHD1_sd, ymax = PHD1_mean + PHD1_sd), width = 0.2) +
    geom_line(aes(group = strain)) +
    labs(x = "Time", y = "Relative Expression", title = "PHD1 Relative Expression PHD1/SCF") +
    theme_bw()

ggsave("figure/qPCR_validation/PHD1_SCF_qPCR_Period1_1point_NA.pdf")

FL_1_N_quant %>%
    ggplot(aes(x = timeString, y = PHD2_mean, color = strain)) +
    geom_point() +
    geom_errorbar(aes(ymin = PHD2_mean - PHD2_sd, ymax = PHD2_mean + PHD2_sd), width = 0.2) +
    geom_line(aes(group = strain)) +
    labs(x = "Time", y = "Relative Expression", title = "PHD2 Relative Expression PHD2/SCF") +
    theme_bw()

ggsave("figure/qPCR_validation/PHD2_SCF_qPCR_Period1_1point_NA.pdf")

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

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
    mutate(deDeltaCq = deltaCq - controlList, 
           quant = 2^(-deDeltaCq)) %>% 
           group_by(Sample, strain, day, phase) %>% 
           summarise(means = mean(quant), sd = sd(quant)) %>% 
           mutate(timeString = str_c(day, "Dpa", phase))
# head(pRR5_qPCR_tbl_quant)


pRR5_qPCR_tbl_quant %>% ggplot(aes(timeString, means, color = strain)) +
geom_point() +
geom_errorbar(aes(ymin = means - sd, ymax = means + sd), width = 0.2) +
geom_line(aes(group = strain))+
labs(x = "Time", y = "Relative Expression", title = "PRR5 qPCR validation")+
theme_bw()

ggsave("figure/qPCR_validation/PRR5_2_5.pdf")

library(patchwork)
# 分别绘制生物学重复，不区分WT，FL
metaInfo_WT_FL_0_2day_TMM_sample_exp <- WT_FL_0_2day_TMM_sample_exp %>% select(1:6)
plotRepCirca <- function(str, alia_name = "") {
    df <- cbind(metaInfo_WT_FL_0_2day_TMM_sample_exp, measure = WT_FL_0_2day_TMM_sample_exp[str])
    colnames(df)[7] <- "measure"


    rects <-
        data.frame(
            xstart = c(14.5, 38.5, 62.5),
            xend = c(23, 47, 71)
        )
    ggplot(data = df, aes(time, measure)) +
        geom_point(aes(col = strain)) +
        geom_smooth(aes(group = interaction(as.factor(replicate), strain), color = strain), span = 0.3) +
        facet_wrap(~replicate, nrow = 2) +
        # labels = df$labs[1:(length(df$labs) / 4)]
        scale_x_continuous(breaks = seq(1, 69, 4), ) +
        geom_rect(
            data = rects,
            aes(
                xmin = xstart,
                xmax = xend,
                ymin = 0,
                ymax = Inf
            ),
            inherit.aes = FALSE,
            alpha = 0.2
        ) +
        labs(title = str, subtitle = alia_name) +
        ylab("relative expression/(TMM)") +
        theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
}

# 合并生物学重复，区分WT FL

WT_FL_0_2day_rep_reduce <- function(st) {
    tempvar <- sym(st)
    WT_FL_0_2day_TMM_sample_exp %>%
        select(c(1:6), !!tempvar, ) %>%
        group_by(strain, period, time, labs) %>%
        summarise(
            "mean_li" = mean(!!tempvar),
            "std_li" = sd(!!tempvar),
            .groups = "keep"
        )
}

my_cir_plot_WT_FL_0_2dpa <- function(geneid, alia_name = NULL) {
    tryCatch(
        error = function(cnd) {
            str_c("No RNA-seq data WT_FL_0_2dpa available! ", geneid)
        },
        {
            # 光照时间为 7am-10:30pm
            rects <-
                data.frame(
                    xstart = c(14.5, 38.5, 62.5),
                    xend = c(23, 47, 71)
                )
            sub_df <- WT_FL_0_2day_rep_reduce(geneid)
            ggplot(sub_df, aes(x = time, y = mean_li, color = strain)) +
                geom_point(size = 3) +
                facet_grid(rows = vars(strain)) +
                geom_errorbar(
                    aes(
                        ymax = std_li + mean_li,
                        ymin = mean_li - std_li
                    ),
                    width = 2,
                    color = "black"
                ) +
                geom_line() +
                ylab("relative expression/(TMM)") +
                scale_x_continuous(breaks = seq(1, 69, 4), labels = sub_df$labs[1:(length(sub_df$labs) /
                    2)]) +
                geom_rect(
                    data = rects,
                    aes(
                        xmin = xstart,
                        xmax = xend,
                        ymin = 0,
                        ymax = Inf
                    ),
                    inherit.aes = FALSE,
                    alpha = 0.2
                ) +
                labs(title = geneid, subtitle = alia_name) +
                theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
        }
    )
}

# plot "Ghir_D12G020190" "Ghir_A12G005280" "Ghir_A13G013150" "Ghir_D07G006230" "Ghir_D02G019940" "Ghir_D08G002570"
p1 <- plotRepCirca("Ghir_D12G020190")
p2 <- plotRepCirca("Ghir_A12G005280", "AT5G26830.1 / The mRNA is cell-to-cell mobile")
p3 <- plotRepCirca("Ghir_A13G013150", "AT5G61170")
p4 <- plotRepCirca("Ghir_D07G006230", "AT4G11630")
p5 <- plotRepCirca("Ghir_D02G019940", "AT2G43090 / The mRNA is cell-to-cell mobile")
p6 <- plotRepCirca("Ghir_D08G002570", "AT5G57330")

p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(guides = 'collect')
ggsave("figure/circ_WT_FL_TMM_FL_specific_mtx_kmeans_seed1_cluster1_Genes_6gene_test.pdf")

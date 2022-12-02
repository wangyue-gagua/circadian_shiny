library(patchwork)
# 分别绘制生物学重复，不区分WT，FL
metaInfo_WT_FL_0_2day_TMM_sample_exp <- WT_FL_0_2day_TMM_sample_exp %>% select(1:6)
plotRepCirca <- function(str, alia_name = NULL) {
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
        dplyr::select(c(1:6), !!tempvar, ) %>%
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

# grid plot
my_grid_plot_WT_FL_0_2dpa <- function(a_chr_list, alias = NULL) {
    ## input a id list and a gene name, output a grid plot
    gobs <- map(a_chr_list, my_cir_plot_WT_FL_0_2dpa, alias)
    # new_gobs <- map(gobs[seq_along(gobs)], `+`, theme(axis.title.x = element_blank(), axis.text.x = element_blank()))
    # new_gobs[-1] <- gobs[-1]

    gridExtra::grid.arrange(grobs = gobs)
}

metaInfo_salmon_WT_FL_0_2day_TMM_sample_exp <- salmon_WT_FL_0_2day_TMM_sample_exp %>% select(1:6)

salmonPlotRepCirca <- function(str, alia_name = "") {
    df <- cbind(metaInfo_salmon_WT_FL_0_2day_TMM_sample_exp, measure = salmon_WT_FL_0_2day_TMM_sample_exp[str])
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
        scale_x_continuous(breaks = seq(1, 69, 4)) +
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
        ylab("expression level (TMM)") +
        theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
}
library(grid)
salmon_grid_plot_WT_FL_0_2dpa <- function(a_chr_list, alias = NULL) {
    ## input a id list and a gene name, output a grid plot
    gobs <- map(a_chr_list, salmonPlotRepCirca, alias)
    new_gobs <- map(gobs[seq_along(gobs)], `+`, theme(axis.title.x = element_blank(), axis.text.x = element_blank()))
    new_gobs[-1] <- gobs[-1]

    gridExtra::grid.arrange(grobs = new_gobs)
}


plotRepCircaDetrended <- function(str, alia_name = NULL) {
    metaInfo_WT_FL_0_2day_TMM_sample_exp_detrend <- WT_FL_0_2day_TMM_sample_exp_detrended %>% select(1:6)
    df <- cbind(metaInfo_WT_FL_0_2day_TMM_sample_exp_detrend, measure = WT_FL_0_2day_TMM_sample_exp_detrended[str])
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
        ylab("Detrended relative expression") +
        theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
}

# 火山图
deResVolcano <- function(deRes, titles) {
    quantiles <- quantile(deRes$log2FoldChange, c(0.0001, 0.5, 0.9999))
    maxPerdown <- quantiles[[1]]
    maxPerup <- quantiles[[3]]
    numberDown <- sum(deRes$change == "down")
    numberUp <- sum(deRes$change == "up")
    deRes %>%
        mutate(textture = ifelse(log2FoldChange <= maxPerdown | log2FoldChange > maxPerup, Geneid, "")) %>%
        ggplot(aes(x = log2FoldChange, y = -log10(padj), color = change)) +
        geom_point(alpha = 0.4) +
        scale_color_manual(values = c("up" = "#ff4757", "down" = "#546de5", "ns" = "#d2dae2")) +
        geom_vline(xintercept = c(-1, 1), lty = 4, col = "black", lwd = 0.8) +
        geom_hline(yintercept = -log10(padjCutOff), lty = 4, col = "black", lwd = 0.8) +
        geom_label_repel(aes(label = textture),
            size = 3, box.padding = unit(0.5, "lines"),
            point.padding = unit(0.8, "lines"),
            segment.color = "black",
            show.legend = FALSE
        ) +
        labs(title = titles, subtitle = str_c("Down: ", numberDown, " Up: ", numberUp)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

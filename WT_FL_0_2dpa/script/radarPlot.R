library(fmsb)
WT_meta2d <- read_csv("WT_meta2d/JTKresult_WT_0_2day_genes_TMM_EXPR_mergeRep_selected.csv")
FL_meta2d <- read_csv("FL_meta2d/JTKresult_FL_0_2day_genes_TMM_EXPR_mergeRep_selected.csv")

WT_dat <- WT_meta2d %>% filter(BH.Q < 0.01)
FL_dat <- FL_meta2d %>% filter(BH.Q < 0.01)

# The default radar chart
radarchart(data, maxmin = F)
head(WT_dat)
# 光照时间为 7am-10:30pm
WT_dat %>%
    group_by(LAG) %>%
    summarise(n = n()) %>%
    mutate(time = (LAG + 9) %% 24) %>%
    ggplot(aes(x = time, y = n)) +
    geom_bar(stat = "identity", fill = "#3366ff") +
    geom_rect(
        data = data.frame(
            xstart = c(0, 22.5),
            xend = c(7, 24)
        ),
        aes(
            xmin = xstart,
            xmax = xend,
            ymin = 0,
            ymax = Inf
        ),
        inherit.aes = FALSE,
        alpha = 0.2
    ) +
    theme_bw() +
    labs(x = "Time (h) ZT (real)", y = "Number of genes", title = "WT 0-2dpa") +
    scale_x_continuous(breaks = seq(0, 24, 2)) +
    coord_polar()

ggsave("figure/phaseDistribute/WT_0_2dpa_phaseDistribute.pdf")

FL_dat %>%
    group_by(LAG) %>%
    summarise(n = n()) %>%
    mutate(time = (LAG + 9) %% 24) %>%
    ggplot(aes(x = time, y = n)) +
    geom_bar(stat = "identity", fill = "#3366ff") +
    geom_rect(
        data = data.frame(
            xstart = c(0, 22.5),
            xend = c(7, 24)
        ),
        aes(
            xmin = xstart,
            xmax = xend,
            ymin = 0,
            ymax = Inf
        ),
        inherit.aes = FALSE,
        alpha = 0.2
    ) +
    theme_bw() +
    labs(x = "Time (h) ZT (real)", y = "Number of genes", title = "FL 0-2dpa") +
    scale_x_continuous(breaks = seq(0, 24, 2)) +
    coord_polar()
ggsave("figure/phaseDistribute/FL_0_2dpa_phaseDistribute.pdf")
times <- WT_FL_0_2day_TMM_sample_exp$time %>%
    sort() %>%
    unique()
ggplot(
    tibble(
        time = times - 1, period = c(rep(0, 6), rep(1, 6), rep(2, 6)),
        labels = rep(c("9am", "1pm", "5pm", "9pm", "1am", "5am"), 3)
    ),
    aes(x = time %% 24, y = 0),
) +
    geom_bar(stat = "identity") +
    geom_rect(
        # 光照时间为 7am-10:30pm
        data = tibble(
            xstart = c(0, 21.5),
            xend = c(6, 24)
        ),
        aes(
            xmin = xstart,
            xmax = xend,
            ymin = 0,
            ymax = Inf
        ),
        inherit.aes = FALSE,
        alpha = 0.2
    ) +
    geom_point(aes(y = period + 1, size = factor((period + 1)), ), position = position_nudge(x = -0.1)) +
    geom_point(aes(y = period + 1, size = factor((period + 1)), ), position = position_nudge(x = 0.1)) +
    scale_x_continuous(limits = c(-0.1, 24), breaks = seq(0, 20, 4), labels = c("1am", "5am", "9am", "1pm", "5pm", "9pm")) +
    coord_polar() +
    labs(x = "Time (h) ZT", y = "Period", title = "WT FL 0-2dpa sampling", size = "Period") +
    theme_bw() +
    theme(panel.border = element_blank())

ggsave("figure/phaseDistribute/WT_FL_0_2dpa_sampling.pdf")


## circadian event phase distribution
eventDistributionPlot <- function(sample, event) {
    event_sample <- stringr::str_c(event, "_", sample)
    filePath <- stringr::str_c("METACYCLE/rMATs_merge_", event_sample, "_META2d/meta2d_rMATs_merge_", event_sample, "_exp.csv")
    df <- read_csv(filePath)
    df %>%
        filter(meta2d_BH.Q < 0.01) %>%
        group_by(meta2d_phase) %>%
        summarise(n = n()) %>%
        mutate(time = (meta2d_phase + 9) %% 24) %>%
        ggplot(aes(x = time)) +
        geom_histogram(binwidth = 1, fill = "#3366ff") +
        geom_rect(
            data = data.frame(
                xstart = c(0, 22.5),
                xend = c(7, 24)
            ),
            aes(
                xmin = xstart,
                xmax = xend,
                ymin = 0,
                ymax = Inf
            ),
            inherit.aes = FALSE,
            alpha = 0.2
        ) +
        theme_bw() +
        labs(x = "Time (h) ZT (real)", y = "Number of events", title = str_c(sample, " ", event, " Phase distribution")) +
        scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 2), expand = expansion(add = 0)) +
        coord_polar()
}
## A3SS WT
distri_WT_A3SS <- eventDistributionPlot("WT", "A3SS")
distri_FL_A3SS <- eventDistributionPlot("FL", "A3SS")
distri_WT_A5SS <- eventDistributionPlot("WT", "A5SS")
distri_FL_A5SS <- eventDistributionPlot("FL", "A5SS")
distri_WT_MXE <- eventDistributionPlot("WT", "MXE")
distri_FL_MXE <- eventDistributionPlot("FL", "MXE")
distri_WT_RI <- eventDistributionPlot("WT", "RI")
distri_FL_RI <- eventDistributionPlot("FL", "RI")
distri_WT_SE <- eventDistributionPlot("WT", "SE")
distri_FL_SE <- eventDistributionPlot("FL", "SE")
library(patchwork)
(distri_WT_A3SS / distri_WT_A5SS / distri_WT_RI / distri_WT_SE) | (distri_FL_A3SS / distri_FL_A5SS / distri_FL_RI / distri_FL_SE)

allEventDistributionPlot <- function() {
    A3SS_WT <- read_csv("METACYCLE/rMATs_merge_A3SS_WT_META2d/meta2d_rMATs_merge_A3SS_WT_exp.csv") %>%
        filter(meta2d_BH.Q < 0.01) %>%
        group_by(meta2d_phase) %>%
        summarise(n = n()) %>%
        mutate(event = "A3SS", sample = "WT")
    A3SS_FL <- read_csv("METACYCLE/rMATs_merge_A3SS_FL_META2d/meta2d_rMATs_merge_A3SS_FL_exp.csv") %>%
        filter(meta2d_BH.Q < 0.01) %>%
        group_by(meta2d_phase) %>%
        summarise(n = n()) %>%
        mutate(event = "A3SS", sample = "FL")
    A5SS_WT <- read_csv("METACYCLE/rMATs_merge_A5SS_WT_META2d/meta2d_rMATs_merge_A5SS_WT_exp.csv") %>%
        filter(meta2d_BH.Q < 0.01) %>%
        group_by(meta2d_phase) %>%
        summarise(n = n()) %>%
        mutate(event = "A5SS", sample = "WT")
    A5SS_FL <- read_csv("METACYCLE/rMATs_merge_A5SS_FL_META2d/meta2d_rMATs_merge_A5SS_FL_exp.csv") %>%
        filter(meta2d_BH.Q < 0.01) %>%
        group_by(meta2d_phase) %>%
        summarise(n = n()) %>%
        mutate(event = "A5SS", sample = "FL")
    RI_WT <- read_csv("METACYCLE/rMATs_merge_RI_WT_META2d/meta2d_rMATs_merge_RI_WT_exp.csv") %>%
        filter(meta2d_BH.Q < 0.01) %>%
        group_by(meta2d_phase) %>%
        summarise(n = n()) %>%
        mutate(event = "RI", sample = "WT")
    RI_FL <- read_csv("METACYCLE/rMATs_merge_RI_FL_META2d/meta2d_rMATs_merge_RI_FL_exp.csv") %>%
        filter(meta2d_BH.Q < 0.01) %>%
        group_by(meta2d_phase) %>%
        summarise(n = n()) %>%
        mutate(event = "RI", sample = "FL")
    SE_WT <- read_csv("METACYCLE/rMATs_merge_SE_WT_META2d/meta2d_rMATs_merge_SE_WT_exp.csv") %>%
        filter(meta2d_BH.Q < 0.01) %>%
        group_by(meta2d_phase) %>%
        summarise(n = n()) %>%
        mutate(event = "SE", sample = "WT")
    SE_FL <- read_csv("METACYCLE/rMATs_merge_SE_FL_META2d/meta2d_rMATs_merge_SE_FL_exp.csv") %>%
        filter(meta2d_BH.Q < 0.01) %>%
        group_by(meta2d_phase) %>%
        summarise(n = n()) %>%
        mutate(event = "SE", sample = "FL")
    bind_rows(A3SS_WT, A3SS_FL, A5SS_WT, A5SS_FL, RI_WT, RI_FL, SE_WT, SE_FL) %>%
        mutate(time = (meta2d_phase + 9) %% 24) %>%
        ggplot(aes(x = time)) +
        geom_histogram(binwidth = 1, fill = "#3366ff") +
        geom_rect(
            data = data.frame(
                xstart = c(0, 22.5),
                xend = c(7, 24)
            ),
            aes(
                xmin = xstart,
                xmax = xend,
                ymin = 0,
                ymax = Inf
            ),
            inherit.aes = FALSE,
            alpha = 0.2
        ) +
        theme_bw() +
        scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 2), expand = expansion(add = 0)) +
        facet_grid(rows = vars(sample), cols = vars(event)) +
        labs(title = "Rhythmic alternative event distribution") +
        coord_polar()
}
allEventDistributionPlot()
ggsave("figure/phaseDistribute/allEventDistribution.pdf")

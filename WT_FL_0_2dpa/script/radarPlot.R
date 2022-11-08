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
    labs(x = "Time (h) ZT (real)", y = "Number of genes", title = "WT 0-2dpa")+
    scale_x_continuous(breaks = seq(0, 24, 2)) + 
    coord_polar()

ggsave("figure/phaseDistribute/WT_0_2dpa_phaseDistribute.pdf")

FL_dat  %>% 
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
    labs(x = "Time (h) ZT (real)", y = "Number of genes", title = "FL 0-2dpa")+
    scale_x_continuous(breaks = seq(0, 24, 2)) + 
    coord_polar()
ggsave("figure/phaseDistribute/FL_0_2dpa_phaseDistribute.pdf")

arabExpData <- read_tsv("arabExpData/GSE137732_LD_ZT2_47_transcript_expr.txt")
head(arabExpData)
arabExpData %>% filter(str_detect(t_name, "^AT5G61380"))
arabPlotZT <- function(t_name, alias) {
    testDf <- tibble(
        times = seq(2, 47, 3),
        value = arabExpData %>% filter(t_name == !!t_name) %>% select(-c(1:9)) %>% as_vector()
    )

    rects <-
        data.frame(
            xstart = c(17, 41),
            xend = c(23, 47)
        )
    testDf %>%
        ggplot(aes(x = times, y = value)) +
        geom_point() +
        geom_line() +
        xlab("time (h) ZT") +
        ylab("relative expression/(FPKM)") +
        scale_x_continuous(breaks = seq(2, 47, 3), labels = seq(2, 47, 3)) +
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
        labs(title = t_name, subtitle = alias) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
}

arabPlotZT("AT1G65060.2", "4CL3")
arabPlotZT("AT1G65060.1", "4CL3")
arabPlotZT("AT2G46830.1", "CCA1")
ggsave("figure/arabPlotZT/AT2G46830_1_CCA1.pdf", arabPlotZT("AT2G46830.1", "CCA1"), width = 10, height = 5)
arabPlotZT("AT2G46830.2", "CCA1")
arabPlotZT("AT2G46830.3", "CCA1")
arabIsoformPlotZT <- function(g_name, alias = NULL) {
    plotDf <- arabExpData %>%
        filter(gene_name == !!g_name) %>%
        select(-c("chr", "strand", "start", "end", "num_exons", "gene_id", "length")) %>%
        pivot_longer(
            cols = -c(t_name, gene_name),
            names_to = "time",
        ) %>%
        mutate(timeZT = str_extract(time, "\\d+") %>% as.numeric())

    rects <-
        data.frame(
            xstart = c(17, 41),
            xend = c(23, 47)
        )
    plotDf %>%
        ggplot(aes(x = timeZT, y = value, group = t_name, color = t_name)) +
        geom_point() +
        geom_line() +
        xlab("time (h) ZT") +
        ylab("relative expression/(FPKM)") +
        scale_x_continuous(breaks = seq(2, 47, 3), labels = seq(2, 47, 3)) +
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
        labs(title = g_name, subtitle = alias, color = "isoform") +
        theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
}
arabIsoformPlotZT("CCA1")
ggsave("figure/arabPlotZT/CCA1.pdf", arabIsoformPlotZT("CCA1"), width = 10, height = 5)
arabIsoformPlotZT("APRR5")
ggsave("figure/arabPlotZT/APRR5.pdf", arabIsoformPlotZT("APRR5"), width = 10, height = 5)
arabIsoformPlotZT("APRR7")
ggsave("figure/arabPlotZT/APRR7.pdf", arabIsoformPlotZT("APRR7"), width = 10, height = 5)
arabIsoformPlotZT("APRR1", "TOC1")
ggsave("figure/arabPlotZT/APRR1.pdf", arabIsoformPlotZT("APRR1", "TOC1"), width = 10, height = 5)
arabIsoformPlotZT("APRR9")
ggsave("figure/arabPlotZT/APRR9.pdf", arabIsoformPlotZT("APRR9"), width = 10, height = 5)

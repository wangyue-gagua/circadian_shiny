# 第一次测序与第二次测序检测基因数对比
first_sample_info_exp <- data.table::fread("../merged_counts/sample_info_exp.csv") %>%
    as_tibble()
first_sample_list <- first_sample_info_exp %>%
    group_by(strain) %>%
    arrange(time, .by_group = T) %>%
    pull(sample)
first_count_matrix <- read_tsv("../merged_counts/genes.counts.matrix") %>%
    rename(Geneid = `...1`) %>%
    select(Geneid, !!first_sample_list)

second_sample_info_exp <- data.table::fread("../merged_counts/WT_FL_0_2day_TMM_sample_exp.csv") %>%
    as_tibble()
second_sample_list <- second_sample_info_exp %>%
    group_by(strain) %>%
    arrange(time, .by_group = T) %>%
    pull(sample)


second_count_matrix <- read_tsv("../merged_counts/WT_FL_0_2day_genes.counts.matrix") %>%
    rename(Geneid = `...1`) %>%
    select(Geneid, !!second_sample_list)


first_count_matrix_longer <- first_count_matrix %>%
    pivot_longer(cols = -Geneid, names_to = "sample", values_to = "count") %>%
    mutate(sample = factor(sample, levels = first_sample_list))

second_count_matrix_longer <- second_count_matrix %>%
    pivot_longer(cols = -Geneid, names_to = "sample", values_to = "count") %>%
    mutate(sample = factor(sample, levels = second_sample_list))

stat_box_data <- function(y, upper_limit = log10(max(first_count_matrix_longer$count)) * 1.15) {
    return(
        data.frame(
            y = 0.95 * upper_limit,
            label = length(y[10^y > 3])
        )
    )
}
first_count_matrix_longer_pdf <- first_count_matrix_longer %>% ggplot(aes(x = sample, y = log10(count))) +
    geom_boxplot() +
    stat_summary(
        fun.data = stat_box_data,
        geom = "text",
        size = 2,
        color = "red",
        hjust = 0.5,
        vjust = 0.9
    ) +
    labs(title = "First sequencing Gene Counts") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("figure/geneExpDiff/first_count_matrix_longer.pdf", first_count_matrix_longer_pdf, width = 10, height = 4)

second_stat_box_data <- function(y, upper_limit = log10(max(second_count_matrix_longer$count)) * 1.15) {
    return(
        data.frame(
            y = 0.95 * upper_limit,
            label = length(y[10^y > 3])
        )
    )
}
second_count_matrix_longer_pdf <- second_count_matrix_longer %>% ggplot(aes(x = sample, y = log10(count))) +
    geom_boxplot() +
    stat_summary(
        fun.data = second_stat_box_data,
        geom = "text",
        size = 2,
        color = "red",
        hjust = 0.5,
        vjust = 0.9
    ) +
    labs(title = "Second sequencing Gene Counts") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("figure/geneExpDiff/second_count_matrix_longer.pdf", second_count_matrix_longer_pdf, width = 20, height = 4)


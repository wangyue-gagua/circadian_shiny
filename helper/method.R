## reduce replication ,set mean of replication as exp value, a new function
rep_reduce <- function(st) {
  tempvar <- sym(st)
  sample_info_exp %>% select(c(1:6), !!tempvar,) %>% group_by(strain, period, time, labs) %>%
    summarise(
      'mean_li' = mean(!!tempvar),
      'std_li' = sd(!!tempvar),
      .groups = "keep"
    )
}
## the plot function
my_cir_plot <- function(geneid, alia_name = NULL) {
  tryCatch(
    error = function(cnd) {
      str_c("No RNA-seq data available! ", geneid)
    },
    {
      rects <-
        data.frame(xstart = seq(-36, 36, 24),
                   xend = seq(-24, 48, 24))
      sub_df <- rep_reduce(geneid)
      ggplot(sub_df, aes(x = time, y = mean_li, color = strain)) +
        geom_point(size = 3) +
        # scale_x_continuous(breaks = seq(-72,72,24))+
        facet_grid(rows = vars(strain)) +
        geom_errorbar(
          aes(ymax = std_li + mean_li,
              ymin = mean_li - std_li),
          width = 2,
          color = 'black'
        ) +
        geom_line() +
        ylab("relative expression/(TMM)") +
        scale_x_continuous(breaks = seq(-48, 48, 12), labels = sub_df$labs[1:(length(sub_df$labs) /
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
        labs(title = geneid, subtitle = alia_name)
      
    }
    
  )
}

## WT_FL_0_2day
WT_FL_0_2day_rep_reduce <- function(st) {
  tempvar <- sym(st)
  WT_FL_0_2day_TMM_sample_exp %>%
    select(c(1:6), !!tempvar,) %>%
    group_by(strain, period, time, labs) %>%
    summarise(
      'mean_li' = mean(!!tempvar),
      'std_li' = sd(!!tempvar),
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
        data.frame(xstart = c(14.5, 38.5, 62.5),
                   xend = c(23, 47, 71))
      sub_df <- WT_FL_0_2day_rep_reduce(geneid)
      ggplot(sub_df, aes(x = time, y = mean_li, color = strain)) +
        geom_point(size = 3) +
        facet_grid(rows = vars(strain)) +
        geom_errorbar(
          aes(ymax = std_li + mean_li,
              ymin = mean_li - std_li),
          width = 2,
          color = 'black'
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


my_tissue_plot <- function(st) {
  tryCatch(
    error = function(cnd) {
      str_c("No RNA-seq data available! ", st)
    },
    {
      temp <-
        iso_exp_tpm %>% filter(str_match(tracking_id, "(.*)\\.\\d")[,2] %in% st) %>% column_to_rownames(var = "tracking_id") %>%
        mutate("leaf_0_" = leaf_0) %>% select(!leaf_0)
      if (nrow(temp) > 0) {
        temp <- temp[1,]
      } else{
        stop("No data available!")
      }
      df <-
        tibble(
          "samples" = factor(names(temp), levels = names(temp)),
          "fpkm" = as.numeric(temp),
          "tissue" = map_chr(names(temp), ~ str_match(.x, "(\\w+?)_.*")[[2]])
        )
      df %>% ggplot(aes(x = samples, y = log(fpkm + 1))) +
        geom_point(aes(color = tissue), size = 3) +
        geom_path(aes(
          group = tissue,
          y = log(fpkm + 1),
          color = tissue
        )) +
        theme(axis.text.x = element_text(angle = 90)) +
        ylab("relative expression/log(fpkm+1)") +
        labs(title = st)
    }
  )
}

my_prot_plot <- function(st) {
  tryCatch(
    error = function(cnd) {
      str_c("No MS data available! ", st)
    },
    
    {
      df <- AD_pro_19920 %>% filter(str_detect(tracking_id, st))
      if (nrow(df) > 0) {
        df <- df[1,]
      }
      df_pl <-
        tibble(
          samples = factor(
            df %>% select(contains("WT")) %>% names(),
            levels = (df %>% select(contains("WT")) %>% names())
          ),
          exps = (df %>% select(contains("WT")) %>% as_vector())
        )
      ggplot(data = df_pl, aes(
        x = samples,
        y = exps,
        color = str_detect(samples, "tmt")
      )) +
        geom_point(size = 3) +
        labs(
          title = st,
          subtitle = str_c("phosphorylation site: ", df$sites),
          col = "treatment"
        ) +
        theme(legend.text = element_blank())
    }
  )
}

my_go_table <- function(str) {
  str_gtl <- genes2Go %>% filter(str_detect(Query, pattern = str))
  if (nrow(str_gtl) == 0)
    map_dfr(str_gtl, function(x)
      x <- "No GO result")
  else
    return(str_gtl)
}

my_ara_homo_tbl <- function(str) {
  str_homo <-
    blastp_AD1_HAU_v1_0_vs_arabidopsis_1 %>% filter(str_detect(Query, pattern = str))
  if (nrow(str_homo) == 0)
    map_dfr(str_homo, function(x)
      x <- "No homolog in tair10")
  else
    return(str_homo)
}

## ID convert
IdConvert <- function(convert_bed, input_ID) {
  return(
    convert_bed %>% filter(str_detect(ID_In, input_ID))
  )
}
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
tablePath <- args[1] # 需要手动额外添加一列列名
outPath <- args[2]
WT_FL_0_2day_genes_TMM_EXPR <- read_delim(tablePath, 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)


WT_FL_0_2day_genes_TMM_EXPR_longer <- t(WT_FL_0_2day_genes_TMM_EXPR %>% as_tibble())
colnames(WT_FL_0_2day_genes_TMM_EXPR_longer) <- WT_FL_0_2day_genes_TMM_EXPR_longer[1,]
WT_FL_0_2day_genes_TMM_EXPR_longer_trans <- WT_FL_0_2day_genes_TMM_EXPR_longer[2:nrow(WT_FL_0_2day_genes_TMM_EXPR_longer),]
trans <- function(sample) {
  if(str_detect(sample, "\\dDPA_1am")) {
    return(17 + parse_integer(str_match(sample, "(\\d)DPA")[[2]]) * 24)
  }else if (str_detect(sample, "\\dDPA_5am")) {
    return(21 + parse_integer(str_match(sample, "(\\d)DPA")[[2]]) * 24)
  }else if (str_detect(sample, "\\dDPA_9am")) {
    return(1 + parse_integer(str_match(sample, "(\\d)DPA")[[2]]) * 24)
  }else if (str_detect(sample, "\\dDPA_1pm")) {
    return(5 + parse_integer(str_match(sample, "(\\d)DPA")[[2]]) * 24)
  }else if (str_detect(sample, "\\dDPA_5pm")) {
    return(9 + parse_integer(str_match(sample, "(\\d)DPA")[[2]]) * 24)
  }else if (str_detect(sample, "\\dDPA_9pm")) {
    return(13 + parse_integer(str_match(sample, "(\\d)DPA")[[2]]) * 24)
  }
}
WT_FL_0_2day_genes_TMM_EXPR_longer_trans_final <- WT_FL_0_2day_genes_TMM_EXPR_longer_trans %>%
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "sample") %>% 
  mutate(strain=map_chr(sample, function(str) str_extract(str,"(FL)|(WT)")),
         period=map_chr(sample, function(str) str_match(str, ".*?_(\\d)DPA_.*")[[2]]),
         time= map_dbl(sample, trans),
         replicate=map_chr(sample, function(str) str_match(str, ".*?_\\dDPA_.*(\\d)")[[2]]),
         labs=map_chr(sample, function(str) str_match(str, ".*?_(\\dDPA_.*)\\d")[[2]]), .after = sample)


write_csv(WT_FL_0_2day_genes_TMM_EXPR_longer_trans_final, outPath)
library(cosinor2)
library(ggplot2)

fit.pa_ext.cosinor <- population.cosinor.lm(data = PA_extraverts, time = PA_time, period = 24)
head(PA_extraverts)
head(PA_time)
nrow(PA_extraverts)
# 不同重复可以用population lm计算
library(circacompare)
set.seed(42)
data_grouped <- make_data(phi1 = 6, noise_sd = 1)
head(data_grouped)

set.seed(42)
phi1_in <- 3
mixed_data <- function(n) {
    counter <- 1
    for (i in 1:n) {
        x <- make_data(k1 = 0, alpha1 = 5, phi1 = rnorm(1, phi1_in, 1), hours = 72, noise_sd = 2)
        x$id <- counter
        counter <- counter + 1
        if (i == 1) {
            res <- x
        } else {
            res <- rbind(res, x)
        }
    }
    return(res)
}
df <- mixed_data(20)
out <- circacompare_mixed(
    x = df,
    col_time = "time",
    col_group = "group",
    col_outcome = "measure",
    col_id = "id",
    control = list(grouped_params = c("alpha", "phi"), random_params = c("phi1")),
    period = 24
)

ggplot(data = df[df$id %in% c(1:6), ], aes(time, measure)) +
    geom_point(aes(col = group)) +
    geom_smooth(aes(group = interaction(as.factor(id), group)), span = 0.3) +
    facet_wrap(~id)





plotRepCirca("Ghir_D11G030190")
plotRepCirca("Ghir_A10G001070")
plotRepCirca("Ghir_D02G008700")
plotRepCirca("Ghir_A05G021450")


plotRepCirca("Ghir_A01G000210")
plotRepCirca("Ghir_D04G000830")
plotRepCirca("Ghir_A10G017040")
ggsave("figure/Ghir_A10G017040.pdf")
plotRepCirca("Ghir_A11G010900", "PRR9")
ggsave("figure/plotRepCirca_Ghir_A11G010900.pdf")





circacompare:::extract_model_coefs(out_Ghir_D12G009360$fit)
out_Ghir_D12G009360$summary %>% filter(parameter == "P-value for mesor difference")
out_Ghir_D12G009360



my_cir_plot_WT_FL_0_2dpa("Ghir_A11G010900", "PRR9")

fit <- ssections(data = PANAS_november, time = PANAS_time, period = 5, interval = 6, increment = 6)

head(PANAS_time)

rownames(PANAS_november)
PANAS_november %>% select(X01, X04)


## run cosinor analysis on the WT FL respectively
### import data
WT_0_2day_genes_TMM_EXPR_mergeRep_selected <- read_csv("WT_0_2day_genes_TMM_EXPR_mergeRep_selected.csv")
FL_0_2day_genes_TMM_EXPR_mergeRep_selected <- read_csv("FL_0_2day_genes_TMM_EXPR_mergeRep_selected.csv")

sampleTime <- seq(1, 69, 4)
WT_0_2day_genes_TMM_EXPR_mergeRep_selected_df <- WT_0_2day_genes_TMM_EXPR_mergeRep_selected %>% column_to_rownames("Geneid")

WT_0_2day_genes_TMM_EXPR_mergeRep_selected_df
filter(WT_0_2day_genes_TMM_EXPR_mergeRep_selected_df, rownames(WT_0_2day_genes_TMM_EXPR_mergeRep_selected_df) == "Ghir_D01G016160")
WT_0_2day_genes_TMM_EXPR_mergeRep_selected_ssections <- ssections(
    data = filter(WT_0_2day_genes_TMM_EXPR_mergeRep_selected_df, rownames(WT_0_2day_genes_TMM_EXPR_mergeRep_selected_df) == "Ghir_D01G016160"),
    time = seq(1, 69, 4), period = 24, interval = 6, increment = 6
)
WT_0_2day_genes_TMM_EXPR_mergeRep_selected_ssections$cosinors[[1]]$single.cos[[1]]$coefficients
head(WT_0_2day_genes_TMM_EXPR_mergeRep_selected_df)
head(PANAS_november)
head(PANAS_time)
sampleTime
PANAS_time
class(PANAS_november)
class(WT_0_2day_genes_TMM_EXPR_mergeRep_selected_df)
seq(1, 18, 1)
p1 <- salmonPlotRepCirca("Ghir_D13G018480.1", "Proline-rich spliceosome-associated (PSP) family protein / zinc knuckle (CCHC-type) family protein")
ggsave(p1, "figure/transcriptionPlot/Ghir_D13G018480.1_PSP.pdf", width = 10, height = 10)

p2 <- salmonPlotRepCirca("Ghir_A08G008200.1", "spliceosomal protein U1A")
p3 <- salmonPlotRepCirca("Ghir_D13G018480.2", "Proline-rich spliceosome-associated (PSP) family protein / zinc knuckle (CCHC-type) family protein")

p4 <- salmonPlotRepCirca("Ghir_D06G014080.1", "proline-rich spliceosome-associated (PSP) family protein")
p5 <- salmonPlotRepCirca("Ghir_A13G017740.1", "proline-rich spliceosome-associated (PSP) family protein")

p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 2, nrow = 3, guides = "collect")
ggsave("figure/transcriptionPlot/PSP_5.pdf", width = 10, height = 10)

# important function
convertToDblAndMean <- function(str) {
    if (is.na(str)) {
        return(0)
    } else {
        a1 <- as.numeric(str_split(str, ",")[[1]][[1]])
        a2 <- as.numeric(str_split(str, ",")[[1]][[2]])
        means <- (a1 + a2) / 2
        return(ifelse(is.na(means), 0, means))
    }
}

# 0dpa SE
rMATs_0DPA_9am_SE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_9am/SE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    select(c(GeneID, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, IncLevel1, IncLevel2))
rMATs_0DPA_1pm_SE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_1pm/SE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    select(c(GeneID, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, IncLevel1, IncLevel2))
rMATs_0DPA_5pm_SE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_5pm/SE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    select(c(GeneID, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, IncLevel1, IncLevel2))
rMATs_0DPA_9pm_SE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_9pm/SE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    select(c(GeneID, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, IncLevel1, IncLevel2))
rMATs_0DPA_1am_SE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_1am/SE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    select(c(GeneID, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, IncLevel1, IncLevel2))
rMATs_0DPA_5am_SE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_5am/SE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    select(c(GeneID, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, IncLevel1, IncLevel2))

# 1dpa SE
rMATs_1DPA_9am_SE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_9am/SE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    select(c(GeneID, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, IncLevel1, IncLevel2))
rMATs_1DPA_1pm_SE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_1pm/SE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    select(c(GeneID, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, IncLevel1, IncLevel2))
rMATs_1DPA_5pm_SE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_5pm/SE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    select(c(GeneID, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, IncLevel1, IncLevel2))
rMATs_1DPA_9pm_SE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_9pm/SE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    select(c(GeneID, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, IncLevel1, IncLevel2))
rMATs_1DPA_1am_SE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_1am/SE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    select(c(GeneID, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, IncLevel1, IncLevel2))
rMATs_1DPA_5am_SE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_5am/SE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    select(c(GeneID, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, IncLevel1, IncLevel2))

# 2dpa SE
rMATs_2DPA_9am_SE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_9am/SE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    select(c(GeneID, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, IncLevel1, IncLevel2))
rMATs_2DPA_1pm_SE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_1pm/SE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    select(c(GeneID, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, IncLevel1, IncLevel2))
rMATs_2DPA_5pm_SE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_5pm/SE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    select(c(GeneID, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, IncLevel1, IncLevel2))
rMATs_2DPA_9pm_SE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_9pm/SE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    select(c(GeneID, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, IncLevel1, IncLevel2))
rMATs_2DPA_1am_SE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_1am/SE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    select(c(GeneID, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, IncLevel1, IncLevel2))
rMATs_2DPA_5am_SE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_5am/SE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    select(c(GeneID, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, IncLevel1, IncLevel2))

# merge SE
rMATs_merge_SE <- full_join(rMATs_0DPA_9am_SE, rMATs_0DPA_1pm_SE,
    by = c("GeneID", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"), suffix = c("0DPA_9am", "0DPA_1pm")
) %>%
    full_join(rMATs_0DPA_5pm_SE,
        by = c("GeneID", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"), suffix = c("0DPA_1pm", "0DPA_5pm")
    ) %>%
    full_join(rMATs_0DPA_9pm_SE,
        by = c("GeneID", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"), suffix = c("0DPA_5pm", "0DPA_9pm")
    ) %>%
    full_join(rMATs_0DPA_1am_SE,
        by = c("GeneID", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"), suffix = c("0DPA_9pm", "0DPA_1am")
    ) %>%
    full_join(rMATs_0DPA_5am_SE,
        by = c("GeneID", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"), suffix = c("0DPA_1am", "0DPA_5am")
    ) %>%
    full_join(rMATs_1DPA_9am_SE,
        by = c("GeneID", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"), suffix = c("0DPA_5am", "1DPA_9am")
    ) %>%
    full_join(rMATs_1DPA_1pm_SE,
        by = c("GeneID", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"), suffix = c("1DPA_9am", "1DPA_1pm")
    ) %>%
    full_join(rMATs_1DPA_5pm_SE,
        by = c("GeneID", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"), suffix = c("1DPA_1pm", "1DPA_5pm")
    ) %>%
    full_join(rMATs_1DPA_9pm_SE,
        by = c("GeneID", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"), suffix = c("1DPA_5pm", "1DPA_9pm")
    ) %>%
    full_join(rMATs_1DPA_1am_SE,
        by = c("GeneID", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"), suffix = c("1DPA_9pm", "1DPA_1am")
    ) %>%
    full_join(rMATs_1DPA_5am_SE,
        by = c("GeneID", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"), suffix = c("1DPA_1am", "1DPA_5am")
    ) %>%
    full_join(rMATs_2DPA_9am_SE,
        by = c("GeneID", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"), suffix = c("1DPA_5am", "2DPA_9am")
    ) %>%
    full_join(rMATs_2DPA_1pm_SE,
        by = c("GeneID", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"), suffix = c("2DPA_9am", "2DPA_1pm")
    ) %>%
    full_join(rMATs_2DPA_5pm_SE,
        by = c("GeneID", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"), suffix = c("2DPA_1pm", "2DPA_5pm")
    ) %>%
    full_join(rMATs_2DPA_9pm_SE,
        by = c("GeneID", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"), suffix = c("2DPA_5pm", "2DPA_9pm")
    ) %>%
    full_join(rMATs_2DPA_1am_SE,
        by = c("GeneID", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"), suffix = c("2DPA_9pm", "2DPA_1am")
    ) %>%
    full_join(rMATs_2DPA_5am_SE,
        by = c("GeneID", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"), suffix = c("2DPA_1am", "2DPA_5am")
    )

## save middle file
# write_csv(rMATs_merge_SE, "mediumDataSave/rMATs_merge_SE.csv")

rMATs_merge_SE_WT_exp <- rMATs_merge_SE %>%
    select(-c(exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE)) %>%
    mutate(
        WT_0DPA_9am = map_dbl(IncLevel10DPA_9am, convertToDblAndMean),
        WT_0DPA_1pm = map_dbl(IncLevel10DPA_1pm, convertToDblAndMean),
        WT_0DPA_5pm = map_dbl(IncLevel10DPA_5pm, convertToDblAndMean),
        WT_0DPA_9pm = map_dbl(IncLevel10DPA_9pm, convertToDblAndMean),
        WT_0DPA_1am = map_dbl(IncLevel10DPA_1am, convertToDblAndMean),
        WT_0DPA_5am = map_dbl(IncLevel10DPA_5am, convertToDblAndMean),
        WT_1DPA_9am = map_dbl(IncLevel11DPA_9am, convertToDblAndMean),
        WT_1DPA_1pm = map_dbl(IncLevel11DPA_1pm, convertToDblAndMean),
        WT_1DPA_5pm = map_dbl(IncLevel11DPA_5pm, convertToDblAndMean),
        WT_1DPA_9pm = map_dbl(IncLevel11DPA_9pm, convertToDblAndMean),
        WT_1DPA_1am = map_dbl(IncLevel11DPA_1am, convertToDblAndMean),
        WT_1DPA_5am = map_dbl(IncLevel11DPA_5am, convertToDblAndMean),
        WT_2DPA_9am = map_dbl(IncLevel12DPA_9am, convertToDblAndMean),
        WT_2DPA_1pm = map_dbl(IncLevel12DPA_1pm, convertToDblAndMean),
        WT_2DPA_5pm = map_dbl(IncLevel12DPA_5pm, convertToDblAndMean),
        WT_2DPA_9pm = map_dbl(IncLevel12DPA_9pm, convertToDblAndMean),
        WT_2DPA_1am = map_dbl(IncLevel12DPA_1am, convertToDblAndMean),
        WT_2DPA_5am = map_dbl(IncLevel12DPA_5am, convertToDblAndMean)
    ) %>%
    select(c(GeneID, starts_with("WT_")))
# write_csv(rMATs_merge_SE_WT_exp, "mediumDataSave/rMATs_merge_SE_WT_exp.csv")
## FL SE
rMATs_merge_SE <- read_csv("mediumDataSave/rMATs_merge_SE.csv")
rMATs_merge_SE_FL_exp <- rMATs_merge_SE %>%
    select(-c(exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE)) %>%
    mutate(
        FL_0DPA_9am = map_dbl(IncLevel20DPA_9am, convertToDblAndMean),
        FL_0DPA_1pm = map_dbl(IncLevel20DPA_1pm, convertToDblAndMean),
        FL_0DPA_5pm = map_dbl(IncLevel20DPA_5pm, convertToDblAndMean),
        FL_0DPA_9pm = map_dbl(IncLevel20DPA_9pm, convertToDblAndMean),
        FL_0DPA_1am = map_dbl(IncLevel20DPA_1am, convertToDblAndMean),
        FL_0DPA_5am = map_dbl(IncLevel20DPA_5am, convertToDblAndMean),
        FL_1DPA_9am = map_dbl(IncLevel21DPA_9am, convertToDblAndMean),
        FL_1DPA_1pm = map_dbl(IncLevel21DPA_1pm, convertToDblAndMean),
        FL_1DPA_5pm = map_dbl(IncLevel21DPA_5pm, convertToDblAndMean),
        FL_1DPA_9pm = map_dbl(IncLevel21DPA_9pm, convertToDblAndMean),
        FL_1DPA_1am = map_dbl(IncLevel21DPA_1am, convertToDblAndMean),
        FL_1DPA_5am = map_dbl(IncLevel21DPA_5am, convertToDblAndMean),
        FL_2DPA_9am = map_dbl(IncLevel22DPA_9am, convertToDblAndMean),
        FL_2DPA_1pm = map_dbl(IncLevel22DPA_1pm, convertToDblAndMean),
        FL_2DPA_5pm = map_dbl(IncLevel22DPA_5pm, convertToDblAndMean),
        FL_2DPA_9pm = map_dbl(IncLevel22DPA_9pm, convertToDblAndMean),
        FL_2DPA_1am = map_dbl(IncLevel22DPA_1am, convertToDblAndMean),
        FL_2DPA_5am = map_dbl(IncLevel22DPA_5am, convertToDblAndMean)
    ) %>%
    select(c(GeneID, starts_with("FL_")))
write_csv(rMATs_merge_SE_FL_exp, "mediumDataSave/rMATs_merge_SE_FL_exp.csv")

# 0DPA A3SS
rMATs_0DPA_9am_A3SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_9am/A3SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_0DPA_9am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_0DPA_9am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_0DPA_9am, FL_0DPA_9am))
rMATs_0DPA_1pm_A3SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_1pm/A3SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_0DPA_1pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_0DPA_1pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_0DPA_1pm, FL_0DPA_1pm))
rMATs_0DPA_5pm_A3SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_5pm/A3SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_0DPA_5pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_0DPA_5pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_0DPA_5pm, FL_0DPA_5pm))
rMATs_0DPA_9pm_A3SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_9pm/A3SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_0DPA_9pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_0DPA_9pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_0DPA_9pm, FL_0DPA_9pm))
rMATs_0DPA_1am_A3SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_1am/A3SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_0DPA_1am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_0DPA_1am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_0DPA_1am, FL_0DPA_1am))
rMATs_0DPA_5am_A3SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_5am/A3SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_0DPA_5am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_0DPA_5am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_0DPA_5am, FL_0DPA_5am))

## 1DPA A3SS
rMATs_1DPA_9am_A3SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_9am/A3SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_1DPA_9am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_1DPA_9am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_1DPA_9am, FL_1DPA_9am))
rMATs_1DPA_1pm_A3SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_1pm/A3SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_1DPA_1pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_1DPA_1pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_1DPA_1pm, FL_1DPA_1pm))
rMATs_1DPA_5pm_A3SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_5pm/A3SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_1DPA_5pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_1DPA_5pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_1DPA_5pm, FL_1DPA_5pm))
rMATs_1DPA_9pm_A3SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_9pm/A3SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_1DPA_9pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_1DPA_9pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_1DPA_9pm, FL_1DPA_9pm))
rMATs_1DPA_1am_A3SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_1am/A3SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_1DPA_1am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_1DPA_1am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_1DPA_1am, FL_1DPA_1am))
rMATs_1DPA_5am_A3SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_5am/A3SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_1DPA_5am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_1DPA_5am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_1DPA_5am, FL_1DPA_5am))

## 2DPA A3SS
rMATs_2DPA_9am_A3SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_9am/A3SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_2DPA_9am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_2DPA_9am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_2DPA_9am, FL_2DPA_9am))
rMATs_2DPA_1pm_A3SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_1pm/A3SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_2DPA_1pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_2DPA_1pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_2DPA_1pm, FL_2DPA_1pm))
rMATs_2DPA_5pm_A3SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_5pm/A3SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_2DPA_5pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_2DPA_5pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_2DPA_5pm, FL_2DPA_5pm))
rMATs_2DPA_9pm_A3SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_9pm/A3SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_2DPA_9pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_2DPA_9pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_2DPA_9pm, FL_2DPA_9pm))
rMATs_2DPA_1am_A3SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_1am/A3SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_2DPA_1am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_2DPA_1am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_2DPA_1am, FL_2DPA_1am))
rMATs_2DPA_5am_A3SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_5am/A3SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_2DPA_5am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_2DPA_5am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_2DPA_5am, FL_2DPA_5am))

## merge
rMATs_merge_A3SS <- rMATs_0DPA_9am_A3SS %>%
    full_join(
        rMATs_0DPA_1pm_A3SS,
        by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")
    ) %>%
    full_join(
        rMATs_0DPA_5pm_A3SS,
        by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")
    ) %>%
    full_join(
        rMATs_0DPA_9pm_A3SS,
        by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")
    ) %>%
    full_join(
        rMATs_0DPA_1am_A3SS,
        by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")
    ) %>%
    full_join(
        rMATs_0DPA_5am_A3SS,
        by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")
    ) %>%
    full_join(
        rMATs_1DPA_9am_A3SS,
        by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")
    ) %>%
    full_join(
        rMATs_1DPA_1pm_A3SS,
        by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")
    ) %>%
    full_join(
        rMATs_1DPA_5pm_A3SS,
        by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")
    ) %>%
    full_join(
        rMATs_1DPA_9pm_A3SS,
        by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")
    ) %>%
    full_join(
        rMATs_1DPA_1am_A3SS,
        by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")
    ) %>%
    full_join(
        rMATs_1DPA_5am_A3SS,
        by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")
    ) %>%
    full_join(
        rMATs_2DPA_9am_A3SS,
        by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")
    ) %>%
    full_join(
        rMATs_2DPA_1pm_A3SS,
        by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")
    ) %>%
    full_join(
        rMATs_2DPA_5pm_A3SS,
        by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")
    ) %>%
    full_join(
        rMATs_2DPA_9pm_A3SS,
        by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")
    ) %>%
    full_join(
        rMATs_2DPA_1am_A3SS,
        by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")
    ) %>%
    full_join(
        rMATs_2DPA_5am_A3SS,
        by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")
    )

# write_csv(rMATs_merge_A3SS, "mediumDataSave/rMATs_merge_A3SS.csv")
rMATs_merge_A3SS_WT_exp <- rMATs_merge_A3SS %>%
    select(c(GeneID, starts_with("WT_")))
rMATs_merge_A3SS_FL_exp <- rMATs_merge_A3SS %>%
    select(c(GeneID, starts_with("FL_")))
# write_csv(rMATs_merge_A3SS_WT_exp, "mediumDataSave/rMATs_merge_A3SS_WT_exp.csv")
# write_csv(rMATs_merge_A3SS_FL_exp, "mediumDataSave/rMATs_merge_A3SS_FL_exp.csv")

# 0DPA A5SS
rMATs_0DPA_9am_A5SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_9am/A5SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_0DPA_9am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_0DPA_9am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_0DPA_9am, FL_0DPA_9am))
rMATs_0DPA_1pm_A5SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_1pm/A5SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_0DPA_1pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_0DPA_1pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_0DPA_1pm, FL_0DPA_1pm))
rMATs_0DPA_5pm_A5SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_5pm/A5SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_0DPA_5pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_0DPA_5pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_0DPA_5pm, FL_0DPA_5pm))
rMATs_0DPA_9pm_A5SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_9pm/A5SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_0DPA_9pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_0DPA_9pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_0DPA_9pm, FL_0DPA_9pm))
rMATs_0DPA_1am_A5SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_1am/A5SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_0DPA_1am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_0DPA_1am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_0DPA_1am, FL_0DPA_1am))
rMATs_0DPA_5am_A5SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_5am/A5SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_0DPA_5am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_0DPA_5am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_0DPA_5am, FL_0DPA_5am))
## 1DPA A5SS
rMATs_1DPA_9am_A5SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_9am/A5SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_1DPA_9am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_1DPA_9am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_1DPA_9am, FL_1DPA_9am))
rMATs_1DPA_1pm_A5SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_1pm/A5SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_1DPA_1pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_1DPA_1pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_1DPA_1pm, FL_1DPA_1pm))
rMATs_1DPA_5pm_A5SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_5pm/A5SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_1DPA_5pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_1DPA_5pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_1DPA_5pm, FL_1DPA_5pm))
rMATs_1DPA_9pm_A5SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_9pm/A5SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_1DPA_9pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_1DPA_9pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_1DPA_9pm, FL_1DPA_9pm))
rMATs_1DPA_1am_A5SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_1am/A5SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_1DPA_1am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_1DPA_1am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_1DPA_1am, FL_1DPA_1am))
rMATs_1DPA_5am_A5SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_5am/A5SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_1DPA_5am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_1DPA_5am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_1DPA_5am, FL_1DPA_5am))
## 2DPA A5SS
rMATs_2DPA_9am_A5SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_9am/A5SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_2DPA_9am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_2DPA_9am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_2DPA_9am, FL_2DPA_9am))
rMATs_2DPA_1pm_A5SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_1pm/A5SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_2DPA_1pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_2DPA_1pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_2DPA_1pm, FL_2DPA_1pm))
rMATs_2DPA_5pm_A5SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_5pm/A5SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_2DPA_5pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_2DPA_5pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_2DPA_5pm, FL_2DPA_5pm))
rMATs_2DPA_9pm_A5SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_9pm/A5SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_2DPA_9pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_2DPA_9pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_2DPA_9pm, FL_2DPA_9pm))
rMATs_2DPA_1am_A5SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_1am/A5SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_2DPA_1am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_2DPA_1am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_2DPA_1am, FL_2DPA_1am))
rMATs_2DPA_5am_A5SS <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_5am/A5SS.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_2DPA_5am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_2DPA_5am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, WT_2DPA_5am, FL_2DPA_5am))
## merge A5SS
rMATs_merge_A5SS <- rMATs_0DPA_9am_A5SS %>%
    full_join(rMATs_0DPA_1pm_A5SS, by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")) %>%
    full_join(rMATs_0DPA_5pm_A5SS, by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")) %>%
    full_join(rMATs_0DPA_9pm_A5SS, by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")) %>%
    full_join(rMATs_0DPA_1am_A5SS, by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")) %>%
    full_join(rMATs_0DPA_5am_A5SS, by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")) %>%
    full_join(rMATs_1DPA_9am_A5SS, by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")) %>%
    full_join(rMATs_1DPA_1pm_A5SS, by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")) %>%
    full_join(rMATs_1DPA_5pm_A5SS, by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")) %>%
    full_join(rMATs_1DPA_9pm_A5SS, by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")) %>%
    full_join(rMATs_1DPA_1am_A5SS, by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")) %>%
    full_join(rMATs_1DPA_5am_A5SS, by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")) %>%
    full_join(rMATs_2DPA_9am_A5SS, by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")) %>%
    full_join(rMATs_2DPA_1pm_A5SS, by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")) %>%
    full_join(rMATs_2DPA_5pm_A5SS, by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")) %>%
    full_join(rMATs_2DPA_9pm_A5SS, by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")) %>%
    full_join(rMATs_2DPA_1am_A5SS, by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")) %>%
    full_join(rMATs_2DPA_5am_A5SS, by = c("GeneID", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"))
## write to file
# write_csv(rMATs_merge_A5SS, "mediumDataSave/rMATs_merge_A5SS.csv")
rMATs_merge_A5SS_WT_exp <- rMATs_merge_A5SS %>%
    select(c(GeneID, starts_with("WT_")))
rMATs_merge_A5SS_FL_exp <- rMATs_merge_A5SS %>%
    select(c(GeneID, starts_with("FL_")))
# write_csv(rMATs_merge_A5SS_WT_exp, "mediumDataSave/rMATs_merge_A5SS_WT_exp.csv")
# write_csv(rMATs_merge_A5SS_FL_exp, "mediumDataSave/rMATs_merge_A5SS_FL_exp.csv")

# MXE
## 0DPA MXE
rMATs_0DPA_9am_MXE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_9am/MXE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_0DPA_9am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_0DPA_9am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, `1stExonStart_0base`, `1stExonEnd`, `2ndExonStart_0base`, `2ndExonEnd`, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_0DPA_9am, FL_0DPA_9am))
rMATs_0DPA_1pm_MXE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_1pm/MXE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_0DPA_1pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_0DPA_1pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, `1stExonStart_0base`, `1stExonEnd`, `2ndExonStart_0base`, `2ndExonEnd`, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_0DPA_1pm, FL_0DPA_1pm))
rMATs_0DPA_5pm_MXE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_5pm/MXE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_0DPA_5pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_0DPA_5pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, `1stExonStart_0base`, `1stExonEnd`, `2ndExonStart_0base`, `2ndExonEnd`, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_0DPA_5pm, FL_0DPA_5pm))
rMATs_0DPA_9pm_MXE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_9pm/MXE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_0DPA_9pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_0DPA_9pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, `1stExonStart_0base`, `1stExonEnd`, `2ndExonStart_0base`, `2ndExonEnd`, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_0DPA_9pm, FL_0DPA_9pm))
rMATs_0DPA_1am_MXE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_1am/MXE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_0DPA_1am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_0DPA_1am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, `1stExonStart_0base`, `1stExonEnd`, `2ndExonStart_0base`, `2ndExonEnd`, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_0DPA_1am, FL_0DPA_1am))
rMATs_0DPA_5am_MXE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_5am/MXE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_0DPA_5am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_0DPA_5am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, `1stExonStart_0base`, `1stExonEnd`, `2ndExonStart_0base`, `2ndExonEnd`, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_0DPA_5am, FL_0DPA_5am))
## 1DPA MXE
rMATs_1DPA_9am_MXE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_9am/MXE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_1DPA_9am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_1DPA_9am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, `1stExonStart_0base`, `1stExonEnd`, `2ndExonStart_0base`, `2ndExonEnd`, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_1DPA_9am, FL_1DPA_9am))
rMATs_1DPA_1pm_MXE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_1pm/MXE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_1DPA_1pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_1DPA_1pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, `1stExonStart_0base`, `1stExonEnd`, `2ndExonStart_0base`, `2ndExonEnd`, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_1DPA_1pm, FL_1DPA_1pm))
rMATs_1DPA_5pm_MXE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_5pm/MXE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_1DPA_5pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_1DPA_5pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, `1stExonStart_0base`, `1stExonEnd`, `2ndExonStart_0base`, `2ndExonEnd`, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_1DPA_5pm, FL_1DPA_5pm))
rMATs_1DPA_9pm_MXE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_9pm/MXE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_1DPA_9pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_1DPA_9pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, `1stExonStart_0base`, `1stExonEnd`, `2ndExonStart_0base`, `2ndExonEnd`, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_1DPA_9pm, FL_1DPA_9pm))
rMATs_1DPA_1am_MXE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_1am/MXE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_1DPA_1am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_1DPA_1am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, `1stExonStart_0base`, `1stExonEnd`, `2ndExonStart_0base`, `2ndExonEnd`, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_1DPA_1am, FL_1DPA_1am))
rMATs_1DPA_5am_MXE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_5am/MXE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_1DPA_5am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_1DPA_5am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, `1stExonStart_0base`, `1stExonEnd`, `2ndExonStart_0base`, `2ndExonEnd`, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_1DPA_5am, FL_1DPA_5am))
## 2DPA MXE
rMATs_2DPA_9am_MXE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_9am/MXE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_2DPA_9am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_2DPA_9am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, `1stExonStart_0base`, `1stExonEnd`, `2ndExonStart_0base`, `2ndExonEnd`, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_2DPA_9am, FL_2DPA_9am))
rMATs_2DPA_1pm_MXE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_1pm/MXE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_2DPA_1pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_2DPA_1pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, `1stExonStart_0base`, `1stExonEnd`, `2ndExonStart_0base`, `2ndExonEnd`, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_2DPA_1pm, FL_2DPA_1pm))
rMATs_2DPA_5pm_MXE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_5pm/MXE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_2DPA_5pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_2DPA_5pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, `1stExonStart_0base`, `1stExonEnd`, `2ndExonStart_0base`, `2ndExonEnd`, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_2DPA_5pm, FL_2DPA_5pm))
rMATs_2DPA_9pm_MXE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_9pm/MXE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_2DPA_9pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_2DPA_9pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, `1stExonStart_0base`, `1stExonEnd`, `2ndExonStart_0base`, `2ndExonEnd`, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_2DPA_9pm, FL_2DPA_9pm))
rMATs_2DPA_1am_MXE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_1am/MXE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_2DPA_1am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_2DPA_1am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, `1stExonStart_0base`, `1stExonEnd`, `2ndExonStart_0base`, `2ndExonEnd`, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_2DPA_1am, FL_2DPA_1am))
rMATs_2DPA_5am_MXE <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_5am/MXE.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_2DPA_5am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_2DPA_5am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, `1stExonStart_0base`, `1stExonEnd`, `2ndExonStart_0base`, `2ndExonEnd`, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_2DPA_5am, FL_2DPA_5am))

## merge MXE
MXE_by_Column <- c("GeneID", "1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE")

rMATs_merge_MXE <- rMATs_0DPA_9am_MXE %>%
    full_join(rMATs_0DPA_1pm_MXE, by = MXE_by_Column) %>%
    full_join(rMATs_0DPA_5pm_MXE, by = MXE_by_Column) %>%
    full_join(rMATs_0DPA_9pm_MXE, by = MXE_by_Column) %>%
    full_join(rMATs_0DPA_1am_MXE, by = MXE_by_Column) %>%
    full_join(rMATs_0DPA_5am_MXE, by = MXE_by_Column) %>%
    full_join(rMATs_1DPA_9am_MXE, by = MXE_by_Column) %>%
    full_join(rMATs_1DPA_1pm_MXE, by = MXE_by_Column) %>%
    full_join(rMATs_1DPA_5pm_MXE, by = MXE_by_Column) %>%
    full_join(rMATs_1DPA_9pm_MXE, by = MXE_by_Column) %>%
    full_join(rMATs_1DPA_1am_MXE, by = MXE_by_Column) %>%
    full_join(rMATs_1DPA_5am_MXE, by = MXE_by_Column) %>%
    full_join(rMATs_2DPA_9am_MXE, by = MXE_by_Column) %>%
    full_join(rMATs_2DPA_1pm_MXE, by = MXE_by_Column) %>%
    full_join(rMATs_2DPA_5pm_MXE, by = MXE_by_Column) %>%
    full_join(rMATs_2DPA_9pm_MXE, by = MXE_by_Column) %>%
    full_join(rMATs_2DPA_1am_MXE, by = MXE_by_Column) %>%
    full_join(rMATs_2DPA_5am_MXE, by = MXE_by_Column)
# write_csv(rMATs_merge_MXE, "mediumDataSave/rMATs_merge_MXE.csv")
rMATs_merge_MXE_WT_exp <- rMATs_merge_MXE %>%
    select(c(GeneID, starts_with("WT_")))
rMATs_merge_MXE_FL_exp <- rMATs_merge_MXE %>%
    select(c(GeneID, starts_with("FL_")))
# write_csv(rMATs_merge_MXE_WT_exp, "mediumDataSave/rMATs_merge_MXE_WT_exp.csv")
# write_csv(rMATs_merge_MXE_FL_exp, "mediumDataSave/rMATs_merge_MXE_FL_exp.csv")

# RI
## 0DPA RI

rMATs_0DPA_9am_RI <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_9am/RI.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_0DPA_9am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_0DPA_9am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, riExonStart_0base, riExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_0DPA_9am, FL_0DPA_9am))
rMATs_0DPA_1pm_RI <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_1pm/RI.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_0DPA_1pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_0DPA_1pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, riExonStart_0base, riExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_0DPA_1pm, FL_0DPA_1pm))
rMATs_0DPA_5pm_RI <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_5pm/RI.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_0DPA_5pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_0DPA_5pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, riExonStart_0base, riExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_0DPA_5pm, FL_0DPA_5pm))
rMATs_0DPA_9pm_RI <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_9pm/RI.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_0DPA_9pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_0DPA_9pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, riExonStart_0base, riExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_0DPA_9pm, FL_0DPA_9pm))
rMATs_0DPA_1am_RI <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_1am/RI.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_0DPA_1am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_0DPA_1am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, riExonStart_0base, riExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_0DPA_1am, FL_0DPA_1am))
rMATs_0DPA_5am_RI <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_0DPA_5am/RI.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_0DPA_5am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_0DPA_5am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, riExonStart_0base, riExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_0DPA_5am, FL_0DPA_5am))
## 1DPA RI
rMATs_1DPA_9am_RI <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_9am/RI.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_1DPA_9am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_1DPA_9am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, riExonStart_0base, riExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_1DPA_9am, FL_1DPA_9am))
rMATs_1DPA_1pm_RI <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_1pm/RI.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_1DPA_1pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_1DPA_1pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, riExonStart_0base, riExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_1DPA_1pm, FL_1DPA_1pm))
rMATs_1DPA_5pm_RI <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_5pm/RI.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_1DPA_5pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_1DPA_5pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, riExonStart_0base, riExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_1DPA_5pm, FL_1DPA_5pm))
rMATs_1DPA_9pm_RI <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_9pm/RI.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_1DPA_9pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_1DPA_9pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, riExonStart_0base, riExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_1DPA_9pm, FL_1DPA_9pm))
rMATs_1DPA_1am_RI <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_1am/RI.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_1DPA_1am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_1DPA_1am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, riExonStart_0base, riExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_1DPA_1am, FL_1DPA_1am))
rMATs_1DPA_5am_RI <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_1DPA_5am/RI.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_1DPA_5am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_1DPA_5am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, riExonStart_0base, riExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_1DPA_5am, FL_1DPA_5am))
## 2DPA RI
rMATs_2DPA_9am_RI <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_9am/RI.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_2DPA_9am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_2DPA_9am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, riExonStart_0base, riExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_2DPA_9am, FL_2DPA_9am))
rMATs_2DPA_1pm_RI <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_1pm/RI.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_2DPA_1pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_2DPA_1pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, riExonStart_0base, riExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_2DPA_1pm, FL_2DPA_1pm))
rMATs_2DPA_5pm_RI <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_5pm/RI.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_2DPA_5pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_2DPA_5pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, riExonStart_0base, riExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_2DPA_5pm, FL_2DPA_5pm))
rMATs_2DPA_9pm_RI <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_9pm/RI.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_2DPA_9pm = map_dbl(IncLevel1, convertToDblAndMean),
        FL_2DPA_9pm = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, riExonStart_0base, riExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_2DPA_9pm, FL_2DPA_9pm))
rMATs_2DPA_1am_RI <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_1am/RI.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_2DPA_1am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_2DPA_1am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, riExonStart_0base, riExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_2DPA_1am, FL_2DPA_1am))
rMATs_2DPA_5am_RI <- data.table::fread("/data1/wangy/circadian/rMATsMetaData/rMATs_output_2DPA_5am/RI.MATS.JC.txt") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(
        WT_2DPA_5am = map_dbl(IncLevel1, convertToDblAndMean),
        FL_2DPA_5am = map_dbl(IncLevel2, convertToDblAndMean)
    ) %>%
    select(c(GeneID, riExonStart_0base, riExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, WT_2DPA_5am, FL_2DPA_5am))
## merge RI
RI_by_Column <- c("GeneID", "riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE")
rMATs_merge_RI <- rMATs_0DPA_9am_RI %>%
    full_join(rMATs_0DPA_1pm_RI, by = RI_by_Column) %>%
    full_join(rMATs_0DPA_5pm_RI, by = RI_by_Column) %>%
    full_join(rMATs_0DPA_9pm_RI, by = RI_by_Column) %>%
    full_join(rMATs_0DPA_1am_RI, by = RI_by_Column) %>%
    full_join(rMATs_0DPA_5am_RI, by = RI_by_Column) %>%
    full_join(rMATs_1DPA_9am_RI, by = RI_by_Column) %>%
    full_join(rMATs_1DPA_1pm_RI, by = RI_by_Column) %>%
    full_join(rMATs_1DPA_5pm_RI, by = RI_by_Column) %>%
    full_join(rMATs_1DPA_9pm_RI, by = RI_by_Column) %>%
    full_join(rMATs_1DPA_1am_RI, by = RI_by_Column) %>%
    full_join(rMATs_1DPA_5am_RI, by = RI_by_Column) %>%
    full_join(rMATs_2DPA_9am_RI, by = RI_by_Column) %>%
    full_join(rMATs_2DPA_1pm_RI, by = RI_by_Column) %>%
    full_join(rMATs_2DPA_5pm_RI, by = RI_by_Column) %>%
    full_join(rMATs_2DPA_9pm_RI, by = RI_by_Column) %>%
    full_join(rMATs_2DPA_1am_RI, by = RI_by_Column) %>%
    full_join(rMATs_2DPA_5am_RI, by = RI_by_Column)
write_csv(rMATs_merge_RI, "mediumDataSave/rMATs_merge_RI.csv")
rMATs_merge_RI_WT_exp <- rMATs_merge_RI %>% 
select(c(GeneID, starts_with("WT_")))
rMATs_merge_RI_FL_exp <- rMATs_merge_RI %>%
select(c(GeneID, starts_with("FL_")))
write_csv(rMATs_merge_RI_WT_exp, "mediumDataSave/rMATs_merge_RI_WT_exp.csv")
write_csv(rMATs_merge_RI_FL_exp, "mediumDataSave/rMATs_merge_RI_FL_exp.csv")

p1 <- salmonPlotRepCirca("Ghir_D13G018480.1", "Proline-rich spliceosome-associated (PSP) family protein / zinc knuckle (CCHC-type) family protein")
ggsave(p1, "figure/transcriptionPlot/Ghir_D13G018480.1_PSP.pdf", width = 10, height = 10)

p2 <- salmonPlotRepCirca("Ghir_A08G008200.1", "spliceosomal protein U1A")
p3 <- salmonPlotRepCirca("Ghir_D13G018480.2", "Proline-rich spliceosome-associated (PSP) family protein / zinc knuckle (CCHC-type) family protein")

p4 <- salmonPlotRepCirca("Ghir_D06G014080.1", "proline-rich spliceosome-associated (PSP) family protein")
p5 <- salmonPlotRepCirca("Ghir_A13G017740.1", "proline-rich spliceosome-associated (PSP) family protein")

p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 2, nrow = 3, guides = "collect")
ggsave("figure/transcriptionPlot/PSP_5.pdf", width = 10, height = 10)

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


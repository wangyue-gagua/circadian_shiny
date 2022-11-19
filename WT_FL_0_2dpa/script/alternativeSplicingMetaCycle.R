library(MetaCycle)
library(tidyverse)
# 在metacycle分析之前，需要将NA值替换为0, GeneID去重以便于下游分析
## A3SS WT
rMATs_merge_A3SS_WT_exp <- read_csv("./mediumDataSave/rMATs_merge_A3SS_WT_exp.csv") %>% 
    mutate(GeneID = str_c(GeneID, "_", seq_along(GeneID)))
rMATs_merge_A3SS_WT_exp[is.na(rMATs_merge_A3SS_WT_exp)] <- 0
write.csv(rMATs_merge_A3SS_WT_exp, "./mediumDataSave/removeNaMeta/rMATs_merge_A3SS_WT_exp.csv", row.names = FALSE)
meta2d(
    infile="mediumDataSave/removeNaMeta/rMATs_merge_A3SS_WT_exp.csv",
    filestyle = "csv",
    timepoints = seq(1, 69, 4),
    outdir = "METACYCLE/rMATs_merge_A3SS_WT_META2d",
    parallelize = TRUE,
    nCores = 10,
    minper = 24,
    maxper = 24,
    ARSdefaultPer = 24,
    cycMethod = c("JTK", "LS"),
)
## A3SS FL
rMATs_merge_A3SS_FL_exp <- read_csv("./mediumDataSave/rMATs_merge_A3SS_FL_exp.csv") %>% 
    mutate(GeneID = str_c(GeneID, "_", seq_along(GeneID)))
rMATs_merge_A3SS_FL_exp[is.na(rMATs_merge_A3SS_FL_exp)] <- 0
write.csv(rMATs_merge_A3SS_FL_exp, "./mediumDataSave/removeNaMeta/rMATs_merge_A3SS_FL_exp.csv", row.names = FALSE)
meta2d(
    infile="mediumDataSave/removeNaMeta/rMATs_merge_A3SS_FL_exp.csv",
    filestyle = "csv",
    timepoints = seq(1, 69, 4),
    outdir = "METACYCLE/rMATs_merge_A3SS_FL_META2d",
    parallelize = TRUE,
    nCores = 10,
    minper = 24,
    maxper = 24,
    ARSdefaultPer = 24,
    cycMethod = c("JTK", "LS"),
)

## A5SS WT
rMATs_merge_A5SS_WT_exp <- read_csv("./mediumDataSave/rMATs_merge_A5SS_WT_exp.csv") %>% 
    mutate(GeneID = str_c(GeneID, "_", seq_along(GeneID)))
rMATs_merge_A5SS_WT_exp[is.na(rMATs_merge_A5SS_WT_exp)] <- 0
write.csv(rMATs_merge_A5SS_WT_exp, "./mediumDataSave/removeNaMeta/rMATs_merge_A5SS_WT_exp.csv", row.names = FALSE)
meta2d(
    infile="mediumDataSave/removeNaMeta/rMATs_merge_A5SS_WT_exp.csv",
    filestyle = "csv",
    timepoints = seq(1, 69, 4),
    outdir = "METACYCLE/rMATs_merge_A5SS_WT_META2d",
    parallelize = TRUE,
    nCores = 10,
    minper = 24,
    maxper = 24,
    ARSdefaultPer = 24,
    cycMethod = c("JTK", "LS"),
)
## A5SS FL
rMATs_merge_A5SS_FL_exp <- read_csv("./mediumDataSave/rMATs_merge_A5SS_FL_exp.csv") %>% 
    mutate(GeneID = str_c(GeneID, "_", seq_along(GeneID)))
rMATs_merge_A5SS_FL_exp[is.na(rMATs_merge_A5SS_FL_exp)] <- 0
write.csv(rMATs_merge_A5SS_FL_exp, "./mediumDataSave/removeNaMeta/rMATs_merge_A5SS_FL_exp.csv", row.names = FALSE)
meta2d(
    infile="mediumDataSave/removeNaMeta/rMATs_merge_A5SS_FL_exp.csv",
    filestyle = "csv",
    timepoints = seq(1, 69, 4),
    outdir = "METACYCLE/rMATs_merge_A5SS_FL_META2d",
    parallelize = TRUE,
    nCores = 10,
    minper = 24,
    maxper = 24,
    ARSdefaultPer = 24,
    cycMethod = c("JTK", "LS"),
)
## MXE WT
rMATs_merge_MXE_WT_exp <- read_csv("./mediumDataSave/rMATs_merge_MXE_WT_exp.csv") %>% 
    mutate(GeneID = str_c(GeneID, "_", seq_along(GeneID)))
rMATs_merge_MXE_WT_exp[is.na(rMATs_merge_MXE_WT_exp)] <- 0
write.csv(rMATs_merge_MXE_WT_exp, "./mediumDataSave/removeNaMeta/rMATs_merge_MXE_WT_exp.csv", row.names = FALSE)
meta2d(
    infile="mediumDataSave/removeNaMeta/rMATs_merge_MXE_WT_exp.csv",
    filestyle = "csv",
    timepoints = seq(1, 69, 4),
    outdir = "METACYCLE/rMATs_merge_MXE_WT_META2d",
    parallelize = TRUE,
    nCores = 10,
    minper = 24,
    maxper = 24,
    ARSdefaultPer = 24,
    cycMethod = c("JTK", "LS"),
)
## MXE FL
rMATs_merge_MXE_FL_exp <- read_csv("./mediumDataSave/rMATs_merge_MXE_FL_exp.csv") %>% 
    mutate(GeneID = str_c(GeneID, "_", seq_along(GeneID)))
rMATs_merge_MXE_FL_exp[is.na(rMATs_merge_MXE_FL_exp)] <- 0
write.csv(rMATs_merge_MXE_FL_exp, "./mediumDataSave/removeNaMeta/rMATs_merge_MXE_FL_exp.csv", row.names = FALSE)
meta2d(
    infile="mediumDataSave/removeNaMeta/rMATs_merge_MXE_FL_exp.csv",
    filestyle = "csv",
    timepoints = seq(1, 69, 4),
    outdir = "METACYCLE/rMATs_merge_MXE_FL_META2d",
    parallelize = TRUE,
    nCores = 10,
    minper = 24,
    maxper = 24,
    ARSdefaultPer = 24,
    cycMethod = c("JTK", "LS"),
)
## RI WT
rMATs_merge_RI_WT_exp <- read_csv("./mediumDataSave/rMATs_merge_RI_WT_exp.csv") %>% 
    mutate(GeneID = str_c(GeneID, "_", seq_along(GeneID)))
rMATs_merge_RI_WT_exp[is.na(rMATs_merge_RI_WT_exp)] <- 0
write.csv(rMATs_merge_RI_WT_exp, "./mediumDataSave/removeNaMeta/rMATs_merge_RI_WT_exp.csv", row.names = FALSE)
meta2d(
    infile="mediumDataSave/removeNaMeta/rMATs_merge_RI_WT_exp.csv",
    filestyle = "csv",
    timepoints = seq(1, 69, 4),
    outdir = "METACYCLE/rMATs_merge_RI_WT_META2d",
    parallelize = TRUE,
    nCores = 10,
    minper = 24,
    maxper = 24,
    ARSdefaultPer = 24,
    cycMethod = c("JTK", "LS"),
)
## RI FL
rMATs_merge_RI_FL_exp <- read_csv("./mediumDataSave/rMATs_merge_RI_FL_exp.csv") %>% 
    mutate(GeneID = str_c(GeneID, "_", seq_along(GeneID)))
rMATs_merge_RI_FL_exp[is.na(rMATs_merge_RI_FL_exp)] <- 0
write.csv(rMATs_merge_RI_FL_exp, "./mediumDataSave/removeNaMeta/rMATs_merge_RI_FL_exp.csv", row.names = FALSE)
meta2d(
    infile="mediumDataSave/removeNaMeta/rMATs_merge_RI_FL_exp.csv",
    filestyle = "csv",
    timepoints = seq(1, 69, 4),
    outdir = "METACYCLE/rMATs_merge_RI_FL_META2d",
    parallelize = TRUE,
    nCores = 10,
    minper = 24,
    maxper = 24,
    ARSdefaultPer = 24,
    cycMethod = c("JTK", "LS"),
)
## SE WT
rMATs_merge_SE_WT_exp <- read_csv("./mediumDataSave/rMATs_merge_SE_WT_exp.csv") %>% 
    mutate(GeneID = str_c(GeneID, "_", seq_along(GeneID)))
rMATs_merge_SE_WT_exp[is.na(rMATs_merge_SE_WT_exp)] <- 0
write.csv(rMATs_merge_SE_WT_exp, "./mediumDataSave/removeNaMeta/rMATs_merge_SE_WT_exp.csv", row.names = FALSE)
meta2d(
    infile="mediumDataSave/removeNaMeta/rMATs_merge_SE_WT_exp.csv",
    filestyle = "csv",
    timepoints = seq(1, 69, 4),
    outdir = "METACYCLE/rMATs_merge_SE_WT_META2d",
    parallelize = TRUE,
    nCores = 10,
    minper = 24,
    maxper = 24,
    ARSdefaultPer = 24,
    cycMethod = c("JTK", "LS"),
)
## SE FL
rMATs_merge_SE_FL_exp <- read_csv("./mediumDataSave/rMATs_merge_SE_FL_exp.csv") %>% 
    mutate(GeneID = str_c(GeneID, "_", seq_along(GeneID)))
rMATs_merge_SE_FL_exp[is.na(rMATs_merge_SE_FL_exp)] <- 0
write.csv(rMATs_merge_SE_FL_exp, "./mediumDataSave/removeNaMeta/rMATs_merge_SE_FL_exp.csv", row.names = FALSE)
meta2d(
    infile="mediumDataSave/removeNaMeta/rMATs_merge_SE_FL_exp.csv",
    filestyle = "csv",
    timepoints = seq(1, 69, 4),
    outdir = "METACYCLE/rMATs_merge_SE_FL_META2d",
    parallelize = TRUE,
    nCores = 10,
    minper = 24,
    maxper = 24,
    ARSdefaultPer = 24,
    cycMethod = c("JTK", "LS"),
)


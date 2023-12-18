library(tidyverse)
library(Seurat)
library(vroom)

args <- commandArgs(trailingOnly = TRUE)

basedir <- args[1]
codedir <- args[2]

# basedir <- "/home/obaranov/projects/scRNA_newPaper/"
# codedir = "/home/obaranov/projects/scRNA_newPaper//code/"

setwd(basedir)

source(paste0(codedir, "/utility/General_Utility.R"))

# source(paste0(basedir, "FunctionArchive/Convert_Rownames.R"))
# source(paste0(basedir, "FunctionArchive/Annotate_Pathways.R"))
# source(paste0(basedir,"/Tabea_Paper/code/utility/General_Utility.R"))

commonlabs <- vroom(paste0(codedir, "/utility/CellTypeMapping.csv"))

# this way of reading it is not working FIX IT

f <- file("RawData/20210220_NasalSwab_NormCounts.txt")
cols <- readLines(f, n = 1)
close(f)
n.cols <- unlist(strsplit(cols, "\t"))
f <- file("RawData/20210220_NasalSwab_RawCounts.txt")
cols <- readLines(f, n = 1)
close(f)
r.cols <- unlist(strsplit(cols, "\t"))


r.counts <- vroom("RawData/20210220_NasalSwab_RawCounts.txt",
    skip = 2,
    col_names = c("GeneSymbol", r.cols)
)

meta <- vroom("RawData/20210701_NasalSwab_MetaData.txt") %>%
    dplyr::slice(-1) %>%
    column_to_rownames("NAME")
umap <- vroom("RawData/20210220_NasalSwab_UMAP.txt") %>% dplyr::slice(-1)

fullziegler <- CreateSeuratObject(r.counts %>% column_to_rownames("GeneSymbol"),
    min.cells = 3, min.genes = 10,
    project = "Ziegler_2021",
    meta.data = meta
)

fullziegler[["percent.mt"]] <- PercentageFeatureSet(fullziegler, pattern = "^MT-|^mt-|^MT.|^mt.")
fullziegler <- subset(fullziegler, subset = nCount_RNA < 20000 & nCount_RNA > 1000 & percent.mt < 20)

# add custom annotation
fullziegler@meta.data <- makeCommonAnno(commonlabs, fullziegler[[]], "Ziegler", "Detailed_Cell_Annotations")

zieglerno19 <- subset(fullziegler, donor_id != "COVID19_Participant19")
saveRDS(zieglerno19, "FilteredData/no19_seurat.rds")

Tcells <- subset(zieglerno19, Coarse_Cell_Annotations == "T Cells")
saveRDS(Tcells, "FilteredData/no19_tcells_seurat.rds")

rm(zieglerno19)
# rm(Tcells)

cases <- meta %>%
    filter(SARSCoV2_PCR_Status == "pos") %>%
    rownames()

r.counts <- vroom("RawData/20210220_NasalSwab_RawCounts.txt",
    skip = 2,
    col_names = c("GeneSymbol", r.cols)
)

r.counts <- r.counts[c("GeneSymbol", cases)]


covidpos <- CreateSeuratObject(r.counts %>% column_to_rownames("GeneSymbol"),
    min.cells = 3, min.genes = 10,
    project = "Ziegler_2021_Covid_infected_cells_only",
    meta.data = meta %>%
        filter(SARSCoV2_PCR_Status == "pos")
)
covidpos[["percent.mt"]] <- PercentageFeatureSet(covidpos, pattern = "^MT-|^mt-|^MT.|^mt.")
covidpos <- subset(covidpos, subset = nCount_RNA < 20000 & nCount_RNA > 500 & percent.mt < 30)

covidposno19 <- subset(covidpos, donor_id != "COVID19_Participant19")
saveRDS(covidposno19, "FilteredData/no19_cases_seurat.rds")

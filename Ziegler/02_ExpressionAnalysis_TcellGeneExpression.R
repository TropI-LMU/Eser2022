library(tidyverse)
library(Seurat)
library(vroom)
library(cluster, quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)

basedir <- args[1]
codedir <- args[2]

# basedir <- "/home/obaranov/projects/"
# setwd(paste0(basedir, "/Ziegler_cell_ImpairedLocalIntrinsicImmunity/"))


setwd(basedir)

source(paste0(codedir, "/utility/General_Utility.R"))

# source(paste0(basedir, "FunctionArchive/Convert_Rownames.R"))
# source(paste0(basedir, "FunctionArchive/Annotate_Pathways.R"))


getDonorN <- function(gene, seuratObj) {
    pos <- GetAssayData(object = seuratObj, slot = "data")[gene, ] > 0
    poscells = which(pos) %>% names
    posdon = seuratObj@meta.data[poscells,'donor_id'] %>% unique()
    return(length(posdon))
}


ziegler = readRDS("FilteredData/no19_cases_seurat.rds")

ziegler <- ziegler %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(pc.genes = ziegler@var.genes, npcs = 20, verbose = FALSE)

ziegler.tcells = subset(ziegler, Coarse_Cell_Annotations == "T Cells")

bt.genes <- vroom(paste0(codedir, "/utility/BTcellgenes.csv"))
impgenes = c('IFNG','TNF','FASLG','CD40LG','IL2','IL10','IL21','PRF1','GZMA','GZMB')

all.genes = c(bt.genes$SYMBOL, impgenes)

counts = GetAssayData(object = ziegler.tcells, slot = "data")

pctexp = all.genes %>% sapply(function(x) {
    sum(counts[x, ] > 0)
}) / dim(counts)[2]

expimp <- all.genes %>% 
            sapply(getDonorN, seuratObj = ziegler.tcells)


names = names(expimp) %>% unique

goi = data.frame(pctexpressed = pctexp[names],
                expressedinN = expimp[names],
                row.names= names )


dir.create("OutputData/markerGenes/")
write_delim(goi %>% rownames_to_column("Gene"), "OutputData/markerGenes/GenesOfInterest_tcells.csv", delim = '\t')

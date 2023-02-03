library(tidyverse)
library(Seurat)
library(vroom)
library(ggpubr)

args <- commandArgs(trailingOnly = TRUE)

basedir_z <- args[1]
basedir_y <- args[2]
codedir <- args[3]

source(paste0(codedir, "/utility/General_Utility.R"))


getDonorN <- function(gene, seuratObj) {
    pos <- GetAssayData(object = seuratObj, slot = "data")[gene, ] > 0
    poscells = which(pos) %>% names
    posdon = seuratObj@meta.data[poscells,'donor_id'] %>% unique()
    return(length(posdon))
}

ziegler = readRDS(paste0(basedir_z,"/FilteredData/no19_cases_seurat.rds"))

ziegler <- ziegler %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(pc.genes = ziegler@var.genes, npcs = 20, verbose = FALSE) %>% 
    RunUMAP(dims = 1:5)

ziegler.tcells = subset(ziegler, Coarse_Cell_Annotations == "T Cells")

genes = c('IFNG','TNF','CD40LG','PRF1','GZMB','FASLG')

pz = FeaturePlot(ziegler.tcells,genes, ncol = 3, order = TRUE)
ggsave(paste0(basedir_z,'/OutputData/Sup6_ziegler.svg'), plot = pz)


airso <- readRDS(paste0(basedir_y,"/FilteredData/airSeurat_adult_cases.rds"))

airso <- airso %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(pc.genes = airso@var.genes, npcs = 20, verbose = FALSE) %>%
    RunUMAP(dims = 1:5)


airso.t = subset(airso, subset = commonLabel == 'T Cells')


genes = c("ENSG00000111537", "ENSG00000232810", "ENSG00000102245", 
        "ENSG00000180644", "ENSG00000100453", "ENSG00000117560")
names(genes) = c('IFNG','TNF','CD40LG','PRF1','GZMB','FASLG')

py = FeaturePlot(airso.t, genes, ncol = 3, order = TRUE)
ggsave(paste0(basedir_y,'/OutputData/Sup6_yoshida.svg'), plot = py)


sup6 = ggarrange(pz,py, ncol = 1)
ggsave(paste0(basedir_z,'/OutputData/Sup6.svg'), 
    plot = sup6, width = 15, height = 20)

library(Seurat)
library(tidyverse)
library(vroom)

args <- commandArgs(trailingOnly = TRUE)

basedir <- args[1]
codedir <- args[2]

source(paste0(codedir, "/utility/General_Utility.R"))


getDonorN <- function(gene, seuratObj) {
    pos <- GetAssayData(object = seuratObj, slot = "data")[gene, ] > 0
    poscells <- which(pos) %>% names()
    posdon <- seuratObj@meta.data[poscells, "donor"] %>% unique()
    return(length(posdon))
}


airso <- readRDS(paste0(basedir, "/FilteredData/airSeurat_adult_cases.rds") )

airso <- airso %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(pc.genes = airso@var.genes, npcs = 20, verbose = FALSE)

airso.t = subset(airso, subset = commonLabel == 'T Cells')


#1. verify the expression of cytotoxic and T cell relevant molecules in case vs control


bt.genes <- vroom(paste0(codedir, "/utility/BTcellgenes.csv"))
impgenes.names <- c("IFNG", "TNF", "FASLG", "CD40LG", "IL2", "IL10", "IL21", "PRF1", "GZMA", "GZMB")
impgenes <- c("ENSG00000111537", "ENSG00000232810", "ENSG00000117560",
             "ENSG00000102245", "ENSG00000109471", "ENSG00000136634",
             "ENSG00000138684", "ENSG00000180644", "ENSG00000145649", 
             "ENSG00000100453")
impdct = impgenes.names
names(impdct) = impgenes

counts <- GetAssayData(object = airso.t, slot = "data")

impgenes.there <- intersect(impgenes, counts %>% rownames())
bt.genes.there <- intersect(bt.genes %>% pull('ENSEMBL'), counts %>% rownames())

pctexp <- c(impgenes.there, bt.genes.there) %>% sapply(function(x) {
    sum(counts[x, ] > 0)
}) / dim(counts)[2]

expimp <-  c(impgenes.there, bt.genes.there)  %>%
    sapply(getDonorN, seuratObj = airso.t)


namedct = get.annotation(names(pctexp), 'ENSEMBL','hsa','SYMBOL', merge = 'first') %>% 
                column_to_rownames('ENSEMBL')

names(expimp) = names(expimp) %>% sapply(function(x){namedct[x,]}) %>% unname
names(pctexp) = names(pctexp) %>% sapply(function(x){namedct[x,]}) %>% unname

names = names(expimp) %>% unique

goi = data.frame(pctexpressed = pctexp[names],
                expressedinN = expimp[names],
                row.names= names )


outdir = paste0(basedir, "/OutputData/markerGenes/")

dir.create(outdir)
write_delim(goi %>% rownames_to_column('Gene'),
            paste0(outdir,"/GenesOfInterest_tcells.csv") , delim = "\t")

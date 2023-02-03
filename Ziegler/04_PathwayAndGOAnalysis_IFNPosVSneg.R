library(tidyverse)
library(Seurat)
library(vroom)
library(cluster, quietly = TRUE)
library(enrichR)
library(KEGGREST)


args <- commandArgs(trailingOnly = TRUE)

basedir <- args[1]
codedir <- args[2]

source(paste0(codedir, "/utility/General_Utility.R"))


setwd(basedir)

# source(paste0(basedir, "FunctionArchive/Convert_Rownames.R"))
# source(paste0(basedir, "FunctionArchive/Annotate_Pathways.R"))

pws <- c()
for (fpath in c("BRITE_Env_SignalingMolecules.txt", "BRITE_Env_SignalTransduction.txt", "BRITE_OrgSys_Immune.txt")) {
    f <- file(paste0(codedir, "/utility/PathwayLists/", fpath))
    pws <- c(pws, readLines(f)[-1])
    close(f)
}

pwlist <- c(
    vroom(paste0(codedir, "/utility/PathwayLists/BRITE_OrgSys_Immune.txt"), delim = "\t") %>% pull(),
    vroom(paste0(codedir, "/utility/PathwayLists/BRITE_Env_SignalTransduction.txt"), delim = "\t") %>% pull(),
    vroom(paste0(codedir, "/utility/PathwayLists/BRITE_Env_SignalingMolecules.txt"), delim = "\t") %>% pull(),
    vroom(paste0(codedir, "/utility/PathwayLists/BRITE_Cell_growth_and_death.txt"), delim = "\t") %>% pull()
) %>% unique()

findInEnrichR <- function(pathway, enrichr.tab) {
    query <- keggGet(pathway)[[1]]$NAME
    name <- strsplit(query, " - Homo sapiens")[[1]][1]
    return(grep(name, enrichr.tab$Term))
}

filterEnrichr <- function(pathway.list, enrichr.tab) {
    found.idx <- sapply(pathway.list, function(x) {
        findInEnrichR(x, enrichr.tab)
    }) %>% unlist()
    enrichr.tab <- enrichr.tab[found.idx, ]
    enrichr.tab["Adjusted.P.value"] <- p.adjust(enrichr.tab["P.value"] %>% pull(), method = "BH")
    enrichr.tab %>% arrange(Adjusted.P.value)
}


markertib = vroom('OutputData/markerGenes/Ziegler_markerGenes_allCellTypes.csv')


setEnrichrSite("Enrichr")
use.db <- c("GO_Biological_Process_2021", "KEGG_2021_Human")


celltypes = markertib %>%
    pull(celltype) %>%
    unique()

ubertib.go = tibble(Term = character(),
                    Overlap = character(),
                    P.value = numeric(),
                    Adjusted.P.value = numeric(),
                    Old.P.value = numeric(),
                    Old.Adjusted.P.value = numeric(),
                    Odds.Ratio = numeric(),
                    Combined.Score = numeric(),
                    Genes = character(),
                    celltype = character())

ubertib.kegg = tibble(Term = character(),
                    Overlap = character(),
                    P.value = numeric(),
                    Adjusted.P.value = numeric(),
                    Old.P.value = numeric(),
                    Old.Adjusted.P.value = numeric(),
                    Odds.Ratio = numeric(),
                    Combined.Score = numeric(),
                    Genes = character(),
                    celltype = character())

for(ct in celltypes){
    print(ct)

    degset <- markertib %>%
        filter(celltype == ct) %>%
        filter(p_val_adj < 0.05) %>%
        filter(avg_log2FC < -0.25 | avg_log2FC > 0.25) %>%
        arrange(avg_log2FC) %>%
        pull(Gene)


    enriched <- enrichr(degset, use.db)
    if (!is.null(enriched)) {
        enriched$KEGG_2021_Human =  filterEnrichr(pwlist, enriched$KEGG_2021_Human)
        enriched$KEGG_2021_Human['celltype'] = ct
        enriched$GO_Biological_Process_2021['celltype'] = ct
    }

    ubertib.kegg = bind_rows(ubertib.kegg, enriched$KEGG_2021_Human)
    ubertib.go = bind_rows(ubertib.go, enriched$GO_Biological_Process_2021)
}



dir.create("OutputData/pathwayEnrichment/")
write.table(ubertib.kegg, "OutputData/pathwayEnrichment/Ziegler_pathwayEnrichment_allCelltypes.csv", sep = "\t", row.names = F)
write.table(ubertib.go, "OutputData/pathwayEnrichment/Ziegler_goEnrichment_allCelltypes.csv", sep = "\t", row.names = F)

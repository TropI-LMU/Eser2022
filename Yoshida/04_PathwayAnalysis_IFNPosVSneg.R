library(tidyverse)
library(Seurat)
library(vroom)
library(cluster, quietly = TRUE)
# library(gage)
# library(clusterProfiler)
library(enrichR)
library(KEGGREST)


args <- commandArgs(trailingOnly = TRUE)

basedir <- args[1]
codedir <- args[2]

source(paste0(codedir, "/utility/General_Utility.R"))


# basedir <- "/home/obaranov/projects/"
# basedir <- paste0(basedir, "/Tabea_Paper/newDataset/Yoshida/OutputData/")


# source(paste0(basedir, "FunctionArchive/Get_Annotation.R"))


markertib <- vroom(paste0(basedir, "/OutputData/markerGenes/Yoshida_markerGenes_allCellTypes.csv"))

print(markertib %>% pull(Gene))

# mapper = get.annotation(markertib %>% pull(Gene), 
#             'ENSEMBL', 'hsa', 'SYMBOL',  merge = 'first') %>%
#             column_to_rownames('ENSEMBL')
# markertib = left_join(markertib,mapper %>% rownames_to_column('Gene'))

setEnrichrSite("Enrichr")
use.db <- c("GO_Biological_Process_2021", "KEGG_2021_Human")


celltypes = markertib %>%
    pull(celltype) %>%
    unique()

ubertib.go <- tibble(
    Term = character(),
    Overlap = character(),
    P.value = numeric(),
    Adjusted.P.value = numeric(),
    Old.P.value = numeric(),
    Old.Adjusted.P.value = numeric(),
    Odds.Ratio = numeric(),
    Combined.Score = numeric(),
    Genes = character(),
    celltype = character()
)

ubertib.kegg <- tibble(
    Term = character(),
    Overlap = character(),
    P.value = numeric(),
    Adjusted.P.value = numeric(),
    Old.P.value = numeric(),
    Old.Adjusted.P.value = numeric(),
    Odds.Ratio = numeric(),
    Combined.Score = numeric(),
    Genes = character(),
    celltype = character()
)

for (ct in celltypes) {
    print(ct)

    degset <- markertib %>%
        filter(celltype == ct) %>%
        filter(p_val_adj < 0.01) %>%
        filter(avg_log2FC < -0.25 | avg_log2FC > 0.25) %>%
        arrange(avg_log2FC) %>%
        pull(Gene)

    print(length(degset))
    if(length(degset) > 0){
        enriched <- enrichr(degset, use.db)
        if (dim(enriched$KEGG_2021_Human)[1] >0 ) {
            enriched$KEGG_2021_Human["celltype"] <- ct
            ubertib.kegg <- bind_rows(ubertib.kegg, enriched$KEGG_2021_Human)
        }
        if (dim(enriched$GO_Biological_Process_2021)[1] >0 ) {
            enriched$GO_Biological_Process_2021["celltype"] <- ct
            ubertib.go <- bind_rows(ubertib.go, enriched$GO_Biological_Process_2021)
        }
    }
}


dir.create(paste0(basedir,"/OutputData/pathwayEnrichment/"))
write.table(ubertib.kegg, paste0(basedir,"/OutputData/pathwayEnrichment/Yoshida_pathwayEnrichment_allCelltypes.csv"), sep = "\t", row.names = F)
write.table(ubertib.go, paste0(basedir, "/OutputData/pathwayEnrichment/Yoshida_goEnrichment_allCelltypes.csv"), sep = "\t", row.names = F)

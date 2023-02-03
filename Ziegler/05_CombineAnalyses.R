library(tidyverse)
library(vroom)


args <- commandArgs(trailingOnly = TRUE)

basedir_z <- args[1]
basedir_y <- args[2]
codedir <- args[3]

source(paste0(codedir, "/utility/General_Utility.R"))

# ziegler tibs
z.tib.deg <- vroom(paste0(basedir_z, "/OutputData/markerGenes/Ziegler_markerGenes_allCellTypes.csv"))
go.z = vroom(paste0(basedir_z,'/OutputData/pathwayEnrichment/Ziegler_goEnrichment_allCelltypes.csv'))

# yosh tibs
y.tib.deg = vroom(paste0(basedir_y, "/OutputData/markerGenes/Yoshida_markerGenes_allCellTypes.csv"))
go.y <- vroom(paste0(basedir_y,"/OutputData/pathwayEnrichment/Yoshida_goEnrichment_allCelltypes.csv"))




go.z.c = go.z %>% 
            filter(celltype == 'Developing Ciliated Cells') %>% 
            filter(Adjusted.P.value < 0.05 )

z.genes <- go.z.c$Genes %>%
    sapply(function(x) {
        str_split(x, ";")
    }) %>%
    unlist() %>%
    unname() %>%
    unique()

z.tib.deg["label"] = ""
z.tib.deg["colour"] <- "black"
impgenes.z = z.tib.deg$Gene %>% 
            sapply(function(x){x %in% z.genes}) %>% 
            unname
z.tib.deg[impgenes.z, "label"] = z.tib.deg[impgenes.z, "Gene"]


go.y.c <- go.y %>%
    filter(celltype == "Ciliated") %>%
    filter(Adjusted.P.value < 0.05)

y.genes = go.y.c$Genes %>% 
    sapply(function(x) {str_split(x, ";")}) %>%
    unlist() %>% 
    unname() %>% 
    unique()

y.tib.deg["label"] <- ""
y.tib.deg["colour"] = "black"
impgenes.y = y.tib.deg$Gene %>% 
            sapply(function(x){x %in% y.genes}) %>% 
            unname
y.tib.deg[impgenes.y, "label"] = y.tib.deg[impgenes.y, "Gene"]

b.genes = intersect(z.genes, y.genes)
impgenes.by = y.tib.deg$Gene %>% 
            sapply(function(x){x %in% b.genes}) %>% 
            unname
y.tib.deg[impgenes.by,"colour"] <- "red"

impgenes.bz = z.tib.deg$Gene %>% 
            sapply(function(x){x %in% b.genes}) %>% 
            unname
z.tib.deg[impgenes.bz, "colour"] <- "red"


# label
# colour
merge(go.z.c, go.y.c, by = "Term") %>% 
        vroom_write(paste0(basedir_z,"/OutputData/Both_goEnrichment_Ciliated.csv"), 
        delim = '\t')

z.tib.deg %>% filter(celltype == 'Developing Ciliated Cells') %>%
    vroom_write(paste0(basedir_z, "/OutputData/markerGenes/Ziegler_markerGenes_ciliated.csv"),
        delim = "\t"
    )

y.tib.deg %>% filter(celltype == 'Ciliated') %>%
    vroom_write(paste0(basedir_y, "/OutputData/markerGenes/Yoshida_markerGenes_ciliated.csv"),
        delim = "\t"
    )

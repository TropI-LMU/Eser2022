library(Seurat)
library(tidyverse)
library(vroom)


args <- commandArgs(trailingOnly = TRUE)

basedir <- args[1]
codedir <- args[2]

source(paste0(codedir, "/utility/General_Utility.R"))


airso <- readRDS(paste0(basedir,"/FilteredData/airSeurat_adult_cases.rds"))

airso <- airso %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(pc.genes = airso@var.genes, npcs = 20, verbose = FALSE)

airso.t <- subset(airso, commonLabel == "T Cells")

#2. take all marker genes detected in Yoshida and run the pathway analysis 
#       but only considering 
#       significant pathways from airso for the respective cell type


ifnpos <- which(airso.t@assays$RNA@counts["ENSG00000111537", ] > 0) %>% names()
ifnpos.donors <- airso@meta.data[ifnpos,] %>% pull(donor) %>% unique()

infpos.donor.cells <- airso[[]] %>%
    filter(donor %in% ifnpos.donors) %>%
    rownames()
airso@meta.data["ifnpos.donor"] <- "neg"
airso@meta.data[infpos.donor.cells, "ifnpos.donor"] <- "pos"


# might make sense to remove thresholds for the final plot
ubertib <- tibble(Gene = character(), p_val = numeric(), avg_log2FC = numeric(), pct.1 = numeric(), pct.2 = numeric(), p_val_adj = numeric(), celltype = character())

for (celltype in airso[[]] %>%
    pull("Cell_type_annotation_level2") %>%
    unique()) {
    fcall <- tryCatch(
        {
            fcall <- FindMarkers(airso %>% subset(subset = Cell_type_annotation_level2 == celltype),
                group.by = "IFNposdonor", ident.1 = TRUE, logfc.threshold = 0
            ) %>%
                rownames_to_column("Gene")
            fcall["celltype"] <- celltype
            fcall
        },
        error = function(cond) {
            print(cond)
            tibble(Gene = character(), p_val = numeric(), avg_log2FC = numeric(), pct.1 = numeric(), pct.2 = numeric(), p_val_adj = numeric(), celltype = character())
        }
    )
    ubertib <- bind_rows(ubertib, fcall)
}

ubertib['GeneEns'] = ubertib['Gene']
tmp = get.annotation(ubertib$Gene, 
                'ENSEMBL', 'hsa', 'SYMBOL',  merge = 'first')
colnames(tmp) = c('GeneEns', 'Gene' )
ubertib = left_join(ubertib %>% dplyr::select(-Gene), tmp, by = "GeneEns")


# volcano plot
ubertib$delabel <- NA
rename <- (ubertib$p_val_adj < 0.01) & (ubertib$avg_log2FC > 0.25 | ubertib$avg_log2FC < -0.25)
ubertib$diffexpressed <- "neg"
ubertib$diffexpressed[rename] <- "pos"
ubertib$delabel[rename] <- ubertib$Gene[rename]

dir.create(paste0(basedir,"/OutputData/markerGenes/"))
write.table(ubertib %>% 
                    dplyr::select(Gene, everything()), 
                            paste0(basedir,"/OutputData/markerGenes/Yoshida_markerGenes_allCellTypes.csv"), 
                            sep = "\t", row.names = F)
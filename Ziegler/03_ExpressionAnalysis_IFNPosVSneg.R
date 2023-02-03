library(Seurat)
library(tidyverse)
library(vroom)


args <- commandArgs(trailingOnly = TRUE)

basedir <- args[1]
codedir <- args[2]

setwd(basedir)
source(paste0(codedir, "/utility/General_Utility.R"))


# basedir <- "/home/obaranov/projects/"

# source(paste0(basedir, "FunctionArchive/Convert_Rownames.R"))
# source(paste0(basedir, "FunctionArchive/Annotate_Pathways.R"))


ziegler.cases = readRDS(paste0(basedir, "/FilteredData/no19_cases_seurat.rds"))
ziegler.tcells = readRDS(paste0(basedir, "/FilteredData/no19_tcells_seurat.rds"))

ziegler.cases <- ziegler.cases %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(pc.genes = ziegler.cases@var.genes, npcs = 20, verbose = FALSE)

# ziegler <- ziegler %>%
#     RunHarmony("donor_id", plot_convergence = TRUE)

ifnpos <- which(ziegler.tcells@assays$RNA@counts["IFNG", ] > 0) %>% names()
ifnpos.donors <- sapply(ifnpos, function(x) {
    paste(str_split(x, "_")[[1]][2:3], collapse = "_")}) %>%
    unname() %>%
    unique()

infpos.donor.cells <- ziegler.cases[[]] %>%
    filter(donor_id %in% ifnpos.donors) %>%
    rownames()
ziegler.cases@meta.data["ifnpos.donor"] <- "neg"
ziegler.cases@meta.data[infpos.donor.cells, "ifnpos.donor"] <- "pos"



# might make sense to remove thresholds for the final plot
ubertib <- tibble(Gene = character(), p_val = numeric(), avg_log2FC = numeric(), pct.1 = numeric(), pct.2 = numeric(), p_val_adj = numeric(), celltype = character())
for (celltype in ziegler.cases[[]] %>%
    pull("Detailed_Cell_Annotations") %>%
    unique()) {
    fcall <- tryCatch(
        {
            fcall <- FindMarkers(ziegler.cases %>% subset(subset = Detailed_Cell_Annotations == celltype),
                group.by = "ifnpos.donor", ident.1 = "pos", logfc.threshold = 0, min.cells.group = 50
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


# volcano plot
ubertib$delabel <- NA
rename <- (ubertib$p_val_adj < 0.01) & (ubertib$avg_log2FC > 0.25 | ubertib$avg_log2FC < -0.25)
ubertib$diffexpressed <- "neg"
ubertib$diffexpressed[rename] <- "pos"
ubertib$delabel[rename] <- ubertib$Gene[rename]


dir.create(paste0(basedir,"/OutputData/markerGenes/"))
write.table(ubertib, 
        paste0(basedir,"/OutputData/markerGenes/Ziegler_markerGenes_allCellTypes.csv"), 
        sep = "\t", row.names = F)

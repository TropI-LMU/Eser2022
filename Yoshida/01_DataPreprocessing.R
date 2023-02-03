library(Seurat)
library(tidyverse)
library(vroom)


args <- commandArgs(trailingOnly = TRUE)

basedir <- args[1]
codedir <- args[2]

source(paste0(codedir, "/utility/General_Utility.R"))


commonlabs = vroom(paste0(codedir, "/utility/CellTypeMapping.csv"))


air <- ReadMtx(paste0(basedir,"/RawData/X.mtx"),
    paste0(basedir,"/RawData/obs.csv"),
    paste0(basedir,"/RawData/var.csv"),
    cell.sep = ",", feature.sep = ",",
    skip.cell = 1, skip.feature = 1, mtx.transpose = TRUE
)
metaa <- vroom(paste0(basedir, "/RawData/obs.csv")) %>% column_to_rownames("index")
airso <- CreateSeuratObject(counts = air, meta.data = metaa)

airso@meta.data <- makeCommonAnno(commonlabs, airso[[]], 
                                "Yoshida", "Cell_type_annotation_level2")

saveRDS(airso, paste0(basedir,"/FilteredData/airSeurat.rds") )
airso_a <- subset(airso, Group == "Adult")
saveRDS(airso_a, paste0(basedir,"/FilteredData/airSeurat_adult.rds") )


airso_a[[]] %>% pull("Cell_type_annotation_level2") %>% table()
tcellids = airso_a[[]] %>%
    filter(grepl("^T CD4|^T CD8*", airso_a$Cell_type_annotation_level2)) %>% rownames()

poscells = GetAssayData(object = airso_a, slot = "count")['ENSG00000111537',tcellids] %>%
            sapply(function(x){x > 0})

airso_a[['IFNposcell']] = FALSE
airso_a@meta.data[which(poscells) %>% names(),'IFNposcell'] = TRUE

posdonor = airso_a@meta.data %>% filter(IFNposcell == TRUE) %>% pull(donor) %>% unique
tmp = sapply(airso_a[['donor']],function(x){x %in% posdonor}) 
airso_a[["IFNposdonor"]] <- FALSE
airso_a@meta.data[tmp,'IFNposdonor'] = TRUE

saveRDS(airso_a,paste0(basedir,"/FilteredData/airSeurat_adult.rds") )

airso_a_nopc <- subset(airso_a, COVID_status != "Post-COVID")
saveRDS(airso_a_nopc, paste0(basedir,"/FilteredData//airSeurat_adult_nopost.rds") )

airso_a_case <- subset(airso_a_nopc, COVID_status == "COVID+")
saveRDS(airso_a_case, paste0(basedir,"/FilteredData/airSeurat_adult_cases.rds") )

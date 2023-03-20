args <- commandArgs(trailingOnly = TRUE)

predir <- ""
basedir <- args[1]
codedir <- args[2]
tissue_type <- args[3]
# predir = '/home/obaranov/projects/'
# basedir = '/Tabea_Paper/debug/'
# codedir = '/Tabea_Paper/code/'
# tissue_type = "CD8"

library(tidyverse)
library(ggplot2)
library(vroom)
library(gage)
library(pheatmap)
library(matrixStats)
library(KEGGREST)


setwd(basedir)
source(paste0(codedir,'/utility/Plot_Utility.R'))
source(paste0(codedir,'/utility/General_Utility.R'))

# import basic tables
cellfreq <- read.csv('data/bulkRNA/rawData/CellFractions.csv', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
cellfreq['Sample_ID'] = cellfreq %>% pull('Sample_ID') %>% sapply(function(x){as.character(x)})
metadata <- read.csv("data/bulkRNA/rawData/Metadata.csv", header = TRUE, stringsAsFactors = FALSE, sep = '\t')
metadata['Sample_ID'] = metadata %>% pull('Sample_ID') %>% sapply(function(x){as.character(x)})

metadata_ct <- metadata[metadata$Tissue  == tissue_type,]
counts <- read.csv(paste0("data/bulkRNA/rawData/",tissue_type,"_counts.csv"), sep = '\t', check.names=FALSE) %>%
                column_to_rownames('Gene')

# select the sample & control

case = metadata_ct %>% filter(TAG == 'Tag4')
week1 = cellfreq %>% right_join(case, by = 'Sample_ID') %>%
                filter(weeks_after_onset_of_symptoms == 1) %>%
                pull('Sample_ID') %>% unique

ctl = metadata_ct %>% filter(TAG == 'Ctrl')
d4idx = 1:length(week1)
ctlidx = (length(week1) + 1):(length(week1) + dim(ctl)[1])
counts = as.matrix(counts[,c(week1,ctl %>% pull('Sample_ID'))])


counts = EnsToEnt(counts)

# case = metadata_ct %>% filter(TAG == 'Tag4')


# ctl = metadata_ct %>% filter(TAG == 'Ctrl')



# caseidx = 1:length(case)
# ctlidx = (length(case) + 1):(length(case) + dim(ctl)[1])
# counts = as.matrix(counts[,c(case %>% pull('Sample_ID'), ctl %>% pull('Sample_ID'))])

# id conversion to entrez id
# counts = EnsToEnt(counts)



list.immu = readList(paste0(codedir,'/utility/PathwayLists/GeneSet_Immune.gmt'))
list.sigtrans = readList(paste0(codedir,'/utility/PathwayLists/GeneSet_SignalingMolecules.gmt'))
list.sigmol = readList(paste0(codedir,'/utility/PathwayLists/GeneSet_SignalTransduction.gmt'))
list.death = readList(paste0(codedir,'/utility/PathwayLists/GeneSet_Cell_growth_and_death.gmt'))

kegg.set = list.immu %>% append(list.sigtrans) %>% append(list.sigmol) %>% append(list.death)


gsea.kegg <- gage(counts, gsets = kegg.set, ref = ctlidx,
                    samp = caseidx, compare='unpaired', set.size = c(2,500))

#export
rnames.up = rownames(gsea.kegg$greater)
rnames.down = rownames(gsea.kegg$less)
gsea.up = gsea.kegg$greater %>% as_tibble
gsea.down = gsea.kegg$less %>% as_tibble
gsea.up['PathwayID'] = rnames.up
gsea.down['PathwayID'] = rnames.down
gsea.up['PathwayName'] = ''
gsea.down['PathwayName'] = ''

for(idx in rownames(gsea.up)){
    gsea.up[idx,'PathwayName'] = keggGet(gsea.up[idx,'PathwayID'])[[1]]$NAME
    gsea.down[idx,'PathwayName'] = keggGet(gsea.down[idx,'PathwayID'])[[1]]$NAME
}

gsea.up.df = gsea.up %>% dplyr::select('PathwayID','PathwayName',everything()) %>%
    as.data.frame
gsea.down.df = gsea.down %>% dplyr::select('PathwayID','PathwayName',everything()) %>%
    as.data.frame

dir.create(paste0('result_tables/',tissue_type), showWarnings = TRUE, recursive = TRUE)
gsea.up.df  %>% write_delim(paste0('result_tables/',tissue_type,'/PathwayEnrichment_Up.csv'), delim = '\t')
gsea.down.df  %>% write_delim(paste0('result_tables/',tissue_type,'/PathwayEnrichment_Down.csv'), delim = '\t')

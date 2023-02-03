args <- commandArgs(trailingOnly = TRUE)

predir = ''
basedir = args[1]
codedir = args[2]
tissue_type = args[3]
# predir = '/run/user/1000/gvfs/smb-share:server=intmedicine2.local,share=home/Drive/'
# basedir = '/Tabea_Paper/debug/'
# codedir = '/Tabea_Paper/code/'
# tissue_type = "NK"

library(tidyverse)
library(KEGGREST)
library(gage)


setwd(paste0(predir, basedir))
source(paste0(predir, codedir,'/utility/General_Utility.R'))


# read precalculated data #
###########################
virload = read.csv("data/bulkRNA/rawData/Virusload.csv",
            header = TRUE, stringsAsFactors = FALSE, sep = '\t') %>% drop_na
cellfreq <- read.csv("data/bulkRNA/rawData/CellFractions.csv",
                            header = TRUE, stringsAsFactors = FALSE,  sep = '\t', check.names = FALSE)
c_counts <- read.csv(paste0("data/bulkRNA/rawData/",tissue_type,"_counts.csv"), sep = '\t', check.names=FALSE) %>%
                column_to_rownames('Gene')
btSet = read.csv("data/additionalData/BTcellgenes.csv", sep = '\t')
# filter metadata
metadata <- read.csv("data/bulkRNA/rawData//Metadata.csv", header = TRUE, stringsAsFactors = FALSE, sep = '\t')
metadata <- metadata[metadata$Tissue  == tissue_type,]
metadata['Sample_ID'] = metadata %>% pull('Sample_ID') %>% sapply(function(x){as.character(x)})
# get gene sets
list.immu = readList(paste0(predir,codedir,'/utility/PathwayListst/GeneSet_Immune.gmt'))
list.sigtrans = readList(paste0(predir,codedir,'/utility/PathwayListst/GeneSet_SignalingMolecules.gmt'))
list.sigmol = readList(paste0(predir,codedir,'/utility/PathwayListst/GeneSet_SignalTransduction.gmt'))
list.death = readList(paste0(predir,codedir,'/utility/PathwayListst/GeneSet_Cell_growth_and_death.gmt'))
kegg.set = list.immu %>% append(list.sigtrans) %>% append(list.sigmol) %>% append(list.death)


#
mieps = rownames(c_counts) %>% sapply(function(x){strsplit(x,'\\.')[[1]][1]}) %>% unname
c_counts = c_counts[!duplicated(mieps),]
rownames(c_counts) = mieps[!duplicated(mieps)]

cd4.ifn = cellfreq %>% filter(stimulant == "nucleocapsid") %>%
    filter(Sample_ID %in% colnames(c_counts)) %>%
    `[`(c('Sample_ID','CD4IFNg | Freq. of Parent')) %>%
    column_to_rownames('Sample_ID')

cd4.any = cellfreq %>% filter(stimulant == "nucleocapsid") %>%
    filter(Sample_ID %in% colnames(c_counts)) %>%
    `[`(c('Sample_ID','CD4any | Freq. of Parent')) %>%
    column_to_rownames('Sample_ID')

virload = cd4.ifn %>% rownames_to_column('Sample_ID') %>% inner_join(virload, by = 'Sample_ID')
virload['log.VirusLoad'] = log(virload %>% pull('VirusLoad'))
vir.vector = virload[c('log.VirusLoad','Sample_ID')] %>% column_to_rownames('Sample_ID')

genecor = data.frame(matrix(ncol = 7, nrow = dim(btSet)[1]))
colnames(genecor) = c('symbol','cd4ifng.spearman_rho','cd4ifng.p_value','cd4any.spearman_rho','cd4any.p_value',
                            'vir.spearman_rho','vir.p_value')
rownames(genecor) = btSet$ENSEMBL
genecor$symbol = btSet$SYMBOL
for(gene in btSet$ENSEMBL){
    if(gene %in% rownames(c_counts)){
        correlation = cor.test(as.data.frame(t(c_counts[gene,])) %>% pull, cd4.ifn %>% pull, method = 'spearman')
        genecor[gene,'cd4ifng.spearman_rho']  = correlation$estimate[['rho']]
        genecor[gene,'cd4ifng.p_value'] = correlation$p.value

        correlation = cor.test(as.data.frame(t(c_counts[gene,])) %>% pull, cd4.any %>% pull, method = 'spearman')
        genecor[gene,'cd4any.spearman_rho']  = correlation$estimate[['rho']]
        genecor[gene,'cd4any.p_value'] = correlation$p.value
        # virus load is missing for many patients --> filter
        correlation = cor.test(as.data.frame(t(c_counts[gene,rownames(vir.vector)])) %>% pull, vir.vector %>% pull, method = 'spearman')
        genecor[gene,'vir.spearman_rho']  = correlation$estimate[['rho']]
        genecor[gene,'vir.p_value'] = correlation$p.value
    }
}

genecor = genecor %>% arrange(cd4ifng.p_value,cd4any.p_value,vir.p_value) %>% drop_na()

dir.create(paste0("result_tables/",tissue_type), showWarnings = FALSE)
genecor  %>% write.table(paste0("result_tables/",tissue_type,"/GeneCorrelation.csv"),
                row.names = TRUE, sep = '\t', quote = FALSE, col.names=NA)

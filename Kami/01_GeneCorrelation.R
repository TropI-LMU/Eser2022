args <- commandArgs(trailingOnly = TRUE)

predir = ''
basedir = args[1]
codedir = args[2]
tissue_type = args[3]
# predir = "/home/obaranov/projects/"
# basedir = '/Tabea_Paper/debug/'
# codedir = '/Tabea_Paper/code/'
# tissue_type = "CD8"

library(tidyverse)
library(KEGGREST)
library(gage)


setwd(paste0(predir, basedir))
source(paste0(predir, codedir,'/utility/General_Utility.R'))


# read precalculated data #
###########################
virload = read.csv("data/bulkRNA/rawData/Figure4_viruslast.csv",
            header = TRUE, stringsAsFactors = FALSE, sep = '\t') %>% drop_na
cellfreq <- read.csv("data/bulkRNA/rawData/CellFractions.csv",
                            header = TRUE, stringsAsFactors = FALSE,  sep = '\t', check.names = FALSE) %>% 
                            rename(subject_ID = 'Patient_ID') 
c_counts <- read.csv(paste0("data/bulkRNA/rawData/",tissue_type,"_counts.csv"), sep = '\t', check.names=FALSE) %>%
                column_to_rownames('Gene')
btSet = read.csv("data/additionalData/BTcellgenes.csv", sep = '\t')
# filter metadata
metadata <- read.csv("data/bulkRNA/rawData//Metadata.csv", header = TRUE, stringsAsFactors = FALSE, sep = '\t')
metadata <- metadata[metadata$Tissue  == tissue_type,]
metadata['Sample_ID'] = metadata %>% pull('Sample_ID') %>% sapply(function(x){as.character(x)})

# filter the count data to contain only the relevant genes 

c_counts["shortENSEMBL"] <- rownames(c_counts) %>%
    sapply(function(x) {
        strsplit(x, "\\.")[[1]][1]
    })


mieps = rownames(c_counts) %>% sapply(function(x){strsplit(x,'\\.')[[1]][1]}) %>% unname
c_counts = c_counts[!duplicated(mieps),]
rownames(c_counts) = mieps[!duplicated(mieps)]

metadata_d4 <- metadata %>% filter(TAG == "Tag4")

keep_tabi <- cellfreq %>%
    filter(stimulant == "nucleocapsid") %>%
    filter(weeks_after_onset_of_symptoms == 1) %>%
    pull("Sample_ID")

keepsamp <- intersect(keep_tabi, metadata_d4 %>% pull("Sample_ID"))
c_counts <- as.matrix(c_counts[, keepsamp])

cd4.ifn = cellfreq %>% filter(stimulant == "nucleocapsid") %>%
    filter(Sample_ID %in% colnames(c_counts)) %>%
    `[`(c('Sample_ID','Patient_ID','CD4IFNg | Freq. of Parent')) %>%
    column_to_rownames('Sample_ID')

cd8.ifn= cellfreq %>% filter(stimulant == "nucleocapsid") %>%
    filter(Sample_ID %in% colnames(c_counts)) %>%
    `[`(c('Sample_ID','CD8IFNg | Freq. of Parent')) %>%
    column_to_rownames('Sample_ID')

virload = cd4.ifn %>% rownames_to_column('Sample_ID') %>% inner_join(virload, by = 'Patient_ID')
vir.vector = virload[c('log.VirusLoad','Sample_ID')] %>% column_to_rownames('Sample_ID')

genecor = data.frame(matrix(ncol = 7, nrow = dim(btSet)[1]))
colnames(genecor) = c('symbol','cd4ifng.spearman_rho','cd4ifng.p_value','cd8ifng.spearman_rho','cd8ifng.p_value',
                            'vir.spearman_rho','vir.p_value')
rownames(genecor) = btSet$ENSEMBL
genecor$symbol = btSet$SYMBOL

c_counts <- as.data.frame(c_counts)[rownames(cd4.ifn)]

for(gene in btSet$ENSEMBL){
    if(gene %in% rownames(c_counts)){
        correlation = cor.test(as.data.frame(t(c_counts[gene,])) %>% pull, cd4.ifn %>% pull, method = 'spearman')
        genecor[gene,'cd4ifng.spearman_rho']  = correlation$estimate[['rho']]
        genecor[gene,'cd4ifng.p_value'] = correlation$p.value

        correlation = cor.test(as.data.frame(t(c_counts[gene,])) %>% pull, cd8.ifn%>% pull, method = 'spearman')
        genecor[gene,'cd8ifng.spearman_rho']  = correlation$estimate[['rho']]
        genecor[gene,'cd8ifng.p_value'] = correlation$p.value
        # virus load is missing for many patients --> filter
        correlation = cor.test(as.data.frame(t(c_counts[gene,rownames(vir.vector)])) %>% pull, vir.vector %>% pull, method = 'spearman')
        genecor[gene,'vir.spearman_rho']  = correlation$estimate[['rho']]
        genecor[gene,'vir.p_value'] = correlation$p.value
    }
}

genecor = genecor %>% arrange(cd4ifng.p_value,cd8ifng.p_value,vir.p_value) %>% drop_na()

dir.create(paste0("result_tables/",tissue_type), showWarnings = FALSE)
genecor  %>% write.table(paste0("result_tables/",tissue_type,"/GeneCorrelation.csv"),
                row.names = TRUE, sep = '\t', quote = FALSE, col.names=NA)

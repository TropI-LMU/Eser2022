args <- commandArgs(trailingOnly = TRUE)

predir = ''
basedir = args[1]
codedir = args[2]
tissue_type = args[3]
# predir = '/home/obaranov/projects/'
# basedir = '/Tabea_Paper/debug/'
# codedir = '/Tabea_Paper/code/'
# tissue_type = "CD4"

library(tidyverse)
library(magrittr)
library(gage)
library(vroom)

# this script needs to run from console, as it needs a lot of time and resources


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

# data preprocessing
mieps <- rownames(c_counts) %>%
    sapply(function(x) {
        strsplit(x, "\\.")[[1]][1]
    }) %>%
    unname()
c_counts <- c_counts[!duplicated(mieps), ]
rownames(c_counts) <- mieps[!duplicated(mieps)]

metadata_d4 <- metadata %>% filter(TAG == "Tag4")

keep_tabi <- cellfreq %>%
    filter(stimulant == "nucleocapsid") %>%
    filter(weeks_after_onset_of_symptoms == 1) %>%
    pull("Sample_ID")

keepsamp <- intersect(keep_tabi, metadata_d4 %>% pull("Sample_ID"))
c_counts <- as.matrix(c_counts[, keepsamp])

cd4.ifn <- cellfreq %>%
    filter(stimulant == "nucleocapsid") %>%
    filter(Sample_ID %in% colnames(c_counts)) %>%
    `[`(c('Sample_ID','Patient_ID','CD4IFNg | Freq. of Parent')) %>%
    column_to_rownames("Sample_ID")

cd8.ifn <- cellfreq %>%
    filter(stimulant == "nucleocapsid") %>%
    filter(Sample_ID %in% colnames(c_counts)) %>%
    `[`(c("Sample_ID", "CD8any | Freq. of Parent")) %>%
    column_to_rownames("Sample_ID")

virload <- cd4.ifn %>%
    rownames_to_column("Sample_ID") %>%
    inner_join(virload, by = "Patient_ID")
vir.vector <- virload[c("log.VirusLoad", "Sample_ID")] %>% column_to_rownames("Sample_ID")

genecor <- data.frame(matrix(ncol = 7, nrow = dim(btSet)[1]))
colnames(genecor) <- c(
    "symbol", "cd4ifng.spearman_rho", "cd4ifng.p_value", "cd8ifng.spearman_rho", "cd8ifng.p_value",
    "vir.spearman_rho", "vir.p_value"
)
rownames(genecor) <- btSet$ENSEMBL
genecor$symbol <- btSet$SYMBOL

c_counts <- as.data.frame(c_counts)[rownames(cd4.ifn)]


virload = cd4.ifn %>% rownames_to_column('Sample_ID') %>% inner_join(virload, by = 'Sample_ID')
vir.vector = virload[c('log.VirusLoad','Sample_ID')] %>% column_to_rownames('Sample_ID')
cd4.ifn = cd4.ifn %>% dplyr::select(-Patient_ID)

# to get the confidence intervals, we need to bootstrap the data
#    process:   1. get a bootstrapped set of sample colnames
#               2. filter / reorder the tables accordingly
#               3. calculate the scores, save them
#               4. get the 95% range

cd8i.tib = tibble(gene = btSet$ENSEMBL, SYMBOL = btSet$SYMBOL)
cd4i.tib = tibble(gene = btSet$ENSEMBL, SYMBOL = btSet$SYMBOL)
vir.tib = tibble(gene = btSet$ENSEMBL, SYMBOL = btSet$SYMBOL)
iter =1
gene = btSet$ENSEMBL[1]
minifun = function(gene, newcols, valB){
    if(gene %in% rownames(c_counts)){
        tmp = cor.test(as.data.frame(t(c_counts[gene,newcols])) %>% pull, valB[newcols,], method = 'spearman')
        return(tmp$estimate[['rho']])
    } else {return(NaN) }
}


for(iter in 1:1000){
    newcols = sample(colnames(c_counts), size = dim(c_counts)[2], replace = TRUE )
    correlation = sapply(btSet$ENSEMBL, function(x){minifun(x, newcols,cd8.ifn)})
    cd8i.tib[paste0('boot',iter)]  = correlation

    correlation = sapply(btSet$ENSEMBL, function(x){minifun(x, newcols,cd4.ifn)})
    cd4i.tib[paste0('boot',iter)]  = correlation

    newcols = sample(rownames(vir.vector), size = dim(vir.vector)[1], replace = TRUE )

    correlation = sapply(btSet$ENSEMBL, function(x){minifun(x, newcols,vir.vector)})
    vir.tib[paste0('boot',iter)]  = correlation

}

dir.create(paste0("result_tables/",tissue_type,"/bootstrap_pw/"), showWarnings = FALSE)
write.table(cd8i.tib, paste0("result_tables/",tissue_type,"/bootstrap_pw/bootstrapped_cd8ifng.gene.csv"),
            row.names = FALSE, sep = '\t')
write.table(cd4i.tib, paste0("result_tables/",tissue_type,"/bootstrap_pw/bootstrapped_cd4ifng.gene.csv"),
            row.names = FALSE, sep = '\t')
write.table(vir.tib, paste0("result_tables/",tissue_type,"/bootstrap_pw/bootstrapped_vir.gene.csv"),
            row.names = FALSE, sep = '\t')


# get confidence interval

cd8i.tib = cd8i.tib %>% column_to_rownames('gene')
cd4i.tib = cd4i.tib %>% column_to_rownames('gene')
vir.tib = vir.tib %>% column_to_rownames('gene')

minifunint = function(gene, tib){
    nums <- unlist(lapply(tib, is.numeric), use.names = FALSE)
    tibsub = tib[gene,nums]
    tibsub = tibsub[! is.na(tibsub)]
    ord = order(tibsub)
    q025 = (length(ord) * 0.025) %>% round
    q975 = (length(ord) * 0.975) %>% round
    if(q975 == length(ord)){
        upper = NaN
    } else {
        upper = tibsub[ord[q975]]
    }
    if(q025 == 0){
        lower = NaN
    } else {
        lower = tibsub[ord[q025]]
    }
    return(c(lower, upper) )
}

cd8i.int= rownames(cd8i.tib) %>% sapply(function(gene){minifunint(gene, cd8i.tib)}) %>%
                t() %>% as.data.frame %>% set_colnames(c("cd8ifng.lower", "cd8ifng.upper"))

cd4i.int = rownames(cd4i.tib) %>% sapply(function(gene){minifunint(gene, cd4i.tib)}) %>%
                t() %>% as.data.frame %>% set_colnames(c("cd4ifng.lower", "cd4ifng.upper"))

vir.int = rownames(vir.tib) %>% sapply(function(gene){minifunint(gene, vir.tib)}) %>%
                t() %>% as.data.frame %>% set_colnames(c("vir.lower", "vir.upper"))


gene.tib = vroom( paste0("result_tables/",tissue_type,'/GeneCorrelation.csv'), show_col_types = FALSE)
tmp = colnames(gene.tib)
colnames(gene.tib) = c('gene',tmp[2:length(tmp)])
gene.tib = gene.tib %>%
    right_join(cd4i.int %>% rownames_to_column('gene'), by = 'gene') %>%
    right_join(cd8i.int%>% rownames_to_column('gene'), by = 'gene') %>%
    right_join(vir.int %>% rownames_to_column('gene'), by = 'gene') %>%
    dplyr::select('gene','symbol',starts_with('cd4ifng'), starts_with('cd8ifng'),starts_with('vir'))

gene.tib %>% mutate(across(where(is.numeric), round, 4)) %>%
          write.table( paste0("result_tables/",tissue_type,"/GeneCorrelation_ci.csv"),
                        row.names = FALSE, sep = '\t')

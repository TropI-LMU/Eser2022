#!/usr/bin/R

args <- commandArgs(trailingOnly = TRUE)

predir = ''
basedir = args[1]
codedir = args[2]
tissue_type = args[3]
# predir = "/home/obaranov/projects/"
# basedir = '/Tabea_Paper/debug/'
# codedir = '/Tabea_Paper/code/'
# tissue_type = "Monocytes"

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
# filter metadata
metadata <- read.csv("data/bulkRNA/rawData//Metadata.csv", header = TRUE, stringsAsFactors = FALSE, sep = '\t')
metadata <- metadata[metadata$Tissue  == tissue_type,]
metadata['Sample_ID'] = metadata %>% pull('Sample_ID') %>% sapply(function(x){as.character(x)})
# get gene sets
list.immu = readList(paste0(predir,codedir,'/utility/PathwayLists/GeneSet_Immune.gmt'))
list.sigtrans = readList(paste0(predir,codedir,'/utility/PathwayLists/GeneSet_SignalingMolecules.gmt'))
list.sigmol = readList(paste0(predir,codedir,'/utility/PathwayLists/GeneSet_SignalTransduction.gmt'))
list.death = readList(paste0(predir,codedir,'/utility/PathwayLists/GeneSet_Cell_growth_and_death.gmt'))
kegg.set = list.immu %>% append(list.sigtrans) %>% append(list.sigmol) %>% append(list.death)


# data processing
counts_sub = EnsToEnt(c_counts)


cd4.ifn <- cellfreq %>%
    filter(stimulant == "nucleocapsid") %>%
    filter(weeks_after_onset_of_symptoms == 1) %>%
    filter(Sample_ID %in% colnames(counts_sub)) %>%
    `[`(c('Sample_ID','Patient_ID','CD4IFNg | Freq. of Parent')) %>%
    column_to_rownames("Sample_ID")

cd4.any <- cellfreq %>%
    filter(stimulant == "nucleocapsid") %>%
    filter(weeks_after_onset_of_symptoms == 1) %>%
    filter(Sample_ID %in% colnames(counts_sub)) %>%
    `[`(c("Sample_ID", "CD4any | Freq. of Parent")) %>%
    column_to_rownames("Sample_ID")


mean.df.full = calc.pwmeans(counts_sub[rownames(cd4.ifn)], kegg.set)
mean.df.full %>% write.table(paste0("result_tables/Pathways/",tissue_type,"_PathwayMeanExp.csv"),
                row.names = FALSE, sep = '\t')

correlation = calc.spearman.pw(mean.df.full, cd4.ifn)
mean.df.full['cd4ifng.spearman.rho']  = correlation$rhos
mean.df.full['cd4ifng.p.value'] = correlation$p.values

correlation = calc.spearman.pw(mean.df.full, cd4.any)
mean.df.full['cd4any.spearman.rho']  = correlation$rhos
mean.df.full['cd4any.p.value'] = correlation$p.values

cd8.ifn = cellfreq %>% filter(stimulant == "nucleocapsid") %>%
    filter(Sample_ID %in% colnames(mean.df.full)) %>%
    `[`(c('Sample_ID','CD8IFNg | Freq. of Parent')) %>%
    column_to_rownames('Sample_ID')

cd8.any = cellfreq %>% filter(stimulant == "nucleocapsid") %>%
    filter(Sample_ID %in% colnames(mean.df.full)) %>%
    `[`(c('Sample_ID','CD8any | Freq. of Parent')) %>%
    column_to_rownames('Sample_ID')

correlation = calc.spearman.pw(mean.df.full, cd8.ifn)
mean.df.full['cd8ifng.spearman.rho']  = correlation$rhos
mean.df.full['cd8ifng.p.value'] = correlation$p.values

correlation = calc.spearman.pw(mean.df.full, cd8.any)
mean.df.full['cd8any.spearman.rho']  = correlation$rhos
mean.df.full['cd8any.p.value'] = correlation$p.values

virload = cd4.ifn %>% rownames_to_column('Sample_ID') %>% inner_join(virload, by = 'Patient_ID')
vir.vector = virload[c('log.VirusLoad','Sample_ID')] %>% column_to_rownames('Sample_ID')
cd4.ifn <- cd4.ifn %>% dplyr::select(-Patient_ID)
# virus load is missing for many patients --> filter
correlation = calc.spearman.pw(mean.df.full, vir.vector)
mean.df.full['vir.spearman.rho']  = correlation$rhos
mean.df.full['vir.p.value'] = correlation$p.values

# export
#annotate pathways
pw.tib = annotatePW(c(paste0(predir,basedir,'data/additionalData/geneSets/GeneSet_SignalingMolecules.gmt'),
                        paste0(predir,basedir,'data/additionalData/geneSets/GeneSet_SignalTransduction.gmt'),
                        paste0(predir,basedir,'data/additionalData/geneSets/GeneSet_Immune.gmt'),
                        paste0(predir,basedir,'data/additionalData/geneSets/GeneSet_Cell_growth_and_death.gmt') ))

mean.df.full = mean.df.full %>% right_join(pw.tib, by = 'PathwayID') %>%
                dplyr::select('PathwayID', 'Description',
                ends_with('spearman.rho'), ends_with('p.value'), everything()) %>%
                mutate(across(where(is.numeric), round, 4))


dir.create(paste0("result_tables/",tissue_type), showWarnings = FALSE)
mean.df.full  %>% arrange(cd4ifng.spearman.rho) %>%
                write.table(paste0("result_tables/",tissue_type,"/PathwayCorrelation.csv"),
                row.names = FALSE, sep = '\t')

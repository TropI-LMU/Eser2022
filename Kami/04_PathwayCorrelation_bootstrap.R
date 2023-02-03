args <- commandArgs(trailingOnly = TRUE)

predir = ''
basedir = args[1]
codedir = args[2]
tissue_type = args[3]
# predir = '/media/olga/40A732655125C6DB/Arbeit/'
# basedir = '/Tabea_Paper/debug/'
# codedir = '/Tabea_Paper/code/'
# tissue_type = "NK"

library(tidyverse)
library(gage)
library(vroom)

# this script needs to run from console, as it kills atom


# source(paste0(basedir,'FunctionArchive/Convert_Rownames.R'))
# source(paste0(basedir,'FunctionArchive/Calculate_correlations.R'))
# source(paste0(basedir,'FunctionArchive/Calculate_set_mean.R'))
# source(paste0(basedir,'FunctionArchive/Annotate_Pathways.R'))

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

# data preprocessing
counts_sub = EnsToEnt(c_counts)


cd4.ifn = cellfreq %>% filter(stimulant == "nucleocapsid") %>%
    filter(Sample_ID %in% colnames(counts_sub)) %>%
    `[`(c('Sample_ID','CD4IFNg | Freq. of Parent')) %>%
    column_to_rownames('Sample_ID')

cd4.any = cellfreq %>% filter(stimulant == "nucleocapsid") %>%
    filter(Sample_ID %in% colnames(counts_sub)) %>%
    `[`(c('Sample_ID','CD4any | Freq. of Parent')) %>%
    column_to_rownames('Sample_ID')



virload = cd4.ifn %>% rownames_to_column('Sample_ID') %>% inner_join(virload, by = 'Sample_ID')
virload['log.VirusLoad'] = log(virload %>% pull('VirusLoad'))
vir.vector = virload[c('log.VirusLoad','Sample_ID')] %>% column_to_rownames('Sample_ID')

n.biggest = kegg.set %>% sapply(length) %>% max
n.smallest = kegg.set %>% sapply(length) %>% min

cor.tib = tibble(a = c('smallest.cd4ifng','smallest.cd4any', 'smallest.vir',
                        'biggest.cd4ifng','biggest.cd4any', 'biggest.vir')) %>% column_to_rownames('a')

# to get the confidence intervals, we need to bootstrap the data
#    process:   1. get a bootstrapped set of sample colnames
#               2. filter / reorder the tables accordingly
#               3. calculate the scores, save them
#               4. get the 95% range

cd4.tib = tibble(PathwayID = names(kegg.set))
cd4i.tib = tibble(PathwayID = names(kegg.set))
vir.tib = tibble(PathwayID = names(kegg.set))
iter =1
for(iter in 1:1000){
    newcols = sample(colnames(counts_sub), size = dim(counts_sub)[2], replace = TRUE )
    mean.df = calc.pwmeans(counts_sub[newcols], kegg.set)

    cd4.ifn.sub = cd4.ifn[newcols,,drop = FALSE]
    correlation = calc.spearman.pw(mean.df, cd4.ifn.sub)
    cd4i.tib[paste0('boot',iter)]  = correlation$rhos

    cd4.any.sub = cd4.any[newcols,,drop = FALSE]
    correlation = calc.spearman.pw(mean.df, cd4.any.sub)
    cd4.tib[paste0('boot',iter)]  = correlation$rhos

    newcols = sample(rownames(vir.vector), size = dim(vir.vector)[1], replace = TRUE )
    mean.df = calc.pwmeans(counts_sub[newcols], kegg.set)
    vir.vector.sub = vir.vector[newcols,,drop = FALSE]
    correlation = calc.spearman.pw(mean.df, vir.vector.sub)
    vir.tib[paste0('boot',iter)]  = correlation$rhos

    rnd.genes = list(biggest = sample(rownames(counts_sub), n.biggest, replace = FALSE ),
                    smallest = sample(rownames(counts_sub), n.smallest, replace = FALSE ))
    pw.mean = calc.pwmeans(counts_sub, rnd.genes) %>% dplyr::select(-PathwayID) %>% t
    correlation.i.b = cor(pw.mean[,1], cd4.ifn, method = "spearman")
    correlation.a.b = cor(pw.mean[,1], cd4.any, method = "spearman")
    correlation.i.s = cor(pw.mean[,2], cd4.ifn, method = "spearman")
    correlation.a.s = cor(pw.mean[,2], cd4.any, method = "spearman")

    pw.mean = calc.pwmeans(counts_sub[rownames(vir.vector)], rnd.genes) %>% dplyr::select(-PathwayID) %>% t
    correlation.v.b = cor(pw.mean[,1], vir.vector, method = "spearman")
    correlation.v.s = cor(pw.mean[,2], vir.vector, method = "spearman")

    cor.tib[paste0('iter_', iter)] = c(correlation.i.s,correlation.a.s, correlation.v.s,
                                        correlation.i.b,correlation.a.b, correlation.v.b)

}

dir.create(paste0("result_tables/",tissue_type,"/bootstrap_pw/"), showWarnings = FALSE)
write.table(cd4.tib, paste0("result_tables/",tissue_type,"/bootstrap_pw/bootstrapped_cd4any.csv"),
            row.names = FALSE, sep = '\t')
write.table(cd4i.tib, paste0("result_tables/",tissue_type,"/bootstrap_pw/bootstrapped_cd4ifng.csv"),
            row.names = FALSE, sep = '\t')
write.table(vir.tib, paste0("result_tables/",tissue_type,"/bootstrap_pw/bootstrapped_vir.csv"),
            row.names = FALSE, sep = '\t')


cor.int = get.interval(cor.tib, remove.na = TRUE)
cor.int['spearman.rho'] = c(0,0,0)

dir.create(paste0("result_tables/Random/"), showWarnings = FALSE)
cor.int %>% dplyr::select('spearman.rho', everything()) %>% mutate(across(where(is.numeric), round, 4)) %>%
        rownames_to_column('Pathway') %>%
        write.table(paste0("result_tables/Random/RandomPWCorrelation_ci.csv"), row.names = FALSE, sep = '\t')

# get confidence interval
cd4.int = get.interval(cd4.tib)
colnames(cd4.int) = c('cd4any.conf.up','cd4any.conf.down')

cd4i.int = get.interval(cd4i.tib)
colnames(cd4i.int) = c('cd4ifng.conf.up','cd4ifng.conf.down')

vir.int = get.interval(vir.tib)
colnames(vir.int) = c('vir.conf.up','vir.conf.down')


pw.tib = vroom( paste0("result_tables/",tissue_type,'/PathwayCorrelation.csv'))
pw.tib = pw.tib %>%
    right_join(cd4i.int %>% rownames_to_column('PathwayID'), by = 'PathwayID') %>%
    right_join(cd4.int %>% rownames_to_column('PathwayID'), by = 'PathwayID') %>%
    right_join(vir.int %>% rownames_to_column('PathwayID'), by = 'PathwayID') %>%
    dplyr::select('PathwayID', 'Description',
    starts_with('cd4ifng'), starts_with('cd4any'),starts_with('vir'))

pw.tib %>% mutate(across(where(is.numeric), round, 4)) %>%
          write.table( paste0("result_tables/",tissue_type,"/PathwayCorrelation_ci.csv"),
                        row.names = FALSE, sep = '\t')








# # export
# #annotate pathways
# pw.tib = annotatePW(c('/home/obaranov/projects/additionalData/geneSets/GeneSet_SignalingMolecules.gmt',
# '/home/obaranov/projects/additionalData/geneSets/GeneSet_SignalTransduction.gmt',
# '/home/obaranov/projects/additionalData/geneSets/GeneSet_Immune.gmt'))
#
# mean.df.full = mean.df.full %>% right_join(pw.tib, by = 'PathwayID') %>%
# dplyr::select('PathwayID', 'Description',
# ends_with('spearman.rho'), ends_with('p.value'), everything()) %>%
# mutate(across(where(is.numeric), round, 4))
#
#
# dir.create(paste0("OutputData/",tissue_type), showWarnings = FALSE)
# mean.df.full  %>% arrange(cd4ifng.spearman.rho) %>%
# write.table(paste0("OutputData/",tissue_type,"/PathwayCorrelation.csv"),
# row.names = FALSE, sep = '\t')

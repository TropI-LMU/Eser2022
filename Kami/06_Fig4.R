args <- commandArgs(trailingOnly = TRUE)

predir = ''
basedir = args[1]
plotdir = basedir
codedir = args[2]
tissue_type = args[3]
# predir = '/home/obaranov/projects/'
# basedir = 'Tabea_Paper/debug/'
# codedir = 'Tabea_Paper/code/'

library(tidyverse)
library(ggplot2)
library(vroom)
library(gage)
library(pheatmap)
library(matrixStats)
library(readxl)
library(patchwork)
library(grid)
library(extrafont)



setwd(paste0(predir,basedir))
source(paste0(predir,codedir,'/utility/Plot_Utility.R'))
source(paste0(predir,codedir,'/utility/General_Utility.R'))


#########################
### correlation plots ###
#########################
#note that multiple objects required for B and C are defined in A
labs = read_excel(paste0(predir,'/',codedir,'/utility/PW_characterisation.xlsx')) %>%
    dplyr::select(c('PathwayID','Label','Primary_group')) %>%
    rename(Primary_group = 'Group')
labs$Label = sapply(labs$Label, trimws) %>% unname

### panel A: CD4 correlation

tissue_type = 'CD4'

gene.tib = vroom(paste0(predir,basedir,'/result_tables/',tissue_type,'/PathwayCorrelation_ci.csv')) %>%
            column_to_rownames('PathwayID')
rand.tib = vroom(paste0(predir,basedir,'/result_tables/Random/RandomPWCorrelation_ci.csv')) %>%
            column_to_rownames('Pathway')

rand.row = list()

new.row = list()
for(row in rownames(rand.tib)){
    tmp.row = rand.tib[row,]
    colnames(tmp.row) = colnames(tmp.row) %>%
            sapply(function(x){paste0(row,'.',x)})
    new.row = new.row %>% append(tmp.row)
}
small.row = new.row[grepl('smallest', names(new.row))]
names(small.row) = str_replace(names(small.row),'smallest\\.','')
plot.tib = gene.tib %>%
    add_row(as_tibble(small.row))

big.row = new.row[grepl('biggest', names(new.row))]
names(big.row) = str_replace(names(big.row),'biggest\\.','')
plot.tib = plot.tib %>%
    add_row(as_tibble(big.row))
rownames(plot.tib) = c(rownames(gene.tib),'30.random', '300.random')

plot.tib['Pathway'] = rownames(plot.tib)
plot.tib['Pathway'] = factor(rownames(plot.tib),levels=rownames(plot.tib))
plot.tib['col'] = 'pathway'
plot.tib['30.random','col'] = 'random'
plot.tib['300.random','col'] = 'random'
sigs = plot.tib %>% filter(cd4ifng.p.value < 0.05) %>% rownames()
plot.tib[sigs,'col'] = 'significant'


plot.tib = plot.tib %>% rownames_to_column('PathwayID') %>%
            left_join(labs, by = 'PathwayID') %>%
            column_to_rownames('PathwayID')
plot.tib['Group'] = factor(plot.tib %>% pull('Group'),
        levels = c('Migration','Chemokines', 'innate immunity','T cell function','other','ctl'))

plot.tib['30.random', 'Label'] = "random 30 genes"
plot.tib['300.random', 'Label'] = "random 300 genes"
plot.tib['30.random', 'Group'] = "ctl"
plot.tib['300.random', 'Group'] = "ctl"
group.pal = get.palette.eqdist(plot.tib %>% pull('Group'),'Set1')
names(group.pal) = c('Migration','Chemokines','innate immunity','T cell function')
group.pal['other'] ='#AAAAAA'
group.pal['ctl'] ='#FFFFFF'
group.pal['pathway'] ='grey'
group.pal['significant'] ='black'
group.pal['random'] ='red'

# following pathways are non-species specific, having its human equivalent in the data
plot.tib = plot.tib[! rownames(plot.tib) %in% c('hsa04392','hsa04672'),]
plot.tib = plot.tib %>% arrange(by = cd4ifng.spearman.rho)
plot.tib['Pathway'] = factor(rownames(plot.tib), levels=rownames(plot.tib))
plot.tib = plot.tib %>% drop_na(Label)

plot.tib =plot.tib %>% arrange(cd4ifng.spearman.rho)
plot.tib['Pathway'] = factor(rownames(plot.tib), levels=rownames(plot.tib))

tmp.tib = plot.tib  %>% filter((cd4ifng.spearman.rho > 0.25 | cd4ifng.spearman.rho < -0.25) | grepl('random',Pathway))

# loadfonts(device = "postscript")

loadfonts()
p = ggplot(tmp.tib, aes(y = cd4ifng.spearman.rho, x = Pathway, colour = col)) +
    geom_point() +
    geom_errorbar(aes(ymax = cd4ifng.conf.up, ymin = cd4ifng.conf.down)) +
    geom_hline(yintercept=0, color = 'black') +
    theme_minimal(base_size = 15) +
    geom_point(aes(y = rep(-1.05, dim(tmp.tib)[1]), colour = Group),
                    shape=15, size=3)  +
    scale_colour_manual(values = group.pal) +
    scale_x_discrete(labels=tmp.tib$Label) + ylim(c(-1.05,1)) +
    ylab('Spearman correlation') + xlab('') +
    theme(text=element_text(family="Arial"),
        axis.text.x = element_text(angle = 90,hjust=1,vjust=0.2),
        legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="grey"))

ggsave(paste0(predir,basedir, '/newplots/A_', tissue_type,'_correlation_poster.svg'), plot = p,
        width = 16.5, height = 18, unit = 'cm')
tmp.tib %>% 
    dplyr::select(c("Description", "cd4ifng.spearman.rho", "cd4ifng.p.value","cd4ifng.conf.up","cd4ifng.conf.down")) %>%
    rownames_to_column('TermID') %>%
    write_delim(paste0(predir, basedir, "/result_tables/export/Fig4PanelB.csv"), delim = "\t")


### panel C: virus correlation
plot.tib <- plot.tib %>% arrange(vir.spearman.rho)
plot.tib["Pathway"] <- factor(rownames(plot.tib), levels = rownames(plot.tib))
plot.tib["col"] <- "pathway"
plot.tib["30.random", "col"] <- "random"
plot.tib["300.random", "col"] <- "random"
sigs <- plot.tib %>%
    filter(vir.p.value < 0.05) %>%
    rownames()
plot.tib[sigs, "col"] <- "significant"


plot.tib["Group"] <- factor(plot.tib %>% pull("Group"),
    levels = c("Migration", "Chemokines", "innate immunity", "T cell function", "other", "ctl")
)

tmp.tib <- plot.tib %>% filter((vir.spearman.rho > 0.25 | vir.spearman.rho < -0.25) | grepl("random", Pathway))
ggplot(tmp.tib, aes(y = vir.spearman.rho, x = Pathway, colour = col)) +
    geom_point() +
    geom_errorbar(aes(ymax = vir.conf.up, ymin = vir.conf.down)) +
    geom_hline(yintercept = 0, color = "black") +
    theme_minimal(base_size = 15) +
    geom_point(aes(y = rep(-1.05, dim(tmp.tib)[1]), colour = Group),
        shape = 15, size = 3
    ) +
    scale_colour_manual(values = group.pal) +
    scale_x_discrete(labels = tmp.tib$Label) +
    ylim(c(-1.05, 1)) +
    ylab("") +
    xlab("") +
    ylim(c(-1.05, 1)) +
    scale_y_continuous(position = "right") +
    theme(
        text = element_text(family = "Arial"), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
        legend.position = "none",
        axis.ticks.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = .1, color = "grey")
    )

ggsave(paste0(predir, basedir, "/newplots/C_virusload_correlation_poster.svg"),
    width = 11.5, height = 18, unit = "cm"
)

tmp.tib %>%
    dplyr::select(c("Description","vir.spearman.rho", "vir.p.value", "vir.conf.up", "vir.conf.down")) %>%
    rownames_to_column('TermID') %>%
    write_delim(paste0(predir, basedir, "/result_tables/export/Fig4PanelD.csv"), delim = "\t")



### panel B: cd8 correlation
tissue_type = 'CD8'

gene.tib = vroom(paste0(predir,basedir,'/result_tables/',tissue_type,'/PathwayCorrelation_ci.csv')) %>%
            column_to_rownames('PathwayID')
rand.tib = vroom(paste0(predir,basedir,'/result_tables/Random/RandomPWCorrelation_ci.csv')) %>%
            column_to_rownames('Pathway')

rand.row = list()

new.row = list()
for(row in rownames(rand.tib)){
    tmp.row = rand.tib[row,]
    colnames(tmp.row) = colnames(tmp.row) %>%
            sapply(function(x){paste0(row,'.',x)})
    new.row = new.row %>% append(tmp.row)
}
small.row = new.row[grepl('smallest', names(new.row))]
names(small.row) = str_replace(names(small.row),'smallest\\.','')
plot.tib = gene.tib %>%
    add_row(as_tibble(small.row))

big.row = new.row[grepl('biggest', names(new.row))]
names(big.row) = str_replace(names(big.row),'biggest\\.','')
plot.tib = plot.tib %>%
    add_row(as_tibble(big.row))
rownames(plot.tib) = c(rownames(gene.tib),'30.random', '300.random')

plot.tib['Pathway'] = rownames(plot.tib)
plot.tib['Pathway'] = factor(rownames(plot.tib),levels=rownames(plot.tib))
plot.tib['col'] = 'pathway'
plot.tib['30.random','col'] = 'random'
plot.tib['300.random','col'] = 'random'
sigs = plot.tib %>% filter(cd8ifng.p.value < 0.05) %>% rownames()
plot.tib[sigs,'col'] = 'significant'

plot.tib = plot.tib %>% rownames_to_column('PathwayID') %>%
            left_join(labs, by = 'PathwayID') %>%
            column_to_rownames('PathwayID')
plot.tib['Group'] = factor(plot.tib %>% pull('Group'),
            levels = c('Migration','Chemokines', 'innate immunity','T cell function','other','ctl'))

plot.tib['30.random', 'Label'] = "random 30 genes"
plot.tib['300.random', 'Label'] = "random 300 genes"
plot.tib['30.random', 'Group'] = "ctl"
plot.tib['300.random', 'Group'] = "ctl"

plot.tib = plot.tib[! rownames(plot.tib) %in% c('hsa04392','hsa04672'),]
plot.tib = plot.tib %>% arrange(by = cd8ifng.spearman.rho)
plot.tib['Pathway'] = factor(rownames(plot.tib), levels=rownames(plot.tib))
plot.tib = plot.tib %>% drop_na(Label)

plot.tib =plot.tib %>% arrange(cd8ifng.spearman.rho)
plot.tib['Pathway'] = factor(rownames(plot.tib), levels=rownames(plot.tib))


tmp.tib = plot.tib  %>% filter((cd8ifng.spearman.rho > 0.25 | cd8ifng.spearman.rho < -0.25) | grepl('random',Pathway))
ggplot(tmp.tib, aes(y = cd8ifng.spearman.rho, x = Pathway, colour = col)) +
    geom_point() +
    geom_errorbar(aes(ymax = cd8ifng.conf.up, ymin = cd8ifng.conf.down)) +
    geom_hline(yintercept=0, color = 'black') +
    theme_minimal(base_size = 15) +
    geom_point(aes(y = rep(-1.05, dim(tmp.tib)[1]), colour = Group),
                    shape=15, size=3) +
    scale_colour_manual(values = group.pal) +
    scale_x_discrete(labels=tmp.tib$Label) +
    ylab('') + xlab('') +
    ylim(c(-1.05,1)) +
    theme(text=element_text(family="Arial"),axis.text.x = element_text(angle = 90,hjust=1,vjust=0.2),
        legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="grey"))

ggsave(paste0(predir, basedir, '/newplots/B_', tissue_type,'_correlation_poster.svg'),
        width = 11.5, height = 18, unit = 'cm')

tmp.tib %>%
    dplyr::select(c("Description", "cd8ifng.spearman.rho", "cd8ifng.p.value", "cd8ifng.conf.up", "cd8ifng.conf.down")) %>%
    rownames_to_column('TermID') %>%
    write_delim(paste0(predir, basedir, "/result_tables/export/Fig4PanelC.csv"), delim = "\t")


### panel D: Gene correlations

# last minutes changes, as usual
rearrange = function(gene, tib, tissue){
    outtib = tibble(cd8ifn = c(tib[gene,'cd8ifng.spearman_rho'], tib[gene, 'cd8ifng.p_value'],
    tib[gene, 'cd8ifng.lower'], tib[gene, 'cd8ifng.upper']),
    cd4ifn = c(tib[gene,'cd4ifng.spearman_rho'], tib[gene, 'cd4ifng.p_value'],
    tib[gene, 'cd4ifng.lower'], tib[gene, 'cd4ifng.upper']),
    vir = c(tib[gene,'vir.spearman_rho'], tib[gene, 'vir.p_value'],
    tib[gene, 'vir.lower'], tib[gene, 'vir.upper']))
    outtib = outtib %>% bind_rows( as.data.frame(outtib[2,] < 0.05) )
    rownames(outtib) = c('rho','pvalue','lower','upper','sig')
    return(outtib %>% t() %>% as.data.frame)
}


plotforgene = function(gene){
    gene.subtib = NULL
    for (tissue in c('NK','CD4','CD8','Monocytes')){
        gene.tib = vroom(
            paste0(predir,basedir,'/result_tables/',tissue,'/GeneCorrelation_ci.csv')) %>%
        column_to_rownames('gene')
        tmp.tib = rearrange(gene, gene.tib, tissue) %>%
        rownames_to_column('correlate') %>%
        mutate(tissue = tissue)
        gene.subtib = gene.subtib %>% bind_rows(tmp.tib)
    }

    symbol = gene.tib[gene,'symbol']
    gene.subtib$tissue = factor(gene.subtib$tissue, levels = c('Monocytes','NK','CD4','CD8'))
    gene.subtib$correlate = factor(gene.subtib$correlate, levels = c('vir','cd8ifn','cd4ifn'))
    gene.subtib$sig = factor(gene.subtib$sig, levels = c(1,0))
    # gene.subtib$size = ifelse(gene.subtib$pvalue > 0.05, 'nonsig','sig')
    Cell_type = c(Monocytes = '#8dd3c7',NK ='#0f804d',CD4 = '#55a0fb',CD8='#5757f9')
    p = ggplot(gene.subtib %>% arrange(by = desc(tissue),correlate ) %>% mutate(dummy = 1:dim(gene.subtib)[1]) ,
            aes(y = rho, x = dummy, colour = tissue, linetype = sig)) +
            geom_point(aes(shape = correlate), size = 3) +
            geom_errorbar(aes(ymax = upper, ymin = lower)) +
            geom_hline(yintercept=0, color = 'black') +
            theme_minimal(base_size = 15) +
            theme(text=element_text(family="Arial")) +
            scale_x_discrete(labels=gene.subtib$correlate) + ylim(c(-1.05,1)) +
            coord_flip() +
            scale_colour_manual(values = Cell_type) +
            ylab('Spearman correlation') + xlab('') +
            ggtitle(symbol)
    gene.subtib %>% arrange(by = desc(tissue),correlate ) %>%
        write_delim(paste0(predir, basedir, "/result_tables/export/PanelA_", symbol ,".csv"), delim = "\t")
    return(p)
}

genes = c( "ENSG00000089127", # oas
"ENSG00000115415", # stat1
"ENSG00000055332") # prk

uberp = (plotforgene(genes[1]) + theme(text=element_text(family="Arial"),legend.position = "none")) +
(plotforgene(genes[2]) + theme(text=element_text(family="Arial"),legend.position = "none")) +
plotforgene(genes[3])

ggsave(paste0(predir, basedir,'/newplots/Z_Gene_correlation.svg'), plot = uberp, width = 30, height = 8, unit = 'cm')

###############
### heatmap ###
###############


#gather significant tissue - pathway combinations
tis.pw.keep = tibble()
for(tissue_type in c('Monocytes','NK','CD4','CD8')){
    pw.up = vroom(paste0(predir, basedir,'/result_tables/',tissue_type,'/PathwayEnrichment_Up.csv')) %>%
            arrange(by= PathwayID)
    pw.down = vroom(paste0(predir, basedir,'/result_tables/',tissue_type,'/PathwayEnrichment_Down.csv'))%>%
            arrange(by= PathwayID)

    subtib = bind_rows(pw.up %>% filter(p.val < 0.05) %>% mutate(direction = 'up'),
              pw.down %>% filter(p.val < 0.05)%>% mutate(direction = 'down')) %>%
              mutate(tissue_type = tissue_type)
    tis.pw.keep = tis.pw.keep %>% bind_rows(subtib)
}

mc_countlist = c('Monocytes','NK','CD4','CD8') %>%
                lapply(function(x){read.csv(paste0(predir,basedir,"/data/bulkRNA/rawData/",x,"_counts.csv"),
                                    sep = '\t', check.names=FALSE) %>% column_to_rownames('Gene')} )
names(mc_countlist) = c('Monocytes','NK','CD4','CD8')


virload = read.csv(paste0(predir,basedir,"/data/bulkRNA/rawData/Virusload.csv"),
                            header = TRUE, stringsAsFactors = FALSE, sep = '\t') %>% drop_na
cd4ifn.tib = vroom(paste0(predir,basedir,'/data/bulkRNA/rawData/CellFractions.csv')) %>%
                filter(stimulant == 'nucleocapsid') %>%
                column_to_rownames('Sample_ID')
colnames(cd4ifn.tib) = colnames(cd4ifn.tib) %>%
                    sapply(function(x){str_replace(x,' \\| Freq. of Parent','') %>% str_replace(' ','_')})
# create metadata subtable
metadata <- read.csv(paste0(predir,basedir,"/data/bulkRNA/rawData//Metadata.csv"), header = TRUE,
                            stringsAsFactors = FALSE, sep = '\t')
metadata_mc <- metadata[metadata$Tissue  == tissue_type,]
metadata_mc['Sample_ID'] = metadata_mc %>% pull('Sample_ID') %>% sapply(function(x){as.character(x)})


#read gene sets
list.immu = readList(paste0(predir,codedir,'/utility/PathwayLists/GeneSet_Immune.gmt'))
list.sigtrans = readList(paste0(predir,codedir,'/utility/PathwayLists/GeneSet_SignalingMolecules.gmt'))
list.sigmol = readList(paste0(predir,codedir,'/utility/PathwayLists/GeneSet_SignalTransduction.gmt'))
kegg.set = list.immu %>% append(list.sigtrans) %>% append(list.sigmol)
pw.clas = read_excel(paste0(predir,codedir,'/utility/PW_characterisation.xlsx')) %>%
            dplyr::select(-'link') %>%
            dplyr::rename(Group = Primary_group)

# get the group
metadata_mc_d4 = metadata_mc %>% filter(TAG == 'Tag4')
mc_counts_sub = c('Monocytes','NK','CD4','CD8') %>%
            lapply(function(x){as.matrix(mc_countlist[[x]][,metadata_mc_d4 %>% pull('Sample_ID')])})
names(mc_counts_sub) = c('Monocytes','NK','CD4','CD8')
keep_tabi = cd4ifn.tib %>% filter(stimulant == 'nucleocapsid') %>%
                filter(weeks_after_onset_of_symptoms == 1) %>%
                rownames_to_column('Sample_ID') %>%
                pull('Sample_ID')
keepsamp = intersect(keep_tabi, metadata_mc_d4 %>% pull('Sample_ID'))
mc_counts_sub = mc_counts_sub %>% lapply(function(x){x[,keepsamp]})

# convert the gene IDs
for(tissue_type in names(mc_counts_sub)){
    # print(tissue_type)
    mc_counts_sub[[tissue_type]] = EnsToEnt(mc_counts_sub[[tissue_type]])
}

mean.df.full = mc_counts_sub %>% lapply( function(x){calc.pwmeans(x, kegg.set) %>%
                                                column_to_rownames('PathwayID')} )


# building the plot tibble
mean.exp.tib = tibble()
idx = 1
for(idx in 1:dim(tis.pw.keep)[1]){
    pw = tis.pw.keep[idx,'PathwayID'] %>% pull
    tt = tis.pw.keep[idx,'tissue_type'] %>% pull
    row = mean.df.full[[tt]][pw,]
    rownames(row) = paste0(pw, '_', tt)
    mean.exp.tib = mean.exp.tib %>% bind_rows(row)
}


# adding PW enrichment p values to the plot

anno.sam = cd4ifn.tib[keepsamp,] %>%
            rownames_to_column('Sample_ID') %>%
            dplyr::select(c('Sample_ID','CD4IFNg','CD8IFNg')) %>%
            mutate(CD4pos = as.numeric(CD4IFNg > 0.01),
                   CD8pos = as.numeric(CD8IFNg > 0.01))
anno.sam['CD4pos'] = as.factor(anno.sam %>% pull('CD4pos'))
anno.sam['CD8pos'] = as.factor(anno.sam %>% pull('CD8pos'))
anno.sam = anno.sam %>% left_join(virload %>% dplyr::select(-Patient_ID), by = 'Sample_ID')
anno.sam['VirusLoad'] = log(anno.sam['VirusLoad'])

anno.pw = tibble(PathwayID = mean.exp.tib %>% rownames()) %>%
            column_to_rownames('PathwayID')
anno.pw['PathwayID'] = rownames(anno.pw) %>% sapply(function(x){str_split(x,'_')[[1]][1]})

anno.pw = anno.pw %>% rownames_to_column('ID') %>%
        left_join((pw.clas %>% dplyr::select(-'Description')), by = 'PathwayID') %>%
        column_to_rownames('ID')
anno.pw['Group'] = factor(anno.pw %>% pull('Group'), levels = c('Migration','Chemokines', 'innate immunity','T cell function','other'))
anno.pw['Secondary_group'] = factor(anno.pw %>% pull('Secondary_group'), levels = c('Migration','Chemokines', 'innate immunity','T cell function','other','none'))
anno.pw['Tertiary_group'] = factor(anno.pw %>% pull('Tertiary_group'), levels = c('Migration','Chemokines', 'innate immunity','T cell function','other','none'))
anno.pw['Cell_type'] = rownames(anno.pw) %>%
                        sapply(function(x){str_split(x, '_')[[1]][2]}) %>%
                        factor(levels = c('Monocytes','NK','CD4','CD8'))
anno.pw = anno.pw %>% drop_na('Label')
roworder = anno.pw %>% arrange(Cell_type,Group) %>% rownames
# excluding the replicate of Hippo and the intestinal pw
roworder = roworder[! grepl(c('(hsa04392)|(hsa04672)'), roworder)]

# z score increases the interpretability of the plot
plot.tib = log(mean.exp.tib)
zscore = function(x){(x - mean(x)) / sd(x)}

tmp.tib = plot.tib %>% t %>% as_tibble %>%
        mutate_all(zscore) %>%
        t %>% as.tibble
colnames(tmp.tib) = colnames(plot.tib)
tmp.tib['PathwayID_ct'] = rownames(mean.exp.tib)
tmp.tib['PathwayID'] = tmp.tib %>% pull('PathwayID_ct') %>%
            sapply(function(x){strsplit(x,'_')[[1]][1]}) %>%
            unname()

tmp.tib = pw.up %>% dplyr::select('PathwayID','PathwayName') %>%
    right_join(tmp.tib, by = 'PathwayID') %>%
    column_to_rownames('PathwayID_ct')
labs = anno.pw[roworder,'Label']

# ann_colors = palette.list(anno.sam, anno.pw, palettes = c(enrichment = 'RdBu'))
ann_colors = list(CD4pos = c('0'='white', '1'='darkorchid4'),
                  CD8pos = c('0'='white', '1'='darkorchid1'),
                    VirusLoad = get.palette.numeric(anno.sam %>% pull('VirusLoad'),'Oranges'),
                    Group = get.palette.eqdist(anno.pw %>% pull('Group'),'Set1'),
                    Secondary_group = get.palette.eqdist(anno.pw %>% pull('Group'),'Set1'),
                    Tertiary_group = get.palette.eqdist(anno.pw %>% pull('Group'),'Set1'),
                    Cell_type = c(get.palette.eqdist(anno.pw %>% pull('Cell_type'), 'Set3')[1],NK ='#0f804d',CD4 = '#55a0fb',CD8='#5757f9')
                )
ann_colors[['Group']]['other'] = "#AAAAAA"

tmp.tib %>% drop_na() %>%
    rownames_to_column("PathwayID_celltype") %>%
    dplyr::select( - PathwayID) %>% 
    write_delim(paste0(predir, basedir, "/result_tables/export/Fig4PanelE_heatmap.csv"), delim = "\t")

tmp.tib  = tmp.tib %>% dplyr::select(- PathwayName) %>% drop_na()

# ugly hack to rotate labels
draw_colnames_90 <- function (coln, gaps, ...) {
    coord = pheatmap:::find_coordinates(length(coln), gaps)
    x = coord$coord - 0.5 * coord$size
    res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
    return(res)}

## 'Overwrite' default draw_colnames with your own version

draw_colnames_45 <- function (coln, gaps, ...) {
    coord = pheatmap:::find_coordinates(length(coln), gaps)
    x = coord$coord - 0.5 * coord$size
    res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
    return(res)}

## 'Overwrite' default draw_colnames with your own version
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
ns=asNamespace("pheatmap"))

horiz = pheatmap(tmp.tib[roworder,] %>% dplyr::select(-PathwayID) %>% t,
                labels_col = labs %>% unname,
                annotation_row = ( anno.sam %>% column_to_rownames('Sample_ID') %>%
                                    dplyr::select(-'CD4IFNg',-'VirusLoad', -'CD8IFNg') ),
                annotation_col =  anno.pw[roworder,] %>% dplyr::select('Group','Cell_type'),
                annotation_colors = ann_colors,
                cluster_cols = FALSE,
                cluster_rows = TRUE,
                show_rownames = FALSE, fontsize = 9,
                cellheight=10, cellwidth = 12)


save_pheatmap_pdf(horiz,paste0(predir, basedir,'/newplots/D_HeatMap.pdf'), width = 15, height = 5)

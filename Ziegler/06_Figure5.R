library(tidyverse)
library(vroom)
require(gridExtra)
library(ggrepel)
library(cowplot)
library(ggalt)
library(stringr)
library("ggpubr")
theme_set(theme_bw())


args <- commandArgs(trailingOnly = TRUE)

basedir_z <- args[1]
basedir_y <- args[2]
codedir <- args[3]

setwd(basedir_z)
source(paste0(codedir,'/utility/Plot_Utility.R'))


fixCellNames = function(nlist){
    outlist = list()
    for (name in nlist){
        name = as.character(name) %>% strsplit(' ') %>% unlist
        if(length(name) >1){
        name = paste(c(name[1], name[2:length(name)] %>%
            sapply(function(x){ifelse(grepl("[[:lower:]]", x), tolower(x), x)})  %>%
            unlist() %>% unname() ),
            collapse = ' ')
        } else {name}
        outlist = append(outlist, name)
    }
    return(outlist)
}

fixPWNames = function(plist){
    for (idx in 1:length(plist)){
        plist[idx]=paste(toupper(substr(plist[idx], 1, 1)), substr(plist[idx], 2, nchar(plist[idx])), sep="")
    }
    return(plist)
}

# for the labeling we use manually selected gene set
labs = vroom(paste0(codedir, 'utility/VolcanoLabels.txt'))
palette = c('#4d1a17','#8DD3C7')
celltypes = c('Developing Ciliated Cells',
            'CCL5 high Squamous Cells',
            'Basal Cells')
#ziegler tibs
z.tib.deg = vroom(paste0(basedir_z,'/OutputData/markerGenes/Ziegler_markerGenes_allCellTypes.csv'))
# yosh tibs
y.tib.deg = vroom(paste0(basedir_y, "/OutputData/markerGenes/Yoshida_markerGenes_allCellTypes.csv"))

#common tib
go <- vroom(paste0(basedir_z, "/OutputData/Both_goEnrichment_Ciliated.csv"))
go["Label"] <- fixPWNames(go %>% pull("Term"))

ctz = "Developing Ciliated Cells"
cty = "Ciliated"

################################################################################
## panel A: go terms 
################################################################################


makeGo <- function(ct, tib.go) {
    tib.go["Combined.Score"] = tib.go[c("Combined.Score.x", "Combined.Score.y")] %>% 
                                    rowMeans()
    tib.go["Adjusted.P.value"] <- tib.go[c("Adjusted.P.value.x", "Adjusted.P.value.y")] %>%
        rowMeans()
    overmax.tmp = tib.go$Overlap.y %>% 
                            sapply(function(x){str_split(x,'/')[[1]][2] %>% 
                            as.numeric()}) %>% 
                            unname
    tib.go['Overlap.y'] = tib.go$Overlap.y %>% 
                            sapply(function(x){str_split(x,'/')[[1]][1] %>% 
                            as.numeric()}) %>% 
                            unname
    tib.go['Overlap.x'] = tib.go$Overlap.x %>% 
                            sapply(function(x){str_split(x,'/')[[1]][1] %>% 
                            as.numeric()}) %>% 
                            unname
    overlap.tmp <- tib.go[c("Overlap.x", "Overlap.y")] %>%
        rowMeans()
    tib.go['Overlap'] = seq_along(overlap.tmp) %>% sapply(function(x){paste0(
                        as.character(overlap.tmp[x]),'/',
                        as.character(overmax.tmp[x])
                        )})
    p <- pw.dotplot(tib.go,
        x.col = "Combined.Score", y.col = "Label", n.terms = 30
    ) +
        theme(text = element_text(family = "Arial", colour = "black")) +
        xlab("Combined score") +
        ylab("") +
        scale_color_gradient(
            low = "#4d1a17", high = "#8DD3C7",
            trans = "reverse"
        ) + theme_bw() +
        theme(
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black")
        )
    return(p)
}


go$Label = go$Term %>% sapply(parseGo) %>% unname
p_go <- makeGo(ctz, go) + ylab('enriched in both datasets')
ggsave(
    plot = p_go,
    filename = paste0(basedir_z, "/C_", cty, "_GOterms.svg"),
    width = 10, height = 7
)

################################################################################
## panel B: volcano plots
################################################################################

# get genes
genes.y = go$Genes.y %>% strsplit(';') %>% unlist %>% unique()
genes.z = go$Genes.x %>% strsplit(';') %>% unlist %>% unique()
genes.b = intersect(genes.y, genes.z)


# # add labels
# z.tib.deg.cil = z.tib.deg %>% filter(celltype == ctz)
# found.z = z.tib.deg.cil$Gene %>% 
#                 sapply(function(x){any(grepl(paste0('^',x,'$'),
#                                   genes.z))}) %>% unname
# found.bz = z.tib.deg.cil$Gene %>% 
#                 sapply(function(x){any(grepl(paste0('^',x,'$'),
#                                   genes.b))}) %>% unname
# z.tib.deg.cil["delabel"] = ''
# z.tib.deg.cil[found.z, "delabel"] = z.tib.deg.cil[found.z, "Gene"]
# z.tib.deg.cil['colour'] = 'black'
# z.tib.deg.cil[found.bz, "colour"] <- 'red'


# y.tib.deg.cil <- y.tib.deg %>% filter(celltype == cty)
# found.y = y.tib.deg.cil$Gene %>% 
#                 sapply(function(x){any(grepl(paste0('^',x,'$'),
#                                   genes.y))}) %>% unname
# found.by = y.tib.deg.cil$Gene %>% 
#                 sapply(function(x){any(grepl(paste0('^',x,'$'),
#                                   genes.b))}) %>% unname
# y.tib.deg.cil["delabel"] = ''
# y.tib.deg.cil[found.y, "delabel"] = y.tib.deg.cil[found.y, "Gene"]
# y.tib.deg.cil['colour'] = 'black'
# y.tib.deg.cil[found.y, "delabel"] <- y.tib.deg.cil[found.y, "Gene"]
# y.tib.deg.cil[found.by, "colour"] <- 'red'

y.tib.deg.cil = read.delim(paste0(basedir_y,"/OutputData/markerGenes/Yoshida_markerGenes_ciliated.csv"),
                            sep = '\t')
z.tib.deg.cil = read.delim(paste0(basedir_z,"/OutputData/markerGenes/Ziegler_markerGenes_ciliated.csv"),
                            sep = '\t')


sigy <- (y.tib.deg.cil$p_val_adj < 0.05) & (y.tib.deg.cil$avg_log2FC > 0.25 | y.tib.deg.cil$avg_log2FC < -0.25)
sigz <- (z.tib.deg.cil$p_val_adj < 0.05) & (z.tib.deg.cil$avg_log2FC > 0.25 | z.tib.deg.cil$avg_log2FC < -0.25)
z.tib.deg.cil[sigz,'diffexpressed'] = 'pos'
y.tib.deg.cil[sigy,'diffexpressed'] = 'pos'


makeVolc = function(ct, tib.deg, lab){
    ctlabs = labs[ct] %>%
        drop_na() %>%
        pull()

    tmptib = tib.deg %>%
        filter(celltype == ct) %>%
        mutate(chlab = ifelse(Gene %in% ctlabs, Gene, ""))

    p = ggplot(data = tmptib, aes(
        x = avg_log2FC, y = -log10(p_val),
        col = diffexpressed
    )) + 
        geom_point(show.legend = FALSE) +
        xlab("LFC") +
        ylab("-log10(p value)") +
        theme_classic(base_size = 12) +
        theme(text = element_text(family = "Arial",colour = 'black',),
        plot.title = element_text(hjust = 0.5),legend.position = "none") +
        geom_text_repel(aes(label = delabel, color = colour), max.overlaps = 1000,size = 3) +
        scale_color_manual(values = c(neg = palette[1], 
                                    pos = palette[2],
                                    black = 'black',
                                    red = 'red')) +
        ggtitle(lab)
    return(p)
}


p_vol_z = makeVolc(ctz, z.tib.deg.cil,'Ziegler - Dev. ciliated cells')
ggsave(plot = p_vol_z,
        filename = paste0(basedir_z,'/OutputData/B_',ctz,'_volcano.svg'),
        width = 10, height = 7)
p_vol_y = makeVolc(cty, y.tib.deg.cil, "Yoshida - Ciliated cells")
ggsave(plot = p_vol_y,
        filename = paste0(basedir_z,'/OutputData/B_',cty,'_volcano.svg'),
        width = 10, height = 7)


# ################################################################################
# ## merge into multipanel (E needs manual assembly)
# ################################################################################

uber = ggarrange(ggarrange(p_vol_z,p_vol_y, ncol = 2, labels = c('A','B')), p_go, nrow = 2, labels = c('','C'))


ggsave(paste0(basedir_z,'/OutputData/Fig5_v3.svg'), plot = uber,
    width = 300, height = 200, dpi = 300, units = "mm")
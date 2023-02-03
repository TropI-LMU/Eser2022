library(ggplot2)
library(tidyverse)
library(ggalt)
library(RColorBrewer)
library(pheatmap)
library(hablar)
library(extrafont)
library(showtext)

pw.dotplot = function(tab, x.col = "Odds.Ratio", y.col = 'Term',n.terms = 20, c.col = "Adjusted.P.value", s.col = "Overlap",
                      x.lab = '', y.lab = '', c.lab = '', s.lab = '', split.s.col = '/', y.len = 70){
  if (x.lab == '') x.lab = x.col
  if (y.lab == '') y.lab = y.col
  if (s.lab == '') s.lab = s.col
  if (c.lab == '') c.lab = c.col


  tab = tab %>% dplyr::filter(Adjusted.P.value < 0.05)
  if(dim(tab)[1] == 0){print('nothing to plot after filtering'); return(NA)}
  if(split.s.col != ''){
    tab['scol.split'] = tab %>% pull(s.col) %>%
                        sapply(function(x){strsplit(x,split.s.col)[[1]][1] %>% as.numeric})
  }
  levs = tab %>% arrange(desc(!! sym(x.col)) ) %>% pull(y.col)
  tab[y.col] = factor(tab %>% pull(y.col), levels = levs)

  p = ggplot(tab %>% arrange(desc(!! sym(x.col)) ) %>% dplyr::slice(1:n.terms)) +
    geom_point(aes_string(x=x.col, y = y.col, color = c.col, size = 'scol.split')) +
    labs(x = x.lab, y = y.lab, color = c.lab,  size = s.lab) +
    scale_y_discrete(limits=rev, label = function(x) stringr::str_trunc(x, y.len)) +
    scale_color_gradient(high = "#4d1a17", low = "#8ac2ff", trans = 'reverse' )

  p
}


parseGo = function(title,l = 35){
    title = strsplit(title, "\\(")[[1]]
    term = title[1] %>% substr(1,l)
    go =  strsplit(title[2],'\\)')[[1]][1]
    paste(term,go, sep ='...')
}


encircle = function(p, tab, x, y, group,filter.group = c(), size = 2){
    if(length(filter.group) > 0){
        tab = tab %>% filter(!!sym(group) %in% filter.group)
    }
    for(subgroup in tab %>% pull(group) %>% unique()){
        subtab = tab %>% filter(!!sym(group) == subgroup)
        p = p +geom_encircle(aes_string(x=x, y=y),
                             data=subtab, size=size,
                             expand=0.05, s_shape=1.8)
    }
    return(p)
}


palette.qual = function(col.uni, pal.to.use){
    pal = brewer.pal(length(col.uni),pal.to.use)[1:length(col.uni)]
    if(any(is.na(pal)) | (length(pal) < length(col.uni)) ){
        pal = pal[!is.na(pal)]
        repeatx=ceiling(length(col.uni) / length(pal))
        pal = rep(pal, repeatx)[1:length(col.uni)]
    }
    names(pal) = col.uni
    return(pal)
}

palette.seq.div = function(col.uni, pal.to.use){
    pal.fct = colorRamp(brewer.pal(3,pal.to.use))
    minv = sort(col.uni)[1]
    col.uni.norm = col.uni + abs(minv)
    maxv = sort(col.uni.norm, decreasing = TRUE)[1]
    col.uni.norm = col.uni.norm / maxv
    pal = pal.fct(col.uni.norm) %>%
            apply(1,convert.255RGB.to.HEX)
    names(pal) = col.uni
    return(pal)
}

get.palette.eqdist = function(vals, pal){
    vals.u = vals %>% sort %>% unique
    pal.col = brewer.pal(length(vals.u),pal)[1:length(vals.u)]
    names(pal.col) = vals.u
    return(pal.col)
}

get.palette.numeric = function(vals, pal){
    pal.fct = colorRamp(brewer.pal(7,pal))
    vals.u = vals %>% sort(na.last = FALSE) %>% unique
    if(any(is.na(vals.u))){
        vals.u = vals.u[-1]
        add.na = TRUE
    }
    vals.u.norm = vals.u + abs(min(vals.u))
    vals.u.norm = vals.u.norm / max(vals.u.norm)
    pal.col = pal.fct(vals.u.norm) %>%
            apply(1,convert.255RGB.to.HEX)
    names(pal.col) = vals.u
    if(add.na){
        nms = names(pal.col)
        pal.col = c(pal.col,'#000000')
        names(pal.col) = c(nms,'NA')
    }
    return(pal.col)
}

convert.255RGB.to.HEX = function(col){
    rgb(col[1] / 255,col[2] / 255,col[3] / 255)
}

#' creates a list with color palettes for every column in the given data frames;
#' assumes every column is a factor or can be treated as such
#'
#' @param ... any number of data frames
#' @param palettes a named list matching columns to specific palettes
#' @return named list of palettes for all columns in the dataframes
#' @examples
#' add(1, 1)
#' add(10, 1)
palette.list = function(... , palettes = NULL, drop.special = FALSE){
    dfs = list(...)
    pal.list = c()
    pals.avail = brewer.pal.info
    # handle cols with custm palettes first
    if(!is.null(palettes)){
        tmp.list = handle.special.cols(dfs, palettes, pal.list, pals.avail)
        pal.list = tmp.list[['pal.list']]
        dfs = tmp.list[['dfs']]
        unused = rownames(pals.avail) %>% setdiff(palettes)
        if(drop.special){
            pals.avail = pals.avail[unused,]
        } else {
            pals.avail = pals.avail[unused,] %>% rbind(pals.avail[palettes,])
        }
    }
    # index below is allways used with +1
    pal.idx = list(qual = 0, seq = 0, div = 0)
    pal.n = list(qual = pals.avail %>% dplyr::filter(category == 'qual') %>% dim() %>% `[`(1),
            seq = pals.avail %>% dplyr::filter(category == 'seq') %>% dim() %>% `[`(1),
            div = pals.avail %>% dplyr::filter(category == 'div') %>% dim() %>% `[`(1))
    for(df in dfs){
        for(col.idx in 1:dim(df)[2]){
            # class == numeric means continious var
            if(class(df[1,col.idx]) == 'factor'){
                col.val = df[col.idx] %>% pull
            } else {col.val = retype(df[col.idx])%>% pull}
            col.class = class(col.val)
            # switches the palette from seq to div if the values cross 0
            if( (col.class == 'integer') | (col.class == 'numeric') ){
                if ( min(col.val, na.rm=TRUE) < 0 ){
                    col.class = 'numeric.negative'
                }
            }
            pal.type = switch(col.class,'factor' = 'qual',
                            'character' = 'qual',
                            'numeric' = 'seq',
                            'integer' = 'seq',
                            'numeric.negative' = 'div')
            col.uni = col.val %>% sort %>% unique
            pal.to.use = pals.avail %>%
                        dplyr::filter(category == pal.type) %>%
                        rownames() %>%
                        `[`(pal.idx[[pal.type]] +1)
            # below, +1 fixes 1 based indexing
            pal.idx[[pal.type]] = (pal.idx[[pal.type]]+1) %% pal.n[[pal.type]]
            if(pal.type == 'seq' | pal.type == 'div'){
                if(length(col.uni %>% unique) > 1){
                    pal = palette.seq.div(col.uni, pal.to.use)
                } else { pal.type = 'qual' }
            }
            if(pal.type == 'qual') {
                pal = palette.qual(col.uni, pal.to.use)
            }
            pal.list[colnames(df)[col.idx]] = list(pal)
        }
    }
    return(pal.list)
}


handle.special.cols = function(dfs, palettes, pal.list, pals.avail){

    found = sapply(dfs, function(x){colnames(x) %>% intersect(names(palettes))})
    for(df.n in seq_along(found)){
        for(col in found[[df.n]]){
            col.uni = dfs[[df.n]] %>% pull(col) %>% unique
            pal.type = pals.avail[palettes[[col]],] %>% pull('category')
            if(pal.type == 'qual'){
                pal = palette.qual(col.uni, palettes[col])
                names(pal) = col.uni
            } else {
                if(any(is.na(col.uni))){
                    col.uni = sort(col.uni, na.last = FALSE)
                    pal = palette.seq.div(col.uni[-1], palettes[col])
                    pal = c('#000000',pal)
                    names(pal) = c('NA',col.uni[-1])
                } else {
                    pal = palette.seq.div(col.uni, palettes[col])
                    names(pal) = col.uni
                }
            }
            pal.list[col] = list(pal)
        }
        warning(paste0('dimension before delete: ', dim(dfs[[df.n]])))
        dfs[[df.n]] = dfs[[df.n]] %>% dplyr::select(-found[[df.n]])
        warning(paste0('dimension after delete: ', dim(dfs[[df.n]])))
    }
    return(list(dfs = dfs, pal.list = pal.list))
}


save_pheatmap_pdf <- function(x, filename, width=5, height=5) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}

# Theme closely resembling the default graphpad output
theme_graphpad <- function(){
    font <- "sans"   #assign font family up front

    theme_classic() %+replace%    #replace elements we want to change

    theme(

      #grid elements
      panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines

      #since theme_minimal() already strips axis lines,
      #we don't need to do that again

      #text elements
      plot.title = element_text(             #title
                   family = font)    ,        #set font family
                   # size = 20,                #set font size
                   # face = 'bold',            #bold typeface
                   # hjust = 0,                #left align
                   # vjust = 2),               #raise slightly

      plot.subtitle = element_text(          #subtitle
                   family = font,            #font family
                   size = 14),               #font size

      plot.caption = element_text(           #caption
                   family = font,            #font family
                   size = 9,                 #font size
                   hjust = 1),               #right align

      axis.title = element_text(             #axis titles
                   family = font,            #font family
                   size = 12),               #font size

      axis.text = element_text(              #axis text
                   family = font,            #axis famuly
                   size = 10,                 #font size
                   color = 'black'),          #font color

      axis.text.x = element_text(            #margin for axis text
                    margin=margin(5, b = 10))

      #since the legend often requires manual tweaking
      #based on plot content, don't define it here
    )
}

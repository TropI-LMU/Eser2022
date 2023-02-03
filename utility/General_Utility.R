
library(org.Hs.eg.db)
library(tidyverse)


get.annotation <- function(genes, keytype, species, get.columns,  merge = 'none'){
  del.ENS = FALSE
  if(keytype == 'SYMBOL'){
      library(limma)
      normed.genes = sapply(genes,function(x){alias2Symbol(x)[1]})
      genes = names(normed.genes) %>%
              sapply(function(x){if (is.na(normed.genes[[x]])) x else normed.genes[[x]][1]})
  }
  if(! 'ENSEMBL' %in% get.columns & keytype!='ENSEMBL'){get.columns <- c(get.columns,'ENSEMBL'); del.ENS <- TRUE}
  if(species == 'mouse' | species == 'mmu'){
    library(org.Mm.eg.db)
    ano.data <- AnnotationDbi::select(org.Mm.eg.db, keys = genes, columns = get.columns, keytype = keytype)
  } else if (species == 'human' | species == 'hsa') {
    library(org.Hs.eg.db)
    ano.data <- AnnotationDbi::select(org.Hs.eg.db, keys = genes, columns = get.columns, keytype = keytype)
  }
  # implement a warning for one to many mapping
  if(keytype == 'SYMBOL'){
      genes.dct = names(normed.genes)
      names(genes.dct) = genes
      ano.data$origSYMBOL = ano.data$SYMBOL %>%
                      sapply(function(x){genes.dct[x]}) %>% unname
  }
  ano.data <- merge.annotation(ano.data, merge)
  if(del.ENS){ano.data <- ano.data %>% dplyr::select(-'ENSEMBL')}
  return(ano.data)
}


EnsToEnt = function(mc_counts){
    ensids = rownames(mc_counts) %>%
            sapply(function(x){strsplit(x, '\\.')[[1]][1]})
    mc_counts = as_tibble(mc_counts)
    mc_counts['ensembl'] = ensids

    ezid = mc_counts %>%
            pull('ensembl') %>%
            get.annotation( 'ENSEMBL','human','ENTREZID')

    mc_counts = mc_counts %>%
            inner_join(ezid, by = c('ensembl'='ENSEMBL')) %>%
            dplyr::select(-'ensembl') %>%
            drop_na(ENTREZID) %>%
            distinct(ENTREZID, .keep_all = TRUE) %>%
            column_to_rownames('ENTREZID')
}

merge.annotation <- function(ano.data, merge){
  if(merge == 'first'){
    dups <- duplicated(ano.data[,1])
    if(sum(dups >0 )){warning(paste0('multiple annotations found for: ', paste(unlist(ano.data[dups,1]), collapse=' ') ))}
    ano.data <- ano.data[! dups, ]
    # ano.data <- ano.data %>% group_by(!! sym(keytype)) %>% filter(row_number()==1) for a dplyrry way
  } else if (merge == 'all'){
    warning(paste0('WARNING: using merge = all joins annotations by semicolon, you might want to review the data before using the table for further analysis;',
                   'multiple annotations found for: ', paste(unlist(ano.data[dups,1]), collapse=' ')))
    merge.by.semicolon <- function(clmn){paste(clmn, collapse=';')}
    ano.data <- ano.data %>% group_by(!! sym(keytype)) %>% mutate_at(vars(-group_cols()), merge.by.semicolon) %>% ungroup() %>%  unique()
  } else if(merge == 'drop') {
    dups <- duplicated(ano.data[,1])
    dup.genes <- unlist(ano.data[which(dups),1])
    if(sum(dups >0 )){warning(paste0('multiple annotations found for: ', paste(unlist(ano.data[dups,1]), collapse=' '), '; dropping the genes', collapse=' ' ))}
    ano.data <- ano.data  %>% filter(! ENSEMBL %in% dup.genes)
  }
  return(as_tibble(ano.data))
}

calc.pwmeans = function(counts, gene.list){
    mean.list = list()
    for(name in names(gene.list)){
            pw.mean = list(counts[gene.list[[name]],] %>%
                    drop_na() %>%
                    dplyr::select(where(is_numeric)) %>%
                    colMeans)
            names(pw.mean) = name
            mean.list = mean.list %>% append(pw.mean)
    }
    # has to be df to get the rownames
    mean.df = mean.list %>%
            as.data.frame %>%
            t %>%
            as.data.frame %>%
            rownames_to_column('PathwayID') %>%
            as_tibble()
            #   ⠀⠀⠀⠀⠀⢀⡔⠋⢉⠩⡉⠛⠛⠛⠉⣉⣉⠒⠒⡦⣄
            # ⠀⠀⠀⠀⠀⢀⠎⠀⠀⠠⢃⣉⣀⡀⠂⠀⠀⠄⠀⠀⠀⠀⢱⠀⠀⠀⠀⠀⠀ ⠀⠀⠀⠀⠀
            #     ⡰⠟⣀⢀⣒⠐⠛⡛⠳⢭⠆⠀⠤⡶⠿⠛⠂⠀⢈⠳⡀
            # ⠀⠀⢸⢈⢘⢠⡶⢬⣉⠉⠀⠀⡤⠄⠀⠀⠣⣄⠐⠚⣍⠁⢘⡇⠀⠀⠀⠀ ⠀⠀⠀⠀
            #   ⠈⢫⡊⠀⠹⡦⢼⣍⠓⢲⠥⢍⣁⣒⣊⣀⡬⢴⢿⠈⡜⠀⠀⠀⠀⠀ you should be mad!
            #     ⠹⡄⠀⠘⢾⡉⠙⡿⠶⢤⣷⣤⣧⣤⣷⣾⣿⠀⡇⠀⠀⠀⠀⠀ ⠀⠀⠀⠀⠀⠀⠀
            #      ⠘⠦⡠⢀⠍⡒⠧⢄⣀⣁⣀⣏⣽⣹⠽⠊⠀⡇⠀⠀⠀⠀⠀ ⠀⠀⠀⠀⠀⠀⠀⠀⠀
            #        ⠈⠑⠪⢔⡁⠦⠀⢀⡤⠤⠤⠄⠀⠠⠀⡇⠀⠀⠀⠀⠀ ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
            #             ⠈⠑⠲⠤⠤⣀⣀⣀⣀⣀⠔⠁
    return(mean.df)
}

calc.spearman.pw = function(pw.ex.tab,second.var,
                    draw.plot = NULL){
    sam.order = rownames(second.var)
    ps = c()
    rhos = c()
    p = list()
    if(! is_null(draw.plot)){pdf(draw.plot)}
    for(name in pw.ex.tab %>% pull('PathwayID')){
        b = pw.ex.tab[pw.ex.tab['PathwayID'] == name,] %>%
            dplyr::select(all_of(sam.order)) %>% t
        cor.res = cor.test(second.var %>% pull, b[sam.order,], method = 'spearman')
        ps = c(ps, cor.res$p.value)
        rhos = c(rhos,cor.res$estimate)
        p = plot.xvy(second.var %>% pull, b[sam.order,], ylabel = name,
                        rho = cor.res$estimate, p = cor.res$p.value)
        if(! is_null(draw.plot)){print(p)}
    }
    if(! is_null(draw.plot)){dev.off()}
    return(list(p.values = ps, rhos = unname(rhos)))
}

plot.xvy = function(x, y, xlabel = '', ylabel = '',
rho = NULL, p = NULL, trend = FALSE){
    plt = ggplot(tibble(x= x, y=y)) +
    geom_point(aes(x = x, y = y)) +
    xlab(xlabel) + ylab(ylabel)
    if(trend){
        tl = predict(lm(y~x))
        plt = plt + geom_line(aes(x = x, y = tl)) +
        xlab(xlabel) + ylab(ylabel)
    }
    return(plt)
}

get.interval = function(boot.tib, remove.na = FALSE){
    conf.int = tibble(PathwayID = rownames(boot.tib), conf.up = -99, conf.down = -99) %>%
                    column_to_rownames('PathwayID')
    for(name in rownames(boot.tib)){
        conf.tmp = boot.tib[name,] %>%
                    unlist(use.names = FALSE) %>%
                    sort
        if(remove.na){
            conf.tmp = conf.tmp[!is.na(conf.tmp)]
            if(length(conf.tmp) < 0.1 * dim(boot.tib)[2]){
                warning(paste0('Gene ', name, ': more than 90% of values are NAs'))
            } else if(length(conf.tmp) < 0.5 * dim(boot.tib)[2]){
                warning(paste0('Gene ', name, ': more than 50% of values are NAs'))
            }
        }
        q025 = (length(conf.tmp) * 0.025) %>% round
        q975 = (length(conf.tmp) * 0.975) %>% round
        if((q025 == 0) | (q975 == q025) | q975 < 10){
            conf.int[name,'conf.down'] = NA
            conf.int[name,'conf.up'] = NA
        } else {
            conf.int[name,'conf.down'] = conf.tmp[q025]
            conf.int[name,'conf.up'] = conf.tmp[q975]
        }
    }
    return(conf.int)
}


getKeggGenes = function(pwset){
    niter = length(pwset)%/%10
    kegg.set = list()
    if (length(pwset)%%10 != 0) {niter = niter +1}
    for(iter in 1:niter){
        range = (((iter - 1) * 10)+1):(iter * 10)
        kegg.tmp = keggGet(pwset[range] %>% na.omit)
        tmplist = kegg.tmp
        names(tmplist) =  sapply(kegg.tmp, function(x){x$ENTRY[[1]]})
        kegg.set = kegg.set  %>% append(tmplist)
    }
    return(kegg.set)
}


parse.name = function(x){
    out = strsplit(x,'\t') %>% unlist %>% `[`(1:2)
    out[2] = out[2] %>% strsplit(' - Homo') %>%
    unlist %>% `[`(1)
    return(out)
}

annotatePW = function(filenames){
    pw.list = list()
    for(name in filenames){
        f = file(name)
        sub.list = base::readLines(f) %>% lapply(parse.name)
        close(f)
        pw.list = append(pw.list, sub.list)
    }
    pw.tib = tibble(PathwayID = pw.list %>% sapply(function(x){x[1]}),
            Description = pw.list %>% sapply(function(x){x[2]}))
    return(pw.tib)
}

makeCommonAnno = function(anno.tib, meta, col.in.anno, col.in.meta){
    listed <- lapply(anno.tib[col.in.anno], function(x) {
        sapply(strsplit(x, ","),str_trim )
    })[[1]]
    names(listed) = anno.tib %>% pull("Label")
    listed = listed %>% unlist
    names(listed) = gsub("\\d$", "", names(listed))
    listed = setNames(names(listed), listed)
    meta['commonLabel'] = sapply(meta[col.in.meta], function(x){listed[x]}) %>% unname
    return(meta)
}
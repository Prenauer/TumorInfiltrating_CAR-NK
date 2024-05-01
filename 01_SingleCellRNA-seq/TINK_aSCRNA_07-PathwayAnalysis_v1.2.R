#!/usr/bin/env Rscript
library(dplyr)
library(stringr)
library(Seurat)
library(SeuratWrappers)
library(patchwork)
library(ggplot2)
library(reticulate)
library(future)
require(igraph)
require(network)
require(sna)
library(org.Mm.eg.db)
require(ggrastr)
require(ggrepel)


setwd('scRNAseq')
options(future.globals.maxSize= 10000*1024^2)
plan('multicore')


## Load data
so <- readRDS(file = 'Data/so_NK_v1.0.RDS')


## Get transcriptional regulators from Between-NK comparison
tr <- AnnotationDbi::select(org.Mm.eg.db, keytype="GOALL", keys="GO:0140110", columns="SYMBOL")
tr <- unique(tr$SYMBOL)
de <- read.delim('Data/06c_de_compare_NK_subtypes_v1.0.txt')
de <- de[with(de, order(celltype, -lfc)),]
topg <- de[which((de$adj_pval < 0.01) & ((de$lfc) > 0.5) & (de$max.detection.rate > 0.2) & (de$name %in% tr)),]
topg <- slice_head(topg, by=celltype, n=10)[,'name'] %>% unique()
d <- data.frame(Celltype=so$Celltype, t(so@assays$RNA@scale.data[topg,]))
d <- dplyr::reframe(d, .by=c('Celltype'), across(where(is.numeric), mean))
d <- data.frame(d[, 2:ncol(d)], row.names=d[,1])
pheatmap::pheatmap(d, cluster.method='ward.D2', cutree_cols=6, cutree_rows=2, cluster_rows=F,
                   treeheight_row=15,treeheight_col=15, angle_col='90', fontsize_col=7,fontsize_row=8,
                   filename='Figures/07a_hm_top5_TF_byCelltype_NK_comparison_v1.1.pdf', height=2.75, width=7)
topg <- c('Rora','Batf3','Tcf7','Klf2','Zeb2','Nr4a1','Rel','C1qbp','Hmgb2','Gmnn','Ezh2','Ifi204','Trim30a','Irf7','Muc1','Nfat5','Irf8', 'Eomes','Zbtb20')
p <- VlnPlot(so, pt.size=0, group.by='Celltype2', fill.by='feature', topg, stack=T, add.noise=F, adjust=2) + theme(legend.position = 'none')
ggsave(plot=p, 'Figures/07b_vln_top_TF_byCelltype_NK_comparison_v1.0.pdf', height=4, width=6)


## Function for pathway analysis
RunPW <- function(de, order.by='f_statistic', thresh.lfc = 1, thresh.q = 0.01, organism = 'mmusculus'){
    require(gprofiler2)
    if(length(order.by) == nrow(de)) de$statistic <- order.by
    if(length(order.by) == 1) de$statistic <- de[,order.by]
    
    glist <- de[with(de, order(statistic, decreasing = T)),]
    glist <- list(up=glist$name[which((glist$lfc > thresh.lfc) & (glist$adj_pval < thresh.q))] %>% head(n=500),
                  dn=glist$name[which((glist$lfc < -thresh.lfc) & (glist$adj_pval < thresh.q))] %>% head(n=500))
    
    gost <- do.call(rbind, lapply(c('up','dn'), function(direction){
        res <- gost(query = glist[[direction]], organism = organism, evcodes = T, ordered_query = T,significant = F,user_threshold = 1,correction_method = 'gSCS',domain_scope = 'known', sources = 'GO:BP')$result
        res$percent <-  res$intersection_size/res$term_size
        res <- res[,c('term_id','term_name','p_value','term_size','query_size','intersection_size','percent', 'precision','recall','intersection')]
        res$term_name <- paste0(toupper(substr(res$term_name,1,1)), substr(res$term_name,2,nchar(res$term_name)))
        if(nrow(res) > 0) return(cbind(direction=direction, res))
        if(nrow(res) == 0) return(NULL)
    }))
    return(gost)
}
Pathway_Network_Analysis <- function(gost, leiden_res = 1.2, filename.results, n_metapw=5){
    ## Calculate similarity matrix
    SimilarityMatrix <- function(res, overlap.pct = 0.5) # must have intersection and term_id columns
    { 
        require(stringr)
        # iterate rows, convert gene string to lists of character-vectors
        glists <- lapply(res$intersection, function(x) x %>% str_split(., ',') %>% unlist())
        gInter <- lapply(glists, function(g1) lapply(glists,function(g2) intersect(g1,g2) %>% length()) %>% unlist())
        gUnion <- lapply(glists, function(g1) lapply(glists,function(g2) union(g1,g2) %>% length()) %>% unlist())
        gMin <- lapply(glists, function(g1) lapply(glists,function(g2) min(length(g1),length(g2))) %>% unlist())
        m <- sapply(1:length(glists), function(i) {
            jaccard <- (gInter[[i]]/gUnion[[i]]) * (1 - overlap.pct)
            overlap <- (gInter[[i]]/gMin[[i]]) * overlap.pct
            return(jaccard + overlap)
        }) %>% as.matrix()
        dimnames(m) <- list(res$term_id,res$term_id)
        return(m)
    }
    ## Make edge matrix
    EdgeMatrix <- function(res, simThresh = 0.375, parallels = F, overlap.pct = 0.5){
        m <- SimilarityMatrix(res, overlap.pct)
        ## filter by similarity score
        edges <- reshape2::melt(m, as.is=T)
        edges <- edges[(edges$Var1 != edges$Var2) & (edges$value >= simThresh),]
        if(!parallels){
            parallels <- apply(edges,1,function(x) sort(as.character(x[1:2])) %>% paste0(.,collapse = ','))
            edges <- edges[!duplicated(parallels),]
        }
        names(edges) <- c('from','to','similarity_coefficient')
        edges$name <- edges$from
        return(edges)
    }
    
    # Run directions separately
    res <- lapply(c('up','dn'), function(direction){
        gost.subset <- gost[which(gost$direction == direction),]
        if(nrow(gost.subset) <= 2) return(NULL)
        if(nrow(gost.subset) > 2){
            ## make iGraph
            e <- EdgeMatrix(gost.subset)
            v <- gost.subset[,c('term_id','term_name','term_size','percent','intersection','p_value','precision','recall')]
            names(v) <- c('name','description','go_term_size','go_percent','genes','pval','precision','recall')
            v <- v[which(v$name %in% union(e$from,e$to)),]
            net <- igraph::graph_from_data_frame(d = e,directed = F, vertices = v)
            
            ## cluster the nodes
            net <- set_edge_attr(net, 'weight', value =  edge_attr(net, 'similarity_coefficient'))
            #leiden_res <- modularity(net, membership(cluster_leading_eigen(net)))
            m <- rbind(e[,1:3], data.frame(from=unique(c(e$from,e$to)), to=unique(c(e$from,e$to)), similarity_coefficient=1)) %>% unique()
            m <- reshape2::dcast(m, to~from, value.var='similarity_coefficient', fill=0)
            m <- data.frame(m[match(v$name,m[,1]),v$name], row.names=v$name)
            opt.resolution <- do.call(rbind, lapply(seq(0.3,1.5,0.3), function(r){
                l <- cluster_leiden(graph = net, objective_function = 'modularity', resolution_parameter = r, n_iterations = 500)
                #CalcWSS <- function(x,y) (length(x)-1)*(var(x)+var(y))
                #wss <- reframe(cbind(v[,c('V1','V2')],cluster=l$membership), .by='cluster', wss=CalcWSS(V1,V2))$wss %>% sum()
                #asw <- mean(cluster::silhouette(as.integer(l$membership), m)[,3])
                return(data.frame(resolution=r, nclust=l$nb_clusters, quality=l$quality))
            }))
            #leiden_res <- opt.resolution$resolution[which(opt.resolution$quality == max(opt.resolution$quality))]
            cat(leiden_res,'\n')
            l <- cluster_leiden(graph = net, objective_function = 'modularity', resolution_parameter = leiden_res, n_iterations = 500)
            v$cluster <- membership(l)[v$name]
            
            
            ## choosing meta-pathways
            v$intersection_size <- v$go_term_size * v$go_percent
            #metaPW <- v[with(v, order(precision, decreasing = T)),] %>% .[(.[,'go_term_size'] > 50),] %>% .[(.[,'intersection_size'] > quantile(v$intersection_size, 0.1)),] %>% .[!duplicated(.[,'cluster']),] %>% .[with(., order(cluster, decreasing = F)),]
            metaPW <- v[with(v, order(precision, decreasing = T)),] %>% .[(.[,'go_term_size'] > 50),] %>% .[(.[,'intersection_size'] > quantile(v$intersection_size, 0.1)),] %>% .[!duplicated(.[,'cluster']),] %>% .[with(., order(cluster, decreasing = F)),]
            metaPW <- structure(metaPW$description, names = metaPW$cluster)
            v$metaPW <- metaPW[as.character(v$cluster)]
            metaPW <- head(metaPW, n_metapw)
            v$label <- ''
            v[which(v$description %in% metaPW),'label'] <- v[which(v$description %in% metaPW),'description']
            v$logFDR <- -log10(v$pval)
            
            ## prepare graphing attributes
            l <- layout_with_fr(net, minx = rep(-100,nrow(v)), maxx = rep(100,nrow(v)), miny = rep(-100,nrow(v)), maxy = rep(100,nrow(v))) 
            if((max(l[,1]) - min(l[,1])) > (max(l[,2]) - min(l[,2]))) l <- l[,c(2,1)]
            v <- cbind(v,as.data.frame(l))
            v$color <- ifelse(v$cluster > n_metapw, NA, v$cluster) %>% factor()
            v$label <- sapply(1:nrow(v), function(x) ifelse(nchar(v$label[x]) > 0, paste0(toupper(substr(v$label[x],1,1)),substr(v$label[x],2,nchar(v$label[x]))), ''))
            v$label <- sapply(1:nrow(v), function(x) ifelse(nchar(v$label[x]) > 0, paste0(c(unlist(str_split(strwrap(v$label[x], width=25), "\n")), paste0('(adj.p = ',format(v$pval[x],scientific = T,digits = 3),')')),collapse = '\n'),''))
            e <- igraph::as_data_frame(net, 'edges')
            e.names <- colnames(e)
            e <- merge(e, v[,c('name','V1','V2')], by.x = 'from', by.y = 'name')
            e <- merge(e, v[,c('name','V1','V2')], by.x = 'to', by.y = 'name')
            colnames(e) <- c(e.names, 'x0','y0','x1','y1')
            e$size <- e$similarity_coefficient *1
            v$size <- (v$logFDR/max(v$logFDR)) * 14
            v$direction <- direction
            return(list(e = e, v = v))
        }
    })
    ## write results
    output <- do.call(rbind, lapply(res, function(x) if(!is.null(x)) return(x[['v']]) ))
    output$label <- stringr::str_replace_all(output$label, '\n', ' ')
    write.table(output, filename.results, sep='\t', quote=F, row.name=F)
    ## Plot results
    p <- lapply(1:2, function(i) {
        ggcolor <- function(n) hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[1:n]
        x <- res[[i]]
        if(!is.null(x)) {
            e <- x[['e']]
            v <- x[['v']]
            v <- v[which(v$cluster <= n_metapw),]
            e <- e[which((e$from %in% v$name & e$to %in% v$name)),]
            v$color <- ggcolor(max(as.integer(v$color),na.rm = T))[v$color]
            tmp <- ggplot(e) +
                ggrastr::rasterize(geom_segment(aes(x = x0, xend = x1, y = y0, yend = y1),linewidth = e$size*0.5, alpha = 0.2, color = '#AD88CE')) +
                geom_point(data = v, aes(x = V1, y = V2, fill = color, size = size), colour = 'gray60',alpha = 0.6, shape = 21) +
                scale_fill_discrete(na.value = 'gray90',guide = 'none') +
                geom_label_repel(data = v, aes(x = V1, y = V2, label = label),point.padding = unit(1.6, 'lines'), size = 2, 
                                 alpha = .7, seed = 42, force = 25,force_pull = 1,  max.iter=1e7, #nudge_y = 0.5,
                                 min.segment.length = 0, color = 'black', max.overlaps = 500, lineheight = 0.8) + 
                theme_void() +
                theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
                      axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                      plot.background = element_rect(colour = "grey50")) +
                scale_alpha(guide = 'none') + 
                annotate('text', x=min(v$V1), y=max(v$V2), hjust=0, vjust=1, label=c('Upreg.','Downreg.')[i]) +
                scale_size_continuous(name = str_wrap('Adj. p (-log10)',8))
            return(tmp)
        }
        if(is.null(x)) return(patchwork::plot_spacer())
    })
    return(cowplot::plot_grid(plotlist=p, align='hv',ncol=2))
}




## Downstream pathway analyses for tumor-infiltration
# pathway analysis by tumor model
l <- list()
p <- list()
for(model in c('B16_tumor','E0771_tumor')){
    de <- read.delim(paste0('Data/06a_de_',model,'_infil_v1.1.txt'))
    gost <- RunPW(de, thresh.lfc=1, thresh.q=0.01)
    l[[model]] <- gost
    gost=l[[model]]
    gost <- gost[(gost$term_size < 750) & (gost$intersection_size > 2) & (gost$p_value < 0.05),]
    p[[model]] <- Pathway_Network_Analysis(gost, 1, paste0('Data/07a_pw_',model,'_infil_v1.1.txt'), n_metapw=8)
}
p <- cowplot::plot_grid(plotlist=p, align='hv',ncol=2)
ggsave(plot = p, paste0('Figures/07a_networkplot_tumor_infil_v1.1.pdf'), dpi = 300, height = 5, width = 12.5, units = 'in')
# pathway analysis by celltype and tumor model
p <- lapply(c('B16_tumor','E0771_tumor'), function(cancer){
    lapply(levels(so$Celltype), function(celltype){
        de <- read.delim(paste0('Data/06b_de_',cancer,'_',celltype,'_v1.1.txt'))
        gost <- RunPW(de, thresh.lfc=1, thresh.q=0.05)
        gost <- gost[(gost$term_size < 750) & (gost$intersection_size > 2) & (gost$p_value < 0.05),]
        p <- Pathway_Network_Analysis(gost, 1, paste0('Data/07b_pw_',model,'_',celltype,'_v1.1.txt'), n_metapw=8)
        return(p)
    })
})
p <- cowplot::plot_grid(plotlist=p, align='hv',ncol=2)
ggsave(plot = p, paste0('Figures/07a_networkplot_tumor_infil_byCelltype_v1.1.pdf'), dpi = 300, height =4*4, width = 10, units = 'in')


















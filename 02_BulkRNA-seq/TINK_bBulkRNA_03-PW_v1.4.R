library(edgeR)
library(stringr)
library(dplyr)
require(ggrepel)
library(ggplot2)
library(ggrastr)
library(patchwork)
library(reshape2)
options(stringsAsFactors = F)

setwd('prRNA')


## Get DE transcripts
enst <- lapply(c('unstim','stim'), function(gp){
    de <- read.delim(paste0('Data/02a_de_primaryNK_CALHM2KO_',gp,'_v1.1.txt'))
    return(de$enst[which(abs(de$lfc) > 1 & de$adj_pval < 0.01)])
}) %>% unlist() %>% unique()
write.table(enst, 'Data/enst.txt', sep='\t', quote=F, row.names=F)


## Make heatmap of DE transcripts
d <- lcpm[which(lcpm$ENST %in% enst),]
g <- d$GeneSymbol
d <- as.matrix(d[,4:ncol(d)])
rownames(d) <- g
pheatmap::pheatmap(d, cluster_cols=F, clustering_method='ward.D2', scale='row', cutree_rows=4,angle_col='90',
                   cutree_cols=2, treeheight_row=20, treeheight_col=15,color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 8, name ="RdYlBu")))(100),
                   filename=paste0('Figures/heatmaps/03d_hm_all_DEG.pdf'), height=8,width=3)


## Pathway analyses of clustered DE transcripts
cl <- hclust(dist(t(scale(t(d)))), method='ward.D2')
l <- cutree(cl,4)
l <- l[!is.na(names(l))]
res.bp <- do.call(rbind, lapply(1:4, function(i){
    g <- names(l)[l==i]
    res <- gost(query = g, organism = 'hsapiens', evcodes = T, ordered_query = F,significant = F,user_threshold = 1,correction_method = 'gSCS',domain_scope = 'known', sources = c('GO:BP'))$result
    res <- res[which((res$p_value < 0.05) & (res$term_size < 750) & (res$intersection_size > 2) & (res$intersection_size >= 4)),]
    if(nrow(res) > 0) {
        res$percent <-  res$intersection_size/res$term_size
        res <- res[,c('term_id','term_name','p_value','term_size','query_size','intersection_size','percent', 'precision','recall','intersection')]
        res$term_name <- paste0(toupper(substr(res$term_name,1,1)), substr(res$term_name,2,nchar(res$term_name)))
        res <- cbind(Cluster=i, res)
    }
    return(res)
}))
write.table(res.bp, 'Data/03a_pw_primaryNK_CALHM2KO_v1.2b.txt', sep='\t', quote=F, row.names=F)




## Barplots for top pathways
res <- read.delim('Data/03a_pw_primaryNK_CALHM2KO_v1.2b.txt')
# create empty data so that each cluster in plot has at least 5 bars length.
empty.data <- data.frame(Cluster=rep(1:4,5), term_id=paste0('empty',rep(1:4,5),1:5), term_name='', p_value=1, term_size=500, query_size=50,
                         intersection_size=10, percent=0.5, precision=0.01, recall=0.01, intersection='')
res <- rbind(res, empty.data)
res$Cluster <- factor(res$Cluster, levels=(c(3,1,4,2)))
res <- res[with(res, order(Cluster, p_value, -precision)),]
pw <- slice_head(res[res$intersection_size >= 4,], by='Cluster', n=5)
pw$term_id <- factor(pw$term_id, levels=rev(unique(pw$term_id)))
p <- ggplot(pw) + geom_bar(aes(x=term_id, y=-log10(p_value), fill=factor(Cluster)), stat='identity', color='gray30') +
    coord_flip() + scale_x_discrete(labels=structure(str_wrap(pw$term_name,20), names=as.character(pw$term_id))) +
    theme_light() + theme(axis.text.y = element_text(lineheight = 0.6), legend.position='none') + labs(y='q value (-log10)',x='')
ggsave(plot=p, 'Figures/03e_barplot_pathway_analysis.pdf', height=8,width=3)


# Make Volcano plots where dots are colored according to pathway
for(gp in c('unstim','stim')){
    res <- read.delim(paste0('Data/02a_de_primaryNK_CALHM2KO_',gp,'_v1.1.txt'))
    res$logq <- -log10(res$adj_pval)
    res.de <- res[(res$adj_pval < 0.01) & (abs(res$lfc) > 1),]
    signif.lines <- data.frame(x = c(min(res$lfc), -1, 1, 1), xend = c(-1, -1, 1, max(res$lfc)),
                               y = c(-log10(0.05), -log10(0.05), max(res$logq), -log10(0.05)), 
                               yend = c(-log10(0.05),  max(res$logq), -log10(0.05), -log10(0.05)))
    genes.label.sig <- res.de[order(res.de$statistic, decreasing = T),]
    genes.label.sig <- c(genes.label.sig[(genes.label.sig$lfc > 0),'enst'] %>% head(8),
                         genes.label.sig[(genes.label.sig$lfc < 0),'enst'] %>% head(8))
    genes.label.fc <- res.de[order(res.de$lfc, decreasing = T),]
    genes.label.fc <- c(genes.label.fc[which(!is.na(genes.label.fc$enst)),'enst'] %>% head(15),
                        genes.label.fc[which(!is.na(genes.label.fc$enst)),'enst'] %>% tail(15))
    genes2plot.small <- union(genes.label.sig, genes.label.fc)
    res.label.small <- res.de[which(res.de$enst %in% genes2plot.small),]
    res.de$pw <- NA
    #ggcolor <- hcl(h = seq(15, 375, length = nrow(pw) + 1), l = 65, c = 100)[1:nrow(pw)]
    for(i in nrow(pw):1){ # reverse the order so that genes will be labeled by the most significant pathway.
        pw.genes <- stringr::str_split(pw$intersection[i],',') %>% unlist()
        res.de$pw[which(res.de$name %in% pw.genes)] <- pw$term_name[i]
    }
    
    p2 <- ggplot(res, aes(x = lfc, y = logq)) + 
        ggrastr::geom_point_rast(alpha = .25, size = 2.5, color = 'gray10', shape = 16) +
        ggrastr::geom_point_rast(data = res.de, alpha = 1, size = 2.5, color = 'gray10', shape = 16) +
        ggrastr::geom_point_rast(data = res.de, alpha = 1, color = 'gray', shape = 16) +
        #ggrastr::geom_point_rast(data = res.de[which(res.de$lfc > 0),], alpha = 1, color = 'firebrick', shape = 16) +
        #ggrastr::geom_point_rast(data = res.de[which(res.de$lfc < 0),], alpha = 1, color = 'steelblue', shape = 16) +
        ggrastr::geom_point_rast(data = res.de[!is.na(res.de$pw),], color='gray10', size=4, alpha = 1, shape = 16) +
        ggrastr::geom_point_rast(data = res.de[!is.na(res.de$pw),], color='white', size=3, alpha = 1, shape = 16) +
        ggrastr::geom_point_rast(data = res.de[!is.na(res.de$pw),], aes(color=factor(pw)), size=3, alpha = 1, shape = 16) +
        geom_segment(data = signif.lines, aes(x = x, xend = xend, y = y, yend = yend), 
                     linetype = 'dashed', color = 'gray20') +
        ggrepel::geom_text_repel(data = res.label.small, aes(label = name), alpha = 0.8, size = 2.5, nudge_y = 2,force=5, segment.size=0.04,
                                 max.iter = 100000,min.segment.length = 0, fontface = 'italic', max.overlaps = 100) +
        labs(x = 'Log2 fold-change', y = 'q value (-log10)') +
        cowplot::theme_cowplot()
    p1 <- p2 + theme(legend.position='none')
    ggsave(plot = cowplot::plot_grid(p1, p2, align='hv', ncol=2), filename = paste0('Figures/03f_volcano_prNK_',gp,'_v1.2.pdf'), height = 3.5, width = 8)
}


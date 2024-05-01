#!/usr/bin/env Rscript
library(dplyr)
library(stringr)
library(Seurat)
library(SeuratWrappers)
library(patchwork)
library(ggplot2)
library(reticulate)
library(future)

setwd('scRNAseq')
options(future.globals.maxSize= 10000*1024^2)
plan('multicore')


## Load data
so <- readRDS(file = 'Data/so_NK_v1.0.RDS')


## Get lists of NK signature genes
glist <- list(NKT=c('Cd3e','Cd3d'),
              ILC1=c('Il7r','Rora','Gpr183','Cxcr6','Sdc4','Tmem176a','Tmem176b','Tcrg-C4'),
              cNK=c('Eomes','Klra4','Klra8'),
              mNK=c('Itgam','Zeb2','Ly6c2','Klrg1','Cma1'),
              iNK=c('Emb','Xcl1','Cd27'),
              tissueResidency=c('Tcf7','Emb','Kit','Ltb','Cxcr3','Cd160','Cd7','Klrb1f'),# c('Cd69','Cxcr3','Rgs1','Rgs2','Xcl1','Crtam') # iNK4, ILC1
              effectorResponse=c('Gzma','Gzmb','Prf1','Ifng','Lgals1'),
              regulatoryResponse=c('Emb','Xcl1','Tgfb1','Ccl3','Ccl4'),
              survival=c('Bax','Bcl2','Bcl2l11','Bcl10','Cycs'),
              proliferation=c('Mki67','Top2a'))
glist.simplified <- list(NKT=c('Cd3e'),
                         ILC1=c('Rora','Gpr183','Cxcr6','Sdc4','Tmem176a','Tmem176b','Tcrg-C4'),
                         mNK=c('Itgam','Klrg1','Cma1','Zeb2'),
                         iNK=c('Emb','Xcl1','Cd27'),
                         tissueResidency=c('Tcf7','Kit','Cxcr3','Cd160'),# c('Cd69','Cxcr3','Rgs1','Rgs2','Xcl1','Crtam') # iNK4, ILC1
                         effectorResponse=c('Gzma','Gzmb','Gzmc', 'Prf1','Ifng','Lgals1'),
                         regulatoryResponse=c('Emb','Xcl1','Tgfb1','Ccl3','Ccl4'),
                         survival=c('Bax','Bcl2','Bcl2l11','Bcl10','Cycs'),
                         proliferation=c('Mki67','Top2a'))


## Function for dotplots
DotPlotFunction <- function(features, group.by, dot.color='blue', cells.to.include=NULL, cluster.by='mean', scale_size=c(0,4), cluster.row=F,cluster.col=F, cluster.method='ward.D2'){
    if(is.null(cells.to.include)) cells.to.include <- colnames(so@assays$RNA@data) 
    data <- data.frame(t(as.matrix((so@assays$RNA@data[which(rownames(so@assays$RNA@data) %in% features),]))), 
                       group=group.by, color=dot.color)[cells.to.include,]
    data <- reshape2::melt(data, id.vars=c('group','color')) %>% data.frame()
    colnames(data)[3:4] <- c('gene','expr')
    data <- reframe(data, .by=c('gene','group','color'), mean=mean(expr), n_expr_cells=sum(expr > 0), n_total_cells=length(expr), exp_rate=100*sum(expr > 0)/length(expr))
    data <- data[with(data, order(group, color, gene)),]
    features <- stringr::str_replace(features, '-', '\\.')
    # cluster data 
    tmp <- reshape2::dcast(data=data, formula = group ~ gene, value.var=cluster.by, drop=T) 
    tmp <- data.frame(tmp[,-1], row.names=tmp[,1])
    hc.col <- hclust(dist(tmp), method=cluster.method)
    hc.row <- hclust(dist(t(tmp)), method=cluster.method) 
    data$group <- factor(data$group, levels=rev(unique(data$group)))
    data$gene <- factor(data$gene, levels=features)
    if(cluster.row) data$group <- factor(data$group, levels=hc.col$labels[hc.col$order])
    if(cluster.col) data$gene <- factor(data$gene, levels=hc.row$labels[hc.row$order])
    # windsorize
    data$mean[which(data$mean > quantile(data$mean, 0.9))] <- quantile(data$mean, 0.9)
    # dot plot
        p <- ggplot(data, aes(x=gene, y=group)) + 
            geom_point(aes(size=exp_rate), alpha=0.8, color='gray') +
            geom_point(aes(size=exp_rate, alpha=mean), color=dot.color) + theme_test() + 
            scale_size(range=scale_size) + scale_alpha(range=c(0,1))  +
            labs(x='Celltypes',y='', title='') + 
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    if(cluster.row) p <- p + ggh4x::scale_y_dendrogram(hclust=hc.col)
    if(cluster.col) p <- p + ggh4x::scale_x_dendrogram(hclust=hc.row)
    return(p)
}


## Dot plots of NK signature genes
so$TimePretty <- paste0(formatC(so$Time, width=2, flag='0'),'dpi')
p <- DotPlotFunction(features=unlist(glist.simplified[c('NKT','ILC1','mNK','iNK','tissueResidency','proliferation')]), 
    group.by=paste0(so$Celltype,'_',so$TimePretty))
ggsave(plot=p, width = 6, height = 6,filename = 'Figures/05a_dot_NK_signatures_byTime_v1.0.pdf')
p <- DotPlotFunction(features=unlist(glist.simplified[c('NKT','ILC1','mNK','iNK','tissueResidency','proliferation')]), 
    group.by=paste0(so$Celltype,'_',so$Source))
ggsave(plot=p, width = 6, height = 6,filename = 'Figures/05b_dot_NK_signatures_byTissue_v1.0.pdf')
p <- DotPlotFunction(features=unlist(glist.simplified[c('effectorResponse','regulatoryResponse')]), 
    group.by=paste0(so$Celltype,'_',so$TimePretty))
ggsave(plot=p, width = 3.75, height = 6,filename = 'Figures/05c_dot_NK_response_signatures_byTime_v1.0.pdf')
p <- DotPlotFunction(features=unlist(glist.simplified[c('effectorResponse','regulatoryResponse')]), 
    group.by=paste0(so$Celltype,'_',so$Source))
ggsave(plot=p, width = 3.75, height = 6,filename = 'Figures/05d_dot_NK_response_signatures_byTissue_v1.0.pdf')


## Featureplots of NK signatures
p <- FeaturePlot(so, paste0(names(glist),'_score'), order=T, min.cutoff='q01',max.cutoff='q90', raster=T)
ggsave(plot=p, width = 12, height = 7.5,filename = 'Figures/05e_fp_NK_signature_v1.0.pdf')


## Split-view umap plot (colored by celltype, split by dataset)
so$Sample_Type <- paste0(so$Cancer, '_', so$Source, '_', so$TimePretty) %>% stringr::str_replace(., 'Control_Culture', 'Culture')
so$Sample_Type <- factor(so$Sample_Type, levels=c('Culture_00dpi','B16_Spleen_07dpi', 'B16_Spleen_15dpi','B16_Tumor_07dpi', 'B16_Tumor_15dpi','E0771_Spleen_07dpi', 'E0771_Spleen_15dpi','E0771_Tumor_07dpi', 'E0771_Tumor_15dpi','Merged_Datasets'))
data <- data.frame(so@reductions$umap@cell.embeddings, Celltype=so$Celltype, Sample_Type=so$Sample_Type)
data <- rbind(data, data.frame(data[,1:3], Sample_Type='Merged_Datasets'))
p.celltype <- ggplot(data, aes(x=UMAP_1, y=UMAP_2, color=Celltype)) + ggrastr::geom_point_rast(scale=0.25, alpha=0.5, raster.dpi = 300) + facet_wrap(vars(Sample_Type)) + theme_test()
ggsave(plot = p.celltype, file = 'Figures/05f_split_umap_NK_celltypes_v1.0.pdf', unit = 'in', width = 7, height = 6)
celltype.proportions <- dplyr::mutate(so@meta.data, .by=c('Sample'), n_cells_in_dataset = length(Sample))
celltype.proportions <- dplyr::reframe(celltype.proportions, .by=c('Sample','Time','Source','Cancer','Celltype'), n_cells_in_dataset=unique(n_cells_in_dataset), n_cells=length(Sample), pct_cells=100*length(Sample)/n_cells_in_dataset) #%>% unique()
write.table(celltype.proportions, 'Data/05a_NK_celltype_proportions_by_dataset_v1.0.txt',sep = '\t', quote=F, row.names=F)


## Make dotplots for NK function-related markers
glist <- list(ActivatingReceptors = c('Cd28','Fcgr3','Klrc2','Klrk1','Ncr1'),
              InhibitoryReceptors = c('Cd244a','Klrg1','Klrc1','Klrb1c','Ly6c2'),
              ChemokineReceptors = c('Ccr2','Ccr5','Cxcr3','Cxcr4','Cxcr6','S1pr5'),
              AdhesionReceptors = c('Cd2','Spn','Cd44','Sell','Itga2','Itgal'))
DefaultAssay(so) <- 'RNA'
plot_grouping <- paste0(so$Celltype,'_',so$TimePretty,'_',so$Source) 
colorpal <- structure(c('tomato3', 'cyan3', 'orange2','seagreen'), names=names(glist))
p <- lapply(names(glist), function(g){
    DotPlot3(cells.to.include=which(!(so$Celltype %in% c('ILC1','NKT'))), features=glist[[g]], 
             group.by=plot_grouping, dot.color=colorpal[g], cluster.row=T,
             cluster.col=F,cluster.by='exp_rate') + labs(title=g)
})
p <- cowplot::plot_grid(plotlist=p, ncol=4, align='hv')
ggsave(plot=p, 'Figures/05g_dot_NK_functional_characterization_v1.0.pdf', dpi=300, height=7, width=24)



## compare time, tumor-model, and tissue across celltypes
so$TimePretty <- paste0(formatC(so$Time, width=2, flag='0'),'dpi')
o <- list(TimePretty = c('00dpi','07dpi','15dpi'), Source = c('Culture','Tumor','Spleen'), Cancer = c('Control','B16','E0771'))
meta <- data.frame(so@meta.data[,c('Cancer','TimePretty','Celltype','Source')], dataset=apply(so@meta.data[,c('Cancer','TimePretty','Source')], 1,paste0,collapse='_'))
meta <- lapply(c('Cancer','TimePretty','Source'), function(x) {
    tmp <- merge(meta, reframe(meta, .by=x,total_cells=length(dataset)), by=x)
    tmp <- reframe(tmp, .by=c('Celltype', x), ct=length(Celltype), pct=100*(length(Celltype)/total_cells)) %>% unique()
    colnames(tmp)[2] <- 'Comparison'
    tmp$Comparison <- factor(tmp$Comparison, levels=o[[x]])
    return(tmp)
})
write.table(do.call(rbind, meta), file = 'Data/05b_celltype_comparison_percents_v1.0.txt', sep = '\t', quote=F, row.names = F)


## Stats of celltype proportions across different experimental conditions (pairwise Fisher exact with FDR correction) 
chisq.multcomp <- function(d, p.method="fdr") {
        x <- as.matrix(d)
        fun.p <- function(i,j) suppressWarnings(chisq.test(c(x[i], x[j])))$p.value
        tab.p <- pairwise.table(fun.p,as.character(x),p.adjust.method=p.method)
        call <- match.call()
        dname.x <- if(length(call$x)==1) {call$x} else {paste(call$x[1],"(",paste(call$x[-1],collapse=","),")",sep="")}
        r <- tab.p
        columns <- sapply(colnames(d), function(x) { rep(x, nrow(d))}) %>% as.character()
        rows <- rep(rownames(d), ncol(d))
        rownames(r) <- paste0(rows[-1], '_', columns[-1])
        colnames(r) <- paste0(rows[-length(rows)], '_', columns[-length(rows)])
        r <- reshape2::melt(as.matrix(r), na.rm = T)
        r <- r[which((stringr::str_split(r$Var1, '_', simplify = T)[,2]) == (stringr::str_split(r$Var2, '_', simplify = T)[,2])),]
        colnames(r) <- c('Group1','Group2','q_value')
        return(r)
}
d <- lapply(1:3, function(i) {
    # run stats
    tmp <- data.table::dcast(data.table::data.table(meta[[i]]), Comparison ~ Celltype, value.var='pct')
    tmp <- data.frame(tmp[,-1], row.names=unlist(tmp[,1])) 
    r <- chisq.multcomp(tmp)
    return(data.frame(Comparison = c('Cancer','TimePretty','Source')[i], r))
})
write.table(do.call(rbind, d), file = 'Data/05c_celltype_comparison_stats_v1.0.txt', sep = '\t', quote=F, row.names = F)


## Bar plots of celltype proportions across different experimental conditions
p <- lapply(1:3, function(i) {
    ggplot(meta[[i]], aes(fill=Comparison, y= pct, x= Celltype)) + 
        labs(x='', y='Percent', title=c('Cancer','TimePretty','Source')[i]) +
        geom_bar(stat ='identity', position = 'dodge', color = 'gray30') + 
        theme_classic() +  theme(axis.text.x = element_text(angle=90, hjust=1)) %>% return()
})
p <- cowplot::plot_grid(plotlist=p, align='hv', nrow=3)
ggsave(plot = p, filename = file.path(DIR_FIG,paste0('Figures/05h_barplot_NK_celltype_pop_summary_v1.0.pdf')), dpi = 300, width = 5.5, height = 4.5, units = 'in')



## scatter plots comparing tumor models across each immune cell
d <- lapply(levels(so$Celltype2), function(celltype){
    cells <- which((so$Celltype2 == celltype) & (so$Time == 7))
    d <- log2(1 + AverageExpression(so[,cells], group.by='Cancer', slot='data')[['RNA']])
    return(data.frame(x=d[,1], y=d[,2], celltype=celltype))
}) %>% data.table::rbindlist()
d <- mutate(d, .by='celltype', Correlation=pcaPP::cor.fk(x=x,y=y))
d$Correlation <- paste0('tau = ',format(d$Correlation, digits=3))
d$Correlation[duplicated(d$celltype)] <- NA
p <- ggplot(d, aes(x=x, y=y)) + ggrastr::geom_point_rast(shape=20, scale=0.6, alpha=0.6, aes(color=celltype)) + 
    geom_smooth(method='lm') + 
    theme_test() + theme(legend.position='none') + 
    labs(x='B16F10 expression (mean log2-norm)', y='E0771 expression (mean log2-norm)') +
    facet_wrap(vars(celltype), ncol=4) 
ggsave(plot=p, 'Figures/05i_xy_cancer_correlation_07dpi_v1.0.pdf', height=4, width=5.5)
d <- lapply(levels(so$Celltype2), function(celltype){
    cells <- which((so$Celltype2 == celltype) & (so$Time == 15))
    d <- log2(1 + AverageExpression(so[,cells], group.by='Cancer', slot='data')[['RNA']])
    return(data.frame(x=d[,1], y=d[,2], celltype=celltype))
}) %>% data.table::rbindlist()
d <- mutate(d, .by='celltype', Correlation=pcaPP::cor.fk(x=x,y=y))
d$Correlation <- paste0('tau = ',format(d$Correlation, digits=3))
d$Correlation[duplicated(d$celltype)] <- NA
p <- ggplot(d, aes(x=x, y=y)) + ggrastr::geom_point_rast(shape=20, scale=0.6, alpha=0.6, aes(color=celltype)) + 
    geom_smooth(method='lm') + 
    theme_test() + theme(legend.position='none') + 
    labs(x='B16F10 expression (mean log2-norm)', y='E0771 expression (mean log2-norm)') +
    facet_wrap(vars(celltype), ncol=4) 
ggsave(plot=p, 'Figures/05j_xy_cancer_correlation_15dpi_v1.0.pdf', height=4, width=5.5)


## Calhm2 heatmap, comparing tumor-model, tissue, and celltype
DefaultAssay(so) <- 'RNA'
so$TissueAndModel <- paste0(so$Cancer, '_', so$Source) %>% stringr::str_replace(., 'Control_Culture', 'Culture')
so$TissueAndModel <- factor(so$TissueAndModel, levels=rev(c('Culture','B16_Spleen','B16_Tumor','E0771_Spleen','E0771_Tumor')))
data <- data.frame(calhm2 = so@assays$RNA@data['Calhm2', ], tissue=so$TissueAndModel, celltype=so$Celltype)
data <- data[which(so$Time != 15),]
data <- reframe(data, .by=c('tissue','celltype'), mean=mean(calhm2), n_expr_cells=sum(calhm2 > 0), 
                n_total_cells=length(calhm2), exp_rate=sum(calhm2 > 0)/length(calhm2))
data <- data[with(data, order(celltype, tissue)),]
write.table(data, 'Figures/RevFigures/hm_calhm2_detection.txt', sep='\t', quote=F, row.names=F)
data <- reshape2::dcast(data=data, formula = celltype ~ tissue, value.var='exp_rate')
data <- data.frame(data[,-1], row.names=data[,1])
anno_col <- data.frame(Tissue=c('Tumor','Tumor','Culture','Spleen','Spleen'), row.names=c('E0771_Tumor','B16_Tumor','Culture','E0771_Spleen','B16_Spleen'))
anno_color <- list(Tissue=c(Tumor='tomato2',Culture='darkcyan',Spleen='orange2'))
labels_col <- c('E0771','E0771','B16F10','B16F10','Culture')
pheatmap::pheatmap(data, clustering_method='ward.D2',display_numbers=T, angle_col='90', cutree_cols=2,labels_col=labels_col,
                   treeheight_row=15, treeheight_col=15, annotation_col=anno_col, annotation_colors=anno_color, fontsize_number=6,
                   filename='Figures/05k_hm_calhm2_detection.pdf', height=3, width=4.25)


## Calhm2 FeaturePlot
p <- FeaturePlot(so, order=T, 'Calhm2', split.by='Source', raster=T, max.cutoff='q75', pt.size=1)
ggsave(plot=p,'Figures/05l_fp_calhm2_detection.pdf', height=4, width=7.5)


# Calhm2 Dotplot
DefaultAssay(so) <- 'RNA'
so$TissueAndModel <- paste0(so$Cancer, '_', so$Source) %>% stringr::str_replace(., 'Control_Culture', 'Culture')
so$TissueAndModel <- factor(so$TissueAndModel, levels=rev(c('Culture','B16_Spleen','B16_Tumor','E0771_Spleen','E0771_Tumor')))
so$TissueModelTime <- paste0(so$Cancer, '_', so$Source, '_', so$Time, 'dpi') %>% stringr::str_replace(., 'Control_Culture', 'Culture')
data <- data.frame(calhm2 = so@assays$RNA@data['Calhm2', ], tissue=so$TissueModelTime, celltype=so$Celltype, category=so$Source)
data <- data[which(so$Time != 15),]
data <- reframe(data, .by=c('tissue','celltype','category'), mean=mean(calhm2), n_expr_cells=sum(calhm2 > 0), n_total_cells=length(calhm2), exp_rate=100*sum(calhm2 > 0)/length(calhm2))
data <- data[with(data, order(celltype, tissue)),]
write.table(data[,-3], 'Data/05d_hm_calhm2_detection.txt', sep='\t', quote=F, row.names=F)
# cluster data by expression rate
tmp <- reshape2::dcast(data=data, formula = celltype ~ tissue, value.var='exp_rate')
rownames(tmp) <- tmp$celltype
tmp <- tmp[,-1]
hc.col <- hclust(dist(tmp), method='ward.D2')
hc.row <- hclust(dist(t(tmp)), method='ward.D2') 
data$celltype <- factor(data$celltype, levels=hc.col$labels[hc.col$order])
data$tissue <- factor(data$tissue, levels=hc.row$labels[hc.row$order])
# dot plot
p <- ggplot(data, aes(x=celltype, y=tissue)) + 
    geom_point(aes(size=exp_rate), alpha=0.8, color='gray') + theme_test() +
    geom_point(aes(size=exp_rate, color=category, alpha=mean)) + theme_test() + 
    scale_size(range=c(0,10)) + scale_alpha(range=c(0,1)) +
    ggh4x::scale_x_dendrogram(hclust=hc.col) + ggh4x::scale_y_dendrogram(hclust=hc.row) +
    labs(x='Celltypes','') + scale_color_manual(values=c(Spleen='tomato2', Tumor='orange2', Culture='steelblue3')) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(plot=p, width = 8, height = 3,filename = paste0('Figures/05m_dot_calhm2_detection.pdf'))



## Decreased mature markers in tumor mNK cells
gene <- glist$mNK
DefaultAssay(so) <- 'RNA'
so$TissueAndModel <- paste0(so$Cancer, '_', so$Source) %>% stringr::str_replace(., 'Control_Culture', 'Culture')
so$TissueAndModel <- factor(so$TissueAndModel, levels=rev(c('Culture','B16_Spleen','B16_Tumor','E0771_Spleen','E0771_Tumor')))
col <- structure(c('tomato2','cyan4', 'orange2','steelblue3', 'green4'), names= c('Culture','B16_Spleen','B16_Tumor','E0771_Spleen','E0771_Tumor'))
p <- VlnPlot(so[, which((so$Celltype %in% c('mNK')))], adjust=2.5,
             pt.size=0, split.by='TissueAndModel', group.by='Celltype', gene, stack=T,add.noise=F) + scale_fill_manual(values=col) +
    geom_vline(xintercept=c(1:5), linetype='longdash', size=0.1, alpha=0.25)
ggsave(plot=p, width = 5, height = 3,filename = paste0('Figures/05n_vln_mNK_markers.pdf'))
# Export data (mean expression and expression rate)
d <- data.frame(t(as.matrix(so@assays$RNA@data[gene , which(so$Celltype %in% c('mNK'))])),
                Group=so@meta.data[which(so$Celltype %in% c('mNK')), 'TissueAndModel'])
ExpressionRate <- function(x) return(sum(x > 0)/length(x))
d <- reframe(d, .by='Group', Itgam.mean=mean(Itgam), Itgam.ExpressionRate=ExpressionRate(Itgam),
             Zeb2.mean=mean(Zeb2), Zeb2.ExpressionRate=ExpressionRate(Zeb2),
             Ly6c2.mean=mean(Ly6c2), Ly6c2.ExpressionRate=ExpressionRate(Ly6c2),
             Klrg1.mean=mean(Klrg1), Klrg1.ExpressionRate=ExpressionRate(Klrg1),
             Cma1.mean=mean(Cma1), Cma1.ExpressionRate=ExpressionRate(Cma1))
d <- d[c(3,1,2,4,5),]
write.table(d, 'Data/05e_mNK_marker_expression_stats.txt', sep='\t', quote=F, row.names=F)








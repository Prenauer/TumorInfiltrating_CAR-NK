#!/usr/bin/env Rscript
library(dplyr)
library(stringr)
library(Seurat)
library(SeuratWrappers)
library(patchwork)
library(ggplot2)
library(ggrastr)
library(ggrepel)
library(reticulate)
library(future)

setwd('scRNAseq')
options(future.globals.maxSize= 10000*1024^2)
plan('multicore')


## Load data
so <- readRDS(file = 'Data/so_NK_v1.0.RDS')


## Function to filter cells and genes
FilterData <- function(data, design, contrast){
    gdr <- 0.05
    cdr <- 0.10
    # remove Gm42418, as it overlaps the rRNA element Rn45s and represent rRNA contamination (PMID: 34795279).
    data <- data[which(rownames(data) != 'Gm42418'),]
    # get max gene expression rate of contrast and non-contrast groups
    cat('Calculate gene detection rate\n')
    groups <- factor(design[, which(colnames(design)==contrast)])
    gene.detection.rate <- do.call(cbind, lapply(levels(groups), function(group){
        cells <- which(groups == group)
        return(apply(data[,cells], 1, function(x) {sum(x > 0)})/length(cells))
    }))
    gene.detection.rate <- apply(gene.detection.rate, 1, max)
    # filter cells that don't express enough of the filtered genes
    cell.detection.rate <- apply(data[which(gene.detection.rate > gdr),], 2, function(x) {sum(x > 0)})/length(which(gene.detection.rate > gdr))
    data <- data[,which(cell.detection.rate > cdr)]
    design <- design[which(cell.detection.rate > cdr),]
    rownames(design) <- colnames(data)
    # recalculate gene expression rate in filtered cells
    groups <- factor(design[, which(colnames(design)==contrast)])
    gene.detection.rate <- do.call(cbind, lapply(levels(groups), function(group){
        cells <- which(groups == group)
        return(apply(data[,cells], 1, function(x) {sum(x > 0)})/length(cells))
    }))
    gene.detection.rate <- apply(gene.detection.rate, 1, max)
    # filter out lowly expressed genes
    data <- data[which(gene.detection.rate > gdr),]
    cat('Filtered data: ',nrow(data),' genes in ',ncol(data),' cells\n')
    return(list(data=data, design=design, gene.detection.rate=gene.detection.rate))
}


## Function to run DE
DE_TumorAnalysis <- function(data, design, contrast='Tumor'){
    filtered.data <- FilterData(data, design, contrast=contrast, cdr=cdr, gdr=gdr)
    fit <- glm_gp(filtered.data$data, design = filtered.data$design, size_factors = 'deconvolution', 
                  subsample = T, on_disk = F, verbose = T, ridge_penalty=0.0001)
    de <- test_de(fit, contrast = contrast, verbose = T) 
    de$max.detection.rate <- filtered.data$gene.detection.rate[de$name]
    return(de)
}


## Run DE for tumor-infiltration in each tumor model
Tumor <- as.integer(so$Source == 'Tumor')
InVivo <- as.integer(so$Source != 'Culture')
Cell.detection.rate <- so$Cell.detection.rate
# B16F10 tumor infiltration
cells <- which((so$Cancer %in% c('B16','Control')) & (so$Time != 15)) 
design <- model.matrix(~ 0 + InVivo + Cell.detection.rate + Tumor)
de <- DE_TumorAnalysis(data=so$RNA@counts[,cells], design=design[cells,])
write.table(de, 'Data/06a_de_B16_tumor_infil_v1.0.txt', sep='\t', quote=F, row.names=F)
# E0771 tumor infiltration
cells <- which((so$Cancer %in% c('E0771','Control')) & (so$Time != 15))
design <- model.matrix(~ 0 + InVivo + Cell.detection.rate + Tumor)
de <- DE_TumorAnalysis(data=so$RNA@counts[,cells], design=design[cells,])
write.table(de, 'Data/06a_de_E0771_tumor_infil_v1.0.txt', sep='\t', quote=F, row.names=F)




## Run DE for tumor-infiltration across all cell subsets and tumor-models
design <- model.matrix(~ 0 + InVivo + Cell.detection.rate + Tumor)
for(cancer in c('B16','E0771')){
    for(celltype in levels(so$Celltype)){
        cells <- which((so$Celltype %in% celltype) & (so$Cancer == cancer))
        de <- DE_TumorAnalysis(data=so$RNA@counts[,cells], design=design[cells,])
        write.table(de, paste0('Data/06b_de_',cancer,'_tumor_',celltype,'_v1.0.txt'), sep='\t', quote=F, row.names=F)
    }
}



## Run DE between NK subtypes
Time <- as.integer(so$Time == 15) # baseline is 0 and anything in vivo is >0dpi, so you only need to list which samples are 15dpi to make sure the design matrix is not co-linear.
celltypes <- c('NK1','NK2','NK3','NK4','NK5','NK6','NK7','trNK')
cells <- which((so$Celltype %in% celltypes))
de <- do.call(rbind, lapply(celltypes, function(celltype){
    Celltype_Coef <- as.integer(so$Celltype == celltype)
    design <- model.matrix(~ Cell.detection.rate + InVivo + Time + Tumor + Celltype_Coef)
    res <- DE_TumorAnalysis(data=so$RNA@counts[,cells], design=design[cells,], contrast='Celltype_Coef')
    return(data.frame(celltype=celltype, res))
}))
write.table(de, paste0('Data/06c_de_compare_NK_subtypes_v1.0.txt'), sep='\t', quote=F, row.names=F)



## Volcano plots of overall tumor-infiltration
de <- do.call(rbind, lapply(c('B16_tumor','E0771_tumor'), function(cancer){
    de.filepath=paste0('Data/06a_de_',cancer,'_tumor_infil_v1.0.txt')
    de <- read.delim(de.filepath) 
    de$logFDR <- -log10(de$adj_pval)
    de$logFDR[is.infinite(de$logFDR)] <- 1 + max(de$logFDR[!is.infinite(de$logFDR)])
    de <- de[order(de$lfc),]
    deg.label <- c(de[abs(de$lfc) > 1 & de$logFDR > 2,'name'] %>% head(., 8),
                   de[abs(de$lfc) > 1 & de$logFDR > 2,'name'] %>% tail(., 8))
    de <- de[order(de$f_statistic, decreasing=T),]
    deg.label <- c(deg.label, c('Itgam','Zeb2','Calhm2','Spn','Ly6c2','S1pr5'),
                   de[de$lfc > 1 & de$logFDR > 2,'name'] %>% head(., 8),
                   de[de$lfc < -1 & de$logFDR > 2,'name'] %>% head(., 8)) %>% unique()
    de$label <- de$name %in% deg.label
    de$cancer <- paste0(cancer,' tumor infil.')
    return(de)
}))
segmentlines <- do.call(rbind, lapply(c('B16_tumor','E0771_tumor'), function(cancer){
    de.filepath=paste0('Data/06a_de_',cancer,'_tumor_infil_v1.0.txt')
    de <- read.delim(de.filepath) 
    de$logFDR <- -log10(de$adj_pval)
    output <- data.frame(x = c(-1,-1,1,1), 
                         xend = c(min(de$lfc), -1,1,max(de$lfc)), 
                         y = c(2,2,2,2), 
                         yend = c(2,max(de$logFDR),max(de$logFDR),2),
                         celltype=paste0(cancer,' tumor infil.'))
    return(output)
}))
p <- lapply(unique(de$cancer), function(cancer){
    de.subset <- de[which(de$cancer==cancer),]
    ggplot(de.subset, aes(x=lfc, y=logFDR)) +
        geom_point_rast(data=de.subset[which(abs(de.subset$lfc) < 1 | de.subset$logFDR > 2),],raster.dpi=100, scale=0.4, color = 'gray60', fill = 'gray60') + 
        geom_segment(data=segmentlines, aes(x=x, xend=xend, y=y, yend=yend), color = 'gray50', linetype = 'longdash', size=0.1) +
        geom_point_rast(data = de.subset[(de.subset$lfc < -1) & (de.subset$logFDR > 2),],raster.dpi=100, color = 'steelblue', scale=0.5) +
        geom_point_rast(data = de.subset[(de.subset$lfc > 1) & (de.subset$logFDR > 2),],raster.dpi=100, color = 'firebrick', scale=0.5) +
        geom_text_repel(data = de.subset[which(de.subset$label),], aes(label = name), size = 2, nudge_y = -0.5, force_pull=0.01,
                        alpha = .7, max.iter=1e7,min.segment.length = 0, segment.size=0.05, segment.alpha=0.8, #label.padding = 0.15,
                        force = 100, color = 'black', max.overlaps = 30, fontface = 'italic') +
        labs(x = 'Log2 fold-change', y = 'q value (-log10)', title=cancer) + theme_test() 
}) 
p <- cowplot::plot_grid(plotlist=p, align='hv', ncol=2)
ggsave(plot=p, width = 5.5, height = 3,filename = paste0('Figures/06a_volcano_tumor_infil_v1.0.pdf'))



## Volcano plots of celltype-specific tumor-infiltration
de <- do.call(rbind, lapply(c('B16','E0771'), function(cancer){
    return(do.call(rbind, lapply(levels(so$Celltype2), function(celltype){
        de <- read.delim(paste0('Data/06b_de_',cancer,'_tumor_',celltype,'_v1.0.txt'))
        de$logFDR <- -log10(de$adj_pval)
        de$logFDR[is.infinite(de$logFDR)] <- 1 + max(de$logFDR[!is.infinite(de$logFDR)])
        de <- de[order(de$lfc),]
        deg.label <- c(de[abs(de$lfc) > 1 & de$logFDR > 2,'name'] %>% head(., 8),
                       de[abs(de$lfc) > 1 & de$logFDR > 2,'name'] %>% tail(., 8))
        de <- de[order(de$f_statistic, decreasing=T),]
        deg.label <- c(deg.label, c('Itgam','Zeb2','Calhm2','Spn','Ly6c2','S1pr5'),
                       de[de$lfc > 1 & de$logFDR > 2,'name'] %>% head(., 8),
                       de[de$lfc < -1 & de$logFDR > 2,'name'] %>% head(., 8)) %>% unique()
        de$label <- de$name %in% deg.label
        de$celltype <- paste0(celltype,': ',cancer,' infil.')
        return(de)
    })))
}))
segmentlines <- lapply(levels(so$Celltype2), function(celltype){
    do.call(rbind, lapply(c('B16','E0771'), function(cancer){
        de.filepath=paste0('Data/06b_de_',cancer,'_tumor_',celltype,'_v1.0.txt')
        de <- read.delim(de.filepath) 
        de$logFDR <- -log10(de$adj_pval)
        output <- data.frame(x = c(-1,-1,1,1), 
                             xend = c(min(de$lfc), -1,1,max(de$lfc)), 
                             y = c(2,2,2,2), 
                             yend = c(2,max(de$logFDR),max(de$logFDR),2),
                             celltype=paste0(celltype,': ',cancer,' infil.'))
        return(output)
    }))
}) %>% data.table::rbindlist() %>% data.frame()
p <- lapply(unique(de$celltype), function(celltype){
    de.subset <- de[which(de$celltype==celltype),]
    pl <- ggplot(de.subset, aes(x=lfc, y=logFDR)) +
        geom_point_rast(data=de.subset[which(abs(de.subset$lfc) < 1 | de.subset$logFDR > 2),],raster.dpi=100, scale=0.4, color = 'gray60', fill = 'gray60') + 
        geom_segment(data=segmentlines, aes(x=x, xend=xend, y=y, yend=yend), color = 'gray50', linetype = 'longdash', size=0.1) +
        geom_point_rast(data = de.subset[(de.subset$lfc < -1) & (de.subset$logFDR > 2),],raster.dpi=100, color = 'steelblue', scale=0.5) +
        geom_point_rast(data = de.subset[(de.subset$lfc > 1) & (de.subset$logFDR > 2),],raster.dpi=100, color = 'firebrick', scale=0.5) +
        geom_text_repel(data = de.subset[which(de.subset$label),], aes(label = name), size = 2, nudge_y = -0.5, force_pull=0.01,
                        alpha = .7, max.iter=1e7,min.segment.length = 0, segment.size=0.05, segment.alpha=0.2, #label.padding = 0.15,
                        force = 100, color = 'black', max.overlaps = 30, fontface = 'italic') +
        labs(x = 'Log2 fold-change', y = 'q value (-log10)', title=celltype) + theme_test() 
    return(pl)
}) 
p <- cowplot::plot_grid(plotlist=p, align='hv', ncol=6)
ggsave(plot=p, width = 15, height = 9,filename = paste0('Figures/06b_volcano_tumor_infil_byCelltype_v1.0.pdf'))



## heatmap of top tumor-infiltration genes across celltypes
de <- de[with(de, order(celltype, -f_statistic)),]
topg <- de[which((de$adj_pval < 0.01) & (abs(de$lfc) > 2) & !grepl('ILC3',de$celltype) & !grepl('NKT',de$celltype)),]
topg <- slice_head(topg, by=celltype, n=20)[,'name'] %>% unique()
d <- reshape2::dcast(de[which(de$name %in% topg & !grepl('ILC3',de$celltype) & !grepl('NKT',de$celltype)),], celltype~name, value.var='lfc',fill=0)
d <- data.frame(d[,-1], row.names=d[,1])
d[(d > 5)] <- 5
d[(d < -5)] <- -5
pheatmap::pheatmap(d, cluster.method='ward.D2', cutree_cols=6, cutree_rows=2, 
                   treeheight_row=15,treeheight_col=15, angle_col='90', fontsize_col=7,fontsize_row=8,
                   filename='Figures/06c_hm_top_deg_byCelltype_v1.0.pdf', height=3.5, width=10)



## heatmap of top NK-characteristic de genes across celltypes
de <- read.delim('Data/06c_de_compare_NK_subtypes_v3.5.txt')
de <- de[with(de, order(celltype, -f_statistic)),]
topg <- de[which((de$adj_pval < 0.01) & (abs(de$lfc) > 2)),]
topg <- slice_head(topg, by=celltype, n=20)[,'name'] %>% unique()
d <- reshape2::dcast(de[which(de$name %in% topg),], celltype~name, value.var='lfc',fill=0)
d <- data.frame(d[,-1], row.names=d[,1])
d[(d > 5)] <- 5
d[(d < -5)] <- -5
pheatmap::pheatmap(d, cluster.method='ward.D2', cutree_cols=6, cutree_rows=2, cluster_rows=F,
                   treeheight_row=15,treeheight_col=15, angle_col='90', fontsize_col=7,fontsize_row=8,
                   filename='Figures/06d_hm_top_deg_byCelltype_lfc_v1.0.pdf', height=3.5, width=12)




## Violin plots, showing decreased mature marker expression in tumor mNK cells
gene <- c('Nkg7','Zeb2','Spn','Ly6c2','S1pr5','Itgam')
DefaultAssay(so) <- 'RNA'
so$TissueAndModel <- paste0(so$Cancer, '_', so$Source) %>% stringr::str_replace(., 'Control_Culture', 'Culture')
so$TissueAndModel <- factor(so$TissueAndModel, levels=rev(c('Culture','B16_Spleen','B16_Tumor','E0771_Spleen','E0771_Tumor')))
col <- structure(c('tomato2','cyan4', 'orange2','steelblue3', 'green4'), names= c('Culture','B16_Spleen','B16_Tumor','E0771_Spleen','E0771_Tumor'))
p <- VlnPlot(so[, which((so$Cancer %in% c('B16', 'E0771', 'Control')) & (so$Time < 15) & (so$Celltype == c('mNK1')))], 
             pt.size=0, split.by='TissueAndModel', group.by='Celltype', gene, stack=T,add.noise=F) + scale_fill_manual(values=col) +
    geom_vline(xintercept=c(1:4), linetype='longdash', size=0.1, alpha=0.25)
ggsave(plot=p, width = 5, height = 5,filename = paste0('Data/06e_vln_mNK_markers_v1.0.pdf'))




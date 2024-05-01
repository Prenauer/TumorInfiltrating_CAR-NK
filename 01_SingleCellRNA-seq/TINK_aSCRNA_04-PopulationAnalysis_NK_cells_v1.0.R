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
so <- readRDS(file = 'Data/so_integrated_v1.1.RDS')


## Subset NK cells and remove myeloid cells
so <- subset(so, subset = (Full.Celltype %in% c('NK cells', 'NK-T cells')))


## Scale data
plan('sequential')
DefaultAssay(so) <- 'RNA'
so <- ScaleData(so)
DefaultAssay(so) <- 'integrated'
so <- ScaleData(so)


## Linear dim reduction
plan('multicore')
so <- RunPCA(so)
ElbowPlot(so, ndims = 30, reduction = 'pca')


## UMAP dim reduction 
ndim <- 20
plan('multicore')
so <- RunUMAP(so, reduction = 'pca', dims = 1:ndim, seed.use = 11012018)


## Find optimal cell clustering 
plan('multicore')
so <- FindNeighbors(so, reduction = 'pca', dims = 1:ndim)
# Calculate clustering metrics: within-cluster sum-of-squares (wss) and average Silhouette width (asw)
plan('sequential')
CalcWSS <- function(x,y) (length(x)-1)*(var(x)+var(y))
xy <- so@reductions$umap@cell.embeddings %>% data.frame()
d <- dist(xy)
optimalResolution <- lapply(c(0.1, seq(0.2,1.8,0.2)), function(res){
    tmp <- FindClusters(so, graph.name = 'integrated_snn', algorithm = 2, resolution = res, random.seed = 11012018)
    wss <- reframe(cbind(xy,cluster=Idents(tmp)), .by='cluster', wss=CalcWSS(UMAP_1,UMAP_2))$wss %>% sum()
    asw <- mean(cluster::silhouette(as.integer(Idents(tmp)), d)[,3])
    cat(res,'\n')
    return(data.frame(res=res, wss=wss, asw=asw))
}) %>% data.table::rbindlist()
p1 <- ggplot(optimalResolution, aes(x=res, y=asw)) + geom_point(color='tomato2') + geom_path(color='tomato2') + 
    geom_vline(xintercept=c(0.4), linetype='longdash', size=0.2, color='gray40') +
    scale_x_continuous(breaks=seq(0.4,1.6,0.4)) + 
    labs(x='Clustering resolution', y=str_wrap('Average silhouette width',20)) + theme_test()
p2 <- ggplot(optimalResolution, aes(x=res, y=wss)) + geom_point(color='darkcyan') + geom_path(color='darkcyan') + 
    geom_vline(xintercept=c(0.4), linetype='longdash', size=0.2, color='gray40') +
    scale_x_continuous(breaks=seq(0.4,1.6,0.4)) +
    labs(x='Clustering resolution', y=str_wrap('Within-cluster sum-of-squares (wss)',20)) + theme_test()
ggsave(plot = p1 + p2, filename = 'Figures/04a_qc_NK_optimal_clustering.pdf', height = 2, width = 6)


## Cluster cells
plan('sequential')
so <- FindClusters(so, graph.name = 'integrated_snn', algorithm = 2, resolution = 0.4, random.seed = 11012018)


## Heatmap of variable genes to check cluster separation (chosen from Wilcoxon test DE results)
plan('multicore')
de <- presto::wilcoxauc(so, group_by='integrated_snn_res.0.4',seurat_assay='integrated', assay='data')
de$group <- factor(de$group, levels=levels(so$integrated_snn_res.0.4))
de <- de[with(de, order(group, pval, -statistic)),]
write.table(de, file = 'Data/04a_wilcox_NK_clusters_v1.0.txt', sep = '\t', col.names = T, row.names = T)
de <- de[with(de, order(group, -auc, -statistic)),]
de <- de[which((de$padj < 0.01) & ((de$logFC) > 0.1)),]
genes <- dplyr::slice_head(de, by=group, n=20)$feature %>% unique()
DefaultAssay(so) <- 'integrated'
p <- DoHeatmap(subset(so, downsample = 200), features = genes, slot = 'scale.data', label = T, angle=90)
ggsave(plot = p, filename = 'Figures/04b_hm_wilcox_NK_clusters_v1.0.pdf', dpi = 300, width = 10, height = 2.6, units = 'in')


## Determine cluster identities
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
for(n in names(glist)) so@meta.data[,paste0(n,'_score')] <- scale(colSums(so@assays$RNA@data[glist[[n]],]))[,1]
p <- VlnPlot(so, pt.size=0, group.by='integrated_snn_res.0.4', fill.by='feature', stack=T, add.noise=F,
             paste0(c('NKT','ILC1','mNK','iNK','tissueResidency','proliferation'),'_score'))
ggsave(plot=p, width = 5, height = 5,filename = 'Figures/04c_vln_NK_signatures_clusters_v1.0.pdf')


## Rename NK populations
Idents(so) <- so[['integrated_snn_res.0.4']]
so <- RenameIdents(so, '1'='NK1','2'='NK2','3'='NK3','4'='NK6','7'='NK4','9'='NK5', '6'='NK7', '0'='mNK', '8'='NKT','10'='trNK','5'='ILC1')
so$Celltype <- Idents(so)


## Heatmap showing celltypes are transcriptionally distinct (variable genes chosen from Wilcoxon test DE results)
plan('multicore')
de <- read.delim('Data/04a_wilcox_NK_clusters_v1.0.txt')
de <- de[with(de, order(group, -auc, -statistic)),]
de <- de[which((de$padj < 0.01) & ((de$logFC) > 0.1)),]
genes <- dplyr::slice_head(de, by=group, n=20)$feature %>% unique()
DefaultAssay(so) <- 'integrated'
Idents(so) <- so$Celltype
p <- DoHeatmap(subset(so, downsample = 200), features = genes, slot = 'scale.data', label = T, angle=90)
ggsave(plot = p, filename = 'Figures/04d_hm_wilcox_NK_celltypes_v1.0.pdf', dpi = 300, width = 10, height = 2.6, units = 'in')


## Violin plot showing that celltypes match expected transcriptional signatures
p <- VlnPlot(so, pt.size=0, group.by='Celltype', fill.by='feature', stack=T, add.noise=F,
             paste0(c('NKT','ILC1','mNK','iNK','tissueResidency','proliferation'),'_score'))
ggsave(plot=p, width = 5, height = 5,filename = 'Figures/04e_vln_NK_signatures_celltypes_v1.0.pdf')


## save data
saveRDS(so, file = 'Data/so_NK_v1.0.RDS')


## Get celltype summary stats
DefaultAssay(so) <- 'RNA'
so$TimePretty <- paste0(formatC(so$Time, width=2, flag='0'),'dpi')
so$TissueModelTime <- paste0(so$Cancer, '_', so$Source, '_', so$TimePretty) %>% stringr::str_replace(., 'Control_Culture', 'Culture')
data <- data.frame(celltype = so$Celltype, tissue=so$Source, model=so$Cancer, time=so$Time)
data <- list(model=reframe(data, .by='celltype', n_control=sum(model == 'Control'), n_spleen=sum(model == 'B16F10'), n_tumor=sum(model == 'E0771') ),
             time=reframe(data, .by='celltype', n_00dpi=sum(time == 0), n_07dpi=sum(time == 7), n_15dpi=sum(time == 15) ),
             tissue=reframe(data, .by='celltype', n_control=sum(tissue == 'Culture'), n_spleen=sum(tissue == 'Spleen'), n_tumor=sum(tissue == 'Tumor') ))
data <- cbind(data$model, data$time[,-1], data$tissue[,-1])
write.table(data, 'Data/04b_summary_stats_NK_cell.txt',sep = '\t', quote=F, row.names=F)








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
so <- readRDS(file = 'Data/so_integrated_v1.0.RDS')


## Dimensional reduction
plan('multicore')
so <- RunPCA(so)
ElbowPlot(so, ndims = 50, reduction = 'pca')
ndim <- 27
so <- RunUMAP(so, reduction = 'pca', dims = 1:ndim, seed.use = 11012018)


## Cluster cells
plan('sequential')
so <- FindClusters(so, graph.name = 'integrated_snn', algorithm = 2, resolution = 0.4, random.seed = 11012018)
Idents(so) <- so@meta.data$integrated_snn_res.0.4
so <- FindSubCluster(so, cluster = '8', graph.name = 'integrated_snn', subcluster.name = 'sub.cluster.int', algorithm = 2, resolution = 0.1)
Idents(so) <- so@meta.data$sub.cluster.int


## Wilcox test between clusters
plan('multicore')
Idents(so) <- so[['sub_cluster_int']]
Idents(so) <- factor(Idents(so), levels = c('0','1','2','3','4','5','6','7','8_0','8_1','9','10','11','12','13','14','15'))
de <- RunPrestoAll(so, logfc.threshold = 0.5, only.pos = F, min.pct = 0.2, max.cells.per.ident = 2000, slot = 'scale.data', assay = 'integrated')
write.table(de, file = 'Data/03a_wilcox_AllCells_clusters_v1.0.txt', sep = '\t', col.names = T, row.names = T)


## Heatmap of variable genes to check cluster separation (chosen from Wilcoxon test DE results)
de <- read.delim(file = 'Data/03a_wilcox_AllCells_clusters_v1.0.txt')
Idents(so) <- factor(so[['sub_cluster_int']], levels = c('10','11','12','13','14','15','0','1','2','3','4','5','6','7','8_0','8_1','9'))
features <- de[(de$p_val_adj < 0.01) & (abs(de$avg_diff) > 1) & (de$pct.1 > 0.1),]
features <- dplyr::slice_head(features, by='cluster', n=100)$gene %>% unique()
DefaultAssay(so) <- 'integrated'
p <- DoHeatmap(subset(so, downsample = 200, group.by='Full.Celltype'), group.by='Full.Celltype', features = features, slot = 'scale.data', label = F)
ggsave(plot = p, filename = 'Figures/03a_hm_wilcox_AllCells_clusters_v1.0.pdf', dpi = 300, width = 9, height = 5, units = 'in')


## Violin and Feature plot for cell type identification
DefaultAssay(so) <- 'RNA'
genes <- c('Ptprc','Ly6g','Cd3e','Cd4','Cd8a','Cd19','Sdc1','Cd14','Adgre1','H2-Aa','Ncr1','Gypa','Pmel','Mlana')
p.vln <- VlnPlot(so, genes, adjust = 40, pt.size = 0, slot = 'data', x, log = T, combine = T, stack = T, fill.by = 'ident') +
        geom_vline(xintercept = 2, colour = "grey50", linetype = 'longdash') + theme(legend.position = 'none')
ggsave(plot = p.vln, file = 'Figures/03b_vln_All_Celltypes_Markers_adjusted_v8.pdf', unit = 'in', width = 8, height = 6)
p.feat <- FeaturePlot(so, order = T, min.cutoff = 'q01',max.cutoff = 'q99',pt.size = 1, genes, raster = T) + plot_layout(ncol = 3)
ggsave(plot = p.feat, file = 'Figures/03c_umap_AllCells_clusters_markers_v1.0.pdf', unit = 'in', width = 3.5*3, height = 3*ceiling((n.vln)/3))


## Rename populations
markers <- c('Cd19','Sdc1','Cd3e','Cd14','Adgre1','Hbb-bs','Gypa','Pmel','Mlana','H2-Aa','Ly6g','Ptprc')
so <- RenameIdents(so, '0' = 'NK cells', '1' = 'NK cells', '2' = 'NK cells', '3' = 'NK cells', 
                   '4' = 'NK cells', '5' = 'Macropphage', '6' = 'B cells', '7' = 'NK cells', '8_0' = 'CD4 T cells', 
                   '8_1' = 'CD8 T cells', '9' = 'NK-T cells', '10' = 'Erythrocyte', '11' = 'Monocyte', '12' = 'NK cells',
                   '13' = 'Cancer cells', '14' = 'Neutrophils', '15' = 'Plasma Cells')
so[['Full_Celltype']] <- Idents(so)


## Export metadata
meta <- cbind(so@meta.data[,c('nCount_RNA','nFeature_RNA','Sample','Time','Source','Cancer','Full_Celltype')],so@reductions$umap@cell.embeddings)
write.table(meta, 'Data/03b_metadata_allCells_v1.0.txt',sep = '\t')


## save data
saveRDS(so, file = 'Data/so_integrated_v1.1.RDS')


## Heatmap of variable genes with celltype labels
de <- read.delim(file = 'Data/03a_wilcox_AllCells_clusters_v1.0.txt')
Idents(so) <- factor(so[['sub_cluster_int']], levels = c('10','11','12','13','14','15','0','1','2','3','4','5','6','7','8_0','8_1','9'))
genes <- de[(de$p_val_adj < 0.01) & (abs(de$avg_diff) > 1) & (de$pct.1 > 0.1),]
genes <- dplyr::slice_head(genes, by='cluster', n=100)$gene %>% unique()
DefaultAssay(so) <- 'integrated'
p <- DoHeatmap(subset(so, downsample = 200, group.by='Full.Celltype'), group.by='Full.Celltype', features = genes, slot = 'scale.data', label = F)
ggsave(plot = p, filename = 'Figures/03d_hm_wilcox_AllCells_celltypes_v1.0.pdf', dpi = 300, width = 9, height = 5, units = 'in')


## Violin and Feature plot with celltype labels
DefaultAssay(so) <- 'RNA'
genes <- c('Ptprc','Ly6g','Cd3e','Cd4','Cd8a','Cd19','Sdc1','Cd14','Adgre1','H2-Aa','Ncr1','Gypa','Pmel','Mlana')
p.vln <- VlnPlot(so, genes, adjust = 40, pt.size = 0, slot = 'data', x, log = T, combine = T, stack = T, fill.by = 'ident') +
    geom_vline(xintercept = 2, colour = "grey50", linetype = 'longdash') + theme(legend.position = 'none')
ggsave(plot = p.vln, file = 'Figures/03e_vln_AllCells_celltypes_markers_v1.0.pdf', unit = 'in', width = 8, height = 6)
p.feat <- FeaturePlot(so, order = T, min.cutoff = 'q01',max.cutoff = 'q99',pt.size = 1, genes, raster = T) + plot_layout(ncol = 3)
ggsave(plot = p.feat, file = 'Figures/03f_umap_AllCells_celltypes_markers_v1.0.pdf', unit = 'in', width = 3.5*3, height = 3*ceiling((n.vln)/3))


## Visualize UMAP by dataset
data <- data.frame(so@reductions$umap@cell.embeddings, Celltype=so$Full_Celltype, Sample_Type=paste0(so$Cancer,'_',so$Source,'_',so$Time))
data <- rbind(data, data.frame(data[,1:3], Sample_Type='Merged_Datasets'))
p.celltype <- ggplot(data, aes(x=UMAP_1, y=UMAP_2, color=Celltype)) + ggrastr::geom_point_rast(scale=0.25, alpha=0.5, raster.dpi = 300) + facet_wrap(vars(Sample_Type)) + theme_test()
ggsave(plot = p.celltype, file = 'Figures/03g_split_umap_AllCells_celltypes_markers_v1.0.pdf', unit = 'in', width = 7, height = 6)






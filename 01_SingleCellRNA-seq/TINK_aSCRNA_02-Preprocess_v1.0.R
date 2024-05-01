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

## Load the PBMC dataset
data <- Read10X(data.dir = "outs/count/raw_feature_bc_matrix")
## Initialize the Seurat object
so <- CreateSeuratObject(counts = data, project = "TINK", min.cells = 3, min.features = 200)


## Annotate dataset
meta <- read.delim('Data/TINK_aSCRNA_02-Preprocess_sample_data.txt')
so[['Sample']] <- meta[stringr::str_split(colnames(so), '-', simplify = T)[,2] %>% as.integer(),'Sample']
so[['Time']] <- meta[stringr::str_split(colnames(so), '-', simplify = T)[,2] %>% as.integer(),'Time']
so[['Source']] <- meta[stringr::str_split(colnames(so), '-', simplify = T)[,2] %>% as.integer(),'Source']
so[['Cancer']] <- c('B16', 'Control', 'E0771')[(factor(substr(so[['Sample']]$Sample,0,5)) %>% as.integer())]
so[['Cell_detection_rate']] <- log10(so@meta.data$nFeature_RNA)
so[['Percent_mt']] <- PercentageFeatureSet(so, pattern = '^mt-')
so[['Low_Quality']] <- PercentageFeatureSet(so, features = 'Kcnq1ot1')
so[['rRNA_Contam']] <- PercentageFeatureSet(so, features = c('Gm26917','Gm42418'))
so[['B16']] <- PercentageFeatureSet(so, features = 'Pmel')


## Visualize QC metrics
p1 <- VlnPlot(so, pt.size=0, features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "Low_Quality", "rRNA_Contam"), stack=T)
p2 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "percent_mt") + 
p3 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(p1 + p2 + p3, height = 3, width = 9, filename = 'Figures/02a_qc_metrics.pdf')


## Filter
so <- subset(so, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & Percent_mt < 5)
so <- subset(so, subset = Low_Quality < 0.1 & rRNA_Contam < 5)


## norm, scale
so <- SplitObject(so, split.by = 'Sample')
so <- lapply(so, function(x){
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000) 
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(so)
so <- lapply(so, function(x){
  x <- ScaleData(x, features = features) 
  x <- RunPCA(x, features = features)
})


## Integrate samples
plan('sequential')
anchors <- FindIntegrationAnchors(object.list = so, anchor.features = features, reduction = 'rpca', k.anchor = 20)
saveRDS(anchors, file = 'Data/integration_anchors_v6.RDS')
so <- IntegrateData(anchorset = anchors)
DefaultAssay(so) <- "RNA"
so <- ScaleData(so)
DefaultAssay(so) <- "integrated"
so <- ScaleData(so)


## save data
saveRDS(so, file = 'Data/so_integrated_v1.0.RDS')





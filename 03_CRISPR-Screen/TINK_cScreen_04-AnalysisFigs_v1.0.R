#!/usr/bin/env Rscript
library(plyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(gprofiler2)
library(RCurl)
library(ggrepel)
library(pheatmap)
library(limma)
library(edgeR)
library(metap)
library(glmGamPoi)
library(RColorBrewer)
library(future.apply)
require(ggrastr)
library(ggplotify)
options(ggrastr.default.dpi=750)


setwd('Screen')


## Make QC plots for screen sample variability
tumor_models <- c('B16F10','E0771','GL261','Pan02')


## Compare screen results to expressed genes (ImmGen NK cells)
immgen <- read.delim('/Users/pren/Dropbox/Private/TINK/BulkRNAaseq/GSE122597_Gene_count_table_GENCODE_vM25.csv',row.names = 1, sep = ',')
immgen <- apply(immgen[,-ncol(immgen)],2, function(x) { 1e6 * x/sum(x)}) %>% rowMeans()
immgen <- log2(immgen + 1)
res <- do.call(rbind, lapply(tumor_models, function(tumor_model) { 
    res <- read.delim(file.path(DIR_DATA,paste0('Mageck/mageck_',tumor_model,'_results.txt.gene_summary.txt'))) 
    res <- res[,c('id','pos.p.value','pos.score')]
    res$zscore <- scale(-res$pos.score)[,1]
    res$lcpm <- immgen[res$id]
    res <- res[!is.na(res$lcpm),]
    res$model <- tumor_model
    return(res)
}))
p <- ggplot(res, aes(x = zscore, y = lcpm, color = zscore)) +
    geom_point(alpha = 0, aes(size=0.5* -log10(res$pos.p.value))) + 
    geom_point(color = 'gray40', aes(size=1* -log10(res$pos.p.value))) + 
    geom_point(color = 'white', aes(size=0.5* -log10(res$pos.p.value))) + 
    geom_point(alpha = .8, aes(size=0.5* -log10(res$pos.p.value))) + 
    scale_colour_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))(100), limits=range(res$zscore)) +
    scale_size(range=c(0.5,4)) + ylim(c(0,14)) + 
    geom_vline(xintercept = 1, color = 'gray50', linetype = 'dashed') + 
    labs(x = 'Gene score', y = 'ImmGen gene expr. (log2)') + 
    guides(size=guide_legend(title=str_wrap('p value (-log10)',15))) +
    cowplot::theme_cowplot() + 
    facet_wrap(vars(model), ncol=2)
ggsave(plot = p, '03a_Immgen_comparison_v1.0.pdf', height = 5, width = 7)









#!/usr/bin/env Rscript
library(plyr)
library(dplyr)
library(stringr)

setwd('Screen')

## Prepare data for analysis
ct <- read.delim('NKscreen_CountTable.txt')[,-2]
meta <- read.delim('TINK_Screen_Metadata.txt')
meta$Sample_ID <- paste0(meta$Sample_ID,'_',meta$Source)
meta$Tumor <- str_split(meta$Sample_ID, '_', simplify = T)[,1]
colnames(ct) <- c('sgRNA','Gene', meta$Sample_ID)
newSgNames <- stringr::str_split(ct$sgRNA, '_', simplify = T)
newSgNames <- paste0(newSgNames[,1], '_sg', newSgNames[,2])
ct$sgRNA <- newSgNames
rownames(ct) <- ct$sgRNA
rownames(meta) <- meta$Sample_ID
# Set up Metadata of filtered data
meta$TotalReads <- colSums(ct[,meta$Sample_ID])
meta$batch <- c(rep(1,46),rep(2,43))


## Analyze enrichment in each tumor model
tumor_models <- c('B16F10','E0771','Pan02','GL261')
lapply(tumor_models, function(tumor_model){ 
	# remove low-count samples
	lowCountSamples <- colnames(ct[,-c(1:2)])[colSums(ct[,-c(1:2)]) < 100000]

	model <- grep(tumor_model, colnames(ct), value = T)
	controls <- grep('Ctrl', colnames(ct), value = T) %>% setdiff(y = lowCountSamples)
	tumor <- grep('Tumor', colnames(ct), value = T) %>% intersect(model) %>% setdiff(y = lowCountSamples)
	spleen <- grep('Spleen', colnames(ct), value = T) %>% intersect(model) %>% setdiff(y = lowCountSamples)
	samples <- c(tumor, controls)

	# remove low-count guides
	filterByCount <- lib$sgRNA[rowSums(ct[,samples]) > 2]
	filterBySingleGuide <- plyr::count(lib[filterByCount,'Gene'])
	filterBySingleGuide <- filterBySingleGuide$x[which(filterBySingleGuide$freq > 1)]
	filterBySingleGuide <- lib$sgRNA[which(lib$Gene %in% filterBySingleGuide)]
	filter.sg <- intersect(filterBySingleGuide, filterByCount)
	
	df <- apply(ct[filter.sg, samples], 2, function(x) 1e6*x/sum(x))
	df <- cbind(lib[filter.sg,], df)
	ntc <- intersect(filter.sg, lib$sgRNA[which(lib$Gene == 'NTC')])
	
	s.test <- tumor %>% paste(collapse = ',')
	s.ctrl <- controls %>% paste(collapse = ',')
	f.input <- paste0('Mageck/mageck_', tumor_model, '_counts.txt')
	f.output <- paste0('Mageck/mageck_', tumor_model, '_results.txt')
	f.ntc <- paste0('Mageck/mageck_', tumor_model, '_ntc.txt')
	cmd <- paste0('mageck test --skip-gene NTC', 
		' -k ', f.input, ' -t ', s.test, ' -c ', s.ctrl, ' --norm-method none --remove-zero-threshold 2',
		' --control-sgrna ', f.ntc, ' -n ', f.output)
	
	write.table(df, f.input, sep = '\t', row.names = F, quote = F)
	write.table(ntc, f.ntc, sep = '\t', row.names = F, quote = F)
	cat(paste0(cmd,'\n','\n'))
})



## Compare screen results to expressed genes (ImmGen NK cells)
immgen <- read.delim('TINK_cScreen_GSE122597_Gene_count_table_GENCODE_vM25.csv',row.names = 1, sep = ',')
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
    scale_colour_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))(100), 
                           limits=range(res$zscore)) +
    scale_size(range=c(0.5,4)) + ylim(c(0,14)) + 
    geom_vline(xintercept = 1, color = 'gray50', linetype = 'dashed') + 
    labs(x = 'Gene score', y = 'ImmGen gene expr. (log2)') + 
    guides(size=guide_legend(title=str_wrap('p value (-log10)',15))) +
    cowplot::theme_cowplot() + 
    facet_wrap(vars(model), ncol=2)
ggsave(plot = p, '03a_Immgen_comparison_v1.0.pdf', height = 5, width = 7)



## Venn of top screen hits
# immgen data
immgen <- read.delim('TINK_cScreen_GSE122597_Gene_count_table_GENCODE_vM25.csv',row.names = 1, sep = ',')
immgen <- apply(immgen[,-ncol(immgen)],2, function(x) { 1e6 * x/sum(x)}) %>% rowMeans()
coding.genes <- read.delim('TINK_cScreen_gencode.vM25.gene_anno.txt')
coding.genes <- coding.genes$gene_name[which(coding.genes$gene_type == 'protein_coding')]
genes.in.library <- unique(read.delim('NKscreen_CountTable.txt')$Gene)
immgen <- immgen[which(names(immgen) %in% coding.genes)] 
immgen <- immgen[immgen > 0] %>% names()
# screen data
screen <- do.call(rbind, lapply(tumor_models, function(tumor_model) { 
    res <- read.delim(paste0('Mageck/mageck_',tumor_model,'_results.txt.gene_summary.txt')) 
    res <- res[,c('id','pos.p.value','pos.score')]
    res$zscore <- scale(-res$pos.score)[,1]
    res$model <- tumor_model
    return(res)
}))
screen <- screen[which((screen$pos.p.value < 0.05) & (screen$zscore > 1)),]
screen <- plyr::count(screen$id)
screen <- screen$x[which(screen$freq >= 2)]
# tumor DE genes
deg.b16 <- read.delim('Data/06a_de_B16_tumor_infil_v1.0.txt')
deg.b16 <- deg.b16$name[which((deg.b16$adj_pval < 0.01) & (deg.b16$lfc < -1))]
deg.e0771 <- read.delim('Data/06a_de_E0771_tumor_infil_v1.0.txt')
deg.e0771 <- deg.e0771$name[which((deg.e0771$adj_pval < 0.01) & (deg.e0771$lfc < -1))]
data.venn <- list(Immgen=immgen, Screen=screen, B16_dep=deg.b16, E0771_dep=deg.e0771)
p.venn <- venn::venn(data.venn, ellipse = T, zcolor = c('tomato2','seagreen3','gold2','gray'), ggplot = T)
data.venn <- venn::venn(data.venn, ellipse = T, zcolor = c('tomato2','seagreen3','gold2','gray'), ggplot = F)
data.venn <- attr(data.venn,'intersections')
data.venn <- data.frame(Comparison=names(data.venn), Genes=unlist(lapply(data.venn, paste0, collapse=',')))
ggsave(plot = p.venn, height = 3.5, width = 3.5, filename = '03b_venn_screenhits_v1.0.pdf')
write.table(data.venn, '03b_venn_screenhits_v1.0.txt', sep='\t', row.names=F,quote=F)



## Make plots for Mageck analysis results
tumor_models <- c('B16F10','E0771','Pan02','GL261')
p <- lapply(tumor_models, function(tumor_model){
    mageck <- read.delim(paste0('Mageck/mageck_',tumor_model,'_results.txt.gene_summary.txt'))
    mageck <- data.frame(gene = mageck$id,
                         logP = -log10(mageck$pos.p.value),
                         lfc = scale(-mageck$pos.lfc),
                         rank = mageck$pos.rank)
    mageck <- mageck[with(mageck, order(rank, decreasing = F)),]
    glist1 <- mageck$gene[which(mageck$logP > -log10(0.05))]
    glist.label <- mageck[with(mageck, order(rank)),'gene'] %>% head(8)
    glist.label <- c(glist.label, 'Calhm2')
    ggplot(mageck, aes(x = rank, y = logP, label = gene)) + 
        ggrastr::geom_point_rast(color = 'gray50') + 
        ggrastr::geom_point_rast(data = mageck[which(mageck$gene %in% glist1),], aes(color = rank)) + 
        scale_colour_viridis_c(guide='none', option ='turbo', direction = -1) + 
        ggrepel::geom_text_repel(data = mageck[which(mageck$gene %in% glist.label),], segment.size=0.2, nudge_x=1,
                                 min.segment.length = 0, segment.color = 'gray70', max.overlaps = 30, fontface='italic') +
        labs(x = 'Gene rank', y = 'p value (-log10)', title = paste0(tumor_model)) +
        cowplot::theme_cowplot() +
        theme(legend.position = 'none') %>% 
        return()
})
p <- cowplot::plot_grid(plotlist = p, ncol = 4)
ggsave(plot = p, filename = '03d_mageck_plots_v1.0.pdf', height = 2.5, width = 10)





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


for(gp in c('unstim','stim')){
    ## Set up DGElist
    if(gp == 'unstim') data <- catchKallisto(paste0('RawData/Sample',1:6))
    if(gp == 'stim') data <- catchKallisto(paste0('RawData/Sample',7:12))
    groups <- c(rep('AAVS1',3), rep('CALHM2',3))
    rownames(data$counts) <- str_split(rownames(data$counts), '\\.', simplify = T)[,1]
    colnames(data$counts) <- paste0(groups,'_',c(1:3,1:3))
    anno <- read.delim('Homo_sapiens.GRCh38.96.t2g.txt', header=F, col.names=c('ENST','ENSG','GeneSymbol'))
    anno <- anno[which(anno$ENST %in% rownames(data$counts)),]
    rownames(data$annotation) <- str_split(rownames(data$annotation), '\\.', simplify = T)[,1]
    anno <- rbind(anno, data.frame(ENST=setdiff(rownames(data$counts), anno$ENST), ENSG=NA, GeneSymbol=NA))
    rownames(anno) <- anno$ENST
    anno <- anno[rownames(data$counts),]
    dge <- DGEList(counts = data$counts, group = factor(groups), genes=anno[,3, drop=F], remove.zeros = T)
    
    
    ## Pre-process
    #filter low-expressed transcripts
    group <- dge$samples$group
    lcpm.raw <- cpm(dge, log = T)
    dge <- dge[filterByExpr(dge, group=dge$samples$group),, keep.lib.sizes=F]
    lcpm.filt <- cpm(dge, log = T)
    #get cutoff
    L <- mean(dge$samples$lib.size) * 1e-6
    M <- median(dge$samples$lib.size) * 1e-6
    lcpm.cutoff <- log2(2/M + 1/L)
    #plot filtered vs unfiltered
    p1 <- ggplot(melt(lcpm.raw)) + geom_density(aes(x = value, color = Var2)) + 
        geom_vline(xintercept = lcpm.cutoff, color = 'gray50', linetype = 2) + 
        labs(title = 'Raw data',x = 'Log-cpm', y = 'Density', color = 'Sample') + theme_classic() + 
        theme(legend.position = c(0.75, 0.75), legend.key.size = unit(0.1, 'in'), 
              legend.title = element_text(size=8),legend.text = element_text(size=6), 
              legend.background = element_rect(fill = "white", color = "black")) 
    p2 <- ggplot(melt(lcpm.filt)) + geom_density(aes(x = value, color = Var2)) + 
        geom_vline(xintercept = lcpm.cutoff, color = 'gray50', linetype = 2) + 
        labs(title = 'Filtered data',x = 'Log-cpm', y = 'Density', color = 'Sample') + theme_classic() + 
        theme(legend.position = c(0.75, 0.75), legend.key.size = unit(0.1, 'in'), 
              legend.title = element_text(size=8),legend.text = element_text(size=6), 
              legend.background = element_rect(fill = "white", color = "black")) 
    ggsave(plot = p1|p2, filename = paste0('Figures/02a_qc_prNK_filtering_',gp,'_v1.0.pdf'), height = 2.8, width = 4.2)
    rm(lcpm.cutoff, L,M, lcpm.raw, lcpm.filt)
    
    
    ## Normalize
    dge <- calcNormFactors(object = dge, method = "TMM")
    
    
    ## Set design
    CALHM2KO <- as.integer(dge$samples$group == 'CALHM2')
    design <- model.matrix(~ CALHM2KO)
    
    
    # Estimate dispersion
    dge <- estimateGLMCommonDisp(dge, design = design)
    dge <- estimateGLMTrendedDisp(dge, method = 'bin.loess', span = 1/3, design = design)
    dge <- estimateGLMTagwiseDisp(dge, design = design)
    
    
    # Fit model and run QLF test
    fit <- glmFit(dge, design = design)
    lrt <- glmLRT(fit, coef = 'CALHM2KO')

    
    # Export results
    res <- topTags(lrt, n = Inf)$table
    res$ENST <- rownames(res)
    res <- res[with(res, order(LR, decreasing = T)),]
    de <- res[,c('ENST','GeneSymbol','logCPM','LR','logFC','PValue','FDR')]
    colnames(de) <- c('enst', 'name','logCPM','statistic','lfc','pval','adj_pval')
    write.table(de, paste0('Data/02a_de_primaryNK_CALHM2KO_',gp,'_v1.1.txt'), sep = '\t', col.names = T, row.names = F, quote = F)
    
    
    ## Check sample correlation
    lcpm <- cpm(dge, log = T, normalized.lib.sizes = T) %>% data.frame()
    corr <- cor(lcpm, method = 'spearman')
    pheatmap::pheatmap(corr, scale = 'none', treeheight_row = 12, cluster_cols = T, 
                       treeheight_col = 12,  clustering_method = 'complete', 
                       angle_col = 90, width = 4, height = 3, 
                        filename = paste0('Figures/02b_hm_prNK_correlation_heatmap_',gp,'_v1.1.pdf'))
    write.table(corr, file = paste0('Data/02b_hm_prNK_correlation_heatmap_',gp,'_v1.1.txt'),sep = '\t', quote = F)
    pdf(width = 4, height = 4,  file = paste0('Figures/02c_mds_prNK_',gp,'_v1.1.pdf'))
    ggcolor <- function(n) hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[1:n]
    col <- ggcolor(length(unique(group)))[as.integer(group)]
    lcpm <- cpm(dge, log = T, normalized.lib.sizes = T, prior.count = 1) %>% data.frame()
    limma::plotMDS(lcpm, top = 1000, gene.selection = 'pairwise', col = col)
    dev.off()
    
    
    ## Export raw counts
    lcpm <- cpm(dge, log = T)
    lcpm <- data.frame(lcpm)
    lcpm$ENST <- rownames(lcpm)
    lcpm <- merge(anno, lcpm, by = 'ENST', all.x = F, all.y = T)
    write.table(lcpm, paste0('Data/02c_lcpm_prNK_',gp,'_v1.1.txt'), sep = '\t', col.names = T, row.names = F, quote = T)
    lcpm <- dge$counts
    lcpm <- data.frame(lcpm)
    lcpm$ENST <- rownames(lcpm)
    lcpm <- merge(anno, lcpm, by = 'ENST', all.x = F, all.y = T)
    write.table(lcpm, paste0('Data/02c_counts_prNK_',gp,'_v1.1.txt'), sep = '\t', col.names = T, row.names = F, quote = T)
    
    
    # Volcano plot
    res <- de
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
    p <- ggplot(res, aes(x = lfc, y = logq)) + 
        ggrastr::geom_point_rast(alpha = .25, size = 2.5, color = 'gray10', shape = 16) +
        ggrastr::geom_point_rast(data = res.de, alpha = 1, size = 2.5, color = 'gray10', shape = 16) +
        ggrastr::geom_point_rast(data = res.de, alpha = 1, color = 'white', shape = 16) +
        ggrastr::geom_point_rast(data = res.de[which(res.de$lfc > 0),], alpha = 1, color = 'firebrick', shape = 16) +
        ggrastr::geom_point_rast(data = res.de[which(res.de$lfc < 0),], alpha = 1, color = 'steelblue', shape = 16) +
        geom_segment(data = signif.lines, aes(x = x, xend = xend, y = y, yend = yend), 
                     linetype = 'dashed', color = 'gray20') +
        ggrepel::geom_text_repel(data = res.label.small, aes(label = name), alpha = 0.8, size = 2.5, nudge_y = 2,force=5, segment.size=0.04,
                                 max.iter = 100000,min.segment.length = 0, fontface = 'italic', max.overlaps = 100) +
        labs(x = 'Log2 fold-change', y = 'q value (-log10)') +
        cowplot::theme_cowplot()
    ggsave(plot = p, filename = paste0('Figures/02g_volcano_prNK_',gp,'_v1.1.pdf'), height = 4, width = 4)
}



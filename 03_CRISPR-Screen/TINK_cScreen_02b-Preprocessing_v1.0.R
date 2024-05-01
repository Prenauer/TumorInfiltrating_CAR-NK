#!/usr/bin/env Rscript

library('stringr')
library('dplyr')

countFileList <- c('NKscreen_S5_trim_align_count.txt', 'NKscreen_S6_trim_align_count.txt', 
                   'NKscreen_S7_trim_align_count.txt', 'NKscreen_S8_trim_align_count.txt', 
                   'NKscreen_S1_trim_align_count.txt', 'NKscreen_S2_trim_align_count.txt', 
                   'NKscreen_S3_trim_align_count.txt', 'NKscreen_S4_trim_align_count.txt') 
countFileNames <- c('NKscreen_F1R1', 'NKscreen_F1R2', 'NKscreen_F1R3', 'NKscreen_F1R4', 
                    'NKscreen_F2R1', 'NKscreen_F2R2','NKscreen_F2R3', 'NKscreen_F2R4')

meta <- read.delim('/Screen/TINK_Screen_Metadata.txt')
meta$Filename <- paste0('/Screen/AlignedFastq/',meta$Filename)

ct <- read.delim(meta$Filename[1],header = F, sep = ' ', stringsAsFactors = F)
for(i in 2:nrow(meta)) 
    ct <- merge(ct, read.delim(meta$Filename[i],header = F, sep = ' ', stringsAsFactors = F),by = 'V1', all = T)
ct <- merge(read.delim('/Screen/TINK_cScreen_mSurfaceomeV2_libInfo_mageckFormat.txt', header = F), ct, by = 'V1', all = T)
for(i in 4:ncol(ct)) ct[is.na(ct[,i]),i] <- 0
colnames(ct) <- c('sgRNA','Spacer','Gene', meta$Sample_ID)

write.table(ct, file = '/Screen/NKscreen_CountTable.txt', sep = '\t', col.names = T, row.names = F, quote = F)


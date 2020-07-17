setwd("~/Documents/NYUAD/Zf_ATF6/Atf6_target_ATAC")
library(ChIPpeakAnno)
library(ChIPseeker)
### input ATAC peaks, mouse and zebrafish
ATAC <- readPeakFile('xxx/ATAC-peaks.broadPeak', head=F)
seqlevels(ATAC) <-  paste('chr', seqlevels(ATAC), sep='')

### atf6 gene list
atf6_tar <- read.delim('Atf6PutTargets_all_1498.csv', sep = ',')
### convert ID to entrez
# library('biomaRt')
# listMarts()    # to see which database options are present
# mart_zf <- useDataset('drerio_gene_ensembl', useMart('ensembl'))
# listFilters(mart_zf) # RefSeq mRNA ID(s) [e.g. NM_001001398]
# listAttributes(mart_zf)
# atf6_tar_RefID <- getBM(attributes=c("ensembl_gene_id", "refseq_mrna", 'zfin_id_symbol'),
#                        filters='ensembl_gene_id',values=atf6_tar$Zebrafish_Ensembl_ID,
#                        mart=mart_zf) 
library(org.Dr.eg.db)
keytypes(org.Dr.eg.db)

mapIds(org.Dr.eg.db, keys = as.character(atf6_tar$Zebrafish_Ensembl_ID), 
                         keytype = "ENSEMBL", column="ENTREZID")
# remove dplyr
map <- select(org.Dr.eg.db, as.character(atf6_tar$Zebrafish_Ensembl_ID), 
              c("SYMBOL", "REFSEQ"), "ENSEMBL")

atf6_tar_RefID <- map$REFSEQ
### get TSS from targe gene list
library(TxDb.Drerio.UCSC.danRer10.refGene)
txdb = TxDb.Drerio.UCSC.danRer10.refGene
tss <- promoters(txdb, upstream = 0, downstream = 1)
PR <- promoters(txdb, upstream = 2000, downstream = 2000)
tss_atf6_tar <- tss[tss$tx_name %in% atf6_tar_RefID, ]
tss_atf6_pr <- PR[PR$tx_name %in% atf6_tar_RefID, ]
### check ATAC +/-
idx <- findOverlaps(tss_atf6_pr, larva_ATAC)
tss_atf6_tar$score <- FALSE
tss_atf6_tar$score[unique(queryHits(idx))] <- TRUE

### calculate matirx
library(EnrichedHeatmap)
col_fun = colorRamp2(c(0, 0.99), c("grey", "red"))
# ma <- normalizeToMatrix(ATAC, tss_atf6_tar, value_column = "V5", smooth = TRUE,
#                        extend=2000, mean_mode="w0", w=50, background=0)
# EnrichedHeatmap(ma, name='ATAC')
# sum(width(ATAC))/1369631918
ma <- normalizeToMatrix(larva_ATAC, tss_atf6_tar, value_column = "V5", smooth = TRUE,
                        extend=2000, mean_mode="w0", w=50, background=0)

# EnrichedHeatmap(ma, name='larva_ATAC')
partition = tss_atf6_tar$score
EnrichedHeatmap(ma, name="ATAC", column_title='ATAC_ZF', col=col_fun, split = -partition)
# sum(width(larva_ATAC))/1674207132
# sum(width(ATAC))/1674207132
###============================================================================================
### checck gene expression
### extract gene list
counts_df <- read.table('./atf6_con_Exp-reads.txt')
atf6_genes_exp <- subset(counts_df, rownames(counts_df) %in% map$ENSEMBL)
head(map)
### add ATAC + or - to map than combine with gene expression
df_refseq_atac <- data.frame(tss_atf6_tar$tx_name, tss_atf6_tar$score)
head(df_refseq_atac)
map_atac <- merge(map, df_refseq_atac, keep.x=all,
                  by.x='REFSEQ', by.y='tss_atf6_tar.tx_name')
atf6_genes_exp$ENSEMBL <- rownames(atf6_genes_exp)
### merge gene exp and atac 
atf6_genes_exp_atac <- merge(atf6_genes_exp, map_atac, keep.x=all, by='ENSEMBL')
colnames(atf6_genes_exp_atac)[10] <- 'atac_score'
atf6_genes_exp_atac$atf6Avg <- apply(atf6_genes_exp_atac[ ,2:4], 1, mean)
atf6_genes_exp_atac$sibAvg <- apply(atf6_genes_exp_atac[ ,5:7], 1, mean)
### remove duplicated rows
library(dplyr)
atf6_genes_exp_atac <- atf6_genes_exp_atac %>% distinct(ENSEMBL, .keep_all = TRUE)
### Check genes filtered out in DESeq2 are not included in expression file 
idx_miss <- which(!(atf6_tar$Zebrafish_Ensembl_ID %in% atf6_genes_exp_atac$ENSEMB))
### there 22, 0 expression, atac -
sup_df <- data.frame(ENSEMBL = atf6_tar[idx_miss, ]$Zebrafish_Ensembl_ID,
                     nATf6livr_1=rep(0, 22), nATf6livr_2=rep(0, 22), nATf6livr_3=rep(0, 22),
                     sibslivr_1=rep(0, 22), sibslivr_2=rep(0, 22), sibslivr_3=rep(0, 22),
                     REFSEQ=rep(0, 22), SYMBOL=rep(0, 22), atac_score=rep(FALSE, 22),
                     atf6Avg=rep(0, 22), sibAvg=rep(0, 22))
new_atf6_genes_exp_atac <- rbind(atf6_genes_exp_atac, sup_df)
### boxplot gene exp with atac
library(ggplot2)
library(ggpubr)
p <- ggplot(new_atf6_genes_exp_atac, aes(x = atac_score, y = log(sibAvg + 1))) +
  geom_violin(trim = FALSE) +
  geom_jitter(color = "tomato",position=position_jitter(0.1))  +
  stat_summary(fun=median, color="black", geom='point')

p + stat_compare_means(method = "t.test")

### plot atac - between wt and atf6
atac_pos_exp <- subset(atf6_genes_exp_atac, atf6_genes_exp_atac$atac_score=='TRUE')
attach(atac_pos_exp)
plot(log(sibAvg+1), log(atf6Avg+1))
plot(log(atf6_genes_exp_atac$sibAvg+1), 
     log(atf6_genes_exp_atac$atf6Avg+1))


lines(lowess(log(atf6_genes_exp_atac$sibAvg+1), 
             log(atf6_genes_exp_atac$atf6Avg+1)), col="red")
  
subset(atac_pos_exp, sibAvg<10 & atf6Avg>20)

subset(atf6_genes_exp_atac, sibAvg<10 & atf6Avg>20)
subset(atf6_genes_exp_atac, sibAvg < atf6Avg)
subset(atac_pos_exp, sibAvg < atf6Avg)


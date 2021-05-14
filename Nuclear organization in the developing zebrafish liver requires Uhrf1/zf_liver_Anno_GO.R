setwd("~/Documents/NYUAD/Zf_ATAC/Zf_Liver_NB551229")
library(genomation)
library(ChIPseeker)
library(GenomicRanges)
gene.parts <- readTranscriptFeatures('~/Documents/iGenome/zf/GRCz10_refseq.bed',
                                     up.flank = 2000, down.flank = 2000)
zf_lvr_atac <- readPeakFile('./macs2/ATAC_merge_LV_peaks.broadPeak', header=F) 
### add 'chr' to chr name
seqlevels(zf_lvr_atac) <-  paste('chr', seqlevels(zf_lvr_atac), sep='')
zf_lvr_atac
anno <- annotateWithGeneParts(zf_lvr_atac, gene.parts, intersect.chr=T)
plotTargetAnnotation(anno, precedence=T)

### 
library("TxDb.Drerio.UCSC.danRer10.refGene")
library("org.Dr.eg.db")
txdb <- TxDb.Drerio.UCSC.danRer10.refGene
promoter <- getPromoters(TxDb=txdb, upstream=2000, downstream=2000)
Ann <- annotatePeak(zf_lvr_atac,TxDb = txdb, 
                       annoDb = "org.Dr.eg.db")
plotAnnoPie(Ann)
plotAnnoPie(Ann, main= 'zf liver Annotation')
### GO analysis
library(clusterProfiler)
library(dplyr)
tab <- as.data.frame(Ann)
promoter_entriz <- tab %>% dplyr::filter(annotation=='Promoter (1-2kb)' | 
                  annotation=='Promoter (<=1kb)' |
                  annotation=='Promoter (2-3kb)') %>% dplyr::select(geneId)

ego_pro <- clusterProfiler::enrichGO(gene     = unique(promoter_entriz$geneId),
                                    OrgDb         = org.Dr.eg.db,
                                    ont           = "BP",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = 0.05,
                                    qvalueCutoff  = 0.05, 
                                    readable      = TRUE)
dotplot(ego_pro, showCategory=20)

Intergenic_entriz <- tab %>% 
  dplyr::filter(annotation=='Distal Intergenic') %>%
  dplyr::select(geneId)

ego_intergenic <- clusterProfiler::enrichGO(gene     = unique(Intergenic_entriz$geneId),
                                     OrgDb         = org.Dr.eg.db,
                                     ont           = "BP",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     qvalueCutoff  = 0.05, 
                                     readable      = TRUE)
dotplot(ego_intergenic, showCategory=20)
### combine promoter and distal intergenic
entriz_pr_inter <- c(promoter_entriz$geneId, Intergenic_entriz$geneId)
ego_pr_inter <- clusterProfiler::enrichGO(gene     = entriz_pr_inter,
                                            OrgDb         = org.Dr.eg.db,
                                            ont           = "BP",
                                            pAdjustMethod = "BH",
                                            pvalueCutoff  = 0.05,
                                            qvalueCutoff  = 0.05, 
                                            readable      = TRUE)
dotplot(ego_pr_inter, showCategory=20)
tab <- as.data.frame(ego_pr_inter)
write.table(tab, 'zf_lvr_atac_GO_summary.txt', quote=F, 
            col.names = T, sep='\t')
# enlarge the promoter pool from annotatePeak()
# entriz_tss <- seq2gene(zf_lvr_atac, tssRegion = c(-1000, 1000), 
# flankDistance = 1000, TxDb=txdb)
#==============================================================================================
### TE annotation
### add 'chr' into atac peak file
# cat ATAC_merge_LV_peaks.broadPeak | sed 's/^/chr/' > ATAC_merge_LV_peaks.bed
### Annotation with *.bed files
# bedtools annotate -i ~/Documents/iGenome/zf/GRCz10_RM_TE.bed -files ./macs2/ATAC_merge_LV_peaks.bed > zf_lvr_ATAC-TEanno.txt
###  getting result into R
Anno_all <- read.table('./zf_lvr_ATAC-TEanno.txt', row.names = NULL)
### filter  non-annotate
Anno_df <- Anno_all %>% filter(V7 != 0)
head(Anno_df)
### Ann_df show % proportion in each state
colnames(Anno_df) <- c('chr', 'start', 'end','name',
                       'class', 'family','r1')
### convert Annotation df to Granges
Anno_gr <- makeGRangesFromDataFrame(Anno_df, keep.extra.columns = T)

### calculate proportion of LTR, LINE and SINE in each state
### bp-wised 
LTR_gr <- Anno_gr[Anno_gr$class=='LTR']
LINE_gr <- Anno_gr[Anno_gr$class=='LINE']
SINE_gr <- Anno_gr[Anno_gr$class=='SINE']
DNA_gr <- Anno_gr[Anno_gr$class=='DNA']
SR_gr <- Anno_gr[Anno_gr$class=='Simple_repeat']

ltr <- sum(round(width(LTR_gr)*LTR_gr$r1)) # LTR
line <- sum(round(width(LINE_gr)*LINE_gr$r1)) # LINE
sine <- sum(round(width(SINE_gr)*SINE_gr$r1)) #  SINE
dna <- sum(round(width(DNA_gr)*DNA_gr$r1)) # DNA
simple_repeats <- sum(round(width(SR_gr)*SR_gr$r1))
### idx of each state i
whole <- sum(width(zf_lvr_atac))


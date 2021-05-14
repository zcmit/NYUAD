setwd("~/Documents/NYUAD/Zf_ATAC/Zf_Liver_NB551229/Bams/")

## load the library
library(ATACseqQC)
## input the bamFile from the ATACseqQC package 
bamfile <- grep('.bam', list.files('./'), value=T)
bamfile.labels <- gsub("_mdup.withrg.csorted.cleaned.aligned.bam",
                       "", basename(bamfile))

estimateLibComplexity(readsDupFreq(bamfile[1]))
fragSize <- fragSizeDist(bamfile[3], bamfile.labels[3])

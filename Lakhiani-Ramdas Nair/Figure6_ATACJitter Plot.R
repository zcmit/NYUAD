library(ggplot2)
library(ggpubr)
library(dplyr)
## filter non-expressed genes out, which are rowSum(fpkm) < 20
merged<-mergedATAC_RNAseq_readcount[-which(duplicated(mergedATAC_RNAseq_readcount$ensembl_id)), ]
library(tidyverse)
merged <- data.frame(merged, row.names = 1)
merged_filtered<- merged[rowSums(merged)>=20, ]
### normalization
merged_upregulated <- merged_filtered + 1
merged_upregulated_norm <- log(merged_upregulated)

merged_UP_open <- merged_upregulated_norm[which(rownames(merged_upregulated_norm) %in% ATAC_open$zebrafish), ]
merged_UP_closed <- merged_upregulated_norm[which(rownames(merged_upregulated_norm) %in% ATAC_closed$zebrafish), ]
### set bivalent gene expression


dat_open <- data.frame(
  name=c(rep("Sibling", nrow(merged_UP_open )), rep("nAtf6",nrow(merged_UP_open )), rep("UT",nrow(merged_UP_open )), 
         rep("iAs",nrow(merged_UP_open )), rep("UT2", nrow(merged_UP_open )), rep("EtOH", nrow(merged_UP_open ))),
        value=c(merged_UP_open $Atf6_sibslivermean, merged_UP_open $Atf6_Atf6livermean, merged_UP_open $iAs_UTliversmean, merged_UP_open $iAs_liversmean,
          merged_UP_open $EtOH_UTliversmean, merged_UP_open$EtOH_liversmean))
dat_open$name <- factor(dat_open$name, levels = c("Sibling", "nAtf6", "UT", "iAs", "UT2", "EtOH"))





dat_closed <- data.frame(
  name=c(rep("Sibling", nrow(merged_UP_closed )), rep("nAtf6",nrow(merged_UP_closed )), rep("UT",nrow(merged_UP_closed )), 
         rep("iAs",nrow(merged_UP_closed)), rep("UT2", nrow(merged_UP_closed )), rep("EtOH", nrow(merged_UP_closed ))),
  value=c(merged_UP_closed $Atf6_sibslivermean, merged_UP_closed $Atf6_Atf6livermean, merged_UP_closed $iAs_UTliversmean, merged_UP_closed $iAs_liversmean,
          merged_UP_closed $EtOH_UTliversmean, merged_UP_closed $EtOH_liversmean))


dat<- data.frame(
  name=c(rep("Untreated-open", nrow(merged_UP_open )), rep("iAs-open", nrow(merged_UP_open )), rep("Untreated- closed",nrow(merged_UP_closed )), rep("iAs-closed",nrow(merged_UP_closed ))),
  value=c(merged_UP_open $iAs_UTliversmean,  merged_UP_open $iAs_liversmean,merged_UP_closed $iAs_UTliversmean, merged_UP_closed $iAs_liversmean))
dat$name <- factor(dat$name, levels = c("Untreated-open", "iAs-open", "Untreated- closed", "iAs-closed"))

dat_closed$name <- factor(dat_closed$name, levels = c("Sibling", "nAtf6", "UT", "iAs", "UT2", "EtOH"))
comparisons = list( merged_UP_open [1:2], merged_UP_closed[1:2], merged_UP_open[3:4] , merged_UP_closed[3:4])
p <- ggplot(dat, aes(x=name, y=value)) + geom_violin(trim = FALSE)
p + geom_jitter(shape=16, position=position_jitter(0.2), color='black')+
  stat_summary(fun=median, color="red", geom='point')+
  theme_classic()+
stat_compare_means(comparisons = comparisons, label = "p.signif", size=15)
wilcox.test(merged_UP_open$EtOH_UTliversmean,merged_UP_open$EtOH_liversmean)


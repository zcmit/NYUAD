### Full Capstone R Script ###

#### Biomart for gene conversions (zebrafish to human genes) ####
# I need to find human genes for my dnmt1 s904 lvrs dataset that includes only zebrafish genes
# The dataset I received already has the ensembl IDs

# Install biomaRt for R
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

library('biomaRt')
ensembl <- useMart('ensembl')  # put the name of the dataset to use. I'm using the ensembl database
listDatasets(ensembl) # to see which datasets are present in ensembl
# Number 60 drerio_gene_ensembl       Zebrafish genes (GRCz11)
# Number 85 hsapiens_gene_ensembl     Human genes (GRCh38.p13)
mart_hg <- useDataset('hsapiens_gene_ensembl', useMart('ensembl'))
mart_zf <- useDataset('drerio_gene_ensembl', useMart('ensembl'))

# list attributes from each dataset
listAttributes(mart_hg)
listAttributes(mart_zf)

# human: 'hgnc_symbol', zebrafish: 'zfin_id_symbol'

listFilters(mart_hg)  # check which filters are available
listFilters(mart_zf)  # check which filters are available

# Convert IDs
# getBM(attributes=Output, filters=Input_filter, values= Gene or Gene list, mart=TARGET_MART) 
# example: convert 'pcna' to get its ensembl id
dnmt1_data <- dnmt1_lvr_zf$ensembl_id
dnmt1_data_IDtoNAMES <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),filters='ensembl_gene_id',
                              values=dnmt1_data, mart=mart_hg)



#### Import, clean and sort datasets: dnmt1 s904 ####
# Import dataset
dnmt1_lvr_zf <- read.table("DNMT1_lvr_zf.txt", sep="\t", header = T)

# Rename columns so that they make sense to me
colnames(dnmt1_lvr_zf) <- c("ensembl_id","dnmt1_baseMean","dnmt1_log2FC", "dnmt1_ifcSE", "dnmt1_stat", "dnmt1_pvalue",
                            "dnmt1_padj", "dnmt1_weight", "dnmt1_Mut1", "dnmt1_Mut2", "dnmt1_Mut3", "dnmt1_Sib1",
                            "dnmt1_Sib2", "dnmt1_Sib3", "zfin_id_symbol")

# Subset for significant, up-regulated, and down-regulated then export
dnmt1_sig <- subset(dnmt1_lvr_zf, dnmt1_padj<0.05)
UP_dnmt1_sig <- subset(dnmt1_sig, dnmt1_log2FC>1.5)
write.table(UP_dnmt1_sig, file="UP_dnmt1_sig.txt", quote = F, sep="\t", row.names=F)
DOWN_dnmt1_sig <- subset(dnmt1_sig, dnmt1_log2FC< (-1.5))
write.table(DOWN_dnmt1_sig, file="DOWN_dnmt1_sig.txt", quote = F, sep="\t", row.names=F)

# Create ensembl ID lists and export
UP_dnmt1_ensembl <- UP_dnmt1_sig[,1]
write.table(UP_dnmt1_ensembl, file="UP_dnmt1_ensembl.txt", quote = F, sep="\t", row.names=F)
DOWN_dnmt1_ensembl <- DOWN_dnmt1_sig[,1]
write.table(DOWN_dnmt1_ensembl, file="DOWN_dnmt1_ensembl.txt", quote = F, sep="\t", row.names=F)

#### Import, clean and sort datasets: uhrf1 hi272 ####
uhrf1_hi272_lvr <- read.csv("MTvsSib_120uhrf1_lvr.csv", sep=",", header = T)

# Rename columns so that they make sense to me
colnames(uhrf1_hi272_lvr) <- c("ensembl_id","hi272_baseMean","hi272_log2FC", "hi272_ifcSE",
                               "hi272_stat", "hi272_pvalue", "hi272_padj", "hi272_weight", "hi272_Mut1",
                               "hi272_Mut2", "hi272_Sib1", "hi272_Sib2", "hi272_Sib4", "zfin_id_symbol")

# Subset significant, up-regulated and down-regulated then export
hi272_sig <- subset(uhrf1_hi272_lvr, hi272_padj<0.05)
UP_hi272_sig <- subset(hi272_sig, hi272_log2FC>1.5)
write.table(UP_hi272_sig, file="UP_hi272_sig.txt", quote = F, sep="\t", row.names=F)
DOWN_hi272_sig <- subset(hi272_sig, hi272_log2FC< (-1.5))
write.table(DOWN_hi272_sig, file="DOWN_hi272_sig.txt", quote = F, sep="\t", row.names=F)

# Create only ensembl id lists and export
UP_hi272_ensembl <- UP_hi272_sig[,1]
write.table(UP_hi272_ensembl, file="UP_hi272_ensembl.txt", quote = F, sep="\t", row.names=F)
DOWN_hi272_ensembl <- DOWN_hi272_sig[,1]
write.table(DOWN_hi272_ensembl, file="DOWN_hi272_ensembl.txt", quote = F, sep="\t", row.names=F)


#### Import, clean and sort datasets: hUHRF1-overexpression (hUHRF1-high) ####
hUHRF1_high <- read.csv("uhrf1High_vs_sib_hgEnsembl.csv", sep=",", header = T)

# Rename columns so that they make sense to me
colnames(hUHRF1_high) <- c("ensembl_id","high_baseMean", "high_log2FC", "high_ifcSE","high_stat","high_pvalue", "high_padj",
                           "high_weight","hUHRF1_high7_counts", "hUHRF1_high8_counts", "hUHRF1_high9_counts",
                           "hUHRF1_sib7_counts", "hUHRF1_sib8_counts", "hUHRF1_sib9_counts", "zfin_id_symbol",
                           "Gene.stable.ID.1","HGNC.symbol")

# Subset significant, up-regulated, and down-regulated then export
hUHRF1_high_sig <- subset(hUHRF1_high, high_padj<0.05)
UP_hUHRF1_high_sig <- subset(hUHRF1_high_sig, high_log2FC>1.5)
write.table(UP_hUHRF1_high_sig, file="UP_hUHRF1_high_sig.txt", quote = F, sep="\t", row.names=F)
DOWN_hUHRF1_high_sig <- subset(hUHRF1_high_sig, high_log2FC< (-1.5))
write.table(DOWN_hUHRF1_high_sig, file="DOWN_hUHRF1_high_sig.txt", quote = F, sep="\t", row.names=F)

# Create  ensembl ID lists and export
UP_hUHRF1_high_ensembl <- UP_hUHRF1_high_sig[,1]
write.table(UP_hUHRF1_high_ensembl, file="UP_hUHRF1_high_ensembl.txt", quote = F, sep="\t", row.names=F)
DOWN_hUHRF1_high_ensembl <- DOWN_hUHRF1_high_sig[,1]
write.table(DOWN_hUHRF1_high_ensembl, file="DOWN_hUHRF1_high_ensembl.txt", quote = F, sep="\t", row.names=F)


#### Create master dataset for later use ####
# Create and export master dataset with dnmt1 s904, uhrf1 hi272, and hUHRF1-high datasets
Merge_dnmt1s904_hUHRF1high <- merge(dnmt1_lvr_zf, hUHRF1_high, by="ensembl_id")
Master_dnmt1s904_hUHRF1high_uhrf1hi272 <- merge(Merge_dnmt1s904_hUHRF1high, uhrf1_hi272_lvr, by="ensembl_id")
# Subset for only columns that I care about
Master_dnmt1s904_hUHRF1high_uhrf1hi272 <- subset(Master_dnmt1s904_hUHRF1high_uhrf1hi272, select=-c(dnmt1_ifcSE,
 dnmt1_stat, dnmt1_pvalue,dnmt1_weight,high_ifcSE,high_stat,high_pvalue,high_weight,hi272_ifcSE,hi272_stat,
 hi272_pvalue,hi272_weight))
# Export master dataset
write.table(Master_dnmt1s904_hUHRF1high_uhrf1hi272,
  file="Master_dnmt1s904_hUHRF1high_uhrf1hi272.txt", quote = F, sep="\t", row.names=F)



### Figure 5 plots: dnmt1 s904 ###
#### dnmt1 volcano plot ####
plot(dnmt1_lvr_zf$dnmt1_log2FC, -log(dnmt1_lvr_zf$dnmt1_padj), pch=20, xlim=c(-8,10),
     xlab="log2FC", ylab="-log(pdaj)", main="Dnmt1 s904 Volcano Plot", col="grey")
with(subset(dnmt1_lvr_zf, dnmt1_padj< 0.05), 
     points(dnmt1_log2FC, -log(dnmt1_padj), pch=20, col="black"))
with(subset(dnmt1_lvr_zf, dnmt1_padj< 0.05 & dnmt1_log2FC >1.5), 
     points(dnmt1_log2FC, -log(dnmt1_padj), pch=20, col="purple"))
with(subset(dnmt1_lvr_zf, dnmt1_padj< 0.05 & dnmt1_log2FC <(-1.5)), 
     points(dnmt1_log2FC, -log(dnmt1_padj), pch=20, col="orange"))
legend(-11, 390, c("Non-significant gene", "Significant gene", "UP-regulated gene","DOWN-regulated gene"),
       cex=0.8, bty = "n",col=c("grey","black","purple","orange"),
       horiz = F, pch=c(16), y.intersp=0.2)

#### dnmt1 MA plot ####
plot(log(dnmt1_lvr_zf$dnmt1_baseMean), dnmt1_lvr_zf$dnmt1_log2FC, pch=20,
     xlab="log10 Base Mean", ylab="log2 Fold Change", main="Dnmt1 s904 MA Plot", col="grey")
with(subset(dnmt1_lvr_zf, dnmt1_padj< 0.05), 
     points(log(dnmt1_baseMean), dnmt1_log2FC, pch=20, col="black"))
with(subset(dnmt1_lvr_zf, dnmt1_padj< 0.05 & dnmt1_log2FC >1.5), 
     points(log(dnmt1_baseMean), dnmt1_log2FC, pch=20, col="purple"))
with(subset(dnmt1_lvr_zf, dnmt1_padj< 0.05 & dnmt1_log2FC <(-1.5)), 
     points(log(dnmt1_baseMean), dnmt1_log2FC, pch=20, col="orange"))

#### dnmt1 MA plot with WT siblings mean ####
# First, add columns for means of read counts for muts and sibs
dnmt1_lvr_zf$dnmt1_mut_mean <- (dnmt1_lvr_zf$dnmt1_Mut1 + dnmt1_lvr_zf$dnmt1_Mut2 + dnmt1_lvr_zf$dnmt1_Mut3)/3
dnmt1_lvr_zf$dnmt1_sib_mean <- (dnmt1_lvr_zf$dnmt1_Sib1 + dnmt1_lvr_zf$dnmt1_Sib2 + dnmt1_lvr_zf$dnmt1_Sib3)/3
# Then, add column for sibs mean plus 1
# I do this because some genes' mean reads=0, and the log of it is infinity, which isn't
# good for plotting the data. So I add 1 to all values, then take the log of that
dnmt1_lvr_zf$dnmt1_sib_mean_plus1 <- dnmt1_lvr_zf$dnmt1_sib_mean + 1
# Finally, plot with sib reads mean
plot(log(dnmt1_lvr_zf$dnmt1_sib_mean_plus1), dnmt1_lvr_zf$dnmt1_log2FC, pch=20,
     xaxt = "n",
     xlab="log10 WT Sib Mean", ylab="log2 Fold Change", main="Dnmt1 s904 WT Sibs Mean MA Plot", col="grey")
with(subset(dnmt1_lvr_zf, dnmt1_padj< 0.05), 
     points(log(dnmt1_sib_mean_plus1), dnmt1_log2FC, pch=20, col="black"))
with(subset(dnmt1_lvr_zf, dnmt1_padj< 0.05 & dnmt1_log2FC >1.5), 
     points(log(dnmt1_sib_mean_plus1), dnmt1_log2FC, pch=20, col="purple"))
with(subset(dnmt1_lvr_zf, dnmt1_padj< 0.05 & dnmt1_log2FC <(-1.5)), 
     points(log(dnmt1_sib_mean_plus1), dnmt1_log2FC, pch=20, col="orange"))


### Figure 5 plots: hUHRF1-high ###
#### hUHRF1-high volcano plot ####
plot(hUHRF1_high$high_log2FC, -log(hUHRF1_high$high_padj), pch=20, xlim=c(-10,10),
     xaxt="n",
     xlab="log2FC", ylab="-log(pdaj)", main="hUHRF1-high Volcano plot", col="grey")
with(subset(hUHRF1_high, high_padj< 0.05), 
     points(high_log2FC, -log(high_padj), pch=20, col="black"))
with(subset(hUHRF1_high, high_padj< 0.05 & high_log2FC >1.5), 
     points(high_log2FC, -log(high_padj), pch=20, col="purple"))
with(subset(hUHRF1_high, high_padj< 0.05 & high_log2FC <(-1.5)), 
     points(high_log2FC, -log(high_padj), pch=20, col="orange"))
legend(-11, 78, c("Non-significant gene", "Significant gene", "UP-regulated gene","DOWN-regulated gene"),
       cex=0.8, bty = "n",col=c("grey","black","purple","orange"),
       horiz = F, pch=c(16), y.intersp=0.2)


#### hUHRF1-high MA plot ####
plot(log(hUHRF1_high$high_baseMean) ,hUHRF1_high$high_log2FC, pch=20,
     xaxt="n", yaxt="n",
     xlab="log10 Base Mean", ylab="log2 Fold Change", main="hUHRF1-high MA plot", col="grey")
with(subset(hUHRF1_high, high_padj< 0.05), 
     points(log(high_baseMean), high_log2FC, pch=20, col="black"))
with(subset(hUHRF1_high, high_padj< 0.05 & high_log2FC >1.5), 
     points(log(high_baseMean), high_log2FC, pch=20, col="purple"))
with(subset(hUHRF1_high, high_padj< 0.05 & high_log2FC <(-1.5)), 
     points(log(high_baseMean), high_log2FC, pch=20, col="orange"))

#### dnmt1 MA plot with WT siblings mean ####
# First, add columns for means of read counts for highs and sibs
hUHRF1_high$hUHRF1_high_mean <- (hUHRF1_high$hUHRF1_high7_counts + hUHRF1_high$hUHRF1_high8_counts + hUHRF1_high$hUHRF1_high9_counts)/3
hUHRF1_high$hUHRF1_sib_mean <- (hUHRF1_high$hUHRF1_sib7_counts + hUHRF1_high$hUHRF1_sib8_counts + hUHRF1_high$hUHRF1_sib9_counts)/3
# Then, add column for sibs mean plus 1
# I do this because some genes' mean reads=0, and the log of it is infinity, which isn't
# good for plotting the data. So I add 1 to all values, then take the log of that
hUHRF1_high$hUHRF1_sib_mean_plus1 <- hUHRF1_high$hUHRF1_sib_mean + 1
# Finally, plot with sib reads mean
plot(log(hUHRF1_high$hUHRF1_sib_mean_plus1), hUHRF1_high$high_log2FC, pch=20,
     xaxt="n", yaxt="n",
     xlab="log10 WT Sib Mean", ylab="log2 Fold Change", main="hUHRF1 High WT Sibs Mean MA Plot", col="grey")
with(subset(hUHRF1_high, high_padj< 0.05), 
     points(log(hUHRF1_sib_mean_plus1), high_log2FC, pch=20, col="black"))
with(subset(hUHRF1_high, high_padj< 0.05 & high_log2FC >1.5), 
     points(log(hUHRF1_sib_mean_plus1), high_log2FC, pch=20, col="purple"))
with(subset(hUHRF1_high, high_padj< 0.05 & high_log2FC <(-1.5)), 
     points(log(hUHRF1_sib_mean_plus1), high_log2FC, pch=20, col="orange"))



#### Figure 6: dnmt1 s904 & hUHRF1-high mean reads ####
# Plot mean reads for WT sibs from dnmt1 s904 and hUHRF1-high using Master dataset
plot(dnmt1_high$dnmt1_sib_mean ,dnmt1_high$hUHRF1_sib_mean, pch=20, xlim=c(0,100000), ylim=c(0,100000),
     xlab="Dnmt1 WT Siblings Mean Reads", ylab="hUHRF1-high WT Siblings Mean Reads", main="Dnmt1 vs. hUHRF1-high Wildtype Siblings Mean Reads", col="black")
# add line
abline(lm(dnmt1_high$dnmt1_sib_mean ~ dnmt1_high$hUHRF1_sib_mean), col="blue")
# add correlation coefficient
r <- cor(dnmt1_high$dnmt1_sib_mean, dnmt1_high$hUHRF1_sib_mean, method = "pearson")
r
text(90000, 90000, labels='r = 0.98')

#### Figure 6 Venn diagram: dnmt1 s904 & hUHRF1-high ####


library(venneuler)
require(venneuler)
# I use Venny online to find the numbers I used in the function below.
v <- venneuler(c(dnmt1_UP=1564, dnmt1_DOWN=677, hUHRF1_high_UP=802, hUHRF1_high_DOWN=587, "dnmt1_UP&hUHRF1_high_UP"=189, "dnmt1_UP&hUHRF1_high_DOWN"=13,
                 "hUHRF1_high_DOWN&dnmt1_DOWN"=77,"dnmt1_DOWN&hUHRF1_high_UP"=57))
plot(v, col=c("#00664a", "#00f2b1","red4", "pink"))
# Add numbers manually on Adobe Illustrator

#### Figure 6 Venn diagram: dnmt1 s904 & uhrf1 hi272 ####
library(venneuler)
require(venneuler)
# I use Venny online to find the numbers I used in the function below.
v <- venneuler(c(dnmt1_UP=1564, dnmt1_DOWN=677, uhrf11_hi272_UP=1146, uhrf11_hi272_DOWN=1081, "dnmt1_UP&uhrf11_hi272_UP"=632, "dnmt1_UP&uhrf11_hi272_DOWN"=4,
                 "uhrf11_hi272_DOWN&dnmt1_DOWN"=186,"dnmt1_DOWN&uhrf11_hi272_UP"=3))
# I make "dnmt1_UP&uhrf11_hi272_DOWN"=4 just to make the overlap clearer, but the actual number is 2 #
plot(v, col=c("#00664a", "#00f2b1","yellow4", "yellow"))
# Add numbers manually on Adobe Illustrator

#### Figure 6: significant genes heatmap ####
#install.packages("pheatmap")
library("pheatmap")
# Create a new dataframe with only significant human gene names and values from each line (dnmt1 muts, highs, uhrf1 muts)
heatmap <- subset(Master_dnmt1s904_hUHRF1high_uhrf1hi272, (dnmt1_padj<0.05 | high_padj<0.05 | hi272_padj<0.05))
Sig_heatmap <- heatmap[,c(5,6,7,8,9,10,18,19,20,21,22,23,26,33,34,35,36,37)]
# reorder columns to have all sibs next to each other, then the muts and highs
Sig_heatmap <- Sig_heatmap[,c(13,10,11,12,4,5,6,16,17,18,7,8,9,1,2,3,14,15)]

NAMES <- Sig_heatmap$HGNC.symbol
print(NAMES)
rownames(Sig_heatmap) <- NAMES

Sig_heatmap_matrix <- subset(Sig_heatmap, select = -HGNC.symbol)

common_matrix <- as.matrix(Sig_heatmap_matrix)

# z score normalization
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
# Apply function and transpose matrix
common_matrix_norm <- t(apply(common_matrix, 1, cal_z_score))
pheatmap(common_matrix_norm, cluster_cols=F, cutree_rows = 5, show_rownames=F)
# Heatmap function takes a while to load


#### Figure 7 log2FC plot: dnmt1 s904 & hUHRF1-high ####
# Use Master dataset. Assign to new column "Significance". Based on this, colors will be different
# 4 options: only dnmt1 significant, only high significant, both dnmt1 and high significant, none are significant
dnmt1_high$Significance <- (with(dnmt1_high, ifelse(dnmt1_padj<0.05 & high_padj>0.05, Significance <- "Only dnmt1",
                                                    ifelse(dnmt1_padj>0.05 & high_padj<0.05, Significance <- "Only hUHRF1-high",
                                                           ifelse(dnmt1_padj<0.05 & high_padj<0.05, Significance <- "Both",
                                                                  Significance <- "None")))))

# Subset for significantly expressed genes only
dnmt1_high_sig_only <- subset(dnmt1_high, dnmt1_padj<0.05 | high_padj<0.05)

# Plot only significant genes (overlapping 3 plots)
library(ggplot2)
plot <- ggplot(dnmt1_high_sig_only, aes(x=dnmt1_log2FC, y=high_log2FC, color=Significance)) +
  geom_point(size=2.1, shape=1)  +
  theme_minimal() +
  coord_fixed() +  
  xlim(-10, 10) +
  ylim (-6,6) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
  scale_color_manual(values=c("#999999", "#cc79a7", "#56b4e9")) + # colors for significance
  geom_vline(xintercept = 1.5, color="#f0e442", linetype="solid", size=0.6) + # lines for up/down regulation
  geom_vline(xintercept = -1.5, color="#f0e442", linetype="solid", size=0.6) +
  geom_hline(yintercept = -1.5, color="#f0e442", linetype="solid", size=0.6) +
  geom_hline(yintercept = 1.5, color="#f0e442", linetype="solid", size=0.6) +
  geom_text(x= 1.5, y=6, label="log2FC=1.5", color="black") +
  geom_text(x= -1.5, y=6, label="log2FC=-1.5", color="black") +
  geom_text(x= -9.6, y=-1.5, label="log2FC=-1.5", color="black") +
  geom_text(x= -9.6, y=1.5, label="log2FC=-1.5", color="black")
plot

# Both dnmt1 s904 and hUHRF1-highs significant
both <- subset(dnmt1_high_sig_only, dnmt1_padj<0.05 & high_padj<0.05)
plot1 <- ggplot(both, aes(x=dnmt1_log2FC, y=high_log2FC, color=Significance)) +
  geom_point(size=2.1, shape=1)  +
  theme_minimal() +
  coord_fixed() +  
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_vline(xintercept = 1.5, color="#f0e442", linetype="solid", size=0.6) +
  geom_vline(xintercept = -1.5, color="#f0e442", linetype="solid", size=0.6) +
  geom_hline(yintercept = -1.5, color="#f0e442", linetype="solid", size=0.6) +
  geom_hline(yintercept = 1.5, color="#f0e442", linetype="solid", size=0.6) +
  scale_color_manual(values=c("#999999")) +
  geom_text(x= 1.5, y=9.5, label="log2FC=1.5", color="black") +
  geom_text(x= -1.5, y=9.5, label="log2FC=-1.5", color="black") +
  geom_text(x= -5.4, y=-1.5, label="log2FC=-1.5", color="black") +
  geom_text(x= -5.4, y=1.5, label="log2FC=-1.5", color="black")
plot1

# Only dnmt1 s904 significant
only_dnmt1 <- subset(dnmt1_high_sig_only, dnmt1_padj<0.05 &! high_padj<0.05)
plot2 <- ggplot(only_dnmt1, aes(x=dnmt1_log2FC, y=high_log2FC, color=Significance)) +
  geom_point(size=2.1, shape=1)  +
  theme_minimal() +
  coord_fixed() +  
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_vline(xintercept = 1.5, color="#f0e442", linetype="solid", size=0.6) +
  geom_vline(xintercept = -1.5, color="#f0e442", linetype="solid", size=0.6) +
  geom_hline(yintercept = -1.5, color="#f0e442", linetype="solid", size=0.6) +
  geom_hline(yintercept = 1.5, color="#f0e442", linetype="solid", size=0.6) +
  scale_color_manual(values=c("#cc79a7")) +
  geom_text(x= 1.5, y=6.5, label="log2FC=1.5", color="black") +
  geom_text(x= -1.5, y=6.5, label="log2FC=-1.5", color="black") +
  geom_text(x= -5, y=-1.5, label="log2FC=-1.5", color="black") +
  geom_text(x= -5, y=1.5, label="log2FC=-1.5", color="black")
plot2

# Only hUHRF1-highs significant
only_highs <- subset(dnmt1_high_sig_only, high_padj<0.05 &! dnmt1_padj<0.05)
plot3 <- ggplot(only_highs, aes(x=dnmt1_log2FC, y=high_log2FC, color=Significance)) +
  geom_point(size=2.1, shape=1)  +
  theme_minimal() +
  coord_fixed() + 
  xlim (-2.5,2.5) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_vline(xintercept = 1.5, color="#f0e442", linetype="solid", size=0.6) +
  geom_vline(xintercept = -1.5, color="#f0e442", linetype="solid", size=0.6) +
  geom_hline(yintercept = -1.5, color="#f0e442", linetype="solid", size=0.6) +
  geom_hline(yintercept = 1.5, color="#f0e442", linetype="solid", size=0.6) +
  scale_color_manual(values=c("#56b4e9")) +
  geom_text(x= 1.5, y=6.5, label="log2FC=1.5", color="black") +
  geom_text(x= -1.5, y=6.5, label="log2FC=-1.5", color="black") +
  geom_text(x=-1.6, y=-1.5, label="log2FC=-1.5", color="black") +
  geom_text(x=-1.6, y=1.5, label="log2FC=-1.5", color="black")
plot3

#### Figure 7 log2FC plot: dnmt1 s904 & uhrf1 hi272 ####
# Use Master dataset. Assign to new column "Significance". Based on this, colors will be different
dnmt1_hi272 <- Master_dnmt1_hi272_highs[,c(-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13)]
#4 options: only dnmt1 significant, only hi272 significant, both dnmt1 and hi272 significant, none are significant
dnmt1_hi272$Significance <- (with(dnmt1_hi272, ifelse(dnmt1_padj<0.05 & hi272_padj>0.05, Significance <- "dnmt1_only",
                                                      ifelse(dnmt1_padj>0.05 & hi272_padj<0.05, Significance <- "uhrf1_hi272_only",
                                                             ifelse(dnmt1_padj<0.05 & hi272_padj<0.05, Significance <- "both",
                                                                    Significance <- "none")))))

# Subset for significant;y expressed genes only
dnmt1_hi272_sig_only <- subset(dnmt1_hi272, dnmt1_padj<0.05 | hi272_padj<0.05)

# Plot only significant genes (overlapping 3 plots)
library(ggplot2)
plot <- ggplot(dnmt1_hi272_sig_only, aes(x=dnmt1_log2FC, y=hi272_log2FC, color=Significance)) +
  geom_point(size=2.1, shape=1)  +
  theme_minimal() +
  coord_fixed() + 
  ylim (-6,6) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  scale_color_manual(values=c("#999999", "#cc79a7", "#009e73")) +
  geom_vline(xintercept = 1.5, color="#f0e442", linetype="solid", size=0.6) +
  geom_vline(xintercept = -1.5, color="#f0e442", linetype="solid", size=0.6) +
  geom_hline(yintercept = -1.5, color="#f0e442", linetype="solid", size=0.6) +
  geom_hline(yintercept = 1.5, color="#f0e442", linetype="solid", size=0.6) +
  geom_text(x= 1.5, y=6, label="log2FC=1.5", color="black") +
  geom_text(x= -1.5, y=6, label="log2FC=-1.5", color="black") +
  geom_text(x= -5.3, y=-1.5, label="log2FC=-1.5", color="black") +
  geom_text(x= -5.3, y=1.5, label="log2FC=-1.5", color="black")
plot

# Both dnmt1 s904 and uhrf1 hi272 significant
library(ggplot2)
both <- subset(dnmt1_hi272_sig_only, dnmt1_padj<0.05 & hi272_padj<0.05)
plot1 <- ggplot(both, aes(x=dnmt1_log2FC, y=hi272_log2FC, color=Significance)) +
  geom_point(size=2.1, shape=1)  +
  theme_minimal() +
  coord_fixed() +  
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_vline(xintercept = 1.5, color="#f0e442", linetype="solid", size=0.6) +
  geom_vline(xintercept = -1.5, color="#f0e442", linetype="solid", size=0.6) +
  geom_hline(yintercept = -1.5, color="#f0e442", linetype="solid", size=0.6) +
  geom_hline(yintercept = 1.5, color="#f0e442", linetype="solid", size=0.6) +
  scale_color_manual(values=c("#999999"))
# +  geom_text(x= 1.5, y=6, label="log2FC=1.5", color="black") +
#geom_text(x= -1.5, y=6, label="log2FC=-1.5", color="black") +
#geom_text(x= -5.3, y=-1.5, label="log2FC=-1.5", color="black") +
#geom_text(x= -5.3, y=1.5, label="log2FC=-1.5", color="black")
plot1

# Only dnmt1 s904 significant
only_dnmt1 <- subset(dnmt1_hi272_sig_only, dnmt1_padj<0.05 &! hi272_padj<0.05)
plot2 <- ggplot(only_dnmt1, aes(x=dnmt1_log2FC, y=hi272_log2FC, color=Significance)) +
  geom_point(size=2.1, shape=1)  +
  theme_minimal() +
  coord_fixed() +  
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_vline(xintercept = 1.5, color="#f0e442", linetype="solid", size=0.6) +
  geom_vline(xintercept = -1.5, color="#f0e442", linetype="solid", size=0.6) +
  geom_hline(yintercept = -1.5, color="#f0e442", linetype="solid", size=0.6) +
  geom_hline(yintercept = 1.5, color="#f0e442", linetype="solid", size=0.6) +
  scale_color_manual(values=c("#cc79a7"))
#+ geom_text(x= 1.5, y=2.8, label="log2FC=1.5", color="black") +
# geom_text(x= -1.5, y=2.8, label="log2FC=-1.5", color="black") +
#  geom_text(x= -5.3, y=-1.5, label="log2FC=-1.5", color="black") +
# geom_text(x= -5.3, y=1.5, label="log2FC=-1.5", color="black")
plot2

# Only uhrf1 hi272 significant
only_hi272 <- subset(dnmt1_hi272_sig_only, hi272_padj<0.05 &! dnmt1_padj<0.05)
plot3 <- ggplot(only_hi272, aes(x=dnmt1_log2FC, y=hi272_log2FC, color=Significance)) +
  geom_point(size=2.1, shape=1)  +
  theme_minimal() +
  coord_fixed() + 
  xlim (-2.5,2.5) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_vline(xintercept = 1.5, color="#f0e442", linetype="solid", size=0.6) +
  geom_vline(xintercept = -1.5, color="#f0e442", linetype="solid", size=0.6) +
  geom_hline(yintercept = -1.5, color="#f0e442", linetype="solid", size=0.6) +
  geom_hline(yintercept = 1.5, color="#f0e442", linetype="solid", size=0.6) +
  scale_color_manual(values=c("#009e73"))
plot3

#### Figure 9: top 100 Up genes in dnmt1 s904 ####
# create list for only top 100 up sig genes to use in GO analysis
# Arrange genes from lowest padj to highest padj
library(tidyverse)
library(tidyr)
UP_dnmt1_sig_ordered <- arrange(UP_dnmt1_sig, dnmt1_padj) #arranges by padj from lowest to highest
UP_100_dnmt1_sig <- slice(UP_dnmt1_sig_ordered, c(1:100)) #take top 100 most significant (lowest padj)
# export file for use in creating the correlation matrix code
write.table(UP_100_dnmt1_sig, file="UP_100_dnmt1_sig.txt", quote = F, sep="\t", row.names=F)

# Create and export ensembl ID lists
UP_100_dnmt1_sig_IDs <- UP_100_dnmt1_sig[,1]
write.table(UP_100_dnmt1_sig_IDs, file="UP_100_dnmt1_sig_IDs.txt", quote = F, sep="\t", row.names=F)


### Figure 10: oncogenes dnmt1 s904 ###
#### find oncogenes in dnmt1 s904 ####
# To see if there are oncogenes in my dnmt1 dataset
# First, I need human gene names -> I use the merged master dataset I have,
# which is: Master_dnmt1s904_hUHRF1high_uhrf1hi272

# Now, import the oncogenes list that I compiled from researching online.
oncogenes_list <- read.csv("Oncogenes list.csv", sep=",", header = T)
oncogenes_list <- oncogenes_list[,c(1,2)]

# Check to see if human genes in dnmt1 dataset have any matches with the oncogenes dataset,
# and save position to a new variable: oncogenes_dnmt1
oncogenes_dnmt1_position <- which(Master_dnmt1s904_hUHRF1high_uhrf1hi272$HGNC.symbol%in%oncogenes_list$Oncogenes)
oncogenes_dnmt1_position
# 138 matches
# Look for the position and extract the data
oncogenes_dnmt1 <- Master_dnmt1s904_hUHRF1high_uhrf1hi272[oncogenes_dnmt1_position, ]
oncogenes_in_dnmt1 <- oncogenes_dnmt1[,26]
write.table(oncogenes_in_dnmt1, file= "oncogenes_in_dnmt1.txt", quote = F, sep="\t", row.names = F)
# remove replicates and blanks manually
# left with 97 unique oncogene matches

# Subset to see how many significant, UP and DOWN
oncogenes_dnmt1_sig <- subset(oncogenes_dnmt1, dnmt1_padj<0.05)
# 42 significant, 55 not significant
UP_oncogenes_dnmt1 <- subset(oncogenes_dnmt1_sig, dnmt1_log2FC>1.5)
# 8 significantly UP for log2FC>1.5
# 12 significantly UP for log2FC>1.2
# 21 significantly UP for log2FC>1

DOWN_oncogenes_dnmt1 <- subset(oncogenes_dnmt1_sig, dnmt1_log2FC<(-1.5))
# 1 significantly DOWN for log2FC<-1.5
# 2 significantly DOWN for log2FC<-1.2
# 3 significantly DOWN for log2FC<-1
#### dnmt1 s904 oncogenes volcano plot ####
plot(oncogenes_dnmt1$dnmt1_log2FC, -log(oncogenes_dnmt1$dnmt1_padj), pch=20,
     xlab="log2FC", ylab="-log(pdaj)", main="Oncogenes in dnmt1 s904 Volcano Plot", col="grey")
with(subset(oncogenes_dnmt1, dnmt1_padj< 0.05), 
     points(dnmt1_log2FC, -log(dnmt1_padj), pch=20, col="black"))
with(subset(oncogenes_dnmt1, dnmt1_padj< 0.05 & dnmt1_log2FC >1.5), 
     points(dnmt1_log2FC, -log(dnmt1_padj), pch=20, col="purple"))
with(subset(oncogenes_dnmt1, dnmt1_padj< 0.05 & dnmt1_log2FC <(-1.5)), 
     points(dnmt1_log2FC, -log(dnmt1_padj), pch=20, col="orange"))

### Figure 10: oncogenes hUHRF1-high ###
#### find oncogenes in hUHRF1-high ####
# I've already imported my compiled oncogenes list above
# to use with the dnmt1 dataset, so no need to do it again

# Check to see if human genes in highs dataset have any matches with the curated oncogenes dataset,
# and save position to a new variable: oncogenes_highs
oncogenes_high_position <- which(hUHRF1_high$HGNC.symbol%in%oncogenes_list$Oncogenes)
oncogenes_high_position
# 126 matches
# Look for the position and extract the data
oncogenes_high <- hUHRF1_high[oncogenes_high_position, ] # this line extracts data from the row numbers
#specified in "oncogenes_high_position" and assigns them to the new dataset "oncogenes_high"
oncogenes_in_high <- oncogenes_high[,17] # get human gene names for the oncogenes
write.table(oncogenes_in_high, file= "oncogenes_in_high.txt", quote = F, sep="\t", row.names = F)
# remove replicates and blanks manually
# left with 97 unique oncogene matches

# Subset to see how many significant, UP and DOWN
oncogenes_high_sig <- subset(oncogenes_high, high_padj<0.05)
# 15 significant, 42 not significant
UP_oncogenes_high <- subset(oncogenes_high_sig, high_log2FC>1.5)
# 5 significantly UP for log2FC>1.5
# 6 significantly UP for log2FC>1.2
# 6 significantly UP for log2FC>1

DOWN_oncogenes_high <- subset(oncogenes_high_sig, high_log2FC<(-1.5))
# 0 significantly DOWN for log2FC<-1.5
# 3 significantly DOWN for log2FC<-1.2
# 5 significantly DOWN for log2FC<-1
#### hUHRF1-high oncogenes volcano plot ####
plot(oncogenes_high$high_log2FC, -log(oncogenes_high$high_padj), pch=20,
     xlab="log2FC", ylab="-log(pdaj)", main="Oncogenes in hUHRF1-high Volcano Plot", col="grey")
with(subset(oncogenes_high, high_padj< 0.05), 
     points(high_log2FC, -log(high_padj), pch=20, col="black"))
with(subset(oncogenes_high, high_padj< 0.05 & high_log2FC >1.5), 
     points(high_log2FC, -log(high_padj), pch=20, col="purple"))
with(subset(oncogenes_high, high_padj< 0.05 & high_log2FC <(-1.5)), 
     points(high_log2FC, -log(high_padj), pch=20, col="orange"))


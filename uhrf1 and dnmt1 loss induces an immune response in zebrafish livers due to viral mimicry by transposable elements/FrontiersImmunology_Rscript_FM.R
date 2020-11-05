# find all genome CpGs in Zebrafish (all the chromosome including mitocondrial)
########################################################################################################################################
library(BSgenome.Drerio.UCSC.danRer10)  
chrs <- names(Drerio)[1:26]
cgs <- lapply(chrs, function(x) start(matchPattern("CG", Drerio[[x]])))
cpgr <- do.call(c, lapply(1:26, function(x) GRanges(names(Drerio)[x], IRanges(cgs[[x]], width = 2))))
#########################################################################################################

# importing RRBS CpGs
#########################################################################################################
CpGs_All <- read.delim("./Input_Data/CpGs_All.txt")
CpGs_All_LowMeth <- subset(CpGs_All, CpGs_All$X5dpf_WT_sib_RRBS<20)
CpGs_All_HighMeth <- subset(CpGs_All, CpGs_All$X5dpf_WT_sib_RRBS>80)

# converting in GRanges
GR_CpGs_All <- as(CpGs_All, "GRanges")
GR_CpGs_All_LowMeth <- as(CpGs_All_LowMeth, "GRanges")
GR_CpGs_All_HighMeth <- as(CpGs_All_HighMeth, "GRanges")
#########################################################################################################

# import and manipulate Repetitive Elements locations
# (TE nomenclature in objects stand for RE)
#########################################################################################################
All_TE <- read.delim("./Input_Data/GRCz10_RM_TE.txt")
head(All_TE)

# reodering columns
All_TE = All_TE[c("genoName", "genoStart", "genoEnd", "repName",
                  "repClass", "repFamily", "swScore")]
head(All_TE)

# renaming columns for conversion in Genomic Ranges
colnames(All_TE) <- c("chr", "start", "end", "repName",
                      "repClass", "repFamily", "swScore")
#########################################################################################################

# subsetting data.frame
#########################################################################################################
All_TE_LTR <- subset(All_TE, repClass==("LTR"))
All_TE_SINE <- subset(All_TE, repClass==("SINE"))
All_TE_LINE <- subset(All_TE, repClass==("LINE"))
All_TE_DNA <- subset(All_TE, repClass==("DNA"))

# converting in Genomic Ranges
#########################################################################################################
All_TE_LTR_GR <- as(All_TE_LTR, "GRanges")
All_TE_SINE_GR <- as(All_TE_SINE, "GRanges")
All_TE_LINE_GR <- as(All_TE_LINE, "GRanges")
All_TE_DNA_GR <- as(All_TE_DNA, "GRanges")
#########################################################################################################

# import and manipulate Repetitive Elements locations and expression
# in RNA-seq of uhrf1 mutant whole embryos
#########################################################################################################
uhrf1_mut_whole_TE <- read.delim("./Input_Data/whole_uhrf1_TE2_gR.txt")
head(uhrf1_mut_whole_TE)

# reodering columns
uhrf1_mut_whole_TE = uhrf1_mut_whole_TE[c("genoName", "genoStart", "genoEnd", "repName", "baseMean",
                                          "log2FoldChange", "padj", "repClass", "repFamily", "swScore")]
head(uhrf1_mut_whole_TE)

# renaming columns for conversion in Genomic Ranges
colnames(uhrf1_mut_whole_TE) <- c("chr", "start", "end", "repName", "baseMean",
                                  "log2FoldChange", "padj", "repClass", "repFamily", "swScore")

# subsetting data.frame
#########################################################################################################
uhrf1_mut_whole_TE_LTR <- subset(uhrf1_mut_whole_TE, repClass==("LTR"))
uhrf1_mut_whole_TE_SINE <- subset(uhrf1_mut_whole_TE, repClass==("SINE"))
uhrf1_mut_whole_TE_LINE <- subset(uhrf1_mut_whole_TE, repClass==("LINE"))
uhrf1_mut_whole_TE_DNA <- subset(uhrf1_mut_whole_TE, repClass==("DNA"))

uhrf1_mut_whole_Transposons_UP_signif <- subset(uhrf1_mut_whole_Transposons, padj< 0.05 & log2FoldChange >0)
uhrf1_mut_whole_Transposons_DOWN_signif <- subset(uhrf1_mut_whole_Transposons, padj< 0.05 & log2FoldChange <0)

# converting in Genomic Ranges
#########################################################################################################
uhrf1_mut_whole_TE_GR <- as(uhrf1_mut_whole_TE, "GRanges")
uhrf1_mut_whole_TE_LTR_GR <- as(uhrf1_mut_whole_TE_LTR, "GRanges")
uhrf1_mut_whole_TE_SINE_GR <- as(uhrf1_mut_whole_TE_SINE, "GRanges")
uhrf1_mut_whole_TE_LINE_GR <- as(uhrf1_mut_whole_TE_LINE, "GRanges")
uhrf1_mut_whole_TE_DNA_GR <- as(uhrf1_mut_whole_TE_DNA, "GRanges")

uhrf1_mut_whole_Transposons_UP_signif_GR <- as(uhrf1_mut_whole_Transposons_UP_signif, "GRanges")
uhrf1_mut_whole_Transposons_DOWN_signif_GR <- as(uhrf1_mut_whole_Transposons_DOWN_signif, "GRanges")

# subsetting Genomic Ranges by overlap between regions (query, subject) the query drive the GR given as result 
#########################################################################################################
CpGs_uhrf1_mut_whole_Transposons_UP_signif_GR <- subsetByOverlaps(GR_CpGs_All, uhrf1_mut_whole_Transposons_UP_signif_GR)
CpGs_uhrf1_mut_whole_Transposons_DOWN_signif_GR <- subsetByOverlaps(GR_CpGs_All, uhrf1_mut_whole_Transposons_DOWN_signif_GR)
########################################################################


# CpGs Annotation by genomation package - Fig 3A
########################################################################################################################################
library(genomation)
bed.file = "./Input_Data/refseq_danRer10.gz"
gene.parts = readTranscriptFeatures(bed.file)

# Whole Genome CpGs
annot_cpgr_All_Genome = annotateWithGeneParts(cpgr, gene.parts, strand=TRUE, intersect.chr=TRUE)
plotTargetAnnotation(annot_cpgr_All_Genome)
# take the numbers (with promoter > exon > intron precedence)
annot_cpgr_All_Genome

# All RRBS CpGs
annot_GR_CpGs_All = annotateWithGeneParts(GR_CpGs_All, gene.parts, strand=TRUE, intersect.chr=TRUE)
plotTargetAnnotation(annot_GR_CpGs_All)
# take the numbers (with promoter > exon > intron precedence)
annot_GR_CpGs_All

# RRBS CpGs - unmethylated (methylation level <20%)
annot_GR_CpGs_All_LowMeth = annotateWithGeneParts(GR_CpGs_All_LowMeth, gene.parts, strand=TRUE, intersect.chr=TRUE)
plotTargetAnnotation(annot_GR_CpGs_All_LowMeth)
# take the numbers (with promoter > exon > intron precedence)
annot_GR_CpGs_All_LowMeth

# RRBS CpGs - methylated (methylation level >80%)
annot_GR_CpGs_All_HighMeth = annotateWithGeneParts(GR_CpGs_All_HighMeth, gene.parts, strand=TRUE, intersect.chr=TRUE)
plotTargetAnnotation(annot_GR_CpGs_All_HighMeth)
# take the numbers (with promoter > exon > intron precedence)
annot_GR_CpGs_All_HighMeth
#########################################################################################################


# All RRBS CpGs density Plot - Fig. 3B
#########################################################################################################
df_perc.meth_All_CpGs <- data.frame(GR_CpGs_All$X5dpf_WT_sib_RRBS, GR_CpGs_All$X5dpf_uhrf1_mut_RRBS)
colnames(df_perc.meth_All_CpGs) <- c("WT_sib","uhrf1_mut")

library("reshape2")
df_perc.meth_All_CpGs_melt <- melt(df_perc.meth_All_CpGs)
head(df_perc.meth_All_CpGs_melt)
tail(df_perc.meth_All_CpGs_melt)

# create a table with median of each variable for intersection of the line in the geom_vline
########################################################################
library(plyr)
perc.meth_median <- ddply(df_perc.meth_All_CpGs_melt,
                          "variable", summarise, grp.median=median(value))
head(perc.meth_median)

library("ggplot2")
ggplot(df_perc.meth_All_CpGs_melt, aes(value)) +
  geom_density(aes(color = variable, fill=variable), size= 1, alpha=0.5) +
  scale_x_continuous(breaks=seq(0,100,25)) +
  labs(y= "Density of CpGs per each sample", x= "Percentage of CpGs methylation values") +
  theme_classic() + scale_fill_manual(values = c("black", "red")) +
  scale_color_manual(values=c("black", "red")) +
  geom_vline(data=perc.meth_median, aes(xintercept=grp.median, color=variable), linetype="dashed", size= 1)
#########################################################################################################


# All RRBS CpGs LowMeth density Plot - Fig. 3C.a
#########################################################################################################
df_perc.meth_All_CpGs_LowMeth <- data.frame(GR_CpGs_All_LowMeth$X5dpf_WT_sib_RRBS, GR_CpGs_All_LowMeth$X5dpf_uhrf1_mut_RRBS)
colnames(df_perc.meth_All_CpGs_LowMeth) <- c("WT_sib","uhrf1_mut")

library("reshape2")
df_perc.meth_All_CpGs_LowMeth_melt <- melt(df_perc.meth_All_CpGs_LowMeth)
head(df_perc.meth_All_CpGs_LowMeth_melt)
tail(df_perc.meth_All_CpGs_LowMeth_melt)

# create a table with median of each variable for intersection of the line in the geom_vline
########################################################################
library(plyr)
perc.meth_median <- ddply(df_perc.meth_All_CpGs_LowMeth_melt,
                          "variable", summarise, grp.median=median(value))
head(perc.meth_median)

library("ggplot2")
ggplot(df_perc.meth_All_CpGs_LowMeth_melt, aes(value)) +
  geom_density(aes(color = variable, fill=variable), size= 1, alpha=0.5) +
  scale_x_continuous(breaks=seq(0,100,25)) +
  labs(y= "Density of CpGs per each sample", x= "Percentage of CpGs methylation values") +
  theme_classic() + scale_fill_manual(values = c("black", "red")) +
  scale_color_manual(values=c("black", "red")) +
  geom_vline(data=perc.meth_median, aes(xintercept=grp.median, color=variable), linetype= c("twodash", "dotted"), size= 1)
#########################################################################################################

# All CpGs HighMeth density Plot - Fig. 3C.b
#########################################################################################################
df_perc.meth_All_CpGs_HighMeth <- data.frame(GR_CpGs_All_HighMeth$X5dpf_WT_sib_RRBS, GR_CpGs_All_HighMeth$X5dpf_uhrf1_mut_RRBS)
colnames(df_perc.meth_All_CpGs_HighMeth) <- c("WT_sib","uhrf1_mut")

library("reshape2")
df_perc.meth_All_CpGs_HighMeth_melt <- melt(df_perc.meth_All_CpGs_HighMeth)
head(df_perc.meth_All_CpGs_HighMeth_melt)
tail(df_perc.meth_All_CpGs_HighMeth_melt)

# create a table with median of each variable for intersection of the line in the geom_vline
########################################################################
library(plyr)
perc.meth_median <- ddply(df_perc.meth_All_CpGs_HighMeth_melt,
                          "variable", summarise, grp.median=median(value))
head(perc.meth_median)

library("ggplot2")
ggplot(df_perc.meth_All_CpGs_HighMeth_melt, aes(value)) +
  geom_density(aes(color = variable, fill=variable), size= 1, alpha=0.5) +
  scale_x_continuous(breaks=seq(0,100,25)) +
  labs(y= "Density of CpGs per each sample", x= "Percentage of CpGs methylation values") +
  theme_classic() + scale_fill_manual(values = c("black", "red")) +
  scale_color_manual(values=c("black", "red")) +
  geom_vline(data=perc.meth_median, aes(xintercept=grp.median, color=variable), linetype= c("twodash", "dotted"), size= 1)
#########################################################################################################

# All Repetitive Elements in uhrf1 mutant RNA-seq Genomation annotation - Fig 3D
# (RE but object contain TE nomenclature)
#########################################################################################################
# all Repetitive Elements
library(genomation)
annot_CpGs_uhrf1_mut_whole_TE_GR = annotateWithGeneParts(CpGs_uhrf1_mut_whole_TE_GR, gene.parts, strand=TRUE, intersect.chr=TRUE)
plotTargetAnnotation(annot_CpGs_uhrf1_mut_whole_TE_GR)
# take the numbers (with promoter > exon > intron precedence)
annot_CpGs_uhrf1_mut_whole_TE_GR
###################################################################

# TE subclasses
###################################################################
# LTR
annot_CpGs_uhrf1_mut_whole_TE_LTR_GR = annotateWithGeneParts(CpGs_uhrf1_mut_whole_TE_LTR_GR, gene.parts, strand=TRUE, intersect.chr=TRUE)
plotTargetAnnotation(annot_CpGs_uhrf1_mut_whole_TE_LTR_GR)
# take the numbers (with promoter > exon > intron precedence)
annot_CpGs_uhrf1_mut_whole_TE_LTR_GR

# SINE
annot_CpGs_uhrf1_mut_whole_TE_SINE_GR = annotateWithGeneParts(CpGs_uhrf1_mut_whole_TE_SINE_GR, gene.parts, strand=TRUE, intersect.chr=TRUE)
plotTargetAnnotation(annot_CpGs_uhrf1_mut_whole_TE_SINE_GR)
# take the numbers (with promoter > exon > intron precedence)
annot_CpGs_uhrf1_mut_whole_TE_SINE_GR

# LINE
annot_CpGs_uhrf1_mut_whole_TE_LINE_GR = annotateWithGeneParts(CpGs_uhrf1_mut_whole_TE_LINE_GR, gene.parts, strand=TRUE, intersect.chr=TRUE)
plotTargetAnnotation(annot_CpGs_uhrf1_mut_whole_TE_LINE_GR)
# take the numbers (with promoter > exon > intron precedence)
annot_CpGs_uhrf1_mut_whole_TE_LINE_GR

# DNA
annot_CpGs_uhrf1_mut_whole_TE_DNA_GR = annotateWithGeneParts(CpGs_uhrf1_mut_whole_TE_DNA_GR, gene.parts, strand=TRUE, intersect.chr=TRUE)
plotTargetAnnotation(annot_CpGs_uhrf1_mut_whole_TE_DNA_GR)
# take the numbers (with promoter > exon > intron precedence)
annot_CpGs_uhrf1_mut_whole_TE_DNA_GR
#########################################################################################################


# All CpGs in WTsibs of TE diff. expressed in RNA-seq uhrf1_mut_whole - Fig 3E
# Transposons UP-regulated and DOWN-regulated - Significant
#########################################################################################################
sq <- seq(max(length(CpGs_uhrf1_mut_whole_Transposons_UP_signif_GR),
              length(CpGs_uhrf1_mut_whole_Transposons_DOWN_signif_GR)))
# this is forcing to keep the same lenght for all columns and put NAs
df_perc.meth_uhrf1_mut_whole_Transposons_UP_DOWN_signif_WT_sib <- data.frame(CpGs_uhrf1_mut_whole_Transposons_UP_signif_GR$X5dpf_WT_sib_RRBS[sq],
                                                                             CpGs_uhrf1_mut_whole_Transposons_DOWN_signif_GR$X5dpf_WT_sib_RRBS[sq])
# rename columns
colnames(df_perc.meth_uhrf1_mut_whole_Transposons_UP_DOWN_signif_WT_sib) <- c("WT-sib_Transposons-UP",
                                                                              "WT-sib_Transposons-DOWN")
library("reshape2")
df_perc.meth_uhrf1_mut_whole_Transposons_UP_DOWN_signif_WT_sib_melt <- melt(df_perc.meth_uhrf1_mut_whole_Transposons_UP_DOWN_signif_WT_sib)
head(df_perc.meth_uhrf1_mut_whole_Transposons_UP_DOWN_signif_WT_sib_melt)
tail(df_perc.meth_uhrf1_mut_whole_Transposons_UP_DOWN_signif_WT_sib_melt)

library("ggplot2")
ggplot(df_perc.meth_uhrf1_mut_whole_Transposons_UP_DOWN_signif_WT_sib_melt, aes(x = variable, y = value, fill=variable)) + 
  geom_boxplot(width=0.5, notch = FALSE) + scale_y_continuous(breaks=seq(0,100,25)) + 
  labs(x= "Transposons significant UP-regulated and DOWN-regulated", y = "Percentage of CpGs methylation values") + 
  theme_classic() + scale_fill_manual(values = c("violetred4", "palegreen4"))
#########################################################################################################


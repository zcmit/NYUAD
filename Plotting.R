### Histogram of four controls (figure1 in methylome comparison paper) ### 
dat <- read.delim('4_mice.txt', header=F, sep='\t')
library(ggplot2)
dat1 <- data.frame(val=c(as.numeric(dat[1,2:ncol(dat)]),
                         as.numeric(dat[2,2:ncol(dat)]),
                         as.numeric(dat[3,2:ncol(dat)]),
                         as.numeric(dat[4,2:ncol(dat)])), 
                   var=c(rep('Con-1',20), rep('Con-2',20),
                         rep('Con-3',20), rep('Con-4',20)),
                   group=rep(seq(2.5,97.5,by=5),2))
qplot(x=group, y=val, fill=var,
      data=dat1, geom="bar", stat="identity",
      position="dodge")

g <- ggplot(dat1, aes(x=group, y=val, fill=var)) + geom_bar(stat="identity",
                                                            position=position_dodge()) +
  scale_fill_manual(values=c("#CC6666", "#9999CC", "#66CC99","orange4"))

g+theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        legend.title=element_blank()) + scale_x_continuous(breaks=seq(0,100, by=10)) + scale_y_continuous(breaks=seq(0,60, by=10))
# g + theme_economist()+scale_colour_economist()


### Volcano Plot ###
### volcano plot
plot(dat$log2FoldChange, -log(dat$padj), pch=20, xlim=c(-6,6), 
     ylim = c(0,80), xlab="log2FC", ylab="-log(padj)", main="Title")
### color significant genes
with(subset(dat, padj<.05), points(log2FoldChange, -log(padj), pch=20, col="red"))
# with(subset(dat_vol, padj<.05 & log2FoldChange >2), points(log2FoldChange, -log(padj), pch=22, col="blue"))
# set cutoff
cutoff=sort(dat$padj)[20] #the 20th smallest value of res$padj
sig_genes <- which(dat$padj <= cutoff)
sig_genes <- which(dat$padj<0.05)
# add gene name next to points
text(x=dat$log2FoldChange[sig_genes], y=-log(dat$pvalue[sig_genes]),label=dat_cytokines$zfin_id_symbol[sig_genes], cex=0.8)

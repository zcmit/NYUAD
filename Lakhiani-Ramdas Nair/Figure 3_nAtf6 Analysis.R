source('Functions.R')

# Find significantly DE genes in nAtf6 ------------------------------------

natf6_rnaseq = rnaseq_analysis('data/RNAseq files/nAtf6_summary.txt')
atf6_targets=subset(Atf6,Atf6$ensembl_id%in%Atf6PutTargets_all_1498$Zebrafish_Ensembl_ID)
atf6_targets=unique(atf6_targets)
atf6_targets_sig=subset(atf6_targets,atf6_targets$padj<=0.05)

natf6_sig = unique(na.omit(natf6_rnaseq[natf6_rnaseq$pValue <= 0.05,]))

# Figure 3 - Make violin plot ---------------------------------------------

natf6_colour = '#eb3a13'
natf6_violin = plot_violin('All DEGs in nAtf6',natf6_colour,natf6_sig,2.5)

png(filename = 'figures/Fig 3/Fig3_violin.png',850,1000)
natf6_violin
dev.off()

# Figure 3 - Venn ---------------------------------------------------------

fig3_venn = make_3way_eulerr(natf6_sig$ensembl_id,targets,upr,'figures/Fig 3/fig3_venn',
                             c(natf6_colour,targets_colour,upr_colour),7)

# Figure 3 - UpSet plot ---------------------------------------------------

library(UpSetR)
upset_natf6 = unique(c(targets,upr))
n = length(upset_natf6)
upset_natf6 = data.frame(gene = upset_natf6,DE=character(n))
upset_natf6$Atf6_Targets = as.numeric(upset_natf6$gene %in% targets)
# upset_natf6$Atf6_Mouse_KO_Targets = as.numeric(upset_natf6$gene %in% unique(rutkowski2008_mouseconsensus_homologs$zebrafish))
upset_natf6$UPR_Genes = as.numeric(upset_natf6$gene %in% upr)

upset_natf6$DE = c('notDE','DE')[as.numeric(upset_natf6$gene %in% natf6_sig$ensembl_id) + 1]
upset_natf6$DE = as.factor(upset_natf6$DE)

Myfunc <- function(row){
  data <- row["DE"]=='DE'
}

sets = colnames(upset_natf6)[3:4]

#labels need to be added manually as a result of this loop and then moved above bars on illustrator
#for each condiition, bottom number is DE and top number is total
for (i in 0:1){
  for (j in 0:1){
    print(paste(sets[1],i,sets[2],j))
    print(summary(upset_natf6$Atf6_Targets==i & upset_natf6$UPR_Genes==j))        
    print(summary(upset_natf6$Atf6_Targets==i & upset_natf6$UPR_Genes==j & upset_natf6$DE =='DE'))        
   }
}

pdf('figures/Fig 3/fig3_upset.pdf',width = 4,height = 4.8,onefile=FALSE)

upset(upset_natf6,
      mb.ratio = c(0.8,0.2),sets = rev(sets),keep.order = T,
      order.by = 'degree',decreasing = F,
      query.legend = 'top',
      sets.bar.color = c(upr_colour,targets_colour),
      queries = list(
        list(query = Myfunc,active=T,color=natf6_colour,query.name='Differentially Expressed in nAtf6')
      ),
      mainbar.y.max = 1450,
      text.scale = c(1.5,1,1,1,1.2,1),
      sets.x.label = 'Gene List Size',
      set_size.show = TRUE,
      set_size.scale_max = 2000,
      show.numbers = TRUE)
grid.text("354/1444",x = 0.62, y=0.92, gp=gpar(fontsize=8))
grid.text("81/244",x = 0.698, y=0.36, gp=gpar(fontsize=8))
grid.text("25/54",x = 0.76, y=0.235, gp=gpar(fontsize=8))

dev.off()

# Figure 3 - GO Enrichment ------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
OrgDb = org.Hs.eg.db

primary_targets = intersect(targets,natf6_sig$ensembl_id)
secondary_targets = targets[!targets %in% natf6_sig$ensembl_id]
tertiary_targets= targets[!(natf6_sig$ensembl_id%in%targets)]

primary_targets_human = unique(targets_homologs$human[targets_homologs$zebrafish %in% primary_targets])
secondary_targets_human = unique(targets_homologs$human[targets_homologs$zebrafish %in% secondary_targets])
tertiary_targets_human = unique(targets_homologs$human[targets_homologs$zebrafish %in% tertiary_targets])
figS5 = list(primarytargets = primary_targets_human, secondarytargets = secondary_targets_human, tertiarytargets = tertiary_targets_human)
figS5_go = compareCluster(fig3,fun = "enrichGO",keyType = "ENSEMBL", OrgDb = OrgDb,ont = "BP",readable = TRUE)

create_revigo_table(figS5_go,'data/REVIGO/input/fig3_goidswithp.txt')
figS5_go = drop_from_revigo(figS5_go,"data/REVIGO/output/REVIGO_figS5.csv")
scale = 1.5
wrapsize = 38

figS5_plottable = for_plotting_comparecluster(figS5,figS5_go,6,wrapsize)

plot_figS5 = ggplot(figS5_plottable, aes(x = Cluster, y = Description)) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_color_continuous(name=str_wrap('Adjusted p-value', width = wrapsize),low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
  ylab(NULL) +
  xlab('')+
  scale_x_discrete(labels= c('Atf6-Sufficient Targets\n(DE in nAtf6)','Atf6-Insufficient Targets','Indirect Atf6 Targets\n(DE in nAtf6; \nnot putative target)'))+
  theme(legend.key.height = unit(30*scale,'points'),
        axis.text.y = element_text(size = 23*scale,angle=0, hjust=1),
        axis.text.x = element_text(angle = 45,hjust=1,size=18*scale),
        text = element_text(size=18*scale),
        axis.title.x = element_text(size=22*scale))+
  scale_size(name='Ratio of Genes\nin Gene List',range=c(6*scale, 20*scale),breaks=c(0.01,0.02,0.03))+
  theme(strip.text.y = element_blank())+
  geom_text(aes(label=Count,size=0.00028,fontface='bold'),colour='white',show.legend = F)

png('figures/figS5_GO.png',width = 800*scale,height = 800*scale)
plot_figS5
dev.off()

save(fig3_go,file = 'data/RData/fig3_go.RData')


# Figure S5 - GO Enrichment of Atf6 and Atf6 remaining targets ------------------------------------------------

library(clusterProfiler)
library(org.Dr.eg.db)
library(AnnotationDbi)
library(ggplot2)
OrgDb = org.Dr.eg.db
natf6_DEG=unique(natf6_sig$ensembl_id)
atf6targets<-unique(atf6_targets$ensembl_id)

Atf6_filter<-filter(natf6_sig,!(natf6_sig$ensembl_id%in%targets))
Atf6_human<- subset.data.frame(targets_homologs, targets_homologs$zebrafish %in% Atf6_filter$ensembl_id)
figS5 = list(Atf6DEGs = natf6_DEG, Atf6Targets=atf6targets, RemainingAtf6 = Atf6_filter)
figS5_go = compareCluster(figS5,fun = "enrichGO",keyType = "ENSEMBL", OrgDb = OrgDb,ont = "BP",readable = TRUE)

create_revigo_table(figS5_go,'figS5_goidswithp.txt')
figS5_go = drop_from_revigo(figS5_go,"REVIGO_figS5.csv")
scale = 1.5
wrapsize = 38

figS5_plottable = for_plotting_comparecluster(figS5,figS5_go,6,wrapsize)

plot_figS5 = ggplot(figS5_plottable, aes(x = Cluster, y = Description)) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_color_continuous(name=str_wrap('Adjusted p-value', width = wrapsize),low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
  ylab(NULL) +
  xlab('')+
  scale_x_discrete(labels= c('DE in Atf6 ','DE Atf6 targets', 'Remaining DEGS'))+
  theme(legend.key.height = unit(30*scale,'points'),
        axis.text.y = element_text(size = 23*scale,angle=0, hjust=1),
        axis.text.x = element_text(angle = 45,hjust=1,size=18*scale),
        text = element_text(size=18*scale),
        axis.title.x = element_text(size=22*scale))+
  scale_size(name='Ratio of Genes\nin Gene List',range=c(6*scale, 20*scale),breaks=c(0.04,0.06,0.08))+
  theme(strip.text.y = element_blank())+
  geom_text(aes(label=Count,size=0.00001,fontface='bold'),colour='white',show.legend = F)

png('figures/figS5_GO.png',width = 800*scale,height = 800*scale)
plot_figS5
dev.off()

save(fig3_go,file = 'data/RData/fig3_go.RData')

















if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Dr.eg.db")
library(org.Dr.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
library(ggplot2)
OrgDb = org.Dr.eg.db
atf6_sig=subset(natf6_rnaseq,natf6_rnaseq$padj<0.05) 
Fig3e<-unique(atf6_sig$ensembl_id)


gene.df <- bitr(Fig3e, fromType = "ENSEMBL", #Upgenes = whatever ENSEMBL list you'd like
                toType = c("ENTREZID", "SYMBOL"),
                OrgDb = OrgDb) #check organism
# Group GO: BP: biological process
ggo_bp <- clusterProfiler::groupGO(gene     = gene.df$ENTREZID,
                                   OrgDb    = OrgDb,
                                   ont      = "BP",
                                   level    = 3,
                                   readable = TRUE)
head(summary(ggo_bp)[,-5])
barplot(ggo_bp, drop=TRUE, showCategory=12,self_contained = TRUE)
# GO over-representation test for BP
ego_bp <- clusterProfiler::enrichGO(gene          = gene.df$ENTREZID,
                                    OrgDb         = OrgDb,
                                    ont           = "BP",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = 0.05,
                                    qvalueCutoff  = 0.05,
                                    readable      = TRUE)
barplot(ego_bp, showCategory=20)
clusterProfiler::dotplot(ego_bp, showCategory=15)
nottargets<-atf6_sig[!targets]
library(dplyr)
Atf6_filter<-filter(atf6_sig,!(atf6_sig$ensembl_id%in%targets))
fig3f=unique(Atf6_filter$ensembl_id)


  gene.df <- bitr(fig3f, fromType = "ENSEMBL", #Upgenes = whatever ENSEMBL list you'd like
                  toType = c("ENTREZID", "SYMBOL"),
                  OrgDb = OrgDb) #check organism
# Group GO: BP: biological process
ggo_bp <- clusterProfiler::groupGO(gene     = gene.df$ENTREZID,
                                   OrgDb    = OrgDb,
                                   ont      = "BP",
                                   level    = 3,
                                   readable = TRUE)
head(summary(ggo_bp)[,-5])
barplot(ggo_bp, drop=TRUE, showCategory=12,self_contained = TRUE)
# GO over-representation test for BP
ego_bp <- clusterProfiler::enrichGO(gene          = gene.df$ENTREZID,
                                    OrgDb         = OrgDb,
                                    ont           = "BP",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = 0.05,
                                    qvalueCutoff  = 0.05,
                                    readable      = TRUE)
barplot(ego_bp, showCategory=20)

natf6_shoulders<-intersect(shoulders2013_targets,natf6_sig$ensembl_id)

shoulders_notnatf6 = shoulders2013_targets[!shoulders2013_targets %in% natf6_sig$ensembl_id] 

targets_notatf6=subset(Atf6PutTargets_all_1498,Atf6PutTargets_all_1498$Zebrafish_Ensembl_ID%in% shoulders_notnatf6)
targets_etoh=subset(targets_notatf6,targets_notatf6$Zebrafish_Ensembl_ID %in% etoh_sig$ensembl_id)
targets_iAs=subset(targets_notatf6,targets_notatf6$Zebrafish_Ensembl_ID %in% ias_sig$ensembl_id)


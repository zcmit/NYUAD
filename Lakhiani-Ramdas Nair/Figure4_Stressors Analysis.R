source('code/Functions.R')

# Find significantly DE genes in EtOH -------------------------------------

etoh_rnaseq = rnaseq_analysis('data/RNAseq files/EtOH_summary.txt')
etoh_sig = unique(na.omit(etoh_rnaseq[etoh_rnaseq$pValue <= 0.05,]))
etoh_upregulated= subset (etoh_sig, log2FoldChange>=0.02)
etoh_downregulated= subset (etoh_sig, log2FoldChange<=0.02)

# Find significantly DE genes in iAs --------------------------------------

ias_rnaseq = rnaseq_analysis('data/RNAseq files/iAs_summary.txt')
ias_sig = unique(na.omit(ias_rnaseq[ias_rnaseq$pValue <= 0.05,]))
iAs_upregulated= subset (ias_sig, log2FoldChange>=0.02)
iAs_downregulated= subset (ias_sig, log2FoldChange<=0.02)
# Figure 4 - Make violin plots --------------------------------------------

etoh_colour = '#ffd21f'
ias_colour = '#A546FD'

etoh_violin = plot_violin('All DEGs in EtOH',etoh_colour,etoh_sig,2.5)
ias_violin = plot_violin('All DEGs in iAs',ias_colour,ias_sig,2.5)

png(filename = 'figures/Fig 4/Fig4_etoh_violin.png',850,1000)
etoh_violin
dev.off()

png(filename = 'figures/Fig 4/Fig4_ias_violin.png',850,1000)
ias_violin
dev.off()

# Figure 4 - Venns --------------------------------------------------------

fig4_etoh_venn = make_3way_eulerr(etoh_sig$ensembl_id,targets,upr,'figures/Fig 4/fig4_etoh_venn',
                             c(etoh_colour,targets_colour,upr_colour),7)

fig4_ias_venn = make_3way_eulerr(ias_sig$ensembl_id,targets,upr,'figures/Fig 4/fig4_ias_venn',
                                  c(ias_colour,targets_colour,upr_colour),7)

# Figure 4 - UpSet plot - EtOH --------------------------------------------------

library(UpSetR)
upset_etoh = unique(c(targets,upr))
n = length(upset_etoh)
upset_etoh = data.frame(gene = upset_etoh,DE=character(n))
upset_etoh$Atf6_Targets = as.numeric(upset_etoh$gene %in% targets)
upset_etoh$Atf6_Sufficient_Targets = as.numeric(upset_etoh$gene %in% primary_targets)
upset_etoh$UPR_Genes = as.numeric(upset_etoh$gene %in% upr)

upset_etoh$DE = c('notDE','DE')[as.numeric(upset_etoh$gene %in% etoh_sig$ensembl_id) + 1]
upset_etoh$DE = as.factor(upset_etoh$DE)

Myfunc <- function(row){
  data <- row["DE"]=='DE'
}

sets = colnames(upset_etoh)[3:5]

#labels need to be added manually as a result of this loop and then moved above bars on illustrator
#for each condiition, bottom number is DE and top number is total
for (i in 0:1){
  for (j in 0:1){
    for (k in 0:1){
      print(paste(sets[1],i,sets[2],j,sets[3],k))
      print(summary(upset_etoh$Atf6_Targets==i & upset_etoh$Atf6_Sufficient_Targets==j & upset_etoh$UPR_Genes==k))        
      print(summary(upset_etoh$Atf6_Targets==i & upset_etoh$Atf6_Sufficient_Targets==j & upset_etoh$UPR_Genes==k & upset_etoh$DE =='DE'))        
    }
  }
}

pdf('figures/Fig 4/fig4_etoh_upset.pdf',width = 4,height = 4.8,onefile=FALSE)

upset(upset_etoh,
      mb.ratio = c(0.8,0.2),sets = rev(sets),keep.order = T,
      order.by = 'degree',decreasing = F,
      query.legend = 'top',
      sets.bar.color = c(upr_colour,natf6_colour,targets_colour),
      queries = list(
        list(query = Myfunc,active=T,color=etoh_colour,query.name='Differentially Expressed in EtOH')
      ),
      mainbar.y.max = 1180,
      text.scale = c(1.5,1,1,1,1.2,1),
      sets.x.label = 'Gene List Size',
      set_size.show = TRUE,
      set_size.scale_max = 2000,
      show.numbers = T)
library(grid)
library(gridtext)
grid.text("142/1090",x = 0.65, y=0.9, gp=gpar(fontsize=8))
grid.text("47/244",x = 0.70, y=0.85, gp=gpar(fontsize=8))
grid.text("84/354",x = 0.75, y=0.80, gp=gpar(fontsize=8))
grid.text("5/29",x = 0.8, y=0.75, gp=gpar(fontsize=8))
grid.text("10/25",x = 0.85, y=0.70, gp=gpar(fontsize=8))

dev.off()

# Figure 4 - UpSet plot - iAs ---------------------------------------------

library(UpSetR)
upset_ias = upset_etoh
upset_ias$DE = c('notDE','DE')[as.numeric(upset_ias$gene %in% ias_sig$ensembl_id) + 1]
upset_ias$DE = as.factor(upset_ias$DE)

Myfunc <- function(row){
  data <- row["DE"]=='DE'
}

sets = colnames(upset_ias)[3:5]

#labels need to be added manually as a result of this loop and then moved above bars on illustrator
#for each condiition, bottom number is DE and top number is total
for (i in 0:1){
  for (j in 0:1){
    for (k in 0:1){
      print(paste(sets[1],i,sets[2],j,sets[3],k))
      print(summary(upset_ias$Atf6_Targets==i & upset_ias$Atf6_Sufficient_Targets==j & upset_ias$UPR_Genes==k))        
      print(summary(upset_ias$Atf6_Targets==i & upset_ias$Atf6_Sufficient_Targets==j & upset_ias$UPR_Genes==k & upset_ias$DE =='DE'))        
    }
  }
}

pdf('figures/Fig 4/fig4_ias_upset.pdf',width = 4,height = 4.8,onefile=FALSE)

upset(upset_ias,
      mb.ratio = c(0.8,0.2),sets = rev(sets),keep.order = T,
      order.by = 'degree',decreasing = F,
      query.legend = 'top',
      sets.bar.color = c(upr_colour,natf6_colour,targets_colour),
      queries = list(
        list(query = Myfunc,active=T,color=ias_colour,query.name='Differentially Expressed in iAs')
      ),
      mainbar.y.max = 1180,
      text.scale = c(1.5,1,1,1,1.2,1),
      sets.x.label = 'Gene List Size',
      set_size.show = TRUE,
      set_size.scale_max = 2000,
      show.numbers = T)

grid.text("244/1090",x = 0.65, y=0.9, gp=gpar(fontsize=8))
grid.text("78/244",x = 0.70, y=0.85, gp=gpar(fontsize=8))
grid.text("141/354",x = 0.75, y=0.80, gp=gpar(fontsize=8))
grid.text("5/29",x = 0.8, y=0.75, gp=gpar(fontsize=8))
grid.text("7/25",x = 0.85, y=0.70, gp=gpar(fontsize=8))

dev.off()

# Figure 4 - GO Enrichment ------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
OrgDb = org.Hs.eg.db

de_etohias = intersect(ias_sig$ensembl_id,intersect(targets,etoh_sig$ensembl_id))
de_etohonly = targets[targets %in% etoh_sig$ensembl_id & !targets %in% ias_sig$ensembl_id]
de_iasonly = targets[targets %in% ias_sig$ensembl_id & !targets %in% etoh_sig$ensembl_id]
de_none = targets[!targets %in% etoh_sig$ensembl_id & !targets %in% ias_sig$ensembl_id]
write.csv(de_none,"DEnone.csv")
fig5 = list(etohias = get_human_aft6targets(de_etohias),
            etohonly = get_human_aft6targets(de_etohonly),
            iasonly = get_human_aft6targets(de_iasonly),
            none = get_human_aft6targets(de_none))
fig5_go = compareCluster(fig5,fun = "enrichGO",keyType = "ENSEMBL", OrgDb = OrgDb,ont = "BP",readable = TRUE)

create_revigo_table(fig5_go,'data/REVIGO/input/fig5_goidswithp.txt')
fig5_go = drop_from_revigo(fig5_go,"data/REVIGO/output/REVIGO_fig5.csv")
scale = 1.5
wrapsize = 38

fig5_plottable = for_plotting_comparecluster(fig5,fig5_go,4,wrapsize)

plot_fig5 = ggplot(fig5_plottable, aes(x = Cluster, y = Description)) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_color_continuous(name=str_wrap('Adjusted p-value', width = wrapsize),low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
  ylab(NULL) +
  xlab('')+
  scale_x_discrete(labels= c('DE in EtOH & iAs','DE in EtOH only','DE in iAs only','Remaining Atf6\nTargets'))+
  theme(legend.key.height = unit(30*scale,'points'),
        axis.text.y = element_text(size = 22*scale,angle=0, hjust=1),
        axis.text.x = element_text(angle = 45,hjust=1,size=18*scale),
        text = element_text(size=18*scale),
        axis.title.x =element_text(size=22*scale))+
  scale_size(name='Ratio of Genes\nin Gene List',range=c(6*scale, 18*scale),breaks=c(0.02,0.04,0.08,0.14))+
  theme(strip.text.y = element_blank())+
  geom_text(aes(label=Count,size=0.002,fontface='bold'),colour='white',show.legend = F)

png('figures/Fig 5/fig5_GO.png',width = 800*scale,height = 800*scale)
plot_fig5
dev.off()

save(fig5_go,file = 'data/RData/fig5_go.RData')



source('code/Functions.R')

# Curate UPR gene list ----------------------------------------------------

reactome_upr = read.delim('data/UPR files/Reactome UPR.txt',header = F)
reactome_upr = reactome_upr[reactome_upr$V4=='Unfolded Protein Response (UPR)',]
reactome_upr_human = reactome_upr$V1[reactome_upr$V6 == 'Homo sapiens']

go_erstress = read.delim('data/UPR files/GO 0034976 ER Stress.txt',header = F)
go_unfoldedprotein = read.delim('data/UPR files/GO 0006986 Unfolded Protein.txt',header = F)

go_erstress_uniprot = gsub('UniProtKB:','',unique(as.character(go_erstress$V1)))
go_unfoldedprotein_uniprot = gsub('UniProtKB:','',unique(as.character(go_unfoldedprotein$V1)))

library(clusterProfiler)
library(org.Hs.eg.db)
OrgDb = org.Hs.eg.db

go_erstress_ids = bitr(go_erstress_uniprot,fromType = "UNIPROT",toType = "ENSEMBL",OrgDb = OrgDb)
go_unfoldedprotein_ids = bitr(go_unfoldedprotein_uniprot,fromType = "UNIPROT",toType = "ENSEMBL",OrgDb = OrgDb)

go_erstress_human = unique(go_erstress_ids$ENSEMBL)
go_unfoldedprotein_human = unique(go_unfoldedprotein_ids$ENSEMBL)

upr_human = unique(c(reactome_upr_human,go_erstress_human,go_unfoldedprotein_human))
upr_human_tozebrafish = id_conv(upr_human,'human','zebrafish')
upr_human_tomouse = id_conv(upr_human,'human','mouse')

upr_homologs = merge(upr_human_tozebrafish,upr_human_tomouse)
upr = unique(upr_homologs$zebrafish)

erstress_zf = unique(upr_homologs$zebrafish[upr_homologs$human %in% go_erstress_human])
unfoldedprotein_zf = unique(upr_homologs$zebrafish[upr_homologs$human %in% go_unfoldedprotein_human])
reactome_zf = unique(upr_homologs$zebrafish[upr_homologs$human %in% reactome_upr_human])

make_3way_eulerr(erstress_zf,unfoldedprotein_zf,reactome_zf,'figures/Fig 2/fig2A-venn',c('#007d80',upr_colour,'#80fdff'),7)

save(upr,upr_homologs, file = 'data/RData/upr_files.RData')

# Figure 2 - Comparing Atf6 targets and UPR genes - GO Enrichment --------------------

targetandupr = intersect(targets,upr)
targetandupr_human = unique(targets_homologs$human[targets_homologs$zebrafish %in% targetandupr])
onlytarget = targets[!targets %in% upr]
onlytarget_human = unique(targets_homologs$human[targets_homologs$zebrafish %in% onlytarget])
onlyupr = upr[!upr %in% targets]
onlyupr_human = unique(upr_homologs$human[upr_homologs$zebrafish %in% onlyupr])

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
OrgDb = org.Hs.eg.db

fig2 = list(only_UPR = onlyupr_human,UPR_and_Target = targetandupr_human,only_Target = onlytarget_human)
fig2_go = compareCluster(fig2, keyType = "ENSEMBL", OrgDb = OrgDb, ont = "CC", readable = TRUE)

create_revigo_table(fig2_go,'data/REVIGO/input/fig2_goidswithp.txt')
fig2_go = drop_from_revigo(fig2_go,"data/REVIGO/output/REVIGO_fig2.csv")
scale = 1.5
wrapsize = 38

fig2_plottable = for_plotting_comparecluster(fig2,fig2_go,4,wrapsize)

plot_fig2 = ggplot(fig2_plottable, aes(x = Cluster, y = Description)) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_color_continuous(name=str_wrap('Adjusted p-value', width = wrapsize),low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
  ylab(NULL) +
  xlab('')+
  scale_x_discrete(labels= c('Only UPR gene','UPR gene and Atf6 target','Only Atf6 target'))+
  theme(legend.key.height = unit(30*scale,'points'),
        axis.text.y = element_text(size = 23*scale,angle=0, hjust=1),
        axis.text.x = element_text(angle = 45,hjust=1,size=18*scale),
        text = element_text(size=18*scale),
        axis.title.x =element_text(size=22*scale))+
  scale_size(name='Ratio of Genes\nin Gene List',range=c(6*scale, 22*scale),breaks=c(0.01,0.05,0.1,0.2))+
  theme(strip.text.y = element_blank())+
  geom_text(aes(label=Count,size=0.00045,fontface='bold'),colour='white',show.legend = F)

png('figures/Fig 2/fig2_GO.png',width = 800*scale,height = 800*scale)
plot_fig2
dev.off()

save(fig2_go,file = 'data/RData/fig2_go.RData')

# Figure 2 - only Atf6 Target - GO Enrichment -------------------------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
OrgDb = org.Hs.eg.db

fig2_onlytarget = list(only_Target = onlytarget_human)
fig2_onlytarget_go = compareCluster(fig2_onlytarget, keyType = "ENSEMBL", OrgDb = OrgDb, ont = "BP", readable = TRUE)

create_revigo_table(fig2_onlytarget_go,'data/REVIGO/input/fig2_onlytarget_goidswithp.txt')
fig2_onlytarget_go = drop_from_revigo(fig2_onlytarget_go,"data/REVIGO/output/REVIGO_fig2_onlytarget.csv")
scale = 1.5
wrapsize = 38

fig2_onlytarget_plottable = for_plotting_comparecluster(fig2_onlytarget,fig2_onlytarget_go,10,wrapsize)

plot_fig2_onlytarget = ggplot(fig2_onlytarget_plottable, aes(x = Cluster, y = Description)) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_color_continuous(name=str_wrap('Adjusted p-value', width = wrapsize),low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
  ylab(NULL) +
  xlab('')+
  scale_x_discrete(labels= 'Only Atf6 target')+
  theme(legend.key.height = unit(30*scale,'points'),
        axis.text.y = element_text(size = 23*scale,angle=0, hjust=1),
        axis.text.x = element_text(angle = 45,hjust=1,size=18*scale),
        text = element_text(size=18*scale),
        axis.title.x =element_text(size=22*scale))+
  scale_size(name='Ratio of Genes\nin Gene List',range=c(8*scale, 26*scale),breaks=c(0.01,0.02,0.03,0.04))+
  theme(strip.text.y = element_blank())+
  geom_text(aes(label=Count,size=0.01,fontface='bold'),colour='white',show.legend = F)

png('figures/Fig 2/fig2_onlytarget_GO.png',width = 700*scale,height = 800*scale)
plot_fig2_onlytarget
dev.off()

save(fig2_onlytarget_go,file = 'data/RData/fig2_onlytarget_go.RData')


# Figure 2 - Comparing Atf6 targets and UPR genes - Venn ------------------

targets_colour = '#17AA5A'
upr_colour = '#00c3c7'

make_2way_eulerr(targets,upr,'figures/Fig 2/fig2_venn',c(targets_colour,upr_colour),7)

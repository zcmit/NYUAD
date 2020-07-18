source('code/Functions.R')

# Get Promoters -----------------------------------------------------------

#get and save all proximal promoters (250bp upstream) from biomart
zebrafish_promoters = get_promoter(250,'zebrafish')
human_promoters = get_promoter(250,'human')
mouse_promoters = get_promoter(250,'mouse')

save(zebrafish_promoters,human_promoters,mouse_promoters,file='data/RData/promoters.RData')

# Clean Promoters ---------------------------------------------------------

zebrafish_clean_genes = get_bed_file_and_clean_genes(250,'zebrafish')
human_clean_genes = get_bed_file_and_clean_genes(250,'human')
mouse_clean_genes = get_bed_file_and_clean_genes(250,'mouse')

zebrafish_clean_promoters = zebrafish_promoters[zebrafish_promoters$`Gene stable ID` %in% zebrafish_clean_genes,]
human_clean_promoters = human_promoters[human_promoters$`Gene stable ID` %in% human_clean_genes,]
mouse_clean_promoters = mouse_promoters[mouse_promoters$`Gene stable ID` %in% mouse_clean_genes,]

save(zebrafish_clean_promoters,human_clean_promoters,mouse_clean_promoters,file='data/RData/cleanpromoters.RData')

# Find Atf6 Consensus Sequence --------------------------------------------

zebrafish_consensus_targets = find_consensus(zebrafish_clean_promoters)
human_consensus_targets = find_consensus(human_clean_promoters)
mouse_consensus_targets = find_consensus(mouse_clean_promoters)

human_tozebrafish_consensus = id_conv(human_consensus_targets$`Gene stable ID`,'human','zebrafish')
human_tomouse_consensus = id_conv(human_consensus_targets$`Gene stable ID`,'human','mouse')

human_consensus_homologs = merge(human_tozebrafish_consensus,human_tomouse_consensus)

consensus_targets_homologs = human_consensus_homologs[human_consensus_homologs$zebrafish %in% zebrafish_consensus_targets$`Gene stable ID`
                                                      & human_consensus_homologs$mouse %in% mouse_consensus_targets$`Gene stable ID`,]
consensus_targets = unique(consensus_targets_homologs$zebrafish)
length(consensus_targets)

save(consensus_targets_homologs,consensus_targets, file = 'data/RData/consensus_targets_files.RData')

# Find 'tftargets' Predicted Targets  -------------------------------------

# devtools::install_github("slowkow/tftargets")
library(tftargets)

# [1] "ENCODE"      "ITFP"        "Marbach2016"
# [4] "Neph2012"    "TRED"        "TRRUST"

#ITFP has 1 atf6 target, Marbach2016, Neph2012, TRED and TRRUST have many

#entrezid
tred = TRED$ATF6
neph = unlist(Neph2012)
nephnames = names(neph)
neph = unique(neph[grep('ATF6',nephnames)])

tftargets_entrez = c(tred,neph)

library(org.Hs.eg.db)
library(clusterProfiler)
tftargets1_ensembl = bitr(tftargets_entrez,fromType = "ENTREZID",toType = "ENSEMBL",OrgDb = "org.Hs.eg.db")
tftargets1_ensembl = unique(tftargets1_ensembl$ENSEMBL)

#symbols
trrust = TRRUST$ATF6
marbach = Marbach2016$ATF6
tftargets_symbols = unique(c(trrust,marbach))

library(biomaRt)
ens = useMart('ensembl')
ens = useDataset("hsapiens_gene_ensembl",mart=ens)
tftargets2_ensembl = getBM(attributes = c('ensembl_gene_id','external_gene_name'), filters = 'external_gene_name',values = tftargets_symbols,mart=ens)
tftargets2_ensembl = unique(tftargets2_ensembl$ensembl_gene_id)

tftargets_atf6 = unique(c(tftargets1_ensembl,tftargets2_ensembl))

human_tozebrafish_tftargets = id_conv(tftargets_atf6,'human','zebrafish')
human_tomouse_tftargets = id_conv(tftargets_atf6,'human','mouse')

tftargets_homologs = merge(human_tozebrafish_tftargets,human_tomouse_tftargets)
tftargets = unique(tftargets_homologs$zebrafish)
save(tftargets_homologs,tftargets, file = 'data/RData/tftargets_files.RData')



# Find Shoulders2013 human Atf6 OE targets --------------------------------

shoulders2013 = na.omit(read.csv('data/spreadsheets/shoulders2013.csv'))
shoulders2013$xbp1_padj = 10^shoulders2013$Q.value..XBP1s.
shoulders2013$atf6_padj = 10^shoulders2013$Q.value..ATF6.
shoulders2013$xbp1_atf6_padj = 10^shoulders2013$Q.value..XBP1s.and.ATF6

shoulders2013 = shoulders2013[,c(1:4,14:16,11:13)]
colnames(shoulders2013)[2:4] = c('xbp1_FC','atf6_FC','xbp1_atf6_FC')

shoulders2013_atf6targets = shoulders2013[shoulders2013$atf6_padj<0.05,]

library(clusterProfiler)
library(org.Hs.eg.db)
library(Biobase)

bitr_shoulders2013 = bitr(shoulders2013_atf6targets$NCBI.ID,"ENTREZID","ENSEMBL",OrgDb = org.Hs.eg.db)

shoulders2013_atf6targets = merge(shoulders2013_atf6targets,bitr_shoulders2013,by.x='NCBI.ID',by.y="ENTREZID")

shoulders2013_human = intersect(shoulders2013_atf6targets$ENSEMBL,human_clean_genes)

human_tozebrafish_shoulders2013 = id_conv(shoulders2013_human,'human','zebrafish')
human_tomouse_shoulders2013 = id_conv(shoulders2013_human,'human','mouse')

shoulders2013_homologs = merge(human_tozebrafish_shoulders2013,human_tomouse_shoulders2013) 
length(unique(shoulders2013_homologs$human))
#38 unique human genes

shoulders2013_humanconsensus_homologs = shoulders2013_homologs[shoulders2013_homologs$human %in% human_consensus_targets$`Gene stable ID`,] #38 unique zf genes
#remove 1 gene that are does not have 1-to-1 homology in human-mouse
mousehuman_shoulders = unique(shoulders2013_humanconsensus_homologs[,c('human','mouse')])
mousehuman_shoulders = mousehuman_shoulders[isUnique(mousehuman_shoulders$human),]
#32 human genes

shoulders2013_humanconsensus_homologs = shoulders2013_humanconsensus_homologs[shoulders2013_humanconsensus_homologs$human %in% mousehuman_shoulders$human,]

shoulders2013_targets = unique(shoulders2013_humanconsensus_homologs$zebrafish)

save(shoulders2013_homologs,shoulders2013_humanconsensus_homologs,shoulders2013_targets, file = 'data/RData/shoulders2013_files.RData')

# Find Rutkowski2008 MEF Atf6-dependent targets ---------------------------

rutkowski2008 = read.csv("data/spreadsheets/rutkowski2008.csv")
rutkowski2008_entrez = na.omit(as.numeric(rutkowski2008$EntrezGene))

library(clusterProfiler)
library(org.Mm.eg.db)
rutkowski2008_ensembl = bitr(rutkowski2008_entrez,"ENTREZID","ENSEMBL",OrgDb = org.Mm.eg.db)
rutkowski2008_ensembl = rutkowski2008_ensembl[rutkowski2008_ensembl$ENSEMBL %in% mouse_clean_genes,]

rutkowski2008 = merge(rutkowski2008,rutkowski2008_ensembl,by.x='EntrezGene',by.y='ENTREZID')
rutkowski2008 = rutkowski2008[,c(1,22,2,3,6:21)]
#remove genes with no change in any condition
rutkowski2008 = rutkowski2008[!(rutkowski2008$comp.8.WT.treated.vs.control == 'NC' & 
                                  rutkowski2008$comp.8.KO.treated.vs.WT.treated == 'NC' & 
                                  rutkowski2008$comp.24.WT.treated.vs.control == 'NC' & 
                                  rutkowski2008$comp.24.KO.treated.vs.WT.treated == 'NC'),]

#knocking out Atf6 should do something, so no NCs for KO vs  WT treated
rutkowski2008 = rutkowski2008[!(rutkowski2008$comp.8.KO.treated.vs.WT.treated == 'NC' |
                                  rutkowski2008$comp.24.KO.treated.vs.WT.treated == 'NC'),]
#Tm induces Atf6 which should regulate downstream targets, so no NCs for WT treated vs control
rutkowski2008 = rutkowski2008[!(rutkowski2008$comp.8.WT.treated.vs.control == 'NC' |
                                  rutkowski2008$comp.24.WT.treated.vs.control == 'NC'),]
                                  
#Atf6 targets have opposite response in KO treated, so are up/down or down/up at 8h and 24h
rutkowski2008 = rutkowski2008[((rutkowski2008$comp.8.WT.treated.vs.control == 'Up TM WT' &
                                  rutkowski2008$comp.8.KO.treated.vs.WT.treated == 'Dn TM KO') |
                                 (rutkowski2008$comp.8.WT.treated.vs.control == 'Dn Tm WT' &
                                    rutkowski2008$comp.8.KO.treated.vs.WT.treated == 'Up TM KO')),]

rutkowski2008 = rutkowski2008[((rutkowski2008$comp.24.WT.treated.vs.control == 'Up TM WT' &
                                  rutkowski2008$comp.24.KO.treated.vs.WT.treated == 'Dn TM KO') |
                                 (rutkowski2008$comp.24.WT.treated.vs.control == 'Dn Tm WT' &
                                    rutkowski2008$comp.24.KO.treated.vs.WT.treated == 'Up TM KO')),]

rutkowski2008_mouse = unique(rutkowski2008$ENSEMBL)
rutkowski2008_mouse_tohuman = id_conv(rutkowski2008_mouse,'mouse','human')
rutkowski2008_human_tozebrafish = id_conv(unique(rutkowski2008_mouse_tohuman$human),'human','zebrafish')

rutkowski2008_homologs = merge(rutkowski2008_human_tozebrafish,rutkowski2008_mouse_tohuman) #39 unique zf genes
#33 unique mouse genes
length(unique(rutkowski2008_homologs$mouse))


rutkowski2008_mouseconsensus_homologs = rutkowski2008_homologs[rutkowski2008_homologs$mouse %in% mouse_consensus_targets$`Gene stable ID`,] #33 unique zf genes

#remove 1 gene that are does not have 1-to-1 homology in human-mouse
mousehuman_rutkowski = unique(rutkowski2008_mouseconsensus_homologs[,c('human','mouse')])
mousehuman_rutkowski = mousehuman_rutkowski[isUnique(mousehuman_rutkowski$mouse),]
#26 mouse genes

rutkowski2008_mouseconsensus_homologs = rutkowski2008_mouseconsensus_homologs[rutkowski2008_mouseconsensus_homologs$mouse %in% mousehuman_rutkowski$mouse,]

rutkowski2008_targets = unique(rutkowski2008_mouseconsensus_homologs$zebrafish)

save(rutkowski2008_homologs,rutkowski2008_mouseconsensus_homologs,rutkowski2008_targets, file = 'data/RData/rutkowski2008_files.RData')

# Combining Shoulders2013 and Rutkowski2008 -------------------------------

length(unique(c(shoulders2013_humanconsensus_homologs$human,rutkowski2008_mouseconsensus_homologs$mouse)))

# Supplementary Figure 1 - Excluding 'tftargets' - GO Enrichment --------------------------

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
OrgDb = org.Hs.eg.db

suppfig1 = list(Consensus_Targets=unique(consensus_targets_homologs$human),Database_Targets = unique(tftargets_homologs$human))
suppfig1_go = compareCluster(suppfig1,fun = "enrichGO",keyType = "ENSEMBL", OrgDb = OrgDb,ont = "BP",readable = TRUE)

create_revigo_table(suppfig1_go,'data/REVIGO/input/suppfig1_goidswithp.txt')
suppfig1_go = drop_from_revigo(suppfig1_go,"data/REVIGO/output/REVIGO_suppfig1.csv")
scale = 1.5
wrapsize = 38

suppfig1_plottable = for_plotting_comparecluster(suppfig1,suppfig1_go,6,wrapsize)
plot_suppfig1 <- ggplot(suppfig1_plottable, aes(x = Cluster, y = Description)) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_color_continuous(name=str_wrap('Adjusted p-value', width = wrapsize),
                         low="red", high="blue", guide=guide_colorbar(reverse=TRUE)) +
  ylab(NULL) +
  xlab('Gene List') +
  scale_x_discrete(labels= c('Consensus','Database\nTargets')) +
  theme(plot.title = element_text(hjust = 0.5,size=14*scale)) +
  theme(axis.text.y = element_text(size = 18*scale,angle=0, hjust=1), 
        axis.text.x = element_text(angle = 45,hjust=1,size=18*scale),
        text = element_text(size=18*scale), 
        axis.title.x = element_text(size=22*scale)) +
  scale_size(name='Ratio of Genes\nin Gene List',range=c(6*scale, 16*scale),breaks=c(0.025,0.03,0.04,0.05)) +
  theme(strip.text.y = element_blank()) +
  geom_text(aes(label=Count,size=0.01,fontface='bold'),colour='white',show.legend = F)

png('figures/Supp fig 1/suppfig1_GO.png',width = 800*scale,height = 800*scale)
plot_suppfig1
dev.off()

save(suppfig1_go,file = 'data/RData/suppfig1_go.RData')

# Supplementary Figure 1 - Excluding 'tftargets' - Venn -------------------

consensus_targets_colour = '#17AA5A'
tftargets_colour = '#5D47CE'

make_2way_eulerr(consensus_targets,tftargets,'figures/Supp fig 1/suppfig1_venn',c(consensus_targets_colour,tftargets_colour),7)

targets_2=unique(Targets_2$Gene)
targets=unique(Atf6PutTargets_all_1498$Gene_Name)
targets_colour = '#17AA5A'
targets2_colour = '#808000'
make_2way_eulerr(targets,targets_2,'figures/Supp fig 1/suppfig1b_venn',c(targets_colour,targets2_colour),7)
# Figure 1 - Curating Atf6 target list - GO Enrichment ----------------------------------------------------------------

targets_homologs = unique(rbind(consensus_targets_homologs,shoulders2013_humanconsensus_homologs,rutkowski2008_mouseconsensus_homologs))
targets = unique(targets_homologs$zebrafish)
save(targets_homologs,targets, file = 'data/RData/targets_files.RData')

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
OrgDb <- org.Hs.eg.db

fig1 = list(Putative=unique(targets_homologs$human),
            Consensus=unique(consensus_targets_homologs$human),
            Human_Microarray=unique(shoulders2013_humanconsensus_homologs$human),
            Mouse_Microarray=unique(rutkowski2008_mouseconsensus_homologs$human))
fig1_go = compareCluster(fig1,fun = "enrichGO",keyType = "ENSEMBL", OrgDb = OrgDb,ont = "BP",readable = TRUE)

create_revigo_table(fig1_go,'data/REVIGO/input/fig1_goidswithp.txt')
fig1_go = drop_from_revigo(fig1_go,"data/REVIGO/output/REVIGO_fig1.csv")
scale = 1.5
wrapsize = 38

fig1_plottable = for_plotting_comparecluster(fig1,fig1_go,6,wrapsize)
mouserow = data.frame(c('Mouse_Microarray',fig1_plottable[1,2:3],numeric(6),''))
colnames(mouserow) = colnames(fig1_plottable)
fig1_plottable = rbind(fig1_plottable,mouserow)
fig1_plottable$Cluster = factor(fig1_plottable$Cluster,unique(fig1_plottable$Cluster),ordered = T)


plot_fig1 = ggplot(fig1_plottable, aes(x = Cluster, y = Description)) + 
  geom_point(data = fig1_plottable[!(fig1_plottable$Count==''),],aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_color_continuous(name=str_wrap('Adjusted p-value', width = wrapsize),low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
  ylab(NULL) +
  xlab('')+
  scale_x_discrete(limits = unique(fig1_plottable$Cluster),labels= c('Putative Atf6 Targets','Conserved Atf6 Targets','Atf6 Human OE Targets','Atf6 Mouse KO Targets'))+
  theme(legend.key.height = unit(30*scale,'points'),
        axis.text.y = element_text(size = 23*scale,angle=0, hjust=1),
        axis.text.x = element_text(angle = 45,hjust=1,size=18*scale),
        text = element_text(size=18*scale),
        axis.title.x =element_text(size=22*scale))+
  scale_size(name='Ratio of Genes\nin Gene List',range=c(6*scale, 26*scale),breaks=c(0,0.001,0.01,0.05,0.1,0.2))+
  theme(strip.text.y = element_blank())+
  geom_text(aes(label=Count,size=0.00075,fontface='bold'),colour='white',show.legend = F)

png('figures/Fig 1/fig1_go.png',width = 800*scale,height = 800*scale)
plot_fig1
dev.off()

save(fig1_go,file = 'data/RData/fig1_go.RData')

# Figure 1 - Curating Atf6 target list - Venn -----------------------------

consensus_targets_colour = '#17AA5A'
shoulders2013_colour = '#f8cf5a'
rutkowski2008_colour = '#f55859'

make_3way_eulerr(consensus_targets,shoulders2013_targets,rutkowski2008_targets,
                 'figures/Fig 1/fig1_venn',c(consensus_targets_colour,shoulders2013_colour,rutkowski2008_colour),7)

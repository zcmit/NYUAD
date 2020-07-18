# Get Promoters -----------------------------------------------------------

#get all promoters for each species from biomart

get_promoter <- function(bp,species) {
  library("biomaRt")
  mart=useMart("ENSEMBL_MART_ENSEMBL")
  if (species=='zebrafish'){
    species = 'drerio_gene_ensembl'
  } else if (species=='human'){
    species = 'hsapiens_gene_ensembl'
  } else if (species=='mouse'){
    species = 'mmusculus_gene_ensembl'
  }
  ens=useDataset(species,mart)
  chroms = getBM('chromosome_name',mart=ens)
  chroms = chroms[,1]
  chroms = chroms[-grep('CHR|K|G|J|MT',chroms)]
  
  res=data.frame()
  for (i in 1:length(chroms)){
    res_current = getBM(attributes=c("ensembl_gene_id",
                                     "external_gene_name",
                                     "transcript_flank"),
                        filters=c("upstream_flank","chromosome_name"),
                        values=list(upstream_flank=bp,chroms[i]),
                        mart=ens, bmHeader = T,checkFilters = FALSE)
    res=rbind(res,res_current)
  }
  return(res)
}

# Find Clean Promoters ----------------------------------------------------

get_bed_file_and_clean_genes = function(bp,species){
  if (species=='zebrafish'){
    species = 'drerio_gene_ensembl'
    genome = 'danRer11'
  } else if (species=='human'){
    species = 'hsapiens_gene_ensembl'
    genome = 'hg38'
  } else if (species=='mouse'){
    species = 'mmusculus_gene_ensembl'
    genome = 'mm10'
  }
  
  library(genomation)
  library(rtracklayer)
  library(GenomicRanges)
  
  mySession <- browserSession()
  genome(mySession) = genome
  track.names = trackNames(ucscTableQuery(mySession))
  
  mySession = browserSession("UCSC")
  genome(mySession) = genome
  grange <- GRangesForUCSCGenome(genome)
  bed <- getTable(ucscTableQuery(mySession, track="refSeqComposite",
                     range=grange,table='refGene'))
  
  starts = as.character(bed$exonStarts)
  ends = as.character(bed$exonEnds)
  
  starts = strsplit(starts,',')
  ends = strsplit(ends,',')
  
  blocklengths = character()
  
  for (i in 1:length(starts)){
    lengths = as.numeric(ends[[i]])-as.numeric(starts[[i]])
    lengths = paste(paste(lengths,collapse = ',',sep=''),',',sep='')
    blocklengths = c(blocklengths,lengths)
  }
  
  blockstarts = character()
  
  for (i in 1:length(starts)){
    bstarts = as.numeric(starts[[i]])-bed$txStart[i]
    bstarts = paste(paste(bstarts,collapse = ',',sep=''),',',sep='')
    blockstarts = c(blockstarts,bstarts)
  }
  
  bed$blocklengths = blocklengths
  bed$blockstarts = blockstarts
  bed$rgb = 0
  
  #BED file columns: chrom,chromstart,chromend,name,score,strand,thickstart,thickend,itemrgb,blockcount,blocksizes,blockstarts
  
  bed = bed[,c('chrom','txStart','txEnd','name','score','strand','cdsStart','cdsEnd','rgb','exonCount','blocklengths','blockstarts')]
  bedfilename = paste('data/BED files/',genome,'.bed',sep='')
  
  write.table(bed,bedfilename,quote=FALSE,sep='\t',row.names = FALSE,col.names = FALSE)
  
  unclean = readTranscriptFeatures(bedfilename,unique.prom = FALSE,up.flank = bp,down.flank = 0)
  unclean = as.data.frame(unclean)
  uncleanpromoters = unclean[unclean$group_name=='promoters',]
  
  clean = readTranscriptFeatures(bedfilename,unique.prom = TRUE,up.flank = bp,down.flank = 0)
  clean = as.data.frame(clean)
  cleanpromoters = clean[clean$group_name=='promoters',]
  
  allcleanpromoters = merge(cleanpromoters,uncleanpromoters,by=c(1:8))
  allcleanpromoters = allcleanpromoters[,c(2:7,10)]
  
  library("biomaRt")
  mart=useMart("ENSEMBL_MART_ENSEMBL")
  ens=useDataset(species,mart)
  
  clean_ensembl_mrna = getBM(attributes="ensembl_gene_id",
                            filters="refseq_mrna",
                            values=grep('NM',allcleanpromoters$name.y,value = T),
                            mart=ens)
  clean_ensembl_ncrna = getBM(attributes="ensembl_gene_id",
                             filters="refseq_ncrna",
                             values=grep('NR',allcleanpromoters$name.y,value = T),
                             mart=ens)

  clean_ensembl = unique(c(clean_ensembl_mrna[,1],clean_ensembl_ncrna[,1]))
  return(clean_ensembl)
}

# Find Atf6 Consensus Sequence-----------------------------------------------------

find_consensus = function(promoterdf) {
  consensus = promoterdf[grep("CCACG",promoterdf$`Flank (Transcript)`),]
  consensus_complement = promoterdf[grep("CGTGG",promoterdf$`Flank (Transcript)`),]
  allconsensus=rbind(consensus,consensus_complement)
  return(unique(allconsensus[,c(1,3)]))
}

# Convert Ensembl IDs for Species -----------------------------------------

id_conv = function(ensembls,speciesfrom,speciesto){
  library(biomaRt)
  ens = useMart('ensembl')
  if (speciesfrom == 'human'){
    ensfrom = 'hsapiens_gene_ensembl'
  } else if (speciesfrom == 'mouse'){
    ensfrom = 'mmusculus_gene_ensembl'
  } else if (speciesfrom == 'zebrafish'){
    ensfrom = 'drerio_gene_ensembl'
  }
  
  if (speciesto == 'human'){
    ensto = 'hsapiens_gene_ensembl'
  } else if (speciesto == 'mouse'){
    ensto = 'mmusculus_gene_ensembl'
  } else if (speciesto == 'zebrafish'){
    ensto = 'drerio_gene_ensembl'
  }
  
  martfrom = useDataset(ensfrom,ens)
  martto = useDataset(ensto,ens)
  
  newids = getLDS(attributes = 'ensembl_gene_id', filters = "ensembl_gene_id", values = ensembls, 
                  mart = martfrom, attributesL = 'ensembl_gene_id', martL = martto)
  
  colnames(newids) = c(speciesfrom,speciesto)
  return(newids)
}

# Create table input for REVIGO -------------------------------------------

create_revigo_table = function(go, filename){
  db = as.data.frame(go)
  db = db[order(db$p.adjust),]
  keepgos = data.frame(ID=vector(),p.adjust=vector())
  
  #create list of unique GO IDs with lowest p-value of all groups
  for (i in 1:nrow(db)){
    if (!db$ID[i] %in% keepgos$ID){
      keepgos = rbind(keepgos,db[i,])
    }
  }
  write.table(keepgos[,c(2,6)],quote = F,sep = ' ',row.names = F,col.names=F,file = filename)
}


# Drop redundant IDs from REVIGO ------------------------------------------

drop_from_revigo = function(current_go,revigo_filename){
  revi = read.csv(revigo_filename)
  db = as.data.frame(current_go)
  ids = db$ID[!db$ID %in% revi$term_ID[revi$eliminated==0]]
  dropped_go = dropGO(current_go,term = ids)
  return(dropped_go)
}

# Create plot table for CompareCluster ------------------------------------

for_plotting_comparecluster = function(list,finalgo,no_top_IDs,text_wrapsize){
  finaldb = as.data.frame(finalgo)
  topids = vector()
  names = names(list)
  
  #extract topids for each group
  for (i in 1:length(names)){
    name = names[i]
    top = finaldb[finaldb$Cluster==name,]
    top = top[order(top$p.adjust),]
    top = top$ID[1:no_top_IDs]
    topids = c(topids,top)
  }
  
  dropids = unique(finaldb$ID)
  dropids = dropids[!dropids %in% topids]
  #subset GO to only include top IDs
  plotdb = dropGO(finalgo,term=dropids)
  
  library(dplyr)
  library(stringr)
  library(forcats)
  
  plot_table = as.data.frame(plotdb)
  plot_table$GeneRatio = as.numeric(sub('/.*','',plot_table$GeneRatio))/as.numeric(sub('.*/','',plot_table$GeneRatio))
  plot_table$BgRatio = as.numeric(sub('/.*','',plot_table$BgRatio))/as.numeric(sub('.*/','',plot_table$BgRatio))
  plot_table$Description = str_wrap(plot_table$Description, width = text_wrapsize)
  plot_table$p.adjust = as.numeric(plot_table$p.adjust)
  plot_table$Description = factor(plot_table$Description,levels=rev(unique(plot_table$Description)))
  
  return(plot_table)
}

# Make venn diagrams -------------------------------------------------

make_2way_eulerr = function(list1,list2,name,colours,scale){
  library(eulerr)
  n1 = length(list1)
  n2 = length(list2)
  n12 = length(intersect(list1,list2))
  
  comb = c('1'=n1,'2'=n2,'1&2'=n12)
  venn = euler(comb, input = 'union', shape = 'ellipse')
  save_venn_plots(venn,name,colours,scale)
}

make_3way_eulerr = function(list1,list2,list3,name,colours,scale){
  n1 = length(list1)
  n2 = length(list2)
  n3 = length(list3)
  
  n12 = length(intersect(list1,list2))
  n23 = length(intersect(list2,list3))
  n13 = length(intersect(list1,list3))
  
  n123 = length(intersect(intersect(list1,list2),list3))
  
  comb = c('1'=n1,'2'=n2,'3'=n3,
           '1&2'=n12,'2&3'=n23,'1&3'=n13,
           '1&2&3'=n123)
  venn = euler(comb, input = 'union', shape = 'ellipse')
  save_venn_plots(venn,name,colours,scale)
}

save_venn_plots = function(eulerr,name,colours,scale){
  labelname = paste(name,'_labels.png',sep='')
  nolabelname = paste(name,'.png',sep='')
  #to see labels to add manually in illustrator
  png(labelname)
  print({plot(eulerr,
              labels=FALSE,
              quantities = TRUE,
              fills = colours)
  })
  dev.off()

  png(nolabelname,width = 400*scale,height = 400*scale)
  print({plot(eulerr,edges = list(lwd=1*scale),
              labels=FALSE,
              quantities = FALSE,
              fills = colours)
  })
  dev.off()
}

# Process RNAseq count files ----------------------------------------------

rnaseq_analysis = function(filename){
  rnaseq = read.delim(filename)
  rnaseq = na.omit(rnaseq)
  rnaseq = unique(rnaseq[,c('ensembl_id','log2FoldChange','padj')])
  colnames(rnaseq)[3] = 'pValue'
  return(rnaseq)
}

# Make violin plot for RNAseq ---------------------------------------------

plot_violin = function(name,rnaseqcolour,rnaseq_table,scale){
  library(ggplot2)
  library(ggpubr)
  load('data/RData/targets_files.RData')
  load('data/RData/upr_files.RData')
  targets_colour = "#17AA5A"
  upr_colour = "#00c3c7"
  
  namevec = c(name,'Atf6 Targets','UPR Genes')
  all = rnaseq_table
  all$type = namevec[1]
  violin_targets = rnaseq_table[rnaseq_table$ensembl_id %in% targets,]
  violin_targets$type = namevec[2]
  violin_upr = rnaseq_table[rnaseq_table$ensembl_id %in% upr,]
  violin_upr$type = namevec[3]
  
  violin_db = rbind(all,violin_targets,violin_upr)
  comparisons = list( namevec[1:2], namevec[2:3], namevec[c(1,3)] )
  
  namevec[1] = paste(namevec[1],'\n(',nrow(all),')',sep='')
  namevec[2] = paste(namevec[2],'\n(',nrow(violin_targets),')',sep='')
  namevec[3] = paste(namevec[3],'\n(',nrow(violin_upr),')',sep='')
  
  g = ggplot(violin_db,aes(x=type,y=log2FoldChange,fill=type))+
    geom_hline(yintercept = 0)+
    geom_violin()+ geom_boxplot(width=0.04,fill='white')+
    stat_compare_means(comparisons = comparisons,size=5*scale,label="p.signif")+
    scale_fill_manual(values=c(rnaseqcolour,targets_colour,upr_colour))+
    scale_x_discrete(labels=namevec)+
    xlab('')+
    ylab(expression(paste('log'[2],' Fold Change',sep='')))+
    scale_y_continuous(breaks=c(-5,0,5,10))+
    theme_classic()+
    theme(axis.text.y = element_text(size=12*scale),
          axis.title = element_text(size = 16*scale),
          plot.title = element_text(hjust=0.5,size=16*scale),
          axis.title.x = element_text(size=16*scale,margin = margin(t=10*scale)),
          axis.text.x = element_text(size = 12*scale),
          legend.position = 'none')
  return(g)
}




# Get human Atf6 targets from zebrafish --------------------------------------------------

get_human_aft6targets = function(zftargets){
  x = unique(targets_homologs$human[targets_homologs$zebrafish %in% zftargets])
  return(x)
}

# Get gene name from Ensembl ID -------------------------------------------

get_genename = function(ensembls,organism){
  library(biomaRt)
  ens = useEnsembl('ensembl', mirror = 'useast')
  if (organism == 'zf'){
    ens = useDataset("drerio_gene_ensembl",mart=ens)
  }
  if (organism == 'human'){
    ens = useDataset("hsapiens_gene_ensembl",mart=ens)
  }
  if (organism == 'mouse'){
    ens = useDataset("mmusculus_gene_ensembl",mart=ens)
  }
  names = getBM(attributes = c('ensembl_gene_id','external_gene_name'), filters = 'ensembl_gene_id',values = ensembls,mart=ens)
  return(names)
}

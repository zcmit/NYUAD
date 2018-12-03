library('biomaRt')
listMarts()    # to see which database options are present
ensembl <- useMart('ensembl')  # using ensembl database data
listDatasets(ensembl)     # function to see which datasets are present in ensembl

mart_mm <- useDataset('mmusculus_gene_ensembl', useMart('ensembl'))
mart_hg <- useDataset('hsapiens_gene_ensembl', useMart('ensembl'))
mart_zf <- useDataset('drerio_gene_ensembl', useMart('ensembl'))
#  mmusculus_gene_ensembl             Mus musculus genes (GRCm38.p4)         GRCm38.p4
#  drerio_gene_ensembl                 Danio rerio genes (GRCz10)            GRCz10
#  hsapiens_gene_ensembl             Homo sapiens genes (GRCh38.p7)         GRCh38.p7
# human: 'hgnc_symbol', mouse: 'mgi_symbol', zebrafish: 'zfin_symbol' ('zfin_id_symbol')
listFilters(mart_mm, 10)  # check which filters are available

getBM(attributes=c("ensembl_gene_id", "external_gene_id"),filters='external_gene_id',values='pcna', mart=mart_mm) 
# fuction to get  gene id's and gene name from data base

setwd("~/NYU Ad/Mouse_RNAseq")
CellCycle_ids <- read.delim('VJ cell cycle gene panel.txt', header=F, sep='\\')
CellCycle_ids <- as.character(CellCycle_ids$V1)
Immune_ids <- read.delim('YC immune gene panel.txt', header=F, sep='\\')
Immune_ids <- as.character(Immune_ids$V1)
getBM(filters= "mgi_symbol", attributes= c( "hgnc_symbol","ensembl_gene_id"),
      values=CellCycle_ids,mart= mart_mm)
### human: 'hgnc_symbol', mouse: 'mgi_symbol', zebrafish: 'zfin_id_symbol'
head(listAttributes(mart_mm))
chr1_genes <- getBM(attributes=c('ensembl_gene_id',
              'ensembl_transcript_id','hgnc_symbol','chromosome_name','start_position','end_position'), 
              filters = 'chromosome_name', values ="1", mart = mart_mm)


### Homologous Gene IDs conversion
setwd("~/NYU Ad/UPR_Arsenic")
upr_hg_ids <- read.table('UPRGenes_All_Human.txt', header = F, sep='\t')
upr_hg_ids <- as.character(upr_hg_ids$V1)
tab <- getLDS(attributes = "hgnc_symbol", filters = "hgnc_symbol", values = upr_hg_ids, 
              mart = mart_hg, attributesL = 'zfin_id_symbol', martL = mart_zf)

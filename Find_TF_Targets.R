

---
title: "Finding TF targeted genes"
output: html_notebook
---


To find targeted genes of TF. 

1. checking the ChIP-chip and ChIP-Seq data from GEO and ENCODE for your TF. This is the direct evidence that the TF binds to the targets. But the TF binding is a dynamic/tissue specific process, so make sure the ChIP-chip and ChIP-Seq experiment in the similar tissue/cell type you are interested in. 

2. checking the motif of that TF, you can use the motif to scan the gene promoters. If the promoter contain the motif (usually this region should also be conserves across multiple species), then the corresponding gene might be the target of that TF. This is just a computational prediction. It usually gives more targets than reality. 

##### Here is a R package, 'tftargets', provides access to query a particular TF and find its targets in human from different curated databases.

```{r}
devtools::install_github("slowkow/tftargets") ### install the package
library(tftargets)    ### load library
length(TRED)          ### returns 133, worked
```


This package contains the following datasets:

```{r}
text_tbl <- data.frame(
  Dataset = c("TRED", "ITFP", "ENCODE", 'Neph2012', 'TRRUST', 'Marbach2016'),
  Gene_Identifier = c("ENTREZ","HGNC Symbol/Alias", "ENTREZ", "HGNC Symbol/Alias" ,"HGNC Symbol/Alias",
                      "HGNC Symbol/Alias"),
  Structure = c('list', 'list', 'list', 'nested list', 'list', 'list'),
  Number_TF = c(133, 1974, 157, 536, 748, 643))
kable(text_tbl) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)
```

Easily collect targeted genes with Human Gene Symbol

```{r}
### grammar: Dataset_name$gene_name. after $, 'tab' is able to check what genes are in the database
TRED$ATF6  # returns a list of ATF6 targeted genes in TRED, Entrez ID
ENCODE$ATF3
Nrf2 <- arbach2016$NFE2L2
```

Convert Entrez ID to gene symbol

```{r}
library('biomaRt')
ensembl <- useMart('ensembl')  #using ensembl database data
mart_hg <- useDataset('hsapiens_gene_ensembl', useMart('ensembl')) #make human genome mart
entrez_atf6_TRED <- TRED$ATF6  #creat entrez IDs of atf6 targeted genes as one variable
getBM(attributes=c('entrezgene', "hgnc_symbol"),filters='entrezgene',
      values=entrez_atf6_TRED, mart=mart_hg)   #ID converstion
tab <- getBM(attributes=c('entrezgene', "hgnc_symbol"),filters='entrezgene',
      values=entrez_atf6_TRED, mart=mart_hg)  #save it as a varialbe for output
```






library(pathview)
library(org.Hs.eg.db)

# Pathway visualisatie (pathview)

## Gene list maken
gene_list <- res$log2FoldChange
names(gene_list) <- rownames(res)

# NA's verwijderen
gene_list <- gene_list[!is.na(gene_list)]


## Pathway kiezen
pathway.id = "hsa05323" #RA pathway


## Pathview   
pathview(
  gene.data = gene_list,
  pathway.id = "hsa05323",
  species = "hsa",
  gene.idtype = "SYMBOL",
  limit = list(gene = 5),
  low = list(gene = "blue"),    # downregulated
  mid = list(gene = "white"),   # geen verschil
  high = list(gene = "red"),    # upregulated
  out.suffix = "RA_pathway"
)
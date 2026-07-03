library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(KEGGREST)
library(clusterProfiler)

# Set working directory
setwd(...)

# KEGG analyse:

# Selecteer significante genen (DESeq2 resultaten voorbereiden)
res <- as.data.frame(resultaten)
res <- res[!is.na(res$padj), ]

# Upregulated genen
up_genes <- rownames(res[(res$log2FoldChange > 1) & (res$padj < 0.05), ])

# Downregulated genen
down_genes <- rownames(res[(res$log2FoldChange < -1) & (res$padj < 0.05), ])

# Gene symbols omzetten naar ENTREZ IDs
## Upregulated
entrez_up <- mapIds(org.Hs.eg.db,
                    keys = up_genes,
                    column = "ENTREZID",
                    keytype = "SYMBOL",  # pas dit aan indien nodig!
                    multiVals = "first")

## verwijder NA's
entrez_up <- na.omit(entrez_up)

## Downregulated
entrez_down <- mapIds(org.Hs.eg.db,
                      keys = down_genes,
                      column = "ENTREZID",
                      keytype = "SYMBOL",  # pas dit aan indien nodig!
                      multiVals = "first")

## verwijder NA's
entrez_down <- na.omit(entrez_down)


# KEGG enrichment:

## KEGG gene-pathway mapping ophalen
kegg_data <- keggLink("pathway", "hsa")

kegg_df <- data.frame(
  pathway = sub("path:", "", kegg_data),
  gene = sub("hsa:", "", names(kegg_data))
)


# Overlap bepalen
## Upregulated
common_up <- intersect(entrez_up, kegg_df$gene)

## Downregulated
common_down <- intersect(entrez_down, kegg_df$gene)


# KEGG enrichment uitvoeren
## Upregulated
kegg_up <- enricher(
  gene = common_up,
  TERM2GENE = kegg_df,
  universe = kegg_df$gene)

## Downregulated
kegg_down <- enricher(
  gene = common_down,
  TERM2GENE = kegg_df,
  universe = kegg_df$gene)


# Resultaten netter maken
## Pathway namen ophalen
path_names <- keggList("pathway", "hsa")

## Upregulated
kegg_up_res <- as.data.frame(kegg_up)

kegg_up_res$Description <- path_names[kegg_up_res$ID]

## Downregulated
kegg_down_res <- as.data.frame(kegg_down)

kegg_down_res$Description <- path_names[kegg_down_res$ID]


# Filteren op significantie
## Upregulated
sig_kegg_up <- kegg_up_res[kegg_up_res$p.adjust < 0.05, ]

## Downregulated
sig_kegg_down <- kegg_down_res[kegg_down_res$p.adjust < 0.05, ]


# Plotten:
## Upregulated
sig_kegg_up %>%  
  top_n(10, wt = -p.adjust) %>%  
  ggplot(aes(x = Count,             
             y = reorder(Description, Count),             
             colour = p.adjust,             
             size = Count)) +  
  geom_point() +  
  scale_colour_gradient(low = "pink", high = "red") +  
  labs(    
    title = "KEGG pathways (upregulated genen in RA)",    
    x = "Aantal genen",    
    y = "Pathway",    
    colour = "Adjusted p-value",    
    size = "Aantal genen"  
  )

## Downregulated
sig_kegg_down %>%  
  top_n(10, wt = -p.adjust) %>%  
  ggplot(aes(x = Count,             
             y = reorder(Description, Count),             
             colour = p.adjust,             
             size = Count)) +  
  geom_point() +  
  scale_colour_gradient(low = "lightblue", high = "blue") +  
  labs(    
    title = "KEGG pathways (downregulated genen in RA)",    
    x = "Aantal genen",    
    y = "Pathway",    
    colour = "Adjusted p-value",    
    size = "Aantal genen"  
  )

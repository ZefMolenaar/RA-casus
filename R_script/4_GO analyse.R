library(goseq)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

# Set working directory
setwd(...)

# Resultaten opslaan als dataframe
res <- as.data.frame(resultaten)

# Verwijder NA's
res <- res[!is.na(res$padj), ]



# Up- en downregulated genen apart definiëren
up_genes <- (res$log2FoldChange > 1) & (res$padj < 0.05)
down_genes <- (res$log2FoldChange < -1) & (res$padj < 0.05)

# Maak vectors voor goseq
up_vec <- as.integer(up_genes)
down_vec <- as.integer(down_genes)

names(up_vec) <- rownames(res)
names(down_vec) <- rownames(res)



#GO voor upregulated genen:
## Corrigeren op lengte voor GO-analyse met model organisme hg19:
pwf_up <- nullp(up_vec,"hg19","geneSymbol")

## Voor elke GO-term kijken of er meer significante genen in deze categorie zijn
## dan verwacht
GO_up <- goseq(pwf_up,"hg19","geneSymbol")

## Resultaten bekijken
### Corrigeren voor multiple testing
GO_up$padj <- p.adjust(GO.up$over_represented_pvalue, method = "BH")

sig_GO_up <- GO_up[GO_up$padj < 0.05, ]


## Plot maken
### 10 meest significante GO-termen kiezen
### Hoeveel % van de genen in die GO-termen zijn significant
sig_GO_up %>% 
  top_n(10, wt = -padj) %>% 
  mutate(hitsPerc = numDEInCat * 100 / numInCat) %>% 
  ggplot(aes(x = hitsPerc, 
             y = reorder(term, hitsPerc), 
             colour = padj, 
             size = numDEInCat)) +
  geom_point() +
  scale_colour_gradient(low = "pink", high = "red") +
  expand_limits(x = 0) +
  labs(title = "GO enrichment (upregulated genen in RA)",
       x = "Percentage DE genen (%)", 
       y = "GO term (biologische processen", 
       colour = "Adjusted p-value", 
       size = "Aantal DE genen")



#GO voor downregulated genen:
## Identiek aan GO voor upregulated, met enkele aanpassingen aan naamgeving
pwf_down <- nullp(down_vec,"hg19","geneSymbol")

GO_down <- goseq(pwf_down,"hg19","geneSymbol")

GO_down$padj <- p.adjust(GO_down$over_represented_pvalue, method = "BH")

sig_GO_down <- GO_down[GO_down$padj < 0.05, ]

sig_GO_down %>% 
  top_n(10, wt = -padj) %>% 
  mutate(hitsPerc = numDEInCat * 100 / numInCat) %>% 
  ggplot(aes(x = hitsPerc, 
             y = reorder(term, hitsPerc), 
             colour = padj, 
             size = numDEInCat)) +
  geom_point() +
  scale_color_gradient(low = "lightblue", high = "blue") +
  expand_limits(x = 0) +
  labs(title = "GO enrichment (downregulated genen in RA)",
       x = "Percentage DE genen (%)", 
       y = "GO term (biologische processen", 
       colour = "Adjusted p-value", 
       size = "Aantal DE genen")


# Top 10 processen in een staafdiagram visualiseren
library(GO.db)

top10_GO_up <- GO_up %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  slice_head(n = 10)


top10_GO_down <- GO_down %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  slice_head(n = 10)


## Data combineren
### Richting toevoegen
top10_GO_up$direction <- "Upregulated"
top10_GO_down$direction <- "Downregulated"

### Combineer datasets
top10_GO_combined <- rbind(top10_GO_up, top10_GO_down)


## Plotten
ggplot(top10_GO_combined,
       aes(x = reorder(term, -log10(padj)),
           y = -log10(padj),
           fill = direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ direction, scales = "free") +
  scale_fill_manual(values = c(
    "Upregulated" = "firebrick",
    "Downregulated" = "steelblue"
  )) +
  labs(
    title = "GO enrichment in RA",
    x = "GO-term",
    y = "-log10(adjusted p-value)"
  ) +
  theme_light()
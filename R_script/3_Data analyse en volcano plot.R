library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)

# Set working directory
setwd(...)


# Metadata:
## Samplenaam linken met conditie
conditie <- c("gezond", "gezond", "gezond", "gezond", 
              "reuma", "reuma", "reuma", "reuma")
conditie_table <- data.frame(conditie)
rownames(conditie_table) <- c("SRR4785819", "SRR4785820", "SRR4785828", "SRR4785831", 
                              "SRR4785979", "SRR4785980", "SRR4785986", "SRR4785988")

# Voor elk gen kijken: hoeveel expressie in RA en in een gezond persoon:
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = conditie_table,
                              design = ~ conditie)

# Voer analyse uit
dds <- DESeq(dds)
resultaten <- results(dds)

write.table(resultaten, file = 'ResultatenReuma.csv', row.names = TRUE, col.names = TRUE)

# Resultaten samenvatten
## Van hoeveel genen stijgt of daalt de expressie significant
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange > 1, na.rm = TRUE)   #2085
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange < -1, na.rm = TRUE)  #2487


# Volcano plot:
EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj')
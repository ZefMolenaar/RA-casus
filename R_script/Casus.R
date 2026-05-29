library(Rsubread)
library(BiocManager)
setwd("~/HBO biotechnologie/leerjaar 2/periode 4/Casus")

buildindex(
  basename = 'ref_mens',
  reference = 'GCF_000001405.40_GRCh38.p14_genomic.fna',
  memory = 6000,
  indexSplit = TRUE)

align.normaal1 <- align(index = "ref_mens", readfile1 = "SRR4785819_1_subset40k.fastq", readfile2 = "SRR4785819_2_subset40k.fastq", output_file = "normaal1.BAM")
align.normaal2 <- align(index = "ref_mens", readfile1 = "SRR4785820_1_subset40k.fastq", readfile2 = "SRR4785820_2_subset40k.fastq", output_file = "normaal2.BAM")
align.nomraal3 <- align(index = "ref_mens", readfile1 = "SRR4785828_1_subset40k.fastq", readfile2 = "SRR4785828_2_subset40k.fastq", output_file = "normaal3.BAM")
align.normaal4 <- align(index = "ref_mens", readfile1 = "SRR4785831_1_subset40k.fastq", readfile2 = "SRR4785831_2_subset40k.fastq", output_file = "normaal4.BAM")
align.RA1 <- align(index = "ref_mens", readfile1 = "SRR4785979_1_subset40k.fastq", readfile2 = "SRR4785979_2_subset40k.fastq", output_file = "RA1.BAM")
align.RA2 <- align(index = "ref_mens", readfile1 = "SRR4785980_1_subset40k.fastq", readfile2 = "SRR4785980_2_subset40k.fastq", output_file = "RA2.BAM")
align.RA3 <- align(index = "ref_mens", readfile1 = "SRR4785986_1_subset40k.fastq", readfile2 = "SRR4785986_2_subset40k.fastq", output_file = "RA3.BAM")
align.RA4 <- align(index = "ref_mens", readfile1 = "SRR4785988_1_subset40k.fastq", readfile2 = "SRR4785988_2_subset40k.fastq", output_file = "RA4.BAM")

library(Rsamtools)

samples_casus <- c('normaal1', 'normaal2', 'normaal3', 'normaal4', 'RA1', 'RA2', 'RA3', 'RA4')

lapply(samples_casus, function(s) {sortBam(file = paste0(s, '.BAM'), destination = paste0(s, '.sorted'))
})
lapply(samples_casus, function(s) {indexBam(file = paste0(s, '.sorted.bam'))
})

Zefsamples_casus<-c("normaal1.BAM", "normaal2.BAM", "normaal3.BAM", "normaal4.BAM", "RA1.BAM", "RA2.BAM", "RA3.BAM", "RA4.BAM")

count_matrix <- featureCounts(
  files = Zefsamples_casus,
  annot.ext = "genomic.gtf",
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE, 
  GTF.featureType = "gene",
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE
)

str(count_matrix)
head(count_matrix$annotation)
counts <- count_matrix$counts
head(counts)
colnames(counts) <- c('normaal1', 'normaal2', 'normaal3', 'normaal4', 'RA1', 'RA2', 'RA3', 'RA4')
head(counts)
write.csv(counts, "Test_RA_countmatrix.csv")

library(DESeq2)
library(KEGGREST)
library(EnhancedVolcano)
library(pathview)

counts <- read.table("count_matrix.txt", row.names = 1)
colnames(counts)<- c('normaal1', 'normaal2', 'normaal3', 'normaal4', 'RA1', 'RA2', 'RA3', 'RA4')

treatment <- c("normaal", "normaal", "normaal", "normaal", "RA", "RA", "RA", "RA")
treatment_table <- data.frame(treatment)
rownames(treatment_table) <- c('normaal1', 'normaal2', 'normaal3', 'normaal4', 'RA1', 'RA2', 'RA3', 'RA4')
head(treatment_table)


dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = treatment_table,
                              design = ~ treatment)
dds <- DESeq(dds)
resultaten <- results(dds)

write.table(resultaten, file = 'resultatenRA.csv', row.names = TRUE, col.names = TRUE)
resultaten <- read.table("resultatenRA.csv")

sum(resultaten$padj < 0.05 & resultaten$log2FoldChange > 1, na.rm = TRUE)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange < -1, na.rm = TRUE)

hoogste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = TRUE), ]
laagste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = FALSE), ]
laagste_p_waarde <- resultaten[order(resultaten$padj, decreasing = FALSE), ]

EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = "log2FoldChange",
                y = "padj")

dev.copy(png, 'VolcanoplotCASUS.png', 
         width = 8,
         height = 10,
         units = 'in',
         res = 500)
dev.off()

resultaten[1] <- NULL
resultaten[2:5] <- NULL

pathview(
  gene.data = resultaten,
  pathway.id = "hsa05323",  
  species = "hsa",          
  gene.idtype = "KEGG",     
  limit = list(gene = 5)    
)

BiocManager::install("goseq")

library("goseq")
library("geneLenDataBase")


all <- rownames(resultaten)
library(dplyr)
deg <- resultaten %>%
  filter(padj < 0.05)
deg <- rownames(deg)



gene.vector=as.integer(all%in%deg)

options(Ncpus = 8)
BiocManager::install('goseq')
BiocManager::install('geneLenDataBase')
library(goseq)
library(geneLenDataBase)
library(tidyverse)

supportedOrganisms() %>% filter(str_detect(Genome, "hg19"))
sigData <- as.integer(!is.na(resultaten$padj) & resultaten$padj < 0.01)
names(sigData) <- rownames(resultaten)

head(sigData)

pwf <- nullp(sigData, "hg19", "geneSymbol", bias.data = resultaten$padj)
goResults <- goseq(pwf, "hg19","geneSymbol", test.cats=c("GO:BP"))

install.packages('GOplot')
library(GOplot)

goResults %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")

if (!requireNamespace("pathview", quietly = TRUE)) {
  BiocManager::install("pathview")
}
library(pathview)

pathview(
  gene.data = counts,
  pathway.id = "hsa05323",  
  species = "hsa",          
  gene.idtype = "SYMBOL",     
  limit = list(gene = 5)    
)
keggLink("pathway", "hsa:05323")
keggLink("hsa", "hsa:05323")

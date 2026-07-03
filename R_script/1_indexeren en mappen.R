library(Rsubread)
library(Rsamtools)

# Set working directory
setwd(...)

# Referentiegenoom indexeren:
buildindex(
  basename = 'ref_humangenome',
  reference = 'Genoom_Mens.fna',
  memory = 6000,
  indexSplit = TRUE)  


# Mapping (alle monsters alignen op het referentiegenoom):
align.normaal1 <- align(index = "ref_humangenome", 
                        readfile1 = "SRR4785819_1_subset40k.fastq", 
                        readfile2 = "SRR4785819_2_subset40k.fastq",
                        output_file = "normaal1.BAM")
align.normaal2 <- align(index = "ref_humangenome", 
                        readfile1 = "SRR4785820_1_subset40k.fastq", 
                        readfile2 = "SRR4785820_2_subset40k.fastq",
                        output_file = "normaal2.BAM")
align.normaal3 <- align(index = "ref_humangenome", 
                        readfile1 = "SRR4785828_1_subset40k.fastq", 
                        readfile2 = "SRR4785828_2_subset40k.fastq",
                        output_file = "normaal3.BAM")
align.normaal4 <- align(index = "ref_humangenome", 
                        readfile1 = "SRR4785831_1_subset40k.fastq", 
                        readfile2 = "SRR4785831_2_subset40k.fastq",
                        output_file = "normaal4.BAM")

align.reuma1 <- align(index = "ref_humangenome", 
                      readfile1 = "SRR4785979_1_subset40k.fastq", 
                      readfile2 = "SRR4785979_2_subset40k.fastq",
                      output_file = "reuma1.BAM")
align.reuma2 <- align(index = "ref_humangenome", 
                      readfile1 = "SRR4785980_1_subset40k.fastq", 
                      readfile2 = "SRR4785980_2_subset40k.fastq",
                      output_file = "reuma2.BAM")
align.reuma3 <- align(index = "ref_humangenome", 
                      readfile1 = "SRR4785986_1_subset40k.fastq", 
                      readfile2 = "SRR4785986_2_subset40k.fastq",
                      output_file = "reuma3.BAM")
align.reuma4 <- align(index = "ref_humangenome", 
                      readfile1 = "SRR4785988_1_subset40k.fastq", 
                      readfile2 = "SRR4785988_2_subset40k.fastq",
                      output_file = "reuma4.BAM")



# Bestandsnamen van de monsters
samples <- c('normaal1', 'normaal2', 'normaal3', 'normaal4', 'reuma1', 'reuma2',
             'reuma3', 'reuma4')

# Voor elk monster: sorteer en indexeer de BAM-file
## Sorteer BAM-bestanden
lapply(samples, function(s) {sortBam(file = paste0(s, '.BAM'), 
                                     destination = paste0(s, '.sorted'))
})

## Indexeer de gesorteerde BAM-file
lapply(samples, function(s) {indexBam(file = paste0(s, '.sorted.bam'))})
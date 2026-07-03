library(Rsubread)

# Set working directory
setwd(...)

# Count matrix opstellen
allsamples <- c("normaal1.BAM", "normaal2.BAM", "normaal3.BAM", "normaal4.BAM", 
                "reuma1.BAM", "reuma2.BAM", "reuma3.bam", "reuma4.bam")
count_matrix <- featureCounts(
  files = allsamples,
  annot.ext = "genomic.gtf",
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE,
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE
)

# Count matrix bekijken
str(count_matrix)


# Specifieke kolom pakken
counts <- count_matrix$counts
head(counts)

colnames(counts) <- c("gezond1", "gezond2", "gezond3", "gezond4", 
                      "reuma1", "reuma2", "reuma3", "reuma4")
head(counts)

# Opslaan:

write.csv(counts, 'reuma_countmatrix.csv')


# Count matrix verkregen van docent inladen:

counts = count_matrix_RA
# Set a default CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install DESeq2 and edgeR using BiocManager
BiocManager::install("DESeq2")
BiocManager::install("edgeR")
BiocManager::install("limma")
BiocManager::install("AnnotationDbi")
BiocManager::install("GenomicRanges")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("clusterProfiler")



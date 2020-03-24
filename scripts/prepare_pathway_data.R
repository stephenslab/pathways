# TO DO: Explain here what this script does, and how to use it.
library(Matrix)
library(readr)
source("../code/read_data.R")

# LOAD DATA
# ---------
# Read data from gene_info file.
cat("Reading gene data from Homo_sapiens.gene_info.gz.\n")
gene_info <- read_gene_info("../data/Homo_sapiens.gene_info.gz")

# Read and process the BioSystems pathway data.
cat("Reading BioSystems data from bsid2info.gz and biosystems_gene.gz.\n")
bsid2info        <- read_bsid2info("../data/bsid2info.gz")
biosys_gene_sets <- read_biosystems_gene_sets("../data/biosystems_gene.gz",
                                              bsid2info,gene_info)

# Read and process the Pathway Commons pathway data.
cat("Reading Pathway Commons data from PathwayCommons12.All.hgnc.gmt.gz.\n")
out <- read_pathway_commons_data("../data/PathwayCommons12.All.hgnc.gmt.gz",
                                 gene_info)
pc_pathways  <- out$pathways
pc_gene_sets <- out$gene_sets
rm(out)

# Combine the BioSystems and Pathway Commons pathway and gene set data.
# TO DO.

# Remove any pathways that do not have any genes (or at least any
# recognizable gene symbols).
# TO DO.
    
# SUMMARIZE DATA
# --------------
print(nrow(genes))
print(summary(pathways$datasource))

# WRITE DATA TO FILE
# ------------------
# TO DO.

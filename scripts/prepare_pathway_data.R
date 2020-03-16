# TO DO: Explain here what this script does, and how to use it.
library(readr)
source("../code/read_data.R")

# LOAD DATA
# ---------
# Read the HGNC (HUGO Gene Nomenclature Committee) gene data from the
# tab-delimited text file.
cat("Reading HGNC gene data from genes.txt.gz.\n")
genes <- read_hgnc_data("../data/genes.txt.gz")

# Read the Pathway Commons pathway meta data (e.g., pathway names,
# data sources) from the tab-delimited text file.
cat("Reading Pathway Commons pathway data from pathways.txt.gz.\n")
pathways <- read_pc_pathway_data("../data/pathways.txt.gz")

# SUMMARIZE DATA
# --------------
nrow(genes)
print(summary(pathways$datasource))

# WRITE DATA TO FILE
# ------------------
# TO DO.

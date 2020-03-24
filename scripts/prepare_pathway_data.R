# Script to digest the BioSystems and Pathway Commons pathway data for
# use in gene-set enrichment analyses.
library(Matrix)
library(readr)
library(tools)
source("../code/read_data.R")

# LOAD DATA
# ---------
# Read data from gene_info file.
cat("Reading gene data from Homo_sapiens.gene_info.gz.\n")
gene_info <- read_gene_info("../data/Homo_sapiens.gene_info.gz")

# Read and process the BioSystems pathway data.
cat("Reading BioSystems data from bsid2info.gz and biosystems_gene.gz.\n")
bsid2info    <- read_bsid2info("../data/bsid2info.gz")
bs_gene_sets <- read_biosystems_gene_sets("../data/biosystems_gene.gz",
                                          bsid2info,gene_info)

# Read and process the Pathway Commons pathway data.
cat("Reading Pathway Commons data from PathwayCommons12.All.hgnc.gmt.gz.\n")
out <- read_pathway_commons_data("../data/PathwayCommons12.All.hgnc.gmt.gz",
                                 gene_info)
pc_pathways  <- out$pathways
pc_gene_sets <- out$gene_sets
rm(out)

# Combine the BioSystems and Pathway Commons pathway and gene set data.
bsid2info   <- cbind(bsid2info,data.frame(database = "BioSystems"))
pc_pathways <- cbind(pc_pathways,data.frame(database = "PC"))
names(bsid2info)[1]   <- "id"
names(pc_pathways)[3] <- "id"
levels(pc_pathways$data_source) <-
  c("humancyc","inoh","KEGG","netpath","panther","pathbank",
    "Pathway Interaction Database","REACTOME")
pathways  <- merge(bsid2info,pc_pathways,all = TRUE,sort = FALSE,
                   by = c("name","id","data_source","database"))
gene_sets <- cbind(bs_gene_sets,pc_gene_sets)
    
# Remove any pathways that do not have any genes (or at least any
# recognizable gene symbols).
i         <- which(colSums(gene_sets > 0) > 0)
pathways  <- pathways[i,]
gene_sets <- gene_sets[,i]
    
# SUMMARIZE DATA
# --------------
cat("Number of genes:                    ",nrow(gene_info),"\n")
cat("Number of BioSystems pathways:      ",
    sum(pathways$database=="BioSystems"),"\n")
cat("Number of Pathway Commons pathways: ",sum(pathways$database == "PC"),"\n")
cat("Total number of pathways:           ",nrow(pathways),"\n")
cat("Number of gene-pathway associations:",nnzero(gene_sets),"\n")

# SAVE PROCESSED DATA
# -------------------
cat("Saving gene and pathway data to pathways.RData.\n")
save(list = c("gene_info","pathways","gene_sets"),file = "pathway.RData")
resaveRdaFiles("pathway.RData")

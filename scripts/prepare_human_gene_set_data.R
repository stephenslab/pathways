# Script to digest the BioSystems, Pathway Commons and MSigDB human
# gene set data for use in gene-set enrichment analyses.
library(Matrix)
library(readr)
library(tools)
library(xml2)
library(msigdbr)
source("../code/read_gene_set_data.R")

# LOAD DATA
# ---------
# Read gene data from the gene_info files.
cat("Reading gene data from Homo_sapiens.gene_info.gz.\n")
gene_info <- read_gene_info("../data/Homo_sapiens.gene_info.gz")

# Read and process the BioSystems pathway data.
cat("Reading BioSystems data from bsid2info.gz and biosystems_gene.gz.\n")
bsid2info    <- read_bsid2info("../data/bsid2info.gz",organism = 9606)
bs_gene_sets <- read_biosystems_gene_sets("../data/biosystems_gene.gz",
                                          bsid2info,gene_info)

# Read and process the Pathway Commons pathway data.
cat("Reading Pathway Commons data from PathwayCommons12.All.hgnc.gmt.gz.\n")
out <- read_pathway_commons_data("../data/PathwayCommons12.All.hgnc.gmt.gz",
                                 gene_info)
pc_pathways  <- out$pathways
pc_gene_sets <- out$gene_sets
rm(out)

# Read and process MSigDB gene set data.
cat("Reading MSigDB gene set data from msigdb_v7.1.xml, and\n")
cat("extracting MSigDB gene sets using msigdbr package.\n")
out <- get_msigdb_gene_sets("../data/msigdb_v7.1.xml",gene_info,"Homo sapiens")
msigdb_info      <- out$info
msigdb_gene_sets <- out$gene_sets
rm(out)

stop()

# Combine the BioSystems, Pathway Commons and MSigDB gene set data.
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
    
# Remove pathways that do not have any genes (or at least any
# recognizable gene symbols).
i         <- which(colSums(gene_sets > 0) > 0)
pathways  <- pathways[i,]
gene_sets <- gene_sets[,i]

# TO DO: Set all nonzeros in gene_sets to be 1.

# SUMMARIZE DATA
# --------------
cat("genes:                    ",nrow(gene_info),"\n")
cat("BioSystems pathways:      ",
    sum(pathways$database=="BioSystems"),"\n")
cat("Pathway Commons: ",sum(pathways$database == "PC"),"\n")
cat("MSigDB: ")
cat("Total number of pathways:           ",nrow(pathways),"\n")

# SAVE PROCESSED DATA
# -------------------
cat("Saving gene set data to gene_sets_human.RData.\n")
save(list = c("gene_info","pathways","gene_sets"),file = "gene_sets_human.RData")
resaveRdaFiles("gene_sets_human.RData")

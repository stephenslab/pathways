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
cat("Reading MSigDB gene set data from msigdb_v7.2.xml, and\n")
cat("extracting MSigDB gene sets using msigdbr package.\n")
out <- get_msigdb_gene_sets("../data/msigdb_v7.2.xml",gene_info,"Homo sapiens")
msigdb_info      <- out$info
msigdb_gene_sets <- out$gene_sets
rm(out)

# Combine BioSystems, Pathway Commons and MSigDB gene set data.
cat("Merging gene set data.\n")
bsid2info   <- cbind(bsid2info,data.frame(database = "BioSystems"))
pc_pathways <- cbind(pc_pathways,data.frame(database = "PC"))
msigdb_info <- cbind(msigdb_info,data.frame(data_source = as.character(NA),
                                            database = "MSigDB"))
names(bsid2info)[1]        <- "id"
names(pc_pathways)[3]      <- "id"
names(msigdb_info)[c(1,2)] <- c("name","id")
levels(bsid2info$data_source) <- tolower(levels(bsid2info$data_source))
levels(bsid2info$data_source)[3] <- "pid"
gene_set_info <- merge(bsid2info,pc_pathways,all = TRUE,sort = FALSE,
                       by = c("name","id","data_source","database"))
gene_set_info <- merge(gene_set_info,msigdb_info,all = TRUE,sort = FALSE,
                       by = c("name","id","data_source","database"))
gene_sets <- cbind(bs_gene_sets,pc_gene_sets,msigdb_gene_sets)

# Set all nonzeros in the "gene_sets" matrix to be 1.
gene_sets <- as(gene_sets > 0,"dgCMatrix")

# Remove pathways that do not have any genes.
i             <- which(colSums(gene_sets) > 0)
gene_set_info <- gene_set_info[i,]
gene_sets     <- gene_sets[,i]

# SUMMARIZE DATA
# --------------
cat("genes:               ",nrow(gene_info),"\n")
cat("BioSystems gene sets:",sum(gene_set_info$database == "BioSystems"),"\n")
cat("PC gene sets:        ",sum(gene_set_info$database == "PC"),"\n")
cat("MSigDB gene sets:    ",sum(gene_set_info$database == "MSigDB"),"\n")
cat("Total gene sets:     ",nrow(gene_set_info),"\n")

# SAVE PROCESSED DATA
# -------------------
cat("Saving gene and gene-set data to gene_sets_human.RData.\n")
save(list = c("gene_info","gene_set_info","gene_sets"),
     file = "gene_sets_human.RData")
resaveRdaFiles("gene_sets_human.RData")

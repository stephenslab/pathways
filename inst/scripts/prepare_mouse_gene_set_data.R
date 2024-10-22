# Script to digest the BioSystems and MSigDB mouse gene set data for
# use in gene-set enrichment analyses.
library(Matrix)
library(readr)
library(tools)
library(xml2)
library(msigdbr)
source("../code/read_gene_set_data.R")

# LOAD DATA
# ---------
# Read gene data from the gene_info files. (The "HGNC" column is
# removed because it is only relevant to human genes.)
cat("Reading mouse gene data from Mus_musculus.gene_info.gz.\n")
gene_info <- read_gene_info("../datafiles/Mus_musculus.gene_info.gz")
gene_info <- gene_info[-6]

# Read and process the BioSystems pathway data.
cat("Reading BioSystems data from bsid2info.gz and biosystems_gene.gz.\n")
bsid2info    <- read_bsid2info("../datafiles/bsid2info.gz",organism = 10090)
bs_gene_sets <- read_biosystems_gene_sets("../datafiles/biosystems_gene.gz",
                                          bsid2info,gene_info)

# Read and process MSigDB gene set data.
cat("Reading MSigDB gene set data from msigdb_v7.2.xml, and\n")
cat("extracting MSigDB gene sets using msigdbr package.\n")
out <- get_msigdb_gene_sets("../datafiles/msigdb_v7.2.xml",gene_info,
                            "Mus musculus")
msigdb_info      <- out$info
msigdb_gene_sets <- out$gene_sets
rm(out)

# Combine BioSystems and MSigDB gene set data.
cat("Merging gene set data.\n")
bsid2info   <- cbind(bsid2info,data.frame(database = "BioSystems"))
msigdb_info <- cbind(msigdb_info,data.frame(data_source = as.character(NA),
                                            database = "MSigDB"))
names(bsid2info)[1]        <- "id"
names(msigdb_info)[c(1,2)] <- c("name","id")
levels(bsid2info$data_source) <- tolower(levels(bsid2info$data_source))
levels(bsid2info$data_source)[3] <- "pid"
gene_set_info <- merge(bsid2info,msigdb_info,all = TRUE,sort = FALSE,
                       by = c("name","id","data_source","database"))
gene_sets <- cbind(bs_gene_sets,msigdb_gene_sets)

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
cat("MSigDB gene sets:    ",sum(gene_set_info$database == "MSigDB"),"\n")
cat("Total gene sets:     ",nrow(gene_set_info),"\n")

# SAVE PROCESSED DATA
# -------------------
cat("Saving gene and gene-set data to gene_sets_mouse.RData.\n")
save(list = c("gene_info","gene_set_info","gene_sets"),
     file = "gene_sets_mouse.RData")
resaveRdaFiles("gene_sets_mouse.RData")

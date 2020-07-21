# Script to digest the BioSystems, Pathway Commons and MSigDB gene set
# data for use in gene-set enrichment analyses.
library(Matrix)
library(readr)
library(tools)
library(xml2)
library(msigdbr)
source("../code/read_gene_set_data.R")

# LOAD DATA
# ---------
# Read human gene data from the gene_info files.
cat("Reading human gene data from Homo_sapiens.gene_info.gz.\n")
gene_info_human <- read_gene_info("../data/Homo_sapiens.gene_info.gz")

# Read mouse gene data from the gene_info files.
cat("Reading mouse gene data from Mus_musculus.gene_info.gz.\n")
gene_info_mouse <- read_gene_info("../data/Mus_musculus.gene_info.gz")

# Combine the gene data into one master table.
gene_info <- rbind(gene_info_human,gene_info_mouse)
rm(gene_info_human,gene_info_mouse)

# Read and process the BioSystems pathway data.
cat("Reading BioSystems data from bsid2info.gz and biosystems_gene.gz.\n")
bsid2info    <- read_bsid2info("../data/bsid2info.gz")
bsid2info    <- subset(bsid2info,tax_id == 9606 | tax_id == 10090)
bsid2info    <- transform(bsid2info,tax_id = factor(tax_id))
bs_gene_sets <- read_biosystems_gene_sets("../data/biosystems_gene.gz",
                                          bsid2info,gene_info)

# Read and process the Pathway Commons (human) pathway data.
cat("Reading Pathway Commons data from PathwayCommons12.All.hgnc.gmt.gz.\n")
out <- read_pathway_commons_data("../data/PathwayCommons12.All.hgnc.gmt.gz",
                                 gene_info)
pc_pathways  <- out$pathways
pc_gene_sets <- out$gene_sets
rm(out)

# Read and process MSigDB gene set data.
cat("Extracting MSigDB gene sets using msigdbr package.\n")
msigdb_gene_sets_human <- get_msigdb_gene_sets(gene_info,"Homo sapiens")
msigdb_gene_sets_mouse <- get_msigdb_gene_sets(gene_info,"Mus musculus")
msigdb_gene_sets       <- cbind(msigdb_gene_sets_human,msigdb_gene_sets_mouse)
cat("Reading MSigDB gene set data from msigdb_v7.1.xml.\n")
msigdb_info <- read_msigdb_xml("../data/msigdb_v7.1.xml")

# Align the rows of msigdb_info with the columns of msigdb_gene_sets.
msigdb_info <- subset(msigdb_info,
                      is.element(systematic_name,colnames(msigdb_gene_sets)))

cols <- match(msigdb_info$systematic_name,colnames(msigdb_gene_sets))
msigdb_gene_sets <- msigdb_gene_sets[,cols]

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
    
# Remove any pathways that do not have any genes (or at least any
# recognizable gene symbols).
i         <- which(colSums(gene_sets > 0) > 0)
pathways  <- pathways[i,]
gene_sets <- gene_sets[,i]

# TO DO: Set all nonzeros to be 1.

# TO DO: Remove *duplicate* gene sets.

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
cat("Saving gene set data to gene_sets.RData.\n")
save(list = c("gene_info","pathways","gene_sets"),file = "gene_sets.RData")
resaveRdaFiles("gene_sets.RData")

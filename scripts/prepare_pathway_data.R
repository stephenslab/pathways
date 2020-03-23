# TO DO: Explain here what this script does, and how to use it.
library(Matrix)
library(readr)
source("../code/read_data.R")

# LOAD DATA
# ---------
# Read data from gene_info file.
cat("Reading gene data from Homo_sapiens.gene_info.gz.\n")
gene_info <- read_gene_info("../data/Homo_sapiens.gene_info.gz")

# TO DO: Explain here what this code does.
cat("Reading BioSystems data from bsid2info.gz and biosystems_gene.gz.\n")
biosys <- read_bsid2info("../data/bsid2info.gz")

stop()

# Read data from biosystems_gene.
biosys_gene <- read_delim("../data/biosystems_gene.gz",delim = "\t",
                          col_names = c("bsid","geneid","score"),
                          col_types = cols("i","i","i"),progress = FALSE)
class(biosys_gene) <- "data.frame"
biosys_gene <- subset(biosys_gene,
                      is.element(bsid,biosys$bsid)  &
                      is.element(geneid,gene_info$GeneID))

# Create an n x m sparse binary matrix, where n is the number of genes
# and m is the number of pathways.
biosys_gene_sets <-
  sparseMatrix(i = match(biosys_gene$geneid,gene_info$GeneID),
               j = match(biosys_gene$bsid,biosys$bsid),
               x = biosys_gene$score,dims = c(nrow(gene_info),nrow(biosys)))
rownames(biosys_gene_sets) <- gene_info$GeneID
colnames(biosys_gene_sets) <- biosys$bsid

# Remove gene sets with no genes.
cols             <- which(colSums(biosys_gene_sets) > 0)
biosys_gene_sets <- biosys_gene_sets[,cols]

stop()

# Read the HGNC (HUGO Gene Nomenclature Committee) gene data from the
# tab-delimited text file.
cat("Reading HGNC data from genes.txt.gz.\n")
genes <- read_hgnc_data("../data/genes.txt.gz")


stop()

# Read the Pathway Commons pathway meta data (e.g., pathway names,
# data sources) from the tab-delimited text file.
cat("Reading Pathway Commons pathway data from pathways.txt.gz.\n")
# pathways <- read_pc_pathway_data("../data/pathways.txt.gz")

# Read the gene set data.
dat       <- readLines("../data/PathwayCommons12.All.hgnc.gmt.gz")
n         <- length(dat)
m         <- nrow(genes)
pathways  <- data.frame(name       = rep("",n),
                        datasource = rep("",n),
                        url        = rep("",n),
                        stringsAsFactors = FALSE)
i <- NULL
j <- NULL
for (k in 1:n) {
  x <- unlist(strsplit(dat[[k]],"\t",fixed = TRUE))
  y <- trimws(unlist(strsplit(unlist(strsplit(x[2],";",fixed = TRUE)),":",
                              fixed = TRUE)))
  pathways[k,"name"]       <- y[which(y == "name") + 1]
  pathways[k,"datasource"] <- y[which(y == "datasource") + 1]
  pathways[k,"url"]        <- x[1]
  x  <- x[-(1:2)]
  ik <- match(x,genes$approved_symbol)
  if (sum(is.na(ik) > 0))
    cat(sprintf("%d:%d/%d\n",k,sum(is.na(ik) > 0),length(ik)))
  ik <- ik[!is.na(ik)]
  i  <- c(i,ik)
  j  <- c(j,rep(k,length(ik)))
}
pathways <- transform(pathways,datasource = factor(datasource))
gene_sets <- sparseMatrix(i,j,x = 1,dims = c(m,n))
    
# SUMMARIZE DATA
# --------------
print(nrow(genes))
print(summary(pathways$datasource))

# WRITE DATA TO FILE
# ------------------
# TO DO.

# TO DO: Explain here what this script does, and how to use it.
library(Matrix)
library(readr)
source("../code/read_data.R")

# LOAD DATA
# ---------
# Read the HGNC (HUGO Gene Nomenclature Committee) gene data from the
# tab-delimited text file.
cat("Reading HGNC gene data from genes.txt.gz.\n")
genes <- read_hgnc_data("../data/genes.txt.gz")

# TO DO: Explain here what this code does.
cat("Reading NCBI BioSystems pathway data from bsid2info.gz and",
    "biosystems_gene.gz.\n")
suppressWarnings(
  biosys <- read_delim("../data/bsid2info.gz",delim = "\t",quote = "",
                       col_names = c("bsid","datasource","accession","name",
                           "type","scope","taxid","description"),
                       col_types = cols("i","c","c","c","c","c","i","c")))
class(biosys) <- "data.frame"
biosys <- subset(biosys,type == "pathway" & taxid == 9606)
biosys <- biosys[c("bsid","datasource","accession","name","description")]
biosys <- transform(biosys,datasource = factor(datasource))

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

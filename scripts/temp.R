# An initial attempt at reading in the x
library(Matrix)
library(xml2)
library(msigdbr)

# Get the human and mouse gene-set annotations.
msigdb_gene_human <- msigdbr(species = "Homo sapiens")
msigdb_gene_mouse <- msigdbr(species = "Mus musculus")
class(msigdb_gene_human) <- "data.frame"
class(msigdb_gene_mouse) <- "data.frame"
msigdb <- rbind(msigdb_gene_human,msigdb_gene_mouse)

# Remove annotations that do not correspond to an entry in gene
# database.
print(nrow(msigdb))
msigdb <- subset(msigdb,is.element(entrez_gene,gene_info$GeneID))
print(nrow(msigdb))

gsid   <- sort(unique(msigdb$gs_id))
n      <- nrow(msigdb)

# Create the n x m sparse (binary) adjacency matrix.
out <- sparseMatrix(i = match(msigdb$entrez_gene,gene_info$GeneID),
                    j = match(msigdb$gs_id,gsid),x = rep(1,n),
                    dims = c(nrow(gene_info),length(gsid)))
rownames(out) <- gene_info$GeneID
colnames(out) <- gsid

dat <- readLines("../data/msigdb.v7.1.entrez.gmt.gz")
n   <- length(dat)
for (t in 1:n) {
  x <- unlist(strsplit(dat[[t]],"\t",fixed = TRUE))
}
# out <- read_xml("../data/msigdb_v7.1.xml")

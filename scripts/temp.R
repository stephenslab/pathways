# An initial attempt at reading in the msigdb gene sets.
library(Matrix)
library(xml2)
library(msigdbr)

# Retrieve the MSigDB gene set meta-data.
out <- read_xml("../data/msigdb_v7.1.xml")
out <- as_list(out)[[1]]
extract_attribute <- function (dat, attribute)
  sapply(dat,function (x) attr(x,attribute))
msigdb_info <-
  data.frame(standard_name     = extract_attribute(out,"STANDARD_NAME"),
             systematic_name   = extract_attribute(out,"SYSTEMATIC_NAME"),
             category_code     = extract_attribute(out,"CATEGORY_CODE"),
             sub_category_code = extract_attribute(out,"SUB_CATEGORY_CODE"),
             organism          = extract_attribute(out,"ORGANISM"),
             description_brief = extract_attribute(out,"DESCRIPTION_BRIEF"),
             stringsAsFactors = FALSE)
msigdb_info <- transform(msigdb_info,
                         category_code     = factor(category_code),
                         sub_category_code = factor(sub_category_code),
                         organism          = factor(organism))

# Get the human and mouse gene-set annotations.
msigdb_gene_human <- msigdbr(species = "Homo sapiens")
msigdb_gene_mouse <- msigdbr(species = "Mus musculus")
class(msigdb_gene_human) <- "data.frame"
class(msigdb_gene_mouse) <- "data.frame"
msigdb_gene <- rbind(msigdb_gene_human,msigdb_gene_mouse)

# Remove annotations that do not correspond to an entry in gene
# database.
msigdb_gene <- subset(msigdb_gene,is.element(entrez_gene,gene_info$GeneID))
gsid        <- sort(unique(msigdb_gene$gs_id))
n           <- nrow(msigdb_gene)

# Create the n x m sparse (binary) adjacency matrix.
msigdb_gene_sets <-
  sparseMatrix(i = match(msigdb_gene$entrez_gene,gene_info$GeneID),
               j = match(msigdb_gene$gs_id,gsid),
               x = rep(1,n),dims = c(nrow(gene_info),length(gsid)))
rownames(msigdb_gene_sets) <- gene_info$GeneID
colnames(msigdb_gene_sets) <- gsid

# Remove msigdb_info records that do not have an annotation.
msigdb_info <- subset(msigdb_info,
                      is.element(systematic_name,colnames(msigdb_gene_sets)))

# Align the rows of msigdb_info with the columns of msigdb_gene_sets.
cols <- match(msigdb_info$systematic_name,colnames(msigdb_gene_sets))
msigdb_gene_sets <- msigdb_gene_sets[,cols]

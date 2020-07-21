# Read the "gene_info" tab-delimited text file downloaded from the
# NCBI FTP site (ftp.ncbi.nih.gov/gene). The return value is a data
# frame with one row per gene, and the following columns: tax_id,
# GeneID, Symbol, Synonyms, chromosome, Ensembl and HGNC.
read_gene_info <- function (file) {

  # Read the data into a data frame.
  out <- suppressMessages(read_delim(file,delim = "\t",col_names = TRUE))
  class(out) <- "data.frame"
  dbXrefs    <- out$dbXrefs
  out        <- out[c("#tax_id","GeneID","Symbol","Synonyms","chromosome")]
  names(out)[1] <- "tax_id"
  
  # Set any entries with a single hyphen to NA, and convert the
  # "chromosome" column to a factor.
  out$chromosome[out$chromosome == "-"] <- NA
  out$Synonyms[out$Synonyms == "-"]     <- NA
  dbXrefs[dbXrefs == "-"]               <- NA
  out <- transform(out,chromosome = factor(chromosome))

  # Extract the Ensembl ids. Note that a small number of genes map to
  # more than one Ensembl id; in those cases, we retain the first
  # Ensembl id only.
  dbXrefs <- strsplit(dbXrefs,"|",fixed = TRUE)
  out$Ensembl <- sapply(dbXrefs,function (x) {
                            i <- which(substr(x,1,8) == "Ensembl:")
                            if (length(i) > 0)
                              return(substr(x[i[1]],9,nchar(x[i[1]])))
                            else
                                return(as.character(NA))
                          })

  # For human genes, extract the HGNC (HUGO Gene Nomenclature
  # Committee) ids.
  out$HGNC <- sapply(dbXrefs,function (x) {
                i <- which(substr(x,1,10) == "HGNC:HGNC:")
                if (length(i) > 0)
                  return(substr(x[i[1]],6,nchar(x[i[1]])))
                else
                  return(as.character(NA))
              })

  # Return the processed gene data.
  return(out)
}

# Read the "bsid2info" tab-delimited file downloaded from the NCBI FTP
# site (ftp.ncbi.nih.gov/pub/biosystems). The return value is a data
# frame with one row per BioSystems pathway, and the following
# columns: bsid, tax_id, data_source, accession and name.
read_bsid2info <- function (file, organism = 9606) {
  out <- suppressWarnings(
    read_delim(file,delim = "\t",quote = "",progress = FALSE,
               col_names = c("bsid","data_source","accession","name",
                             "type","scope","tax_id","description"),
               col_types = cols("i","c","c","c","c","c","i","c")))
  class(out) <- "data.frame"
  out <- subset(out,type == "pathway" & tax_id == organism)
  out <- out[c("bsid","data_source","accession","name")]
  out <- transform(out,data_source = factor(data_source))
  rownames(out) <- NULL
  return(out)
}

# Read the "biosystems_gene" tab-delimited file downloaded from the
# NCBI FTP site (ftp.ncbi.nih.gov/pub/biosystems). The return value is
# a data frame with one row per gene-pathway association, and the
# following columns: bsid, geneid and score.
read_biosystems_gene <- function (file) {
  out <- read_delim(file,delim = "\t",progress = FALSE,
                    col_names = c("bsid","geneid","score"),
                    col_types = cols("i","i","i"))
  class(out) <- "data.frame"
  return(out)
}

# Read the "biosystems_gene" tab-delimited file downloaded from the
# NCBI FTP site (ftp.ncbi.nih.gov/pub/biosystems), and return an n x m
# sparse adjacency matrix, where n is the number of genes and m is the
# number of BioSystems pathways. The entries of this sparse matrix are
# the "scores" assigned the gene-pathway associations.
read_biosystems_gene_sets <- function (file, bsid2info, gene_info) {

  # Read the "biosystems_gene" file, and retain only gene-pathway
  # associations that have corresponding entries in the gene and
  # pathway tables.
  biosys_gene <- read_biosystems_gene(file)
  biosys_gene <- subset(biosys_gene,
                        is.element(bsid,bsid2info$bsid) &
                        is.element(geneid,gene_info$GeneID))

  # Create the n x m sparse adjacency matrix, where n is the number of
  # genes and m is the number of pathways.
  out <- with(biosys_gene,
              sparseMatrix(i = match(geneid,gene_info$GeneID),
                           j = match(bsid,bsid2info$bsid),x = score,
                           dims = c(nrow(gene_info),nrow(bsid2info))))
  rownames(out) <- gene_info$GeneID
  colnames(out) <- bsid2info$bsid
  return(out)
}

# Read pathway names, gene sets, etc, from the tab-delimited
# ".hgnc.gmt" text file downloaded from the Pathway Commons website
# (https://www.pathwaycommons.org). The return value is a list with
# two list elements, (1) "pathways", the pathway meta-data (columns:
# name, data_source, url), and (2) "gene_sets", an n x m sparse
# adjacency matrix, where n is the number of genes and m is the number
# of Pathway Commons gene sets.
read_pathway_commons_data <- function (file, gene_info) {
  dat <- readLines(file)
  n   <- length(dat)

  # Set up the data frame containing the pathway names and sources.
  pathways <- data.frame(name        = rep("",n),
                         data_source = rep("",n),
                         url         = rep("",n),
                         stringsAsFactors = FALSE)

  # Repeat for each Pathway Commons pathway.
  i <- vector("numeric")
  j <- vector("numeric")
  for (t in 1:n) {

    # Get the pathway URL.
    x <- unlist(strsplit(dat[[t]],"\t",fixed = TRUE))
    pathways[t,"url"] <- x[1]

    # Get the pathway name and source.
    y <- unlist(strsplit(x[2],";",fixed = TRUE))
    y <- trimws(unlist(strsplit(y,":",fixed = TRUE)))
    pathways[t,"name"]        <- y[which(y == "name") + 1]
    pathways[t,"data_source"] <- y[which(y == "datasource") + 1]

    # Get the pathway-gene associations.
    x  <- x[-(1:2)]
    it <- match(x,gene_info$Symbol)
    it <- it[!is.na(it)]
    i  <- c(i,it)
    j  <- c(j,rep(t,length(it)))
  }
  
  # Create the n x m sparse (binary) adjacency matrix, where n is the
  # number of genes and m is the number of Pathway Commons pathways.
  gene_sets <- sparseMatrix(i = i,j = j,x = rep(1,length(i)),
                            dims = c(nrow(gene_info),n))
  rownames(gene_sets) <- gene_info$GeneID
  colnames(gene_sets) <- pathways$url
  
  # Return the pathway meta-data and the n x m adjacency matrix.
  pathways <- transform(pathways,data_source = factor(data_source))
  return(list(pathways = pathways,gene_sets = gene_sets))
}

# Read the MSigDB gene set meta-data from the XML file. The return
# value is a data frame with one row per gene set, and with the
# following columns: standard_name, systematic_name, category_code,
# sub_category_code, organism and description_brief.
read_msigdb_xml <- function (file) {
  out <- read_xml(file)
  out <- as_list(out)[[1]]
  extract_attr <- function (dat, attribute)
    sapply(dat,function (x) attr(x,attribute))
  out <- data.frame(standard_name     = extract_attr(out,"STANDARD_NAME"),
                    systematic_name   = extract_attr(out,"SYSTEMATIC_NAME"),
                    category_code     = extract_attr(out,"CATEGORY_CODE"),
                    sub_category_code = extract_attr(out,"SUB_CATEGORY_CODE"),
                    organism          = extract_attr(out,"ORGANISM"),
                    description_brief = extract_attr(out,"DESCRIPTION_BRIEF"),
                    stringsAsFactors = FALSE)
  return(transform(out,
                   category_code     = factor(category_code),
                   sub_category_code = factor(sub_category_code),
                   organism          = factor(organism)))
}

# Read the MSigDB gene set meta-data from the XML data, and extract
# the MSigDB gene sets using the msigdbr package. The return value is
# a list with two list elements: "gene_sets", the n x m sparse
# adjacency matrix of gene sets, where n is the number of genes and m
# is the number of gene sets; and "info", a data frame with one row
# per gene set, and with columns "standard_name", "systematic_name",
# "category_code", "sub_category_code", "organism" and
# "description_brief".
get_msigdb_gene_sets <- function (file, gene_info, species) {

  # Read the MSigDB gene set meta-data.
  info <- read_msigdb_xml(file)
    
  # Get the gene-set annotations.
  x <- msigdbr(species = species)
  class(x) <- "data.frame"

  # Remove gene-set annotations that do not have a corresponding
  # gene_info entry.
  x   <- subset(x,is.element(entrez_gene,gene_info$GeneID))
  ids <- sort(unique(x$gs_id))

  # Create an n x m sparse (binary) adjacency matrix from the gene-set
  # data.
  gene_sets <- sparseMatrix(i = match(x$entrez_gene,gene_info$GeneID),
                            j = match(x$gs_id,ids),x = rep(1,nrow(x)),
                            dims = c(nrow(gene_info),length(ids)))
  rownames(gene_sets) <- gene_info$GeneID
  colnames(gene_sets) <- ids

  # Align the rows of "info" with the columns of "gene_sets".
  info      <- subset(info,is.element(systematic_name,colnames(gene_sets)))
  cols      <- match(info$systematic_name,colnames(gene_sets))
  gene_sets <- gene_sets[,cols]

  # Return the gene sets (represented as a sparse matrix) and the
  # accompanying meta-data ("info").
  return(list(info = info,gene_sets = gene_sets))
}

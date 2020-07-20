# Read the "gene_info" tab-delimited text file downloaded from the
# NCBI FTP site (ftp.ncbi.nih.gov/gene). The return value is a data
# frame with one row per gene, and the following columns: GeneID,
# Symbol, Synonyms, chromosome, Ensembl and HGNC.
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
# frame with one row per BioSystems pathway---limited to pathways for
# the specified organisms---and the following columns: bsid, tax_id,
# data_source, accession and name.
read_bsid2info <- function (file, organisms = c(9606, 10090)) {
  out <- suppressWarnings(
    read_delim(file,delim = "\t",quote = "",progress = FALSE,
               col_names = c("bsid","data_source","accession","name",
                             "type","scope","tax_id","description"),
               col_types = cols("i","c","c","c","c","c","i","c")))
  class(out) <- "data.frame"
  out <- subset(out,type == "pathway" & is.element(tax_id,organisms))
  out <- out[c("bsid","tax_id","data_source","accession","name")]
  out <- transform(out,
                   tax_id      = factor(tax_id),
                   data_source = factor(data_source))
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
# (https://www.pathwaycommons.org).
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
  
  # Return the pathway meta-data and the n x m gene-set matrix.
  pathways <- transform(pathways,data_source = factor(data_source))
  return(list(pathways = pathways,gene_sets = gene_sets))
}

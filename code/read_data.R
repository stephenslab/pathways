# Read the "gene_info" tab-delimited text file downloaded from the
# NCBI FTP site (ftp.ncbi.nih.gov/gene). The return value is a data
# frame with one row per gene, and the following columns: GeneID,
# Symbol, Synonyms, chromosome, Ensembl and HGNC.
read_gene_info <- function (file) {

  # Read the data into a data frame.
  out <- suppressMessages(read_delim(file,delim = "\t",col_names = TRUE))
  class(out) <- "data.frame"
  dbXrefs    <- out$dbXrefs
  out        <- out[c("GeneID","Symbol","Synonyms","chromosome")]

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

  # Extract the HGNC (HUGO Gene Nomenclature Committee) ids.
  out$HGNC <- sapply(dbXrefs,function (x) {
                i <- which(substr(x,1,10) == "HGNC:HGNC:")
                if (length(i) > 0)
                  return(substr(x[i[1]],6,nchar(x[i[1]])))
                else
                  return(NA)
              })

  # Return the processed gene data.
  return(out)
}

# Read the "bsid2info" tab-delimited file downloaded from the NCBI FTP
# site (ftp.ncbi.nih.gov/pub/biosystems). The return value is a data
# frame with one row per BioSystems pathway, and the following columns:
# bsid, data_source, accession and name.
read_bsid2info <- function (file) {
  out <- suppressWarnings(
    read_delim(file,delim = "\t",quote = "",progress = FALSE,
               col_names = c("bsid","data_source","accession","name",
                             "type","scope","tax_id","description"),
               col_types = cols("i","c","c","c","c","c","i","c")))
  class(out) <- "data.frame"
  out <- subset(out,type == "pathway" & tax_id == 9606)
  out <- transform(out,data_source = factor(data_source))
  out <- out[c("bsid","data_source","accession","name")]
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
# sparse matrix, where n is the number of genes and m is the number of
# BioSystems pathways. The entries of this sparse matrix are the "scores"
# assigned the gene-pathway associations.
read_biosystems_gene_sets <- function (file, bsid2info, gene_info) {

  # Read the "biosystems_gene" file, and retain only gene-pathway
  # associations that have corresponding entries in the gene and
  # pathway tables.
  biosys_gene <- read_biosystems_gene(file)
  biosys_gene <- subset(biosys_gene,
                        is.element(bsid,bsid2info$bsid) &
                        is.element(geneid,gene_info$GeneID))

  # Create an n x m sparse matrix, where n is the number of genes and
  # m is the number of pathways.
  out <- with(biosys_gene,
              sparseMatrix(i = match(geneid,gene_info$GeneID),
                           j = match(bsid,bsid2info$bsid),x = score,
                           dims = c(nrow(gene_info),nrow(bsid2info))))
  rownames(out) <- gene_info$GeneID
  colnames(out) <- bsid2info$bsid
  return(out)
}
  
# Read the HGNC gene data from the tab-delimited text file downloaded
# from the HGNC website (www.genenames.org). Only entries marked as
# being "approved" are outputted.
read_hgnc_data <- function (file) {
  out <- suppressMessages(read_delim(file,delim = "\t"))
  class(out) <- "data.frame"
  names(out) <- c("hgnc_id","approved_symbol","approved_name","status",
                  "previous_symbols","alias_symbols","chr","accession_numbers",
                  "refseq_ids")
  out <-
    transform(out,
      hgnc_id = as.numeric(sapply(strsplit(hgnc_id,":",fixed = TRUE),"[",2)),
      status  = factor(status))
  return(subset(out,status == "Approved"))
}

# Read the pathway meta data (e.g., pathway names, data sources) from
# the tab-delimited text file downloaded from the Pathway Commons
# website (https://www.pathwaycommons.org).
read_pc_pathway_data <- function (file) {
  out <- suppressMessages(read_delim(file,delim = "\t",quote = ""))
  class(out) <- "data.frame"
  names(out) <- tolower(names(out))
  out <- out[c("display_name","datasource","pathway_uri")]
  return(transform(out,datasource = factor(datasource)))
}

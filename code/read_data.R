# Read the HGNC gene data from the tab-delimited text file downloaded
# from the HGNC website (www.genenames.org). Only "approved" entries
# are retained.
read_hgnc_data <- function (file) {
  out <- suppressMessages(read_delim(file,delim = "\t"))
  class(out) <- "data.frame"
  names(out) <- c("hgnc_id","approved_symbol","approved_name","status",
                  "previous_symbols","alias_symbols","chr","accession_numbers",
                  "refseq_ids")
  out <- transform(out,
                   hgnc_id = as.numeric(sapply(strsplit(hgnc_id,":"),"[",2)),
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

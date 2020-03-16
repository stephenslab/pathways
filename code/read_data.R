# Read the pathway meta daata (pathway names and data sources) from
# the tab-delimited text file downloaded from the Pathway Commons
# website.
read_pc_pathways <- function (file) {
  out <- suppressMessages(
           read_delim(file,delim = "\t",col_names = TRUE,quote = ""))
  class(out) <- "data.frame"
  names(out) <- tolower(names(out))
  out <- out[c("display_name","datasource","pathway_uri")]
  return(transform(out,datasource = factor(datasource)))
}

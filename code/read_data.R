# TO DO: Explain here what this function does, and how to use it.
read_pathways <- function (file) {
  out <- suppressMessages(read_delim(file,delim = "\t",col_names = TRUE,
                                     quote = ""))
  class(out) <- "data.frame"
  names(out) <- tolower(names(out))
  out <- out[c("display_name","datasource","pathway_uri")]
  return(transform(out,datasource = factor(datasource)))
}

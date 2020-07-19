# An initial attempt at reading in the Gene Ontology data.
library(readr)
dat <- read_delim("../data/goa_human.gaf.gz",delim = "\t",
                  col_names = FALSE,comment = "!")
class(dat) <- "data.frame"

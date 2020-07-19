# An initial attempt at reading in the x
library(xml2)
library(msigdbr)
out <- msigdbr(species = "Mus musculus")
class(out) <- "data.frame"

dat <- readLines("../data/msigdb.v7.1.entrez.gmt.gz")
n   <- length(dat)
for (t in 1:n) {
  x <- unlist(strsplit(dat[[t]],"\t",fixed = TRUE))
}
                                        # out <- read_xml("../data/msigdb_v7.1.xml")

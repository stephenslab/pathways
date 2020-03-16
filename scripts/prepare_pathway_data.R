# TO DO: Explain here what this script does, and how to use it.
library(readr)
source("../code/read_data.R")

# LOAD DATA
# ---------
# TO DO: Explain here what these few lines of code do.
cat("Reading pathway data from pathways.txt.gz.\n")
out <- read_pathways("../data/pathways.txt.gz")
print(summary(out$datasource))

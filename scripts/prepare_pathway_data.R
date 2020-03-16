# TO DO: Explain here what this script does, and how to use it.
library(readr)
source("../code/read_data.R")

# LOAD DATA
# ---------
# Read the Pathway Commons pathway meta data (e.g., pathway names,
# data sources) from the tab-delimited text file.
cat("Reading Pathway Commons pathway data from pathways.txt.gz.\n")
out <- read_pc_pathways("../data/pathways.txt.gz")

# SUMMARIZE DATA
# --------------
print(summary(out$datasource))

# WRITE DATA TO FILE
# ------------------
# TO DO.

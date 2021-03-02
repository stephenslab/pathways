# Integrate the "data source" or MSigDB "category code" into the the
# "database" column.
gene_set_info <- transform(gene_set_info,
                           database = factor(paste(database,
                                                   ifelse(database == "MSigDB",
                                                          as.character(category_code),
                                                          as.character(data_source)),
                                                   sep = "-")))
gene_set_info <- gene_set_info[c("name","id","database","accession","sub_category_code",
                                 "organism","description_brief")]

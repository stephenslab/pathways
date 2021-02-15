# pathways

Human and mouse gene-set data compiled for gene-set enrichment
analyses. The gene sets are compiled from three sources:
[NCBI BioSystems][biosystems], [Pathway Commons][pc] and
[MSigDB][msigdb].

## Quick Start

Load the gene set data into R, e.g.,

```R
library(Matrix)
load("gene_sets_human.RData")
```

The gene-set data are stored as an n x m sparse matrix, where n is the
number of genes, and m is the number of gene sets. For the human gene
sets, n = 61,676 and m = 37,856.

```R
dim(gene_sets)
```

An entry `gene_sets[i,j]` of sparse matrix is 1 if gene i is annotated
to gene set j; otherwise, the entry is zero. For example, to retrieve
the gene set for the IL12-mediated signaling events pathway from the
BioSystems database, run:

```R
id <- subset(gene_set_info,
             name == "IL12-mediated signaling events" &
             database == "BioSystems")$id
genes <- which(gene_sets[,id] > 0)
```

This will give you the numbers of the rows in the `gene_info`
table. To look up information about these genes, such as their
official symbols, and the Ensembl gene ids, you would do

```R
gene_info[genes,c("GeneID","Symbol","Ensembl")]
```

## Other notes

+ `inst/scripts/prepare_human_gene_set_data.R` and
  `inst/scripts/prepare_mouse_gene_set_data.R` are the R scripts used
  to generate the gene set data files.

+ `inst/code/read_gene_set_data.R` contains some functions
  used in the data preparation scripts.

+ See the `inst/datafiles` directory for data files downloaded from
  their original sources as well as somethe processed data files. Some 
  of the data files are not actually included in the git repository
  because they are large; see below for instructions on downloading the
  data.

+ **Homo_sapiens.gene_info.gz** and **Mus_musculus.gene_info.gz** are
  tab-delimited files containing gene information. These files were
  downloaded from the [NCBI FTP site][ncbi-ftp-gene] on October
  15, 2020.

+ **bsid2info.gz** is a tab-delimited text file containing information
  about the NCBI BioSystems pathways. This file was downloaded from
  the [NCBI FTP site][ncbi-ftp-biosystems] on March 22, 2020.

+ **biosystems_gene.gz** is a tab-delimited file containing gene set
  data for the NCBI BioSystems pathways. This file was downloaded from
  the [NCBI FTP site][ncbi-ftp-biosystems] on March 22, 2020.

+ **PathwayCommons12.All.hgnc.gmt.gz** is a tab-delimited text file
  containing pathway data, including gene sets. This file was
  downloaded from [Pathway Commons][pc-12-downloads] on March
  20, 2020.

+ **msigdb_v7.2.xml.gz** is an XML file containing information about
  the MSigDB gene sets. This file was downloaded from
  [here][msigdb-download] on October 15, 2020. The MSigDB XML format
  is described [here][msigdb-xml-format].

[biosystems]: https://www.ncbi.nlm.nih.gov/biosystems
[pc]: https://www.pathwaycommons.org
[ncbi-ftp-gene]: https://ftp.ncbi.nih.gov/gene
[hgnc]: https://www.genenames.org/download/custom
[ncbi-ftp-biosystems]: https://ftp.ncbi.nih.gov/pub/biosystems
[pc-12-downloads]: https://www.pathwaycommons.org/archives/PC2/v12
[gaf]: http://geneontology.org/docs/go-annotation-file-gaf-format-2.1
[msigdb]: https://www.gsea-msigdb.org/gsea/msigdb
[msigdb-download]: https://www.gsea-msigdb.org/gsea/downloads.jsp
[msigdb-xml-format]: https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/MSigDB_XML_description

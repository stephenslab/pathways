# pathways

Data compiled for gene-set enrichment analyses. The pathway data are
compiled from two sources: [NCBI BioSystems][biosystems] and
[Pathway Commons][pc].

## Quick Start

Load the pathway data into R:

```R
library(Matrix)
load("pathways.RData")
```

The gene-set data are stored as an n x m sparse matrix, where n =
61,630 is the number of genes, and m = 6,737 is the number of
pathways.

```R
dim(gene_sets)
```

For example, to retrieve the gene set for the IL12-mediated signaling
events pathway from Pathway Commons, run:

```R
id <- subset(pathways,
             name == "IL12-mediated signaling events" &
             database == "BioSystems")$id
genes <- which(gene_sets[,id] > 0)
```

This will give you the row numbers of the `gene_info` table. To look
up information about these genes, such as the official gene symbols,
you would do

```R
gene_info[genes,"Symbol"]
```

## Source code

+ **prepare_pathway_data.R** is the R script used to generate the
  **pathways.RData** file.

+ **read_data.R** contains some function definitions used in the data
  preparation script.

## Source data

See the **data** directory for data files downloaded from the original
sources, and the processed data files. Some of the data files are not
actually included in the git repository because they are large; see
below for instructions on downloading the data.

+ **Homo_sapiens.gene_info.gz** and **Mus_musculus.gene_info.gz** are
  tab-delimited files containing gene information. These files were
  downloaded from the [NCBI FTP site][ncbi-ftp-gene] on July 17, 2020.

+ **bsid2info.gz** is a tab-delimited text file containing information
  about the NCBI BioSystems pathways. This file was downloaded from
  the [NCBI FTP site][ncbi-ftp-biosystems] on March 22, 2020.

+ **biosystems_gene.gz** is a tab-delimited file containing gene set
  data for the NCBI BioSystems pathways. This file was downloaded from
  the [NCBI FTP site][ncbi-ftp-biosystems] on March 22, 2020.

+ **PathwayCommons12.All.hgnc.gmt.gz** is a text file containing 
  pathway data, including gene sets. This file was downloaded from
  [Pathway Commons][pc-12-downloads] on March 20, 2020.

+ **go.obo.gz** is a text file containing the full Gene Ontology (GO)
  ontology data in [OBO format][obo]. This file was downloaded from
  [Gene Ontology][go-ontology] on July 17, 2020.

+ **goa_human.gaf.gz** and **mgi.gaf.gz** are text files containing
  Gene Ontology (GO) annotations in [GAF format][gaf]. These files
  were downloaded from [Gene Ontology][go-annotations] on July
  17, 2020.
 
+ **msigdb_v7.1.xml.gz** and **msigdb.v7.1.entrez.gmt.gz** are text
  files containing information about the MSigDB gene sets. These files
  were downloaded from the [GSEA website][gsea] on July 18, 2020.

[biosystems]: https://www.ncbi.nlm.nih.gov/biosystems
[pc]: https://www.pathwaycommons.org
[ncbi-ftp-gene]: https://ftp.ncbi.nih.gov/gene
[hgnc]: https://www.genenames.org/download/custom
[ncbi-ftp-biosystems]: https://ftp.ncbi.nih.gov/pub/biosystems
[pc-12-downloads]: https://www.pathwaycommons.org/archives/PC2/v12
[go-ontology]: http://geneontology.org/docs/download-ontology
[go-annotations]: http://current.geneontology.org/products/pages/downloads.html
[obo]: http://owlcollab.github.io/oboformat/doc/obo-syntax.html
[gaf]: http://geneontology.org/docs/go-annotation-file-gaf-format-2.1
[gsea]: https://www.gsea-msigdb.org/gsea/downloads.jsp

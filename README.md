# pathways

Human and mouse gene sets compiled for gene set enrichment analysis,
and a simple interface for performing gene set enrichment analysis
using fgsea. The gene sets are compiled from NCBI BioSystems, Pathway
Commons and MSigDB. The gene sets are compiled from three sources:
[NCBI BioSystems][biosystems], [Pathway Commons][pc] and 
[MSigDB][msigdb].

See [here](https://stephenslab.github.io/pathways/gsea_b_cells.html)
for an example of an interactive plot for exploring the results of a
gene set enrichment analysis on differential expression in B cells
vs. other immune cell populations.

## Quick Start

Install the package:

```R
library(remotes)
install_github("stephenslab/pathways")
```

Load the package:

```R
library(pathways}
```

Try running the gene set enrichment analysis example:

```R
example("perform_gsea")
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

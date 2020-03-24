# pathways

Pathway data compiled for gene-set enrichment analyses.

See the **data** directory for data files downloaded from the original
sources, and the processed data files. Some of the data files are not
actually included in the git repository because they are large; see
below for instructions on downloading the data.

## Data

+ **Homo_sapiens.gene_info.gz** is a tab-delimitd file containing gene
  information. This file was downloaded from the
  [NCBI FTP site][ncbi-ftp-gene] on March 22, 2020.

+ **bsid2info.gz** is a tab-delimited text file containing information
  about the NCBI BioSystems pathways. This file was downloaded from
  the [NCBI FTP site][ncbi-ftp-biosystems] on March 22, 2020.

+ **biosystems_gene.gz** is a tab-delimited file containing gene set
  data for the NCBI BioSystems pathways. This file was downloaded from
  the [NCBI FTP site][ncbi-ftp-biosystems] on March 22, 2020.

+ **PathwayCommons12.All.hgnc.gmt.gz** is a text file containing 
  pathway data, including gene sets. This file was downloaded from
  [Pathway Commons][pc-12-downloads] on March 20, 2020.

## Data 

[ncbi-ftp-gene]: https://ftp.ncbi.nih.gov/gene
[hgnc]: https://www.genenames.org/download/custom
[ncbi-ftp-biosystems]: https://ftp.ncbi.nih.gov/pub/biosystems
[pc-12-downloads]: https://www.pathwaycommons.org/archives/PC2/v12

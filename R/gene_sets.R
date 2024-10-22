#' @name gene_sets_human
#'
#' @title Human Gene Sets
#'
#' @docType data
#'
#' @description 37,856 gene sets compiled from the NCBI BioSystems,
#'   Pathway Commons and MSigDB databases.
#'
#' @format \code{gene_sets_human} is a list with the following
#'   list elements:
#' 
#' \describe{
#'
#'   \item{gene_info}{Data frame containing information on human
#'     genes, including gene symbols, HGNC and Ensembl ids. This table was
#'     compiled from the \code{Homo_sapiens.gene_info.gz} file downloaded
#'     from the NCBI FTP site.}
#'
#'   \item{gene_set_info}{Data frame containing information on the
#'     gene sets, including gene set name, id and database of origin. The
#'     columns "sub_category_code", "organism" and "description_brief" are
#'     only used for the MSigDB gene sets.}
#'
#'   \item{gene_sets}{Gene sets encoded as a 61,676 x 37,856 sparse
#'     binary matrix, in which \code{gene_sets[i,j] = 1} if and only if
#'     gene i is included in gene set j; otherwise, \code{gene_sets[i,j] =
#'     0}. The row names are the gene Ensembl ids and the column names
#'     are the gene set ids.}}
#'
#' @keywords data
#'
NULL

#' @name gene_sets_mouse
#'
#' @title Mouse Gene Sets
#'
#' @docType data
#'
#' @description 33,380 gene sets compiled from NCBI BioSystems,
#'   Pathway Commons and MSigDB databases.
#'
#' @format \code{gene_sets_mouse} is a list with the following
#'   list elements:
#' 
#' \describe{
#'
#'   \item{gene_info}{Data frame containing information on mouse
#'     genes, including gene symbols and Ensembl ids. This table was
#'     compiled from the \code{Mus_musculus.gene_info.gz} file downloaded
#'     from the NCBI FTP site.}
#'
#'   \item{gene_set_info}{Data frame containing information on the
#'     gene sets, including gene set name, id and database of origin. The
#'     columns "sub_category_code", "organism" and "description_brief" are
#'     only used for the MSigDB gene sets.}
#'
#'   \item{gene_sets}{Gene sets encoded as a 73,202 x 33,380 sparse
#'     binary matrix, in which \code{gene_sets[i,j] = 1} if and only if
#'     gene i is included in gene set j; otherwise, \code{gene_sets[i,j] =
#'     0}. The row names are the gene Ensembl ids and the column names are
#'     the gene set ids.}}
#'
#' @keywords data
#'
NULL

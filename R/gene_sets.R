#' @name gene_sets_human
#'
#' @title Human Gene Sets
#'
#' @docType data
#'
#' @description 37,856 gene sets compiled from NCBI BioSystems,
#'   Pathway Commons and MSigDB databases.
#'
#' @format \code{gene_sets_human} is a list with the following
#'   list elements:
#' 
#' \describe{
#'
#'   \item{gene_info}{Data frame containing information on human
#'     genes, including gene symbols and Ensembl ids.}
#'
#'   \item{gene_set_info}{Data frame containing information on gene
#'     sets, including gene set name, id and database of origin.}
#'
#'   \item{gene_sets}{Gene sets encoded as a 61,676 x 37,856 sparse
#'     matrix, in which gene_sets[i,j] = 1 if gene j belongs to gene set
#'     i; otherwise, gene_sets[i,j] = 0.}}
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
#'     genes, including gene symbols and Ensembl ids.}
#'
#'   \item{gene_set_info}{Data frame containing information on gene
#'     sets, including gene set name, id and database of origin.}
#'
#'   \item{gene_sets}{Gene sets encoded as a 73,202 33,380 sparse
#'     matrix, in which gene_sets[i,j] = 1 if gene j belongs to gene set
#'     i; otherwise, gene_sets[i,j] = 0.}}
#'
#' @keywords data
#'
NULL

#' @title Add Title Here
#'
#' @description Add description here
#'
#' @return Describe return value here.
#'
#' @examples
#' # Add an example here.
#' 
#' @export
#' 
perform_gsea <- function () {

}

#' @importFrom fgsea fgsea
perform_gsea_helper <- function (gene_sets, z, eps = 1e-32, nproc = 1, ...) {

  # Convert the gene sets adjacency matrix into the format accepted by
  # fgsea.
  pathways <- matrix2pathways(gene_sets)

  # Perform gene set enrichment analysis using fgsea.
  out <- suppressWarnings(fgsea(pathways,z,eps = eps,nproc = nproc,...))
  class(out) <- "data.frame"

  # TO DO: Explain here what these lines of code go.
  rownames(out) <- out$pathway
  out <- out[c("pval","log2err","ES","NES")]
  out <- out[colnames(gene_sets),]
  out[is.na(out$log2err),] <- NA
  return(out)
}

# Recover a list of gene sets from an n x m adjacency matrix, A, in
# which n is the number of genes and m is the number of gene sets;
# A[i,j] = 1 if and only if gene i is included in gene set j. This
# function is used to prepare the "pathways" fgsea input from a
# collection of gene sets encoded as a (sparse) matrix. For this
# function to work, the rows and columns of A must be named.
matrix2pathways <- function (A) {
  n          <- ncol(A)
  out        <- vector("list",n)
  names(out) <- colnames(A)
  genes      <- rownames(A)
  for (i in 1:n)
    out[[i]] <- genes[A[,i] == 1]
  return(out)
}

#' @title Add Title Here
#'
#' @description Add description here.
#'
#' @param gene_sets Describe input argument "gene_sets" here.
#'
#' @param Z Describe input argument "Z" here.
#'
#' @param verbose Describe input argument "verbose" here.
#' 
#' @param eps Describe input argument "eps" here.
#'
#' @param nproc Describe input argument "nproc" here.
#'
#' @param \dots Describe dots (...) here.
#' 
#' @return Describe return value here.
#'
#' @examples
#' # Add an example here.
#' 
#' @export 
#' 
perform_gsea <- function (gene_sets, Z, verbose = TRUE, eps = 1e-32,
                          nproc = 1, ...) {

  # Get the number of gene sets (n) and the number of columns in Z (k).
  n <- ncol(gene_sets)
  k <- ncol(Z)
  if (verbose) {
    cat("Computing enrichment statistics for",n,"gene sets and")
    cat(k,"collections of gene-wise statistics:")
  }
  if (is.null(colnames(Z)))
    colnames(Z) <- paste0("k",1:k)

  browser()
  
  # Align the gene and gene set data.
  out       <- align_by_rownames(gene_sets,Z)
  gene_sets <- out$A
  Z         <- out$Z
  
  # Initialize the outputs.
  out <- list(pval    = matrix(0,n,k),
              log2err = matrix(0,n,k),
              ES      = matrix(0,n,k),
              NES     = matrix(0,n,k))
  rownames(out$pval)    <- colnames(gene_sets)
  rownames(out$log2err) <- colnames(gene_sets)
  rownames(out$ES)      <- colnames(gene_sets)
  rownames(out$NES)     <- colnames(gene_sets)
  colnames(out$pval)    <- colnames(Z)
  colnames(out$log2err) <- colnames(Z)
  colnames(out$ES)      <- colnames(Z)
  colnames(out$NES)     <- colnames(Z)
  
  # Run the gene set enrichment analysis for each column of Z.
  for (i in 1:k) {
    if (verbose) {
      if (i > 1)
        cat(", ")
      cat(rownames(Z)[i])
    }
    ans             <- perform_gsea_helper(gene_sets,Z[,i],eps,nproc,...)
    out$pval[,i]    <- ans$pval
    out$log2err[,i] <- ans$log2err
    out$ES[,i]      <- ans$ES
    out$NES[,i]     <- ans$NES
  }
  if (verbose)
    cat("\n")

  # Output the results of the gene set enrichment analyses.
  return(out)
}

# TO DO: Explain what this function does, and how to use it.
#  
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

# This function aligns the gene set data (A) with the gene statistics
# (Z) by their row names to prepare these data for a gene set
# enrichment analysis.
align_by_rownames <- function (A, Z) {
  x   <- rownames(A)
  y   <- rownames(Z)
  ids <- intersect(x,y)
  i   <- match(ids,x)
  j   <- match(ids,y)
  A   <- A[i,]
  Z   <- Z[j,]
  return(list(A = A,Z = Z))
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

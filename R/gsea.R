#' @title Perform Gene Set Enrichment Analysis  
#'
#' @description A simple interface for performing gene set enrichment
#'   analysis on multiple collections of gene-level statistics using
#'   fgsea.
#'
#' @param gene_sets Gene set data encoded as an n x m binary matrix,
#'   where n is the number of gene sets and m is the number of genes:
#'   \code{gene_sets[i,j] = 1} if and only if gene i is included in gene
#'   set j, and otherwise \code{gene_sets[i,j] = 0}. The rows and
#'   columns should be named. Only the rows of \code{gene_sets} matching
#'   the rows of \code{Z} will be used in the enrichment analysis.
#'
#' @param Z Matrix of gene-level statistics such as z-scores, with
#'   rows corresponding to genes. An enrichment analysis is performed
#'   for each column of \code{Z}. The rows and columns should be named;
#'   rows of \code{Z} matching the columns of \code{gene_sets} will
#'   be used in the enrichment analysis.
#'
#' @param verbose When \code{verbose = TRUE}, information about the
#'   progress of the gene set enrichment analysis is printed to the
#'   console.
#' 
#' @param eps The lower bound for calculating p-values; smaller values
#'   of \code{eps} may give more accurate p-values at the possible cost
#'   of slightly longer computation. Passed as the \code{eps} input to
#'   \code{\link[fgsea]{fgsea}}.
#'
#' @param nproc Number of workers used to run gene set enrichment
#'   analysis. Passed as the \code{nproc} input to
#'   \code{\link[fgsea]{fgsea}}.
#'
#' @param \dots Additional arguments passed to
#'   \code{\link[fgsea]{fgsea}}.
#' 
#' @return The return value is a list containing four n x k matrices
#' of gene set enrichment analysis results, where n is the number of
#' gene sets and k is the number of columns in \code{Z}. The matrices
#' give the p-values (pval), enrichment scores (ES), normalized
#' enrichment scores (NES), and expected errors (log2err). See
#' \code{\link[fgsea]{fgseaMultilevel}} for more details about these
#' outputs.
#'
#' @examples
#' data(gene_sets_human)
#' data(pbmc_facs_z)
#' gsea_res <- perform_gsea(gene_sets_human$gene_sets,pbmc_facs_z,
#'                          nproc = 4)
#'
#' @importFrom fgsea fgsea
#' 
#' @export 
#' 
perform_gsea <- function (gene_sets, Z, verbose = TRUE, eps = 1e-32,
                          nproc = 0, ...) {

  # Check the "gene_sets" and "Z" inputs.
  if (!(is.matrix(gene_sets) | inherits(gene_sets,"Matrix")) & is.matrix(Z) &
        !is.null(rownames(gene_sets)) & !is.null(colnames(gene_sets)) &
        !is.null(rownames(Z)) & !is.null(colnames(Z)))
    stop("Input arguments \"gene_sets\" and \"Z\" should be matrices with ",
         "named rows and columns")
    
  # Get the number of gene sets (n) and the number of columns in Z (k).
  n <- ncol(gene_sets)
  k <- ncol(Z)

  # Align the gene and gene set data.
  out       <- align_by_rownames(gene_sets,Z)
  gene_sets <- out$A
  Z         <- out$Z
  if (nrow(Z) == 0)
    stop("Failed to align gene_sets and Z by row names.")
  if (verbose)
    cat("Using statistics from",nrow(Z),"genes for gene set enrichment",
        "analysis.\n")
  
  # Convert the gene sets binary adjacency matrix into the format
  # accepted by fgsea.
  if (verbose)
    cat("Converting gene set data for fgsea.\n")
  pathways <- matrix2pathways(gene_sets)
  
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
  if (verbose) {
    cat("Computing enrichment statistics for",n,"gene sets and\n")
    cat(k,"collections of gene-level statistics:\n")
  }
  for (i in 1:k) {
    if (verbose) {
      if (i > 1)
        cat(", ")
      cat(colnames(Z)[i])
    }

    # Perform gene set enrichment analysis using fgsea.
    ans <- suppressWarnings(fgsea(pathways,Z[,i],eps = eps,nproc = nproc,...))
    class(ans) <- "data.frame"

    # Refactor the outputs as a data frame with one row per gene set,
    # and the following columns: pval, log2err, ES and NES. The
    # outputs are set to mising (NA) whenever the expected error
    # (log2err) is also NA.
    rownames(ans) <- ans$pathway
    ans <- ans[c("pval","log2err","ES","NES")]
    ans <- ans[colnames(gene_sets),]
    ans[is.na(ans$log2err),] <- NA

    # Store the fgsea results.
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

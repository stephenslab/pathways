gsea_plot_shapes <- rep(c("circle","square"),10)
gsea_plot_colors <- c("gold",         # biocyc
                      "yellowgreen",  # C1
                      "dodgerblue",   # C2
                      "olivedrab",    # C3
                      "firebrick",    # C4
                      "darkorange",   # C5
                      "magenta",      # C6
                      "darkblue",     # C7
                      "gold",         # H
                      "tomato",       # humancyc
                      "olivedrab",    # inoh
                      "darkblue",     # kegg
                      "royalblue",    # netpath
                      "darkorange",   # panther
                      "yellowgreen",  # pathbank
                      "magenta",      # pid
                      "dodgerblue")   # reactome

#' @title Add Title Here
#' 
#' @description Create an interactive scatterplot using plotly to
#' explore the gene-set enrichment results for a given topic k. Input
#' argument "gene_set_info" is a data frame containing information
#' about the gene sets; "gsea_res" is an output from function
#' "perform_gsea"; and "label_gene_sets" is a vector of gene set ids
#' to be labeled in the plot.
#'
#' @param gsea_res Describe input argument "gsea_res" here.
#'
#' @param gene_set_info Describe input argument "gene_set_info" here.
#'
#' @param k Describe input argument "k" here.
#'
#' @param file Describe input argument "file" here.
#'
#' @param height Describe input argument "height" here.
#'
#' @param width Describe input argument "width" here.
#'
#' @param title Describe input argument "title" here.
#'
#' @param max_name_len Describe input argument "max_name_len" here.
#' 
#' @return Describe the return value here.
#'
#' @seealso \code{\link{perform_gsea}}
#' 
#' @examples
#' # See perform_gsea for examples.
#' 
#' @importFrom htmlwidgets saveWidget
#' @importFrom plotly plot_ly
#' @importFrom plotly add_trace
#' @importFrom plotly layout
#'
#' @export
#'
gsea_plotly <- function (gsea_res, gene_set_info, k, file, 
                         height = 600, width = 800, title = NULL,
                         max_name_len = 44) {

  # Compile the data for plotting.
  dat <- compile_data_for_gsea_plot(gene_set_info,gsea_res,k,min_pval,
                                    max_name_len)
  
  # Create the plotly plot.
  p <- plot_ly(data = dat,x = ~p0,y = ~p1,color = ~database,symbol = ~database,
               text = ~sprintf("%s\nid: %s\np-value: %0.2e",name,id,10^(-p1)),
               type = "scatter",mode = "markers",hoverinfo = "text",
               symbols = gsea_plot_shapes,colors = gsea_plot_colors,
               marker = list(size = 8,line = list(color = "white",width = 1)),
               height = height,width = width)
  p <- add_trace(p,data = data.frame(x = range(c(dat$p0,dat$p1)),
                                     y = range(c(dat$p0,dat$p1))),
                 x = ~x,y = ~y,mode = "lines",type = "scatter",
                 inherit = FALSE,showlegend = FALSE,
                 line = list(color = "lightgray",dash = "dash",size = 0.3))
  p <- layout(p,
              legend = list(font = list(size = 10)),
              xaxis = list(title = paste("most extreme signed -log10",
                                         "p-value in other columns"),
                           showgrid = FALSE,showline = FALSE,zeroline = FALSE),
              yaxis = list(title = "signed -log10 p-value in selected column",
                           showgrid = FALSE,showline = FALSE,zeroline = FALSE),
              hoverlabel = list(bgcolor = "white",bordercolor = "black",
                                font = list(color = "black",family = "arial",
                                            size = 12)),
              font = list(family = "arial",size = 12),
              title = title)
  if (!missing(file))
    saveWidget(p,file,selfcontained = TRUE,title = title)
  return(p)
}

# Compile the data frame used for create_gsea_plotly.
compile_data_for_gsea_plot <- function (gene_set_info, gsea_res, k,
                                        min_pval, max_name_len) {
  ES   <- gsea_res$ES
  pval <- gsea_res$pval
  
  # If necessary, get the column number.
  if (is.character(k))
    k <- which(colnames(ES) == k)

  # Take care of NAs.
  ES[is.na(ES)]     <- 0
  pval[is.na(pval)] <- 1
  
  # Compute the signed -log10 p-values.
  P <- -sign(ES) * log10(pval)

  # Compute the "most extreme" signed p-value among other columns.
  n  <- nrow(P)
  p0 <- rep(0,n)
  for (i in 1:n) {
    if (sign(ES[i,k]) > 0)
      p0[i] <- max(P[i,-k])
    else
      p0[i] <- min(P[i,-k])
  }

  # Compile the data for plotting.
  return(data.frame(p1         = P[,k],
                    p0         = p0,
                    name       = substr(gene_set_info$name,1,max_name_len),
                    id         = gene_set_info$id,
                    database   = gene_set_info$database,
                    stringsAsFactors = FALSE))
}

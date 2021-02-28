#' @title Interactive Gene Set Enrichment Analysis Plot
#' 
#' @description Create an interactive scatterplot to explore the
#'   results of a gene set enrichment analysis. In the plot, one point
#'   is drawn for each gene set. The color and shape of the points are
#'   varied the gene set database. To compare results across multiple
#'   analyses, the p-value for the selected analysis (or column) is
#'   plotted on the vertical axis, and the \dQuote{most extreme} p-value
#'   from the other analyses (or columns) is plotted along the
#'   horizontal axis.
#'
#' @param gsea_res An output from \code{\link{perform_gsea}}.
#'
#' @param gene_set_info A data frame containing information about the
#'   gene sets. It should have, at a minimum, columns \dQuote{name},
#'   \dQuote{id} and \dQuote{database}.
#'
#' @param k The column of the gene set enrichment results used for the
#'   plot.
#'
#' @param file Save the interactive plot to this HTML file using
#'   \code{\link[htmlwidgets]{saveWidget}}.
#'
#' @param height Height of the plot in pixels. Passed as argument
#'   \dQuote{height} to \code{\link[plotly]{plot_ly}}.
#'
#' @param width Width of the plot in pixels. Passed as argument
#'   \dQuote{width} to \code{\link[plotly]{plot_ly}}.
#'
#' @param title The text used for the plot title.
#'
#' @param max_name_len Pathway names longer than this are shortened to
#'   this number of characters.
#'
#' @param colors Colors used to draw the points in the plot; passed as
#'   argument \dQuote{colors} to \code{plot_ly}.
#'
#' @param shapes Shapes used to draw the points in the plot; passed as
#'   argument \dQuote{shapes} to \code{plot_ly}.
#' 
#' @return A \code{plotly} object.
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
                         height = 600, width = 800,
                         title = "gsea plot", max_name_len = 44,
                         colors = c("gold","darkblue","tomato","dodgerblue",
                                    "darkmagenta","yellowgreen","dodgerblue",
                                    "olivedrab","firebrick","darkorange",
                                    "magenta","darkblue","gold","tomato",
                                    "gold","olivedrab","darkblue",
                                    "darkorange","magenta","yellowgreen",
                                    "tomato","dodgerblue"),
                         shapes = c(rep("circle",4),rep("x",10),
                                    rep("diamond",8))) {
    
  # Compile the data for plotting.
  dat <- compile_data_for_gsea_plot(gene_set_info,gsea_res,k,min_pval,
                                    max_name_len)
  
  # Create the plotly plot.
  p <- plot_ly(data = dat,x = ~p0,y = ~p1,color = ~database,symbol = ~database,
               text = ~sprintf("%s\nid: %s\np-value: %0.2e",name,id,10^(-p1)),
               type = "scatter",mode = "markers",hoverinfo = "text",
               symbols = shapes,colors = colors,
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

gsea_plot_shapes <- c(23,24,21)
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
#' @param margin Describe input argument "margin" here.
#'
#' @param height Describe input argument "height" here.
#'
#' @param width Describe input argument "width" here.
#'
#' @param title Describe input argument "title" here.
#' 
#' @return Describe return value here.
#'
#' @seealso \code{\link{perform_gsea}}
#' 
#' @examples
#' # See perform_gsea for examples.
#' 
#' @importFrom plotly plot_ly
#' @importFrom plotly add_trace
#' @importFrom plotly layout
#'
#' @export
#'
gsea_plotly <- function (gsea_res, gene_set_info, k, margin = 1,
                         height = 675, width = 800, title = NULL) {

  # Compile the data for plotting.
  dat <- compile_data_for_gsea_plot(gene_set_info,gsea_res,k)
  
  # Create the plotly plot.
  pmin <- min(c(dat$p0,dat$p1))
  pmax <- max(c(dat$p0,dat$p1))
  p <- plot_ly(data = dat,x = ~p0,y = ~p1,symbol = ~database,
               color = ~collection,
               text = ~sprintf("%s\nid: %s\np-value: %0.2e",name,id,10^(-p1)),
               type = "scatter",mode = "markers",hoverinfo = "text",
               symbols = gsea_plot_shapes,colors = gsea_plot_colors,
               marker = list(size = 8,line = list(color = "white",width = 1)),
               height = height,width = width)
  p <- add_trace(p,data = data.frame(x = c(pmin,pmax),y = c(pmin,pmax)),
                 x = ~x,y = ~y,mode = "lines",type = "scatter",
                 inherit = FALSE,showlegend = FALSE,
                 line = list(color = "lightgray",dash = "dash",size = 0.3))
  p <- layout(p,
              legend = list(font = list(size = 10)),
              xaxis = list(title="most extreme signed p-value in other topics",
                           showgrid = FALSE,showline = FALSE,zeroline = FALSE),
              yaxis = list(title = "signed p-value in topic",
                           showgrid = FALSE,showline = FALSE,zeroline = FALSE),
              hoverlabel = list(bgcolor = "white",bordercolor = "black",
                                font = list(color = "black",family = "arial",
                                            size = 12)),
              font = list(family = "arial",size = 12),
              title = title)
  return(p)
}

# Compile the data frame used for create_gsea_plotly.
compile_data_for_gsea_plot <- function (gene_set_info, gsea_res, k,
                                        max_name_len = 44) {

  # Compute the "signed p-values".
  P <- with(gsea_res,-sign(ES) * log10(pval))

  # Compute the "most extreme" signed p-value among other topics.
  n  <- nrow(P)
  p0 <- rep(0,n)
  rows <- which(!is.na(gsea_res$ES[,k]))
  for (i in rows) {
    if (sign(gsea_res$ES[i,k]) > 0)
      p0[i] <- max(P[i,-k],na.rm = TRUE)
    else
      p0[i] <- min(P[i,-k],na.rm = TRUE)
  }

  # Compile the data for plotting.
  dat <- data.frame(p1         = P[,k],
                    p0         = p0,
                    name       = substr(gene_set_info$name,1,max_name_len),
                    id         = gene_set_info$id,
                    database   = gene_set_info$database,
                    collection = with(gene_set_info,
                                      ifelse(database == "MSigDB",
                                             as.character(category_code),
                                             as.character(data_source))),
                    stringsAsFactors = FALSE)
  
  # Subsample the gene sets with p-values close to zero because those
  # gene sets are not particularly interesting, and there are many of
  # them, which slows down the plotting.
  i   <- which(with(dat,!(abs(p0) < 3 & abs(p1) < 3)))
  j   <- which(with(dat,abs(p0) < 3 & abs(p1) < 3))
  n   <- length(j)
  j   <- sample(j,ceiling(n/10))
  dat <- dat[c(i,j),]

  # Set any missing signed p-values to zero.
  dat[is.na(dat$p1),"p1"] <- 0
  dat[is.na(dat$p0),"p0"] <- 0
  return(dat)
}

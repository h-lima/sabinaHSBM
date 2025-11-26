#' @name plot_interaction_matrix
#'
#' @title Plot binary Matrix
#'
#' @description This function plots a (recontructed) binary matrix with options to reorder the matrix and customize the colors and axis titles.
#'
#' @param adj_mat A \code{matrix} representing the interaction data to be plotted.
#' @param order_mat (\emph{optional, default} \code{TRUE}) \cr
#' A \code{logical} value indicating whether the matrix should be order by row frequencies.
#' @param y_title (\emph{optional, default} \code{"Rows"}) \cr
#' A \code{character} string specifying the title for the y-axis.
#' @param x_title (\emph{optional, default} \code{"Columns"}) \cr
#' A \code{character} string specifying the title for the x-axis.
#' @param ... Further arguments for \code{image} R native function
#'
#' @details
#' The \code{order_mat} argument indicates whether the interaction matrix should
#' be iteratively reordered to enhance the grouping of links for improved
#' visualization.
#'
#' The default colors for the matrix are "white" for 0s and "red4" for 1s.
#' If a single col is provided it will change the value for the 1s color. If two or more colors
#' are provided in a vector, it will fallback to the behavior of the \code{image} function
#' \code{col} argument.
#'
#' @return
#' A plot of the interaction matrix.
#'
#' @examples
#' # Load example HSBM prediction results
#' data(dat, package = "sabinaHSBM")
#' 
#' # Plot the Known matrix
#' plot_interaction_matrix(dat)
#' 
#' # Plot the Reconstructed matrix
#' # data(myReconst, package = "sabinaHSBM")
#' # plot_interaction_matrix(myReconst$new_mat, order_mat = FALSE)
#'
#' @seealso [graphics::image()]
#'
#' @importFrom graphics box image
#'
#' @export
plot_interaction_matrix <- function(adj_mat, order_mat = TRUE,
                                    y_title = "Rows", x_title = "Columns", ...){

    args <- list(...)
    if(!("col" %in% names(args))){
        col = c("white", "red4")
    }else{
        col <- args$col
    }
    if(length(col) == 1){
        col <- c("white", col)
    }
    if(order_mat){
        order_ij <- order_com(adj_mat, get_orders = TRUE)
        adj_mat <- adj_mat[order(order_ij$i),
                           order(order_ij$j)]
    }
    adj_mat <- t(adj_mat[nrow(adj_mat):1, ])

    args$x <- 1:nrow(adj_mat)
    args$y <- 1:ncol(adj_mat)
    args$z <- adj_mat
    args$ylab <- y_title
    args$xlab <- x_title
    args$col <- col

    do.call(image, args)
    box()

}

# This is a almost straight copy of a section
# of nestedtemp non-exported function in vegan package (GPL2 Licensed)
order_com <- function(comm, get_orders = FALSE){
    colpack <- function(x, rr)
    {
        ind <- matrix(rep(rr, ncol(x)), nrow=nrow(x))
        s <- -colSums((x*ind)^2)
        t <- -colSums((nrow(x) - (1-x)*ind + 1)^2)
        st <- rank(s+t, ties.method = "random")
        st
    }
    rowpack <- function(x, cr)
    {
        ind <- matrix(rep(cr, each=nrow(x)), nrow=nrow(x))
        s <- -rowSums((x*ind)^2)
        t <- -rowSums((ncol(x) - (1-x)*ind + 1)^2)
        st <- rank(s+t, ties.method = "random")
        st
    }
    comm <- ifelse(comm > 0, 1, 0)
    ## Start with columns, expect if nrow > ncol
    if (ncol(comm) >= nrow(comm)) {
        i <- rank(-rowSums(comm), ties.method = "average")
    } else {
        j <- rank(-colSums(comm), ties.method = "average")
        i <- rowpack(comm, j)
    }
    ## Improve eight times
    for (k in seq_len(8)) {
        j <- colpack(comm, i)
        i <- rowpack(comm, j)
    }
    if (ncol(comm) < nrow(comm))
        j <- colpack(comm, i)
    comm <- comm[order(i), order(j)]

    if(get_orders){
        return(list(i = i, j = j))
    }
    return(comm)
}


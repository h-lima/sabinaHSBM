#' @name plot_interaction_matrix
#'
#' @title Plot binary bipartite Matrix
#'
#' @description This function plots a (recontructed) binary bipartite matrix with options to reorder the matrix and customize the colors and axis titles.
#'
#' @param adj_mat A \code{matrix} or \code{data.frame} representing the interaction data to be plotted.
#' @param col (\emph{optional, default} \code{"red4"}) \cr
#' A \code{character} string specifying the color to be used for the linked node pairs.
#' @param order_mat (\emph{optional, default} \code{TRUE}) \cr
#' A \code{logical} value indicating whether the matrix should be reordered.
#' @param y_title (\emph{optional, default} \code{"Rows"}) \cr
#' A \code{character} string specifying the title for the y-axis.
#' @param x_title (\emph{optional, default} \code{"Columns"}) \cr
#' A \code{character} string specifying the title for the x-axis.
#'
#' @details
#' The \code{order_mat} argument indicates whether the interaction matrix should be iteratively reordered to enhance the grouping of links for improved visualization.
#'
#' @return
#' A plot of the interaction matrix.
#'
#' @examples
#' # Load example HSBM prediction results
#' data(dat, package = "sabinaHSBM")  #@@@JMB pensar ejemplo para data directory
#' 
#' # Plot the Known matrix
#' plot_interaction_matrix(dat)
#' 
#' # Plot the Reconstructed matrix
#' # data(myReconst, package = "sabinaHSBM")  #@@@JMB pensar ejemplo para data directory
#' # plot_interaction_matrix(myReconst$new_mat, order_mat = FALSE)
#'
#' @export
plot_interaction_matrix <- function(adj_mat, col = "red4", order_mat = TRUE,
                                    y_title = "Rows", x_title = "Columns"){

    if(order_mat){
        adj_mat <- order_com(adj_mat)
    }
    adj_mat <- t(adj_mat[nrow(adj_mat):1, ])

    image(1:nrow(adj_mat), 1:ncol(adj_mat), adj_mat,
          ylab = y_title, xlab = x_title,
          col = c("white", col))
    box()

}

#' @export
plot.hsbm.input <- function(hsbm_input, col = "red4", order_mat = TRUE){

    mat <- hsbm_input$data
    plot_interaction_matrix(mat, col = col, order_mat = order_mat)

}

# This is a almost straight copy of a section
# of nestedtemp function in vegan package (GPL2 Licensed)
order_com <- function(comm){
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

    return(comm)
}

# This is a almost straight copy of a section
# of nestedtemp function in vegan package (GPL2 Licensed)
order_com <- function(comm){
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

    return(comm)
}

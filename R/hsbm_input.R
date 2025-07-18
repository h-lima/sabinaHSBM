#' @name hsbm.input
#'
#' @title Prepare input data for Hierarchical Stochastic Block Model (HSBM) analyses.
#'
#' @description This function prepares the input data necessary for running the \bold{Hierarchical Stochastic Block Model (HSBM)} analysis by creating cross-validation folds and corresponding edgelists.
#'
#' @param data A binary bipartite \code{matrix} containing the input data where rows and columns represent nodes of two distinct sets of elements, and where the entries (i, j) represent the interactions, edges or links: 1 signifies a relationship/link, and 0 signifies no relationship/no link.
#' @param n_folds (\emph{optional, default} \code{5}) \cr
#' A \code{numeric} value specifying the number of cross-validation folds to generate. If only 1 fold is specified, held-out links will not be computed. This is not recommended, as the main objective of the package is to predict missing links through cross-validation.
#' @param folds (\emph{optional, default} \code{NULL}) \cr
#' An optional \code{matrix} of precomputed cross-validation folds for the input data \code{data}. If \code{NULL} (default), folds are generated internally using \code{create_cv_folds()}.
#' @param min_per_col A \code{numeric} value specifying the minimum number of non-zero entries required per column to be considered for removal. Default is 2.
#' @param min_per_row A \code{numeric} value specifying the minimum number of non-zero entries required per row to be considered for removal. Default is 2.
#' #@@@JMB pendiente param add_spurious?? y no_heldout
#'
#' @return
#' An object of class \code{hsbm.input} containing the organized input data, the cross-validation fold assignment for each held-out edge/link, and a list of the corresponding edge lists generated for each fold for HSBM analysis:
#' - \code{$data} The binary bipartite \code{matrix} of input data.
#' - \code{$folds} A \code{matrix} of cross-validation fold assignments for each held-out edge/link.
#' - \code{$edgelist} A \code{list} of edge lists generated for each fold. Each edge list is a \code{data.frame} with the following columns:
#'   - \code{v1} The index of the first (in rows) type of nodes.
#'   - \code{v2} The index of the second (in columns) type of nodes.
#'   - \code{value} The measurement value for the edge. It is 1 if there is a edge/link between the node pair, and 0 if the edge is held out for cross-validation.
#'   - \code{row_names} Names of the nodes corresponding to the rows in the original bipartite matrix, representing the first type of nodes.
#'   - \code{col_names} Names of the nodes corresponding to the columns in the original bipartite matrix, representing the second type of nodes.
#'   - \code{edge_type} Specifies whether each edge/link is 'documented' (present in the original data) or 'held_out' (excluded for cross-validation).
#'   - \code{n} The number of times each node pair (i, j) was measured. In this model context, it is assumed that every node pair (i, j) was measured exactly once (\( n_{ij} = 1 \)).
#'   - \code{x} The observed values for each node pair (i, j), where \( x_{ij} \in \{0, 1\} \) represents the reported matrix \( A \). Here, \( x_{ij} = 1 \) indicates that there is an edge/link between nodes i and j, while \( x_{ij} = 0 \) indicates that there is no edge/link.
#'
#' @details
#' - The \code{data} parameter should be a \code{matrix} or \code{data.frame} where each entry represents the presence of interaction/link between the nodes represented by rows and columns.
#' - The \code{folds} parameter, if provided, should be a \code{matrix} containing pairs of indices (row, column) representing held-out edges/links during cross-validation.
#'   It typically contains columns:
#'   \describe{
#'     \item{row}{Indices of rows where edges/links are held out.}
#'     \item{col}{Indices of columns where edges/links are held out.}
#'     \item{gr}{Group identifier indicating which cross-validation fold the edge/link pair belongs to.}
#'   }
#' - If \code{folds} is \code{NULL}, \code{create_cv_folds} is called internally to generate folds based on \code{data}. Each fold ensures that during cross-validation, every column in the input matrix \code{data} has at least \code{min_per_col} non-zero entries, and every row has at least \code{min_per_row} non-zero entries.
#' - The \code{edgelists} are created using the \code{hsbm_edgelist} function for each fold, incorporating information on the number of observations per row and column.
#'
#' @examples
#' # Load example data
#' data(dat, package = "sabinaHSBM")
#'
#' # Prepare input data
#' myInput <- hsbm.input(data = dat, 
#'                       n_folds = 10) 
#'
#' # Summary of dataset
#' summary(myInput)
#'
#' # Plot Reconstructed matrix 
#' plot_interaction_matrix(myInput$data, order_mat = FALSE)
#'
#' @export
hsbm.input <- function(data, n_folds = 5, folds = NULL,
                       min_per_col = 2, min_per_row = 2, 
                       add_spurious = FALSE, no_heldout = FALSE){

    if(!is.matrix(data)){
        warning("Argument data should be of type matrix. Converting to matrix.")
        data <- as.matrix(data)
    }
    # data checks, rowsums == 0, etc
    data <- data[rowSums(data) != 0, colSums(data) != 0]

    edgelists <- list()
    if(n_folds == 1){
        edgelists[[1]] <- hsbm_edgelist(data, folds = NULL, 
                                        add_spurious = add_spurious,
                                        no_heldout = no_heldout)
        value <- list(data = data, folds = folds, edgelists = edgelists)
        attr(value, "class") <- "hsbm.input"
        return(value)
    }

    if(is.null(folds)){
        folds <- create_cv_folds(com = data, n = n_folds,
                                 min_per_col = min_per_col, 
                                 min_per_row = min_per_row)
    }

    nr_folds <- length(unique(folds[, 3]))
    for(i in 1:nr_folds){
        edgelists[[i]] <- hsbm_edgelist(data, folds, fold_id = i, 
                                        add_spurious = add_spurious,
                                        no_heldout = no_heldout)

    }

    value <- list(data = data, folds = folds, edgelists = edgelists)

    attr(value, "class") <- "hsbm.input"

    return(value)

}

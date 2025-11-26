#' @name hsbm.input
#'
#' @title Prepare input data for Hierarchical Stochastic Block Model (HSBM) analyses.
#'
#' @description This function prepares the input data necessary for running the \bold{Hierarchical Stochastic Block Model (HSBM)} analysis by creating cross-validation folds and corresponding edgelists.
#'
#' @param data A binary \code{matrix} containing the input data where rows and columns represent nodes and where the entries (i, j) represent the interactions, edges or links: 1 signifies a relationship/link, and 0 signifies no relationship/no link. The matrix should be an undirected network.
#' @param n_folds (\emph{optional, default} \code{10}) \cr
#' A \code{numeric} value specifying the number of cross-validation folds to generate. If only 1 fold is specified, held-out links will not be computed. This is not recommended, as the main objective of the package is to predict missing links through cross-validation.
#' @param folds (\emph{optional, default} \code{NULL}) \cr
#' An optional \code{matrix} of precomputed cross-validation folds for the input data \code{data}. If \code{NULL} (default), folds are generated internally using \code{create_cv_folds()}.
#' @param min_per_row A code{numeric} specifying the minimum degree a row node must have for any of its
#'   links to be considered for removal. Default is 2.
#' @param min_per_col A code {numeric} specifying the minimum degree a column node must have for any of its
#'   links to be considered for removal. Default is 2.
#' @param max_held_per_fold (\emph{optional, default} \code{NULL}) \cr
#'   An optional code{numeric} specifying the maximum number
#'   of held-out links to place in a fold. If a fold grows larger, its links will be
#'   randomly subsampled to this size. Default is NULL (no limit).
#' @param no_heldout (\emph{optional, default} \code{FALSE}) \cr
#'   A \code{logical} value, that when \code{TRUE} does not encode removed links in
#'   cross-validation as held-outs.
#' @param is_bipartite (\emph{optional, default} \code{TRUE}) \cr
#' A \code{logical} indicating if the \code{data} argument is a bipartite matrix (i.e. rows and cols correspond to different types of nodes).
#' If \code{FALSE} the matrix is taken as unipartite.
#' If the matrix is unipartite, only the lower triangular (diagonal included) is used to generate the edgelist and subsequent predictions.
#'
#' @return
#' An object of class \code{hsbm.input} containing the organized input data, the cross-validation fold assignment for each held-out edge/link, and a list of the corresponding edge lists generated for each fold for HSBM analysis:
#' - \code{$data} The binary \code{matrix} of input data.
#' - \code{$folds} A \code{matrix} of cross-validation fold assignments for each held-out edge/link.
#' - \code{$edgelist} A \code{list} of edge lists generated for each fold. Each edge list is a \code{data.frame} with the following columns:
#'   - \code{v1} The index of the first (in rows) type of nodes.
#'   - \code{v2} The index of the second (in columns) type of nodes.
#'   - \code{value} The measurement value for the edge. It is 1 if there is a edge/link between the node pair, and 0 if the edge is held out for cross-validation.
#'   - \code{row_names} Names of the nodes corresponding to the rows in the original matrix.
#'   - \code{col_names} Names of the nodes corresponding to the columns in the original matrix.
#'   - \code{edge_type} Specifies whether each edge/link is 'documented' (present in the original data) or 'held_out' (excluded for cross-validation).
#'   - \code{n} The number of times each node pair (i, j) was measured. In this model context, it is assumed that every node pair (i, j) was measured exactly once \eqn{n_{ij} = 1}.
#'   - \code{x} The observed values for each node pair (i, j), where \eqn{x_{ij} \in \{0, 1\}} represents the reported matrix \eqn{D}. Here, \eqn{x_{ij} = 1} indicates that there is an edge/link between nodes i and j, while \eqn{x_{ij} = 0} indicates that there is no edge/link.
#'
#' @details
#' - The \code{data} parameter should be a \code{matrix} or \code{data.frame} where each entry represents the presence of interaction/link between the nodes represented by rows and columns.
#' - The \code{folds} parameter, if provided, should be a \code{matrix} containing three columns. The first two columns should be pairs of indices (row, column) representing held-out edges/links during cross-validation.
#'   It typically contains columns:
#'   \describe{
#'     \item{row}{Indices of rows where edges/links are held out.}
#'     \item{col}{Indices of columns where edges/links are held out.}
#'     \item{fold_id}{Fold identifier indicating which cross-validation fold the edge/link pair belongs to.}
#'   }
#' - If \code{folds} is \code{NULL}, \code{create_cv_folds} is called internally to generate folds based on \code{data}. Each fold ensures that during cross-validation, every column in the input matrix \code{data} has at least \code{min_per_col} non-zero entries, and every row has at least \code{min_per_row} non-zero entries.
#' - The \code{edgelists} are created for each fold, incorporating information on the number of observations per row and column.
#' - The \code{is_bipartite} is a logical indicator for bipartite/unipartite matrices.
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
hsbm.input <- function(data, n_folds = 10, folds = NULL,
                       min_per_col = 2, min_per_row = 2,
                       max_held_per_fold = NULL,
                       no_heldout = FALSE,
                       is_bipartite = TRUE){

    if(!is.matrix(data)){
        warning("Argument data should be of type matrix. Converting to matrix.")
        data <- as.matrix(data)
    }
    # Remove possible upper diagonal entries to avoid edge multiplicities
    if(!is_bipartite){
        deg <- rowSums(data) + colSums(data) - diag(data)
        data <- data[deg > 0, deg > 0]
        data[upper.tri(data)] <- 0
    }else{
        data <- data[rowSums(data) != 0, colSums(data) != 0]
    }

    if(!(is.matrix(data))) stop("Too many zero columns or rows in adj_mat.")

    if(is.null(colnames(data))) colnames(data) <- paste0("col", 1:ncol(data))
    if(is.null(rownames(data))) rownames(data) <- paste0("row", 1:nrow(data))

    edgelists <- list()
    if(n_folds == 1){
        edgelists[[1]] <- hsbm_edgelist(data, folds = NULL, 
                                        no_heldout = no_heldout,
                                        is_bipartite = is_bipartite)
        value <- list(data = data, folds = folds, edgelists = edgelists)
        attr(value, "class") <- "hsbm.input"
        return(value)
    }

    if(is.null(folds)){
        folds <- create_cv_folds(com = data, n_folds = n_folds,
                                 min_per_col = min_per_col, 
                                 min_per_row = min_per_row,
                                 max_held_per_fold = max_held_per_fold,
                                 is_bipartite = is_bipartite)
    }

    nr_folds <- length(unique(folds[, 3]))
    for(i in 1:nr_folds){
        edgelists[[i]] <- hsbm_edgelist(data, folds, fold_id = i, 
                                        no_heldout = no_heldout,
                                        is_bipartite = is_bipartite)

    }

    value <- list(data = data, folds = folds, edgelists = edgelists,
                  is_bipartite = is_bipartite)

    attr(value, "class") <- "hsbm.input"

    return(value)

}

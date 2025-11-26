#' Example HSBM input object
#'
#' An example of an object of class \code{hsbm.input}, precomputed using the \code{hsbm.input()} function. This object is included to facilitate the demonstration of functionality in examples and vignettes without requiring computational overhead.
#'
#' @format A \code{list} with the following components:
#' \describe{
#'   \item{data}{A \code{matrix} of binary values representing links in an original bipartite network.}
#'   \item{folds}{A \code{matrix} specifying cross-validation fold assignments for held-out edges/links.}
#'   \item{edgelists}{A \code{list} of edge lists generated for each fold, detailing the held-out links for validation.}
#'   \item{is_bipartite}{A code{logical} indicator for bipartite/unipartite matrices.}
#' }
#'
#' @source Generated using the \code{hsbm.input()} function in the \code{sabinaHSBM} package.
#' @usage data(myInput)
"myInput"

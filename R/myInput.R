##' An \code{R} object of class \code{hsbm.input}, resulting from the application of the \code{hsbm.input()} function.
#'
#' @format ## `myInput`
#' A \code{list} with the following components:
#' \describe{
#'   \item{data}{A \code{matrix} of binary values representing links in a original bipartite network.}
#'   \item{folds}{A \code{matrix} of cross-validation fold assignments for each held-out edge/link.}
#'   \item{edgelist}{A \code{list} of edge lists generated for each fold.}
#'   \item{method}{The method used for the HSBM analysis, here "binary_classifier".}
#'   \item{iter}{The number of iterations for the HSBM analysis.}
#' }
#'
#' @source Data is of own creation.
"myInput"

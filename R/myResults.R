##' An \code{R} object of class \code{hsbm.predict}, resulting from the application of the \code{get_hsbm_results()} function.
#'
#' @format A \code{list} with the following components:
#' \describe{
#'   \item{data}{A \code{matrix} representing the binary bipartite input data.}
#'   \item{folds}{A \code{matrix} detailing the cross-validation fold assignments for each held-out edge/link.}
#'   \item{edgelist}{A \code{list} of different \code{data.frame} containing the edge lists generated for each fold.}
#'   \item{predictions}{A \code{list} with:
#'     \describe{
#'       \item{probs}{A \code{list} of different \code{data.frame} with the predicted probabilities for the edges/links in the corresponding edge list.}
#'       \item{groups}{A \code{list} of different \code{data.frame} containing the group assignments for each node in the network for the corresponding edge list.}
#'       \item{res_folds}{A \code{list} of different \code{data.frame} containing the results of the predictions for each fold, including predicted probabilities and edge types.}
#'       \item{res_averaged}{A \code{data.frame} summarizing the averaged results across all folds.}
#'     }}
#' }
#'
#' @source Data is of own creation.
"myResults"

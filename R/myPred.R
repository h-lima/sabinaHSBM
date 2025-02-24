#' Example HSBM prediction object
#'
#' An example of an object of class \code{hsbm.predict}, precomputed using the \code{hsbm.predict()} function. This object is included to facilitate the demonstration of functionality in examples and vignettes without requiring computational overhead.
#'
#' @format A \code{list} with the following components:
#' \describe{
#'   \item{data}{A \code{matrix} representing the binary bipartite input data.}
#'   \item{folds}{A \code{matrix} specifying cross-validation fold assignments for held-out edges/links.}
#'   \item{probs}{Predicted probabilities for edges/links, stored as \code{data.frame}s, one for each fold.}
#'   \item{groups}{Group assignments for nodes across hierarchical levels, stored as \code{data.frame}s.}
#'   \item{min_dl}{A \code{list} of minimum description length values for each fold.}
#' }
#'
#' @source Generated using the \code{hsbm.predict()} function in the \code{sabinaHSBM} package.
#' @usage data(myPred)
"myPred"

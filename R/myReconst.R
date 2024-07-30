##' An \code{R} object of class \code{hsbm.reconstructed}, resulting from the application of the \code{hsbm.reconstructed()} function.
#'
#' @format A \code{list} with the following components:
#' \describe{
#'   \item{data}{A \code{matrix} representing the binary bipartite input data.}
#'   \item{pred_mats}{A \code{list} of matrices where ...........#@@JMB temrinar....}
#'   \item{reconstructed_bin_mats}{A \code{list} of binary matrices representing the reconstructed data for each fold based on the applied threshold.}
#'   \item{reconstructed_stats}{A \code{list} of numeric vectors where each vector contains various statistics for the reconstructed data of each fold.}
#'   \item{tb}{A \code{data.frame} with 10 rows and 16 columns summarizing performance metrics across folds.}
#'   \item{new_mat}{A \code{numeric} matrix representing the final reconstructed binary matrix.}
#'   \item{threshold}{A \code{character} string indicating the thresholding method used for reconstruction.}
#' }
#'
#' @source Data is of own creation.
"myReconst"

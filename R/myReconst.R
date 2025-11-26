#' Example HSBM reconstructed object
#'
#' An example of an object of class \code{hsbm.reconstructed}, precomputed using the \code{hsbm.reconstructed()} function.
#' This object is included to facilitate the demonstration of functionality
#' in examples and vignettes without requiring computational overhead.
#'
#' @format A \code{list} with the following components:
#' \describe{
#'   \item{data}{A \code{matrix} representing the binary input data.}
#'   \item{pred_mats}{A \code{list} of matrices containing predicted probabilities for links in each cross-validation fold.}
#'   \item{reconstructed_df}{A \code{list} with the following components:
#'       \describe{
#'          \item{res_folds}{A \code{list} of binary matrices representing reconstructed data for each fold, stored as \code{data.frame}s.}
#'          \item{res_averaged}{A \code{data.frame} summarizing predictions averaged across all folds}
#'      }
#'   }
#'   \item{stats}{A \code{data.frame} summarizing performance metrics across folds.}
#'   \item{new_mat}{A \code{matrix} representing the final reconstructed binary matrix.}
#'   \item{threshold}{A \code{character} string indicating the thresholding method used.}
#' }
#'
#' @source Generated using the \code{hsbm.reconstructed()} function in the \code{sabinaHSBM} package.
#' @usage data(myReconst)
"myReconst"

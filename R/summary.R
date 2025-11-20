#' @title Summary for hsbm.input Objects
#' @description Provides a summary for objects of class \code{hsbm.input}.
#' @param object An object of class \code{hsbm.input}.
#' @param ... Additional arguments (not used).
#' @return A data frame containing:
#' \describe{
#'   \item{n_rows}{Number of rows (first set of nodes).}
#'   \item{n_cols}{Number of columns (second set of nodes).}
#'   \item{n_links}{Number of observed links in the binary input matrix.}
#'   \item{rows_single_link}{Number of rows with only one link.}
#'   \item{cols_single_link}{Number of columns with only one link.}
#'   \item{possible_links}{Total number of possible links in the input matrix.}
#' }
#' @seealso \code{\link{hsbm.input}}
#' @export
#' @method summary hsbm.input
summary.hsbm.input <- function(object, ...){
    mat <- as.matrix(object$data)

    summary_mat <- data.frame(n_rows = nrow(mat),
                            n_cols = ncol(mat),
                            n_links = sum(mat),
                            rows_single_link = sum(rowSums(mat) == 1),
                            cols_single_link = sum(colSums(mat) == 1),
                            possible_links = length(mat))
    return(summary_mat)
}

#' @title Summary for hsbm.reconstructed Objects
#' @description Provides a summary for objects of class \code{hsbm.reconstructed}.
#' @param object An object of class \code{hsbm.reconstructed}.
#' @param ... Additional arguments (not used).
#' @return A data frame containing the following summary statistics:
#' \describe{
#'   \item{obs_links}{Number of observed links in the original input matrix.}
#'   \item{unobs_links}{Number of unobserved links in the original input matrix.}
#'   \item{pred_links}{Number of predicted links in the reconstructed binary matrix.}
#'   \item{kept_links}{Number of links retained (unchanged) in the reconstructed matrix.}
#'   \item{spurious_links}{Number of links removed (false positives) in the reconstructed matrix.}
#'   \item{missing_links}{Number of new links added (false negatives) in the reconstructed matrix.}
#'   \item{mean_RLF}{Mean Recovery Link Fraction (proportion of held-out links correctly predicted as missing links) across all folds.}
#'   \item{mean_auc}{Mean AUC across all folds.}
#'   \item{mean_aucpr}{Mean AUC-PR across all folds.}
#'   \item{mean_yPRC}{Mean Precision-Recall Curve baseline across all folds.}
#'   \item{mean_prec}{Mean precision across all folds.}
#'   \item{mean_sens}{Mean sensitivity across all folds.}
#'   \item{mean_spec}{Mean specificity across all folds.}
#'   \item{mean_ACC}{Mean accuracy across all folds.}
#'   \item{mean_ERR}{Mean error rate across all folds.}
#'   \item{mean_tss}{Mean True Skill Statistic (TSS) across all folds.}
#' }
#' @seealso \code{\link{hsbm.reconstructed}}
#' @export
#' @method summary hsbm.reconstructed
summary.hsbm.reconstructed <- function(object, ...){
    matrix1 <- as.matrix(object$data)
    matrix2 <- as.matrix(object$new_mat)

   # Number of observed links
   obs_links <- sum(matrix1, na.rm = TRUE)
   # Number of unobserved links
   unobs_links <- length(matrix1) - obs_links
   # Number of predicted links
   pred_links <- sum(matrix2, na.rm = TRUE)
   # Number of links kept the same
   kept_links <- sum(matrix1 == matrix2 & matrix1 == 1, na.rm = TRUE)
   # Number of lost links in matrix2
   spurious_links <- sum(matrix1 == 1 & matrix2 == 0, na.rm = TRUE)
   # Number of new links in matrix2
   missing_links <- sum(matrix1 == 0 & matrix2 == 1, na.rm = TRUE)

   summary_mat <- data.frame(
                       obs_links = obs_links,
                       unobs_links = unobs_links,
                       pred_links = pred_links,
                       kept_links = kept_links,
                       spurious_links = spurious_links,
                       missing_links = missing_links
   )

   summary_eval <- data.frame(
                       mean_RLF = mean(object$stats$pred_held_ones),
                       mean_auc = mean(object$stats$auc),
                       mean_aucpr = mean(object$stats$aucpr),
                       mean_yPRC = mean(object$stats$yPRC),
                       mean_prec = mean(object$stats$precision),
                       mean_sens = mean(object$stats$sens),
                       mean_spec = mean(object$stats$spec),
                       mean_ACC = mean(object$stats$ACC),
                       mean_ERR = mean(object$stats$ERR),
                       mean_tss = mean(object$stats$tss)
   )

  return(list(
        `Reconstructed Network Metrics` = summary_mat,
        `Evaluation Metrics` = summary_eval
    ))
}

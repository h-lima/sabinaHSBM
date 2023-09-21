#' @export
summary.hsbm.input <- function(hsbm_input){
    mat <- as.matrix(hsbm_input$data)

    summary_mat <- data.frame(n_rows = nrow(mat),
                            n_cols = ncol(mat),
                            n_links = sum(mat),
                            rows_single_link = sum(rowSums(mat) == 1),
                            cols_single_link = sum(colSums(mat) == 1),
                            possible_links = length(mat))
    return(summary_mat)
}

#' @export
summary.hsbm.predict <- function(hsbm_predict){
    mat <- as.matrix(hsbm_predict$data)

    summary_mat <- data.frame(n_rows = nrow(mat),
                            n_cols = ncol(mat),
                            n_links = sum(mat),
                            rows_single_link = sum(rowSums(mat) == 1),
                            cols_single_link = sum(colSums(mat) == 1),
                            possible_links = length(mat))

    return(summary_mat)
}

#' @export 								
summary.hsbm.reconstructed <- function(hsbm_new){
    matrix1 <- as.matrix(hsbm_new$data)
    matrix2 <- as.matrix(hsbm_new$new_mat)

   # Number of observed links
   obs_links <- sum(matrix1)
   # Number of unobserved links
   unobs_links <- length(matrix1) - obs_links
   # Number of predicted links
   pred_links <- sum(matrix2)
   # Number of links kept the same
   kept_links <- sum(matrix1 == matrix2 & matrix1 == 1)
   # Number of lost links in matrix2
   spurious_links <- sum(matrix1 == 1 & matrix2 == 0)
   # Number of new links in matrix2
   missing_links <- sum(matrix1 == 0 & matrix2 == 1)

   summary_mat <- data.frame(obs_links = obs_links,
                       unobs_links = unobs_links,
                       pred_links = pred_links,
                       kept_links = kept_links,
                       spurious_links = spurious_links,
                       missing_links = missing_links, 
                       mean_auc = mean(hsbm_new$tb$auc),
                       mean_aucpr = mean(hsbm_new$tb$aucpr),
                       mean_prec = mean(hsbm_new$tb$precision),
		               mean_sens = mean(hsbm_new$tb$sens),
		               mean_spec = mean(hsbm_new$tb$spec),
		               mean_ACC = mean(hsbm_new$tb$ACC),
		               mean_ERR = mean(hsbm_new$tb$ERR),
		               mean_tss = mean(hsbm_new$tb$tss))

  return(summary_mat)
}

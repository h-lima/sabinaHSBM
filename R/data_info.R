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

#' @export
hsbm.input <- function(data, folds = NULL, n_folds = 5, method = "binary_classifier", iter = 10000,
                       add_n_x = TRUE, min_per_col = 2, min_per_row = 2){
    # data checks, rowsums == 0, etc
    data <- data[rowSums(data) != 0, colSums(data) != 0]

    if(is.null(folds)){
        folds <- create_cv_folds(Z = data, n = n_folds,
                                 min_per_col = min_per_col, min_per_row = min_per_row)
    }

    edgelists <- list()

    for(i in 1:n_folds){
        edgelists[[i]] <- hsbm_edgelist(data, folds, fold_id = i, add_n_x = add_n_x)

    }

    value <- list(data = data, folds = folds, edgelist = edgelists,
                  method = method, iter = iter)

    attr(value, "class") <- "hsbm.input"

    return(value)

}

hsbm_edgelist <- function(adj_mat, folds, fold_id = NULL, no_heldout = FALSE,
                          is_bipartite = TRUE){
    col_names <- c("row_names", "col_names", "value")

    if(!(is.matrix(adj_mat))) stop("adj_mat argument must be of type 'matrix'.")

    adj_mat_train <- adj_mat
    if(!is_bipartite) adj_mat_train[upper.tri(adj_mat_train)] <- 0
    if(!is.null(fold_id)){
        folds <- as.matrix(folds[which(folds[, 3] == fold_id), ])
        adj_mat_train[folds[, c('row', 'col')]] <- 0
    }else if(!is.null(folds)){
        adj_mat_train[folds[, c('row', 'col')]] <- 0
    }

    long_mat <- reshape2::melt(adj_mat_train)
    colnames(long_mat) <- col_names
    cond <- (long_mat$value == 1)
    long_mat <- dplyr::filter(long_mat, cond)
    long_mat$edge_type <- "documented"

    if(!is.null(folds)){
        colnames(folds) <- col_names
        folds[, 'value'] <- 0
        folds <- as.data.frame(folds)
        folds$edge_type <- "held_out"
        folds$row_names <- rownames(adj_mat)[folds$row_names]
        folds$col_names <- colnames(adj_mat)[folds$col_names]

        long_mat <- rbind(long_mat, folds)
    }

    long_mat$n <- 1
    long_mat$x <- 1
    if(!no_heldout){
        long_mat$x <- ifelse(long_mat$edge_type == "held_out", 0, long_mat$x)
    }else{
        long_mat$edge_type <- ifelse(long_mat$edge_type == "held_out", "documented",
                                     long_mat$edge_type)
    }

    long_mat$v1 <- as.numeric(long_mat[[col_names[1]]]) - 1
    if(is_bipartite){
        long_mat$v2 <- as.numeric(long_mat[[col_names[2]]]) + max(long_mat$v1)
    }else{
        long_mat$v2 <- as.numeric(long_mat[[col_names[2]]]) - 1
    }

    long_mat <- dplyr::select(long_mat, "v1", "v2", "value", col_names[1],
                              col_names[2], "edge_type", "n", "x")

    return(long_mat)
}


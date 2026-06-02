hsbm_edgelist <- function(adj_mat, folds, fold_id = NULL, no_heldout = FALSE,
                          is_bipartite = TRUE, custom_nx = NULL){
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
    long_mat$is_custom <- FALSE
    
    if(!is.null(custom_nx)){
        long_mat <- dplyr::left_join(long_mat, custom_nx, 
                                     by = c("row_names" = "row", 
                                            "col_names" = "col"),
                                     suffix = c("", "_custom"))
        
        long_mat$is_custom <- !is.na(long_mat$n_custom) | !is.na(long_mat$x_custom)
        
        long_mat$n <- ifelse(!is.na(long_mat$n_custom), 
                             long_mat$n_custom, long_mat$n)
        long_mat$x <- ifelse(!is.na(long_mat$x_custom), 
                             long_mat$x_custom, long_mat$x)
    }
    if(!no_heldout){
        long_mat$x <- ifelse(long_mat$edge_type == "held_out", 
                             0, long_mat$x)
        long_mat$n <- ifelse(long_mat$edge_type == "held_out", 
                             1, long_mat$n)
    }else{
        long_mat$edge_type <- ifelse(long_mat$edge_type == "held_out", 
                                     "documented", long_mat$edge_type)
    }

    long_mat$v1 <- match(long_mat[[col_names[1]]], rownames(adj_mat)) - 1
    if(is_bipartite){
        long_mat$v2 <- match(long_mat[[col_names[2]]], 
                             colnames(adj_mat)) + max(long_mat$v1, na.rm = TRUE)
    }else{
        long_mat$v2 <- match(long_mat[[col_names[2]]], 
                             colnames(adj_mat)) - 1
    }
    long_mat <- dplyr::select(long_mat, "v1", "v2", "value", col_names[1],
                              col_names[2], "edge_type", "n", "x")

    return(long_mat)
}

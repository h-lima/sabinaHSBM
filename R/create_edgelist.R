hsbm_edgelist <- function(adj_mat, folds, fold_id = NULL, add_spurious = FALSE){
    col_names <- c("row_names", "col_names", "value")

    if(!(is.matrix(adj_mat))) stop("adj_mat argument must be of type 'matrix'.")
    # check values 0 and 1
    # check any string
    # check names

    # Remove empty columns and rows
    adj_mat <- adj_mat[rowSums(adj_mat) != 0,
                       colSums(adj_mat) != 0]
    if(!(is.matrix(adj_mat))) stop("Too many zero columns or rows in adj_mat.")

    if(is.null(colnames(adj_mat))) colnames(adj_mat) <- paste0("col", 1:ncol(adj_mat))
    if(is.null(rownames(adj_mat))) rownames(adj_mat) <- paste0("row", 1:nrow(adj_mat))

    if(!is.null(fold_id)){
        folds <- as.matrix(folds[which(folds[, 'gr'] == fold_id), ])
        adj_mat_train <- adj_mat
        adj_mat_train[folds[, c('row', 'col')]] <- 0
    }

    long_mat <- reshape2::melt(adj_mat_train)
    colnames(long_mat) <- col_names
    long_mat <- dplyr::filter(long_mat, value == 1)
    long_mat$edge_type <- "documented"

    colnames(folds) <- col_names
    folds[, 'value'] <- 0
    folds <- as.data.frame(folds)
    folds$edge_type <- "held_out"
    folds$row_names <- rownames(adj_mat)[folds$row_names]
    folds$col_names <- colnames(adj_mat)[folds$col_names]

    long_mat <- rbind(long_mat, folds)

    long_mat$n <- 1
    long_mat$x <- 1
    long_mat[long_mat$edge_type == "held_out", ]$n <- 1
    long_mat[long_mat$edge_type == "held_out", ]$x <- 0

    long_mat$v1 <- as.numeric(factor(long_mat[[col_names[1]]])) - 1
    long_mat$v2 <- as.numeric(factor(long_mat[[col_names[2]]])) + max(long_mat$v1)

    long_mat <- dplyr::select(long_mat, v1, v2, value, col_names[1], col_names[2], edge_type, n, x)

    if(add_spurious){
        long_mat <- add_spurious_edges(long_mat, nrow(adj_mat), col_names[1:2])
    }

    return(long_mat)
}

add_spurious_edges <- function(long_mat, n_row, col_names, add_n_x = TRUE){

    v1_v2 <- paste0(long_mat$v1, "_", long_mat$v2)
    spurious_v1 <- sample(unique(long_mat$v1),
                          size = length(which(long_mat$edge_type == "held_out")),
                          replace = TRUE)
    spurious_v2 <- sample(unique(long_mat$v2),
                          size = length(which(long_mat$edge_type == "held_out")),
                          replace = TRUE)
    spurious_v1_v2 <- paste0(spurious_v1, "_", spurious_v2)
    for(i in 1:length(spurious_v1_v2)){
        if(!(spurious_v1_v2[i] %in% v1_v2)){
            next
        }
        while(spurious_v1_v2[i] %in% v1_v2){
            spurious_v1_v2[i] <- paste0(sample(unique(long_mat$v1), size = 1),
                                        "_",
                                        sample(unique(long_mat$v2), size = 1))
        }
    }
    spurious_mat <- stringr::str_match(spurious_v1_v2, "(\\d+)_(\\d+)")[, 2:3]
    x_names <- paste0("row", as.numeric(spurious_mat[, 1]) - 1)
    y_names <- paste0("col", as.numeric(spurious_mat[, 2]) - n_row + 1)
    colnames(spurious_mat) <- c("v1", "v2")
    spurious_mat <- cbind(spurious_mat, value = 1, x_names,
                          y_names, edge_type = "spurious_edge")
    colnames(spurious_mat)[4:5] <- col_names
    if(add_n_x){
        spurious_mat <- cbind(spurious_mat, n = 1, x = 0)
    }

    long_mat <- rbind(long_mat, spurious_mat)
    long_mat <- dplyr::mutate_if(long_mat, is.factor, as.character)
    long_mat <- dplyr::mutate_at(long_mat, c("v1", "v2", "value", "n", "x"),
                                 as.numeric)

    return(long_mat)

}

#' @export
create_cv_folds <-function(com, n = 10, min_per_col = 2, min_per_row = 2,
                           get_zeros = FALSE){

    if(!(is.matrix(com))) stop("com argument must be of type matrix")
    if(max(range(com))>1) com[com>0]<-1
    pairs = which(com==1, arr.ind=T)
    colnames(pairs)<-c('row', 'col')

    colm = colSums(com)
    rowm <- rowSums(com)

    # Exclude pairs that cannot be removed
    pairs_dim <- 'col'
    has_mins <- TRUE
    while(has_mins){
        min_val <- get(paste0("min_per_", pairs_dim))
        pairs <- rm_pairs_rows(pairs, row_or_col = pairs_dim,
                               min_val = min_val)
        if(pairs_dim == 'row') pairs_dim <- 'col' else pairs_dim <- 'row'
        freq_count <- data.frame(table(pairs[, pairs_dim]))
        has_mins <- sum(freq_count$Freq < min_val) > 0
    }

    size = floor(nrow(pairs)/n)
    gr = rep(size, n)
    if(nrow(pairs) %% size!=0) gr[n] =  gr[n] + nrow(pairs) %% size

    # Distribute pairs over folds
    folds_lst <- create_folds(com, pairs, gr, min_per_col, min_per_row)

    # Correct possible cases when distribution leads to zero dim sums
    held_res <- NULL
    stop_while <- 0
    while(any(has_0_sum_folds(folds_lst$com, folds_lst$held_pairs, n, 'col')) ||
          any(has_0_sum_folds(folds_lst$com, folds_lst$held_pairs, n, 'row'))){

        if(stop_while > 1000){
            stop("Failed to redistribute 0 sum folds.",
                 "Run create_cv_folds again and if persists file an issue in github.")
        }

        fold_to_correct <- which(has_0_sum_folds(folds_lst$com, folds_lst$held_pairs, n, 'col'))[1]
        held_res <- distribute_zero_sum_folds(folds_lst$com, folds_lst$held_pairs,
                                              folds_lst$pairs_sums, pairs_dim = "col",
                                              min_per_col = min_per_col, min_per_row = min_per_row,
                                              fold_i = fold_to_correct, n = n)
        folds_lst$held_pairs <- held_res
        if(any(has_0_sum_folds(folds_lst$com,folds_lst$held_pairs, n, 'row'))){
            fold_to_correct <- which(has_0_sum_folds(folds_lst$com, folds_lst$held_pairs, n, 'row'))[1]
            held_res <- distribute_zero_sum_folds(folds_lst$com, folds_lst$held_pairs,
                                                  folds_lst$pairs_sums, pairs_dim = "row",
                                                  min_per_col = min_per_col, min_per_row = min_per_row,
                                                  fold_i = fold_to_correct, n = n)
        }

        folds_lst$held_pairs <- held_res

        stop_while <- stop_while + 1

    }

    if(any(has_0_sum_folds(folds_lst$com, folds_lst$held_pairs, n, 'col'))){
        stop("Zero sum column detected after removing fold pairs.")
    }else if(any(has_0_sum_folds(folds_lst$com, folds_lst$held_pairs, n, 'row'))){
        stop("Zero sum row detected after removing fold pairs.")
    }

    held_pairs <- folds_lst$held_pairs

    print(sprintf("Actual cross-validation rate is %0.3f" , table(held_pairs[,'gr'])/sum(1*(com>0))))

    return(held_pairs[order(held_pairs[,'gr']),])

}

rm_pairs_rows <- function(pairs, row_or_col, min_val){
    i_count <- data.frame(table(pairs[, row_or_col]))
    colnames(i_count) <- c("i", "freq")
    rm_i <- i_count$i[i_count$freq < min_val]
    if(any(pairs[, row_or_col] %in% rm_i)){
        pairs <- pairs[-which(pairs[, row_or_col] %in% rm_i), ]
    }

    return(pairs)
}

get_pairs_sums <- function(pairs){
    row_count <- data.frame(table(pairs[, 'row']))
    colnames(row_count) <- c("row", "rowsums")
    row_count$row <- as.numeric(levels(row_count$row))
    pairs_sums <- dplyr::left_join(as.data.frame(pairs), row_count, by = "row")
    col_count <- data.frame(table(pairs[, 'col']))
    colnames(col_count) <- c("col", "colsums")
    col_count$col <- as.numeric(levels(col_count$col))
    pairs_sums <- as.matrix(dplyr::left_join(pairs_sums, col_count, by = "col"))

    return(pairs_sums)
}

subtract_pairs_sums <- function(pairs_sums, pairs_dim, dim_i, min_val){
    forbidden = c()
    if(pairs_dim == "row"){
        j <- 3
    }else if(pairs_dim == "col"){
        j <- 4
    }

    if(any(pairs_sums[, pairs_dim] == dim_i)){
        minus_s <- which(pairs_sums[, pairs_dim] == dim_i)
        pairs_sums[minus_s, j] <- pairs_sums[minus_s, j] - 1
        if(unique(pairs_sums[minus_s, j]) < min_val){
            forbidden <- c(forbidden, unique(pairs_sums[minus_s, pairs_dim]))
        }
    }
    return(list(pairs_sums = pairs_sums, forbidden = forbidden))
}

has_0_sum_folds <- function(com, held_pairs, n, pairs_dim){
    test <- logical(n)
    if(pairs_dim == "row") sum_func <- rowSums else sum_func <- colSums
    for(i in 1:n){
        com_test <- com
        row_col <- held_pairs[held_pairs[,3] == i, c(1,2)]
        com_test[row_col] <- 0
        test[i] <- any(sum_func(com_test) == 0)
    }

    return(test)
}

create_folds <- function(com, pairs, gr, min_per_col, min_per_row){
    fold_sz <- gr
    n <- length(gr)
    pairs_sums <- get_pairs_sums(pairs)
    held_pairs <- matrix(nrow = nrow(pairs_sums), ncol = 3)
    colnames(held_pairs) <- c("row", "col", "gr")
    held_pairs_i <- 1
    for(i in 1:(n - 1)){
        forbidden_rows <- c()
        forbidden_cols <- c()
        pairs_sums <- pairs_sums[pairs_sums[, 'rowsums'] != 0, ]
        pairs_sums <- pairs_sums[pairs_sums[, 'colsums'] != 0, ]
        while(fold_sz[i] > 0){
            held_i <- sample.int(nrow(pairs_sums), 1)
            pairs_sums_i <- pairs_sums[held_i, ]
            row_i <- pairs_sums_i[1]
            col_i <- pairs_sums_i[2]
            while(row_i %in% forbidden_rows || col_i %in% forbidden_cols){
                held_i <- sample.int(nrow(pairs_sums), 1)
                pairs_sums_i <- pairs_sums[held_i, ]
                row_i <- pairs_sums_i[1]
                col_i <- pairs_sums_i[2]
            }
            held_pairs[held_pairs_i, ] <- c(row_i, col_i, i)

            pairs_sums <- pairs_sums[-held_i, ]
            pairs_sums_lst <- subtract_pairs_sums(pairs_sums, "row", row_i, min_per_row)
            forbidden_rows <- c(forbidden_rows, pairs_sums_lst$forbidden)
            pairs_sums_lst <- subtract_pairs_sums(pairs_sums_lst$pairs_sums, "col", col_i, min_per_col)
            forbidden_cols <- c(forbidden_cols, pairs_sums_lst$forbidden)
            pairs_sums <- pairs_sums_lst$pairs_sums

            held_pairs_i <- held_pairs_i + 1
            fold_sz[i] <- fold_sz[i] - 1
        }
    }
    # Last fold gets the remaining pairs
    last_fold_i <- (nrow(held_pairs) - fold_sz[n] + 1):nrow(held_pairs)
    held_pairs[last_fold_i, ] <- cbind(pairs_sums[, c(1, 2)], n)

    return(list(com = com, held_pairs = held_pairs, pairs_sums = pairs_sums))
}


# Get zero sum dimensions and distribute to other folds
distribute_zero_sum_folds <- function(com, held_pairs, pairs_sums, pairs_dim,
                                      min_per_col, min_per_row, fold_i, n, debug = FALSE){
    com_test <- com
    row_col <- held_pairs[held_pairs[,3] == fold_i, c(1,2)]
    com_test[row_col] <- 0
    if(pairs_dim == "row") sum_func <- rowSums else sum_func <- colSums
    if(all(sum_func(com_test) != 0) & debug){
        warning("Nothing to do in distribute_zero_sum_folds. Return original held_pairs.")
        return(held_pairs)
    }
    dim_zero <- which(sum_func(com_test) == 0)
    for(dim_j in dim_zero){
        pairs_i <- pairs_sums[pairs_sums[, pairs_dim] == dim_j, , drop = FALSE]
        pairs_rsums <- pairs_i[, 3]
        pairs_csums <- pairs_i[, 4]
        if(any(pairs_rsums < min_per_row)){
            if(all(pairs_rsums < min_per_row)){
                pairs_i <- pairs_i
            }else{
                pairs_i <- pairs_i[-(which(pairs_rsums < min_per_row)), , drop = FALSE]
            }
        }
        pairs_csums <- pairs_i[, 4]
        if(any(pairs_csums < min_per_col)){
            if(all(pairs_csums < min_per_col)){
                pairs_i <- pairs_i
            }else{
                pairs_i <- pairs_i[-(which(pairs_csums < min_per_col)), , drop = FALSE]
            }
        }

        move_row <- sample(nrow(pairs_i), pmax(nrow(pairs_i) - 1, 1))
        if(length(move_row) > length(setdiff(1:n, fold_i))){
            to_fold <- sample(setdiff(1:n, fold_i), length(move_row), replace = TRUE)
        }else{
            to_fold <- sample(setdiff(1:n, fold_i), length(move_row))
        }
        for(i in 1:length(move_row)){
            pairs_move <- c(pairs_i[move_row[i], c(1, 2)], gr = to_fold[i])
            pairs_row <- which(held_pairs[, 1] == pairs_move[1] & held_pairs[, 2] == pairs_move[2])
            held_pairs[pairs_row, 3] <- pairs_move[3]
        }
    }

    return(held_pairs)

}


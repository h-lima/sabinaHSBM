#' Code adapted from HPprediction package
#' @export
create_cv_folds <-function(Z, n= 10, min_per_col = 2, min_per_row = 2,
                           get_zeros = FALSE){

    if(!(is.matrix(Z))) stop("Z argument must be of type matrix")
    ## n-fold cross validation
    if(max(range(Z))>1) Z[Z>0]<-1
    pairs = which(Z==1, arr.ind=T)
    colnames(pairs)<-c('row', 'col')

    colm = colSums(Z)
    rowm <- rowSums(Z)

    cat("nrow pairs first: ", nrow(pairs), '\n')
    cat("sum(colm): ", sum(colm), '\n')
    cat("sum(rowm): ", sum(rowm), '\n')

    # Exclude rows that cannot be removed
    pairs_dim <- 'col'
    freq_count <- data.frame(table(pairs[, pairs_dim]))
    min_val <- get(paste0("min_per_", pairs_dim))
    has_mins <- sum(freq_count$Freq < min_val) > 0
    while(has_mins){
        min_val <- get(paste0("min_per_", pairs_dim))
        pairs <- rm_pairs_rows(pairs, row_or_col = pairs_dim,
                            min_val = min_val)
        cat("nrow pairs after func ", nrow(pairs), '\n')
        if(pairs_dim == 'row') pairs_dim <- 'col' else pairs_dim <- 'row'
        freq_count <- data.frame(table(pairs[, pairs_dim]))
        #print(freq_count)
        #print(sum(freq_count$Freq < min_val))
        has_mins <- sum(freq_count$Freq < min_val) > 0
        #print(has_mins)
    }



    size = floor(nrow(pairs)/n)
    gr = rep(size, n)
    if(nrow(pairs) %% size!=0) gr[n] =  gr[n] + nrow(pairs) %% size
    cat("gr: ", gr, "\n")

    # Distribute pairs over folds
    folds_lst <- create_folds(Z, pairs, gr, min_per_col, min_per_row)

    # Correct possible cases when distribution leads to zero dim sums
    held_res <- NULL
    while(any(has_0_sum_folds(folds_lst$com, folds_lst$held_pairs, n, 'col'))){
        #print("fail col")
        fold_to_correct <- which(has_0_sum_folds(folds_lst$com, folds_lst$held_pairs, n, 'col'))[1]
        held_res <- distribute_last_held(folds_lst$com, folds_lst$held_pairs,
                                         folds_lst$pairs_sums, pairs_dim = "col",
                                         min_per_col = min_per_col, min_per_row = min_per_row,
                                         fold_i = fold_to_correct, n = n)
        #print("done for col")
        if(any(has_0_sum_folds(folds_lst$com, held_res, n, 'row'))){
            #print("fail row")
            fold_to_correct <- which(has_0_sum_folds(folds_lst$com, folds_lst$held_pairs, n, 'row'))[1]
            held_res <- distribute_last_held(folds_lst$com, held_res,
                                             folds_lst$pairs_sums, pairs_dim = "row",
                                             min_per_col = min_per_col, min_per_row = min_per_row,
                                             fold_i = fold_to_correct, n = n)
            #print("done for row")
        }

        folds_lst$held_pairs <- held_res
    }
    if(!is.null(held_res)){
        #print("#########CREATE FOLDS CORRECTED########")
        if(any(has_0_sum_folds(folds_lst$com, held_res, n, 'col'))){
            #print("fail col")
            stop("fail col")
        }else if(any(has_0_sum_folds(folds_lst$com, held_res, n, 'row'))){
            #print("fail row1")
            stop("fail row")
        }
    }else{
        #print("PERFECT CREATE FOLDS")
    }


    held_pairs <- folds_lst$held_pairs

    print(sprintf("Actual cross-validation rate is %0.3f" , table(held_pairs[,'gr'])/sum(1*(Z>0))))

    return(held_pairs[order(held_pairs[,'gr']),])

}

rm_pairs_rows <- function(pairs, row_or_col, min_val){
    i_count <- data.frame(table(pairs[, row_or_col]))
    colnames(i_count) <- c("i", "freq")
    rm_i <- i_count$i[i_count$freq < min_val]
    if(any(pairs[, row_or_col] %in% rm_i)){
        pairs <- pairs[-which(pairs[, row_or_col] %in% rm_i), ]
    }
    cat("number of rm_i: ", length(rm_i), '\n')
    cat("nrow pairs after rm i: ", nrow(pairs), '\n')
    cat("any rm i in pairs: ", any(rm_i %in% pairs[, row_or_col]), '\n')
    return(pairs)
}

get_pairs_sums <- function(pairs){
    row_count <- data.frame(table(pairs[, 'row']))
    colnames(row_count) <- c("row", "rowsums")
    row_count$row <- as.numeric(levels(row_count$row))
    pairs_sums <- left_join(as.data.frame(pairs), row_count, by = "row")
    col_count <- data.frame(table(pairs[, 'col']))
    colnames(col_count) <- c("col", "colsums")
    col_count$col <- as.numeric(levels(col_count$col))
    pairs_sums <- as.matrix(left_join(pairs_sums, col_count, by = "col"))

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
        #print(pairs_sums[minus_s, ])
        pairs_sums[minus_s, j] <- pairs_sums[minus_s, j] - 1
        if(unique(pairs_sums[minus_s, j]) < min_val){
            forbidden <- c(forbidden, unique(pairs_sums[minus_s, pairs_dim]))
        }
        #print(pairs_sums[minus_s, ])
    }
    return(list(pairs_sums = pairs_sums, forbidden = forbidden))
}

has_0_sum_folds <- function(com, held_pairs, n, pairs_dim){
    test <- logical(n)
    if(pairs_dim == "row") sum_func <- rowSums else sum_func <- colSums
    for(i in 1:n){
        Z_test <- com
        row_col <- held_pairs[held_pairs[,3] == i, c(1,2)]
        Z_test[row_col] <- 0
        test[i] <- any(sum_func(Z_test) == 0)
        #cat("\tfold i zero sum FALSE: ", test[i], '\n')
        if(!test[i]){


        }
    }

    return(test)
}

create_folds <- function(Z, pairs, gr, min_per_col, min_per_row){
    fold_sz <- gr
    n <- length(gr)
    pairs_sums <- get_pairs_sums(pairs)
    held_pairs <- matrix(nrow = nrow(pairs_sums), ncol = 3)
    colnames(held_pairs) <- c("row", "col", "gr")
    held_pairs_i <- 1
    cat("nrow pairs sums: ", nrow(pairs_sums), '\n')
    for(i in 1:(n - 1)){
        forbidden_rows <- c()
        forbidden_cols <- c()
        pairs_sums <- pairs_sums[pairs_sums[, 'rowsums'] != 0, ]
        pairs_sums <- pairs_sums[pairs_sums[, 'colsums'] != 0, ]
        cat("nrow(pairs_sums) before fold ", i, ":", nrow(pairs_sums), "\n")
        while(fold_sz[i] > 0){
            held_i <- sample.int(nrow(pairs_sums), 1)
            pairs_sums_i <- pairs_sums[held_i, ]
            row_i <- pairs_sums_i[1]
            col_i <- pairs_sums_i[2]
            while(row_i %in% forbidden_rows || col_i %in% forbidden_cols){
                #print("HERE")
                #print(pairs_sums_i)
                held_i <- sample.int(nrow(pairs_sums), 1)
                pairs_sums_i <- pairs_sums[held_i, ]
                row_i <- pairs_sums_i[1]
                col_i <- pairs_sums_i[2]
            }
            if(is.na(row_i)) stop("na row_i")
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
        cat("nrow(pairs_sums) after fold ", i, ":", nrow(pairs_sums), "\n")
    }
    last_fold_i <- (nrow(held_pairs) - fold_sz[n] + 1):nrow(held_pairs)
    held_pairs[last_fold_i, ] <- cbind(pairs_sums[, c(1, 2)], n)

    if(any(has_0_sum_folds(Z, held_pairs, n, 'col'))){
        print("fail col")

    }
    if(any(has_0_sum_folds(Z, held_pairs, n, 'row'))){
        print("fail row")
    }

    return(list(com = Z, held_pairs = held_pairs, pairs_sums = pairs_sums))
}


distribute_last_held <- function(com, held_pairs, pairs_sums, pairs_dim,
                                 min_per_col, min_per_row, fold_i, n){
    Z_test <- com
    row_col <- held_pairs[held_pairs[,3] == fold_i, c(1,2)]
    cat("fold_i: ", fold_i, '\n')
    Z_test[row_col] <- 0
    if(pairs_dim == "row") sum_func <- rowSums else sum_func <- colSums
    cat("\tfold i zero sum FALSE: ", any(sum_func(Z_test) == 0), '\n')
    if(all(sum_func(Z_test) != 0)){
        print("Nothing to do")
        return(NULL)
    }
    held_pairs_df <- as.data.frame(held_pairs)
    dim_zero <- which(sum_func(Z_test) == 0)
    print(dim_zero)
    for(dim_j in dim_zero){
        forbidden_fold <- c()
        cat("NEXT\n")
        #print(dplyr::filter(held_pairs_df, gr == i, col == col_j))
        pairs_i <- pairs_sums[pairs_sums[, pairs_dim] == dim_j, , drop = FALSE]
        print(class(pairs_i))
        print(pairs_i)
        pairs_rsums <- pairs_i[, 3]
        pairs_csums <- pairs_i[, 4]
        if(any(pairs_rsums < min_per_row)){
            print("Here row")
            if(all(pairs_rsums < min_per_row)){
                print("HERE ALL")
                pairs_i <- pairs_i
            }else{
                pairs_i <- pairs_i[-(which(pairs_rsums < min_per_row)), , drop = FALSE]
            }
            print(!pairs_rsums < min_per_row)
                print(pairs_i)

        }
        pairs_csums <- pairs_i[, 4]
        if(any(pairs_csums < min_per_col)){
            print("Here col")
            if(all(pairs_csums < min_per_col)){
                print("HERE ALL")
                pairs_i <- pairs_i
            }else{
                pairs_i <- pairs_i[-(which(pairs_csums < min_per_col)), , drop = FALSE]
            }
            print(!pairs_rsums < min_per_row)
                print(pairs_i)
        }
        cat("pairs_i\n")
        print(pairs_i)
        print(nrow(pairs_i))

        move_row <- sample(nrow(pairs_i), pmax(nrow(pairs_i) - 1, 1))
        cat("move row: ", move_row, '\n')
        to_fold <- sample(1:(n-1), length(move_row))
        cat("to_fold: ", to_fold, '\n')
        while(to_fold %in% forbidden_fold){
            to_fold <- sample(1:(n-1), 1)
            cat("to_fold: ", to_fold, '\n')
        }
        forbidden_fold <- c(forbidden_fold, to_fold)
        pairs_move <- c(pairs_i[move_row, c(1, 2)], gr = to_fold)
        cat("pairs_move: ", pairs_move, '\n')
        pairs_row <- which(held_pairs[, 1] == pairs_move[1] & held_pairs[, 2] == pairs_move[2])
        cat("held_pairs_df[pairs_row] before: ", as.matrix(held_pairs_df[pairs_row, ]), '\n')
        held_pairs_df[pairs_row, 3] <- pairs_move[3]
        cat("held_pairs_df[pairs_row] after: ", as.matrix(held_pairs_df[pairs_row, ]), '\n')

        #print(colSums(com)[col_j])
    }

    return(as.matrix(held_pairs_df))


}

#' Code adapted from HPprediction package
#' @export
create_cv_folds_old <-function(Z, n= 10, min_per_col = 2, min_per_row = 2,
                           get_zeros = FALSE){

    # TODO error number columns must match with bat data

    if(!(is.matrix(Z))) stop("Z argument must be of type matrix")
    ## n-fold cross validation
    if(max(range(Z))>1) Z[Z>0]<-1
    pairs = which(Z==1, arr.ind=T)
    colnames(pairs)<-c('row', 'col')

    colm = colSums(Z)
    min_rows <- which(rowSums(Z) < min_per_row)
    # Discount min_rows elements
    for(i in min_rows){
        pairs <- pairs[-which(pairs[, 'row'] == i), ]
        colm[which(Z[i, ] == 1)] <- colm[which(Z[i, ] == 1)] - 1
    }
    # Place zero to ignore columns
    for(i in 1:length(colm)){
        if(colm[i] < min_per_col){
            colm[i] <- 0
        }
    }
    size = floor(sum(colm)/n)
    gr = rep(size, n)
    if(sum(colm) %% size!=0) gr[n] =  gr[n] + sum(colm) %% size

    group_colm = rep(1:n,times = gr)[sample.int(sum(colm), sum(colm))]
    pair_list = numeric(sum(colm))
    for(i in 1:sum(colm)){
        a = which(colm>0)
        b = a[sample.int(length(a),1)]
        colm[b] = colm[b]-1
        pair_list[i]<-b
    }
    pair_list= tapply(pair_list, group_colm,identity)

    gr_list= list()
    bank= c()
    for(i in 1:n){
        a= table(pair_list[[i]])
        gr_rows = unlist(sapply(1:length(a), function(r){
            b = which(pairs[,'col']== as.numeric(names(a[r])))
            b =setdiff(b, bank)
            b[sample.int(length(b), a[r])]
        }))
        bank = c(bank, gr_rows)
        gr_list[[i]]<-cbind(gr_rows, i)
    }

    aux = do.call('rbind', gr_list)
    held_pairs = cbind(pairs[aux[,1], ], gr = aux[,2])

    if(get_zeros){
        pairs_0 = which(Z==0, arr.ind=T)
        held_zeros <- cbind(pairs_0[sample.int(nrow(pairs_0), sum(gr)), ], 
                            gr = rep(1:10, times = gr), value = 0)
        held_pairs <- cbind(held_pairs, value = 1)
        held_pairs <- rbind(held_pairs, held_zeros)
    }

    print(sprintf("Actual cross-validation rate is %0.3f" , table(held_pairs[,'gr'])/sum(1*(Z>0))))

    return(held_pairs[order(held_pairs[,'gr']),])

}

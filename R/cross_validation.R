#' Code adapted from HPprediction package
#' @export
create_cv_folds <-function(Z, n= 10, min_per_col = 2){

    # TODO error number columns must match with bat data

    # Remove empty columns and rows
    Z <- Z[rowSums(Z) != 0,
            colSums(Z) != 0]
    if(!(is.matrix(Z))) stop("Too many zero columns or rows in adj_mat.")
    ## n-fold cross validation
    if(max(range(Z))>1) Z[Z>0]<-1
    pairs = which(Z==1, arr.ind=T)
    colnames(pairs)<-c('row', 'col')

    colm = colSums(Z)
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

    print(sprintf("Actual cross-validation rate is %0.3f" , table(held_pairs[,'gr'])/sum(1*(Z>0))))

    return(held_pairs[order(held_pairs[,'gr']),])

}

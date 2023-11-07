#' @export
get_hsbm_results <- function(hsbm_output, input_names = TRUE){

    if(!inherits(hsbm_output, "hsbm.predict")){
        stop("hsbm_output argument must be an hsbm.predict object.")
    }

    com <- hsbm_output$data
    n_v1 <- nrow(hsbm_output$data)
    folds_res_list <- lapply(hsbm_output$predictions$probs,
                             function(x){
                                 df <- tidy_hsbm_results(x, n_v1)
                                 return(df)
                              })

    # Join all fold results
    folds_res_all <- Reduce(function(...)
                                dplyr::full_join(..., by = c("v1", "v2")),
                             folds_res_list)
    p_cols <- stringr::str_detect(colnames(folds_res_all), "p\\.")
    folds_res_all <- dplyr::mutate_if(folds_res_all, p_cols, ~tidyr::replace_na(.,0))
    folds_res_all$p <- rowMeans(folds_res_all[, p_cols])
    folds_res_all$sd <- apply(folds_res_all[, p_cols], 1, sd)
    folds_res_all$range <- apply(folds_res_all[, p_cols], 1, function(x) max(x) - min(x))
    edge_type_cols <- stringr::str_detect(colnames(folds_res_all),
                                 "edge_type\\.")

    et <- apply(folds_res_all[, edge_type_cols], 1,
                  function(x){
                          edge_types <- unique(x)
                          if("documented" %in% edge_types){
                              return("documented")
                          }else if("spurious_edge" %in% edge_types){
                              return("reconstructed")
                          }else{
                              return(paste0(stringr::str_replace_na(edge_types, replacement=""),
                                            collapse = ""))
                          }
                  })

    folds_res_all$edge_type <- et
    folds_select <- dplyr::select(folds_res_all, v1, v2,
                                  v1_names = v1_names.x,
                                  v2_names = v2_names.x,
                                  p, sd, range, edge_type)
    if(input_names){
        v1_row <- folds_select$v1 + 1
        v2_col <- folds_select$v2 - n_v1 + 1
        folds_select$v1_names <- rownames(com)[v1_row]
        folds_select$v2_names <- colnames(com)[v2_col]
    }

    hsbm_output$predictions$res_folds <- folds_res_list
    hsbm_output$predictions$res_averaged <- folds_select

    return(hsbm_output)

}


#' @export
top_links <- function(hsbm_output, n = 10){

    res_averaged <-  hsbm_output$predictions$res_averaged
    reconstructed <- dplyr::filter(res_averaged, edge_type == "reconstructed")

    cols <- c("v1_names", "v2_names", "p", "sd")
    reconstructed <- reconstructed[order(reconstructed$p, decreasing = TRUE), ][1:n, cols]
    row.names(reconstructed) <- NULL

    return(reconstructed)
}



#' @export
hsbm.reconstructed <- function(hsbm_out, pred_all = FALSE, rm_documented = FALSE, threshold = "prc_closest_topright"){

    hsbm_reconstructed <- list()
    hsbm_reconstructed$data <- hsbm_out$data
    hsbm_reconstructed$reconstructed_mats <- list()
    hsbm_reconstructed$reconstructed_stats <- list()
    n_folds <- length(hsbm_out$predictions$probs)
    for(i in 1:n_folds){
        reconstruction_i <- get_reconstruction(hsbm_out, fold_id = i, pred_all = pred_all,
                                               rm_documented, threshold = threshold)
        hsbm_reconstructed$reconstructed_mats[[i]] <- reconstruction_i$new_mat
        hsbm_reconstructed$reconstructed_stats[[i]] <- reconstruction_i$stats
    }

    tb_all <- do.call(rbind, hsbm_reconstructed$reconstructed_stats)
    tb_all <- data.frame(cbind(1:n_folds, tb_all))
    #yPRC is the baseline of Precision-Recall Curve
    colnames(tb_all) <- c("folds", "auc",
                          "aucpr", "yPRC", "thresh",    
                          "n_heldout", "pred_held_ones",
                          "n_ones", "pred_tot_ones", "total_pred_ones",
                          "precision", "sens", "spec", "ACC", "ERR","tss")

    hsbm_reconstructed$tb <- tb_all
    hsbm_reconstructed$new_mat <- avg_mat(hsbm_reconstructed$reconstructed_mats,
                                          thresh = mean(tb_all$thresh))

    hsbm_reconstructed$threshold <- threshold
    
    attr(hsbm_reconstructed, "class") <- "hsbm.reconstructed"
    
    return(hsbm_reconstructed)
}


get_reconstruction <- function(hsbm_out, fold_id, threshold, pred_all = FALSE, rm_documented = FALSE){
    df <- hsbm_out$predictions$probs[[fold_id]]
    com <- hsbm_out$data
    folds <- as.data.frame(hsbm_out$folds)
    gt_data <- tidy_hsbm_results(df, n_v1 = nrow(com))
    rows <- as.numeric(gt_data$v1) + 1
    cols <- as.numeric(gt_data$v2) - nrow(com) + 1
    ps <- gt_data$p
    com_train <- com
    row_col <- as.matrix(dplyr::filter(folds, gr == fold_id)[, c('row', 'col')])
    com_train[row_col] <- 0
    com_fit <- com_train
    com_fit[cbind(rows, cols)] <- ps
    if(!pred_all){
        com_tmp <- com
        com_tmp[] <- 0
        com_tmp[row_col] <- 1
        com_i <- com_tmp
    }else{
        com_i <- com
    }
    if(rm_documented){
        com_fit[com_train == 1] <- -1
        com_i[com_train == 1] <- -1
        yPRC <- sum(com_i[com_i != -1])/length(com_i[com_i != -1]) 
    }else{
	yPRC <- sum(com_i)/length(com_i) 
    }
    com_fit_long <- reshape2::melt(com_fit)$value
    com_i_long <- reshape2::melt(com_i)$value
    if(rm_documented){
        com_fit_long <- com_fit_long[com_fit_long != -1]
        com_i_long <- com_i_long[com_i_long != -1]
    }

    pred <- ROCR::prediction(predictions = com_fit_long, labels = com_i_long)
    aucpr <- as.numeric(ROCR::performance(pred, 'aucpr')@y.values)
    auc <- as.numeric(ROCR::performance(pred, 'auc')@y.values)
    f <- ROCR::performance(pred, "f")
    perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
    perf2 <- ROCR::performance(pred, "prec", "rec")

    thresh <-sel_thresh(threshold, perf, perf2, f)
	
    com_fit_bin <- ifelse(com_fit > thresh, 1, 0)

    n_heldout <- nrow(row_col)
    documented_ones <- sum(com)
    pred_tot_ones <- sum((com + com_fit_bin) == 2)/ sum(com)
    pred_held_ones <- sum((com[row_col] + com_fit_bin[row_col]) == 2)/sum(com[row_col])
    total_pred_ones <- sum(com_fit_bin)

    precision <- get_precision(com_i, com_train, com_fit_bin, pred_all = pred_all)

    return(list(new_mat = com_fit_bin,
                stats = c(auc, 
                          aucpr,
                          yPRC,
                          thresh,
                          n_heldout,
                          pred_held_ones,
                          documented_ones,
                          pred_tot_ones,
                          total_pred_ones,
                          precision$precision,
                          precision$sens,
                          precision$spec,
                          precision$ACC,
                          precision$ERR,
                          precision$tss)))
}

get_precision <- function(com, com_train, com_fit_bin, pred_all = FALSE){

    com= 1*(com>0)
    if(pred_all){
        Z1 = com==1
    }else{
        Z1 = (com - com_train)==1
    }
    m = sum(1*(Z1))
    Z2 = com==0
    n = sum(1*Z2)
    TP = sum(com_fit_bin[Z1])
    FN = m - TP
    FP = sum(com_fit_bin[Z2])
    TN = n - FP

    precision <- TP/(TP + FP)
    sens <- TP/(TP + FN)
    spec <- TN/(TN + FP) 
    ACC <- (TP + TN) / (TP + TN + FN + FP)
    ERR <- (FP + FN) / (TP + TN + FN + FP)
    tss <- (TP/(TP+FN))+(TN/(TN+FP))-1

    return(list(precision = precision, sens = sens, spec= spec, ACC = ACC, ERR = ERR, tss = tss))
}


avg_mat <- function(reconstructed_mats_list, thresh){

    mat_avg <- Reduce("+", reconstructed_mats_list)/length(reconstructed_mats_list)
    mat_avg_bin <- 1*(mat_avg > thresh)
    return(mat_avg_bin)
}

# First column is v1 and second is v2
tidy_hsbm_results <- function(gt_df, n_v1 = 447){
    last_v1_v <- n_v1 - 1
    # Remove v2 to v2 and v1 to v1 links
    # (this can be aranged by flag the graph as bipartite)
    v2_v2_link <- apply(gt_df[, 1:2], 1,
                          function(x) x[1] > last_v1_v & x[2] > last_v1_v)
    gt_df <- gt_df[!v2_v2_link, ]
    v1_v1_link <- apply(gt_df[, 1:2], 1,
                           function(x) x[1] <= last_v1_v & x[2] <= last_v1_v)
    gt_df <- gt_df[!v1_v1_link, ]

    #Swap links where v1 is v2 and v2 v1
    gt_df[, 1:2] <- t(apply(gt_df[, 1:2], 1,
                      function(x){
                          if(x[1] > last_v1_v & x[2] <= last_v1_v){
                              return(c(x[2],x[1]))
                          }else{
                              return(x)
                          }
                         }
                        )
                       )

    return(gt_df)
}



sel_thresh<- function(threshold, perf, perf2, f) {
	 THRESH_names<-c("roc_youden",
			  "roc_closest_topleft", 
			  "roc_equal_sens_spec", 
			  "roc_no_omission", 
			  "prc_min_rec_prec", 
			  "prc_equal_rec_prec", 
			  "prc_closest_topright",
			  "prc_max_F1")
  if (threshold %in% THRESH_names) {
	  if (threshold == "roc_youden") {
	    df <- data.frame(cut = perf@alpha.values[[1]], fpr = perf@x.values[[1]], tpr = perf@y.values[[1]])
	    roc_youden <- df[which.max(df$tpr - df$fpr), "cut"] 
	    thresh<-roc_youden
	  } else if (threshold == "roc_closest_topleft") {
	    df <- data.frame(cut = perf@alpha.values[[1]], fpr = perf@x.values[[1]], tpr = perf@y.values[[1]])
	    roc_closest_topleft <- df[which.min((1-df$tpr)^2 + (df$fpr)^2), "cut"]
	    thresh<-roc_closest_topleft 
	  } else if (threshold == "roc_equal_sens_spec") {
	    df <- data.frame(cut = perf@alpha.values[[1]], fpr = perf@x.values[[1]], tpr = perf@y.values[[1]])
 	    roc_equal_sens_spec <- df[which.min(abs(df$tpr - (1-df$fpr))), "cut"]
	    thresh<-roc_equal_sens_spec
	  } else if (threshold == "roc_no_omission") {
	    df <- data.frame(cut = perf@alpha.values[[1]], fpr = perf@x.values[[1]], tpr = perf@y.values[[1]])
 	    roc_no_omission <- max(pred@cutoffs[[1]][pred@fn[[1]] == 0])
	    thresh<-roc_no_omission
	  } else if (threshold == "prc_min_rec_prec") {
	    df2 <- data.frame(cut = perf2@alpha.values[[1]], Recall = perf2@x.values[[1]], Precision = perf2@y.values[[1]])
 	    prc_min_rec_prec <- df2[which.min(df2$Recall + df2$Precision), "cut"]
	    thresh<-prc_min_rec_prec
	  } else if (threshold == "prc_equal_rec_prec") {
	    df2 <- data.frame(cut = perf2@alpha.values[[1]], Recall = perf2@x.values[[1]], Precision = perf2@y.values[[1]])
 	    prc_equal_rec_prec <- df2[which.min(abs(df2$Recall - df2$Precision)), "cut"]
	    thresh<-prc_equal_rec_prec
	  } else if (threshold == "prc_closest_topright") {
	    df2 <- data.frame(cut = perf2@alpha.values[[1]], Recall = perf2@x.values[[1]], Precision = perf2@y.values[[1]])
	    prc_closest_topright <- df2[which.min((1-df2$Recall)^2 + (1-df2$Precision)^2), "cut"]
	    thresh<-prc_closest_topright
	  } else if (threshold == "prc_max_F1") {
	    df3 <- data.frame(Cutoff = f@x.values[[1]], PrecisionRecallFmeasure = f@y.values[[1]])
	    prc_max_F1 <- df3[which.max(df3$PrecisionRecallFmeasure), "Cutoff"]
	    thresh<-prc_max_F1
	  } 
  } else if (is.numeric(threshold)) {
        if (threshold >= 0 && threshold <= 1) {
          thresh <- threshold
        } else {
          stop("Threshold must be a numeric value between 0 and 1")
        }
  } else {
    stop("Invalid threshold option")
  }	
return(thresh)
}

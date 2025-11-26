#' @name hsbm.reconstructed
#'
#' @title Generate a reconstructed binary matrix from
#' Hierarchical Stochastic Block Model (HSBM)
#'
#' @description This function generates a reconstructed binary
#' matrix based on the HSBM model.
#'
#' @param hsbm_out An object of class \code{hsbm.output} containing the output
#' from the HSBM analysis, including predictions and input data.
#' @param rm_documented (\emph{optional, default} \code{TRUE}) \cr
#' A \code{logical} value indicating whether to remove documented entries (1s) from the evaluation.
#' @param na_treatment (\emph{optional, default} \code{"na_to_0"}) \cr
#' A \code{character} string specifying how to handle \code{NA} values derived from HSBM predictions.
#' Options include \code{"na_to_0"}, \code{"ignore_na"}, and \code{"keep_na"}. See Details for available options.
#' @param threshold (\emph{optional, default} \code{"prc_closest_topright"}) \cr
#' A \code{character} string or \code{numeric} value specifying the method to determine the threshold for
#' binary classification of predictions. See Details for available options.
#' @param consistency_matrix (\emph{optional, default} \code{"average_thresholded"}) \cr
#' A \code{character} string specifying the method for creating the new reconstructed matrix.
#' Options include \code{"average_thresholded"} and \code{"ensemble_binary"}. See Details for available options.
#' @param ensemble_threshold (\emph{optional, default} \code{NULL}) \cr
#' A \code{numeric} value \code{(0-1)} for \code{consistency_matrix = "ensemble_binary"} specifying the minimum proportion
#' of folds that must predict a link as present \code{(1)} for it to be included in the final binary matrix. Default
#' is \code{0.1} (for 10 folds it means at least one fold must predict presence).
#' The behavior of \code{ensemble_threshold} depends on the value of \code{consistency_matrix}:
#' - If \code{"average_thresholded"}: \code{ensemble_threshold} is applied directly to the averaged probabilities across all folds to binarize the final matrix. If \code{NULL}, the mean of the thresholds computed for each fold is used as the default.
#' - If \code{"ensemble_binary"}: \code{ensemble_threshold} specifies the proportion of folds in which a link must be predicted as 1 to be classified as 1 in the final matrix. If \code{NULL}, the default value is \code{0.1} (i.e., the link must be predicted as 1 in at least 10\% of the folds).
#'
#' @return
#' An object of class \code{hsbm.reconstructed} containing:
#' - \code{$data} The original binary input matrix.
#' - \code{$pred_mats} A \code{list} of matrices with predicted probabilities for each fold.
#' - \code{$reconstructed_df}: A \code{list} with:
#'   - \code{res_folds}: A \code{list} of \code{data.frame} objects summarizing predictions for each fold. It is a cleaned version of the \code{$predictions$probs} component in \code{hsbm.predict}.
#'     - \code{v1}, \code{v2}: Indices of the first (in rows) and second (in columns) of nodes.
#'     - \code{v1_names}, \code{v2_names}: Node names for the first (in rows) and second (in columns) of nodes.
#'     - \code{p}: Predicted probabilities for each edge/link.
#'     - \code{edge_type}: Type of edge (e.g., "documented" or "reconstructed").
#'   - \code{res_averaged}: A \code{data.frame} summarizing predictions averaged across all folds, with columns:
#'     - \code{v1}, \code{v2}: Indices of the first (in rows) and second (in columns) of nodes.
#'     - \code{v1_names}, \code{v2_names}: Node names for the first (in rows) and second (in columns) of nodes.
#'     - \code{p}: Average predicted probabilities for each edge/link.
#'     - \code{sd}: Standard deviation of predicted probabilities across folds.
#'     - \code{range}: Range of predicted probabilities across folds.
#'     - \code{edge_type}: Type of edge (e.g., "documented", "reconstructed").
#'     - \code{nr_na}: Number of folds where a given edge/link had \code{NA} predictions.
#' - \code{$stats} A \code{data.frame} summarizing evaluation metrics for each fold. It includes the following columns:
#'   - \code{folds}: The index of the cross-validation fold.
#'   - \code{auc}: The Area Under the Curve (AUC) of the Receiver Operating Characteristic (ROC) curve for the fold.
#'   - \code{aucpr}: The Area Under the Curve (AUC) of the Precision-Recall Curve (PRC) for the fold.
#'   - \code{yPRC}: The baseline of Precision-Recall Curve for the fold.
#'   - \code{thresh}: The binary classification threshold value applied to the predicted probabilities to convert them into binary classifications (0 or 1) for the fold.
#'   - \code{n_heldout}: The number of edges/links in the held-out set of the cross-validation fold, i.e., the edges that were excluded from the model training and used for evaluation.
#'   - \code{pred_held_ones}: The proportion of held-out edges/links that were correctly predicted as 1s by the model.
#'   - \code{n_ones}: The total number of positive edges/links in the original data.
#'   - \code{pred_tot_ones}: The proportion of positive edges/links in the original data that were correctly predicted as positive by the model.
#'   - \code{total_pred_ones}: The total number of positive edges/links predicted by the model in the reconstructed matrix.
#'   - \code{precision}: The precision of the model, calculated as the ratio of true positives to the sum of true positives and false positives.
#'   - \code{sens}: The sensitivity (or recall) of the model.
#'   - \code{spec}: The specificity of the model.
#'   - \code{ACC}: The overall accuracy of the model.
#'   - \code{ERR}: The error rate of the model.
#'   - \code{tss}: The True Skill Statistic (TSS).
#' - \code{$new_mat}: The final reconstructed binary matrix, combining predictions across folds using the specified \code{consistency_matrix}.
#' - \code{$threshold}: The method or value used to determine the binary classification threshold.
#'
#' @details
#' - The \code{rm_documented} parameter determines whether observed/documented edges/links (1s) are considered in the evaluation and
#' reconstruction process. When using \code{method = "conditional_missing"}, it is recommended to have \code{rm_documented = TRUE}
#' to exclude observed/documented links (1s) from evaluation and network reconstruction, as no probabilities are computed for them.
#' Conversely, with \code{"method = "marginal_all""}, we recommend setting
#' \code{rm_documented = FALSE} to include observed/documented edges/links in the evaluation and reconstruction process.
#'
#' - The \code{na_treatment} parameter specifies how to handle \code{NA} values in the predictions. Available options are:
#'   - \code{"na_to_0"}: Interprets \code{NA} values as no evidence for an existing link and assigns them a value of zero.
#' This assumes that \code{NA} indicates the absence of evidence for a link.
#'   - \code{"ignore_na"}: Ignores \code{NA} values when calculating averages and other metrics, treating them as missing data.
#'   - \code{"keep_na"}: Retains \code{NA} values in the output when at least one fold predicts an \code{NA}, preserving uncertainty.
#'
#' - The \code{threshold} parameter defines the method for determining the binary classification threshold. Available options include:
#'   - \code{"roc_youden"}: Maximizes the Youden's J statistic (sensitivity + specificity - 1).
#'   - \code{"roc_closest_topleft"}: Minimizes the distance to the top-left corner in the ROC space.
#'   - \code{"roc_equal_sens_spec"}: Equalizes sensitivity and specificity.
#'   - \code{"roc_no_omission"}: Maximizes the threshold with no false negatives.
#'   - \code{"prc_min_rec_prec"}: Minimizes the sum of recall and precision.
#'   - \code{"prc_equal_rec_prec"}: Equalizes recall and precision.
#'   - \code{"prc_closest_topright"}: Minimizes the distance to the top-right corner in the PRC space.
#'   - \code{"prc_max_F1"}: Maximizes the F1 score.
#'   - Alternatively, a numeric value between 0 and 1 can be provided as a custom threshold.
#'
#' - The \code{consistency_matrix} parameter specifies the method used to generate the new reconstructed matrix. Valid options are:
#'   - \code{"average_thresholded"}: This method averages the predicted probability matrices across all folds and applies a threshold to
#' transform the averaged probabilities into a final binary matrix. The threshold used is the average of the thresholds calculated for each fold.
#'   - \code{"ensemble_binary"}: This method first transform predicted probability matrices into binary matrices using fold-specific thresholds.
#' Then, an \code{ensemble_threshold} is applied to aggregate these binary matrices into a final binary matrix, setting an entry to 1 if it exceeds
#' the `ensemble_threshold`. The default `ensemble_threshold` is `0.1` if not specified.
#'
#'
#' @seealso \code{\link{hsbm.predict}}
#'
#' @examples
#' ## Not run:
#' # Example workflow to generate `myPred`:
#' data(dat, package = "sabinaHSBM")
#'
#' # Prepare input for HSBM
#' myInput <- hsbm.input(data = dat,
#'                       n_folds = 10)
#'
#' # Run HSBM predictions
#' myPred <- hsbm.predict(hsbm_input = myInput,
#'                       method = "conditional_missing",
#'                       iter = 1000,
#'                       wait = 1000)
#' ## End(Not run)
#'
#' # Load example HSBM reconstructed results
#' data(myPred, package = "sabinaHSBM")
#'
#' myReconst <- hsbm.reconstructed(myPred)
#'
#' # Get the final reconstructed binary matrix
#' reconstructed_matrix <- myReconst$new_mat
#'
#' # Evaluation metrics
#' eval_metrics <- myReconst$stats
#' print(eval_metrics)
#'
#' # Final averaged matrix
#' averaged_matrix <- myReconst$reconstructed_df$res_averaged
#'
#' # Plot Reconstructed matrix
#' plot_interaction_matrix(myReconst$new_mat, order_mat = FALSE)
#'
#' @export
hsbm.reconstructed <- function(hsbm_out, rm_documented = TRUE,
                               na_treatment = "na_to_0",
                               threshold = "roc_youden",
                               consistency_matrix = "average_thresholded",
                               ensemble_threshold = NULL){

    if(!inherits(hsbm_out, "hsbm.predict")) {
        stop("Error: hsbm_out must be an object of hsbm.predict class. Consider running hsbm.predict() function.")
    }

    if(!(consistency_matrix %in% c("average_thresholded", "ensemble_binary"))){
        stop("\nError: Invalid value for `consistency_matrix`.",
             " It must be 'ensemble_binary' or 'average_thresholded'.\n")
    }
    if(na_treatment != "na_to_0" & hsbm_out$method == "conditional_missing"){
        stop("\nError: `na_treatment` cannot be set for `conditional_missing` method since there are no NAs.\n")
    }
    if(rm_documented & hsbm_out$method == "marginal_all") {
        warning("\`rm_documented` was set to TRUE for the `marginal_all` method. ",
                "Documented links (original 1s) were removed for evaluation. ",
                "Ensure this setting is appropriate for your analysis.\n")
    }
    if(!rm_documented & hsbm_out$method == "conditional_missing") {
        warning("\n`rm_documented` was set to FALSE for the `conditional_missing` method. ",
                "Documented links (1s) were included in evaluation, even though no probabilities ",
                "were computed for them. Ensure this setting is appropriate for your analysis.\n")
    }
    if (consistency_matrix == "ensemble_binary" & is.null(ensemble_threshold)) {
        message("\nNo `ensemble_threshold` provided; defaulting to 0.1.\n")
    }

    hsbm_reconstructed <- list()
    hsbm_reconstructed$data <- hsbm_out$data
    hsbm_reconstructed$pred_mats <- list()
    hsbm_reconstructed$reconstructed_df <- list()
    binary_mats <- list()
    reconstructed_stats <- list()

    is_bipartite <- hsbm_out$is_bipartite

    hsbm_reconstructed$reconstructed_df <- get_hsbm_results(hsbm_out,
                                                            input_names = TRUE,
                                                            na_treatment = na_treatment,
                                                            is_bipartite = is_bipartite)
    n_folds <- length(hsbm_out$probs)
    for(i in 1:n_folds){
        reconstruction_i <- get_reconstruction(hsbm_reconstructed$reconstructed_df$res_folds,
                                               fold_id = i,
                                               com = hsbm_out$data,
                                               folds = hsbm_out$folds,
                                               method = hsbm_out$method,
                                               rm_documented = rm_documented,
                                               na_treatment = na_treatment,
                                               threshold = threshold,
                                               is_bipartite = is_bipartite)
        hsbm_reconstructed$pred_mats[[i]] <- reconstruction_i$pred_mat
        reconstructed_stats[[i]] <- reconstruction_i$stats
        binary_mats[[i]] <- reconstruction_i$new_bin_mat
    }

    tb_all <- do.call(rbind, reconstructed_stats)
    tb_all <- data.frame(cbind(1:n_folds, tb_all))
    #yPRC is the baseline of Precision-Recall Curve
    colnames(tb_all) <- c("folds", "auc",
                          "aucpr", "yPRC", "thresh",
                          "n_heldout", "pred_held_ones",
                          "n_ones", "pred_tot_ones", "total_pred_ones",
                          "precision", "sens", "spec", "ACC", "ERR","tss")

    hsbm_reconstructed$stats <- tb_all
    if(consistency_matrix == "ensemble_binary"){
        ens_thresh <- if(is.null(ensemble_threshold)) 0.1 else ensemble_threshold
        hsbm_reconstructed$new_mat <- avg_mat(binary_mats,
                                              thresh = ens_thresh,
                                              na_treatment = na_treatment)
    }else{
        thresh <- mean(tb_all$thresh)
        hsbm_reconstructed$new_mat <- avg_mat(hsbm_reconstructed$pred_mats,
                                              thresh = thresh,
                                              na_treatment = na_treatment)
    }

    hsbm_reconstructed$threshold <- threshold

    attr(hsbm_reconstructed, "class") <- "hsbm.reconstructed"

    return(hsbm_reconstructed)
}

get_hsbm_results <- function(hsbm_output, input_names = TRUE,
                             na_treatment = "na_to_0",
                             is_bipartite = TRUE){

    if(!inherits(hsbm_output, "hsbm.predict")){
        stop("Error: hsbm_output argument must be an hsbm.predict object.")
    }

    com <- hsbm_output$data
    n_v1 <- nrow(hsbm_output$data)
    folds_res_list <- lapply(hsbm_output$probs,
                            function(x){
                                df <- tidy_hsbm_results(gt_df = x,
                                                        n_v1 = n_v1,
                                                        is_bipartite = is_bipartite)
                                return(df)
                              })

    # Join all fold results
    folds_res_all <- Reduce(function(...)
                                dplyr::full_join(..., by = c("v1", "v2")),
                             folds_res_list)
    p_cols <- stringr::str_detect(colnames(folds_res_all), "^p\\.?")
    nr_na <- apply(folds_res_all[, p_cols], 1, function(x) sum(is.na(x)))
    if(na_treatment == "na_to_0"){
        na_rm <- TRUE
        folds_res_all <- dplyr::mutate_if(folds_res_all, p_cols,
                                          ~tidyr::replace_na(.,0))
    }else if(na_treatment == "ignore_na"){
        na_rm <- TRUE
    }else if(na_treatment == "keep_na"){
        na_rm <- FALSE
    }else{
        stop("Error: Unkown na_treatment argument.")
    }
    all_p <- rowMeans(folds_res_all[, p_cols], na.rm = na_rm)
    all_sd <- apply(folds_res_all[, p_cols], 1, function(x) stats::sd(x, na.rm = na_rm))
    all_range <- apply(folds_res_all[, p_cols], 1,
                       function(x) max(x, na.rm = na_rm) - min(x, na.rm = na_rm))
    folds_res_all$p <- all_p
    folds_res_all$sd <- all_sd
    folds_res_all$range <- all_range
    folds_res_all$nr_na <- nr_na
    edge_type_cols <- stringr::str_detect(colnames(folds_res_all), "edge_type\\.")

    et <- apply(folds_res_all[, edge_type_cols], 1,
                  function(x){
                          edge_types <- unique(x)
                          if("documented" %in% edge_types){
                              return("documented")
                          }else if("spurious_edge" %in% edge_types){
                              return("spurious_edge")
                         }else{
                              return(paste0(stringr::str_replace_na(edge_types,
                                                                    replacement=""),
                                            collapse = ""))
                          }
                  })

    folds_res_all$edge_type <- et
    folds_select <- dplyr::select(folds_res_all, "v1", "v2",
                                  v1_names = "v1_names.x",
                                  v2_names = "v2_names.x",
                                  "p", "sd", "range", "edge_type", "nr_na")
    if(input_names){
        v1_row <- folds_select$v1 + 1
        folds_select$v1_names <- rownames(com)[v1_row]
        if(is_bipartite){
            v2_col <- folds_select$v2 - n_v1 + 1
            folds_select$v2_names <- colnames(com)[v2_col]
        }else{
            v2_col <- folds_select$v2 + 1
            folds_select$v2_names <- colnames(com)[v2_col]
        }
    }

    min_size_warning <- 0.25 * prod(dim(com))
    if((hsbm_output$method == "marginal_all") && (nrow(folds_select) < min_size_warning)) {
        percentage <- 100 * (nrow(folds_select) / prod(dim(com)))
        warning(sprintf("Predictions obtained for %.2f%% of the links. Consider increasing the number of iterations (iter argument in hsbm.predict()).", percentage))
    }

    return(list(res_folds = folds_res_list, res_averaged = folds_select))

}

get_reconstruction <- function(res_folds, fold_id, com, folds, method, threshold,
                               na_treatment, rm_documented = FALSE,
                               spurious_edges = FALSE, is_bipartite = TRUE){
    df <- res_folds[[fold_id]]
    rows <- as.numeric(df$v1) + 1
    if(is_bipartite){
        cols <- as.numeric(df$v2) - nrow(com) + 1
    }else{
        cols <- as.numeric(df$v2) + 1
    }
    ps <- df$p
    com_train <- com
    com_i <- com
    cond <- folds[, 3] == fold_id
    row_col <- as.matrix(dplyr::filter(as.data.frame(folds), cond)[, c('row', 'col')])
    com_train[row_col] <- 0
    com_fit <- com_train
    com_fit[cbind(rows, cols)] <- ps

    if(rm_documented){
        com_fit[com_train == 1] <- NA
        com_i[com_train == 1] <- NA
    }
    com_fit_vec <- reshape2::melt(com_fit, na.rm = TRUE)$value
    com_i_vec <- reshape2::melt(com_i, na.rm = TRUE)$value

    yPRC <- sum(com_i_vec)/length(com_i_vec)
    pred <- ROCR::prediction(predictions = com_fit_vec, labels = com_i_vec)
    aucpr <- as.numeric(ROCR::performance(pred, 'aucpr')@y.values)
    auc <- as.numeric(ROCR::performance(pred, 'auc')@y.values)
    f <- ROCR::performance(pred, "f")
    perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
    perf2 <- ROCR::performance(pred, "prec", "rec")

    thresh <-sel_thresh(threshold, perf, perf2, f, pred)

    perf_prec <- ROCR::performance(pred, measure = "prec")
    perf_sens <- ROCR::performance(pred, measure = "sens")
    perf_spec <- ROCR::performance(pred, measure = "spec")
    perf_acc <- ROCR::performance(pred, measure = "acc")

    cutoffs <- perf_prec@x.values[[1]]
    precisions <- perf_prec@y.values[[1]]
    sensitivities <- perf_sens@y.values[[1]]
    specificities <- perf_spec@y.values[[1]]
    accuracies <- perf_acc@y.values[[1]]

    closest_index <- which.min(abs(cutoffs - thresh))

    precision <- precisions[closest_index]
    sens <- sensitivities[closest_index]
    spec <- specificities[closest_index]
    ACC <- accuracies[closest_index]
    ERR <- 1 - ACC
    tss <- sens + spec - 1

    # Get documented values back
    if(rm_documented & method == "conditional_missing"){
        com_fit[com_train == 1] <- 1
    }else if(rm_documented & method == "marginal_all"){
        com_fit[cbind(rows, cols)] <- ps
        com_fit[com_train == 1] <- 1
    }

    if(na_treatment != "na_to_0"){
        com_fit[] <- NA
        com_fit[cbind(rows, cols)] <- ps
    }

    # Make sure when p == 0.5 does not give missing links
    if(thresh == 0.5){
        com_fit_bin <- ifelse(com_fit > thresh, 1, 0)
    }else{
        com_fit_bin <- ifelse(com_fit >= thresh, 1, 0)
    }
    n_heldout <- nrow(row_col)
    documented_ones <- sum(com)
    pred_tot_ones <- sum((com + com_fit_bin) == 2, na.rm = TRUE)/ sum(com)
    pred_held_ones <- sum((com[row_col] + com_fit_bin[row_col]) == 2, na.rm = TRUE)/sum(com[row_col])
    total_pred_ones <- sum(com_fit_bin, na.rm = TRUE)

    spurious_ones <- NULL
    if(spurious_edges){
        cond <- (df$edge_type == "spurious_edge")
        spurious <- dplyr::filter(df, cond)
        spurious_ones <- sum(spurious$p > thresh)
    }

    return(list(pred_mat = com_fit,
                new_bin_mat = com_fit_bin,
                stats = c(auc,
                          aucpr,
                          yPRC,
                          thresh,
                          n_heldout,
                          pred_held_ones,
                          spurious_ones,
                          documented_ones,
                          pred_tot_ones,
                          total_pred_ones,
                          precision,
                          sens,
                          spec,
                          ACC,
                          ERR,
                          tss)))

}


avg_mat <- function(reconstructed_mats_list, thresh, na_treatment, method = "average_thresholded"){

    if(na_treatment == "keep_na"){
        na_rm <- FALSE
    }else{
        na_rm <- TRUE
    }
    mat_avg <- apply(simplify2array(reconstructed_mats_list), 1:2, mean,
                     na.rm = na_rm)

    mat_avg_bin <- 1*(mat_avg >= thresh)

    return(mat_avg_bin)

}


#' @importFrom rlang .data
tidy_hsbm_results <- function(gt_df, n_v1 = 447, is_bipartite = TRUE){

    if(!is_bipartite){
        gt_df <- dplyr::mutate(gt_df,
                              new_v1 = pmax(.data$v1, .data$v2),
                              new_v2 = pmin(.data$v1, .data$v2))
        gt_df <- dplyr::group_by(gt_df,
                                .data$new_v1,
                                .data$new_v2)
        gt_df <- dplyr::summarise(gt_df,
                                  p = mean(.data$p),
                                  v1_names = dplyr::first(.data$v1_names),
                                  v2_names = dplyr::first(.data$v2_names),
                                  edge_type = dplyr::first(.data$edge_type),
                                  .groups = 'drop')
        gt_df <- dplyr::mutate_at(gt_df, c("new_v1", "new_v2"), as.numeric)
        colnames(gt_df)[1:2] <- c("v1", "v2")
        return(gt_df)
    }

    last_v1_v <- n_v1 - 1
    # Remove v2 to v2 and v1 to v1 links
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

sel_thresh<- function(threshold, perf, perf2, f, pred) {
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
          stop("Error: Threshold must be a numeric value between 0 and 1")
        }
  } else {
    stop("Error: Invalid threshold option")
  }
  return(thresh)
}

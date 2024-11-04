#' @name hsbm.predict
#'
#' @title Generate edge/link predictions using Hierarchical Stochastic Block Model (HSBM).
#'
#' @description This function performs predictions based on the HSBM for a specified edge list (or all edge lists) from an `hsbm.input` object.
#'
#' @param hsbm_input An object of class \code{hsbm.input} containing the necessary data and configurations for running the HSBM analysis.
#' A \code{logical} value indicating whether to add the number of observations per row and column to the edge list.
#' @param elist_i (\emph{optional, default} \code{NULL}) \cr
#' A \code{numeric} value specifying the index of the edge list (fold) from an `hsbm.input`to run predictions on. If \code{NULL}, predictions are run for all edge lists.
#' @param method (\emph{optional, default} \code{"binary_classifier"}) \cr
#' A \code{character} string specifying the method used for the HSBM prediction. Options include \code{"binary_classifier"} and \code{"full_reconstruction"}.
#' @param verbose (\emph{optional, default} \code{TRUE}) \cr
#' A \code{logical} value indicating whether to print progress messages during prediction.
#'
#' @return
#' An object of class \code{hsbm.predict} containing the edge/link predictions and group assignments for the specified edge lists (fold):
#' - \code{$data} The binary bipartite \code{matrix} of input data.
#' - \code{$folds} A \code{matrix} of cross-validation fold assignments for each held-out edge/link.
#' - \code{$edgelist} A \code{list} of edge lists generated for each fold.
#' - \code{$predictions$probs} A \code{list} where each element is a \code{data.frame} with the predicted probabilities (p) for the edges/links in the corresponding edge list, according to the HSBM model.
#' - \code{$predictions$groups} A \code{list} where each element is a \code{data.frame} that contains the group assignments for each node in the network for the corresponding edge list.
#'
#' @details
#' - The \code{hsbm_input} parameter should be an object of class \code{hsbm.input}, which includes the input data, the cross-validation folds, and corresponding edge lists.
#' - The \code{elist_i} parameter allows you to specify a particular edge list to run predictions on. If not specified, predictions are run on all edge lists.
#' - The \code{method} parameter defines the prediction method to be used. Currently, both \code{"binary_classifier"} (see \url{https:// #@@@JMB info binary}) and \code{"full_reconstruction"} are supported.
#' - The \code{verbose} parameter, when set to \code{TRUE}, enables the display of progress messages, which is useful for tracking the computation process.
#'
#' @seealso \code{\link{hsbm.input}}
#'
#' @examples
#' # Load example data and prepare input
#' data(dat, package = "sabinaHSBM")
#' myInput <- hsbm.input(data = dat, 
#'                       n_folds = 10, 
#'                       method = "binary_classifier", 
#'                       iter = 1000)
#'
#' # Generate predictions for the prepared input
#' myPred <- hsbm.predict(hsbm_input = myInput)
#'
#' @export
hsbm.predict <- function(hsbm_input, elist_i = NULL, method = "binary_classifier",
                         verbose = TRUE){
    hsbm_name <- as.list(match.call())$hsbm_input
    if(!inherits(hsbm_input, "hsbm.input")){
        stop("hsbm must be an object of hsbm.input class. Consider running hsbm_input() function.")
    }

    reticulate::py_run_string(import_modules())
    reticulate::py_run_string(add_taxa_vertex_prop())
    reticulate::py_run_string(create_graph())
    reticulate::py_run_string(get_missing_edges())
    reticulate::py_run_string(hsbm_link_prediction())
    reticulate::py_run_string(get_groups())
    reticulate::py_run_string(get_predicted_network())

    hsbm_output <- hsbm_input

    predictions <- list()
    predictions$probs <- list()
    predictions$groups <- list()
    if(is.null(elist_i)){
        elist_predict <- 1:length(hsbm_input$edgelist)
    }else{
        elist_predict <- elist_i
    }
    for(i in elist_predict){
        if(verbose){
            cat("Computing predictions for fold ", i, "\n")
        }
        i_py <- i - 1
        reticulate::py_run_string(stringr::str_glue("elists = r.{hsbm_name}['edgelist']"))
        reticulate::py_run_string(paste0("elist = elists[", i_py, "]"))
        reticulate::py_run_string("g = create_graph(elist)")
        reticulate::py_run_string(stringr::str_glue("res_dict = hsbm_predict(g, elist, ",
                                                    "method = r.{hsbm_name}['method'], ",
                                                    "force_niter = r.{hsbm_name}['iter'], ",
                                                    "wait = r.{hsbm_name}['wait'])"))
        reticulate::py_run_string("groups_df = get_groups(res_dict['state_min_dl'], res_dict['graph'])")

        hsbm_output$predictions$probs[[i]] <- py$res_dict$pred_probs
        hsbm_output$predictions$groups[[i]] <- py$groups_df

        #reticulate::py_save_object(py$res_dict, str_glue("res_dict{i}.pkl"))

        reticulate::py_run_string("del res_dict; del groups_df; del g; gc.collect()")

    }

    attr(hsbm_output, "class") <- "hsbm.predict"

    return(hsbm_output)
}


hsbm.predict2 <- function(hsbm_input, add_n_x = TRUE, verbose = TRUE){
    hsbm_name <- as.list(match.call())$hsbm_input
    if(!inherits(hsbm_input, "hsbm.input")){
        stop("hsbm must be an object of hsbm.input class. Consider running hsbm_input() function.")
    }

    cat("\tpy object begin: ", format(object.size(py), units = "MB"), "\n")

    reticulate::py_run_string(import_modules())
    reticulate::py_run_string(add_taxa_vertex_prop())
    reticulate::py_run_string(create_graph())
    reticulate::py_run_string(get_missing_edges())
    reticulate::py_run_string(hsbm_link_prediction())
    reticulate::py_run_string(get_groups())

    hsbm_output <- hsbm_input
    cat("\thsbm_output begin: ", format(
                                     object.size(hsbm_output),
                                     units = "MB"),
        "\n")
    rm(hsbm_input); gc()

    hsbm_output$predictions <- list()
    hsbm_output$predictions$probs <- list()
    hsbm_output$predictions$groups <- list()
    for(i in 1:length(hsbm_output$edgelist)){
        if(verbose){
            cat("Computing predictions for fold ", i, "\n")
        }
        i_py <- i - 1
        reticulate::py_run_string(stringr::str_glue("elists = r.{hsbm_name}['edgelist']"))
        reticulate::py_run_string(paste0("elist = elists[", i_py, "]"))
        reticulate::py_run_string("g = create_graph(elist)")
        reticulate::py_run_string(stringr::str_glue("res_dict = hsbm_predict(g, elist, ",
                                                    "force_niter = r.{hsbm_name}['iter'])"))
        reticulate::py_run_string("groups_df = get_groups(res_dict['state_min_dl'], res_dict['graph'])")

        hsbm_output$predictions$probs[[i]] <- py$res_dict$pred_probs
        hsbm_output$predictions$groups[[i]] <- py$groups_df

        reticulate::py_run_string("del res_dict; del groups_df; del g; gc.collect()")

        cat("\tpy object: ", format(object.size(py), units = "MB"), "\n")
        cat("\thsbm_output: ", format(object.size(hsbm_output), units = "MB"), "\n")
    }

    attr(hsbm_output, "class") <- "hsbm.predict"

    return(hsbm_output)
}

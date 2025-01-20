#' @name hsbm.predict
#'
#' @title Generate edge/link predictions using Hierarchical Stochastic Block Model (HSBM).
#'
#' @description This function performs predictions based on the HSBM for a specified edge list (or all edge lists) from an `hsbm.input` object.
#'
#' @param hsbm_input An object of class \code{hsbm.input} containing the necessary data and configurations for running the HSBM analysis.
#' @param elist_i (\emph{optional, default} \code{NULL}) \cr
#' A \code{numeric} value specifying the index of the edge list (fold) from an `hsbm.input`to run predictions on. If \code{NULL}, predictions are run for all edge lists.
#' @param method (\emph{optional, default} \code{"binary_classifier"}) \cr
#' A \code{character} string specifying the method used for the HSBM prediction. Options include \code{"binary_classifier"} and \code{"full_reconstruction"}.
#' @param verbose (\emph{optional, default} \code{TRUE}) \cr
#' A \code{logical} value indicating whether to print progress messages during prediction.
#' @param save_blocks (\emph{optional, default} \code{TRUE}) \cr
#' A \code{logical} value indicating whether to save group assignments (blocks) for nodes in the network during the prediction process.
#' @param save_pickle (\emph{optional, default} \code{FALSE}) \cr
#' A \code{logical} value indicating whether to save the model results as Python pickle files for each fold.
#'
#' @return
#' An object of class \code{hsbm.predict} containing the edge/link predictions and group assignments for the specified edge lists (fold):
#' - \code{$data} The binary bipartite \code{matrix} of input data.
#' - \code{$folds} A \code{matrix} of cross-validation fold assignments for each held-out edge/link.
#' - \code{$predictions$probs} A \code{list} where each element is a \code{data.frame} with the predicted probabilities for the edges/links in the corresponding edge list (fold), according to the HSBM model. Each \code{data.frame} contains:
#'   - \code{v1}: The index of the first type of node (rows in the original matrix).
#'   - \code{v2}: The index of the second type of node (columns in the original matrix).
#'   - \code{p}: Predicted probabilities of link between the nodes \code{v1} and \code{v2}.
#'   - \code{v1_names}: Names of the nodes corresponding to \code{v1}, derived from the row names of the original input matrix.
#'   - \code{v2_names}: Names of the nodes corresponding to \code{v2}, derived from the column names of the original input matrix.
#'   - \code{edge_type}: The type of edge/link, such as \code{"documented"} for edges/links observed/documented in the original data or \code{"reconstructed"} for undomented edges/links.
#' - \code{$predictions$groups} (if \code{save_blocks = TRUE}) A \code{list} where each element is a \code{data.frame} containing the group assignments for each node in the network for the corresponding edge list (fold). Each \code{data.frame} includes:
#'   - \code{nodes}: Indices of the nodes in the network.
#'   - \code{G1, G2, G3, ...}: Group assignments for each node across different hierarchical levels, where each column represents a specific level of the hierarchy and the values indicate the group to which the node belongs.
#'   - \code{taxa}: Names of the nodes, derived from the row and column names of the input matrix.   #@@@JMB la columna taxa puede llamrse de una forma más genérica????
#' - \code{$min_dl} A \code{list} of minimum description length values for each fold.
#'
#' @details
#' - The \code{hsbm_input} parameter should be an object of class \code{hsbm.input}, which includes the input data, the cross-validation folds, and corresponding edge lists.
#' - The \code{elist_i} parameter allows you to specify a particular edge list to run predictions on. If not specified, predictions are run on all edge lists.
#' - The \code{method} parameter defines the prediction method to be used. Options are:
#'   \describe{
#'     \item{\code{"binary_classifier"}}{Focuses on predicting probabilities for unobserved links (\code{0s}). Use this to identify missing links in partially incomplete networks.}
#'     \item{\code{"full_reconstruction"}}{Estimates probabilities for all links (\code{0s} and \code{1s}), enabling the identification of both missing and spurious links in incomplete or error-prone networks.}
#'   }
#' - The \code{save_blocks} parameter determines whether group assignments (blocks) for nodes are saved during prediction. Set this to \code{FALSE} to skip saving block information.
#' - The \code{save_pickle} parameter, when \code{TRUE}, saves the model results as Python pickle file. Files are saved in the working directory named as \code{hsbm_res_fold<i>.pkl}, where \code{i} corresponds to the fold index. The pickle object is a dictionary with 5 elements. enum..... #@@@
#'
#' @seealso \code{\link{hsbm.input}}
#'
#' @examples
#' ## Not run:
#' # Example workflow to generate `myInput`:
#' data(dat, package = "sabinaHSBM")
#' 
#' # Prepare input for HSBM
#' myInput <- hsbm.input(data = dat, 
#'                       n_folds = 10, 
#'                       method = "binary_classifier", 
#'                       iter = 1000)
#' ## End(Not run)
#'
#' # Load example HSBM predicted results
#' data(myInput, package = "sabinaHSBM")
#'
#' myPred <- hsbm.predict(hsbm_input = myInput)
#'
#' # View link probabilities of fold 1
#' myPred$predictions$probs[[1]]
#'
#' # View groups of fold 1
#' myPred$predictions$groups[[1]]
#'
#' @export
hsbm.predict <- function(hsbm_input, elist_i = NULL, 
                         verbose = TRUE, save_blocks = TRUE, save_pickle = FALSE){

    if(!inherits(hsbm_input, "hsbm.input")){
        stop("hsbm must be an object of hsbm.input class. Consider running hsbm_input() function.")
    }
    hsbm_name <- as.list(match.call())$hsbm_input

    reticulate::py_run_string(import_modules())
    reticulate::py_run_string(add_taxa_vertex_prop())
    reticulate::py_run_string(create_graph())
    reticulate::py_run_string(get_missing_edges())
    reticulate::py_run_string(hsbm_link_prediction())
    reticulate::py_run_string(get_groups())
    reticulate::py_run_string(get_predicted_network())
    reticulate::py_run_string(save_pickle())

    hsbm_output <- list()
    hsbm_output$data <- hsbm_input$data
    hsbm_output$folds <- hsbm_input$folds
    hsbm_output$method <- hsbm_input$method
    hsbm_output$iter <- hsbm_input$iter
    hsbm_output$wait <- hsbm_input$wait
    hsbm_output$min_dl <- list()

    predictions <- list()
    predictions$probs <- list()
    predictions$groups <- list()
    if(is.null(elist_i)){
        elist_predict <- 1:length(hsbm_input$edgelists)
    }else{
        elist_predict <- elist_i
    }

    for(i in elist_predict){
        if(verbose){
            cat("Computing predictions for fold ", i, "\n")
        }
        i_py <- i - 1
        reticulate::py_run_string(stringr::str_glue("elists = r.{hsbm_name}['edgelists']"))
        reticulate::py_run_string(paste0("elist = elists[", i_py, "]"))
        reticulate::py_run_string("g = create_graph(elist)")
        reticulate::py_run_string(stringr::str_glue("res_dict = hsbm_predict(g, elist, ",
                                                    "method = r.{hsbm_name}['method'], ",
                                                    "force_niter = r.{hsbm_name}['iter'], ",
                                                    "wait = r.{hsbm_name}['wait'])"))
        reticulate::py_run_string("groups_df = get_groups(res_dict['state_min_dl'], res_dict['graph'])")
        if(save_pickle){
            #reticulate::py_save_object(py$res_dict, str_glue("res_dict{i}.pkl"))
            reticulate::py_run_string(paste0("save_pickle(res_dict,", i,")"))
        }

        hsbm_output$predictions$probs[[i]] <- reticulate::py$res_dict$pred_probs
        if(save_blocks){
            hsbm_output$predictions$groups[[i]] <- reticulate::py$groups_df
        }
        hsbm_output$min_dl[[i]] <- reticulate::py$res_dict$min_dl

        reticulate::py_run_string("del res_dict; del groups_df; del g; gc.collect()")

    }

    attr(hsbm_output, "class") <- "hsbm.predict"

    return(hsbm_output)
}



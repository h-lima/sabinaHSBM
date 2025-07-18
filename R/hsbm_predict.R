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
#' @param iter A \code{numeric} value specifying the number of iterations for the HSBM analysis. Default is 10000.
#' @param wait A \code{numeric} value specifying the number of iterations needed for MCMC equilibration. Default is 1000.
#' @param rnd_seed A \code{numeric} value specifying a seed for the random number generator in graph-tool and numpy python modules. Default is \code{NULL} which does not set a random seed.
#' @param verbose (\emph{optional, default} \code{TRUE}) \cr
#' A \code{logical} value indicating whether to print progress messages during prediction.
#' @param save_blocks (\emph{optional, default} \code{TRUE}) \cr
#' A \code{logical} value indicating whether to save group assignments (blocks) for nodes in the network during the prediction process.
#' @param save_pickle (\emph{optional, default} \code{FALSE}) \cr
#' A \code{logical} value indicating whether to save the model results as Python pickle files for each fold.
#' @param save_plots (\emph{optional, default} \code{FALSE}) \cr
#' A \code{logical} value indicating whether to save a hierarchical edge bundling plot for each fold. When \code{TRUE}, a figure is saved in the working directory as a visual representation of the inferred structure.
#' @param n_cores (\emph{optional, default} \code{1}) \cr
#' A \code{numeric} value indicating the number of CPU cores to use for parallel computation; use \code{parallel::detectCores()} to check the maximum available on your system.
#'
#' @return
#' An object of class \code{hsbm.predict} containing the edge/link predictions and group assignments for the specified edge lists (fold):
#' - \code{$data} The binary bipartite \code{matrix} of input data.
#' - \code{$folds} A \code{matrix} of cross-validation fold assignments for each held-out edge/link.
#' - \code{$method} The method used for the HSBM analysis, as specified by the user.
#' - \code{$iter} The number of iterations used to extract link probabilities, as specified by the user. Default is 10000.
#' - \code{$wait} The number of iterations needed to test for equilibration in the \code{mcmc_equilibrate} function from \code{graph-tool}.
#' - \code{$min_dl} A \code{list} of minimum description length values for each fold.
#' - \code{$probs} A \code{list} where each element is a \code{data.frame} with the predicted probabilities for the edges/links in the corresponding edge list (fold), according to the HSBM model. Each \code{data.frame} contains:
#'   - \code{v1}: The index of the first type of node (rows in the original matrix).
#'   - \code{v2}: The index of the second type of node (columns in the original matrix).
#'   - \code{p}: Predicted probabilities of link between the nodes \code{v1} and \code{v2}.
#'   - \code{v1_names}: Names of the nodes corresponding to \code{v1}, derived from the row names of the original input matrix.
#'   - \code{v2_names}: Names of the nodes corresponding to \code{v2}, derived from the column names of the original input matrix.
#'   - \code{edge_type}: The type of edge/link, such as \code{"documented"} for edges/links observed/documented in the original data, \code{"held_out"} for edges/links observed/documented that is retained/masked during cross-validation to assess the model's ability to predict known edges/links, or \code{"reconstructed"} for undocumented edges/links.
#' - \code{$groups} (if \code{save_blocks = TRUE}) A \code{list} where each element is a \code{data.frame} containing the group assignments for each node in the network for the corresponding edge list (fold). Each \code{data.frame} includes:
#'   - \code{nodes}: Indices of the nodes in the network.
#'   - \code{G1, G2, G3, ...}: Group assignments for each node across different hierarchical levels, where each column represents a specific level of the hierarchy and the values indicate the group to which the node belongs.
#'   - \code{names}: Names of the nodes, derived from the row and column names of the input matrix.
#'
#' @details
#' - The \code{hsbm_input} parameter should be an object of class \code{hsbm.input}, which includes the input data, the cross-validation folds, and corresponding edge lists.
#' - The \code{elist_i} parameter allows you to specify a particular edge list to run predictions on. If not specified, predictions are run on all edge lists.
#' - The \code{method} parameter determines the link prediction approach:
#'   \describe{
#'     \item{\code{"binary_classifier"}}{
#'       Computes conditional probabilities for each unobserved/undocumented edges/links (\code{0s}), given the inferred block structure. 
#'       The probability is computed as: \eqn{p_{ij} = \frac{1}{N} \sum_{k=1}^{N} P(A_{ij} = 1 \mid b^{(k)})}, where \eqn{b^{(k)}} is the block partition sampled at iteration \eqn{k}.
#'       This method is particularly useful for identifying missing links—unobserved links that are likely to exist in partially incomplete networks.
#'     }
#'     \item{\code{"full_reconstruction"}}{
#'       Estimates the marginal posterior probabilities for each links (\code{0s} and \code{1s}), by averaging over sampled network configurations. 
#'       The probability is given by: \eqn{p_{ij} = \frac{1}{N} \sum_{k=1}^{N} A_{ij}^{(k)}}, where \eqn{A_{ij}^{(k)} = 1} if the link is present in sample \eqn{k}, and \eqn{0} otherwise.
#'       This method allows identifying both missing links (unobserved/undocumented edges/links likely to exist) and potentially spurious links (observed/documented edges/links that might be erroneous). It is suitable for networks that may contain errors or require a more comprehensive assessment of link validity.
#'     }
#'   }
#' - The \code{save_blocks} parameter determines whether group assignments (blocks) for nodes are saved during prediction. Set this to \code{FALSE} to skip saving block information.
#' - The \code{save_pickle} parameter, when \code{TRUE}, saves the model results as Python pickle file. Files are saved in the working directory named as \code{hsbm_res_fold<i>.pkl}, where \code{i} corresponds to the fold index. The pickle object is a dictionary with 5 elements. enum..... #@@@
#' - If \code{save_plots = TRUE}, hierarchical edge bundling plots are saved in the working directory as PDF files named \code{hsbm_plot_foldi.pdf}, where \code{i} refers to the fold index.
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
#'                       n_folds = 10)
#' ## End(Not run)
#'
#' # Load example HSBM predicted results
#' data(myInput, package = "sabinaHSBM")
#'
#' myPred <- hsbm.predict(hsbm_input = myInput,
#'                       method = "binary_classifier", 
#'                       iter = 1000,
#'                       wait=1000)
#'
#' # View link probabilities of fold 1
#' myPred$probs[[1]]
#'
#' # View groups of fold 1
#' myPred$groups[[1]]
#'
#' @export
hsbm.predict <- function(hsbm_input, elist_i = NULL, 
                         method = "binary_classifier",
                         iter = 10000, wait = 1000,
                         rnd_seed = NULL,
                         verbose = TRUE, 
                         save_blocks = TRUE,
                         save_pickle = FALSE,
                         save_plots = FALSE,
                         n_cores = 1){

    if(!inherits(hsbm_input, "hsbm.input")){
        stop("Error: hsbm_input must be an object of hsbm.input class. Consider running hsbm.input() function.")
    }
    if(!(method %in% c("binary_classifier", "full_reconstruction"))){
        stop("Error: Unknown method. Method should be binary_classifier or",
             " full_reconstruction")
    }

    hsbm_output <- list()
    hsbm_output$data <- hsbm_input$data
    hsbm_output$folds <- hsbm_input$folds
    hsbm_output$method <- method
    hsbm_output$iter <- iter
    hsbm_output$wait <- wait
    hsbm_output$min_dl <- list()
    hsbm_output$probs <- list()
    if(save_blocks){
        hsbm_output$groups <- list()
    }

    if(is.null(elist_i)){
        elist_predict <- 1:length(hsbm_input$edgelists)
    }else{
        elist_predict <- elist_i
    }

    # For n_cores = 1 mclapply calls lapply
    out_lst <- parallel::mclapply(elist_predict,
                                  function(i){
                                     run_hsbm(hsbm_input = hsbm_input,
                                              hsbm_output = hsbm_output,
                                              method = method,
                                              iter = iter,
                                              wait = wait,
                                              rnd_seed = rnd_seed,
                                              save_pickle = save_pickle,
                                              save_blocks = save_blocks,
                                              save_plots = save_plots,
                                              verbose = verbose,
                                              i = i)
                                    }, mc.cores = n_cores)

    hsbm_output <- merge_hsbm.out(out_lst)

    attr(hsbm_output, "class") <- "hsbm.predict"

    return(hsbm_output)
}


merge_hsbm.out <- function(out_lst){ 

    merged <- out_lst[[1]]

    for(el in c("min_dl", "probs", "groups")){
        if(!(el %in% names(merged))) next
        preds <- lapply(out_lst,
                        function(x){
                            res <- unlist(x[[el]],
                                          recursive = FALSE)
                            if(length(res) == 1){
                                return(res)
                            }else{
                               return(data.frame(res))
                            }
                        })

        merged[[el]] <- preds

    }

    return(merged)
}


run_hsbm <- function(hsbm_input, hsbm_output, method, iter, wait, rnd_seed,
                     save_pickle, save_blocks, save_plots, verbose, i){

    reticulate::py_run_string(import_modules())
    reticulate::py_run_string(add_names_vertex_prop())
    reticulate::py_run_string(create_graph())
    reticulate::py_run_string(get_missing_edges())
    reticulate::py_run_string(hsbm_link_prediction())
    reticulate::py_run_string(get_groups())
    reticulate::py_run_string(get_predicted_network())
    reticulate::py_run_string(save_pickle())
    reticulate::py_run_string(save_plot())

    options("reticulate.engine.environment" = environment())

    if(is.null(rnd_seed)){
        rnd_seed <- reticulate::py_none()
    }

    if(verbose){
        cat("Computing predictions for fold ", i, "\n")
    }
    i_py <- i - 1
    reticulate::py_run_string(stringr::str_glue("elists = r.hsbm_input['edgelists']"))
    reticulate::py_run_string(paste0("elist = elists[", i_py, "]"))
    reticulate::py_run_string("g = create_graph(elist)")
    reticulate::py_run_string(stringr::str_glue("res_dict = hsbm_predict(g, elist, ",
                                                "method = r.method, ",
                                                "force_niter = r.iter, ",
                                                "wait = r.wait, ",
                                                "rnd_seed = r.rnd_seed)"))
    reticulate::py_run_string("groups_df = get_groups(res_dict['state_min_dl'], res_dict['graph'])")
    if(save_pickle){
        #reticulate::py_save_object(py$res_dict, str_glue("res_dict{i}.pkl"))
        reticulate::py_run_string(paste0("save_pickle(res_dict,", i,")"))
    }
    if(save_plots){
        reticulate::py_run_string(
                        paste0("save_plot(res_dict['state_min_dl'],", i,")"))
    }

    hsbm_output$probs[[i]] <- reticulate::py$res_dict$pred_probs
    if(save_blocks){
        hsbm_output$groups[[i]] <- reticulate::py$groups_df
    }
    hsbm_output$min_dl[[i]] <- reticulate::py$res_dict$min_dl

    reticulate::py_run_string("del res_dict; del groups_df; del g; gc.collect()")

    options("reticulate.engine.environment" = NULL)

    return(hsbm_output)

}


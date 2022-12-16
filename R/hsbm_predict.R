#' @export
hsbm.predict <- function(hsbm_input, add_n_x = TRUE, verbose = TRUE){
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

    hsbm_output <- hsbm_input

    predictions <- list()
    predictions$probs <- list()
    predictions$groups <- list()
    for(i in 1:length(hsbm_input$edgelist)){
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

        #reticulate::py_save_object(py$res_dict, str_glue("res_dict{i}_bat.pkl"))

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

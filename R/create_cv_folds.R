#' @title Create cross-validation folds with guaranteed connectivity and "maximum" link removal.
#'
#' @description This function identifies the "maximum" number of removable links and partitions
#' them into n_folds. It assigns links to folds to ensure that no
#' node becomes disconnected (i.e., has a degree of zero) in the training
#' matrix of any fold.
#'
#' @param com A binary (0/1) \code{matrix}.
#' @param n_folds (\emph{optional, default} \code{5}) \cr
#' A \code{numeric} value specifying the number of cross-validation folds to generate. If only 1 fold is specified,
#' held-out links will not be computed. This is not recommended, as the main objective of the package is to predict
#' missing links through cross-validation.
#' @param min_per_row A code{numeric} specifying the minimum degree a row node must have for any of its
#'   links to be considered for removal. Default is 2.
#' @param min_per_col A code {numeric} specifying the minimum degree a column node
#' must have for any of its links to be considered for removal. Default is 2.
#' @param max_held_per_fold A code{numeric} specifying the maximum number
#'   of held-out links to place in a fold. If a fold grows larger, its links will be
#'   randomly subsampled to this size. Default is NULL (no limit).
#' @param is_bipartite (\emph{optional, default} \code{TRUE}) \cr
#' A \code{logical} indicating if interactions of node of same type (i.e. nodes in rows and nodes in columns) are to be removed from the inferred probabilities.
#'
#' @return A matrix with columns `row`, `col`, and `fold_id`, indicating
#'   the removed links and their assigned fold for testing.
#'
#' @export
create_cv_folds <- function(com, n_folds = 5, min_per_row = 2, min_per_col = 2,
                            max_held_per_fold = NULL, is_bipartite = TRUE) {

    if(min_per_row < 2 || min_per_col < 2){
        stop("Error: min_per_row and min_per_col must be equal or greater than 2",
             " to guarantee only connected nodes.")
    }
    if(!is.matrix(com)) stop("Error: 'com' must be of type matrix.")
    if(sum(com) ==  0) stop("Error: The provided matrix contains no links.")
    if(!all(com %in% c(0, 1))){
        com[com > 0] <- 1
        warning("Matrix was not binary and has been coerced to 0/1.")
    }
    if(is_bipartite) com[upper.tri(com)] <- 0

    all_links <- which(com == 1, arr.ind = TRUE)

    row_degrees <- rowSums(com)
    col_degrees <- colSums(com)

    # Create list of possible of links to remove
    is_eligible <- apply(all_links, 1,
                         function(link){
                             return(row_degrees[link[1]] >= min_per_row &&
                                    col_degrees[link[2]] >= min_per_col)
                         })
    removable_links <- all_links[is_eligible, , drop = FALSE]

    if (nrow(removable_links) < n_folds){
        stop(paste("Error: Not enough removable links found to create n_folds.",
                   "Found", nrow(removable_links), "but need at least", n_folds))
    }

    shuffled_links <- removable_links[sample(nrow(removable_links)), , drop = FALSE]

    # Track link assignment per fold
    row_node_counts <- matrix(0, nrow = nrow(com), ncol = n_folds)
    col_node_counts <- matrix(0, nrow = ncol(com), ncol = n_folds)
    fold_assignments <- vector(mode = "numeric",
                               length = nrow(shuffled_links))
    fold_sizes <- vector(mode = "numeric",
                         length = n_folds)

    # Assign each link to best possible fold
    for(i in 1:nrow(shuffled_links)){
        link <- shuffled_links[i, ]
        r_node <- link[1]
        c_node <- link[2]

        # Find folds that guarantee node connectivity
        is_valid_fold <- (row_node_counts[r_node, ] < (row_degrees[r_node] - 1)) &
                         (col_node_counts[c_node, ] < (col_degrees[c_node] - 1))
        valid_folds <- which(is_valid_fold)

        if (length(valid_folds) == 0) {
            stop(paste("Could not find a valid fold for link", i,
                       "(row=", r_node, ", col=", c_node, ").",
                       " Re-run function or adjust arguments."))
        }

        # Choose smallest fold to keep them balanced
        best_fold <- valid_folds[which.min(fold_sizes[valid_folds])]

        # Assign link and update trackers
        fold_assignments[i] <- best_fold
        fold_sizes[best_fold] <- fold_sizes[best_fold] + 1
        row_node_counts[r_node, best_fold] <- row_node_counts[r_node, best_fold] + 1
        col_node_counts[c_node, best_fold] <- col_node_counts[c_node, best_fold] + 1
    }

    held_pairs <- cbind(shuffled_links, fold_id = fold_assignments)
    held_pairs <- held_pairs[order(held_pairs[, 3]), ]

    if(!is.null(max_held_per_fold)){
        all_row_ind <- 1:nrow(held_pairs)
        ind_by_fold <- split(all_row_ind, held_pairs[, 3])

        subset_ind_lst <-
            lapply(ind_by_fold,
                function(row_ind) {
                    if (length(row_ind) > max_held_per_fold) {
                        return(sample(row_ind, max_held_per_fold))
                    } else {
                        return(row_ind)
                    }
                }
            )

        subset_ind <- unlist(subset_ind_lst, use.names = FALSE)
        held_pairs <- held_pairs[subset_ind, , drop = FALSE]

    }

    return(held_pairs)
}

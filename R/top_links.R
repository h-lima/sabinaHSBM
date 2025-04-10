#' @name top_links
#'
#' @title Extract Top Predicted Links
#'
#' @description
#' Extracts the top \code{n} predicted links from the HSBM model, allowing for prioritization of the links most likely to be false negatives or false positives. It ranks "undocumented" links (originally 0s) by descending probability, with higher probabilities indicating the most likely false negatives (links that were not observed but are likely true/present). Ranks "documented" links (originally 1s) by ascending probability, with lower probabilities suggesting the most likely false positives (links that were observed but are likely false/absent). Similarly, ranks identified "missing" links and "spurious" links.
#'
#' @param hsbm_reconstructed An object of class \code{hsbm.reconstructed} containing the results of the HSBM model reconstruction.
#' @param n (\emph{optional, default} \code{10}) \cr
#' An \code{integer} specifying the number of top links to extract.
#' @param edge_type (\emph{optional, default} \code{"undocumented"}) \cr
#' A \code{character} indicating the type of edge/link to extract and rank: "undocumented" or "documented", "spurious" or "missing".
#'
#' @return
#' A \code{data.frame} with the top \code{n} predicted links (links with the highest probabilities), including the node names, predicted probabilities (\code{p}), and their standard deviations (\code{sd}) as a uncertainty measure.
#'
#' @seealso \code{\link{hsbm.predict}}, \code{\link{hsbm.reconstructed}}
#'
#' @examples
#' # Load example HSBM prediction results
#' # Assuming `hsbm_reconstructed` is an object of class `hsbm.reconstructed`
#' data(myReconst, package = "sabinaHSBM")
#'
#' top_10_links <- top_links(myPred, n = 10, edge_type = "undocumented")
#' print(top_10_links)
#'
#' @export
top_links <- function(hsbm_reconstructed, n = 10, edge_type = "undocumented") {

    if(!inherits(hsbm_reconstructed, "hsbm.reconstructed")) {
        stop("hsbm_reconstructed must be of class hsbm.reconstructed. Consider running hsbm.reconstructed() function.")
    }

    mat1 <- as.matrix(hsbm_reconstructed$data)
    mat2 <- as.matrix(hsbm_reconstructed$new_mat)

    if(edge_type == "spurious") {
      link_condition <- (mat1 == 1 & mat2 == 0)
    } else if(edge_type == "missing") {
      link_condition <- (mat1 == 0 & mat2 == 1)
    } else if(edge_type == "undocumented") {
      link_condition <- (mat1 == 0)
    } else if(edge_type == "documented") {
      link_condition <- (mat1 == 1)
    } else {
      stop("The 'edge_type' argument must be either 'spurious', 'missing', 'reconstructed' or 'documented'.")
    }

    links <- which(link_condition, arr.ind = TRUE)

    rows <- links[, 1]
    cols <- links[, 2]

    selected_rows <- rownames(mat1)[rows]
    selected_cols <- colnames(mat2)[cols]
    selected_links <- paste(selected_rows, selected_cols, sep = "_")

    res_averaged <- hsbm_reconstructed$reconstructed_df$res_averaged
    res_averaged$combined_names <- paste(res_averaged$v1_names, res_averaged$v2_names, sep = "_")

    reconstructed <- res_averaged[res_averaged$combined_names %in% selected_links, c("v1_names", "v2_names", "p", "sd")]

    if(edge_type %in% c("documented", "spurious")) {
      reconstructed <- reconstructed[order(reconstructed$p), ]
    } else {
      reconstructed <- reconstructed[order(reconstructed$p, decreasing = TRUE), ]
    }

    reconstructed <- reconstructed[1:n, ]

    row.names(reconstructed) <- NULL

    return(reconstructed)
}

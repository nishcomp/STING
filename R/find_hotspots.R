#' find_hotspots
#'
#' Identifies spatial hotspots of gene expression in a Seurat object by classifying spots as hotspots or non-hotspots based on gene expression levels and local neighborhood agreement.
#'
#' @param seurat_object A Seurat object containing spatial transcriptomics data.
#' @param gene_name A character string specifying the gene to analyze for hotspot detection.
#' @param num_neighbors An integer specifying the number of nearest neighbors to consider when assessing local agreement. Default is \code{20}.
#' @param neighbor_fraction A numeric value between 0 and 1 indicating the fraction of neighboring spots that must also be hotspots for a spot to be classified as a hotspot. Default is \code{0.05}.
#' @param percentile A numeric value between 0 and 1 specifying the percentile threshold for defining initial hotspots based on gene expression. Default is \code{0.8}.
#'
#' @return The modified Seurat object with an additional metadata column named \code{hotspot_<gene_name>} indicating whether each spot is classified as a "Hotspot" or "Non-hotspot".
#'
#' @details
#' This function identifies gene expression hotspots by:
#' \itemize{
#'   \item Selecting initial hotspot spots where the gene expression is above the specified percentile threshold.
#'   \item Using k-nearest neighbors (\code{num_neighbors}) to assess whether a spot remains a hotspot based on its neighbors.
#'   \item Classifying each spot as either "Hotspot" or "Non-hotspot" based on local neighborhood agreement.
#' }
#' The classification results are stored in the Seurat object under a new metadata column named \code{hotspot_<gene_name>}.
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(RANN)
#'
#' # Example usage
#' seurat_obj <- Load10XSpatial(data.dir)
#' seurat_obj <- find_hotspots(seurat_object = seurat_obj, gene_name = "GeneX", num_neighbors = 15, neighbor_fraction = 0.1, percentile = 0.85)
#' }
#'
#' @importFrom Seurat FetchData GetTissueCoordinates
#' @importFrom RANN nn2
#' @importFrom stats quantile
#' @export
#'

find_hotspots <- function(seurat_object, gene_name, num_neighbors = 20, neighbor_fraction = 0.05, percentile=0.8){

  gene_expression <- FetchData(seurat_object, vars = gene_name)[, 1]


  expression_threshold <- quantile(gene_expression[gene_expression >0], probs = percentile)

  #initial_hotspot_spots <- which(gene_expression >= mean(gene_expression[gene_expression >0]))
  initial_hotspot_spots <- which(gene_expression >= expression_threshold)
  coordinates <- GetTissueCoordinates(seurat_object, scale = NULL)[, c(1,2)]

  # k-nearest neighbors for each spot using k-d tree
  kdtree <- nn2(coordinates, k = num_neighbors + 1)
  # classification
  regions <- rep("Non-hotspot", length(gene_expression))
  regions[initial_hotspot_spots] <- "Hotspot"

  # checck if nearest neighbors are also hotspots
  for (i in initial_hotspot_spots) {
    nearest_neighbors <- kdtree$nn.idx[i, -1]

    hotspot_neighbors <- sum(regions[nearest_neighbors] == "Hotspot")

    if (hotspot_neighbors / num_neighbors >= neighbor_fraction) {
      # If true, keep this spot as a hotspot
      regions[i] <- "Hotspot"
    } else {
      regions[i] <- "Non-hotspot"
    }
  }
  column_name <- paste0("hotspot_", gene_name)
  seurat_object[[column_name]] <- regions

  return(seurat_object)
}


#' hotspot_neighbors
#'
#'
#' Identifies the neighbors of hotspots by extending upto a specified radius from the hotspot
#' @param seurat_obj A Seurat object containing spatial transcriptomics data.
#' @param distance Distance to extend from hotspots to define neighbors
#' @param metada_column A metadta column of the Seurat Object which has the "hotspot" classification.
#' export
#'
#'
hotspot_neighbors <- function(seurat_obj, distance_threshold, metadata_column) {
  coords <- GetTissueCoordinates(seurat_obj, scale=NULL)

  # hotspot and non-hotspot cells
  regions_data <- FetchData(seurat_obj, vars = metadata_column)
  hotspots <- rownames(regions_data[regions_data[1] == "Hotspot", , drop = FALSE])
  non_hotspots <-  rownames(regions_data[regions_data[1] == "Non-hotspot", , drop = FALSE])

  #  nearest neighbor search
  nn <- nn2(coords[hotspots,], coords[non_hotspots,], k=1)

  # classify non-hotspots based on distance threshold
  neighbors <- non_hotspots[nn$nn.dists <= distance_threshold]
  non_neighbors <- non_hotspots[nn$nn.dists > distance_threshold]

  #  new category labels
  new_categories <- rep("non-neighbors", ncol(seurat_obj))
  names(new_categories) <- colnames(seurat_obj)
  new_categories[hotspots] <- "hotspots"
  new_categories[neighbors] <- "neighbors"

  # Add new category to Seurat object
  column_name <- paste0(metadata_column, "_neighbors")
  seurat_obj[[column_name]] <- new_categories

  return(seurat_obj)
}

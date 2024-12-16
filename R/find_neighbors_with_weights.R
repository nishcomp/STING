#' Calculate Neighbors with Weights
#'
#' @param seurat_obj A Seurat object.
#' @param method Method for finding neighbors ("knns" or "distance").
#' @param k Number of nearest neighbors for "knns".
#' @param dist Maximum distance for "distance".
#' @param bin Resolution bin size.
#' @return A list containing neighbors and weights.
#' @export
find_neighbors_with_weights <- function(seurat_obj, method = c("knns", "distance"), k = 20, dist = 10, bin=8) {

  if (!inherits(seurat_obj, "Seurat")){
    stop("Invalid Seurat object provided")
  }
  if (missing(seurat_obj)) {
    stop("Seurat object is required.")
  }
  if (method == "distance"){
    if (length(Images(seurat_obj))  == 0){
      stop("No image information found in the Seurat object.")
    }
    microns_per_pixel <- bin / seurat_obj@images$slice1@scale.factors$spot
    dist <- dist * microns_per_pixel
  }

  # Get spot coordinates
  spot_coordinates <- GetTissueCoordinates(seurat_obj, scale = NULL)
  sp_points <- sf::st_as_sf(spot_coordinates, coords = c("imagecol", "imagerow"))

  # Determine neighbors
  neighbors_list <- if (method == "knns") {
    neighbors <- knearneigh(spot_coordinates, k = k)
    knn2nb(neighbors)
  } else if (method == "distance") {
    neighbors_list <- dnearneigh(sp_points, d1 = 0, d2 = dist)
  } else {
    stop("Invalid method. Please choose 'knns' or 'distance'.")
  }

  # Compute distance matrix
  dist_matrix <- nbdists(neighbors_list, sf::st_coordinates(sp_points))

  # Inverse distance decay weights
  decay_weights <- lapply(dist_matrix, function(dist) {
    1 / (dist + 1e-16)
  })

  # Adjust for empty neighbor sets
  decay_weights <- lapply(decay_weights, function(weights) {
    if (length(weights) == 0) {
      return(NULL)
    } else {
      return(weights)
    }
  })

  weights <- nb2listw(neighbors_list, glist = decay_weights, style = "B", zero.policy = TRUE)

  return(list(neighbors_list = neighbors_list, weights = weights))
}

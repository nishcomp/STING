#' calculate_moranSVGs
#'
#' Computes spatially variable genes (SVGs) using Moran's I statistic, saves the results to a CSV file if an output path is provided, and returns the results as a data frame.
#'
#' @param sample_name A character string specifying the name of the sample being processed. This will be used in messages and output filenames. Default seurat_obj@project.name.
#' @param seurat_obj A Seurat object containing the spatial transcriptomics data for the sample.
#' @param output_path A character string specifying the directory where the output CSV file will be saved. If \code{NULL}, the results will not be saved to a file. Default is \code{NULL}.
#' @param layer A character string indicating the assay layer to use for expression data. Default is \code{"data"}.
#' @param min_count The minimum number of cells a gene has to be present in (count greater than 0 in a given cell).
#' @param method A character string specifying the method to calculate spatial neighbors. Options are:
#'   \itemize{
#'     \item \code{"knns"}: Uses k-nearest neighbors.
#'     \item \code{"distance"}: Uses a distance threshold to determine neighbors.
#'   }
#'   Default is \code{"knns"}.
#' @param k An integer specifying the number of nearest neighbors to consider when \code{method = "knns"}. Default is \code{20}.
#' @param dist A numeric value specifying the distance threshold (in microns) for neighbor determination when \code{method = "distance"}. Default is \code{10}.
#' @param bin A numeric value representing the binning size in microns per pixel. Default is \code{8}.
#'
#' @return A data frame containing the spatially variable genes and their Moran's I statistics. The data frame includes:
#' \itemize{
#'   \item \code{gene}: Gene name or identifier.
#'   \item \code{I}: Moran's I statistic.
#'   \item \code{p.value}: P-value for the Moran's I test.
#'   \item \code{Adjusted_P}: Adjusted p-value (FDR corrected).
#' }
#'
#' @details
#' This function identifies spatially variable genes by:
#' \itemize{
#'   \item Calculating spatial neighbors based on either k-nearest neighbors (\code{"knns"}) or distance thresholds (\code{"distance"}).
#'   \item Computing Moran's I statistic for each gene's expression values to assess spatial autocorrelation.
#'   \item Adjusting p-values using the False Discovery Rate (FDR) method.
#' }
#' If an \code{output_path} is provided, the results are saved in a CSV file named \code{"<sample_name>_spatially_variable_genes.csv"} in the specified directory.
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(spdep)
#' library(doParallel)
#' library(foreach)
#'
#' # Example usage
#' sample_name <- "Sample_1"
#' sample <- CreateSeuratObject(counts = matrix(rpois(2000, lambda = 10), nrow = 100, ncol = 20))
#' output_path <- "~/output/"
#'
#' results <- calculate_moranSVGs(sample_name = sample_name,
#'                                sample = sample,
#'                                output_path = output_path,
#'                                layer = "data",
#'                                method = "knns",
#'                                k = 15)
#' }
#'
#' @importFrom Seurat GetAssayData GetTissueCoordinates
#' @importFrom parallel detectCores
#' @importFrom spdep knn2nb knearneigh dnearneigh nbdists nb2listw moran.test
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach %dopar% foreach
#' @importFrom stats p.adjust
#' @export
calculate_moranSVGs <- function(seurat_obj, sample_name = seurat_obj@project.name, output_path = NULL, layer = "counts", k = 50, min_count = 10) {

  cat("Processing sample:", sample_name, "\n")

  # Filter genes
  genes_to_keep <- function(obj, min_count) {
    expr_data <- GetAssayData(obj, layer = layer)
    expr_data <- as(expr_data, "dgCMatrix")
    expr_data <- expr_data[!grepl("^MT-", rownames(expr_data)), ]
    gene_totals <- rowSums(expr_data > 0)
    names(gene_totals[gene_totals > min_count])
  }

  genes_to_keep <- genes_to_keep(seurat_obj, min_count)

  # Get spot coordinates
  spot_coordinates <- GetTissueCoordinates(seurat_obj, scale = NULL)
  coord_names <- intersect(c("x", "y", "imagerow", "imagecol"), colnames(spot_coordinates))
  if (length(coord_names) != 2) {
    stop("Expected columns ('x', 'y') or ('imagerow', 'imagecol') not found in spot_coordinates.")
  }
  spot_coordinates <- spot_coordinates[, coord_names]

  # Create spatial points
  sp_points <- sf::st_as_sf(spot_coordinates, coords = coord_names)

  # Find neighbors and calculate weights
  neighbors <- knearneigh(spot_coordinates, k = k)
  neighbors_list <- knn2nb(neighbors)
  dist_matrix <- nbdists(neighbors_list, sp_points)
  decay_weights <- lapply(dist_matrix, function(dist) 1 / (dist + 1e-16))
  weights <- nb2listw(neighbors_list, decay_weights, style = "B")

  cat("Calculated weights for:", seurat_obj@project.name, "\n")

  # Get expression data for filtered genes
  expr_data <- GetAssayData(seurat_obj, layer = layer)[genes_to_keep, ]

  # Parallel computation
  num_cores <- detectCores() - 1
  registerDoParallel(num_cores)

  results <- foreach(i = 1:nrow(expr_data), .combine = 'rbind', .packages = c("spdep")) %dopar% {
    gene_name <- rownames(expr_data)[i]
    moran <- moran.test(expr_data[i, ], weights, alternative = "greater")
    data.frame(row.names = gene_name,
               gene = gene_name,
               I = moran$estimate["Moran I statistic"],
               p.value = moran$p.value,
               stringsAsFactors = FALSE
    )
  }

  stopImplicitCluster()
  gc()

  # Process results
  results <- results %>%
    mutate(Adjusted_P = p.adjust(p.value, method = "fdr")) %>%
    arrange(desc(I))

  # Save results if output_path is provided
  if (!is.null(output_path)) {
    if (!dir.exists(output_path)) {
      dir.create(output_path, recursive = TRUE)
    }
    output_file <- file.path(output_path, paste0(sample_name, "_spatially_variable_genes.csv"))
    write.csv(results, output_file, row.names = FALSE)
    cat("Results saved for:", sample_name, "at", output_file, "\n")
  }

  return(results)
}


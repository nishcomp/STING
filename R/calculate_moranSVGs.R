#' calculate_moranSVGs
#'
#' Computes spatially variable genes (SVGs) using Moran's I statistic, saves the results to a CSV file if an output path is provided, and returns the results as a data frame.
#'
#' @param seurat_obj A Seurat object containing the spatial transcriptomics data for the sample.
#' @param sample_name A character string specifying the name of the sample being processed. This will be used in messages and output filenames. Default seurat_obj@project.name.
#' @param output_path A character string specifying the directory where the output CSV file will be saved. If \code{NULL}, the results will not be saved to a file. Default is \code{NULL}.
#' @param layer A character string indicating the assay layer to use for expression data. Default is \code{"counts"}.
#' @param k An integer specifying the number of nearest neighbors to consider. Default is \code{50}.
#' @param min_count The minimum number of cells a gene has to be present in (count greater than 0 in a given cell). Default is \code{10}.
#'
#' @return A data frame containing the spatially variable genes and their Moran's I statistics.
#'
#' @importFrom Seurat GetAssayData GetTissueCoordinates
#' @importFrom parallel detectCores
#' @importFrom spdep knn2nb knearneigh nbdists nb2listw moran.test
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach %dopar% foreach
#' @importFrom stats p.adjust
#' @export
calculate_moranSVGs <- function(seurat_obj, sample_name = seurat_obj@project.name, output_path = NULL, layer = "counts", k = 50, min_count = 10) {

  cat("Processing sample:", sample_name, "\n")
  cat("Layer used for expression data:", layer, "\n")

  # Filter genes
  genes_to_keep <- function(obj, min_count) {
    expr_data <- GetAssayData(obj, layer = layer)
    expr_data <- as(expr_data, "dgCMatrix")
    expr_data <- expr_data[!grepl("^MT-", rownames(expr_data)), ]
    if (nrow(expr_data) == 0) stop("No genes remaining after filtering.")
    gene_totals <- rowSums(expr_data > 0)
    names(gene_totals[gene_totals > min_count])
  }

  genes_to_keep <- genes_to_keep(seurat_obj, min_count)
  cat("Number of genes retained:", length(genes_to_keep), "\n")

  # Get spot coordinates
  spot_coordinates <- GetTissueCoordinates(seurat_obj, scale = NULL)
  if (is.null(spot_coordinates) || ncol(spot_coordinates) < 2) {
    stop("Spot coordinates could not be retrieved or are insufficient.")
  }
  coord_names <- intersect(c("x", "y", "imagerow", "imagecol"), colnames(spot_coordinates))
  spot_coordinates <- spot_coordinates[, coord_names]

  # Create spatial points
  sp_points <- sf::st_as_sf(spot_coordinates, coords = coord_names)

  # Find neighbors and calculate weights
  neighbors <- knearneigh(spot_coordinates, k = k)
  neighbors_list <- knn2nb(neighbors)
  dist_matrix <- nbdists(neighbors_list, sp_points)
  decay_weights <- lapply(dist_matrix, function(dist) 1 / (dist + 1e-16))
  weights <- nb2listw(neighbors_list, decay_weights, style = "B")
  cat("Calculated spatial weights.\n")

  # Get expression data
  expr_data <- GetAssayData(seurat_obj, layer = layer)[genes_to_keep, ]

  # Parallel computation
  num_cores <- detectCores() - 1
  registerDoParallel(num_cores)
  on.exit(stopImplicitCluster())

  results <- foreach(i = 1:nrow(expr_data), .combine = 'rbind', .packages = c("spdep")) %dopar% {
    tryCatch({
      gene_name <- rownames(expr_data)[i]
      moran <- moran.test(expr_data[i, ], weights, alternative = "greater")
      data.frame(
        row.names = gene_name,
        gene = gene_name,
        I = moran$estimate["Moran I statistic"],
        p.value = moran$p.value,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      message("Error with gene ", rownames(expr_data)[i], ": ", e)
      NULL
    })
  }

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

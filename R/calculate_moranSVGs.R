#' calculate_moranSVGs
#'
#' Computes Moran's I spatial autocorrelation statistic for genes in spatial transcriptomics data.
#'
#' @param expr_data A numeric matrix or data frame where rows represent genes and columns represent spatial spots or cells.
#' @param weights A spatial weights matrix or listw object created using `spdep` for spatial autocorrelation calculations.
#' @param chunk_size Integer, optional (default = 150). The number of genes to process in each chunk. Helps manage memory usage and parallel processing.
#' @param cores Integer, optional (default = `future::availableCores() - 1`). The number of cores to use for parallel processing. The default utilizes all but one available core.
#'
#' @details
#' This function computes Moran's I statistic for each gene in the spatial transcriptomics dataset, testing for spatial autocorrelation.
#' It uses parallel processing via the `future` and `foreach` packages to optimize computations across multiple cores.
#'
#' The results include raw p-values for Moran's I statistic and adjusted p-values corrected using the false discovery rate (FDR) method.
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{I}: The Moran's I statistic for each gene.
#'   \item \code{p.value}: The raw p-value associated with the test of spatial autocorrelation.
#'   \item \code{Adjusted_P}: The FDR-adjusted p-value for each gene.
#' }
#' The rows of the data frame are ordered by Moran's I statistic in descending order.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' library(spdep)
#' library(future)
#' library(foreach)
#'
#' # Simulated spatial transcriptomics data (100 genes, 20 spots)
#' set.seed(123)
#' expr_data <- matrix(rnorm(2000), nrow = 100, ncol = 20)
#' rownames(expr_data) <- paste0("Gene_", 1:100)
#'
#' # Create spatial weights
#' coords <- cbind(runif(20), runif(20)) # Simulated coordinates
#' nb <- knn2nb(knearneigh(coords, k = 5)) # Neighbors based on k=5
#' weights <- nb2listw(nb, style = "W")
#'
#' # Calculate Moran's I
#' results <- calculate_moranSVGs(expr_data, weights)
#' head(results)
#' }

#' @export
calculate_moranSVGs <- function(expr_data, weights, chunk_size = 150, cores = future::availableCores() - 1) {
  plan(multisession, workers = cores)
  chunks <- split(1:nrow(expr_data), ceiling(seq_along(1:nrow(expr_data)) / chunk_size))

  results <- foreach(chunk = chunks, .combine = 'rbind', .options.future = list(seed = TRUE)) %dofuture% {
    do.call(rbind, lapply(chunk, function(i) {
      gene_name <- rownames(expr_data)[i]
      moran <- spdep::moran.test(expr_data[i, ], weights, alternative = "greater")
      data.frame(
        I = moran$estimate["Moran I statistic"],
        p.value = moran$p.value,
        row.names = gene_name
      )
    }))
  }
  results$Adjusted_P <- p.adjust(results$p.value, method = "fdr")
  results[order(results$I, decreasing = TRUE), ]
}

#' plot_svgs
#'
#' Visualizes spatially variable genes (SVGs) as either a circular or linear graph, showing the relationships between genes,
#' their Moran's I values, and their categories.
#'
#' @param df A data frame containing at least the following columns:
#'   \itemize{
#'     \item \code{genes}: Gene names or identifiers.
#'     \item \code{I}: Moran's I statistic values for spatial autocorrelation.
#'     \item \code{category}: Categories or clusters to which genes belong.
#'   }
#' @param title A character string specifying the title of the plot. Default is "Spatially Variable Genes".
#' @param type A character string indicating the type of plot to create. Options are:
#'   \itemize{
#'     \item \code{"circular"}: A circular graph layout.
#'     \item \code{"linear"}: A linear graph layout.
#'   }
#'   Default is \code{"circular"}.
#' @param output_file A character string specifying the path to save the plot as a file (e.g., "plot.png").
#'   If \code{NULL} (the default), the plot is not saved to a file.
#'
#' @return A ggplot object visualizing the spatially variable genes in the chosen graph layout.
#'
#' @details
#' This function creates a graph layout where:
#' \itemize{
#'   \item Nodes represent genes, with sizes corresponding to Moran's I values.
#'   \item Edges connect all genes to each other, except for self-connections.
#'   \item Nodes are colored based on their categories.
#' }
#' Labels are adjusted to prevent overlap and ensure readability.
#'
#' For the circular layout:
#' \itemize{
#'   \item Nodes are arranged in a circular pattern.
#'   \item Labels are placed around the circumference of the circle.
#'   \item Node sizes represent Moran's I values, and colors correspond to gene categories.
#' }
#'
#' For the linear layout:
#' \itemize{
#'   \item Nodes are arranged along the x-axis according to their category and Moran's I values.
#'   \item Labels are placed along the graph, with adjustments to prevent overlap.
#'   \item Node sizes represent Moran's I values, and colors correspond to gene categories.
#' }
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' library(ggraph)
#' library(dplyr)
#' library(RColorBrewer)
#'
#' # Example data
#' df <- data.frame(
#'   genes = paste0("Gene_", 1:10),
#'   I = runif(10, 0, 1),
#'   category = sample(letters[1:3], 10, replace = TRUE)
#' )
#'
#' # Create circular plot
#' p_circular <- plot_svgs(df, title = "Circular SVGs Plot", type = "circular")
#' print(p_circular)
#'
#' # Create linear plot
#' p_linear <- plot_svgs(df, title = "Linear SVGs Plot", type = "linear")
#' print(p_linear)
#' }
#'
#' @import igraph
#' @import ggraph
#' @import dplyr
#' @import RColorBrewer
#' @import ggplot2
#' @import ggrepel
#' @export
plot_svgs <- function(df, title = "Spatially Variable Genes", type = c("circular", "linear"), output_file = NULL) {
  type <- match.arg(type)

  # Validate input data frame
  required_cols <- c("genes", "I", "category")
  if (!all(required_cols %in% colnames(df))) {
    stop("Input data frame must contain columns: 'genes', 'I', and 'category'.")
  }

  if (type == "circular") {
    # Sort data frame
    df <- df %>% arrange(category, desc(I))

    # Create edges and vertices
    edges <- expand.grid(from = df$genes, to = df$genes) %>% filter(from != to)
    vertices <- df %>% mutate(angle = 90 - 360 * (row_number() - 1) / n(),
                              hjust = ifelse(angle < -90, 1, 0),
                              angle = ifelse(angle < -90, angle + 180, angle))

    mygraph <- graph_from_data_frame(edges, vertices = vertices)

    # Define dynamic color palette
    category_colors <- setNames(
      colorRampPalette(brewer.pal(min(9, length(unique(vertices$category))), "Set1"))(length(unique(vertices$category))),
      unique(vertices$category)
    )

    # Generate circular plot
    p <- ggraph(mygraph, layout = 'linear', circular = TRUE) +
      geom_node_point(aes(x = x * 1.05, y = y * 1.05, size = value, color = category), alpha = 0.8) +
      geom_node_text(aes(x = x * 1.2, y = y * 1.2, label = name, angle = angle, hjust = hjust),
                     size = 3, vjust = 0.5) +
      scale_size_continuous(range = c(3, 10)) +
      scale_color_manual(values = category_colors) +
      theme_void() +
      labs(color = "Gene Category", size = "I Value", title = title)
  } else if (type == "linear") {
    df <- df %>%
      group_by(category) %>%
      mutate(position = row_number()) %>%
      ungroup()

    p <- ggplot(df, aes(x = position, y = category)) +
      geom_point(aes(size = I, fill = category), shape = 21, stroke = 1, color = "black") +
      geom_text_repel(aes(label = genes), size = 3, box.padding = 0.5) +
      scale_size_continuous(range = c(2, 8), name = "I Value") +
      scale_fill_brewer(palette = "Set2") +
      theme_bw() +
      labs(title = title)
  }

  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 8, height = 6)
  }

  return(p)
}

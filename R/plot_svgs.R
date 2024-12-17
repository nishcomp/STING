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
#' library(dplyr)
#' library(ggraph)
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
#' import::from(dplyr, arrange, desc, mutate, group_by, row_number, ungroup)
#' import::from(igraph, graph_from_data_frame)
#' import::from(ggraph, ggraph, geom_node_point, geom_node_text)
#' import::from(ggplot2, ggplot, aes, geom_point, scale_size_continuous, scale_fill_brewer, theme_bw, element_blank, labs)
#' import::from(RColorBrewer, brewer.pal)
#' import::from(ggrepel, geom_text_repel)
#' @export
plot_svgs <- function(df, title = "Spatially Variable Genes", type = c("circular", "linear"), output_file = NULL) {
  type <- match.arg(type)

  # Validate input data frame
  required_cols <- c("genes", "I", "category")
  if (!all(required_cols %in% colnames(df))) {
    stop("Input data frame must contain columns: 'genes', 'I', and 'category'.")
  }

  if (type == "circular") {
    # Sort the dataframe by category and I value
    df <- dplyr::arrange(df, category, dplyr::desc(I))

    # Create edges (connecting all genes to each other)
    edges <- data.frame(
      from = rep(df$genes, each = length(df$genes)),
      to = rep(df$genes, times = length(df$genes))
    )

    # Remove self-edges (genes connected to themselves)
    edges <- edges[edges$from != edges$to, ]

    # Create vertices (nodes) without the "center" node
    vertices <- data.frame(
      name = df$genes,
      value = df$I,  # Use 'I' values for node sizes
      category = as.character(df$category)
    )

    # Calculate angles for labels
    n <- nrow(vertices)
    vertices$angle <- 90 - 360 * (seq_len(n) - 1) / n

    # Adjust label alignment
    vertices$hjust <- ifelse(vertices$angle < -90, 1, 0)
    vertices$angle <- ifelse(vertices$angle < -90, vertices$angle + 180, vertices$angle)

    # Create graph object (without the "center" node)
    mygraph <- igraph::graph_from_data_frame(edges, vertices = vertices)

    # Define dynamic color palette
    category_colors <- setNames(
      colorRampPalette(RColorBrewer::brewer.pal(min(9, length(unique(vertices$category))), "Set1"))(length(unique(vertices$category))),
      unique(vertices$category)
    )

    # Generate circular plot
    # Make the plot
    p <- ggraph::ggraph(mygraph, layout = 'linear', circular = TRUE) +
      ggraph::geom_node_point(aes(x = x * 1.05, y = y * 1.05, size = value, color = category), alpha = 0.8) +
      ggraph::geom_node_text(aes(x = x * 1.2, y = y * 1.2, label = name, angle = angle, hjust = hjust),
                             size = 3, vjust = 0.5) +
      ggplot2::scale_size_continuous(range = c(3, 10)) +
      ggplot2::scale_color_manual(values = category_colors) +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.position = "right",
        plot.margin = unit(c(1,1,1,1), "cm"),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14)
      ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(title = NULL),
        size = ggplot2::guide_legend(title = "Moran's I Value")
      ) +
      ggplot2::expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5)) +
      ggplot2::coord_fixed() +
      ggplot2::labs(color = "Gene Category", size = "I Value", title = title)

  } else if (type == "linear") {
    df <- dplyr::group_by(df, category) %>%
      dplyr::mutate(position = dplyr::row_number()) %>%
      dplyr::ungroup()

    p <- ggplot2::ggplot(df, ggplot2::aes(x = position, y = category)) +
      ggplot2::geom_point(ggplot2::aes(size = I, fill = category), shape = 21, stroke = 1, color="black") + # Filled circles with black border
      ggrepel::geom_text_repel(
        ggplot2::aes(label = genes),
        size = 3,
        box.padding = 0.5, # Adjust padding around labels
        point.padding = 0.3,
        segment.color = "grey50",
        min.segment.length = 0,
        max.overlaps = Inf
      ) +
      ggplot2::scale_size_continuous(range = c(2, 8), name = "I Value") +  # Adjust size range and add a legend title
      ggplot2::scale_fill_brewer(palette = "Set2", name = "Category") +
      ggplot2::theme_bw() +
      ggplot2::theme(title = ggplot2::element_text(face = "bold", hjust = "center"),
                     panel.grid = ggplot2::element_blank(),
                     axis.title = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(face = "bold", size = 12),
                     legend.position = "top",
                     panel.border = ggplot2::element_blank(),
                     axis.line.y = ggplot2::element_line(color = "black")
      ) +
      ggplot2::guides(
        fill = "none",
        size = ggplot2::guide_legend(title = "Moran's I Value")
      ) +
      ggplot2::labs(title = title)
  }

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = 8, height = 6)
  }

  return(p)
}


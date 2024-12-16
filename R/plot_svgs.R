#' plot_spatial_genes
#'
#' Visualizes spatially variable genes (SVGs) as a circular graph, showing the relationships between genes,
#' their Moran's I values, and their categories.
#'
#' @param df A data frame containing at least the following columns:
#'   \itemize{
#'     \item \code{genes}: Gene names or identifiers.
#'     \item \code{I}: Moran's I statistic values for spatial autocorrelation.
#'     \item \code{category}: Categories or clusters to which genes belong.
#'   }
#' @param title A character string specifying the title of the plot. Default is "Spatially Variable Genes".
#'
#' @return A ggplot object visualizing the spatially variable genes in a circular graph layout.
#'
#' @details
#' This function creates a circular graph layout where:
#' \itemize{
#'   \item Nodes represent genes, with sizes corresponding to Moran's I values.
#'   \item Edges connect all genes to each other, except for self-connections.
#'   \item Nodes are colored based on their categories.
#' }
#' Labels are adjusted to prevent overlap and ensure readability in the circular layout.
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
#' # Create plot
#' p <- plot_spatial_genes(df, title = "Example SVGs Plot")
#' print(p)
#' }
#'
#' @import igraph
#' @import ggraph
#' @import dplyr
#' @import RColorBrewer
#' @import ggplot2
#' @import ggrepl
#' @export
plot_svgs <- function(df, title = "Spatially Variable Genes", type = c("circular", "linear")) {
  if (type == "circular"){
    # Sort the data frame by category and Moran's I value
    df <- df %>% arrange(category, desc(I))

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
    mygraph <- graph_from_data_frame(edges, vertices = vertices)

    # Define a color palette for categories
    category_colors <- setNames(brewer.pal(length(unique(vertices$category)), "Set1"),
                                unique(vertices$category))

    # Create the plot
    p <- ggraph(mygraph, layout = 'linear', circular = TRUE) +
      geom_node_point(aes(x = x * 1.05, y = y * 1.05, size = value, color = category), alpha = 0.8) +
      geom_node_text(aes(x = x * 1.2, y = y * 1.2, label = name, angle = angle, hjust = hjust),
                     size = 3, vjust = 0.5) +
      scale_size_continuous(range = c(3, 10)) +
      scale_color_manual(values = category_colors) +
      theme_void() +
      theme(
        legend.position = "right",
        plot.margin = unit(c(1,1,1,1), "cm"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
      ) +
      guides(
        color = guide_legend(title = NULL),
        size = guide_legend(title = "Moran's I Value")
      ) +
      expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5)) +
      coord_fixed() +
      labs(color = "Gene Category", size = "I Value", title = title)
  }
  if (type == "linear"){
    df <- df %>%
      group_by(category) %>%
      mutate(position = row_number()) %>%
      ungroup()

    p <- ggplot(df, aes(x = position, y = category)) +
      geom_point(aes(size = I, fill = category), shape = 21, stroke = 1, color="black") + # Filled circles with black border
      geom_text_repel(
        aes(label = genes),
        size = 3,
        box.padding = 0.5, # Adjust padding around labels
        point.padding = 0.3,
        segment.color = "grey50",
        min.segment.length = 0,
        max.overlaps = Inf
      ) +
      scale_size_continuous(range = c(2, 8), name = "I Value") +  # Adjust size range and add a legend title
      scale_fill_brewer(palette = "Set2", name = "Category") +
      theme_bw() +
      theme(title = element_text(face = "bold", hjust = "center"),
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_text(face = "bold", size = 12),
            legend.position = "top",
            panel.border = element_blank(),
            axis.line.y = element_line(color = "black")

      ) +
      guides(
        fill = "none",
        size = guide_legend(title = "Moran's I Value")
      ) +
      labs(title = "Gene Expression by Category")
  }

  return(p)
}

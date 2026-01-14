#' Network Graph Visualization for Pleiotropic SNPs
#'
#' Creates an interactive network graph showing connections between pleiotropic SNPs
#' and their associated traits. Nodes represent SNPs and traits, edges represent
#' significant associations.
#'
#' @param pleio_data A data.frame containing pleiotropy analysis results
#' @param min_n_traits Integer. Minimum number of traits per SNP to include (default: 2)
#' @param top_n_snps Integer. Number of top pleiotropic SNPs to highlight (default: 10)
#' @param node_size_snp Numeric. Base size for SNP nodes (default: 5)
#' @param node_size_trait Numeric. Base size for trait nodes (default: 8)
#' @param edge_width_max Numeric. Maximum edge width (default: 3)
#' @param layout Character. Layout algorithm: "fr", "kk", "dh", "grid", "circle" (default: "fr")
#' @param show_labels Logical. Show node labels (default: TRUE)
#' @param label_size Numeric. Label font size (default: 3)
#' @param color_palette Character. Color palette: "okabe_ito", "viridis", "default" (default: "okabe_ito")
#'
#' @return A ggraph ggplot2 object
#'
#' @import ggplot2
#' @importFrom dplyr group_by summarise mutate filter arrange slice_head ungroup rename
#' @importFrom igraph graph_from_data_frame set_vertex_attr
#' @import ggraph
#'
#' @examples
#' data(gwas_subset)
#' pleio_results <- detect_pleiotropy(gwas_subset)
#' if (nrow(pleio_results) > 0) {
#'     p <- plot_pleiotropy_network(pleio_results, top_n_snps = 5)
#'     print(p)
#' }
#'
#' @export
plot_pleiotropy_network <- function(pleio_data,
                                    min_n_traits = 2,
                                    top_n_snps = 10,
                                    node_size_snp = 5,
                                    node_size_trait = 8,
                                    edge_width_max = 3,
                                    layout = "fr",
                                    show_labels = TRUE,
                                    label_size = 3,
                                    color_palette = "okabe_ito") {
    if (!is.data.frame(pleio_data)) {
        stop("Input must be a data.frame or data.table")
    }

    required_cols <- c("SNPS", "MAPPED_TRAIT", "PVALUE_MLOG")
    missing_cols <- setdiff(required_cols, names(pleio_data))
    if (length(missing_cols) > 0) {
        stop("Required columns missing: ", paste(missing_cols, collapse = ", "))
    }

    if (nrow(pleio_data) == 0) {
        stop("Input data is empty")
    }

    snp_traits <- pleio_data |>
        group_by(SNPS, MAPPED_TRAIT) |>
        summarise(
            MAX_PVALUE = max(PVALUE_MLOG),
            .groups = "drop"
        )

    snp_counts <- snp_traits |>
        group_by(SNPS) |>
        summarise(N_TRAITS = n(), .groups = "drop") |>
        filter(N_TRAITS >= min_n_traits) |>
        arrange(desc(N_TRAITS)) |>
        slice_head(n = top_n_snps)

    if (nrow(snp_counts) == 0) {
        warning("No SNPs meet the minimum trait criteria")
        return(ggplot() +
            theme_void() +
            annotate("text",
                x = 0.5, y = 0.5,
                label = "No data to display\nAdjust filtering criteria"
            ))
    }

    top_snps <- snp_counts$SNPS
    plot_edges <- snp_traits |>
        filter(SNPS %in% top_snps) |>
        rename(from = SNPS, to = MAPPED_TRAIT, weight = MAX_PVALUE)

    g <- igraph::graph_from_data_frame(plot_edges, directed = FALSE)

    if ("CHR_ID" %in% names(pleio_data)) {
        snp_chrom <- pleio_data |>
            group_by(SNPS) |>
            summarise(CHR = dplyr::first(CHR_ID), .groups = "drop") |>
            filter(SNPS %in% top_snps)

        g <- igraph::set_vertex_attr(g, "type",
            value = ifelse(igraph::V(g)$name %in% snp_chrom$SNP, "SNP", "Trait")
        )
        g <- igraph::set_vertex_attr(g, "n_traits",
            value = ifelse(igraph::V(g)$name %in% snp_counts$SNP,
                snp_counts$N_TRAITS[match(igraph::V(g)$name, snp_counts$SNP)], 1
            )
        )
    } else {
        g <- igraph::set_vertex_attr(g, "type",
            value = ifelse(igraph::V(g)$name %in% top_snps, "SNP", "Trait")
        )
        g <- igraph::set_vertex_attr(g, "n_traits",
            value = ifelse(igraph::V(g)$name %in% top_snps,
                snp_counts$N_TRAITS[match(igraph::V(g)$name, snp_counts$SNP)], 1
            )
        )
    }

    colors <- get_pleiotropy_colors(color_palette, n = 10)

    p <- ggraph(g, layout = layout) +
        geom_edge_link(aes(width = weight),
            alpha = 0.6,
            color = "gray50"
        ) +
        scale_edge_width_continuous(
            range = c(0.5, edge_width_max),
            guide = "none"
        ) +
        geom_node_point(aes(color = type, size = n_traits), alpha = 0.8) +
        scale_color_manual(values = c("SNP" = colors[1], "Trait" = colors[2])) +
        scale_size_continuous(
            range = c(node_size_snp, node_size_trait),
            guide = "none"
        ) +
        labs(
            title = "Pleiotropic SNP-Trait Network",
            subtitle = paste("Top", top_n_snps, "pleiotropic SNPs and their associated traits"),
            color = "Node Type"
        ) +
        theme_void() +
        theme(
            legend.position = "bottom",
            legend.title = element_text(size = 10, face = "bold"),
            legend.text = element_text(size = 9),
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 10, hjust = 0.5)
        )

    if (show_labels) {
        p <- p + geom_node_text(aes(label = name),
            repel = TRUE,
            size = label_size,
            color = "black",
            bg.color = "white",
            bg.r = 0.15,
            segment.color = "gray50",
            min.segment.length = 0
        )
    }

    return(p)
}

#' Trait Co-occurrence Network Visualization
#'
#' Creates a network showing how often pairs of traits share pleiotropic SNPs.
#' This helps identify trait clusters and shared biological mechanisms.
#'
#' @param pleio_data A data.frame containing pleiotropy analysis results
#' @param min_shared_snps Integer. Minimum shared SNPs to show edge (default: 1)
#' @param top_n_traits Integer. Number of top traits to include (default: 15)
#' @param node_size Numeric. Base node size (default: 8)
#' @param layout Character. Layout algorithm: "fr", "kk", "dh", "grid", "circle" (default: "fr")
#' @param show_edge_labels Logical. Show shared SNP counts on edges (default: TRUE)
#' @param color_by Character. Color nodes by: "n_snps", "degree", "none" (default: "n_snps")
#'
#' @return A ggraph ggplot2 object
#'
#' @import ggplot2
#' @importFrom dplyr group_by summarise mutate filter arrange slice_head ungroup rename
#' @importFrom igraph graph_from_data_frame
#' @import ggraph
#' @importFrom ggrepel geom_text_repel
#'
#' @examples
#' data(gwas_subset)
#' pleio_results <- detect_pleiotropy(gwas_subset)
#' if (nrow(pleio_results) > 0) {
#'     p <- plot_trait_cooccurrence_network(pleio_results, top_n_traits = 5)
#'     print(p)
#' }
#'
#' @export
plot_trait_cooccurrence_network <- function(pleio_data,
                                            min_shared_snps = 1,
                                            top_n_traits = 15,
                                            node_size = 8,
                                            layout = "fr",
                                            show_edge_labels = TRUE,
                                            color_by = "n_snps") {
    if (!is.data.frame(pleio_data)) {
        stop("Input must be a data.frame or data.table")
    }

    if (nrow(pleio_data) == 0) {
        stop("Input data is empty")
    }

    trait_counts <- pleio_data |>
        group_by(MAPPED_TRAIT) |>
        summarise(N_SNPS = dplyr::n_distinct(SNPS), .groups = "drop") |>
        arrange(desc(N_SNPS)) |>
        slice_head(n = top_n_traits)

    if (nrow(trait_counts) < 2) {
        stop("Need at least 2 traits to create network")
    }

    plot_data <- pleio_data |>
        filter(MAPPED_TRAIT %in% trait_counts$MAPPED_TRAIT)

    trait_pairs <- plot_data |>
        group_by(MAPPED_TRAIT) |>
        summarise(SNPS = list(unique(SNPS)), .groups = "drop") |>
        mutate(N_SNPS = lengths(SNPS))

    cooccurrence <- expand.grid(
        trait1 = trait_pairs$MAPPED_TRAIT,
        trait2 = trait_pairs$MAPPED_TRAIT,
        stringsAsFactors = FALSE
    ) |>
        filter(trait1 < trait2) |>
        mutate(
            shared_snps = mapply(function(t1, t2) {
                length(intersect(
                    trait_pairs$SNPS[trait_pairs$MAPPED_TRAIT == t1],
                    trait_pairs$SNPS[trait_pairs$MAPPED_TRAIT == t2]
                ))
            }, trait1, trait2)
        ) |>
        filter(shared_snps >= min_shared_snps)

    if (nrow(cooccurrence) == 0) {
        warning("No trait pairs meet minimum shared SNP criteria")
        return(ggplot() +
            theme_void() +
            annotate("text",
                x = 0.5, y = 0.5,
                label = "No trait co-occurrences found\nAdjust filtering criteria"
            ))
    }

    edges <- cooccurrence |>
        rename(from = trait1, to = trait2, weight = shared_snps)

    g <- igraph::graph_from_data_frame(edges, directed = FALSE)

    node_data <- trait_counts |>
        rename(name = MAPPED_TRAIT, n_snps = N_SNPS)

    colors <- get_pleiotropy_colors("viridis", n = 10)

    p <- ggraph(g, layout = layout) +
        geom_edge_link(aes(width = weight),
            alpha = 0.6,
            color = "gray50"
        ) +
        scale_edge_width_continuous(
            range = c(0.5, 4),
            guide = guide_legend(title = "Shared SNPs")
        ) +
        geom_node_point(aes(size = n_snps, color = n_snps), alpha = 0.8) +
        scale_size_continuous(
            range = c(node_size, node_size * 2),
            guide = guide_legend(title = "Total SNPs")
        ) +
        scale_color_gradientn(
            colors = colors,
            guide = guide_legend(title = "Total SNPs")
        ) +
        labs(
            title = "Trait Co-occurrence Network",
            subtitle = "Edges represent number of shared pleiotropic SNPs"
        ) +
        theme_void() +
        theme(
            legend.position = "right",
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 10, hjust = 0.5)
        )

    if (show_edge_labels) {
        p <- p + geom_edge_link(aes(label = weight),
            label_dodge = 0.1,
            vjust = -0.5,
            size = 2.5,
            color = "black",
            angle_calc = "along"
        )
    }

    return(p)
}

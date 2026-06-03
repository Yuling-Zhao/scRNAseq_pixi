mean_by_column <- function(matrix) {
  as.numeric(Matrix::colMeans(matrix))
}

fraction_positive_by_column <- function(matrix) {
  as.numeric(Matrix::colMeans(matrix > 0))
}

cluster_order <- function(metadata, cluster_key) {
  clusters <- metadata[[cluster_key]]
  if (is.factor(clusters)) {
    return(as.character(levels(clusters)))
  }

  values <- unique(as.character(clusters))
  numeric_values <- suppressWarnings(as.integer(values))
  if (all(!is.na(numeric_values))) {
    values[order(numeric_values)]
  } else {
    sort(values)
  }
}

get_expression_matrix <- function(object, assay = NULL, slot = "data") {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("The R implementation expects a Seurat object and requires the Seurat package.")
  }

  assay <- assay %||% Seurat::DefaultAssay(object)
  Seurat::GetAssayData(object = object, assay = assay, slot = slot)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

semi_automated_marker_annotation <- function(
    object,
    cluster_key,
    marker_gene_sets_json,
    annotation_dir = "results/annotation",
    marker_gene_dir = "data/marker_gene",
    max_global_fraction = 0.60,
    min_global_fraction = 0.01,
    min_specificity_score = 1.5,
    min_retained_markers = 2,
    high_confidence_min_score = 0.75,
    high_confidence_min_margin = 0.20,
    medium_confidence_min_score = 0.50,
    manual_review_margin_threshold = 0.15,
    embedding_key = "umap_harmony_pcs30_nn30",
    max_markers_per_cell_type = 5,
    dotplot_standard_scale = "var",
    dotplot_swap_axes = TRUE,
    dotplot_dendrogram = FALSE,
    dotplot_dpi = 400,
    umap_dpi = 300,
    output_prefix = "",
    assay = NULL,
    slot = "data") {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required.")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.")
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required.")
  }

  dir.create(annotation_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(marker_gene_dir, recursive = TRUE, showWarnings = FALSE)

  metadata <- object[[]]
  if (!cluster_key %in% colnames(metadata)) {
    stop(sprintf("'%s' is not present in object metadata", cluster_key))
  }

  marker_gene_sets <- jsonlite::fromJSON(marker_gene_sets_json, simplifyVector = FALSE)
  if (!is.list(marker_gene_sets) || length(marker_gene_sets) == 0) {
    stop("marker_gene_sets_json must contain a non-empty JSON object")
  }
  marker_gene_sets <- lapply(marker_gene_sets, function(markers) unique(as.character(markers)))

  prefix <- if (nzchar(output_prefix)) paste0(output_prefix, "_") else ""
  filtered_marker_path <- file.path(marker_gene_dir, paste0(prefix, "filtered_marker_gene_sets_core.json"))
  marker_presence_path <- file.path(annotation_dir, paste0(prefix, "marker_presence_summary.csv"))
  filtration_metrics_path <- file.path(annotation_dir, paste0(prefix, "marker_filtration_metrics.csv"))
  cluster_scores_path <- file.path(annotation_dir, paste0(prefix, "cluster_celltype_marker_scores.csv"))
  suggested_annotation_path <- file.path(annotation_dir, paste0(prefix, "suggested_cluster_annotation.csv"))
  validation_dotplot_path <- file.path(annotation_dir, paste0(prefix, "filtered_marker_validation_dotplot.png"))
  validation_dotplot_svg_path <- file.path(annotation_dir, paste0(prefix, "filtered_marker_validation_dotplot.svg"))
  suggested_umap_path <- file.path(annotation_dir, paste0(prefix, "suggested_annotation_umap.png"))
  suggested_umap_svg_path <- file.path(annotation_dir, paste0(prefix, "suggested_annotation_umap.svg"))
  manual_review_path <- file.path(annotation_dir, paste0(prefix, "clusters_for_manual_review.csv"))

  clusters <- as.character(metadata[[cluster_key]])
  names(clusters) <- rownames(metadata)
  ordered_clusters <- cluster_order(metadata, cluster_key)
  expression_matrix <- get_expression_matrix(object, assay = assay, slot = slot)

  marker_gene_sets_present <- lapply(marker_gene_sets, function(genes) genes[genes %in% rownames(expression_matrix)])
  marker_genes_missing <- lapply(marker_gene_sets, function(genes) genes[!genes %in% rownames(expression_matrix)])

  marker_presence_summary <- do.call(rbind, lapply(names(marker_gene_sets), function(cell_type) {
    data.frame(
      cell_type = cell_type,
      total_markers = length(marker_gene_sets[[cell_type]]),
      present_markers = length(marker_gene_sets_present[[cell_type]]),
      missing_markers = length(marker_genes_missing[[cell_type]]),
      missing_marker_genes = paste(marker_genes_missing[[cell_type]], collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(marker_presence_summary, marker_presence_path, row.names = FALSE)

  present_markers <- unique(unlist(marker_gene_sets_present, use.names = FALSE))
  if (length(present_markers) == 0) {
    stop("No marker genes from marker_gene_sets_json were found in the expression matrix")
  }

  x_markers <- expression_matrix[present_markers, , drop = FALSE]
  global_mean_expression <- setNames(mean_by_column(t(x_markers)), present_markers)
  global_fraction_expressing <- setNames(fraction_positive_by_column(t(x_markers)), present_markers)

  cluster_mean_expression <- matrix(
    NA_real_,
    nrow = length(ordered_clusters),
    ncol = length(present_markers),
    dimnames = list(ordered_clusters, present_markers)
  )
  cluster_fraction_expressing <- cluster_mean_expression
  for (cluster in ordered_clusters) {
    cell_names <- names(clusters)[clusters == cluster]
    x_cluster <- expression_matrix[present_markers, cell_names, drop = FALSE]
    cluster_mean_expression[cluster, ] <- mean_by_column(t(x_cluster))
    cluster_fraction_expressing[cluster, ] <- fraction_positive_by_column(t(x_cluster))
  }

  metric_rows <- list()
  for (cell_type in names(marker_gene_sets)) {
    for (marker in marker_gene_sets[[cell_type]]) {
      if (!marker %in% present_markers) {
        metric_rows[[length(metric_rows) + 1]] <- data.frame(
          cell_type = cell_type,
          marker_gene = marker,
          marker_present = FALSE,
          global_mean_expression = NA_real_,
          global_fraction_expressing = NA_real_,
          max_cluster_mean = NA_real_,
          second_cluster_mean = NA_real_,
          specificity_score = NA_real_,
          cluster_separation = NA_real_,
          top_cluster = NA_character_,
          broadly_expressed = NA,
          low_detection = NA,
          low_specificity = NA,
          keep_marker = FALSE,
          stringsAsFactors = FALSE
        )
        next
      }

      cluster_means <- as.numeric(cluster_mean_expression[, marker])
      names(cluster_means) <- ordered_clusters
      sorted_means <- sort(cluster_means, decreasing = TRUE)
      max_cluster_mean <- sorted_means[[1]]
      second_cluster_mean <- if (length(sorted_means) > 1) sorted_means[[2]] else 0
      specificity_score <- max_cluster_mean / (global_mean_expression[[marker]] + 1e-6)
      cluster_separation <- max_cluster_mean - second_cluster_mean
      broadly_expressed <- global_fraction_expressing[[marker]] > max_global_fraction
      low_detection <- global_fraction_expressing[[marker]] < min_global_fraction
      low_specificity <- specificity_score < min_specificity_score
      keep_marker <- !(broadly_expressed || low_detection || low_specificity)

      row <- data.frame(
        cell_type = cell_type,
        marker_gene = marker,
        marker_present = TRUE,
        global_mean_expression = global_mean_expression[[marker]],
        global_fraction_expressing = global_fraction_expressing[[marker]],
        max_cluster_mean = max_cluster_mean,
        second_cluster_mean = second_cluster_mean,
        specificity_score = specificity_score,
        cluster_separation = cluster_separation,
        top_cluster = names(sorted_means)[[1]],
        broadly_expressed = broadly_expressed,
        low_detection = low_detection,
        low_specificity = low_specificity,
        keep_marker = keep_marker,
        stringsAsFactors = FALSE
      )
      for (cluster in ordered_clusters) {
        row[[paste0("cluster_", cluster, "_mean_expression")]] <- cluster_mean_expression[cluster, marker]
        row[[paste0("cluster_", cluster, "_fraction_expressing")]] <- cluster_fraction_expressing[cluster, marker]
      }
      metric_rows[[length(metric_rows) + 1]] <- row
    }
  }

  marker_filtration_metrics <- dplyr::bind_rows(metric_rows)
  write.csv(marker_filtration_metrics, filtration_metrics_path, row.names = FALSE)

  keep_idx <- marker_filtration_metrics$keep_marker %in% TRUE
  kept_marker_keys <- paste(
    marker_filtration_metrics$cell_type[keep_idx],
    marker_filtration_metrics$marker_gene[keep_idx],
    sep = "\r"
  )
  marker_gene_sets_filtered <- lapply(names(marker_gene_sets_present), function(cell_type) {
    markers <- marker_gene_sets_present[[cell_type]]
    markers[paste(cell_type, markers, sep = "\r") %in% kept_marker_keys]
  })
  names(marker_gene_sets_filtered) <- names(marker_gene_sets_present)
  marker_gene_sets_filtered <- marker_gene_sets_filtered[lengths(marker_gene_sets_filtered) >= min_retained_markers]
  jsonlite::write_json(marker_gene_sets_filtered, filtered_marker_path, pretty = TRUE, auto_unbox = TRUE)

  retained_markers <- unique(unlist(marker_gene_sets_filtered, use.names = FALSE))
  score_rows <- list()
  if (length(retained_markers) > 0) {
    cluster_average_expression <- cluster_mean_expression[ordered_clusters, retained_markers, drop = FALSE]
    gene_means <- colMeans(cluster_average_expression)
    gene_stds <- apply(
      cluster_average_expression,
      2,
      function(values) sqrt(mean((values - mean(values))^2))
    )
    gene_stds[gene_stds == 0 | is.na(gene_stds)] <- 1
    marker_zscores <- sweep(sweep(cluster_average_expression, 2, gene_means, "-"), 2, gene_stds, "/")

    for (cluster in ordered_clusters) {
      raw_cluster_rows <- lapply(names(marker_gene_sets_filtered), function(cell_type) {
        markers <- marker_gene_sets_filtered[[cell_type]]
        marker_scores <- as.numeric(marker_zscores[cluster, markers])
        mean_scaled_marker_expression <- mean(marker_scores)
        fraction_positive_markers <- mean(marker_scores > 0)
        specificity_adjusted_score <- mean_scaled_marker_expression * fraction_positive_markers
        data.frame(
          cluster = cluster,
          cell_type = cell_type,
          n_markers_used = length(markers),
          mean_scaled_marker_expression = mean_scaled_marker_expression,
          fraction_positive_markers = fraction_positive_markers,
          specificity_adjusted_score = specificity_adjusted_score,
          markers_used = paste(markers, collapse = ";"),
          stringsAsFactors = FALSE
        )
      })
      raw_cluster_rows <- dplyr::bind_rows(raw_cluster_rows)
      raw_scores <- raw_cluster_rows$specificity_adjusted_score
      score_min <- if (length(raw_scores)) min(raw_scores) else 0
      score_max <- if (length(raw_scores)) max(raw_scores) else 0
      score_range <- score_max - score_min
      raw_cluster_rows$possibility_score <- if (score_range > 0) {
        (raw_scores - score_min) / score_range
      } else if (score_max > 0) {
        1
      } else {
        0
      }
      score_rows[[length(score_rows) + 1]] <- raw_cluster_rows
    }
  }

  cluster_celltype_marker_scores <- dplyr::bind_rows(score_rows)
  write.csv(cluster_celltype_marker_scores, cluster_scores_path, row.names = FALSE)

  annotation_rows <- list()
  if (nrow(cluster_celltype_marker_scores) > 0) {
    for (cluster in ordered_clusters) {
      cluster_scores <- cluster_celltype_marker_scores[cluster_celltype_marker_scores$cluster == cluster, , drop = FALSE]
      cluster_scores <- cluster_scores[order(cluster_scores$possibility_score, cluster_scores$specificity_adjusted_score, decreasing = TRUE), , drop = FALSE]
      if (nrow(cluster_scores) == 0) next

      best <- cluster_scores[1, , drop = FALSE]
      second <- if (nrow(cluster_scores) > 1) cluster_scores[2, , drop = FALSE] else NULL
      second_best_cell_type <- if (is.null(second)) NA_character_ else second$cell_type
      second_best_score <- if (is.null(second)) 0 else second$possibility_score
      confidence_margin <- best$possibility_score - second_best_score

      confidence_label <- if (
        best$possibility_score >= high_confidence_min_score &&
          confidence_margin >= high_confidence_min_margin
      ) {
        "high"
      } else if (best$possibility_score >= medium_confidence_min_score) {
        "medium"
      } else {
        "low"
      }

      explanation <- sprintf(
        "Cluster %s is suggested as %s based on retained markers %s with marker score %.3f and possibility score %.3f. The competing annotation is %s with score %.3f; confidence margin is %.3f.",
        cluster,
        best$cell_type,
        best$markers_used,
        best$specificity_adjusted_score,
        best$possibility_score,
        second_best_cell_type,
        second_best_score,
        confidence_margin
      )

      annotation_rows[[length(annotation_rows) + 1]] <- data.frame(
        cluster = cluster,
        suggested_cell_type = best$cell_type,
        possibility_score = best$possibility_score,
        confidence_label = confidence_label,
        confidence_margin = confidence_margin,
        markers_used = best$markers_used,
        n_markers_used = best$n_markers_used,
        mean_scaled_marker_expression = best$mean_scaled_marker_expression,
        fraction_positive_markers = best$fraction_positive_markers,
        specificity_adjusted_score = best$specificity_adjusted_score,
        second_best_cell_type = second_best_cell_type,
        second_best_score = second_best_score,
        explanation = explanation,
        stringsAsFactors = FALSE
      )
    }
  }

  suggested_cluster_annotation <- dplyr::bind_rows(annotation_rows)
  write.csv(suggested_cluster_annotation, suggested_annotation_path, row.names = FALSE)

  if (nrow(suggested_cluster_annotation) == 0) {
    object$semi_manual_anno <- "unassigned"
  } else {
    cluster_to_cell_type <- setNames(
      suggested_cluster_annotation$suggested_cell_type,
      suggested_cluster_annotation$cluster
    )
    object$semi_manual_anno <- unname(cluster_to_cell_type[clusters])
    object$semi_manual_anno[is.na(object$semi_manual_anno)] <- "unassigned"
  }

  reduction_name <- sub("^X_", "", embedding_key)
  if (reduction_name %in% names(object@reductions)) {
    umap_plot <- patchwork::wrap_plots(
      Seurat::DimPlot(object, reduction = reduction_name, group.by = cluster_key, label = TRUE),
      Seurat::DimPlot(object, reduction = reduction_name, group.by = "semi_manual_anno"),
      ncol = 2
    )
    png(suggested_umap_path, width = 15, height = 6, units = "in", res = umap_dpi)
    print(umap_plot)
    dev.off()
    ggplot2::ggsave(
      suggested_umap_svg_path,
      plot = umap_plot,
      width = 15,
      height = 6
    )
  }

  marker_gene_sets_plotting <- list()
  metrics_for_kept <- marker_filtration_metrics[marker_filtration_metrics$keep_marker %in% TRUE, c("cell_type", "marker_gene", "specificity_score"), drop = FALSE]
  for (cell_type in names(marker_gene_sets_filtered)) {
    markers <- marker_gene_sets_filtered[[cell_type]]
    markers <- markers[markers %in% rownames(expression_matrix)]
    if (length(markers) <= max_markers_per_cell_type) {
      marker_gene_sets_plotting[[cell_type]] <- markers
    } else {
      ranked <- metrics_for_kept[
        metrics_for_kept$cell_type == cell_type & metrics_for_kept$marker_gene %in% markers,
        ,
        drop = FALSE
      ]
      ranked <- ranked[order(ranked$specificity_score, decreasing = TRUE), , drop = FALSE]
      marker_gene_sets_plotting[[cell_type]] <- head(ranked$marker_gene, max_markers_per_cell_type)
    }
  }
  marker_gene_sets_plotting <- marker_gene_sets_plotting[lengths(marker_gene_sets_plotting) > 0]

  if (length(marker_gene_sets_plotting) > 0) {
    features <- unique(unlist(marker_gene_sets_plotting, use.names = FALSE))
    dotplot <- Seurat::DotPlot(object, features = features, group.by = cluster_key) +
      Seurat::RotatedAxis()
    if (isTRUE(dotplot_swap_axes)) {
      dotplot <- dotplot + ggplot2::coord_flip()
    }
    ggplot2::ggsave(validation_dotplot_path, plot = dotplot, dpi = dotplot_dpi, width = 12, height = 8)
    ggplot2::ggsave(validation_dotplot_svg_path, plot = dotplot, width = 12, height = 8)
  }

  clusters_for_manual_review <- if (nrow(suggested_cluster_annotation) == 0) {
    suggested_cluster_annotation
  } else {
    suggested_cluster_annotation[
      suggested_cluster_annotation$confidence_label == "low" |
        suggested_cluster_annotation$confidence_margin < manual_review_margin_threshold,
      ,
      drop = FALSE
    ]
  }
  write.csv(clusters_for_manual_review, manual_review_path, row.names = FALSE)

  parameters <- list(
    cluster_key = cluster_key,
    marker_gene_sets_json = marker_gene_sets_json,
    annotation_dir = annotation_dir,
    marker_gene_dir = marker_gene_dir,
    max_global_fraction = max_global_fraction,
    min_global_fraction = min_global_fraction,
    min_specificity_score = min_specificity_score,
    min_retained_markers = min_retained_markers,
    high_confidence_min_score = high_confidence_min_score,
    high_confidence_min_margin = high_confidence_min_margin,
    medium_confidence_min_score = medium_confidence_min_score,
    manual_review_margin_threshold = manual_review_margin_threshold,
    embedding_key = embedding_key,
    max_markers_per_cell_type = max_markers_per_cell_type,
    dotplot_standard_scale = dotplot_standard_scale,
    dotplot_swap_axes = dotplot_swap_axes,
    dotplot_dendrogram = dotplot_dendrogram,
    dotplot_dpi = dotplot_dpi,
    umap_dpi = umap_dpi,
    output_prefix = output_prefix,
    assay = assay %||% Seurat::DefaultAssay(object),
    slot = slot,
    object_metadata_annotation_key = "semi_manual_anno"
  )

  message("Semi-automated marker annotation completed.")
  message(sprintf("Original cell types: %s", length(marker_gene_sets)))
  message(sprintf("Retained cell types: %s", length(marker_gene_sets_filtered)))
  message(sprintf("Original marker count: %s", sum(lengths(marker_gene_sets))))
  message(sprintf("Retained marker count: %s", sum(lengths(marker_gene_sets_filtered))))

  list(
    object = object,
    marker_presence_summary = marker_presence_summary,
    marker_filtration_metrics = marker_filtration_metrics,
    marker_gene_sets_filtered = marker_gene_sets_filtered,
    cluster_celltype_marker_scores = cluster_celltype_marker_scores,
    suggested_cluster_annotation = suggested_cluster_annotation,
    clusters_for_manual_review = clusters_for_manual_review,
    marker_gene_sets_plotting = marker_gene_sets_plotting,
    parameters = parameters,
    output_paths = list(
      marker_presence_summary = marker_presence_path,
      marker_filtration_metrics = filtration_metrics_path,
      filtered_marker_gene_sets = filtered_marker_path,
      cluster_celltype_marker_scores = cluster_scores_path,
      suggested_cluster_annotation = suggested_annotation_path,
      clusters_for_manual_review = manual_review_path,
      suggested_annotation_umap_png = suggested_umap_path,
      suggested_annotation_umap_svg = suggested_umap_svg_path,
      filtered_marker_validation_dotplot_png = validation_dotplot_path,
      filtered_marker_validation_dotplot_svg = validation_dotplot_svg_path
    )
  )
}

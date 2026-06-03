from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def _mean_by_column(matrix: Any) -> np.ndarray:
    return np.asarray(matrix.mean(axis=0)).ravel()


def _fraction_positive_by_column(matrix: Any) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray((matrix > 0).mean(axis=0)).ravel()
    return (np.asarray(matrix) > 0).mean(axis=0)


def _cluster_order(adata, cluster_key: str) -> list[str]:
    clusters = adata.obs[cluster_key]
    if isinstance(clusters.dtype, pd.CategoricalDtype):
        return [str(cluster) for cluster in clusters.cat.categories]
    return sorted(
        clusters.astype(str).unique(),
        key=lambda value: int(value) if value.isdigit() else value,
    )


def semi_automated_marker_annotation(
    adata,
    cluster_key: str,
    marker_gene_sets_json: str | Path,
    annotation_dir: str | Path = "results/annotation",
    marker_gene_dir: str | Path = "data/marker_gene",
    *,
    max_global_fraction: float = 0.60,
    min_global_fraction: float = 0.01,
    min_specificity_score: float = 1.5,
    min_retained_markers: int = 2,
    high_confidence_min_score: float = 0.75,
    high_confidence_min_margin: float = 0.20,
    medium_confidence_min_score: float = 0.50,
    manual_review_margin_threshold: float = 0.15,
    embedding_key: str = "X_umap_harmony_pcs30_nn30",
    max_markers_per_cell_type: int = 5,
    dotplot_standard_scale: str | None = "var",
    dotplot_swap_axes: bool = True,
    dotplot_dendrogram: bool = False,
    dotplot_dpi: int = 400,
    umap_dpi: int = 300,
    output_prefix: str = "",
) -> dict[str, Any]:
    """Run semi-automated marker-based cluster annotation for a Scanpy AnnData.

    The marker dictionary is read from JSON, marker filtering and cluster scoring
    metrics are written to disk, and suggested cell types are stored in
    ``adata.obs["semi_manual_anno"]`` for the current in-memory AnnData object.
    """

    marker_gene_sets_json = Path(marker_gene_sets_json)
    annotation_dir = Path(annotation_dir)
    marker_gene_dir = Path(marker_gene_dir)
    annotation_dir.mkdir(parents=True, exist_ok=True)
    marker_gene_dir.mkdir(parents=True, exist_ok=True)

    if cluster_key not in adata.obs:
        raise KeyError(f"{cluster_key!r} is not present in adata.obs")

    with marker_gene_sets_json.open() as handle:
        marker_gene_sets = json.load(handle)

    if not isinstance(marker_gene_sets, dict) or not marker_gene_sets:
        raise ValueError("marker_gene_sets_json must contain a non-empty JSON object")

    marker_gene_sets = {
        str(cell_type): list(dict.fromkeys(markers))
        for cell_type, markers in marker_gene_sets.items()
    }

    prefix = f"{output_prefix}_" if output_prefix else ""
    filtered_marker_path = marker_gene_dir / f"{prefix}filtered_marker_gene_sets_core.json"
    marker_presence_path = annotation_dir / f"{prefix}marker_presence_summary.csv"
    filtration_metrics_path = annotation_dir / f"{prefix}marker_filtration_metrics.csv"
    cluster_scores_path = annotation_dir / f"{prefix}cluster_celltype_marker_scores.csv"
    suggested_annotation_path = annotation_dir / f"{prefix}suggested_cluster_annotation.csv"
    validation_dotplot_path = annotation_dir / f"{prefix}filtered_marker_validation_dotplot.png"
    validation_dotplot_svg_path = annotation_dir / f"{prefix}filtered_marker_validation_dotplot.svg"
    suggested_umap_path = annotation_dir / f"{prefix}suggested_annotation_umap.png"
    suggested_umap_svg_path = annotation_dir / f"{prefix}suggested_annotation_umap.svg"
    manual_review_path = annotation_dir / f"{prefix}clusters_for_manual_review.csv"

    clusters = adata.obs[cluster_key].astype(str)
    cluster_order = _cluster_order(adata, cluster_key)

    marker_gene_sets_present = {
        cell_type: [gene for gene in genes if gene in adata.var_names]
        for cell_type, genes in marker_gene_sets.items()
    }
    marker_genes_missing = {
        cell_type: [gene for gene in genes if gene not in adata.var_names]
        for cell_type, genes in marker_gene_sets.items()
    }
    marker_presence_summary = pd.DataFrame(
        [
            {
                "cell_type": cell_type,
                "total_markers": len(marker_gene_sets[cell_type]),
                "present_markers": len(marker_gene_sets_present[cell_type]),
                "missing_markers": len(marker_genes_missing[cell_type]),
                "missing_marker_genes": ";".join(marker_genes_missing[cell_type]),
            }
            for cell_type in marker_gene_sets
        ]
    )
    marker_presence_summary.to_csv(marker_presence_path, index=False)

    present_markers = list(
        dict.fromkeys(gene for genes in marker_gene_sets_present.values() for gene in genes)
    )
    if not present_markers:
        raise ValueError("No marker genes from marker_gene_sets_json were found in adata.var_names")

    marker_indices = adata.var_names.get_indexer(present_markers)
    X_markers = adata.X[:, marker_indices]

    global_mean_expression = pd.Series(_mean_by_column(X_markers), index=present_markers)
    global_fraction_expressing = pd.Series(
        _fraction_positive_by_column(X_markers),
        index=present_markers,
    )

    cluster_mean_expression = pd.DataFrame(index=cluster_order, columns=present_markers, dtype=float)
    cluster_fraction_expressing = pd.DataFrame(index=cluster_order, columns=present_markers, dtype=float)
    for cluster in cluster_order:
        cluster_mask = (clusters == cluster).to_numpy()
        X_cluster = X_markers[cluster_mask, :]
        cluster_mean_expression.loc[cluster] = _mean_by_column(X_cluster)
        cluster_fraction_expressing.loc[cluster] = _fraction_positive_by_column(X_cluster)

    metric_rows = []
    for cell_type, markers in marker_gene_sets.items():
        for marker in markers:
            if marker not in present_markers:
                metric_rows.append(
                    {
                        "cell_type": cell_type,
                        "marker_gene": marker,
                        "marker_present": False,
                        "global_mean_expression": np.nan,
                        "global_fraction_expressing": np.nan,
                        "max_cluster_mean": np.nan,
                        "second_cluster_mean": np.nan,
                        "specificity_score": np.nan,
                        "cluster_separation": np.nan,
                        "top_cluster": pd.NA,
                        "broadly_expressed": pd.NA,
                        "low_detection": pd.NA,
                        "low_specificity": pd.NA,
                        "keep_marker": False,
                    }
                )
                continue

            cluster_means = cluster_mean_expression[marker].astype(float)
            sorted_means = cluster_means.sort_values(ascending=False)
            max_cluster_mean = float(sorted_means.iloc[0])
            second_cluster_mean = float(sorted_means.iloc[1]) if len(sorted_means) > 1 else 0.0
            specificity_score = max_cluster_mean / (float(global_mean_expression[marker]) + 1e-6)
            cluster_separation = max_cluster_mean - second_cluster_mean
            broadly_expressed = float(global_fraction_expressing[marker]) > max_global_fraction
            low_detection = float(global_fraction_expressing[marker]) < min_global_fraction
            low_specificity = specificity_score < min_specificity_score
            keep_marker = not (broadly_expressed or low_detection or low_specificity)

            row = {
                "cell_type": cell_type,
                "marker_gene": marker,
                "marker_present": True,
                "global_mean_expression": float(global_mean_expression[marker]),
                "global_fraction_expressing": float(global_fraction_expressing[marker]),
                "max_cluster_mean": max_cluster_mean,
                "second_cluster_mean": second_cluster_mean,
                "specificity_score": float(specificity_score),
                "cluster_separation": float(cluster_separation),
                "top_cluster": sorted_means.index[0],
                "broadly_expressed": broadly_expressed,
                "low_detection": low_detection,
                "low_specificity": low_specificity,
                "keep_marker": keep_marker,
            }
            for cluster in cluster_order:
                row[f"cluster_{cluster}_mean_expression"] = float(
                    cluster_mean_expression.loc[cluster, marker]
                )
                row[f"cluster_{cluster}_fraction_expressing"] = float(
                    cluster_fraction_expressing.loc[cluster, marker]
                )
            metric_rows.append(row)

    marker_filtration_metrics = pd.DataFrame(metric_rows)
    marker_filtration_metrics.to_csv(filtration_metrics_path, index=False)

    kept_marker_pairs = set(
        marker_filtration_metrics.loc[
            marker_filtration_metrics["keep_marker"].fillna(False),
            ["cell_type", "marker_gene"],
        ].itertuples(index=False, name=None)
    )
    marker_gene_sets_filtered = {
        cell_type: [
            marker
            for marker in marker_gene_sets_present[cell_type]
            if (cell_type, marker) in kept_marker_pairs
        ]
        for cell_type in marker_gene_sets_present
    }
    marker_gene_sets_filtered = {
        cell_type: markers
        for cell_type, markers in marker_gene_sets_filtered.items()
        if len(markers) >= min_retained_markers
    }

    with filtered_marker_path.open("w") as handle:
        json.dump(marker_gene_sets_filtered, handle, indent=2)

    retained_markers = list(
        dict.fromkeys(gene for genes in marker_gene_sets_filtered.values() for gene in genes)
    )
    score_rows = []
    if retained_markers:
        cluster_average_expression = cluster_mean_expression.loc[cluster_order, retained_markers].astype(float)
        gene_means = cluster_average_expression.mean(axis=0)
        gene_stds = cluster_average_expression.std(axis=0, ddof=0).replace(0, 1.0)
        marker_zscores = (cluster_average_expression - gene_means) / gene_stds

        for cluster in cluster_order:
            raw_cluster_rows = []
            for cell_type, markers in marker_gene_sets_filtered.items():
                marker_scores = marker_zscores.loc[cluster, markers].astype(float)
                mean_scaled_marker_expression = float(marker_scores.mean())
                fraction_positive_markers = float((marker_scores > 0).mean())
                specificity_adjusted_score = mean_scaled_marker_expression * fraction_positive_markers
                raw_cluster_rows.append(
                    {
                        "cluster": cluster,
                        "cell_type": cell_type,
                        "n_markers_used": len(markers),
                        "mean_scaled_marker_expression": mean_scaled_marker_expression,
                        "fraction_positive_markers": fraction_positive_markers,
                        "specificity_adjusted_score": float(specificity_adjusted_score),
                        "markers_used": ";".join(markers),
                    }
                )

            raw_scores = np.array([row["specificity_adjusted_score"] for row in raw_cluster_rows], dtype=float)
            score_min = float(raw_scores.min()) if len(raw_scores) else 0.0
            score_max = float(raw_scores.max()) if len(raw_scores) else 0.0
            score_range = score_max - score_min
            for row in raw_cluster_rows:
                if score_range > 0:
                    row["possibility_score"] = (
                        row["specificity_adjusted_score"] - score_min
                    ) / score_range
                else:
                    row["possibility_score"] = 1.0 if score_max > 0 else 0.0
                score_rows.append(row)

    cluster_celltype_marker_scores = pd.DataFrame(score_rows)
    cluster_celltype_marker_scores.to_csv(cluster_scores_path, index=False)

    annotation_rows = []
    if not cluster_celltype_marker_scores.empty:
        for cluster in cluster_order:
            cluster_scores = cluster_celltype_marker_scores.loc[
                cluster_celltype_marker_scores["cluster"] == cluster
            ].sort_values(["possibility_score", "specificity_adjusted_score"], ascending=False)
            if cluster_scores.empty:
                continue

            best = cluster_scores.iloc[0]
            second = cluster_scores.iloc[1] if len(cluster_scores) > 1 else None
            second_best_cell_type = second["cell_type"] if second is not None else pd.NA
            second_best_score = float(second["possibility_score"]) if second is not None else 0.0
            confidence_margin = float(best["possibility_score"] - second_best_score)

            if best["possibility_score"] >= high_confidence_min_score and confidence_margin >= high_confidence_min_margin:
                confidence_label = "high"
            elif best["possibility_score"] >= medium_confidence_min_score:
                confidence_label = "medium"
            else:
                confidence_label = "low"

            explanation = (
                f"Cluster {cluster} is suggested as {best['cell_type']} based on retained markers "
                f"{best['markers_used']} with marker score {best['specificity_adjusted_score']:.3f} "
                f"and possibility score {best['possibility_score']:.3f}. "
                f"The competing annotation is {second_best_cell_type} with score {second_best_score:.3f}; "
                f"confidence margin is {confidence_margin:.3f}."
            )

            annotation_rows.append(
                {
                    "cluster": cluster,
                    "suggested_cell_type": best["cell_type"],
                    "possibility_score": float(best["possibility_score"]),
                    "confidence_label": confidence_label,
                    "confidence_margin": confidence_margin,
                    "markers_used": best["markers_used"],
                    "n_markers_used": int(best["n_markers_used"]),
                    "mean_scaled_marker_expression": float(best["mean_scaled_marker_expression"]),
                    "fraction_positive_markers": float(best["fraction_positive_markers"]),
                    "specificity_adjusted_score": float(best["specificity_adjusted_score"]),
                    "second_best_cell_type": second_best_cell_type,
                    "second_best_score": second_best_score,
                    "explanation": explanation,
                }
            )

    suggested_cluster_annotation = pd.DataFrame(annotation_rows)
    suggested_cluster_annotation.to_csv(suggested_annotation_path, index=False)

    if suggested_cluster_annotation.empty:
        adata.obs["semi_manual_anno"] = "unassigned"
    else:
        cluster_to_cell_type = suggested_cluster_annotation.set_index("cluster")["suggested_cell_type"].to_dict()
        cluster_to_cell_type = {str(cluster): cell_type for cluster, cell_type in cluster_to_cell_type.items()}
        adata.obs["semi_manual_anno"] = clusters.map(cluster_to_cell_type).fillna("unassigned")

    if embedding_key in adata.obsm:
        embedding_basis = embedding_key.removeprefix("X_")
        fig, axes = plt.subplots(1, 2, figsize=(15, 6), constrained_layout=True)
        sc.pl.embedding(
            adata,
            basis=embedding_basis,
            color=cluster_key,
            legend_loc="on data",
            legend_fontoutline=2,
            title=cluster_key,
            ax=axes[0],
            show=False,
        )
        sc.pl.embedding(
            adata,
            basis=embedding_basis,
            color="semi_manual_anno",
            legend_loc="right margin",
            title="semi_manual_anno",
            ax=axes[1],
            show=False,
        )
        fig.savefig(suggested_umap_path, dpi=umap_dpi, bbox_inches="tight")
        fig.savefig(suggested_umap_svg_path, bbox_inches="tight")
        plt.close(fig)

    marker_gene_sets_plotting: dict[str, list[str]] = {}
    metrics_for_kept = marker_filtration_metrics.loc[
        marker_filtration_metrics["keep_marker"].fillna(False),
        ["cell_type", "marker_gene", "specificity_score"],
    ]
    for cell_type, markers in marker_gene_sets_filtered.items():
        markers = [marker for marker in markers if marker in adata.var_names]
        if len(markers) <= max_markers_per_cell_type:
            marker_gene_sets_plotting[cell_type] = markers
        else:
            ranked = metrics_for_kept.loc[
                (metrics_for_kept["cell_type"] == cell_type)
                & (metrics_for_kept["marker_gene"].isin(markers))
            ].sort_values("specificity_score", ascending=False)
            marker_gene_sets_plotting[cell_type] = ranked["marker_gene"].head(max_markers_per_cell_type).tolist()

    marker_gene_sets_plotting = {
        cell_type: markers
        for cell_type, markers in marker_gene_sets_plotting.items()
        if markers
    }
    if marker_gene_sets_plotting:
        dotplot = sc.pl.dotplot(
            adata,
            marker_gene_sets_plotting,
            groupby=cluster_key,
            standard_scale=dotplot_standard_scale,
            swap_axes=dotplot_swap_axes,
            dendrogram=dotplot_dendrogram,
            show=False,
            return_fig=True,
        )
        dotplot.savefig(validation_dotplot_path, dpi=dotplot_dpi, bbox_inches="tight")
        dotplot.savefig(validation_dotplot_svg_path, bbox_inches="tight")
        plt.close("all")

    if suggested_cluster_annotation.empty:
        clusters_for_manual_review = suggested_cluster_annotation.copy()
    else:
        clusters_for_manual_review = suggested_cluster_annotation.loc[
            (suggested_cluster_annotation["confidence_label"] == "low")
            | (suggested_cluster_annotation["confidence_margin"] < manual_review_margin_threshold)
        ].copy()
    clusters_for_manual_review.to_csv(manual_review_path, index=False)

    parameters = {
        "cluster_key": cluster_key,
        "marker_gene_sets_json": str(marker_gene_sets_json),
        "annotation_dir": str(annotation_dir),
        "marker_gene_dir": str(marker_gene_dir),
        "max_global_fraction": max_global_fraction,
        "min_global_fraction": min_global_fraction,
        "min_specificity_score": min_specificity_score,
        "min_retained_markers": min_retained_markers,
        "high_confidence_min_score": high_confidence_min_score,
        "high_confidence_min_margin": high_confidence_min_margin,
        "medium_confidence_min_score": medium_confidence_min_score,
        "manual_review_margin_threshold": manual_review_margin_threshold,
        "embedding_key": embedding_key,
        "max_markers_per_cell_type": max_markers_per_cell_type,
        "dotplot_standard_scale": dotplot_standard_scale,
        "dotplot_swap_axes": dotplot_swap_axes,
        "dotplot_dendrogram": dotplot_dendrogram,
        "dotplot_dpi": dotplot_dpi,
        "umap_dpi": umap_dpi,
        "output_prefix": output_prefix,
        "adata_obs_annotation_key": "semi_manual_anno",
    }

    print("Semi-automated marker annotation completed.")
    print(f"Original cell types: {len(marker_gene_sets)}")
    print(f"Retained cell types: {len(marker_gene_sets_filtered)}")
    print(f"Original marker count: {sum(len(markers) for markers in marker_gene_sets.values())}")
    print(f"Retained marker count: {sum(len(markers) for markers in marker_gene_sets_filtered.values())}")
    print("Output files:")
    for path in [
        marker_presence_path,
        filtration_metrics_path,
        filtered_marker_path,
        cluster_scores_path,
        suggested_annotation_path,
        manual_review_path,
        suggested_umap_path,
        suggested_umap_svg_path,
        validation_dotplot_path,
        validation_dotplot_svg_path,
    ]:
        print(f"  {path}")
    print("Parameters:")
    for key, value in parameters.items():
        print(f"  {key}: {value}")

    return {
        "marker_presence_summary": marker_presence_summary,
        "marker_filtration_metrics": marker_filtration_metrics,
        "marker_gene_sets_filtered": marker_gene_sets_filtered,
        "cluster_celltype_marker_scores": cluster_celltype_marker_scores,
        "suggested_cluster_annotation": suggested_cluster_annotation,
        "clusters_for_manual_review": clusters_for_manual_review,
        "marker_gene_sets_plotting": marker_gene_sets_plotting,
        "parameters": parameters,
        "output_paths": {
            "marker_presence_summary": marker_presence_path,
            "marker_filtration_metrics": filtration_metrics_path,
            "filtered_marker_gene_sets": filtered_marker_path,
            "cluster_celltype_marker_scores": cluster_scores_path,
            "suggested_cluster_annotation": suggested_annotation_path,
            "clusters_for_manual_review": manual_review_path,
            "suggested_annotation_umap_png": suggested_umap_path,
            "suggested_annotation_umap_svg": suggested_umap_svg_path,
            "filtered_marker_validation_dotplot_png": validation_dotplot_path,
            "filtered_marker_validation_dotplot_svg": validation_dotplot_svg_path,
        },
    }

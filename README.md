# scRNAseq Pixi Environment

This repository provides a reproducible `pixi` environment for scRNA-seq analysis and a small collection of reusable analysis functions. It is intended to make environment setup easier on HPC systems, especially when some R/Bioconductor packages are difficult to install reliably with `pixi install` or `BiocManager::install()`.

The repository currently includes a semi-automated cell type annotation function for both Scanpy and Seurat workflows. The function labels clusters from user-provided marker gene sets and writes marker filtering, scoring, annotation, and review outputs.

## 1. Install Pixi

Install `pixi` in a terminal. On Linux or macOS, the official installer is:

```bash
curl -fsSL https://pixi.sh/install.sh | sh
```

If `curl` is not available:

```bash
wget -qO- https://pixi.sh/install.sh | sh
```

Restart the shell or add Pixi to `PATH`:

```bash
export PATH="$HOME/.pixi/bin:$PATH"
pixi --version
```

Pixi installation reference: https://pixi.prefix.dev/latest/installation/

## 2. Configure Pixi On HPC

On HPC systems, avoid placing large caches in a small or quota-limited home directory. Put Pixi environments and package caches on a project, work, or scratch filesystem.

Example:

```bash
mkdir -p /path/to/project/pixi-home
mkdir -p /path/to/project/pixi-cache

export PIXI_HOME=/path/to/project/pixi-home
export PIXI_CACHE_DIR=/path/to/project/pixi-cache
export PATH="$HOME/.pixi/bin:$PATH"
```

To make this persistent, add those lines to `~/.bashrc` or your HPC shell startup file, adjusting the paths for your allocation.

This setup is expected to work from a normal terminal or an Open OnDemand VS Code terminal. Open OnDemand RStudio can have trouble inheriting the same shell environment, so setting up and validating this Pixi environment from RStudio is not currently recommended.

If you want to code in R with Pixi from an LRZ Open OnDemand VS Code session, follow exaclty [OnDemand_VScode.md](OnDemand_VScode.md). It uses `part2` as the example environment for running R notebooks through Jupyter and VS Code.

To build up a reproducible workflow as documentation for the project, follow GOOD_CONDUCT.md. It is very important to have a clear structure for publication.

## 3. Clone This Repository

Clone the repository and enter it:

```bash
git clone https://github.com/Yuling-Zhao/scRNAseq_pixi.git
cd scRNAseq_pixi
```

The repository contains:

- `pixi.toml`: environment definitions
- `pixi.lock`: locked dependency versions
- `bio_pkg/GenomeInfoDbData_1.2.13.tar.gz`: local Bioconductor data package
- `bio_pkg/celldex_1.16.0.tar.gz`: local Bioconductor annotation reference package
- `install_bioc_packages.sh`: installs the bundled packages into the Pixi R environment
- `validate_environment.sh`: checks that the environments and required tools load correctly

## 4. Install All Pixi Environments

Install all environments from the lockfile:

```bash
pixi install --all
```

The environments are split by analysis stage:

- `part1`: command-line tools, including SRA Toolkit, FastQC, MultiQC, and SAMtools
- `part2`: R and core QC packages, including Seurat, scater, SingleCellExperiment, scDblFinder, DropletUtils, SoupX, and supporting R packages
- `part3`: integration and visualization packages, including harmony, batchelor, ggplot2, patchwork, ggalluvial, viridis, RColorBrewer, and ggrepel
- `part4`: annotation packages, including SingleR, HGNChelper, openxlsx, and celldex

## 5. Install Bundled Bioconductor Packages

Two packages are bundled because they can be difficult to install reproducibly through `pixi install` or `BiocManager::install()` in some HPC environments:

- `GenomeInfoDbData_1.2.13.tar.gz`
- `celldex_1.16.0.tar.gz`

Install them with:

```bash
bash install_bioc_packages.sh
```

The script installs both packages into the `part4` Pixi environment and then tests:

```r
library(GenomeInfoDbData)
library(celldex)
```

## 6. Validate The Environment

Run the validation script:

```bash
bash validate_environment.sh
```

The script checks:

- `part1`: SRA Toolkit, FastQC, MultiQC, SAMtools
- `part2`: R, Seurat, DropletUtils, scater, SingleCellExperiment, scDblFinder
- `part3`: harmony, batchelor, ggplot2, patchwork, ggalluvial
- `part4`: SingleR, celldex, HGNChelper, openxlsx

If all checks pass, the environment is ready for scRNA-seq analysis.

## 7. Use The Callable Functions

### Python / Scanpy

Import the Scanpy implementation:

```python
from py_functions import semi_automated_marker_annotation

result = semi_automated_marker_annotation(
    adata,
    cluster_key="leiden",
    marker_gene_sets_json="marker_gene_sets.json",
    annotation_dir="results/annotation",
    marker_gene_dir="data/marker_gene",
)
```

The function updates the in-memory `AnnData` object with:

```python
adata.obs["semi_manual_anno"]
```

It also returns tables and output paths in `result`.

### R / Seurat

Source the R functions:

```r
source("r_functions/load_functions.R")

result <- semi_automated_marker_annotation(
  object = seurat_obj,
  cluster_key = "seurat_clusters",
  marker_gene_sets_json = "marker_gene_sets.json",
  annotation_dir = "results/annotation",
  marker_gene_dir = "data/marker_gene"
)
```

The returned Seurat object is available as:

```r
seurat_obj <- result$object
```

The annotation is stored in:

```r
seurat_obj$semi_manual_anno
```

## Marker Gene JSON Format

The marker gene file should be a JSON object where each key is a cell type and each value is a list of marker genes:

```json
{
  "T cell": ["CD3D", "CD3E", "TRAC"],
  "B cell": ["MS4A1", "CD79A", "CD79B"],
  "Myeloid": ["LYZ", "S100A8", "S100A9"]
}
```

## Notes

- Run commands from the repository root unless noted otherwise.
- Prefer terminal or Open OnDemand VS Code terminal workflows.
- Open OnDemand RStudio may not inherit the expected Pixi shell configuration, so package discovery can fail there even when terminal validation passes.
- The helper functions are designed for repeatable work and can be extended with additional scRNA-seq analysis functions over time.

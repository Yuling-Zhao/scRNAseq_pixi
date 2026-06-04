# Protocol: Reproducible Omics Analysis with Pixi, VS Code, Notebooks, and Quarto

## Philosophy

The workflow consists of two phases:

### Phase 1: Interactive Exploration

Goal:

* Understand the data
* Test parameters
* Visualize results
* Make analysis decisions

Tools:

* VS Code (OnDemand)
* Jupyter Notebooks (.ipynb)
* Pixi environments
* Git

Outputs:

* Notes
* Figures
* Parameter choices
* Prototype code

### Phase 2: Reproducible Analysis

Goal:

* Convert exploratory work into a reproducible pipeline
* Generate figures and tables automatically
* Enable rerunning the entire analysis on new data

Tools:

* Quarto (.qmd)
* R scripts / Python scripts
* Pixi environments
* Git

Outputs:

* HTML/PDF reports
* Publication figures
* Processed data objects
* Reproducible workflow

---

# Project Structure

```text
project/
├── pixi.toml
├── notebooks/
│   ├── 01_qc.ipynb
│   ├── 02_integration.ipynb
│   └── 03_annotation.ipynb
│
├── scripts/
│   ├── 01_qc.py
│   ├── 02_integration.py
│   └── 03_annotation.R
│
├── reports/
│   ├── qc_report.qmd
│   ├── integration_report.qmd
│   └── annotation_report.qmd
│
├── data/
│   ├── raw/
│   └── processed/
│
├── results/
│   ├── figures/
│   ├── tables/
│   └── objects/
│
└── README.md
```

---

# Phase 1: Interactive Exploration

## Step 1. Create Pixi Environments

Example:

```toml
[feature.python.dependencies]
python = "3.11.*"
scanpy = "*"
scvi-tools = "*"
jupyterlab = "*"
ipykernel = "*"

[feature.r.dependencies]
r-base = "*"
r-irkernel = "*"
r-seurat = "*"
bioconductor-singler = "*"
bioconductor-celldex = "*"

[environments]
python = ["python"]
r = ["r"]
```

Install:

```bash
pixi install
```

---

## Step 2. Register Notebook Kernels (Once Only)

Python:

```bash
pixi run -e python python -m ipykernel install \
  --user \
  --name python \
  --display-name "Python (pixi)"
```

R:

```bash
pixi run -e r R
```

```r
IRkernel::installspec(
    name="r",
    displayname="R (pixi)",
    user=TRUE
)
```

This only needs to be done once.

---

## Step 3. Start Jupyter From Pixi

Start Jupyter Lab from the same VS Code terminal:

```bash
pixi run -e r jupyter lab --no-browser --ip=127.0.0.1 --port=8888
```

Copy the local URL printed by Jupyter. It will look similar to:

```text
http://127.0.0.1:8888/lab?token=...
```

Keep this terminal running while you work.

Connect VS Code To Jupyter:

Open an `.ipynb` notebook in VS Code.

Select the kernel picker, choose `Existing Jupyter Server`, and paste the Jupyter URL from the previous step. Then select:

```text
R (pixi)
```

Notebook cells now run R code from the `R (pixi)` Pixi environment.

---

## Step 4. Explore in Notebooks

Examples:

### QC

```python
sc.pp.calculate_qc_metrics(...)
sc.pl.violin(...)
```

### Integration

```python
sc.external.pp.harmony_integrate(...)
```

### Annotation

```python
sc.tl.rank_genes_groups(...)
```

### Visualization

```python
sc.pl.umap(...)
```

Plots appear directly beneath notebook cells.

---

## Step 5. Record Decisions

Every notebook should contain markdown sections documenting:

### Example

```text
Decision:
Use 30 PCs.

Reason:
PC30 and PC50 produced similar biological structure.
PC50 showed slightly increased sample-driven separation.
```

Record:

* QC thresholds
* Number of PCs
* Integration method
* Clustering resolution
* Annotation rationale

---

## Step 6. Save Intermediate Objects

Never rely on notebook memory.

Example:

```python
adata.write(
    "data/processed/20260604_integrated.h5ad"
)
```

or

```r
saveRDS(
    seurat_obj,
    "data/processed/seurat_integrated.rds"
)
```

---

# Phase 2: Convert Exploration into Reproducible Analysis

## Step 7. Move Stable Code into Scripts

Notebook:

```python
sc.pp.normalize_total(...)
sc.pp.log1p(...)
sc.tl.pca(...)
```

becomes:

```python
# scripts/02_integration.py
```

All parameter choices are fixed.

No manual intervention required.

---

## Step 8. Create Quarto Reports

Example:

```yaml
---
title: "QC Report"
format: html
---
```

```python
#| label: qc-plots

sc.pl.violin(...)
```

or

```r
#| label: qc-plots

VlnPlot(...)
```

The report combines:

* text
* code
* tables
* figures

into one document.

---

## Step 9. Execute Reports from the Terminal

Render a report:

```bash
pixi run quarto render reports/qc_report.qmd
```

Output:

```text
reports/qc_report.html
```

---

## Step 10. Run the Entire Analysis

Create a master script:

```bash
#!/bin/bash

pixi run python scripts/01_qc.py
pixi run python scripts/02_integration.py
pixi run Rscript scripts/03_annotation.R

pixi run quarto render reports/qc_report.qmd
pixi run quarto render reports/integration_report.qmd
pixi run quarto render reports/annotation_report.qmd
```

Run:

```bash
bash run_analysis.sh
```

Everything is regenerated automatically.

---

# Git Workflow

During exploration:

```bash
git add notebooks/
git commit -m "Test harmony integration"
```

After finalizing:

```bash
git add scripts/
git add reports/
git commit -m "Finalize integration workflow"
```

Avoid committing:

```text
results/
data/raw/
```

unless intentionally versioned.

---

# Recommended Workflow for scRNA-seq / scATAC-seq Projects

1. Explore interactively in notebooks.
2. Save figures and intermediate objects.
3. Document decisions in notebook markdown.
4. Move finalized code into scripts.
5. Build Quarto reports from finalized scripts.
6. Render reports through Pixi.
7. Version-control everything with Git.
8. Re-run the complete analysis from a single command.

The notebook is the laboratory notebook.

The script is the reproducible method.

The Quarto report is the final documented result.

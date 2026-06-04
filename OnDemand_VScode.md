# 

# Open OnDemand VS Code With Pixi R Environments

This guide shows how to use the `part2` Pixi environment as an R notebook kernel from an LRZ Open OnDemand VS Code session. The same pattern can be reused for other Pixi environments that include R and Jupyter.

## 1. Start VS Code On LRZ

Open an interactive VS Code session from LRZ Open OnDemand. In the VS Code terminal, move to this repository:

```bash
cd /path/to/scRNAseq_pixi
```

If your Pixi cache should live outside your home directory, set the cache paths before installing:

```bash
mkdir -p /path/to/project/pixi-home
mkdir -p /path/to/project/pixi-cache

export PIXI_HOME=/path/to/project/pixi-home
export PIXI_CACHE_DIR=/path/to/project/pixi-cache
export PATH="$HOME/.pixi/bin:$PATH"
```

Use paths from your own LRZ project, work, or scratch allocation.

## 2. Install The Example Environment

Install `part2`, which contains R, Seurat, Bioconductor QC packages, Jupyter, and `IRkernel`:

```bash
pixi install -e part2
```

You can check that R starts inside the environment:

```bash
pixi run -e part2 R --version
```

## 3. Register The R Kernel

Register an R kernel named `part2` once per Pixi installation:

```bash
pixi run -e part2 R -e 'IRkernel::installspec(name = "part2", displayname = "R: pixi part2", user = TRUE)'
```

Verify that Jupyter can see it:

```bash
pixi run -e part2 jupyter kernelspec list
```

The output should include `part2`.

## 4. Start Jupyter From Pixi

Start Jupyter Lab from the same VS Code terminal:

```bash
pixi run -e part2 jupyter lab --no-browser --ip=127.0.0.1 --port=8888
```

Copy the local URL printed by Jupyter. It will look similar to:

```text
http://127.0.0.1:8888/lab?token=...
```

Keep this terminal running while you work.

## 5. Connect VS Code To Jupyter

Open an `.ipynb` notebook in VS Code.

Select the kernel picker, choose `Existing Jupyter Server`, and paste the Jupyter URL from the previous step. Then select:

```text
R: pixi part2
```

Notebook cells now run R code from the `part2` Pixi environment.

## 6. Work Reproducibly

Use notebooks for interactive decisions such as QC thresholds, number of PCs, clustering resolution, or annotation choices. Save important intermediate objects to files instead of relying on notebook memory:

```r
saveRDS(seurat_obj, "data/processed/part2_qc.rds")
```

When the workflow is stable, move repeated code into scripts such as:

```text
scripts/01_qc.R
scripts/02_integration.R
scripts/03_annotation.R
```

Run scripts through Pixi so they use the same environment:

```bash
pixi run -e part2 Rscript scripts/01_qc.R
```

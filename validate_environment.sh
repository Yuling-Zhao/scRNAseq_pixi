#!/bin/bash
 
echo "=========================================="
echo "  scRNA-seq Environment Validation"
echo "  Testing all tools from Parts 1-4"
echo "=========================================="
 
#---------------------------------------
# Part 1: Command-Line Tools
#---------------------------------------
echo -e "\n=== Part 1: Command-Line Tools ===\n"
 
echo "[1/4] Testing SRA Toolkit..."
if pixi run -e part1 fastq-dump --version > /tmp/test_sra.out 2>&1; then
    echo "  ✓ SRA Toolkit working"
else
    echo "  ✗ SRA Toolkit FAILED"
    cat /tmp/test_sra.out
    exit 1
fi
 
echo "[2/4] Testing FastQC..."
if pixi run -e part1 fastqc --version > /tmp/test_fastqc.out 2>&1; then
    echo "  ✓ FastQC working"
else
    echo "  ✗ FastQC FAILED"
    cat /tmp/test_fastqc.out
    exit 1
fi
 
echo "[3/4] Testing MultiQC..."
if pixi run -e part1 multiqc --version > /tmp/test_multiqc.out 2>&1; then
    echo "  ✓ MultiQC working"
else
    echo "  ✗ MultiQC FAILED"
    cat /tmp/test_multiqc.out
    exit 1
fi
 
echo "[4/4] Testing SAMtools..."
if pixi run -e part1 samtools --version > /tmp/test_samtools.out 2>&1; then
    echo "  ✓ SAMtools working"
else
    echo "  ✗ SAMtools FAILED"
    cat /tmp/test_samtools.out
    exit 1
fi
 
#---------------------------------------
# Part 2: R QC Packages
#---------------------------------------
echo -e "\n=== Part 2: Quality Control Packages ===\n"
 
echo "[1/3] Testing R installation..."
if pixi run -e part2 Rscript -e "cat('R working\n')" > /tmp/test_r.out 2>&1; then
    echo "  ✓ R installation working"
else
    echo "  ✗ R FAILED"
    cat /tmp/test_r.out
    exit 1
fi
 
echo "[2/3] Testing Seurat..."
if pixi run -e part2 Rscript -e "suppressMessages(library(Seurat)); cat('Seurat OK\n')" > /tmp/test_seurat.out 2>&1; then
    echo "  ✓ Seurat loaded successfully"
else
    echo "  ✗ Seurat FAILED"
    cat /tmp/test_seurat.out
    exit 1
fi
 
echo "[3/3] Testing Bioconductor packages..."
cat > /tmp/test_bioc.R << 'RCODE'
suppressMessages({
  library(DropletUtils)
  library(scater)
  library(SingleCellExperiment)
  library(scDblFinder)
})
cat('Bioconductor OK\n')
RCODE
 
if pixi run -e part2 Rscript /tmp/test_bioc.R > /tmp/test_bioc.out 2>&1; then
    echo "  ✓ Bioconductor QC packages loaded successfully"
else
    echo "  ✗ Bioconductor packages FAILED"
    cat /tmp/test_bioc.out
    exit 1
fi
 
#---------------------------------------
# Part 3: Integration Packages
#---------------------------------------
echo -e "\n=== Part 3: Integration & Visualization ===\n"
 
echo "[1/2] Testing integration packages..."
cat > /tmp/test_integration.R << 'RCODE'
suppressMessages({
  library(harmony)
  library(batchelor)
})
cat('Integration OK\n')
RCODE
 
if pixi run -e part3 Rscript /tmp/test_integration.R > /tmp/test_int.out 2>&1; then
    echo "  ✓ Integration packages loaded successfully"
else
    echo "  ✗ Integration packages FAILED"
    cat /tmp/test_int.out
    exit 1
fi
 
echo "[2/2] Testing visualization packages..."
cat > /tmp/test_viz.R << 'RCODE'
suppressMessages({
  library(ggplot2)
  library(patchwork)
  library(ggalluvial)
})
cat('Visualization OK\n')
RCODE
 
if pixi run -e part3 Rscript /tmp/test_viz.R > /tmp/test_viz.out 2>&1; then
    echo "  ✓ Visualization packages loaded successfully"
else
    echo "  ✗ Visualization packages FAILED"
    cat /tmp/test_viz.out
    exit 1
fi
 
#---------------------------------------
# Part 4: Annotation Packages
#---------------------------------------
echo -e "\n=== Part 4: Cell Type Annotation ===\n"
 
echo "[1/1] Testing annotation packages..."
cat > /tmp/test_annotation.R << 'RCODE'
suppressMessages({
  library(SingleR)
  library(celldex)
  library(HGNChelper)
  library(openxlsx)
})
cat('Annotation OK\n')
RCODE
 
if pixi run -e part4 Rscript /tmp/test_annotation.R > /tmp/test_annot.out 2>&1; then
    echo "  ✓ Annotation packages loaded successfully"
else
    echo "  ✗ Annotation packages FAILED"
    cat /tmp/test_annot.out
    exit 1
fi
 
#---------------------------------------
# Summary
#---------------------------------------
echo -e "\n=========================================="
echo "  ✓✓✓ All validations passed!"
echo "=========================================="
echo ""
echo "Your environment is ready for scRNA-seq analysis!"
echo ""
echo "Next steps:"
echo "  • Part 1: pixi run -e part1 <command>"
echo "  • Part 2: pixi run -e part2 R"
echo "  • Part 3: pixi run -e part3 R"
echo "  • Part 4: pixi run -e part4 R"
echo ""
echo "Follow the NGS101 tutorial series:"
echo "  https://ngs101.com"
echo ""
 
# Cleanup
rm -f /tmp/test_*.out /tmp/test_*.R

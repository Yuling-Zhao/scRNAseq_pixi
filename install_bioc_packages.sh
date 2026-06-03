#!/bin/bash

set -e

PKG_DIR="/dss/dssfs05/lwp-dss-0003/pn46mo/pn46mo-dss-0000/github/scRNAseq_pixi/bio_pkg"

export LC_ALL=C.UTF-8
export LANG=C.UTF-8
export LC_CTYPE=C.UTF-8

echo "Installing GenomeInfoDbData..."
pixi run -e part4 env \
    R_LIBS_USER="" \
    R_LIBS_SITE="" \
    R CMD INSTALL "${PKG_DIR}/GenomeInfoDbData_1.2.13.tar.gz"

echo "Installing celldex..."
pixi run -e part4 env \
    R_LIBS_USER="" \
    R_LIBS_SITE="" \
    R CMD INSTALL "${PKG_DIR}/celldex_1.16.0.tar.gz"

echo "Testing installation..."

pixi run -e part4 env \
    R_LIBS_USER="" \
    R_LIBS_SITE="" \
    Rscript -e '
        library(GenomeInfoDbData)
        library(celldex)
        cat("GenomeInfoDbData OK\n")
        cat("celldex OK\n")
    '

echo "Done."

#!/bin/bash

set -e

pixi run -e part4 \
    jupyter lab \
    --no-browser \
    --ip=127.0.0.1 \
    --port=8888

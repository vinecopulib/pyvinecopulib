#!/bin/bash
# Create and resize the pyvinecopulib logo

Rscript --vanilla scripts/mklogo.R docs/_static/pyvinecopulib.png
convert -resize 100 docs/_static/pyvinecopulib.png docs/_static/pyvinecopulib.png

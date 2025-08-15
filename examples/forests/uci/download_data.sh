#!/bin/bash

# Download datasets for density estimation
# Papamakarios, George. (2018). Preprocessed datasets for MAF experiments [Data set]. Zenodo. http://doi.org/10.5281/zenodo.1161203

wget https://zenodo.org/record/1161203/files/data.tar.gz
tar -zxvf data.tar.gz
rm -rf data/mnist/
rm -rf data/cifar10/
rm -rf data/BSDS300/
rm data.tar.gz

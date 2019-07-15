#!/usr/bin/env bash

# Script to quickly restore project folder for quick testing.

git add rename.sh
git add restore.sh
git clean -fd

git checkout .

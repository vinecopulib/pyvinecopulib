#!/usr/bin/env bash
# Build the JMLR MLOSS manuscript into pyvinecopulib.pdf.
#
# JMLR's style file is fetched rather than committed, because its repository
# (JmlrOrg/jmlr-style-file) declares no license. The commit is pinned so that
# every build typesets the same way.
set -euo pipefail
cd "$(dirname "$0")"

sty_commit="f413f638b407af76074813f8f88a82a7a5a81e9d"
sty_url="https://raw.githubusercontent.com/JmlrOrg/jmlr-style-file/${sty_commit}/jmlr2e.sty"

if [ ! -f jmlr2e.sty ]; then
  curl -fsSL -o jmlr2e.sty "$sty_url"
fi
latexmk -pdf -interaction=nonstopmode -halt-on-error pyvinecopulib.tex

#!/bin/bash
# Extracts the doxygen documentation from the vinecopulib sources

python3 scripts/mkdocs.py -I'lib/vinecopulib/include' -I'lib/eigen/' $(find lib/vinecopulib/include -regextype awk -regex ".*(class|controls|family|structure|stats)\.(hpp|ipp)" -print) -output='src/docstr.hpp'

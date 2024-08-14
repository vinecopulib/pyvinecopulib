#!/bin/bash
# Extracts the doxygen documentation from the vinecopulib sources

python3 scripts/mkdocs.py -I'lib/vinecopulib/include' -I'/usr/include' $(find lib/vinecopulib/include -regextype awk -regex ".*(class|controls|family|structure|stats)\.(hpp|ipp)" -print) -output='src/include/docstr.hpp' -library_file='/home/tvatter/mambaforge/envs/pykde1d311/lib/libclang.so'
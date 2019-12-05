python3 scripts/mkdocs.py -I'lib/vinecopulib/include' `find lib/vinecopulib/include -regextype awk -regex ".*(class|controls|family|structure|stats)\.(hpp|ipp)" -print` -output='src/docstr.hpp'

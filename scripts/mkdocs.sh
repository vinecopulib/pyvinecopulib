python3 scripts/mkdocs.py -I'lib/vinecopulib/include' `find lib/vinecopulib/include -regextype awk -regex ".*(class|structure)\.(hpp|ipp)" -print` -output='src/docstr.hpp'

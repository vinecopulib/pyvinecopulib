#!/bin/bash
##
## create_boost.sh
##
## Extract and repackage boost.
##
## Derived from CreateBoost.sh from the BH R package by
## Jay Emerson and Dirk Eddelbuettel,  2012 - 2016
##
## Modified by Thibault Vatter,  2019.


## (1) Adjust these variables as needed
## -- repo where the boost will be packaged
pkgdir="${HOME}/Dropbox/github/pyvinecopulib"
## -- repo where the boost sources are currently 
srcdir="${HOME}/Downloads"
## -- current boost sources, placed eg in ${pkgdir}
boosttargz="boost_1_71_0.tar.gz"
## -- current package version and date (and other metadata as needed)
date="2019-08-30"

## (2) Additional resources we require and need to test for
## 'progs' lists the programs we need
progs="bcp"

## (3) Some internal constants and variables
## local libraries directory in git repo
libdir="${pkgdir}/lib"
## same for boost 
boostdir="${libdir}/boost"
## Derive the 'bootroot' name from the tarball, using basename(1)
boostver=$(basename ${boosttargz} ".tar.gz")
## create boost root directory name
boostroot="${libdir}/${boostver}"
## create boost tarball file name with full path
boostsources="${srcdir}/${boosttargz}"

## (4) Display current settings
echo "Using these settings:
Date:          ${date}
Version:       ${boostver}
PkgDir:        ${pkgdir}
BoostDir:      ${boostdir}
Boostsources:  ${boostsources}
"
# exit -1

## (5) Some sanity check here before continuing
for prog in ${progs}; do
    if [ ! -x /usr/bin/${prog} ] && [ ! -x /usr/local/bin/${prog} ]; then
	echo "Program '${prog}' not found, exiting"
	exit 1
    fi
done

if [ ! -f "${boostsources}" ]; then
    echo "Boost input file ${boostsources} missing, exiting."
    exit 1
fi

if [ -d "${boostdir}" ]; then
    echo "Old boost library exists (${boostdir}), removing it ."
    rm -rf "${boostdir}"
fi

for i in $(find "${libdir}" -iname *.tar.gz); do
  echo "Old archive ($i) exists, removing it."
  rm -rf "$i"
done

## (6) Unpack boost -- note that for tarballs straight from Debian we need a rename step
echo "Unpacking ${boosttargz} into ${libdir}."
#(cd ${libdir} && tar xfz ${boostsources} && mv boost*-*.orig/ ${boostver})
(cd "${libdir}" && tar xfz "${boostsources}")
# exit -1

## (7) Install Boost libraries
echo "Copying Boost libraries into ${libdir}"

boostlibs="bind concept config container date_time detail exception functional integer \
           interprocess intrusive io iterator math move mpl numeric pending preprocessor \
           random range smart_ptr tuple typeof type_traits unordered utility uuid \
           filesystem spirit foreach algorithm iostreams \
           dynamic_bitset heap any circular_buffer geometry fusion graph \
           multiprecision phoenix bimap icl flyweight property_tree \
           scope_exit atomic align sort compute"

## this copies the Boost libraries listed in ${boostlibs} from the
## Boost sources in ${boostroot} into the target directory ${libdir}
bcp --boost="${boostroot}"  "${boostlibs}"  "${libdir}"   > /dev/null  2>&1

## (8) Some post processing and cleanup
echo "Post-processing and cleanup"
rm -rf "${libdir}/libs" \
       "${libdir}/Jamroot" \
       "${libdir}/boost.png" \
       "${libdir}/doc" \
       "${libdir}/boost.css" \
       "${libdir}/rst.css" \
       "${boostroot}"

## (9) Create archive to be committed
(cd "${libdir}" && tar -zcf "${libdir}/${boosttargz}" boost)

## (10) And done
echo "Now check with 'git status' and add and commit as needed."

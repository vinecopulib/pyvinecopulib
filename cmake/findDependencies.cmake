include(cmake/findEigen3.cmake            REQUIRED)
message(STATUS "Eigen3 version: ${EIGEN3_VERSION}")
find_package(Boost 1.56                   REQUIRED)
find_package(Threads                      REQUIRED)

set(wdm_INCLUDE_DIRS "lib/wdm/include")
set(external_includes ${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${wdm_INCLUDE_DIRS})

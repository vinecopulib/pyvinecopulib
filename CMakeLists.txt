cmake_minimum_required(VERSION 3.10)

project(FindLibraries)

find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)

# Output the include directories
message(STATUS "EIGEN3_INCLUDE_DIRS=${EIGEN3_INCLUDE_DIRS}")
message(STATUS "Boost_INCLUDE_DIRS=${Boost_INCLUDE_DIRS}")
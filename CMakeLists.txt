cmake_minimum_required(VERSION 2.8.12)
project(pyvinecopulib)

add_subdirectory(lib/pybind11)
pybind11_add_module(pyvinecopulib src/main.cpp)

find_package(Threads                      REQUIRED)
set(eigen3_INCLUDE_DIRS "lib/eigen")
set(boost_INCLUDE_DIRS "lib/boost")
set(wdm_INCLUDE_DIRS "lib/wdm/include")
set(vinecopulib_INCLUDE_DIRS "lib/vinecopulib/include")
set(external_INCLUDE_DIRS ${eigen3_INCLUDE_DIRS} ${boost_INCLUDE_DIRS} 
  ${wdm_INCLUDE_DIRS} ${vinecopulib_INCLUDE_DIRS})

include_directories(src)
include_directories(SYSTEM ${external_INCLUDE_DIRS})


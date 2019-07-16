cmake_minimum_required(VERSION 2.8.12)
project(pyvinecopulib)

add_subdirectory(lib/pybind11)
pybind11_add_module(pyvinecopulib src/main.cpp)

include(cmake/findDependencies.cmake      REQUIRED)

include_directories(src)
include_directories(SYSTEM ${external_includes})


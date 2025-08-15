#pragma once

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>

#include <kde1d.hpp>

#include "docstr.hpp"

namespace nb = nanobind;
using namespace nb::literals;

inline void init_kde1d(nb::module_& module) {}

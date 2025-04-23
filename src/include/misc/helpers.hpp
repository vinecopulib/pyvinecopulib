#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <string>

namespace nb = nanobind;

inline std::string
get_helper_doc(const std::string& module,
               const std::string& attr,
               const std::string& fallback)
{
  try {
    auto mod = nb::module_::import_(module.c_str());
    return nb::cast<std::string>(mod.attr(attr.c_str()));
  } catch (...) {
    return fallback;
  }
}
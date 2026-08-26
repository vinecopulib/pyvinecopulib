#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <string>
#include <vinecopulib.hpp>

namespace nb = nanobind;
using namespace vinecopulib;

inline std::string python_doc_helper(const std::string& module,
                                     const std::string& attr,
                                     const std::string& fallback) {
  try {
    auto mod = nb::module_::import_(module.c_str());
    return nb::cast<std::string>(mod.attr(attr.c_str()));
  } catch (...) {
    return fallback;
  }
}

// Appends a note to the generated docstring of a method that takes or returns
// a vinecopulib ``TriangularArray``: the Python binding represents it as a
// nested list.
inline std::string struct_array_list_doc(const char* base_doc) {
  return std::string(base_doc) +
         R"""(

Notes
-----
The triangular array is represented in Python as a nested list indexed
``[tree][edge]``, where tree ``i`` holds ``d - 1 - i`` entries.
)""";
}

// template<typename T>
// inline std::string
// python_str_helper(const T& obj, const std::string& label)
// {
//   std::string full = obj.str();
//   auto pos = full.find('\n');
//   return pos != std::string::npos ? "<" + label + ">\n" + full.substr(pos +
//   1)
//                                   : "<" + label + ">";
// }

// Defaults for the `controls` arguments. A *constructed* default argument is
// materialized by nanobind as a Python object owned by the function record,
// which then outlives module teardown and is reported as a leaked instance;
// a nullable pointer costs neither that nor the copy an `std::optional` by
// value would force on every call that does pass controls.
inline const FitControlsBicop& default_bicop_controls() {
  static const FitControlsBicop kDefaults{};
  return kDefaults;
}

inline const FitControlsVinecop& default_vinecop_controls() {
  static const FitControlsVinecop kDefaults{};
  return kDefaults;
}

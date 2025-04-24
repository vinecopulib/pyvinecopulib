#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <string>

namespace nb = nanobind;

inline std::string
python_doc_helper(const std::string& module,
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
#include "pyvinecopulib.hpp"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vinecopulib.hpp>

namespace py = pybind11;
using namespace vinecopulib;

PYBIND11_MODULE(_pyvinecopulib, pv)
{

  pv.doc() = R"pbdoc(
  The pyvinecopulib package
  -------------------------
  )pbdoc";

  init_stats(pv);

  init_bicop_family(pv);
  init_bicop_fit_controls(pv);
  init_bicop_class(pv);

  init_vinecop_rvine_structure(pv);
  init_vinecop_fit_controls(pv);
  init_vinecop_class(pv);

#ifdef VERSION_INFO
  pv.attr("__version__") = VERSION_INFO;
#else
  pv.attr("__version__") = "dev";
#endif
}

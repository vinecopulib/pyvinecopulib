#include "bicop/class.hpp"
#include "bicop/family.hpp"
#include "bicop/fit_controls.hpp"
#include "misc/stats.hpp"
#include "vinecop/class.hpp"
#include "vinecop/fit_controls.hpp"
#include "vinecop/rvine_structure.hpp"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vinecopulib.hpp>

namespace py = pybind11;
using namespace vinecopulib;

PYBIND11_MODULE(pyvinecopulib, pv)
{
  pv.doc() = R"pbdoc(
  The pyvinecopulib package
  -------------------------
  )pbdoc";

  py::module_ stats =
    pv.def_submodule("stats", "Misc statistics tools for vine copulas");
  init_stats(stats);

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

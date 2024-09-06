#include "pyvinecopulib.hpp"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vinecopulib.hpp>

namespace py = pybind11;
using namespace vinecopulib;

/// Overrides the `__name__` of a module.  Classes defined by pybind11 use the
/// `__name__` of the module as of the time they are defined, which affects the
/// `__repr__` of the class type objects.
/// See
/// https://github.com/google/tensorstore/blob/94be4f2e8715511bb60fc0a0eaf07335881673b3/python/tensorstore/tensorstore.cc#L75
class ScopedModuleNameOverride
{
public:
  explicit ScopedModuleNameOverride(const py::module& m,
                                    const std::string& name)
    : module_(std::move(m))
  {
    original_name_ = module_.attr("__name__");
    module_.attr("__name__") = name;
  }
  ~ScopedModuleNameOverride() { module_.attr("__name__") = original_name_; }

private:
  py::module module_;
  py::object original_name_;
};

PYBIND11_MODULE(_pyvinecopulib, pv)
{

  // Ensure that members of this module display as `pyvinecopulib.X` rather than
  // `pyvinecopulib.pyvinecopulib.X`.
  ScopedModuleNameOverride name_override(pv, "pyvinecopulib");

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

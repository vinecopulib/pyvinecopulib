#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>

#include "pyvinecopulib.hpp"
// #include <nanobind/stl.h>
#include <vinecopulib.hpp>

namespace nb = nanobind;
using namespace vinecopulib;

/// Overrides the `__name__` of a module.  Classes defined by pybind11 use the
/// `__name__` of the module as of the time they are defined, which affects the
/// `__repr__` of the class type objects.
/// See
/// https://github.com/google/tensorstore/blob/94be4f2e8715511bb60fc0a0eaf07335881673b3/python/tensorstore/tensorstore.cc#L75
class ScopedModuleNameOverride {
 public:
  explicit ScopedModuleNameOverride(const nb::module_& m,
                                    const std::string& name)
      : module_(std::move(m)), original_name_(module_.attr("__name__")) {
    module_.attr("__name__") = name;
  }
  ~ScopedModuleNameOverride() { module_.attr("__name__") = original_name_; }

 private:
  nb::module_ module_;
  nb::object original_name_;
};

NB_MODULE(pyvinecopulib_ext, pv) {
  pv.doc() = R"pbdoc(
  The pyvinecopulib package
  -------------------------
  )pbdoc";

  // Bind each class with __module__ set to the subpackage where the symbol
  // is canonically located (mirrors the pure-Python pyvinecopulib.{core,
  // families,utils} layout). Drives the module path embedded in pickles —
  // the corresponding repr strings are hardcoded per-class to match.
  //
  // Order matters: init_bicop_family (binds BicopFamily) must come before
  // init_bicop_class / init_bicop_fit_controls, which reference it.
  {
    ScopedModuleNameOverride n(pv, "pyvinecopulib.families");
    init_bicop_family(pv);
  }

  {
    ScopedModuleNameOverride n(pv, "pyvinecopulib.core");
    init_bicop_fit_controls(pv);
    init_bicop_class(pv);
    init_vinecop_rvine_structure(pv);
    init_vinecop_fit_controls(pv);
    init_vinecop_class(pv);
  }

  {
    ScopedModuleNameOverride n(pv, "pyvinecopulib.utils");
    init_stats(pv);
    init_benchmark(pv);
    init_kde1d(pv);
  }

#ifdef VERSION_INFO
  pv.attr("__version__") = VERSION_INFO;
#else
  pv.attr("__version__") = "dev";
#endif
}

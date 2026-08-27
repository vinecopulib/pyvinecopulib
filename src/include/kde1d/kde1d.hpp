#pragma once

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <cmath>
#include <kde1d.hpp>
#include <limits>
#include <tuple>

#include "docstr.hpp"
#include "misc/helpers.hpp"

namespace nb = nanobind;
using namespace nb::literals;
using namespace kde1d;

constexpr auto& kde1d_doc = pyvinecopulib_doc.kde1d.Kde1d;

// Python-binding-only docstring for the unified `__init__` factory — the
// upstream C++ class has four constructor overloads that libclang cannot
// disambiguate (so the auto-extracted `kde1d_doc.ctor.doc_*` falls back to
// `doc_was_unable_to_choose_unambiguous_names`). Keep this hand-written
// since the Python signature does not match any single C++ overload.
constexpr const char* kde1d_constructor_doc = R"""(
Constructs a `Kde1d` instance.

Parameters
----------
xmin :
    Lower bound for the support of the density. `NaN` means no
    boundary.
xmax :
    Upper bound for the support of the density. `NaN` means no
    boundary.
type :
    Variable type. One of ``"continuous"``, ``"discrete"``, or
    ``"zero_inflated"``.
multiplier :
    Bandwidth multiplier. The actual bandwidth used is
    ``bandwidth * multiplier``.
bandwidth :
    Bandwidth parameter. `None` / `NaN` selects the bandwidth
    automatically via the plug-in methodology.
degree :
    Degree of the local polynomial (0, 1, or 2 for log-constant,
    log-linear, or log-quadratic fitting).
grid_size :
    Number of grid points for the interpolation grid (must be
    >= 4).
boundary_repair :
    Whether a finite bound may be fitted with a dedicated boundary
    estimator instead of the transformed bulk fit. Eligibility, not a
    guarantee: it needs a finite bound and at least 16 observations,
    and an endpoint whose local behavior is unclear keeps the bulk
    estimate. Has no effect when neither bound is set.
)""";

// Factory function to create a Kde1d from xmin, xmax, type string, multiplier,
// bandwidth, degree
inline Kde1d kde1d_from_params(std::optional<double> xmin = std::nullopt,
                               std::optional<double> xmax = std::nullopt,
                               const std::string& type = "continuous",
                               double multiplier = 1.0,
                               std::optional<double> bandwidth = std::nullopt,
                               size_t degree = 2, size_t grid_size = 400,
                               bool boundary_repair = true) {
  return Kde1d(xmin.value_or(NAN), xmax.value_or(NAN), type, multiplier,
               bandwidth.value_or(NAN), degree, grid_size, boundary_repair);
}

// Factory function to create a Kde1d from grid, xmin, xmax, type string, prob0
inline Kde1d kde1d_from_grid(const Eigen::VectorXd& grid_points,
                             const Eigen::VectorXd& values,
                             std::optional<double> xmin = std::nullopt,
                             std::optional<double> xmax = std::nullopt,
                             const std::string& type = "continuous",
                             double prob0 = 0.0) {
  interp::InterpolationGrid grid(grid_points, values, 0);
  return Kde1d(grid, xmin.value_or(NAN), xmax.value_or(NAN), type, prob0);
}

// Whether `fit` has populated the estimator. The upstream C++ class has no
// such accessor, so it is derived from the interpolation grid here; a
// `Kde1d::is_fitted()` upstream would retire this.
inline bool kde1d_is_fitted(const Kde1d& kde) {
  return kde.get_grid_points().size() > 0;
}

// `Kde1d` satisfies the `MarginLike` contract directly, the way `Bicop` and
// `Vinecop` satisfy theirs. These helpers translate what the estimator already
// knows into the margin layer's vocabulary; none of it is new modeling.

// Variable type in the margin layer's spelling.
inline std::string kde1d_var_type(const Kde1d& kde) {
  const std::string t = kde.get_type_str();
  if (t == "discrete") return "d";
  // The enum is `zero_inflated`; `as_str` renders it hyphenated.
  if (t == "zero_inflated" || t == "zero-inflated") return "zi";
  return "c";
}

// Support as a pair, with an absent bound as an infinity rather than `NaN`, so
// comparisons against it are total.
inline std::tuple<double, double> kde1d_support(const Kde1d& kde) {
  const double lo = kde.get_xmin();
  const double hi = kde.get_xmax();
  constexpr double inf = std::numeric_limits<double>::infinity();
  return {std::isnan(lo) ? -inf : lo, std::isnan(hi) ? inf : hi};
}

// Left limit `F(x^-)`, which a variable with atoms needs and a continuous one
// does not. Computed here because the estimator already knows its own type and
// zero-inflation probability.
inline Eigen::VectorXd kde1d_cdf_left(const Kde1d& kde,
                                      const Eigen::VectorXd& x) {
  const std::string t = kde.get_type_str();
  if (t == "discrete") {
    const Eigen::VectorXd below = x.array() - 1.0;
    return kde.cdf(below);
  }
  if (t == "zero_inflated" || t == "zero-inflated") {
    Eigen::VectorXd out = kde.cdf(x);
    const Eigen::VectorXd mass = kde.pdf(x);
    for (Eigen::Index i = 0; i < x.size(); ++i) {
      if (x[i] == 0.0) out[i] -= mass[i];
    }
    return out;
  }
  return kde.cdf(x);
}

// Log-density. A margin's consumers read `logpdf` when it exists rather than
// taking a log themselves, and the estimator can do it without a Python round
// trip. Zero density gives `-inf`, which is the honest answer.
inline Eigen::VectorXd kde1d_logpdf(const Kde1d& kde,
                                    const Eigen::VectorXd& x) {
  return kde.pdf(x).array().log();
}

// Log-likelihood: the value attained at the fit when called without data, the
// log-density sum otherwise. Same shape as `Bicop.loglik` and `Vinecop.loglik`.
inline double kde1d_loglik(const Kde1d& kde, const Eigen::VectorXd& x) {
  if (x.size() == 0) return kde.get_loglik();
  return kde.pdf(x).array().log().sum();
}

// Wrapper function for set_xmin_xmax with optional parameters
inline void kde1d_set_xmin_xmax(Kde1d& self,
                                std::optional<double> xmin = std::nullopt,
                                std::optional<double> xmax = std::nullopt) {
  self.set_xmin_xmax(xmin.value_or(NAN), xmax.value_or(NAN));
}

// Wrapper function to call the Python kde1d_plot function
inline void kde1d_plot_wrapper(const Kde1d& kde, nb::object xlim,
                               nb::object ylim, int grid_size,
                               bool show_zero_mass) {
  auto mod = nb::module_::import_("pyvinecopulib._python_helpers.kde1d");
  auto kde1d_plot = mod.attr("kde1d_plot");
  kde1d_plot(nb::cast(kde), xlim, ylim, grid_size, show_zero_mass);
}

inline void init_kde1d(nb::module_& module) {
  auto cls =
      nb::class_<Kde1d>(module, "Kde1d", kde1d_doc.doc)
          // Default constructor
          .def(
              "__init__",
              [](Kde1d* self, std::optional<double> xmin,
                 std::optional<double> xmax, const std::string& type,
                 double multiplier, std::optional<double> bandwidth,
                 size_t degree, size_t grid_size, bool boundary_repair) {
                new (self) Kde1d(xmin.value_or(NAN), xmax.value_or(NAN), type,
                                 multiplier, bandwidth.value_or(NAN), degree,
                                 grid_size, boundary_repair);
              },
              "xmin"_a = std::nullopt, "xmax"_a = std::nullopt,
              "type"_a = "continuous", "multiplier"_a = 1.0,
              "bandwidth"_a = std::nullopt, "degree"_a = 2, "grid_size"_a = 400,
              "boundary_repair"_a = true, kde1d_constructor_doc,
              nb::call_guard<nb::gil_scoped_release>())
          .def_static("from_params", &kde1d_from_params,
                      "xmin"_a = std::nullopt, "xmax"_a = std::nullopt,
                      "type"_a = "continuous", "multiplier"_a = 1.0,
                      "bandwidth"_a = std::nullopt, "degree"_a = 2,
                      "grid_size"_a = 400, "boundary_repair"_a = true,
                      kde1d_constructor_doc,
                      nb::call_guard<nb::gil_scoped_release>())
          .def_static("from_grid", &kde1d_from_grid, "grid_points"_a,
                      "values"_a, "xmin"_a = std::nullopt,
                      "xmax"_a = std::nullopt, "type"_a = "continuous",
                      "prob0"_a = 0.0,
                      "Constructs a `Kde1d` from a pre-computed interpolation "
                      "grid (skipping the kernel-density fit).",
                      nb::call_guard<nb::gil_scoped_release>())

          // Properties (getters) — auto-extracted from `lib/kde1d` upstream
          // `//!` comments.
          .def_prop_ro("xmin", &Kde1d::get_xmin, kde1d_doc.get_xmin.doc)
          .def_prop_ro("xmax", &Kde1d::get_xmax, kde1d_doc.get_xmax.doc)
          .def_prop_ro("type", &Kde1d::get_type_str, kde1d_doc.get_type_str.doc)
          .def_prop_ro("prob0", &Kde1d::get_prob0, kde1d_doc.get_prob0.doc)
          .def_prop_ro("multiplier", &Kde1d::get_multiplier,
                       kde1d_doc.get_multiplier.doc)
          .def_prop_ro("bandwidth", &Kde1d::get_bandwidth,
                       kde1d_doc.get_bandwidth.doc)
          .def_prop_ro("degree", &Kde1d::get_degree, kde1d_doc.get_degree.doc)
          .def_prop_ro("grid_size", &Kde1d::get_grid_size,
                       kde1d_doc.get_grid_size.doc)
          .def_prop_ro("actual_grid_size", &Kde1d::get_actual_grid_size,
                       kde1d_doc.get_actual_grid_size.doc)
          .def_prop_ro("boundary_repair", &Kde1d::get_boundary_repair,
                       kde1d_doc.get_boundary_repair.doc)
          .def_prop_ro("edf", &Kde1d::get_edf, kde1d_doc.get_edf.doc)
          .def_prop_ro("grid_points", &Kde1d::get_grid_points,
                       kde1d_doc.get_grid_points.doc,
                       nb::call_guard<nb::gil_scoped_release>())
          .def_prop_ro("values", &Kde1d::get_values, kde1d_doc.get_values.doc,
                       nb::call_guard<nb::gil_scoped_release>())
          .def_prop_ro("is_fitted", &kde1d_is_fitted,
                       "Whether the estimator has been fitted to data.")
          .def_prop_ro(
              "var_type", &kde1d_var_type,
              "Variable type as the margin layer spells it: ``\"c\"``, "
              "``\"d\"``, or ``\"zi\"``.")
          .def_prop_ro("support", &kde1d_support,
                       "Lower and upper bounds of the support, as a pair. An "
                       "absent bound is an infinity rather than `NaN`.")
          .def_prop_ro(
              "family_name", [](const Kde1d&) { return std::string("kde1d"); },
              "Family name, as a selection report spells it.")
          .def_prop_ro("n_parameters", &Kde1d::get_edf,
                       "Effective degrees of freedom, under the margin layer's "
                       "name for a parameter count. See also ``edf``.")

          // Methods — auto-extracted from `lib/kde1d` upstream `//!` comments.
          // Returns `self` so a fit chains, as it does on every Python
          // estimator. The GIL is released around the fit itself rather than by
          // a call guard, because handing back the object needs it.
          .def(
              "fit",
              [](Kde1d& self, const Eigen::VectorXd& x,
                 const Eigen::VectorXd& weights) -> Kde1d& {
                {
                  nb::gil_scoped_release release;
                  self.fit(x, weights);
                }
                return self;
              },
              "x"_a, "weights"_a = Eigen::VectorXd(), kde1d_doc.fit.doc,
              nb::rv_policy::reference_internal)
          .def("pdf", &Kde1d::pdf, "x"_a, "check_fitted"_a = true,
               kde1d_doc.pdf.doc, nb::call_guard<nb::gil_scoped_release>())
          .def("cdf", &Kde1d::cdf, "x"_a, "check_fitted"_a = true,
               kde1d_doc.cdf.doc, nb::call_guard<nb::gil_scoped_release>())
          .def("logpdf", &kde1d_logpdf, "x"_a,
               "Log of the density, or of the probability mass at an atom.\n"
               "\n"
               "Parameters\n"
               "----------\n"
               "x : array, shape (n,), dtype float\n"
               "    Evaluation points.\n"
               "\n"
               "Returns\n"
               "-------\n"
               "array, shape (n,), dtype float\n"
               "    Log-densities; ``-inf`` where the density is zero.",
               nb::call_guard<nb::gil_scoped_release>())
          .def("cdf_left", &kde1d_cdf_left, "x"_a,
               "Left limit ``F(x^-)`` of the distribution function.\n"
               "\n"
               "Parameters\n"
               "----------\n"
               "x : array, shape (n,), dtype float\n"
               "    Evaluation points.\n"
               "\n"
               "Returns\n"
               "-------\n"
               "array, shape (n,), dtype float\n"
               "    ``F(x)`` for a continuous variable, ``F(x - 1)`` for a "
               "discrete one, and\n"
               "    ``F(x) - f(x)`` at the atom of a zero-inflated one.",
               nb::call_guard<nb::gil_scoped_release>())
          .def("loglik", &kde1d_loglik, "x"_a = Eigen::VectorXd(),
               "Log-likelihood attained at the fit, or of given data.\n"
               "\n"
               "Parameters\n"
               "----------\n"
               "x : array, shape (n,), dtype float, optional\n"
               "    Observations. Omitted, the value attained at the fit is "
               "returned.\n"
               "\n"
               "Returns\n"
               "-------\n"
               "float\n"
               "    The log-likelihood.",
               nb::call_guard<nb::gil_scoped_release>())
          .def("icdf", &Kde1d::quantile, "x"_a, "check_fitted"_a = true,
               kde1d_doc.quantile.doc, nb::call_guard<nb::gil_scoped_release>())
          .def("sample", &Kde1d::simulate, "n"_a,
               "seeds"_a = std::vector<int>(), "check_fitted"_a = true,
               kde1d_doc.simulate.doc, nb::call_guard<nb::gil_scoped_release>())
          .def("set_xmin_xmax", &kde1d_set_xmin_xmax, "xmin"_a = std::nullopt,
               "xmax"_a = std::nullopt, kde1d_doc.set_xmin_xmax.doc)
          .def("plot", &kde1d_plot_wrapper, "xlim"_a = nb::none(),
               "ylim"_a = nb::none(), "grid_size"_a = 200,
               "show_zero_mass"_a = true,
               python_doc_helper("pyvinecopulib._python_helpers.kde1d",
                                 "KDE1D_PLOT_DOC",
                                 "Plot the KDE (extended doc unavailable) ")
                   .c_str())

          // String representation
          .def(
              "__repr__",
              [](const Kde1d& kde) {
                return "<pyvinecopulib.core.Kde1d> " + kde.str();
              },
              "Return string representation of the Kde1d object.")
          .def(
              "__str__",
              [](const Kde1d& kde) {
                return "<pyvinecopulib.core.Kde1d> " + kde.str();
              },
              "Return string representation of the Kde1d object.")

          // Serialization support
          .def("__getstate__",
               [](const Kde1d& kde) {
                 nb::dict s;
                 const bool fitted = kde1d_is_fitted(kde);
                 s["fitted"] = fitted;
                 s["xmin"] = kde.get_xmin();
                 s["xmax"] = kde.get_xmax();
                 s["type"] = kde.get_type_str();
                 if (fitted) {
                   // For fitted models: save all data needed to reconstruct
                   s["prob0"] = kde.get_prob0();
                   s["grid_points"] = kde.get_grid_points();
                   s["values"] = kde.get_values();
                 } else {
                   // For unfitted models: save parameters only
                   s["multiplier"] = kde.get_multiplier();
                   s["bandwidth"] = kde.get_bandwidth();
                   s["degree"] = static_cast<std::size_t>(kde.get_degree());
                   s["grid_size"] =
                       static_cast<std::size_t>(kde.get_grid_size());
                   s["boundary_repair"] = kde.get_boundary_repair();
                 }
                 return s;
               })

          .def("__setstate__", [](Kde1d& kde, nb::dict s) {
            const bool fitted = nb::cast<bool>(s["fitted"]);
            const double xmin = nb::cast<double>(s["xmin"]);
            const double xmax = nb::cast<double>(s["xmax"]);
            const std::string type = nb::cast<std::string>(s["type"]);

            if (fitted) {
              const double prob0 = nb::cast<double>(s["prob0"]);
              const Eigen::VectorXd grid_points =
                  nb::cast<Eigen::VectorXd>(s["grid_points"]);
              const Eigen::VectorXd values =
                  nb::cast<Eigen::VectorXd>(s["values"]);
              // Create interpolation grid and construct object
              interp::InterpolationGrid grid(grid_points, values, 0);
              new (&kde) Kde1d(grid, xmin, xmax, type, prob0);
            } else {
              // For unfitted models, construct from parameters
              const double multiplier = nb::cast<double>(s["multiplier"]);
              const double bandwidth = nb::cast<double>(s["bandwidth"]);
              const std::size_t degree = nb::cast<std::size_t>(s["degree"]);
              const std::size_t grid_size =
                  nb::cast<std::size_t>(s["grid_size"]);
              const bool boundary_repair =
                  nb::cast<bool>(s["boundary_repair"]);
              new (&kde) Kde1d(xmin, xmax, type, multiplier, bandwidth, degree,
                               grid_size, boundary_repair);
            }
          });

  // Declared on the class, matching `MarginBase.supports_weights`, so the
  // margin layer's capability read works off the class or an instance.
  cls.attr("supports_weights") = true;
}

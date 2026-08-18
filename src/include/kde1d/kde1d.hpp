#pragma once

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <kde1d.hpp>

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
)""";

// Factory function to create a Kde1d from xmin, xmax, type string, multiplier,
// bandwidth, degree
inline Kde1d kde1d_from_params(std::optional<double> xmin = std::nullopt,
                               std::optional<double> xmax = std::nullopt,
                               const std::string& type = "continuous",
                               double multiplier = 1.0,
                               std::optional<double> bandwidth = std::nullopt,
                               size_t degree = 2, size_t grid_size = 400) {
  return Kde1d(xmin.value_or(NAN), xmax.value_or(NAN), type, multiplier,
               bandwidth.value_or(NAN), degree, grid_size);
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
  nb::class_<Kde1d>(module, "Kde1d", kde1d_doc.doc)
      // Default constructor
      .def(
          "__init__",
          [](Kde1d* self, std::optional<double> xmin,
             std::optional<double> xmax, const std::string& type,
             double multiplier, std::optional<double> bandwidth, size_t degree,
             size_t grid_size) {
            new (self)
                Kde1d(xmin.value_or(NAN), xmax.value_or(NAN), type, multiplier,
                      bandwidth.value_or(NAN), degree, grid_size);
          },
          "xmin"_a = std::nullopt, "xmax"_a = std::nullopt,
          "type"_a = "continuous", "multiplier"_a = 1.0,
          "bandwidth"_a = std::nullopt, "degree"_a = 2, "grid_size"_a = 400,
          kde1d_constructor_doc, nb::call_guard<nb::gil_scoped_release>())
      .def_static("from_params", &kde1d_from_params, "xmin"_a = std::nullopt,
                  "xmax"_a = std::nullopt, "type"_a = "continuous",
                  "multiplier"_a = 1.0, "bandwidth"_a = std::nullopt,
                  "degree"_a = 2, "grid_size"_a = 400, kde1d_constructor_doc,
                  nb::call_guard<nb::gil_scoped_release>())
      .def_static("from_grid", &kde1d_from_grid, "grid_points"_a, "values"_a,
                  "xmin"_a = std::nullopt, "xmax"_a = std::nullopt,
                  "type"_a = "continuous", "prob0"_a = 0.0,
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
      .def_prop_ro("loglik", &Kde1d::get_loglik, kde1d_doc.get_loglik.doc)
      .def_prop_ro("edf", &Kde1d::get_edf, kde1d_doc.get_edf.doc)
      .def_prop_ro("grid_points", &Kde1d::get_grid_points,
                   kde1d_doc.get_grid_points.doc,
                   nb::call_guard<nb::gil_scoped_release>())
      .def_prop_ro("values", &Kde1d::get_values, kde1d_doc.get_values.doc,
                   nb::call_guard<nb::gil_scoped_release>())
      .def_prop_ro("is_fitted", &kde1d_is_fitted,
                   "Whether the estimator has been fitted to data.")

      // Methods — auto-extracted from `lib/kde1d` upstream `//!` comments.
      .def("fit", &Kde1d::fit, "x"_a, "weights"_a = Eigen::VectorXd(),
           kde1d_doc.fit.doc, nb::call_guard<nb::gil_scoped_release>())
      .def("pdf", &Kde1d::pdf, "x"_a, "check_fitted"_a = true,
           kde1d_doc.pdf.doc, nb::call_guard<nb::gil_scoped_release>())
      .def("cdf", &Kde1d::cdf, "x"_a, "check_fitted"_a = true,
           kde1d_doc.cdf.doc, nb::call_guard<nb::gil_scoped_release>())
      .def("icdf", &Kde1d::quantile, "x"_a, "check_fitted"_a = true,
           kde1d_doc.quantile.doc, nb::call_guard<nb::gil_scoped_release>())
      .def("simulate", &Kde1d::simulate, "n"_a, "seeds"_a = std::vector<int>(),
           "check_fitted"_a = true, kde1d_doc.simulate.doc,
           nb::call_guard<nb::gil_scoped_release>())
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
            return "<pyvinecopulib.utils.Kde1d> " + kde.str();
          },
          "Return string representation of the Kde1d object.")
      .def(
          "__str__",
          [](const Kde1d& kde) {
            return "<pyvinecopulib.utils.Kde1d> " + kde.str();
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
               s["grid_size"] = static_cast<std::size_t>(kde.get_grid_size());
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
          const Eigen::VectorXd values = nb::cast<Eigen::VectorXd>(s["values"]);
          // Create interpolation grid and construct object
          interp::InterpolationGrid grid(grid_points, values, 0);
          new (&kde) Kde1d(grid, xmin, xmax, type, prob0);
        } else {
          // For unfitted models, construct from parameters
          const double multiplier = nb::cast<double>(s["multiplier"]);
          const double bandwidth = nb::cast<double>(s["bandwidth"]);
          const std::size_t degree = nb::cast<std::size_t>(s["degree"]);
          const std::size_t grid_size = nb::cast<std::size_t>(s["grid_size"]);
          new (&kde)
              Kde1d(xmin, xmax, type, multiplier, bandwidth, degree, grid_size);
        }
      });
}

#pragma once

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <kde1d.hpp>

#include "kde1d/docstr.hpp"

namespace nb = nanobind;
using namespace nb::literals;
using namespace kde1d;

// Factory function to create a Kde1d from xmin, xmax, type string, multiplier,
// bandwidth, degree
inline Kde1d kde1d_from_params(std::optional<double> xmin = std::nullopt,
                               std::optional<double> xmax = std::nullopt,
                               const std::string& type = "continuous",
                               double multiplier = 1.0,
                               std::optional<double> bandwidth = std::nullopt,
                               size_t degree = 2) {
  return Kde1d(xmin.value_or(NAN), xmax.value_or(NAN), type, multiplier,
               bandwidth.value_or(NAN), degree);
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

// Wrapper function for set_xmin_xmax with optional parameters
inline void kde1d_set_xmin_xmax(Kde1d& self,
                                std::optional<double> xmin = std::nullopt,
                                std::optional<double> xmax = std::nullopt) {
  self.set_xmin_xmax(xmin.value_or(NAN), xmax.value_or(NAN));
}

inline void init_kde1d(nb::module_& module) {
  nb::class_<Kde1d>(module, "Kde1d", kde1d_docstrings::kde1d_class_doc)
      // Default constructor
      .def(
          "__init__",
          [](Kde1d* self, std::optional<double> xmin,
             std::optional<double> xmax, const std::string& type,
             double multiplier, std::optional<double> bandwidth,
             size_t degree) {
            new (self) Kde1d(xmin.value_or(NAN), xmax.value_or(NAN), type,
                             multiplier, bandwidth.value_or(NAN), degree);
          },
          "xmin"_a = std::nullopt, "xmax"_a = std::nullopt,
          "type"_a = "continuous", "multiplier"_a = 1.0,
          "bandwidth"_a = std::nullopt, "degree"_a = 2,
          kde1d_docstrings::kde1d_constructor_doc)
      .def_static("from_params", &kde1d_from_params, "xmin"_a = std::nullopt,
                  "xmax"_a = std::nullopt, "type"_a = "continuous",
                  "multiplier"_a = 1.0, "bandwidth"_a = std::nullopt,
                  "degree"_a = 2, kde1d_docstrings::kde1d_from_params_doc)
      .def_static("from_grid", &kde1d_from_grid, "grid_points"_a, "values"_a,
                  "xmin"_a = std::nullopt, "xmax"_a = std::nullopt,
                  "type"_a = "continuous", "prob0"_a = 0.0,
                  kde1d_docstrings::kde1d_from_grid_doc)

      // Properties (getters)
      .def_prop_ro("xmin", &Kde1d::get_xmin, kde1d_docstrings::xmin_doc)
      .def_prop_ro("xmax", &Kde1d::get_xmax, kde1d_docstrings::xmax_doc)
      .def_prop_ro("type", &Kde1d::get_type_str, kde1d_docstrings::type_doc)
      .def_prop_ro("prob0", &Kde1d::get_prob0, kde1d_docstrings::prob0_doc)
      .def_prop_ro("multiplier", &Kde1d::get_multiplier,
                   kde1d_docstrings::multiplier_doc)
      .def_prop_ro("bandwidth", &Kde1d::get_bandwidth,
                   kde1d_docstrings::bandwidth_doc)
      .def_prop_ro("degree", &Kde1d::get_degree, kde1d_docstrings::degree_doc)
      .def_prop_ro("loglik", &Kde1d::get_loglik, kde1d_docstrings::loglik_doc)
      .def_prop_ro("edf", &Kde1d::get_edf, kde1d_docstrings::edf_doc)
      .def_prop_ro("grid_points", &Kde1d::get_grid_points,
                   kde1d_docstrings::grid_points_doc)
      .def_prop_ro("values", &Kde1d::get_values, kde1d_docstrings::values_doc)

      // Methods
      .def("fit", &Kde1d::fit, "x"_a, "weights"_a = Eigen::VectorXd(),
           kde1d_docstrings::fit_doc)
      .def("pdf", &Kde1d::pdf, "x"_a, "check_fitted"_a = true,
           kde1d_docstrings::pdf_doc)
      .def("cdf", &Kde1d::cdf, "x"_a, "check_fitted"_a = true,
           kde1d_docstrings::cdf_doc)
      .def("quantile", &Kde1d::quantile, "x"_a, "check_fitted"_a = true,
           kde1d_docstrings::quantile_doc)
      .def("simulate", &Kde1d::simulate, "n"_a, "seeds"_a = std::vector<int>(),
           "check_fitted"_a = true, kde1d_docstrings::simulate_doc)
      .def("set_xmin_xmax", &kde1d_set_xmin_xmax, "xmin"_a = std::nullopt,
           "xmax"_a = std::nullopt, kde1d_docstrings::set_xmin_xmax_doc)

      // String representation
      .def(
          "__repr__",
          [](const Kde1d& kde) { return "<pyvinecopulib.Kde1d> " + kde.str(); },
          "Return string representation of the Kde1d object.")
      .def(
          "__str__",
          [](const Kde1d& kde) { return "<pyvinecopulib.Kde1d> " + kde.str(); },
          "Return string representation of the Kde1d object.")

      // Serialization support
      .def("__getstate__",
           [](const Kde1d& kde) -> nb::object {
             // Check if model is fitted by testing if loglik is not NaN
             bool is_fitted = !(kde.get_loglik() != kde.get_loglik());

             if (is_fitted) {
               // For fitted models: save all data needed to reconstruct
               return nb::make_tuple(kde.get_xmin(), kde.get_xmax(),
                                     kde.get_type_str(), kde.get_prob0(),
                                     kde.get_grid_points(), kde.get_values(),
                                     kde.get_loglik(), kde.get_edf(),
                                     true);  // fitted flag
             } else {
               // For unfitted models: save parameters only
               return nb::make_tuple(kde.get_xmin(), kde.get_xmax(),
                                     kde.get_type_str(), kde.get_multiplier(),
                                     kde.get_bandwidth(),
                                     static_cast<double>(kde.get_degree()),
                                     false);  // fitted flag
             }
           })
      .def("__setstate__", [](Kde1d& kde, nb::tuple state) {
        bool is_fitted = nb::cast<bool>(state[state.size() - 1]);

        if (is_fitted) {
          // For fitted models, construct from grid
          double xmin = nb::cast<double>(state[0]);
          double xmax = nb::cast<double>(state[1]);
          std::string type_str = nb::cast<std::string>(state[2]);
          double prob0 = nb::cast<double>(state[3]);
          Eigen::VectorXd grid_points = nb::cast<Eigen::VectorXd>(state[4]);
          Eigen::VectorXd values = nb::cast<Eigen::VectorXd>(state[5]);

          // Create interpolation grid and construct object
          interp::InterpolationGrid grid(grid_points, values, 0);
          new (&kde) Kde1d(grid, xmin, xmax, type_str, prob0);

          // Note: loglik and edf are read-only calculated properties,
          // they will be recalculated based on the grid
        } else {
          // For unfitted models, construct from parameters
          double xmin = nb::cast<double>(state[0]);
          double xmax = nb::cast<double>(state[1]);
          std::string type_str = nb::cast<std::string>(state[2]);
          double multiplier = nb::cast<double>(state[3]);
          double bandwidth = nb::cast<double>(state[4]);
          size_t degree = static_cast<size_t>(nb::cast<double>(state[5]));

          new (&kde) Kde1d(xmin, xmax, type_str, multiplier, bandwidth, degree);
        }
      });
}

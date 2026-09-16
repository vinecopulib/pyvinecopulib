#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>

#include <cmath>
#include <optional>
#include <sstream>
#include <string>

namespace nb = nanobind;
using namespace nb::literals;

//! Configuration for a kernel-density fit.
//!
//! Defined here rather than in `lib/kde1d` because it is API shape for the
//! binding rather than behavior: the library's own C++ interface takes these
//! as constructor arguments, and needs no struct to do so. What it buys on
//! this side is that `Kde1d` takes controls where every other estimator in the
//! package does -- the observations, then `controls`.
//!
//! It carries the fit knobs only. A variable's type and its bounds are a
//! *declaration*, which travels keyword-only, exactly as `var_types` does on
//! `Bicop.from_data`.
struct FitControlsKde1d {
  double multiplier = 1.0;
  std::optional<double> bandwidth = std::nullopt;
  size_t degree = 2;
  size_t grid_size = 400;
  bool boundary_repair = true;

  //! The bandwidth as `Kde1d`'s constructor wants it: NaN for "choose one".
  double bandwidth_or_nan() const {
    return bandwidth.has_value() ? *bandwidth : std::nan("");
  }

  std::string str() const {
    std::stringstream ss;
    ss << "Multiplier: " << multiplier << "\n";
    ss << "Bandwidth: ";
    if (bandwidth.has_value()) {
      ss << *bandwidth << "\n";
    } else {
      ss << "selected\n";
    }
    ss << "Degree: " << degree << "\n";
    ss << "Grid size: " << grid_size << "\n";
    ss << "Boundary repair: " << (boundary_repair ? "yes" : "no") << "\n";
    return ss.str();
  }
};

inline void check_kde1d_controls(const FitControlsKde1d& controls) {
  if (!(controls.multiplier > 0.0)) {
    throw std::invalid_argument("multiplier must be positive");
  }
  if (controls.bandwidth.has_value() && !(*controls.bandwidth > 0.0)) {
    throw std::invalid_argument("bandwidth must be positive when given");
  }
  if (controls.degree > 2) {
    throw std::invalid_argument("degree must be 0, 1 or 2");
  }
  if (controls.grid_size < 4) {
    throw std::invalid_argument("grid_size must be at least 4");
  }
}

inline nb::dict kde1d_controls_to_dict(const FitControlsKde1d& controls) {
  nb::dict state;
  state["multiplier"] = controls.multiplier;
  state["bandwidth"] = controls.bandwidth;
  state["degree"] = controls.degree;
  state["grid_size"] = controls.grid_size;
  state["boundary_repair"] = controls.boundary_repair;
  return state;
}

inline void init_kde1d_fit_controls(nb::module_& module) {
  nb::class_<FitControlsKde1d>(module, "FitControlsKde1d",
                               R"(Controls for a ``Kde1d`` fit.

The knobs a kernel-density estimate is fitted with. A variable's type and its
bounds are not among them: those are a declaration about the variable, passed
keyword-only to ``fit`` / ``select`` / ``from_data``, exactly as ``var_types``
is on ``Bicop.from_data``.

Parameters
----------
multiplier : float, default=1.0
    Scales the selected bandwidth; larger values smooth more.
bandwidth : float, or None, optional
    The bandwidth to use. `None` selects one from the data.
degree : {0, 1, 2}, default=2
    Degree of the local polynomial: log-constant, log-linear or
    log-quadratic.
grid_size : int, default=400
    Number of interpolation grid points; at least 4.
boundary_repair : bool, default=True
    Whether a finite bound is fitted with a boundary correction.

Raises
------
ValueError
    If ``multiplier`` or ``bandwidth`` is not positive, ``degree`` exceeds 2,
    or ``grid_size`` is below 4.
)")
      .def(
          "__init__",
          [](FitControlsKde1d* self, double multiplier,
             std::optional<double> bandwidth, size_t degree, size_t grid_size,
             bool boundary_repair) {
            FitControlsKde1d controls{multiplier, bandwidth, degree, grid_size,
                                      boundary_repair};
            check_kde1d_controls(controls);
            new (self) FitControlsKde1d(std::move(controls));
          },
          "multiplier"_a = 1.0, "bandwidth"_a = nb::none(), "degree"_a = 2,
          "grid_size"_a = 400, "boundary_repair"_a = true)
      // Properties rather than plain fields, as `FitControlsBicop` binds
      // them: each setter re-runs the constructor's checks, so an assignment
      // is refused where the same value at construction would be. Assigning a
      // `grid_size` of 2 otherwise failed later, from inside the fit.
      .def_prop_rw(
          "multiplier",
          [](const FitControlsKde1d& self) { return self.multiplier; },
          [](FitControlsKde1d& self, double value) {
            FitControlsKde1d updated = self;
            updated.multiplier = value;
            check_kde1d_controls(updated);
            self = updated;
          },
          "Scales the selected bandwidth.")
      .def_prop_rw(
          "bandwidth",
          [](const FitControlsKde1d& self) { return self.bandwidth; },
          [](FitControlsKde1d& self, std::optional<double> value) {
            FitControlsKde1d updated = self;
            updated.bandwidth = value;
            check_kde1d_controls(updated);
            self = updated;
          },
          "The bandwidth to use, or `None` to select one.")
      .def_prop_rw(
          "degree", [](const FitControlsKde1d& self) { return self.degree; },
          [](FitControlsKde1d& self, size_t value) {
            FitControlsKde1d updated = self;
            updated.degree = value;
            check_kde1d_controls(updated);
            self = updated;
          },
          "Degree of the local polynomial.")
      .def_prop_rw(
          "grid_size",
          [](const FitControlsKde1d& self) { return self.grid_size; },
          [](FitControlsKde1d& self, size_t value) {
            FitControlsKde1d updated = self;
            updated.grid_size = value;
            check_kde1d_controls(updated);
            self = updated;
          },
          "Number of interpolation grid points.")
      .def_prop_rw(
          "boundary_repair",
          [](const FitControlsKde1d& self) { return self.boundary_repair; },
          [](FitControlsKde1d& self, bool value) {
            self.boundary_repair = value;
          },
          "Whether a finite bound gets a boundary correction.")
      .def("__repr__",
           [](const FitControlsKde1d& controls) {
             return "<pyvinecopulib.core.FitControlsKde1d>\n" + controls.str();
           })
      .def("__str__",
           [](const FitControlsKde1d& controls) {
             return "<pyvinecopulib.core.FitControlsKde1d>\n" + controls.str();
           })
      .def("__getstate__", &kde1d_controls_to_dict)
      .def("to_dict", &kde1d_controls_to_dict,
           R"(Return the settings as a plain dictionary.

Returns
-------
dict
    One entry per setting, keyed by the attribute name.

See Also
--------
pyvinecopulib.core.ControlsLike : The contract this satisfies.
)")
      .def("__setstate__", [](FitControlsKde1d& controls, nb::dict state) {
        new (&controls) FitControlsKde1d{
            nb::cast<double>(state["multiplier"]),
            nb::cast<std::optional<double>>(state["bandwidth"]),
            nb::cast<size_t>(state["degree"]),
            nb::cast<size_t>(state["grid_size"]),
            nb::cast<bool>(state["boundary_repair"])};
      });
}

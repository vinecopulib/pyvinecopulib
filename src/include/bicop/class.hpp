#pragma once

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/tuple.h>

#include <stdexcept>  // For std::invalid_argument
#include <vinecopulib.hpp>

#include "docstr.hpp"
#include "misc/helpers.hpp"

namespace nb = nanobind;
using namespace nb::literals;
using namespace vinecopulib;

// Wrapper function to call the Python bicop_plot function
inline void bicop_plot_wrapper(const Bicop& cop, const std::string& type,
                               const std::string& margin_type, nb::object xylim,
                               nb::object grid_size) {
  auto mod = nb::module_::import_("pyvinecopulib._python_helpers.bicop");
  auto bicop_plot = mod.attr("bicop_plot");
  bicop_plot(nb::cast(cop), type, margin_type, xylim, grid_size);
}

// Factory function to create a Bicop from family, rotation, parameters, and
// variable types
inline Bicop bc_from_family(const BicopFamily& family, int rotation,
                            const nb::DRef<Eigen::MatrixXd>& parameters,
                            const std::vector<std::string>& var_types = {"c",
                                                                         "c"}) {
  return Bicop(family, rotation, parameters, var_types);
}

// Factory function to create a Bicop from data, controls, and variable types
// `data` is dynamically sized, not statically two-column: a discrete variable
// needs a left-limit column, so the accepted shapes are `n x 2`, `n x (2 + k)`
// and `n x 4` (vinecopulib#729).
inline Bicop bc_from_data(const Eigen::MatrixXd& data,
                          const FitControlsBicop& controls = FitControlsBicop(),
                          const std::vector<std::string>& var_types = {"c",
                                                                       "c"}) {
  return Bicop(data, controls, var_types);
}

// Factory function to create a Bicop from a filename
inline Bicop bc_from_file(const std::string& filename) {
  return Bicop(filename);
}

// Factory function to create a Bicop from json
inline Bicop bc_from_json(const std::string& json) {
  nlohmann::json json_obj = nlohmann::json::parse(json);
  return Bicop(json_obj);
}

// Dispatch helper shared by the evaluation and score methods: pick the
// stored-parameter or per-row overload depending on whether `parameters` was
// supplied. Templated on the return type (VectorXd for pdf / hfunc / gradient,
// MatrixXd for scores / hessian / scores_cov). Argument conversion happens
// before the GIL is released, so the returned lambda touches only C++ types.
template <class R>
auto make_opt_dispatch(R (Bicop::*one)(const Eigen::MatrixXd&) const,
                       R (Bicop::*many)(const Eigen::MatrixXd&,
                                        const Eigen::MatrixXd&, size_t) const) {
  return [one, many](const Bicop& self, const Eigen::MatrixXd& u,
                     const std::optional<Eigen::MatrixXd>& parameters,
                     size_t num_threads) -> R {
    if (parameters) return (self.*many)(u, *parameters, num_threads);
    return (self.*one)(u);
  };
}

inline void init_bicop_class(nb::module_& module) {
  constexpr auto& bicop_doc = pyvinecopulib_doc.vinecopulib.Bicop;

  const char* default_constructor_doc =
      R"""(Default constructor for the ``Bicop`` class.

The default constructor uses ``Bicop.from_family()`` to instantiate an
independent bivariate copula. It can then be used to select a model from data using ``Bicop.select()``. Or if a ``BicopFamily`` is passed to the constructor, then the method ``Bicop.fit()`` can be used to fit the copula to data.

Alternatives to instantiate bivariate copulas are:

- ``Bicop.from_family()``: Instantiate from a family, rotation, parameters, and variable types.
- ``Bicop.from_data()``: Instantiate from data, as well as optional controls and variable types.
- ``Bicop.from_file()``: Instantiate from a JSON file, or a CBOR file when the filename ends in ``.cbor``.
- ``Bicop.from_json()``: Instantiate from a JSON string.
)""";

  // `simulate` takes either a sample size or one parameter set per drawn
  // observation (vinecopulib#719), which fixes the sample size by its row
  // count. `parameters` is keyword-only so that the long-standing positional
  // `sample(n, qrng, seeds)` keeps its meaning.
  const std::string sample_doc = std::string(bicop_doc.simulate.doc_3args) +
                                 R"""(

Notes
-----
Pass ``parameters`` instead of ``n`` to draw each observation from a different
parameter set: an ``(n, p)`` array with ``p == len(self.parameters)`` columns in
the family's natural (unrotated) parameterization, drawing one observation per
row and parallelized over ``num_threads``. Exactly one of ``n`` and
``parameters`` may be given. ``parameters`` must be two-dimensional, in either
memory order. This is supported for parametric families only; a nonparametric
family, a wrong shape, or non-finite or out-of-bounds values raise
``RuntimeError``. The draws are continuous even when ``var_types`` marks a
variable discrete. ``num_threads`` applies only to this form; drawing ``n``
observations from the stored parameters is single-threaded.
)""";

  // `pdf` / `cdf` / `hfunc*` / `hinv*` / `loglik` each expose an optional
  // `parameters` argument (vinecopulib#675): when omitted the copula's stored
  // parameters are used (the original behavior); when given an (n, p) matrix
  // they are evaluated with a different parameter set per row of `u`. They are
  // bound as a single Python method (rather than two overloads) so the
  // generated `.pyi` stub carries the full signature. The base docstring is the
  // libclang-extracted one for the stored-parameter form, extended with a note
  // describing the per-row form.
  const std::string per_row_note =
      R"""(

Notes
-----
If ``parameters`` is given, the copula is evaluated with a different parameter
set per row of ``u`` instead of the stored parameters. ``parameters`` is then an
``(n, p)`` array with one row per row of ``u`` and ``p == len(self.parameters)``
columns, in the family's natural (unrotated) parameterization, and the
evaluation may be parallelized over ``num_threads``. This is supported for
parametric families only; a nonparametric family, a wrong shape, or non-finite
or out-of-bounds values raise ``RuntimeError``.
)""";
  const std::string pdf_doc =
      std::string(bicop_doc.pdf.doc_1args) + per_row_note;
  const std::string cdf_doc =
      std::string(bicop_doc.cdf.doc_1args) + per_row_note;
  const std::string hfunc1_doc =
      std::string(bicop_doc.hfunc1.doc_1args) + per_row_note;
  const std::string hfunc2_doc =
      std::string(bicop_doc.hfunc2.doc_1args) + per_row_note;
  const std::string hinv1_doc =
      std::string(bicop_doc.hinv1.doc_1args) + per_row_note;
  const std::string hinv2_doc =
      std::string(bicop_doc.hinv2.doc_1args) + per_row_note;
  const std::string loglik_doc =
      std::string(bicop_doc.loglik.doc_1args) + per_row_note;

  // The eight derivative methods (vinecopulib#683/#687/#694) follow the same
  // one-Python-method-per-name design: the 2-argument (stored-parameter) and
  // 4-argument (per-row-parameter) C++ overloads are dispatched on whether
  // `parameters` was supplied.
  const std::string pdf_deriv_doc =
      std::string(bicop_doc.pdf_deriv.doc_2args) + per_row_note;
  const std::string pdf_deriv2_doc =
      std::string(bicop_doc.pdf_deriv2.doc_2args) + per_row_note;
  const std::string hfunc1_deriv_doc =
      std::string(bicop_doc.hfunc1_deriv.doc_2args) + per_row_note;
  const std::string hfunc1_deriv2_doc =
      std::string(bicop_doc.hfunc1_deriv2.doc_2args) + per_row_note;
  const std::string hfunc2_deriv_doc =
      std::string(bicop_doc.hfunc2_deriv.doc_2args) + per_row_note;
  const std::string hfunc2_deriv2_doc =
      std::string(bicop_doc.hfunc2_deriv2.doc_2args) + per_row_note;
  const std::string logpdf_deriv_doc =
      std::string(bicop_doc.logpdf_deriv.doc_2args) + per_row_note;
  const std::string logpdf_deriv2_doc =
      std::string(bicop_doc.logpdf_deriv2.doc_2args) + per_row_note;

  // The log-likelihood score family (vinecopulib#699) also takes the optional
  // per-row `parameters`, so it shares `per_row_note`. `scores_full` returns a
  // dict on the Python side, so its docstring is hand-written (the generated
  // one describes the C++ `ScoresResult`).
  const std::string scores_doc =
      std::string(bicop_doc.scores.doc_1args) + per_row_note;
  const std::string gradient_doc =
      std::string(bicop_doc.gradient.doc_1args) + per_row_note;
  const std::string hessian_doc =
      std::string(bicop_doc.hessian.doc_1args) + per_row_note;
  const std::string hessian_full_doc =
      std::string(bicop_doc.hessian_full.doc_1args) + per_row_note;
  const std::string scores_cov_doc =
      std::string(bicop_doc.scores_cov.doc_1args) + per_row_note;
  const std::string scores_full_doc =
      std::string(
          R"""(Evaluates the log-likelihood scores, bundled in a dict.

Provided for parity with ``Vinecop.scores_full``; a single pair copula has no
cascade caches, so the returned dict carries only the score matrix.

Returns
-------
dict
    A dict with a single key ``"scores"`` -- the per-observation score matrix
    (shape ``(n, p)``), identical to ``Bicop.scores``.
)""") +
      per_row_note;

  // `flip` binds the libclang-extracted doc, extended with the return
  // description (the Python method returns a flipped copy; the underlying C++
  // method mutates in place, which users need not know about).
  const std::string flip_doc = std::string(bicop_doc.flip.doc) +
                               R"""(

The flipped copula satisfies ``c'(u1, u2) = c(u2, u1)``, with the two
h-functions (and their inverses) exchanged.

Returns
-------
Bicop
    The argument-swapped copula; the object itself is left unchanged.
)""";

  // Dispatch for the derivative methods, which additionally thread the `deriv`
  // selector string through both overloads (the plain evaluation / score
  // dispatch lives in the namespace-scope `make_opt_dispatch` above).
  auto deriv_with_optional_params =
      [](Eigen::VectorXd (Bicop::*one)(const Eigen::MatrixXd&,
                                       const std::string&) const,
         Eigen::VectorXd (Bicop::*many)(const Eigen::MatrixXd&,
                                        const std::string&,
                                        const Eigen::MatrixXd&, size_t) const) {
        return [one, many](const Bicop& self, const Eigen::MatrixXd& u,
                           const std::string& deriv,
                           const std::optional<Eigen::MatrixXd>& parameters,
                           size_t num_threads) -> Eigen::VectorXd {
          if (parameters)
            return (self.*many)(u, deriv, *parameters, num_threads);
          return (self.*one)(u, deriv);
        };
      };

  nb::class_<Bicop>(module, "Bicop", bicop_doc.doc)
      // Default constructor
      .def(nb::init<const BicopFamily, const int,
                    const nb::DRef<Eigen::MatrixXd>&,
                    const std::vector<std::string>&>(),
           "family"_a = BicopFamily::indep, "rotation"_a = 0,
           "parameters"_a = Eigen::MatrixXd(),
           "var_types"_a = std::vector<std::string>(2, "c"),
           default_constructor_doc, nb::call_guard<nb::gil_scoped_release>())
      .def_static("from_family", &bc_from_family,
                  "family"_a = BicopFamily::indep, "rotation"_a = 0,
                  "parameters"_a = Eigen::MatrixXd(),
                  "var_types"_a = std::vector<std::string>(2, "c"),
                  bicop_doc.ctor.doc_4args_family_rotation_parameters_var_types,
                  nb::call_guard<nb::gil_scoped_release>())
      .def_static("from_data", &bc_from_data, "data"_a,
                  "controls"_a.sig("FitControlsBicop()") = FitControlsBicop(),
                  "var_types"_a = std::vector<std::string>(2, "c"),
                  bicop_doc.ctor.doc_3args_data_controls_var_types,
                  nb::call_guard<nb::gil_scoped_release>())
      .def_static("from_file", &bc_from_file, "filename"_a,
                  bicop_doc.ctor.doc_1args_filename,
                  nb::call_guard<nb::gil_scoped_release>())
      .def_static("from_json", &bc_from_json, "json"_a,
                  bicop_doc.ctor.doc_1args_input,
                  nb::call_guard<nb::gil_scoped_release>())
      .def("to_file", &Bicop::to_file, "filename"_a, bicop_doc.to_file.doc,
           nb::call_guard<nb::gil_scoped_release>())
      .def(
          "to_json",
          [](Bicop& self) -> std::string { return self.to_json().dump(); },
          bicop_doc.to_json.doc, nb::call_guard<nb::gil_scoped_release>())
      .def_prop_rw("rotation", &Bicop::get_rotation, &Bicop::set_rotation,
                   bicop_doc.get_rotation.doc)
      .def_prop_rw(
          "parameters", &Bicop::get_parameters,
          [](Bicop& self, const nb::DRef<Eigen::MatrixXd>& parameters) {
            self.set_parameters(parameters);
          },
          bicop_doc.get_parameters.doc,
          nb::call_guard<nb::gil_scoped_release>())
      .def_prop_rw("var_types", &Bicop::get_var_types, &Bicop::set_var_types,
                   bicop_doc.get_var_types.doc)
      .def_prop_ro("family", &Bicop::get_family, bicop_doc.get_family.doc)
      .def_prop_ro("family_name", &Bicop::get_family_name,
                   bicop_doc.get_family_name.doc)
      .def_prop_ro("tau", &Bicop::get_tau, bicop_doc.get_tau.doc)
      .def_prop_ro("taildep", &Bicop::get_taildep, bicop_doc.get_taildep.doc,
                   nb::call_guard<nb::gil_scoped_release>())
      .def_prop_ro("beta", &Bicop::get_beta, bicop_doc.get_beta.doc,
                   nb::call_guard<nb::gil_scoped_release>())
      .def_prop_ro("npars", &Bicop::get_npars, bicop_doc.get_npars.doc)
      .def(
          "loglik",
          [](const Bicop& self, const Eigen::MatrixXd& u,
             const std::optional<Eigen::MatrixXd>& parameters,
             size_t num_threads) -> double {
            if (parameters) return self.loglik(u, *parameters, num_threads);
            return self.loglik(u);
          },
          "u"_a = Eigen::MatrixXd(), "parameters"_a = nb::none(),
          "num_threads"_a = static_cast<size_t>(1), loglik_doc.c_str(),
          nb::call_guard<nb::gil_scoped_release>())
      .def_prop_ro("nobs", &Bicop::get_nobs, bicop_doc.get_nobs.doc)
      .def("aic", &Bicop::aic, "u"_a = Eigen::MatrixXd(), bicop_doc.aic.doc,
           nb::call_guard<nb::gil_scoped_release>())
      .def("bic", &Bicop::bic, "u"_a = Eigen::MatrixXd(), bicop_doc.bic.doc,
           nb::call_guard<nb::gil_scoped_release>())
      .def("mbic", &Bicop::mbic, "u"_a = Eigen::MatrixXd(), "psi0"_a = 0.9,
           bicop_doc.mbic.doc, nb::call_guard<nb::gil_scoped_release>())
      .def(
          "__repr__",
          [](const Bicop& cop) {
            return "<pyvinecopulib.core.Bicop> " + cop.str();
          },
          bicop_doc.str.doc)
      .def(
          "__str__",
          [](const Bicop& cop) {
            return "<pyvinecopulib.core.Bicop> " + cop.str();
          },
          bicop_doc.str.doc)
      .def("parameters_to_tau", &Bicop::parameters_to_tau, "parameters"_a,
           bicop_doc.parameters_to_tau.doc,
           nb::call_guard<nb::gil_scoped_release>())
      .def("tau_to_parameters", &Bicop::tau_to_parameters, "tau"_a,
           bicop_doc.tau_to_parameters.doc,
           nb::call_guard<nb::gil_scoped_release>())
      .def("parameters_to_taildep", &Bicop::parameters_to_taildep,
           "parameters"_a, bicop_doc.parameters_to_taildep.doc,
           nb::call_guard<nb::gil_scoped_release>())
      .def("parameters_to_beta", &Bicop::parameters_to_beta, "parameters"_a,
           bicop_doc.parameters_to_beta.doc,
           nb::call_guard<nb::gil_scoped_release>())
      .def_prop_ro("parameters_lower_bounds",
                   &Bicop::get_parameters_lower_bounds,
                   bicop_doc.get_parameters_lower_bounds.doc,
                   nb::call_guard<nb::gil_scoped_release>())
      .def_prop_ro("parameters_upper_bounds",
                   &Bicop::get_parameters_upper_bounds,
                   bicop_doc.get_parameters_upper_bounds.doc,
                   nb::call_guard<nb::gil_scoped_release>())
      .def("pdf", make_opt_dispatch(&Bicop::pdf, &Bicop::pdf), "u"_a,
           "parameters"_a = nb::none(),
           "num_threads"_a = static_cast<size_t>(1), pdf_doc.c_str(),
           nb::call_guard<nb::gil_scoped_release>())
      .def("cdf", make_opt_dispatch(&Bicop::cdf, &Bicop::cdf), "u"_a,
           "parameters"_a = nb::none(),
           "num_threads"_a = static_cast<size_t>(1), cdf_doc.c_str(),
           nb::call_guard<nb::gil_scoped_release>())
      .def("hfunc1", make_opt_dispatch(&Bicop::hfunc1, &Bicop::hfunc1), "u"_a,
           "parameters"_a = nb::none(),
           "num_threads"_a = static_cast<size_t>(1), hfunc1_doc.c_str(),
           nb::call_guard<nb::gil_scoped_release>())
      .def("hfunc2", make_opt_dispatch(&Bicop::hfunc2, &Bicop::hfunc2), "u"_a,
           "parameters"_a = nb::none(),
           "num_threads"_a = static_cast<size_t>(1), hfunc2_doc.c_str(),
           nb::call_guard<nb::gil_scoped_release>())
      .def("hinv1", make_opt_dispatch(&Bicop::hinv1, &Bicop::hinv1), "u"_a,
           "parameters"_a = nb::none(),
           "num_threads"_a = static_cast<size_t>(1), hinv1_doc.c_str(),
           nb::call_guard<nb::gil_scoped_release>())
      .def("hinv2", make_opt_dispatch(&Bicop::hinv2, &Bicop::hinv2), "u"_a,
           "parameters"_a = nb::none(),
           "num_threads"_a = static_cast<size_t>(1), hinv2_doc.c_str(),
           nb::call_guard<nb::gil_scoped_release>())
      .def("pdf_deriv",
           deriv_with_optional_params(&Bicop::pdf_deriv, &Bicop::pdf_deriv),
           "u"_a, "deriv"_a, "parameters"_a = nb::none(),
           "num_threads"_a = static_cast<size_t>(1), pdf_deriv_doc.c_str(),
           nb::call_guard<nb::gil_scoped_release>())
      .def("pdf_deriv2",
           deriv_with_optional_params(&Bicop::pdf_deriv2, &Bicop::pdf_deriv2),
           "u"_a, "deriv"_a, "parameters"_a = nb::none(),
           "num_threads"_a = static_cast<size_t>(1), pdf_deriv2_doc.c_str(),
           nb::call_guard<nb::gil_scoped_release>())
      .def("hfunc1_deriv",
           deriv_with_optional_params(&Bicop::hfunc1_deriv,
                                      &Bicop::hfunc1_deriv),
           "u"_a, "deriv"_a, "parameters"_a = nb::none(),
           "num_threads"_a = static_cast<size_t>(1), hfunc1_deriv_doc.c_str(),
           nb::call_guard<nb::gil_scoped_release>())
      .def("hfunc1_deriv2",
           deriv_with_optional_params(&Bicop::hfunc1_deriv2,
                                      &Bicop::hfunc1_deriv2),
           "u"_a, "deriv"_a, "parameters"_a = nb::none(),
           "num_threads"_a = static_cast<size_t>(1), hfunc1_deriv2_doc.c_str(),
           nb::call_guard<nb::gil_scoped_release>())
      .def("hfunc2_deriv",
           deriv_with_optional_params(&Bicop::hfunc2_deriv,
                                      &Bicop::hfunc2_deriv),
           "u"_a, "deriv"_a, "parameters"_a = nb::none(),
           "num_threads"_a = static_cast<size_t>(1), hfunc2_deriv_doc.c_str(),
           nb::call_guard<nb::gil_scoped_release>())
      .def("hfunc2_deriv2",
           deriv_with_optional_params(&Bicop::hfunc2_deriv2,
                                      &Bicop::hfunc2_deriv2),
           "u"_a, "deriv"_a, "parameters"_a = nb::none(),
           "num_threads"_a = static_cast<size_t>(1), hfunc2_deriv2_doc.c_str(),
           nb::call_guard<nb::gil_scoped_release>())
      .def("logpdf_deriv",
           deriv_with_optional_params(&Bicop::logpdf_deriv,
                                      &Bicop::logpdf_deriv),
           "u"_a, "deriv"_a, "parameters"_a = nb::none(),
           "num_threads"_a = static_cast<size_t>(1), logpdf_deriv_doc.c_str(),
           nb::call_guard<nb::gil_scoped_release>())
      .def("logpdf_deriv2",
           deriv_with_optional_params(&Bicop::logpdf_deriv2,
                                      &Bicop::logpdf_deriv2),
           "u"_a, "deriv"_a, "parameters"_a = nb::none(),
           "num_threads"_a = static_cast<size_t>(1), logpdf_deriv2_doc.c_str(),
           nb::call_guard<nb::gil_scoped_release>())
      .def("scores", make_opt_dispatch(&Bicop::scores, &Bicop::scores), "u"_a,
           "parameters"_a = nb::none(),
           "num_threads"_a = static_cast<size_t>(1), scores_doc.c_str(),
           nb::call_guard<nb::gil_scoped_release>())
      .def("gradient", make_opt_dispatch(&Bicop::gradient, &Bicop::gradient),
           "u"_a, "parameters"_a = nb::none(),
           "num_threads"_a = static_cast<size_t>(1), gradient_doc.c_str(),
           nb::call_guard<nb::gil_scoped_release>())
      .def("hessian", make_opt_dispatch(&Bicop::hessian, &Bicop::hessian),
           "u"_a, "parameters"_a = nb::none(),
           "num_threads"_a = static_cast<size_t>(1), hessian_doc.c_str(),
           nb::call_guard<nb::gil_scoped_release>())
      .def("scores_cov",
           make_opt_dispatch(&Bicop::scores_cov, &Bicop::scores_cov), "u"_a,
           "parameters"_a = nb::none(),
           "num_threads"_a = static_cast<size_t>(1), scores_cov_doc.c_str(),
           nb::call_guard<nb::gil_scoped_release>())
      .def(
          "hessian_full",
          [](const Bicop& self, const Eigen::MatrixXd& u,
             const std::optional<Eigen::MatrixXd>& parameters,
             size_t num_threads) -> nb::list {
            std::vector<Eigen::MatrixXd> hess;
            {
              // Release the GIL only around the C++ computation; the list
              // construction below must hold it.
              nb::gil_scoped_release release;
              hess = parameters ? self.hessian_full(u, *parameters, num_threads)
                                : self.hessian_full(u);
            }
            nb::list out;
            for (const auto& h : hess) out.append(nb::cast(h));
            return out;
          },
          "u"_a, "parameters"_a = nb::none(),
          "num_threads"_a = static_cast<size_t>(1), hessian_full_doc.c_str())
      .def(
          "scores_full",
          [](const Bicop& self, const Eigen::MatrixXd& u,
             const std::optional<Eigen::MatrixXd>& parameters,
             size_t num_threads) -> nb::dict {
            Bicop::ScoresResult res;
            {
              // Release the GIL only around the C++ computation; the dict
              // construction below must hold it.
              nb::gil_scoped_release release;
              res = parameters ? self.scores_full(u, *parameters, num_threads)
                               : self.scores_full(u);
            }
            nb::dict out;
            out["scores"] = nb::cast(res.scores);
            return out;
          },
          "u"_a, "parameters"_a = nb::none(),
          "num_threads"_a = static_cast<size_t>(1), scores_full_doc.c_str())
      .def(
          "sample",
          [](const Bicop& self, std::optional<size_t> n, bool qrng,
             const std::vector<int>& seeds,
             const std::optional<Eigen::MatrixXd>& parameters,
             size_t num_threads) -> Eigen::MatrixXd {
            if (n.has_value() == parameters.has_value()) {
              throw std::invalid_argument(
                  "give exactly one of `n` and `parameters`; `parameters` "
                  "already fixes the number of observations by its row "
                  "count");
            }
            nb::gil_scoped_release release;
            if (parameters) {
              return self.simulate(*parameters, qrng, seeds, num_threads);
            }
            return self.simulate(*n, qrng, seeds);
          },
          "n"_a = nb::none(), "qrng"_a = false, "seeds"_a = std::vector<int>(),
          nb::kw_only(), "parameters"_a = nb::none(),
          "num_threads"_a = static_cast<size_t>(1), sample_doc.c_str())
      .def(
          "flip",
          [](const Bicop& cop) {
            Bicop flipped = cop;  // deep copy (Bicop's copy ctor clones)
            flipped.flip();
            return flipped;
          },
          flip_doc.c_str(), nb::call_guard<nb::gil_scoped_release>())
      .def("as_continuous", &Bicop::as_continuous, bicop_doc.as_continuous.doc,
           nb::call_guard<nb::gil_scoped_release>())
      .def("fit", &Bicop::fit, "data"_a,
           "controls"_a.sig("FitControlsBicop()") = FitControlsBicop(),
           bicop_doc.fit.doc, nb::call_guard<nb::gil_scoped_release>())
      .def("select", &Bicop::select, "data"_a,
           "controls"_a.sig("FitControlsBicop()") = FitControlsBicop(),
           bicop_doc.select.doc, nb::call_guard<nb::gil_scoped_release>())
      .def("plot", &bicop_plot_wrapper, "type"_a = "surface",
           "margin_type"_a = "unif", "xylim"_a = nb::none(),
           "grid_size"_a = nb::none(),
           python_doc_helper(
               "pyvinecopulib._python_helpers.bicop", "BICOP_PLOT_DOC",
               "Plot the bivariate copula (extended doc unavailable) ")
               .c_str())
      .def("__getstate__",
           [](const Bicop& cop) {
             nb::dict state;
             state["family"] = cop.get_family();
             state["rotation"] = cop.get_rotation();
             state["parameters"] = cop.get_parameters();
             state["var_types"] = cop.get_var_types();
             return state;
           })

      .def("__setstate__", [](Bicop& cop, nb::dict state) {
        new (&cop)
            Bicop(nb::cast<BicopFamily>(state["family"]),
                  nb::cast<int>(state["rotation"]),
                  nb::cast<Eigen::MatrixXd>(state["parameters"]),
                  nb::cast<std::vector<std::string>>(state["var_types"]));
      });
}

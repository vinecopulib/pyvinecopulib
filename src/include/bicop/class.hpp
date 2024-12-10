#pragma once

#include "docstr.hpp"
#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <vinecopulib.hpp>

namespace nb = nanobind;
using namespace nb::literals;
using namespace vinecopulib;

// Factory function to create a Bicop from family, rotation, parameters, and
// variable types
inline Bicop
bc_from_family(const BicopFamily& family,
               int rotation,
               const nb::DRef<Eigen::MatrixXd>& parameters,
               const std::vector<std::string>& var_types = { "c", "c" })
{
  return Bicop(family, rotation, parameters, var_types);
}

// Factory function to create a Bicop from data, controls, and variable types
inline Bicop
bc_from_data(const Eigen::Matrix<double, Eigen::Dynamic, 2>& data,
             const FitControlsBicop& controls = FitControlsBicop(),
             const std::vector<std::string>& var_types = { "c", "c" })
{
  return Bicop(data, controls, var_types);
}

// Factory function to create a Bicop from a filename
inline Bicop
bc_from_file(const std::string& filename)
{
  return Bicop(filename);
}

inline void
init_bicop_class(nb::module_& module)
{

  constexpr auto& bicop_doc = pyvinecopulib_doc.vinecopulib.Bicop;

  const char* default_constructor_doc =
    R"""(Default constructor for the ``Bicop`` class.

The default constructor uses ``Bicop.from_family()`` to instantiate an
independent bivariate copula. It can then be used to select a model from data using ``Bicop.select()``. Or if a ``BicopFamily`` is passed to the constructor, then the method ``Bicop.fit()`` can be used to fit the copula to data.
To instantiate directly from data or from a file, use ``Bicop.from_data()``
and ``Bicop.from_file()`` respectively.)""";

  nb::class_<Bicop>(module, "Bicop", bicop_doc.doc)
    .def(nb::init<const BicopFamily,
                  const int,
                  const nb::DRef<Eigen::MatrixXd>&,
                  const std::vector<std::string>&>(),
         "family"_a = BicopFamily::indep,
         "rotation"_a = 0,
         "parameters"_a = Eigen::MatrixXd(),
         "var_types"_a = std::vector<std::string>(2, "c"),
         default_constructor_doc) // Default constructor
    .def_static("from_family",
                &bc_from_family,
                "family"_a = BicopFamily::indep,
                "rotation"_a = 0,
                "parameters"_a = Eigen::MatrixXd(),
                "var_types"_a = std::vector<std::string>(2, "c"),
                bicop_doc.ctor.doc_4args_family_rotation_parameters_var_types)
    .def_static("from_data",
                &bc_from_data,
                "data"_a,
                "controls"_a.sig("FitControlsBicop()") = FitControlsBicop(),
                "var_types"_a = std::vector<std::string>(2, "c"),
                bicop_doc.ctor.doc_3args_data_controls_var_types)
    .def_static("from_file",
                &bc_from_file,
                "filename"_a,
                bicop_doc.ctor.doc_1args_filename)
    .def("to_file", &Bicop::to_file, "filename"_a, bicop_doc.to_file.doc)
    .def("to_json", &Bicop::to_json, bicop_doc.to_json.doc)
    .def_prop_rw("rotation",
                 &Bicop::get_rotation,
                 &Bicop::set_rotation,
                 "The copula rotation.")
    .def_prop_rw(
      "parameters",
      &Bicop::get_parameters,
      [](Bicop& self, const nb::DRef<Eigen::MatrixXd>& parameters) {
        self.set_parameters(parameters);
      },
      "The copula parameter(s).")
    .def_prop_rw("var_types",
                 &Bicop::get_var_types,
                 &Bicop::set_var_types,
                 "The type of the two variables.")
    .def_prop_ro("family", &Bicop::get_family, "The copula family.")
    .def_prop_ro("tau", &Bicop::get_tau, "The Kendall's tau.")
    .def_prop_ro("npars",
                 &Bicop::get_npars,
                 "The number of parameters (for nonparametric "
                 "families, a conceptually similar definition).")
    .def("loglik",
         &Bicop::loglik,
         "u"_a = Eigen::Matrix<double, Eigen::Dynamic, 2>(),
         bicop_doc.loglik.doc)
    .def_prop_ro("nobs",
                 &Bicop::get_nobs,
                 "The number of observations (only for fitted objects).")
    .def("aic",
         &Bicop::aic,
         "u"_a = Eigen::Matrix<double, Eigen::Dynamic, 2>(),
         bicop_doc.aic.doc)
    .def("bic",
         &Bicop::bic,
         "u"_a = Eigen::Matrix<double, Eigen::Dynamic, 2>(),
         bicop_doc.bic.doc)
    .def("mbic",
         &Bicop::mbic,
         "u"_a = Eigen::Matrix<double, Eigen::Dynamic, 2>(),
         "psi0"_a = 0.9,
         bicop_doc.mbic.doc)
    .def("__repr__",
         [](const Bicop& cop) { return "<pyvinecopulib.Bicop>\n" + cop.str(); })
    .def("str", &Bicop::str, bicop_doc.str.doc)
    .def("parameters_to_tau",
         &Bicop::parameters_to_tau,
         "parameters"_a,
         bicop_doc.parameters_to_tau.doc)
    .def("tau_to_parameters",
         &Bicop::tau_to_parameters,
         "tau"_a,
         bicop_doc.tau_to_parameters.doc)
    .def_prop_ro("parameters_lower_bounds",
                 &Bicop::get_parameters_lower_bounds,
                 bicop_doc.get_parameters_lower_bounds.doc)
    .def_prop_ro("parameters_upper_bounds",
                 &Bicop::get_parameters_upper_bounds,
                 bicop_doc.get_parameters_upper_bounds.doc)
    .def("pdf", &Bicop::pdf, "u"_a, bicop_doc.pdf.doc)
    .def("cdf", &Bicop::cdf, "u"_a, bicop_doc.cdf.doc)
    .def("hfunc1", &Bicop::hfunc1, "u"_a, bicop_doc.hfunc1.doc)
    .def("hfunc2", &Bicop::hfunc2, "u"_a, bicop_doc.hfunc2.doc)
    .def("hinv1", &Bicop::hinv1, "u"_a, bicop_doc.hinv1.doc)
    .def("hinv2", &Bicop::hinv2, "u"_a, bicop_doc.hinv2.doc)
    .def("simulate",
         &Bicop::simulate,
         "n"_a,
         "qrng"_a = false,
         "seeds"_a = std::vector<int>(),
         bicop_doc.simulate.doc)
    .def("fit",
         &Bicop::fit,
         "data"_a,
         "controls"_a.sig("FitControlsBicop()") = FitControlsBicop(),
         bicop_doc.fit.doc)
    .def("select",
         &Bicop::select,
         "data"_a,
         "controls"_a.sig("FitControlsBicop()") = FitControlsBicop(),
         bicop_doc.select.doc)
    .def(
      "plot",
      [](const Bicop& cop,
         const std::string& type = "surface",
         const std::string& margin_type = "unif",
         nb::object xylim = nb::none(),
         nb::object grid_size = nb::none()) {
        auto python_helpers_plotting =
          nb::module_::import_("pyvinecopulib._python_helpers.bicop");

        // Import the Python plotting function
        nb::object bicop_plot = python_helpers_plotting.attr("bicop_plot");

        // Call the Python function with the provided arguments
        bicop_plot(nb::cast(cop), type, margin_type, xylim, grid_size);
      },
      "type"_a = "surface",
      "margin_type"_a = "unif",
      "xylim"_a = nb::none(),
      "grid_size"_a = nb::none(),
      nb::cast<std::string>(
        nb::module_::import_("pyvinecopulib._python_helpers.bicop")
          .attr("BICOP_PLOT_DOC"))
        .c_str());
}
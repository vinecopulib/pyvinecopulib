#pragma once

#include "docstr.hpp"
#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
// #include <nanobind/stl.h>
#include <vinecopulib.hpp>

namespace nb = nanobind;
using namespace vinecopulib;

inline void
init_bicop_class(nb::module_& module)
{

  constexpr auto& bicop_doc = pyvinecopulib_doc.vinecopulib.Bicop;

  nb::class_<Bicop>(module, "Bicop", bicop_doc.doc)
    .def(nb::init<const BicopFamily,
                  const int,
                  const Eigen::MatrixXd&,
                  const std::vector<std::string>&>(),
         nb::arg("family") = BicopFamily::indep,
         nb::arg("rotation") = 0,
         nb::arg("parameters") = Eigen::MatrixXd(),
         nb::arg("var_types") = std::vector<std::string>(2, "c"),
         bicop_doc.ctor.doc_4args_family_rotation_parameters_var_types)
    .def(nb::init<const Eigen::Matrix<double, Eigen::Dynamic, 2>&,
                  const FitControlsBicop&,
                  const std::vector<std::string>&>(),
         nb::arg("data"),
         nb::arg("controls") = FitControlsBicop(),
         nb::arg("var_types") = std::vector<std::string>(2, "c"),
         bicop_doc.ctor.doc_3args_data_controls_var_types)
    .def(nb::init<const std::string>(),
         nb::arg("filename"),
         bicop_doc.ctor.doc_1args_filename)
    .def("to_json", &Bicop::to_file, nb::arg("filename"), bicop_doc.to_file.doc)
    .def_prop_rw("rotation",
                 &Bicop::get_rotation,
                 &Bicop::set_rotation,
                 "The copula rotation.")
    .def_prop_rw("parameters",
                 &Bicop::get_parameters,
                 &Bicop::set_parameters,
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
         nb::arg("u") = Eigen::Matrix<double, Eigen::Dynamic, 2>(),
         bicop_doc.loglik.doc)
    .def_prop_ro("nobs",
                 &Bicop::get_nobs,
                 "The number of observations (only for fitted objects).")
    .def("aic",
         &Bicop::aic,
         nb::arg("u") = Eigen::Matrix<double, Eigen::Dynamic, 2>(),
         bicop_doc.aic.doc)
    .def("bic",
         &Bicop::bic,
         nb::arg("u") = Eigen::Matrix<double, Eigen::Dynamic, 2>(),
         bicop_doc.bic.doc)
    .def("mbic",
         &Bicop::mbic,
         nb::arg("u") = Eigen::Matrix<double, Eigen::Dynamic, 2>(),
         nb::arg("psi0") = 0.9,
         bicop_doc.mbic.doc)
    .def("__repr__",
         [](const Bicop& cop) { return "<pyvinecopulib.Bicop>\n" + cop.str(); })
    .def("str", &Bicop::str, bicop_doc.str.doc)
    .def("parameters_to_tau",
         &Bicop::parameters_to_tau,
         nb::arg("parameters"),
         bicop_doc.parameters_to_tau.doc)
    .def("tau_to_parameters",
         &Bicop::tau_to_parameters,
         nb::arg("tau"),
         bicop_doc.tau_to_parameters.doc)
    .def("parameters_lower_bounds",
         &Bicop::get_parameters_lower_bounds,
         bicop_doc.get_parameters_lower_bounds.doc)
    .def("parameters_upper_bounds",
         &Bicop::get_parameters_upper_bounds,
         bicop_doc.get_parameters_upper_bounds.doc)
    .def("pdf", &Bicop::pdf, nb::arg("u"), bicop_doc.pdf.doc)
    .def("cdf", &Bicop::cdf, nb::arg("u"), bicop_doc.cdf.doc)
    .def("hfunc1", &Bicop::hfunc1, nb::arg("u"), bicop_doc.hfunc1.doc)
    .def("hfunc2", &Bicop::hfunc2, nb::arg("u"), bicop_doc.hfunc2.doc)
    .def("hinv1", &Bicop::hinv1, nb::arg("u"), bicop_doc.hinv1.doc)
    .def("hinv2", &Bicop::hinv2, nb::arg("u"), bicop_doc.hinv2.doc)
    .def("simulate",
         &Bicop::simulate,
         nb::arg("n"),
         nb::arg("qrng") = false,
         nb::arg("seeds") = std::vector<int>(),
         bicop_doc.simulate.doc)
    .def("fit",
         &Bicop::fit,
         nb::arg("data"),
         nb::arg("controls") = FitControlsBicop(),
         bicop_doc.fit.doc)
    .def("select",
         &Bicop::select,
         nb::arg("data"),
         nb::arg("controls") = FitControlsBicop(),
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
      nb::arg("type") = "surface",
      nb::arg("margin_type") = "unif",
      nb::arg("xylim") = nb::none(),
      nb::arg("grid_size") = nb::none(),
      nb::cast<std::string>(
        nb::module_::import_("pyvinecopulib._python_helpers.bicop")
          .attr("BICOP_PLOT_DOC"))
        .c_str());
}
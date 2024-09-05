#pragma once

#include "docstr.hpp"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vinecopulib.hpp>

namespace py = pybind11;
using namespace vinecopulib;

inline void
init_bicop_class(py::module_& module)
{

  constexpr auto& doc = pyvinecopulib_doc;
  constexpr auto& bicop_doc = doc.vinecopulib.Bicop;

  py::class_<Bicop>(module, "Bicop", bicop_doc.doc)
    .def(py::init<const BicopFamily,
                  const int,
                  const Eigen::MatrixXd&,
                  const std::vector<std::string>&>(),
         py::arg("family") = BicopFamily::indep,
         py::arg("rotation") = 0,
         py::arg("parameters") = Eigen::MatrixXd(),
         py::arg("var_types") = std::vector<std::string>(2, "c"),
         bicop_doc.ctor.doc_4args_family_rotation_parameters_var_types)
    .def(py::init<const Eigen::Matrix<double, Eigen::Dynamic, 2>&,
                  const FitControlsBicop&,
                  const std::vector<std::string>&>(),
         py::arg("data"),
         py::arg_v("controls", FitControlsBicop(), "FitControlsBicop()"),
         py::arg("var_types") = std::vector<std::string>(2, "c"),
         bicop_doc.ctor.doc_3args_data_controls_var_types)
    .def(py::init<const std::string>(),
         py::arg("filename"),
         bicop_doc.ctor.doc_1args_filename)
    .def("to_json", &Bicop::to_file, py::arg("filename"), bicop_doc.to_file.doc)
    .def_property("rotation",
                  &Bicop::get_rotation,
                  &Bicop::set_rotation,
                  "The copula rotation.")
    .def_property("parameters",
                  &Bicop::get_parameters,
                  &Bicop::set_parameters,
                  "The copula parameter(s).")
    .def_property("var_types",
                  &Bicop::get_var_types,
                  &Bicop::set_var_types,
                  "The type of the two variables.")
    .def_property_readonly("family", &Bicop::get_family, "The copula family.")
    .def_property_readonly("tau", &Bicop::get_tau, "The Kendall's tau.")
    .def_property_readonly("npars",
                           &Bicop::get_npars,
                           "The number of parameters (for nonparametric "
                           "families, a conceptually similar definition).")
    .def("loglik",
         &Bicop::loglik,
         py::arg("u") = Eigen::Matrix<double, Eigen::Dynamic, 2>(),
         bicop_doc.loglik.doc)
    .def_property_readonly(
      "nobs",
      &Bicop::get_nobs,
      "The number of observations (only for fitted objects).")
    .def("aic",
         &Bicop::aic,
         py::arg("u") = Eigen::Matrix<double, Eigen::Dynamic, 2>(),
         bicop_doc.aic.doc)
    .def("bic",
         &Bicop::bic,
         py::arg("u") = Eigen::Matrix<double, Eigen::Dynamic, 2>(),
         bicop_doc.bic.doc)
    .def("mbic",
         &Bicop::mbic,
         py::arg("u") = Eigen::Matrix<double, Eigen::Dynamic, 2>(),
         py::arg("psi0") = 0.9,
         bicop_doc.mbic.doc)
    .def("__repr__",
         [](const Bicop& cop) { return "<pyvinecopulib.Bicop>\n" + cop.str(); })
    .def("str", &Bicop::str, bicop_doc.str.doc)
    .def("parameters_to_tau",
         &Bicop::parameters_to_tau,
         py::arg("parameters"),
         bicop_doc.parameters_to_tau.doc)
    .def("tau_to_parameters",
         &Bicop::tau_to_parameters,
         py::arg("tau"),
         bicop_doc.tau_to_parameters.doc)
    .def("parameters_lower_bounds",
         &Bicop::get_parameters_lower_bounds,
         bicop_doc.get_parameters_lower_bounds.doc)
    .def("parameters_upper_bounds",
         &Bicop::get_parameters_upper_bounds,
         bicop_doc.get_parameters_upper_bounds.doc)
    .def("pdf", &Bicop::pdf, py::arg("u"), bicop_doc.pdf.doc)
    .def("cdf", &Bicop::cdf, py::arg("u"), bicop_doc.cdf.doc)
    .def("hfunc1", &Bicop::hfunc1, py::arg("u"), bicop_doc.hfunc1.doc)
    .def("hfunc2", &Bicop::hfunc2, py::arg("u"), bicop_doc.hfunc2.doc)
    .def("hinv1", &Bicop::hinv1, py::arg("u"), bicop_doc.hinv1.doc)
    .def("hinv2", &Bicop::hinv2, py::arg("u"), bicop_doc.hinv2.doc)
    .def("simulate",
         &Bicop::simulate,
         py::arg("n"),
         py::arg("qrng") = false,
         py::arg("seeds") = std::vector<int>(),
         bicop_doc.simulate.doc)
    .def("fit",
         &Bicop::fit,
         py::arg("data"),
         py::arg_v("controls", FitControlsBicop(), "FitControlsBicop()"),
         bicop_doc.fit.doc)
    .def("select",
         &Bicop::select,
         py::arg("data"),
         py::arg_v("controls", FitControlsBicop(), "FitControlsBicop()"),
         bicop_doc.select.doc);
}
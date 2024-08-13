#pragma once

#include "docstr.hpp"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vinecopulib.hpp>

namespace py = pybind11;
using namespace vinecopulib;

inline void
init_vinecop_class(py::module_& module)
{

  constexpr auto& doc = pyvinecopulib_doc;
  constexpr auto& vinecop_doc = doc.vinecopulib.Vinecop;

  py::class_<Vinecop>(module, "VinecopCpp", vinecop_doc.doc)
    .def(py::init<const size_t>(), vinecop_doc.ctor.doc_1args_d, py::arg("d"))
    .def(py::init<const RVineStructure&,
                  const std::vector<std::vector<Bicop>>&,
                  const std::vector<std::string>&>(),
         py::arg("structure"),
         py::arg("pair_copulas") = std::vector<size_t>(),
         py::arg("var_types") = std::vector<std::string>(),
         vinecop_doc.ctor.doc_2args_structure_constint)
    .def(py::init<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>&,
                  const std::vector<std::vector<Bicop>>&,
                  const std::vector<std::string>&>(),
         py::arg("matrix"),
         py::arg("pair_copulas") = std::vector<size_t>(),
         py::arg("var_types") = std::vector<std::string>(),
         vinecop_doc.ctor.doc_2args_matrix_constint)
    .def(py::init<const Eigen::MatrixXd&,
                  const RVineStructure&,
                  const std::vector<std::string>&,
                  const FitControlsVinecop&>(),
         py::arg("data"),
         py::arg("structure") = RVineStructure(),
         py::arg("var_types") = std::vector<std::string>(),
         py::arg_v("controls", FitControlsVinecop(), "FitControlsVinecop()"),
         vinecop_doc.ctor.doc_4args_data_structure_var_types_controls)
    .def(py::init<const Eigen::MatrixXd&,
                  const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>&,
                  const std::vector<std::string>&,
                  const FitControlsVinecop&>(),
         py::arg("data"),
         py::arg("matrix") =
           Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>(),
         py::arg("var_types") = std::vector<std::string>(),
         py::arg_v("controls", FitControlsVinecop(), "FitControlsVinecop()"),
         vinecop_doc.ctor.doc_4args_data_matrix_var_types_controls)
    .def(py::init<const std::string, bool>(),
         py::arg("filename"),
         py::arg("check") = true,
         vinecop_doc.ctor.doc_2args_filename_check)
    .def("to_json",
         &Vinecop::to_file,
         py::arg("filename"),
         vinecop_doc.to_file.doc)
    .def_property("var_types",
                  &Vinecop::get_var_types,
                  &Vinecop::set_var_types,
                  "The types of each variables.")
    .def_property_readonly(
      "trunc_lvl", &Vinecop::get_trunc_lvl, "The truncation level.")
    .def_property_readonly("dim", &Vinecop::get_dim, "The dimension.")
    .def("get_pair_copula",
         &Vinecop::get_pair_copula,
         "Gets a pair-copula.",
         py::arg("tree"),
         py::arg("edge"))
    .def("get_family",
         &Vinecop::get_family,
         "Gets the family of a pair-copula.",
         py::arg("tree"),
         py::arg("edge"))
    .def("get_rotation",
         &Vinecop::get_rotation,
         "Gets the rotation of a pair-copula.",
         py::arg("tree"),
         py::arg("edge"))
    .def("get_parameters",
         &Vinecop::get_parameters,
         "Gets the parameters of a pair-copula.",
         py::arg("tree"),
         py::arg("edge"))
    .def("get_tau",
         &Vinecop::get_tau,
         "Gets the kendall's tau of a pair-copula.",
         py::arg("tree"),
         py::arg("edge"))
    .def_property_readonly(
      "pair_copulas", &Vinecop::get_all_pair_copulas, "All pair-copulas.")
    .def_property_readonly(
      "families", &Vinecop::get_all_families, "Families of all pair-copulas.")
    .def_property_readonly("rotations",
                           &Vinecop::get_all_rotations,
                           "The rotations of all pair-copulas.")
    .def_property_readonly("parameters",
                           &Vinecop::get_all_parameters,
                           "The parameters of all pair-copulas.")
    .def_property_readonly(
      "taus", &Vinecop::get_all_taus, "The Kendall's taus of all pair-copulas.")
    .def_property_readonly(
      "order", &Vinecop::get_order, "The R-vine structure's order.")
    .def_property_readonly(
      "matrix", &Vinecop::get_matrix, "The R-vine structure's matrix.")
    .def_property_readonly(
      "structure", &Vinecop::get_rvine_structure, "The R-vine structure.")
    .def_property_readonly(
      "npars", &Vinecop::get_npars, "The total number of parameters.")
    .def_property_readonly(
      "nobs",
      &Vinecop::get_nobs,
      "The number of observations (for fitted objects only).")
    .def_property_readonly("threshold",
                           &Vinecop::get_threshold,
                           "The threshold (for thresholded copulas only).")
    .def("select",
         &Vinecop::select,
         py::arg("data"),
         py::arg_v("controls", FitControlsVinecop(), "FitControlsVinecop()"),
         vinecop_doc.select.doc)
    .def("pdf",
         &Vinecop::pdf,
         py::arg("u"),
         py::arg("num_threads") = 1,
         vinecop_doc.pdf.doc)
    .def("cdf",
         &Vinecop::cdf,
         py::arg("u"),
         py::arg("N") = 10000,
         py::arg("num_threads") = 1,
         py::arg("seeds") = std::vector<int>(),
         vinecop_doc.cdf.doc)
    .def("simulate",
         &Vinecop::simulate,
         py::arg("n"),
         py::arg("qrng") = false,
         py::arg("num_threads") = 1,
         py::arg("seeds") = std::vector<int>(),
         vinecop_doc.simulate.doc)
    .def("rosenblatt",
         &Vinecop::rosenblatt,
         py::arg("u"),
         py::arg("num_threads") = 1,
         vinecop_doc.rosenblatt.doc)
    .def("inverse_rosenblatt",
         &Vinecop::inverse_rosenblatt,
         py::arg("u"),
         py::arg("num_threads") = 1,
         vinecop_doc.inverse_rosenblatt.doc)
    .def("loglik",
         &Vinecop::loglik,
         py::arg("u") = Eigen::MatrixXd(),
         py::arg("num_threads") = 1,
         vinecop_doc.loglik.doc)
    .def("aic",
         &Vinecop::aic,
         py::arg("u") = Eigen::MatrixXd(),
         py::arg("num_threads") = 1,
         vinecop_doc.aic.doc)
    .def("bic",
         &Vinecop::bic,
         py::arg("u") = Eigen::MatrixXd(),
         py::arg("num_threads") = 1,
         vinecop_doc.bic.doc)
    .def("mbicv",
         &Vinecop::mbicv,
         py::arg("u") = Eigen::MatrixXd(),
         py::arg("psi0") = 0.9,
         py::arg("num_threads") = 1,
         vinecop_doc.mbicv.doc)
    .def("__repr__",
         [](const Vinecop& cop) {
           return "<pyvinecopulib.Vinecop>\n" + cop.str();
         })
    .def("str", &Vinecop::str, vinecop_doc.str.doc)
    .def("truncate",
         &Vinecop::truncate,
         py::arg("trunc_lvl"),
         vinecop_doc.truncate.doc);
}
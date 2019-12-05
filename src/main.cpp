#include "docstr.hpp"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vinecopulib.hpp>

namespace py = pybind11;
using namespace vinecopulib;

PYBIND11_MODULE(pyvinecopulib, pv)
{

  constexpr auto& doc = pyvinecopulib_doc;
  constexpr auto& tools_stat_doc = doc.vinecopulib.tools_stats;
  constexpr auto& bicop_doc = doc.vinecopulib.Bicop;
  constexpr auto& bicopfamily_doc = doc.vinecopulib.BicopFamily;
  constexpr auto& fitcontrolsbicop_doc = doc.vinecopulib.FitControlsBicop;
  constexpr auto& rvinestructure_doc = doc.vinecopulib.RVineStructure;
  constexpr auto& dvinestructure_doc = doc.vinecopulib.DVineStructure;
  constexpr auto& cvinestructure_doc = doc.vinecopulib.CVineStructure;
  constexpr auto& vinecop_doc = doc.vinecopulib.Vinecop;
  constexpr auto& fitcontrolsvinecop_doc = doc.vinecopulib.FitControlsVinecop;

  pv.doc() = R"pbdoc(
  The pyvinecopulib package
  -------------------------
  )pbdoc";

  py::enum_<BicopFamily>(pv, "BicopFamily", py::arithmetic(), R"pbdoc(
   A bivariate copula family identifier.

   The following convenient sets of families are also provided:

   - ``all`` contains all the families,
   - ``parametric`` contains the parametric families (all except ``tll``),
   - ``nonparametric`` contains the nonparametric families
     (``indep`` and ``tll``)
   - ``onepar`` contains the parametric families with a single parameter,
     (``gaussian``, ``clayton``, ``gumbel``, ``frank``, and ``joe``),
   - ``twopar`` contains the parametric families with two parameters
     (``student``, ``bb1``, ``bb6``, ``bb7``, and ``bb8``),
   - ``elliptical`` contains the elliptical families,
   - ``archimedean`` contains the archimedean families,
   - ``BB`` contains the BB families,
   - ``itau`` families for which estimation by Kendall's tau inversion is
     available (``indep``, ``gaussian``, ``student``, ``clayton``,
     ``gumbel``, ``frank``, ``joe``),
   - ``lt`` contains the families that are lower-tail dependent,
   - ``ut`` contains the families that are upper-tail dependent.
   )pbdoc")
    .value("indep", BicopFamily::indep, bicopfamily_doc.indep.doc)
    .value("gaussian", BicopFamily::gaussian, bicopfamily_doc.gaussian.doc)
    .value("student", BicopFamily::student, bicopfamily_doc.student.doc)
    .value("clayton", BicopFamily::clayton, bicopfamily_doc.clayton.doc)
    .value("gumbel", BicopFamily::gumbel, bicopfamily_doc.gumbel.doc)
    .value("frank", BicopFamily::frank, bicopfamily_doc.frank.doc)
    .value("joe", BicopFamily::joe, bicopfamily_doc.joe.doc)
    .value("bb1", BicopFamily::bb1, bicopfamily_doc.bb1.doc)
    .value("bb6", BicopFamily::bb6, bicopfamily_doc.bb6.doc)
    .value("bb7", BicopFamily::bb7, bicopfamily_doc.bb7.doc)
    .value("bb8", BicopFamily::bb8, bicopfamily_doc.bb8.doc)
    .value("tll", BicopFamily::tll, bicopfamily_doc.tll.doc)
    .export_values();

  pv.attr("all") = bicop_families::all;
  pv.attr("parametric") = bicop_families::parametric;
  pv.attr("nonparametric") = bicop_families::nonparametric;
  pv.attr("one_par") = bicop_families::one_par;
  pv.attr("two_par") = bicop_families::two_par;
  pv.attr("elliptical") = bicop_families::elliptical;
  pv.attr("archimedean") = bicop_families::archimedean;
  pv.attr("bb") = bicop_families::bb;
  pv.attr("lt") = bicop_families::lt;
  pv.attr("ut") = bicop_families::ut;
  pv.attr("itau") = bicop_families::itau;

  py::class_<FitControlsBicop>(pv, "FitControlsBicop", fitcontrolsbicop_doc.doc)
    .def(py::init<std::vector<BicopFamily>,
                  std::string,
                  std::string,
                  double,
                  std::string,
                  const Eigen::VectorXd&,
                  double,
                  bool,
                  size_t>(),
         fitcontrolsbicop_doc.ctor.doc_9args,
         py::arg("family_set") = bicop_families::all,
         py::arg("parametric_method") = "mle",
         py::arg("nonparametric_method") = "quadratic",
         py::arg("nonparametric_mult") = 1.0,
         py::arg("selection_criterion") = "bic",
         py::arg("weights") = Eigen::VectorXd(),
         py::arg("psi0") = 0.9,
         py::arg("preselect_families") = true,
         py::arg("num_threads") = 1)
    /*      .def(py::init<std::string>(), */
    //      "creates default controls except for the parameteric method.",
    //      py::arg("parametric_method"))
    // .def(py::init<std::string, double>(),
    //      "creates default controls except for the nonparametric method.",
    /* py::arg("nonparametric_method"), py::arg("mult") = 1.0) */
    .def_property("family_set",
                  &FitControlsBicop::get_family_set,
                  &FitControlsBicop::set_family_set,
                  "The family set.")
    .def_property("parametric_method",
                  &FitControlsBicop::get_parametric_method,
                  &FitControlsBicop::set_parametric_method,
                  "The fit method for parametric families.")
    .def_property("nonparametric_method",
                  &FitControlsBicop::get_nonparametric_method,
                  &FitControlsBicop::set_nonparametric_method,
                  "The fit method for nonparametric families.")
    .def_property("nonparametric_mult",
                  &FitControlsBicop::get_nonparametric_mult,
                  &FitControlsBicop::set_nonparametric_method,
                  "The multiplier for the smoothing parameters.")
    .def_property("selection_criterion",
                  &FitControlsBicop::get_selection_criterion,
                  &FitControlsBicop::set_selection_criterion,
                  "The selection criterion.")
    .def_property("weights",
                  &FitControlsBicop::get_weights,
                  &FitControlsBicop::set_weights,
                  "The weights for the observations.")
    .def_property("psi0",
                  &FitControlsBicop::get_psi0,
                  &FitControlsBicop::set_psi0,
                  "The prior probability of non-independence.")
    .def_property("preselect_families",
                  &FitControlsBicop::get_preselect_families,
                  &FitControlsBicop::set_preselect_families,
                  "Whether to exclude families based on symmetry properties "
                  "(see ``FitControlsBicop``)")
    .def_property("num_threads",
                  &FitControlsBicop::get_num_threads,
                  &FitControlsBicop::set_num_threads,
                  "The number of threads.");
  /* .def("__repr__", */
  //      [](const FitControlsBicop& ctrls) {
  //        return "<pyvinecopulib.FitControlsBicop>\n" + ctrls.str();
  //      })
  // .def("str",
  //      &FitControlsBicop::str,
  /* "summarizes the controls into a string (can be used for printing)."); */

  py::class_<Bicop>(pv, "Bicop", bicop_doc.doc)
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
    .def("to_json", &Bicop::to_json, py::arg("filename"), bicop_doc.to_json.doc)
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

  py::class_<RVineStructure>(pv, "RVineStructure", rvinestructure_doc.doc)
    .def(py::init<const size_t&, const size_t&>(),
         py::arg("d") = static_cast<size_t>(1),
         py::arg("trunc_lvl") = std::numeric_limits<size_t>::max(),
         rvinestructure_doc.ctor.doc_2args_d_trunc_lvl)
    .def(py::init<const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>&,
                  bool>(),
         py::arg("mat"),
         py::arg("check") = true,
         rvinestructure_doc.ctor.doc_2args_mat_check)
    .def(py::init<const std::vector<size_t>&, const size_t&, bool>(),
         py::arg("order"),
         py::arg("trunc_lvl") = std::numeric_limits<size_t>::max(),
         py::arg("check") = true,
         rvinestructure_doc.ctor.doc_3args_order_trunc_lvl_check)
    .def(py::init<const std::string, bool>(),
         py::arg("filename"),
         py::arg("check") = true,
         rvinestructure_doc.ctor.doc_2args_filename_check)
    .def("to_json",
         &RVineStructure::to_json,
         py::arg("filename"),
         rvinestructure_doc.to_json.doc)
    .def_property_readonly("dim", &RVineStructure::get_dim, "The dimension.")
    .def_property_readonly(
      "trunc_lvl", &RVineStructure::get_trunc_lvl, "The truncation level.")
    .def_property_readonly("order",
                           (std::vector<size_t>(RVineStructure::*)() const) &
                             RVineStructure::get_order,
                           "The variable order.")
    .def("struct_array",
         &RVineStructure::struct_array,
         py::arg("tree"),
         py::arg("edge"),
         py::arg("natural_order") = false,
         rvinestructure_doc.struct_array.doc)
    .def("truncate",
         &RVineStructure::truncate,
         py::arg("trunc_lvl"),
         rvinestructure_doc.truncate.doc)
    .def_static("simulate",
                &RVineStructure::simulate,
                py::arg("d"),
                py::arg("natural order") = false,
                py::arg("seeds") = std::vector<size_t>(),
                rvinestructure_doc.simulate.doc)
    .def("__repr__",
         [](const RVineStructure& rvs) {
           return "<pyvinecopulib.RVineStructure>\n" + rvs.str();
         })
    .def("str", &RVineStructure::str, rvinestructure_doc.str.doc);

  py::class_<DVineStructure, RVineStructure>(
    pv, "DVineStructure", dvinestructure_doc.doc)
    .def(py::init<const std::vector<size_t>&>(),
         py::arg("order"),
         dvinestructure_doc.ctor.doc_1args)
    .def(py::init<const std::vector<size_t>&, size_t>(),
         py::arg("order"),
         py::arg("trunc_lvl"),
         dvinestructure_doc.ctor.doc_2args)
    .def("__repr__", [](const DVineStructure& rvs) {
      return "<pyvinecopulib.DVineStructure>\n" + rvs.str();
    });

  py::class_<CVineStructure, RVineStructure>(
    pv, "CVineStructure", cvinestructure_doc.doc)
    .def(py::init<const std::vector<size_t>&>(),
         cvinestructure_doc.ctor.doc_1args,
         py::arg("order"))
    .def(py::init<const std::vector<size_t>&, size_t>(),
         py::arg("order"),
         py::arg("trunc_lvl"),
         cvinestructure_doc.ctor.doc_2args)
    .def("__repr__", [](const CVineStructure& rvs) {
      return "<pyvinecopulib.CVineStructure>\n" + rvs.str();
    });

  py::class_<FitControlsVinecop>(
    pv, "FitControlsVinecop", fitcontrolsvinecop_doc.doc)
    .def(py::init<std::vector<BicopFamily>,
                  std::string,
                  std::string,
                  double,
                  size_t,
                  std::string,
                  double,
                  std::string,
                  const Eigen::VectorXd&,
                  double,
                  bool,
                  bool,
                  bool,
                  bool,
                  size_t>(),
         py::arg("family_set") = bicop_families::all,
         py::arg("parametric_method") = "mle",
         py::arg("nonparametric_method") = "quadratic",
         py::arg("nonparametric_mult") = 1.0,
         py::arg("trunc_lvl") = std::numeric_limits<size_t>::max(),
         py::arg("tree_criterion") = "tau",
         py::arg("threshold") = 0.0,
         py::arg("selection_criterion") = "bic",
         py::arg("weights") = Eigen::VectorXd(),
         py::arg("psi0") = 0.9,
         py::arg("preselect_families") = true,
         py::arg("select_trunc_lvl") = false,
         py::arg("select_threshold") = false,
         py::arg("show_trace") = false,
         py::arg("num_threads") = 1,
         fitcontrolsvinecop_doc.ctor.doc_15args)
    .def_property("family_set",
                  &FitControlsVinecop::get_family_set,
                  &FitControlsVinecop::set_family_set,
                  "The family set.")
    .def_property("parametric_method",
                  &FitControlsVinecop::get_parametric_method,
                  &FitControlsVinecop::set_parametric_method,
                  "The fit method for parametric families.")
    .def_property("nonparametric_method",
                  &FitControlsVinecop::get_nonparametric_method,
                  &FitControlsVinecop::set_nonparametric_method,
                  "The fit method for nonparametric families.")
    .def_property("nonparametric_mult",
                  &FitControlsVinecop::get_nonparametric_mult,
                  &FitControlsVinecop::set_nonparametric_method,
                  "The multiplier for the smoothing parameters.")
    .def_property("trunc_lvl",
                  &FitControlsVinecop::get_trunc_lvl,
                  &FitControlsVinecop::set_trunc_lvl,
                  "The truncation level.")
    .def_property("tree_criterion",
                  &FitControlsVinecop::get_tree_criterion,
                  &FitControlsVinecop::set_tree_criterion,
                  "The tree criterion.")
    .def_property("threshold",
                  &FitControlsVinecop::get_threshold,
                  &FitControlsVinecop::set_threshold,
                  "The threshold.")
    .def_property("selection_criterion",
                  &FitControlsVinecop::get_selection_criterion,
                  &FitControlsVinecop::set_selection_criterion,
                  "The selection criterion.")
    .def_property("weights",
                  &FitControlsVinecop::get_weights,
                  &FitControlsVinecop::set_weights,
                  "The weights for the observations.")
    .def_property("psi0",
                  &FitControlsVinecop::get_psi0,
                  &FitControlsVinecop::set_psi0,
                  "The prior probability of non-independence.")
    .def_property(
      "preselect_families",
      &FitControlsVinecop::get_preselect_families,
      &FitControlsVinecop::set_preselect_families,
      "Preselection based on symmetry properties (see ``__init__``).")
    .def_property("select_trunc_lvl",
                  &FitControlsVinecop::get_select_trunc_lvl,
                  &FitControlsVinecop::set_select_trunc_lvl,
                  "Whether to select the truncation level.")
    .def_property("select_threshold",
                  &FitControlsVinecop::get_select_threshold,
                  &FitControlsVinecop::set_select_threshold,
                  "Whether to select the threshold.")
    .def_property("show_trace",
                  &FitControlsVinecop::get_show_trace,
                  &FitControlsVinecop::set_show_trace,
                  "Whether to show the trace.")
    .def_property("num_threads",
                  &FitControlsVinecop::get_num_threads,
                  &FitControlsVinecop::set_num_threads,
                  "The number of threads.");
  /*   .def("__repr__", */
  //    [](const FitControlsVinecop& ctrls) {
  //      return "<pyvinecopulib.FitControlsRinecop>\n" + ctrls.str();
  //    })
  // .def("str",
  //      &FitControlsVinecop::str,
  /* "summarizes the controls into a string (can be used for printing)."); */

  py::class_<Vinecop>(pv, "Vinecop", vinecop_doc.doc)
    .def(py::init<const size_t>(), vinecop_doc.ctor.doc_1args_d, py::arg("d"))
    .def(py::init<const RVineStructure&,
                  const std::vector<std::vector<Bicop>>&,
                  const std::vector<std::string>&>(),
         py::arg("structure"),
         py::arg("pair_copulas") = std::vector<size_t>(),
         py::arg("var_types") = std::vector<std::string>(),
         vinecop_doc.ctor.doc_3args_structure_pair_copulas_var_types)
    .def(py::init<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>&,
                  const std::vector<std::vector<Bicop>>&,
                  const std::vector<std::string>&>(),
         py::arg("matrix"),
         py::arg("pair_copulas") = std::vector<size_t>(),
         py::arg("var_types") = std::vector<std::string>(),
         vinecop_doc.ctor.doc_3args_matrix_pair_copulas_var_types)
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
         &Vinecop::to_json,
         py::arg("filename"),
         vinecop_doc.to_json.doc)
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
         py::arg("qrn") = false,
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

  pv.def("simulate_uniform",
         &tools_stats::simulate_uniform,
         tools_stat_doc.simulate_uniform.doc,
         py::arg("n"),
         py::arg("d"),
         py::arg("qrng") = false,
         py::arg("seeds") = std::vector<int>());

  pv.def("sobol",
         &tools_stats::sobol,
         tools_stat_doc.sobol.doc,
         py::arg("n"),
         py::arg("d"),
         py::arg("seeds") = std::vector<int>());

  pv.def("ghalton",
         &tools_stats::ghalton,
         tools_stat_doc.ghalton.doc,
         py::arg("n"),
         py::arg("d"),
         py::arg("seeds") = std::vector<int>());

  pv.def("to_pseudo_obs",
         &tools_stats::to_pseudo_obs,
         tools_stat_doc.to_pseudo_obs.doc,
         py::arg("x"),
         py::arg("ties_method") = "average");

#ifdef VERSION_INFO
  pv.attr("__version__") = VERSION_INFO;
#else
  pv.attr("__version__") = "dev";
#endif
}

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vinecopulib.hpp>

namespace py = pybind11;
using namespace vinecopulib;

PYBIND11_MODULE(pyvinecopulib, pv)
{

  /*   py::options options; */
  /* options.disable_function_signatures(); */

  pv.doc() = R"pbdoc(
        pyvinecopulib library
        -----------------------

        .. currentmodule:: pyvinecopulib

        .. autosummary::
           :toctree: _generate

           Bicop
           BicopFamily
           FitControlsBicop
           Vinecop
           FitControlsVinecop
    )pbdoc";

  py::module pv_bicop_families =
    pv.def_submodule("bicop_families",
                     "A submodule of 'pyvinecopulib' with convenience "
                     "definitions for bivariate families");

  py::enum_<BicopFamily>(pv, "BicopFamily", py::arithmetic())
    .value("indep", BicopFamily::indep)
    .value("gaussian", BicopFamily::gaussian)
    .value("student", BicopFamily::student)
    .value("clayton", BicopFamily::clayton)
    .value("gumbel", BicopFamily::gumbel)
    .value("frank", BicopFamily::frank)
    .value("joe", BicopFamily::joe)
    .value("bb1", BicopFamily::bb1)
    .value("bb6", BicopFamily::bb6)
    .value("bb7", BicopFamily::bb7)
    .value("bb8", BicopFamily::bb8)
    .value("tll", BicopFamily::tll);

  pv_bicop_families.attr("all") = bicop_families::all;
  pv_bicop_families.attr("parametric") = bicop_families::parametric;
  pv_bicop_families.attr("nonparametric") = bicop_families::nonparametric;
  pv_bicop_families.attr("one_par") = bicop_families::one_par;
  pv_bicop_families.attr("two_par") = bicop_families::two_par;
  pv_bicop_families.attr("elliptical") = bicop_families::elliptical;
  pv_bicop_families.attr("archimedean") = bicop_families::archimedean;
  pv_bicop_families.attr("bb") = bicop_families::bb;
  pv_bicop_families.attr("rotationless") = bicop_families::rotationless;
  pv_bicop_families.attr("lt") = bicop_families::lt;
  pv_bicop_families.attr("ut") = bicop_families::ut;
  pv_bicop_families.attr("itau") = bicop_families::itau;
  pv_bicop_families.attr("flip_by_rotation") = bicop_families::flip_by_rotation;

  py::class_<FitControlsBicop>(pv, "FitControlsBicop")
    .def(py::init<std::vector<BicopFamily>,
                  std::string,
                  std::string,
                  double,
                  std::string,
                  const Eigen::VectorXd&,
                  double,
                  bool,
                  size_t>(),
         "creates the controls for fitting bivariate copula models.",
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
                  &FitControlsBicop::set_family_set)
    .def_property("parametric_method",
                  &FitControlsBicop::get_parametric_method,
                  &FitControlsBicop::set_parametric_method)
    .def_property("nonparametric_method",
                  &FitControlsBicop::get_nonparametric_method,
                  &FitControlsBicop::set_nonparametric_method)
    .def_property("nonparametric_mult",
                  &FitControlsBicop::get_nonparametric_mult,
                  &FitControlsBicop::set_nonparametric_method)
    .def_property("selection_criterion",
                  &FitControlsBicop::get_selection_criterion,
                  &FitControlsBicop::set_selection_criterion)
    .def_property(
      "weights", &FitControlsBicop::get_weights, &FitControlsBicop::set_weights)
    .def_property(
      "psi0", &FitControlsBicop::get_psi0, &FitControlsBicop::set_psi0)
    .def_property("preselect_families",
                  &FitControlsBicop::get_preselect_families,
                  &FitControlsBicop::set_preselect_families)
    .def_property("num_threads",
                  &FitControlsBicop::get_num_threads,
                  &FitControlsBicop::set_num_threads);

  py::class_<Bicop>(pv, "Bicop")
    .def(py::init<const BicopFamily, const int, const Eigen::MatrixXd&>(),
         "creates a specific bivariate copula model.",
         py::arg("family") = BicopFamily::indep,
         py::arg("rotation") = 0,
         py::arg("parameters") = Eigen::MatrixXd())
    .def(py::init<const Eigen::Matrix<double, Eigen::Dynamic, 2>&,
                  const FitControlsBicop&>(),
         "create a copula model from the data, equivalent to cop = Bicop(); "
         "cop.select(data, controls).",
         py::arg("data"),
         py::arg("controls") = FitControlsBicop())
    .def(py::init<const std::string>(),
         "creates from a JSON file.",
         py::arg("filename"))
    .def("to_json",
         &Bicop::to_json,
         "writes the copula object into a JSON file.",
         py::arg("filename"))
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
    .def_property_readonly(
      "npars",
      &Bicop::get_npars,
      "The number of parameters. For nonparametric families, there is a "
      "conceptually similar definition in the sense that it can be used in "
      "the calculation of fit statistics.")
    .def("loglik",
         &Bicop::loglik,
         "computes the log-likelihood (for fitted objects, passing an "
         "empty 'u' returns the fitted criterion).",
         py::arg("u") = Eigen::Matrix<double, Eigen::Dynamic, 2>())
    .def_property_readonly(
      "nobs",
      &Bicop::get_nobs,
      "The number of observations (for fitted objects only).")
    .def("aic",
         &Bicop::aic,
         "computes the Akaike Information Criterion (for fitted objects, "
         "passing an empty 'u' returns the fitted criterion).",
         py::arg("u") = Eigen::Matrix<double, Eigen::Dynamic, 2>())
    .def("bic",
         &Bicop::bic,
         "computes the Bayesian Information Criterion (for fitted objects, "
         "passing an empty 'u' returns the fitted criterion).",
         py::arg("u") = Eigen::Matrix<double, Eigen::Dynamic, 2>())
    .def("mbic",
         &Bicop::mbic,
         "computes the Modified Bayesian Information Criterion (for "
         "fitted objects, passing an "
         "empty 'u' returns the fitted criterion).",
         py::arg("u") = Eigen::Matrix<double, Eigen::Dynamic, 2>(),
         py::arg("psi0") = 0.9)
    .def("__repr__",
         [](const Bicop& cop) { return "<pyvinecopulib.Bicop>\n" + cop.str(); })
    .def("str",
         &Bicop::str,
         "summarizes the model into a string (can be used for printing).")
    .def("parameters_to_tau",
         &Bicop::parameters_to_tau,
         "returns the Kendall's tau corresponding to the parameters passed "
         "as arguments.",
         py::arg("parameters"))
    .def("tau_to_parameters",
         &Bicop::tau_to_parameters,
         "returns the parameters corresponding to the Kendall's tau passed "
         "as arguments.",
         py::arg("tau"))
    .def("flip",
         &Bicop::flip,
         "adjust's the copula model to a change in the variable order.")
    .def("parameters_lower_bounds",
         &Bicop::get_parameters_lower_bounds,
         "returns the lower bounds for the copula's parameters.")
    .def("parameters_upper_bounds",
         &Bicop::get_parameters_upper_bounds,
         "returns the upper bounds for the copula's parameters.")
    .def("pdf", &Bicop::pdf, "evaluates the copula density.", py::arg("u"))
    .def("cdf", &Bicop::cdf, "evaluates the copula distribution.", py::arg("u"))
    .def("hfunc1",
         &Bicop::hfunc1,
         "evaluates the first h-function, that is the partial derivative of "
         "the copula distribution w.r.t. the first argument.",
         py::arg("u"))
    .def("hfunc2",
         &Bicop::hfunc2,
         "evaluates the second h-function, that is the partial derivative of "
         "the copula distribution w.r.t. the second argument.",
         py::arg("u"))
    .def("hinv1",
         &Bicop::hinv1,
         "evaluates the inverse of the first h-function (hfunc1) w.r.t. the "
         "second argument.",
         py::arg("u"))
    .def("hinv2",
         &Bicop::hinv2,
         "evaluates the inverse of the second h-function (hfunc2) w.r.t. the "
         "first argument.",
         py::arg("u"))
    .def("simulate",
         &Bicop::simulate,
         "simulates from the bivariate model.",
         py::arg("n"),
         py::arg("qrng") = false,
         py::arg("seeds") = std::vector<int>())
    .def("fit",
         &Bicop::fit,
         "fits a bivariate copula (with fixed family) to data.",
         py::arg("data"),
         py::arg("controls") = FitControlsBicop())
    .def("select",
         &Bicop::select,
         "selects the best fitting model.",
         py::arg("data"),
         py::arg("controls") = FitControlsBicop());

  py::class_<RVineStructure>(pv, "RVineStructure")
    .def(py::init<const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>&,
                  bool>(),
         "creates an array representation of regular vine structures.",
         py::arg("mat"),
         py::arg("check") = true)
    .def(py::init<const std::vector<size_t>&, bool>(),
         "creates a D-vine with given ordering of variables.",
         py::arg("order"),
         py::arg("check") = true)
    .def(py::init<const std::string, bool>(),
         "creates a structure from a JSON file.",
         py::arg("filename"),
         py::arg("check") = true)
    .def("to_json",
         &RVineStructure::to_json,
         "writes the structure in a JSON file.",
         py::arg("filename"))
    .def_property_readonly("dim", &RVineStructure::get_dim, "The dimension.")
    .def_property_readonly(
      "trunc_lvl", &RVineStructure::get_trunc_lvl, "The truncation level.")
    .def_property_readonly("order",
                           (std::vector<size_t>(RVineStructure::*)() const) &
                             RVineStructure::get_order,
                           "The variable order.")
    .def("struct_array",
         &RVineStructure::struct_array,
         "accesses elements of the structure array.",
         py::arg("tree"),
         py::arg("edge"),
         py::arg("natural_order") = false)
    .def("truncate",
         &RVineStructure::truncate,
         "truncates the R-vine structure.",
         py::arg("trunc_lvl"))
    .def_static("simulate",
                &RVineStructure::simulate,
                "simulates a random R-vine array.",
                py::arg("d"),
                py::arg("natural order") = false,
                py::arg("seeds") = std::vector<size_t>())
    .def("__repr__",
         [](const RVineStructure& rvs) {
           return "<pyvinecopulib.RVineStructure>\n" + rvs.str();
         })
    .def("str",
         &RVineStructure::str,
         "summarizes the model into a string (can be used for printing).");

  py::class_<DVineStructure>(pv, "DVineStructure")
    .def(py::init<const std::vector<size_t>&>(),
         "creates a D-vine with given ordering of variables.",
         py::arg("order"))
    .def(
      py::init<const std::vector<size_t>&, size_t>(),
      "creates a D-vine with given ordering of variables and truncation level.",
      py::arg("order"),
      py::arg("trunc_lvl"));

  py::class_<CVineStructure>(pv, "CVineStructure")
    .def(py::init<const std::vector<size_t>&>(),
         "creates a C-vine with given ordering of variables.",
         py::arg("order"))
    .def(
      py::init<const std::vector<size_t>&, size_t>(),
      "creates a C-vine with given ordering of variables and truncation level.",
      py::arg("order"),
      py::arg("trunc_lvl"));

  py::class_<FitControlsVinecop>(pv, "FitControlsVinecop")
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
         "creates the controls for fitting bivariate copula models.",
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
         py::arg("num_threads") = 1)
    .def_property("family_set",
                  &FitControlsVinecop::get_family_set,
                  &FitControlsVinecop::set_family_set)
    .def_property("parametric_method",
                  &FitControlsVinecop::get_parametric_method,
                  &FitControlsVinecop::set_parametric_method)
    .def_property("nonparametric_method",
                  &FitControlsVinecop::get_nonparametric_method,
                  &FitControlsVinecop::set_nonparametric_method)
    .def_property("nonparametric_mult",
                  &FitControlsVinecop::get_nonparametric_mult,
                  &FitControlsVinecop::set_nonparametric_method)
    .def_property("trunc_lvl",
                  &FitControlsVinecop::get_trunc_lvl,
                  &FitControlsVinecop::set_trunc_lvl)
    .def_property("tree_criterion",
                  &FitControlsVinecop::get_tree_criterion,
                  &FitControlsVinecop::set_tree_criterion)
    .def_property("threshold",
                  &FitControlsVinecop::get_threshold,
                  &FitControlsVinecop::set_threshold)
    .def_property("selection_criterion",
                  &FitControlsVinecop::get_selection_criterion,
                  &FitControlsVinecop::set_selection_criterion)
    .def_property("weights",
                  &FitControlsVinecop::get_weights,
                  &FitControlsVinecop::set_weights)
    .def_property(
      "psi0", &FitControlsVinecop::get_psi0, &FitControlsVinecop::set_psi0)
    .def_property("preselect_families",
                  &FitControlsVinecop::get_preselect_families,
                  &FitControlsVinecop::set_preselect_families)
    .def_property("select_trunc_lvl",
                  &FitControlsVinecop::get_select_trunc_lvl,
                  &FitControlsVinecop::set_select_trunc_lvl)
    .def_property("select_threshold",
                  &FitControlsVinecop::get_select_threshold,
                  &FitControlsVinecop::set_select_threshold)
    .def_property("show_trace",
                  &FitControlsVinecop::get_show_trace,
                  &FitControlsVinecop::set_show_trace)
    .def_property("num_threads",
                  &FitControlsVinecop::get_num_threads,
                  &FitControlsVinecop::set_num_threads);

  py::class_<Vinecop>(pv, "Vinecop")
    .def(py::init<const size_t>(),
         "creates a D-vine with independence copulas.",
         py::arg("d"))
    .def(py::init<const RVineStructure&>(),
         "creates a vine copula with structure specified by an RVineStructure "
         "object; all pair-copulas are set to independence.",
         py::arg("structure"))
    .def(py::init<const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>&,
                  const bool>(),
         "creates a vine copula with structure specified by an R-vine matrix; "
         "all pair-copulas are set to independence.",
         py::arg("matrix"),
         py::arg("check") = true)
    .def(
      py::init<const std::vector<std::vector<Bicop>>&, const RVineStructure&>(),
      "creates an arbitrary vine copula model.",
      py::arg("pair_copulas"),
      py::arg("structure"))
    .def(py::init<const std::vector<std::vector<Bicop>>&,
                  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>&,
                  const bool>(),
         "creates an arbitrary vine copula model.",
         py::arg("pair_copulas"),
         py::arg("matrix"),
         py::arg("check") = true)
    .def(py::init<const std::string, bool>(),
         "creates a vine copula from a JSON file.",
         py::arg("filename"),
         py::arg("check") = true)
    .def(py::init<const Eigen::MatrixXd&, const FitControlsVinecop&>(),
         "constructs a vine copula model from data by creating a model and "
         "calling select() to select the structure, families and parameters.",
         py::arg("data"),
         py::arg("controls") = FitControlsVinecop())
    .def(py::init<const Eigen::MatrixXd&,
                  const RVineStructure&,
                  const FitControlsVinecop&>(),
         "constructs a vine copula model from data by creating a model and "
         "calling select() to select the families and parameters.",
         py::arg("data"),
         py::arg("structure"),
         py::arg("controls") = FitControlsVinecop())
    .def(py::init<const Eigen::MatrixXd&,
                  const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>&,
                  const FitControlsVinecop&,
                  const bool>(),
         "constructs a vine copula model from data by creating a model and "
         "calling select() to select the families and parameters.",
         py::arg("data"),
         py::arg("matrix"),
         py::arg("controls") = FitControlsVinecop(),
         py::arg("check") = true)
    .def("to_json",
         &Vinecop::to_json,
         "write a vine copula to a JSON file.",
         py::arg("filename"))
    .def_property("var_types",
                  &Vinecop::get_var_types,
                  &Vinecop::set_var_types,
                  "The types of each variables.")
    .def_property_readonly(
      "trunc_lvl", &Vinecop::get_trunc_lvl, "The truncation level.")
    .def_property_readonly("dim", &Vinecop::get_dim, "The dimension.")
    /* .def("get_pair_copula", */
    //      &Vinecop::get_pair_copula,
    //      "extracts a pair-copula.",
    //      py::arg("tree"),
    //      py::arg("edge"))
    // .def("get_family",
    //      &Vinecop::get_family,
    //      "extracts the family of a pair-copula.",
    //      py::arg("tree"),
    //      py::arg("edge"))
    // .def("get_rotation",
    //      &Vinecop::get_rotation,
    //      "extracts the rotation of a pair-copula.",
    //      py::arg("tree"),
    //      py::arg("edge"))
    // .def("get_parameters",
    //      &Vinecop::get_parameters,
    //      "extracts the parameters of a pair-copula.",
    //      py::arg("tree"),
    //      py::arg("edge"))
    // .def("get_tau",
    //      &Vinecop::get_tau,
    //      "extracts the Kendall's tau of a pair-copula.",
    //      py::arg("tree"),
    /* py::arg("edge")) */
    .def_property_readonly("get_all_pair_copulas",
                           &Vinecop::get_all_pair_copulas,
                           "extracts all pair-copulas.")
    .def_property_readonly("get_all_families",
                           &Vinecop::get_all_families,
                           "extracts the families of all pair-copulas.")
    .def_property_readonly("get_all_rotations",
                           &Vinecop::get_all_rotations,
                           "extracts the rotations of all pair-copulas.")
    .def_property_readonly("get_all_parameters",
                           &Vinecop::get_all_parameters,
                           "extracts the parameters of all pair-copulas.")
    .def_property_readonly("get_all_taus",
                           &Vinecop::get_all_taus,
                           "extracts the Kendall's taus of all pair-copulas.")
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
         "automatically fits and selects a vine copula model.",
         py::arg("data"),
         py::arg("controls") = FitControlsVinecop())
    .def("pdf",
         &Vinecop::pdf,
         "returns the probability density function.",
         py::arg("u"),
         py::arg("num_threads") = 1)
    .def("cdf",
         &Vinecop::cdf,
         "returns the cumulative distribution.",
         py::arg("u"),
         py::arg("N") = 10000,
         py::arg("num_threads") = 1,
         py::arg("seeds") = std::vector<int>())
    .def("simulate",
         &Vinecop::simulate,
         "sample (quasi-)random numbers from the model.",
         py::arg("n"),
         py::arg("qrn") = false,
         py::arg("num_threads") = 1,
         py::arg("seeds") = std::vector<int>())
    .def("rosenblatt",
         &Vinecop::rosenblatt,
         "computes the Rosenblatt transform.",
         py::arg("u"),
         py::arg("num_threads") = 1)
    .def("inverse_rosenblatt",
         &Vinecop::inverse_rosenblatt,
         "computes the inverse Rosenblatt transform.",
         py::arg("u"),
         py::arg("num_threads") = 1)
    .def("loglik",
         &Vinecop::loglik,
         "computes the log-likelihood (for fitted objects, passing an "
         "empty 'u' returns the fitted criterion).",
         py::arg("u") = Eigen::MatrixXd(),
         py::arg("num_threads") = 1)
    .def("aic",
         &Vinecop::aic,
         "computes the Akaike Information Criterion (for fitted objects, "
         "passing an empty 'u' returns the fitted criterion).",
         py::arg("u") = Eigen::MatrixXd(),
         py::arg("num_threads") = 1)
    .def("bic",
         &Vinecop::bic,
         "computes the Bayesian Information Criterion (for fitted objects, "
         "passing an empty 'u' returns the fitted criterion).",
         py::arg("u") = Eigen::MatrixXd(),
         py::arg("num_threads") = 1)
    .def("mbicv",
         &Vinecop::mbicv,
         "computes the modified Bayesiand Information Criterion for Vines (for "
         "fitted objects, passing an empty 'u' returns the fitted criterion).",
         py::arg("u") = Eigen::MatrixXd(),
         py::arg("psi0") = 0.9,
         py::arg("num_threads") = 1)
    .def("__repr__",
         [](const Vinecop& cop) {
           return "<pyvinecopulib.Vinecop>\n" + cop.str();
         })
    .def("str",
         &Vinecop::str,
         "summarizes the model into a string (can be used for printing).")
    .def("truncate",
         &Vinecop::truncate,
         "truncates the vine copula model.",
         py::arg("trunc_lvl"));

  pv.def("simulate_uniform",
         &tools_stats::simulate_uniform,
         "simulate uniform random numbers.",
         py::arg("n"),
         py::arg("d"),
         py::arg("qrng") = false,
         py::arg("seeds") = std::vector<int>());

  pv.def(
    "to_pseudo_obs",
    &tools_stats::to_pseudo_obs,
    "applies the empirical probability integral transform to a data matrix.",
    py::arg("x"),
    py::arg("ties_method") = "average");

#ifdef VERSION_INFO
  pv.attr("__version__") = VERSION_INFO;
#else
  pv.attr("__version__") = "dev";
#endif
}

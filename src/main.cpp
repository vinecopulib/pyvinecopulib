#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vinecopulib.hpp>

namespace py = pybind11;
using namespace vinecopulib;

PYBIND11_MODULE(pyvinecopulib, pv) {
  pv.doc() = R"pbdoc(
        Pyvinecopulib library
        -----------------------

        .. currentmodule:: pyvinecopulib

        .. autosummary::
           :toctree: _generate

           simulate_uniform
    )pbdoc";

  py::module pv_bicop_families = pv.def_submodule(
      "bicop_families", "A submodule of 'pyvinecopulib' with convenience "
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
      .def(py::init<std::vector<BicopFamily>, std::string, std::string, double,
                    std::string, const Eigen::VectorXd &, double, bool,
                    size_t>(),
           "creates the controls for fitting bivariate copula models.",
           py::arg("family_set") = bicop_families::all,
           py::arg("parametric_method") = "mle",
           py::arg("nonparametric_method") = "quadratic",
           py::arg("nonparametric_mult") = 1.0,
           py::arg("selection_criterion") = "bic",
           py::arg("weights") = Eigen::VectorXd(), py::arg("psi0") = 0.9,
           py::arg("preselect_families") = true, py::arg("num_threads") = 1)
      /*      .def(py::init<std::string>(), */
      //      "creates default controls except for the parameteric method.",
      //      py::arg("parametric_method"))
      // .def(py::init<std::string, double>(),
      //      "creates default controls except for the nonparametric method.",
      /* py::arg("nonparametric_method"), py::arg("mult") = 1.0) */
      .def_property("family_set", &FitControlsBicop::get_family_set,
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
      .def_property("weights", &FitControlsBicop::get_weights,
                    &FitControlsBicop::set_weights)
      .def_property("psi0", &FitControlsBicop::get_psi0,
                    &FitControlsBicop::set_psi0)
      .def_property("preselect_families",
                    &FitControlsBicop::get_preselect_families,
                    &FitControlsBicop::set_preselect_families)
      .def_property("num_threads", &FitControlsBicop::get_num_threads,
                    &FitControlsBicop::set_num_threads);

  py::class_<Bicop>(pv, "Bicop")
      .def(py::init<const BicopFamily, const int, const Eigen::MatrixXd &>(),
           "creates a specific bivariate copula model.",
           py::arg("family") = BicopFamily::indep, py::arg("rotation") = 0,
           py::arg("parameters") = Eigen::MatrixXd())
      .def(py::init<const Eigen::Matrix<double, Eigen::Dynamic, 2> &,
                    const FitControlsBicop &>(),
           "create a copula model from the data, equivalent to cop = Bicop(); "
           "cop.select(data, controls).",
           py::arg("data"), py::arg("controls") = FitControlsBicop())
      .def(py::init<const std::string>(), "creates from a JSON file.",
           py::arg("filename"))
      .def("to_json", &Bicop::to_json,
           "writes the copula object into a JSON file.", py::arg("filename"))
      .def_property("rotation", &Bicop::get_rotation, &Bicop::set_rotation,
                    "The copula rotation.")
      .def_property("parameters", &Bicop::get_parameters,
                    &Bicop::set_parameters, "The copula parameter(s).")
      .def_property_readonly("family", &Bicop::get_family, "The copula family.")
      .def_property_readonly("tau", &Bicop::get_tau, "The Kendall's tau.")
      .def_property_readonly(
          "npars", &Bicop::get_npars,
          "The number of parameters. For nonparametric families, there is a "
          "conceptually similar definition in the sense that it can be used in "
          "the calculation of fit statistics.")
      .def("loglik", &Bicop::loglik,
           "computes the log-likelihood (for fitted objects, passing an "
           "empty 'u' returns the fitted criterion).",
           py::arg("u") = Eigen::Matrix<double, Eigen::Dynamic, 2>())
      .def("nobs", &Bicop::get_nobs,
           "returns the number of observations (for fitted objects only).")
      .def("aic", &Bicop::aic,
           "computes the Akaike Information Criterion (for fitted objects, "
           "passing an empty 'u' returns the fitted criterion).",
           py::arg("u") = Eigen::Matrix<double, Eigen::Dynamic, 2>())
      .def("bic", &Bicop::bic,
           "computes the Bayesian Information Criterion (for fitted objects, "
           "passing an empty 'u' returns the fitted criterion).",
           py::arg("u") = Eigen::Matrix<double, Eigen::Dynamic, 2>())
      .def("mbic", &Bicop::mbic,
           "computes the Modified Bayesian Information Criterion (for "
           "fitted objects, passing an "
           "empty 'u' returns the fitted criterion).",
           py::arg("u") = Eigen::Matrix<double, Eigen::Dynamic, 2>(),
           py::arg("psi0") = 0.9)
      .def("__repr__",
           [](const Bicop &cop) {
             return "<Bicop, family = " + cop.str() + ">";
           })
      .def("str", &Bicop::str,
           "summarizes the model into a string (can be used for printing).")
      .def("parameters_to_tau", &Bicop::parameters_to_tau,
           "returns the Kendall's tau corresponding to the parameters passed "
           "as arguments.",
           py::arg("parameters"))
      .def("tau_to_parameters", &Bicop::tau_to_parameters,
           "returns the parameters corresponding to the Kendall's tau passed "
           "as arguments.",
           py::arg("tau"))
      .def("flip", &Bicop::flip,
           "adjust's the copula model to a change in the variable order.")
      .def("parameters_lower_bounds", &Bicop::get_parameters_lower_bounds,
           "returns the lower bounds for the copula's parameters.")
      .def("parameters_upper_bounds", &Bicop::get_parameters_upper_bounds,
           "returns the upper bounds for the copula's parameters.")
      .def("pdf", &Bicop::pdf, "evaluates the copula density.", py::arg("u"))
      .def("cdf", &Bicop::cdf, "evaluates the copula distribution.",
           py::arg("u"))
      .def("hfunc1", &Bicop::hfunc1,
           "evaluates the first h-function, that is the partial derivative of "
           "the copula distribution w.r.t. the first argument.",
           py::arg("u"))
      .def("hfunc2", &Bicop::hfunc2,
           "evaluates the second h-function, that is the partial derivative of "
           "the copula distribution w.r.t. the second argument.",
           py::arg("u"))
      .def("hinv1", &Bicop::hinv1,
           "evaluates the inverse of the first h-function (hfunc1) w.r.t. the "
           "second argument.",
           py::arg("u"))
      .def("hinv2", &Bicop::hinv2,
           "evaluates the inverse of the second h-function (hfunc2) w.r.t. the "
           "first argument.",
           py::arg("u"))
      .def("simulate", &Bicop::simulate, "simulates from the bivariate model.",
           py::arg("n"), py::arg("qrng") = false,
           py::arg("seeds") = std::vector<int>())
      .def("fit", &Bicop::fit,
           "fits a bivariate copula (with fixed family) to data.",
           py::arg("data"), py::arg("controls") = FitControlsBicop())
      .def("select", &Bicop::select, "selects the best fitting model.",
           py::arg("data"), py::arg("controls") = FitControlsBicop());

  pv.def("simulate_uniform", &tools_stats::simulate_uniform,
         R"pbdoc(
        Simulate uniform random numbers.

        Simulate a matrix of random numbers, with an option to get quasi random
        numbers or a seed.
    )pbdoc");

  py::class_<RVineStructure>(pv, "RVineStructure")
      .def(py::init<
               const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> &,
               bool>(),
           "creates an array representation of regular vine structures.",
           py::arg("mat"), py::arg("check") = true)
      .def(py::init<const std::vector<size_t> &, bool>(),
           "creates a D-vine with given ordering of variables.",
           py::arg("order"), py::arg("check") = true)
      .def(py::init<const std::string, bool>(),
           "creates a structure from a JSON file.", py::arg("filename"),
           py::arg("check") = true)
      .def("to_json", &RVineStructure::to_json,
           "writes the file in a JSON file.", py::arg("filename"))
      .def_property_readonly("dim", &RVineStructure::get_dim, "The dimension.")
      .def_property_readonly("trunc_lvl", &RVineStructure::get_trunc_lvl,
                             "The truncation level.")
      .def_property_readonly("order",
                             (std::vector<size_t>(RVineStructure::*)() const) &
                                 RVineStructure::get_order,
                             "The variable order.")
      .def("struct_array", &RVineStructure::struct_array, py::arg("tree"),
           py::arg("edge"), py::arg("natural_order") = false)
      .def("truncate", &RVineStructure::truncate, "The truncation level.",
           py::arg("trunc_lvl"))
      .def_static("simulate", &RVineStructure::simulate,
                  "simulates a random R-vine array.", py::arg("d"),
                  py::arg("natural order") = false,
                  py::arg("seeds") = std::vector<size_t>())
      .def("__repr__", [](const RVineStructure &rvs) {
        return "<pyvinecopulib.RVineStructure>\n" + rvs.str();
      });

#ifdef VERSION_INFO
  pv.attr("__version__") = VERSION_INFO;
#else
  pv.attr("__version__") = "dev";
#endif
}

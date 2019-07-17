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

  py::class_<Bicop> bicop(pv, "Bicop");
  bicop.def(py::init<const BicopFamily, const int, const Eigen::MatrixXd &>(),
            "creates a specific bivariate copula model.",
            py::arg("family") = BicopFamily::indep, py::arg("rotation") = 0,
            py::arg("parameters") = Eigen::MatrixXd());
  bicop.def(py::init<const Eigen::Matrix<double, Eigen::Dynamic, 2> &,
                     const FitControlsBicop &>(),
            "create a copula model from the data, equivalent to cop = Bicop(); "
            "cop.select(data, controls).",
            py::arg("data"), py::arg("controls") = FitControlsBicop());
  bicop.def_property_readonly("rotation", &Bicop::get_rotation);
  bicop.def_property_readonly("parameters", &Bicop::get_parameters);
  bicop.def_property_readonly("family", &Bicop::get_family);
  bicop.def("__repr__", [](const Bicop &cop) {
    return "<Bicop, family = " + cop.str() + ">";
  });
  pv.def("simulate_uniform", &tools_stats::simulate_uniform,
         R"pbdoc(
        Simulate uniform random numbers.

        Simulate a matrix of random numbers, with an option to get quasi random
        numbers or a seed.
    )pbdoc");

#ifdef VERSION_INFO
  pv.attr("__version__") = VERSION_INFO;
#else
  pv.attr("__version__") = "dev";
#endif
}

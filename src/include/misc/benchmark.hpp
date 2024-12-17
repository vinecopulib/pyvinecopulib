#include <chrono>
#include <iostream>
#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <string>
#include <vinecopulib.hpp>

using namespace vinecopulib;

inline std::vector<double>
benchmark(const Eigen::MatrixXd& data)
{

  // Define different FitControls configurations
  std::vector<FitControlsVinecop> controls = {
    FitControlsVinecop(bicop_families::itau),
    FitControlsVinecop(bicop_families::itau, "itau"),
    FitControlsVinecop({ BicopFamily::tll }),
  };

  std::vector<double> times; // Store times
  for (auto& control : controls) {
    auto func = [&control](const Eigen::MatrixXd& u) {
      Vinecop vc(u, RVineStructure(), {}, control);
    };
    auto data_copy = data.block(0, 0, data.rows(), data.cols());
    auto start = std::chrono::high_resolution_clock::now();
    func(data_copy);
    auto end = std::chrono::high_resolution_clock::now();
    times.push_back(std::chrono::duration<double>(end - start).count());
  }

  return times;
}

inline void
init_benchmark(nb::module_& m)
{
  m.def("benchmark", &benchmark, "data"_a);
}
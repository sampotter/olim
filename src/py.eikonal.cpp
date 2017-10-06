#include <pybind11/pybind11.h>
#include <pybind11/functional.h>

#include "fmm.h"

#include "basic_marcher.hpp"

namespace py = pybind11;
using namespace py::literals;

PYBIND11_MODULE(eikonal, m) {
  m.doc() = "Testing testing";

  m.def("fmm", &fmm,
        "testing",
        "out"_a.none(false),
        "in"_a.none(false),
        "M"_a,
        "N"_a,
        "h"_a,
        "S"_a,
        "type"_a);

  m.def("fmm3d", &fmm, "testing");

  py::class_<basic_marcher>(m, "BasicMarcher")
    .def(
      py::init<int, int, double, std::function<double(double, double)>, double, double>(),
      "height"_a,
      "width"_a,
      "h"_a,
      "s"_a = py::cpp_function(default_speed_func),
      "x0"_a = 0.0,
      "y0"_a = 0.0)
    .def("run", &basic_marcher::run)
    .def(
      "addBoundaryNode", &basic_marcher::add_boundary_node,
      "i"_a,
      "j"_a,
      "value"_a = 0.0)
    .def(
      "getValue", &basic_marcher::get_value,
      "i"_a,
      "j"_a);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:

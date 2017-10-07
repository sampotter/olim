#include <pybind11/pybind11.h>
#include <pybind11/functional.h>

#include "fmm.h"

#include "basic_marcher.hpp"
#include "basic_marcher_3d.hpp"

namespace py = pybind11;
using namespace py::literals;

using speed_function = std::function<double(double, double)>;
using speed_function_3d = std::function<double(double, double, double)>;

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

  py::class_<basic_marcher>(m, "BasicMarcher", py::buffer_protocol())
    .def_buffer([] (basic_marcher & m_) -> py::buffer_info {
        auto const format =
          py::format_descriptor<basic_marcher::float_type>::format();
        return {
          m_.get_node_pointer(),                       // pointer to buffer
          sizeof(basic_marcher::float_type),           // size of one scalar
          format,                                      // format
          basic_marcher::ndims,                        // number of dimensions
          {m_.get_height(), m_.get_width()},           // buffer dimensions
          {sizeof(node)*m_.get_height(), sizeof(node)} // stride (in bytes)
        };
      })
    .def(
      py::init<int, int, double, speed_function, double, double>(),
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
    .def("getValue", &basic_marcher::get_value, "i"_a, "j"_a);

  py::class_<basic_marcher_3d>(m, "BasicMarcher3D")
    .def(
      py::init<int, int, int, double, speed_function_3d, double, double, double>(),
      "height"_a,
      "width"_a,
      "depth"_a,
      "h"_a,
      "s"_a = py::cpp_function(default_speed_func_3d),
      "x0"_a = 0.0,
      "y0"_a = 0.0,
      "z0"_a = 0.0)
    .def("run", &basic_marcher_3d::run)
    .def(
      "addBoundaryNode", &basic_marcher_3d::add_boundary_node,
      "i"_a,
      "j"_a,
      "k"_a,
      "value"_a = 0.0)
    .def("getValue", &basic_marcher_3d::get_value, "i"_a, "j"_a, "k"_a);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:

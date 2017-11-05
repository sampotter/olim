#!/usr/bin/env python

from string import Template

marchers = {
    'basic_marcher': 'BasicMarcher',
    'olim4_mp0': 'Olim4Mid0',
    'olim4_rhr': 'Olim4Rect',
    'olim8_mp0': 'Olim8Mid0',
    'olim8_mp1': 'Olim8Mid1',
    'olim8_rhr': 'Olim8Rect',
    'olim8hu_mp0': 'Olim8HierMid0',
    'olim8hu_mp1': 'Olim8HierMid1',
    'olim8hu_rhr': 'Olim8HierRect',
    'olim8lut_mp0': 'Olim8LutMid0',
    'olim8lut_mp1': 'Olim8LutMid1',
    'olim8lut_rhr': 'Olim8LutRect',
}

marcher_template = Template('''
  py::class_<${cpp_class_name}>(m, "${py_class_name}", py::buffer_protocol())
    .def_buffer([] (${cpp_class_name} & m_) -> py::buffer_info {
        auto const format =
          py::format_descriptor<${cpp_class_name}::float_type>::format();
        return {
          m_.get_node_pointer(),
          sizeof(${cpp_class_name}::float_type),
          format,
          ${cpp_class_name}::ndims,
          {m_.get_height(), m_.get_width()},
          {sizeof(node)*m_.get_height(), sizeof(node)}
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
    .def("run", &${cpp_class_name}::run)
    .def(
      "addBoundaryNode",
      py::overload_cast<int, int, double>(
        &${cpp_class_name}::add_boundary_node),
      "i"_a,
      "j"_a,
      "value"_a = 0.0)
    .def("getValue", &${cpp_class_name}::get_value, "i"_a, "j"_a);
''')

marchers3d = {
    'basic_marcher_3d': 'BasicMarcher3D',
    'olim6_mp0': 'Olim6Mid0',
    'olim6_mp1': 'Olim6Mid1',
    'olim6_rhr': 'Olim6Rect',
    'olim18_mp0': 'Olim18Mid0',
    'olim18_mp1': 'Olim18Mid1',
    'olim18_rhr': 'Olim18Rect',
    'olim26_mp0': 'Olim26Mid0',
    'olim26_mp1': 'Olim26Mid1',
    'olim26_rhr': 'Olim26Rect'}

marcher3d_template = Template('''
  py::class_<${cpp_class_name}>(m, "${py_class_name}")
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
    .def("run", &${cpp_class_name}::run)
    .def(
      "addBoundaryNode", &${cpp_class_name}::add_boundary_node,
      "i"_a,
      "j"_a,
      "k"_a,
      "value"_a = 0.0)
    .def("getValue", &${cpp_class_name}::get_value, "i"_a, "j"_a, "k"_a);
''')

def build_src_txt():
    src_txt = '''
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>

#include "fmm.h"

#include "basic_marcher.hpp"
#include "basic_marcher_3d.hpp"
#include "olim4.hpp"
#include "olim6.hpp"
#include "olim8.hpp"
#include "olim18.hpp"
#include "olim26.hpp"

namespace py = pybind11;
using namespace py::literals;

using speed_function = std::function<double(double, double)>;
using speed_function_3d = std::function<double(double, double, double)>;

PYBIND11_MODULE(eikonal, m) {
  m.doc() = "Testing testing";
'''
    for k in sorted(marchers.keys()):
        v = marchers[k]
        src_txt += marcher_template.substitute(
            cpp_class_name=k, py_class_name=v)
    for k in sorted(marchers3d.keys()):
        v = marchers3d[k]
        src_txt += marcher3d_template.substitute(
            cpp_class_name=k, py_class_name=v)
    src_txt += '''}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
'''
    return src_txt

if __name__ == '__main__':
    src_txt = build_src_txt()
    print(src_txt.strip())

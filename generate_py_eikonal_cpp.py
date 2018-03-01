#!/usr/bin/env python3

import argparse

from string import Template

marchers = {
    'basic_marcher': 'BasicMarcher',
    'olim4_mp0': 'Olim4Mid0',
    'olim4_rhr': 'Olim4Rect',
    'olim8_mp0': 'Olim8Mid0',
    'olim8_mp1': 'Olim8Mid1',
    'olim8_rhr': 'Olim8Rect'
}

# TODO: we should probably make the format that's computed here
# something that's computed and returned by the marcher on the C++
# side. This would make things a little simpler and cleaner here.
marcher_template = Template('''
  py::class_<${cpp_class_name}>(m, "${py_class_name}", py::buffer_protocol())
    .def_buffer([] (${cpp_class_name} & m_) -> py::buffer_info {
      auto const format =
        py::format_descriptor<${cpp_class_name}::float_type>::format();
      return {
        m_.get_node_pointer(),
        sizeof(${cpp_class_name}::float_type),
        format,
        ${cpp_class_name}::ndim,
        {m_.get_height(), m_.get_width()},
        {
          sizeof(${cpp_class_name}::node_type)*m_.get_width(),
          sizeof(${cpp_class_name}::node_type)
        }
      };
    })
    .def(
      py::init([] (py::array_t<double, py::array::c_style | py::array::forcecast> arr, double h) {
        py::buffer_info info = arr.request();
        if (info.format != py::format_descriptor<double>::format()) {
          throw std::runtime_error("Bad format: expected double array");
        }
        if (info.ndim != 2) {
          throw std::runtime_error("Expected `ndim == 2'");
        }
        if (info.shape[0] > std::numeric_limits<int>::max()) {
          throw std::runtime_error(
            "Error: size of first dimension is larger than INT_MAX");
        }
        if (info.shape[1] > std::numeric_limits<int>::max()) {
          throw std::runtime_error(
            "Error: size of second dimension is larger than INT_MAX");
        }
        int height = static_cast<int>(info.shape[0]); // height
        int width = static_cast<int>(info.shape[1]); // width
        auto m_ptr = new ${cpp_class_name} {height, width, h, no_speed_func};
        memcpy(
          (double *) m_ptr->get_s_cache_data(),
          info.ptr,
          sizeof(double)*height*width);
        return m_ptr;
      }),
      "s_cache"_a,
      "h"_a = 1.0)
    .def("run", &${cpp_class_name}::run)
    .def("step", &${cpp_class_name}::step)
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

# TODO: see comment above for `marcher_template' variable.
marcher3d_template = Template('''
  py::class_<${cpp_class_name}>(m, "${py_class_name}", py::buffer_protocol())
    .def_buffer([] (${cpp_class_name} & m_) -> py::buffer_info {
        auto const format =
          py::format_descriptor<${cpp_class_name}::float_type>::format();
        return {
          m_.get_node_pointer(),
          sizeof(${cpp_class_name}::float_type),
          format,
          ${cpp_class_name}::ndim,
          {m_.get_height(), m_.get_width(), m_.get_depth()},
          {
            sizeof(${cpp_class_name}::node_type)*m_.get_width(),
            sizeof(${cpp_class_name}::node_type),
            sizeof(${cpp_class_name}::node_type)*m_.get_width()*m_.get_height(),
          }
        };
      })
    .def(
      py::init([] (py::array_t<double, py::array::c_style | py::array::forcecast> arr, double h) {
        py::buffer_info info = arr.request();
        if (info.format != py::format_descriptor<double>::format()) {
          throw std::runtime_error("Bad format: expected double array");
        }
        if (info.ndim != 3) {
          throw std::runtime_error("Expected `ndim == 3'");
        }
        if (info.shape[0] > std::numeric_limits<int>::max()) {
          throw std::runtime_error(
            "Error: size of first dimension is larger than INT_MAX");
        }
        if (info.shape[1] > std::numeric_limits<int>::max()) {
          throw std::runtime_error(
            "Error: size of second dimension is larger than INT_MAX");
        }
        if (info.shape[2] > std::numeric_limits<int>::max()) {
          throw std::runtime_error(
            "Error: size of third dimension is larger than INT_MAX");
        }
        int height = static_cast<int>(info.shape[0]); // height
        int width = static_cast<int>(info.shape[1]); // width
        int depth = static_cast<int>(info.shape[2]); // depth
        auto m_ptr = new ${cpp_class_name} {height, width, depth, h, no_speed_func};
        memcpy(
          (double *) m_ptr->get_s_cache_data(),
          info.ptr,
          sizeof(double)*height*width*depth);
        return m_ptr;
      }),
      "s_cache"_a,
      "h"_a = 1.0)
    .def("run", &${cpp_class_name}::run)
    .def("step", &${cpp_class_name}::step)
    .def(
      "addBoundaryNode",
      py::overload_cast<int, int, int, double>(
        &${cpp_class_name}::add_boundary_node),
      "i"_a,
      "j"_a,
      "k"_a,
      "value"_a = 0.0)
    .def("getValue", &${cpp_class_name}::get_value, "i"_a, "j"_a, "k"_a);
''')

def build_src_txt(compile_all_olim3d):
    src_txt = '''
#include <limits>

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <basic_marcher.hpp>
#include <basic_marcher_3d.hpp>
#include <olim.hpp>
#include <olim3d.hpp>
'''
    if compile_all_olim3d:
        src_txt += '#include "py.olim3d.hpp"\n'
    src_txt += '''
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
    if compile_all_olim3d:
        src_txt += '''  m.def("olim3d", &olim3d_group_spec, "TODO", "s_cache"_a, "h"_a, "marcher"_a, 
    "cost_func"_a, "bd_points"_a);
}'''
    else:
        src_txt += '''
}'''
    return src_txt

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--all_olim3d', action='store_true')
    args = parser.parse_args()
    src_txt = build_src_txt(args.all_olim3d)
    print(src_txt.strip())

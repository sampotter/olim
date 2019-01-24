#!/usr/bin/env python3

import argparse

from string import Template

marchers = {
    'basic_marcher': 'BasicMarcher',
    'olim4_mp0': 'Olim4Mid0',
    'olim4_mp1': 'Olim4Mid1',
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
    .def(
      py::init([] (py::array_t<double, py::array::c_style> arr, double h) {
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
    .def(
      "add_boundary_node",
      py::overload_cast<int, int, double>(&${cpp_class_name}::add_boundary_node),
      "i"_a,
      "j"_a,
      "value"_a = 0.0)
    .def(
      "add_boundary_node",
      py::overload_cast<double, double, double, double>(
        &${cpp_class_name}::add_boundary_node),
      "i"_a,
      "j"_a,
      "s"_a,
      "value"_a = 0.0)
    .def("add_boundary_nodes", [&] (
        ${cpp_class_name} & m,
        std::vector<std::tuple<int, int, double>> const & nodes)
      {
        m.add_boundary_nodes(nodes.data(), nodes.size());
      })
    .def("set_fac_src", &${cpp_class_name}::set_fac_src, "i"_a, "j"_a, "fc"_a)
    .def("get_speed", &${cpp_class_name}::get_speed, "i"_a, "j"_a)
    .def("get_value", &${cpp_class_name}::get_value, "i"_a, "j"_a)
    .def("get_height", &${cpp_class_name}::get_height)
    .def("get_width", &${cpp_class_name}::get_width);
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
    'olim26_rhr': 'Olim26Rect',
    'olim3d_hu_mp0': 'Olim3dHuMid0',
    'olim3d_hu_mp1': 'Olim3dHuMid1',
    'olim3d_hu_rhr': 'Olim3dHuRect'}

# TODO: see comment above for `marcher_template' variable.
marcher3d_template = Template('''
py::class_<${cpp_class_name}>(m, "${py_class_name}")
    .def(
      py::init([] (py::array_t<double, py::array::f_style | py::array::forcecast> arr, double h) {
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
    .def(
      "add_boundary_node",
      py::overload_cast<int, int, int, double>(
        &${cpp_class_name}::add_boundary_node),
      "i"_a,
      "j"_a,
      "k"_a,
      "value"_a = 0.0)
    .def(
      "add_boundary_node",
      py::overload_cast<double, double, double, double, double>(
        &${cpp_class_name}::add_boundary_node),
      "i"_a,
      "j"_a,
      "k"_a,
      "s"_a,
      "value"_a = 0.0)
    .def("get_speed", &${cpp_class_name}::get_speed, "i"_a, "j"_a, "k"_a)
    .def("get_value", &${cpp_class_name}::get_value, "i"_a, "j"_a, "k"_a)
    .def("set_fac_src", &${cpp_class_name}::set_fac_src, "i"_a, "j"_a, "k"_a, "fc"_a)
    .def("get_height", &${cpp_class_name}::get_height)
    .def("get_width", &${cpp_class_name}::get_width)
    .def("get_depth", &${cpp_class_name}::get_depth);
''')

def build_src_txt(args):
    src_txt = '''
#include <limits>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <src/config.hpp>

#include <basic_marcher.hpp>
#include <basic_marcher_3d.hpp>
#include <olim.hpp>
#include <olim3d.hpp>
'''
    if args.all_olim3d:
        src_txt += '#include "py.olim3d.hpp"\n'
    src_txt += '''
namespace py = pybind11;
using namespace py::literals;

using speed_function = std::function<double(double, double)>;
using speed_function_3d = std::function<double(double, double, double)>;

PYBIND11_MODULE(pyolim, m) {
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

    if args.all_olim3d:
        src_txt += '''  m.def("olim3d", &olim3d_group_spec, "TODO", "s_cache"_a, "h"_a, "marcher"_a, 
    "cost_func"_a, "bd_points"_a);'''

    with open('pyolim_extra.cpp') as f:
        src_txt += f.read()

    return src_txt + '}'

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--all_olim3d', action='store_true')
    args = parser.parse_args()
    src_txt = build_src_txt(args)
    print(src_txt.strip())

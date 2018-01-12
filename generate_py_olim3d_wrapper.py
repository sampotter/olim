#!/usr/bin/env python3

import sys

from string import Template

cost_func_types = ['rhr', 'mp0', 'mp1']

def get_groups_type_str(i):
    return ', '.join([str((i & (1 << j)) >> j) for j in range(8)])

hpp_file_str = '''#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

pybind11::array_t<double> olim3d_group_spec(
  pybind11::array_t<
    double, 
    pybind11::array::c_style | pybind11::array::forcecast> s_cache, 
  double h, int marcher_num, std::string const & cost_func_type,
  std::vector<std::tuple<int, int, int>> const & pts);
'''

cpp_file_template = Template('''#include "py.olim3d.hpp"

#include <olim3d.hpp>

namespace py = pybind11;
using namespace py::literals;

py::array_t<double> olim3d_group_spec(
  py::array_t<double, py::array::c_style | py::array::forcecast> s_cache, 
  double h, int marcher_num, std::string const & cost_func_type,
  std::vector<std::tuple<int, int, int>> const & pts)
{
  py::buffer_info info = s_cache.request();
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
  using marcher_t = marcher_3d<node_3d>;
  marcher_t * m = nullptr;
  if (cost_func_type == "rhr") {
    switch (marcher_num) {
${rhr_cases}
    default:
      throw std::runtime_error("bad marcher: should have 0 <= num < 256");
    }
  } else if (cost_func_type == "mp0") {
    switch (marcher_num) {
${mp0_cases}
    default:
      throw std::runtime_error("bad marcher: should have 0 <= num < 256");
    }
  } else if (cost_func_type == "mp1") {
    switch (marcher_num) {
${mp1_cases}
    default:
      throw std::runtime_error("bad marcher: should have 0 <= num < 256");
    }
  } else {
    throw std::runtime_error("bad cost func: should be rhr, mp0, or mp1");
  }
  for (auto const & pt: pts) {
    int i, j, k;
    std::tie(i, j, k) = pt;
    m->add_boundary_node(i, j, k);
  }
  m->run();
  py::array_t<double> arr({height, width, depth});
  // TODO: this sucks---should replace this with a more low level memory
  // copying approach
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      for (int k = 0; k < depth; ++k) {
        arr.mutable_at(i, j, k) = m->get_value(i, j, k);
      }
    }
  }
  delete m;
  return arr;
}''')

case_template = Template('''    case ${marcher_num}:
      m = dynamic_cast<marcher_t *>(
        new olim3d_${cost_func_type}<groups_t<${groups_type_str}>>(
          height, width, depth, h, no_speed_func));
      memcpy((void *) m->get_s_cache_data(), (void *) info.ptr, sizeof(double)*height*width*depth);
      break;''')

def get_cases_for_cost_func_type(cost_func_type):
    return '\n'.join([case_template.substitute(
        marcher_num=i,
        cost_func_type=cost_func_type,
        groups_type_str=get_groups_type_str(i)) for i in range(256)])

if __name__ == '__main__':
    assert(len(sys.argv) == 3)

    hpp_file_name = sys.argv[1]
    with open(hpp_file_name, 'w') as f:
        print(hpp_file_str, file=f)

    cpp_file_name = sys.argv[2]
    with open(cpp_file_name, 'w') as f:
        cpp_file_str = cpp_file_template.substitute(
            rhr_cases=get_cases_for_cost_func_type('rhr'),
            mp0_cases=get_cases_for_cost_func_type('mp0'),
            mp1_cases=get_cases_for_cost_func_type('mp1')).strip()
        print(cpp_file_str, file=f)

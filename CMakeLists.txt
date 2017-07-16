project (fastmarching)

cmake_minimum_required (VERSION 3.0)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++1z")

find_package (gsl)

add_library (fmm
  src/abstract_marcher.cpp
  src/basic_marcher.cpp
  src/basic_marcher_3d.cpp
  src/fast_marcher.cpp
  src/fast_marcher_3d.cpp
  src/fmm_.cpp
  src/moore_marcher.cpp
  src/neumann_marcher.cpp
  src/neumann_marcher_3d.cpp
  src/node.cpp
  src/node_3d.cpp
  src/olim_util.cpp
  src/olim4_mp0.cpp
  src/olim4_rhr.cpp
  src/olim4_rhr_lut.cpp
  src/olim8_mp0.cpp
  src/olim8_mp1.cpp
  src/olim8_rhr.cpp
  src/speed_funcs.cpp)

target_link_libraries (fmm ${GSL_LIBRARIES})

set (tests
  basic_marcher
  heap
  olim_util
  olim4_rhr_lut
  olim8_mp0
  olim8_mp1
  olim8_rhr)

foreach (test ${tests})
  add_executable (${test}.test src/${test}.test.cpp)
  target_link_libraries (${test}.test fmm)
  add_test (${test} ${test}.test)
endforeach ()

add_executable (basic_marcher.speed src/basic_marcher.speed.cpp)
target_link_libraries (basic_marcher.speed fmm)

add_executable (olim4_rhr_lut.speed src/olim4_rhr_lut.speed.cpp)
target_link_libraries (olim4_rhr_lut.speed fmm)

# See this link about maybe improving the way we're doing this:
# http://public.kitware.com/pipermail/cmake/2015-April/060350.html
set (LIBFMM_PATH "${PROJECT_BINARY_DIR}/libfmm.so")
configure_file (src/fmm.in.py ${PROJECT_BINARY_DIR}/fmm.py)

add_custom_command (
  TARGET fmm
  POST_BUILD
  COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/run_mex.sh"
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_BUILD_TYPE})

enable_testing ()

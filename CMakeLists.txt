project (fastmarching)

cmake_minimum_required (VERSION 3.0)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++1z")

find_package (gsl)

add_library (fmm
  src/basic_marcher.cpp
  src/fast_marcher.cpp
  src/fmm_mex.cpp
  src/heap.cpp
  src/moore_marcher.cpp
  src/neumann_marcher.cpp
  src/node.cpp
  src/olim8_mp0.cpp
  src/olim8_mp1.cpp
  src/olim8_rhr.cpp
  src/olim8_util.cpp
  src/speed_funcs.cpp)

target_link_libraries (fmm ${GSL_LIBRARIES})

set (tests
  basic_marcher
  olim8_mp0
  olim8_mp1
  olim8_rhr
  olim8_util)

foreach (test ${tests})
  add_executable (${test}.test src/${test}.test.cpp)
  target_link_libraries (${test}.test fmm)
  add_test (${test} ${test}.test)
endforeach ()

add_custom_command (
  TARGET fmm
  POST_BUILD
  COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/run_mex.sh"
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_BUILD_TYPE})

enable_testing ()

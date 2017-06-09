project (fastmarching)

cmake_minimum_required (VERSION 3.0)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++1z")

add_library (fmm
  src/basic_marcher.cpp
  src/fast_marcher.cpp
  src/fmm_mex.cpp
  src/heap.cpp
  src/moore_marcher.cpp
  src/neumann_marcher.cpp
  src/node.cpp
  src/olim8pt_marcher.cpp)

add_executable (basic_marcher.test src/basic_marcher.test.cpp)
target_link_libraries (basic_marcher.test fmm)
add_test (basic_marcher basic_marcher.test)

add_executable (olim8pt_marcher.test src/olim8pt_marcher.test.cpp)
target_link_libraries (olim8pt_marcher.test fmm)
add_test (olim8pt_marcher olim8pt_marcher.test)

add_custom_command (
  TARGET fmm
  POST_BUILD
  COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/run_mex.sh"
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_BUILD_TYPE})

enable_testing ()

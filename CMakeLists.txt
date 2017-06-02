project (fastmarching)

cmake_minimum_required (VERSION 3.0)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++1z")

add_library (fmm src/basic_marcher.cpp src/fast_marcher.cpp
  src/fmm_mex.cpp src/heap.cpp src/node.cpp)

add_executable (main src/main.cpp)
target_link_libraries (main fmm)

add_custom_command (
  TARGET fmm
  POST_BUILD
  COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/run_mex.sh" ${CMAKE_CURRENT_SOURCE_DIR})

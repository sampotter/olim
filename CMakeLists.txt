project (fastmarching)

cmake_minimum_required (VERSION 3.0)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++1z")

add_library (grid src/grid.cpp src/heap.cpp src/node.cpp)

add_executable (main src/main.cpp)
target_link_libraries (main grid)

add_custom_command (
  TARGET grid
  POST_BUILD
  COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/run_mex.sh" ${CMAKE_CURRENT_SOURCE_DIR})

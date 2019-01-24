project (olim)

cmake_minimum_required (VERSION 3.9)

include (CheckIPOSupported)
check_ipo_supported (RESULT IPO_RESULT OUTPUT IPO_OUTPUT)
if (IPO_RESULT)
  message ("-- IPO is supported")
else ()
  if (VERBOSE)
    message (WARNING "IPO/LTO is not supported: ${IPO_OUTPUT}")
  else ()
    message ("-- IPO/LTO is not supported")
  endif ()
endif ()

option (COLLECT_STATS "Collect statistics about olimlib's performance." OFF)
option (CHECK_HEAP_PROP_IN_DEBUG "Validate heap continuously in debug mode." OFF)
option (BUILD_PYTHON_BINDINGS "Build Python bindings to libolim." ON)
option (BUILD_SCRATCH "Build a scratch executable (for testing)." OFF)
option (BUILD_TIMINGS "Build the timings executable." OFF)
option (BUILD_GEN_PLOT_DATA "Build the gen_plot_data executables." OFF)
option (BUILD_GEN_STATS "Build the gen_stats executables." OFF)
option (LINK_PROFILE "Link Google's CPU profiler." OFF)

include (CMakeDependentOption)
include (GoogleTest)
include (gtest.cmake)

list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Wall -Wextra -Werror -std=c++1z")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DOLIM_DEBUG")
set (CMAKE_CXX_FLAGS_RELWITHDEBINFO
  "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DOLIM_DEBUG -DRELWITHDEBINFO")
set (CMAKE_CXX_VISIBILITY_PRESET hidden)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  string (REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_RELWITHDEBINFO
    "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
else ()
  string (REPLACE "/DNDEBUG" "" CMAKE_CXX_FLAGS_RELWITHDEBINFO
    "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
endif ()

if (LINK_PROFILE)
  find_package (Gperftools REQUIRED)
endif ()

configure_file (src/config.hpp.in src/config.hpp)

set (OLIM_SRC_FILES
  src/basic_marcher.cpp
  src/basic_marcher_3d.cpp
  src/cost_funcs.cpp
  src/numopt.cpp
  src/slow.cpp)

add_library (olim STATIC ${OLIM_SRC_FILES})

if (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  if (IPO_RESULT)
    set_target_properties (olim PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
  endif ()
endif ()
target_include_directories (olim PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
if (LINK_PROFILE)
  target_link_libraries (olim ${GPERFTOOLS_PROFILER})
endif ()

# This executable is just for convenience. If you want to write a
# quick test and have it build with the rest of the library, you can
# do this here.
if (BUILD_SCRATCH)
  add_executable (scratch misc/cpp/scratch.cpp)
  target_include_directories (scratch PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_link_libraries (scratch olim)
endif ()

if (BUILD_TIMINGS)
  add_executable (timings misc/cpp/timings.cpp)
  target_include_directories (timings PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_link_libraries (timings olim)
endif ()

if (BUILD_GEN_PLOT_DATA)
  if (COLLECT_STATS)
    message (WARNING "Building gen_plot_data while collecting "
      "statistics! This wiill give unreliable timings.")
  endif ()

  add_executable (gen_plot_data_2d misc/cpp/gen_plot_data_2d.cpp)
  target_include_directories (gen_plot_data_2d
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_link_libraries (gen_plot_data_2d olim)

  add_executable (gen_plot_data_3d misc/cpp/gen_plot_data_3d.cpp)
  target_include_directories (gen_plot_data_3d
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_link_libraries (gen_plot_data_3d olim)

  add_executable (gen_qv_plot_data misc/cpp/gen_qv_plot_data.cpp)
  target_include_directories (gen_qv_plot_data
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_link_libraries (gen_qv_plot_data olim)
endif ()

if (BUILD_GEN_STATS)
  if (NOT COLLECT_STATS)
    message (WARNING "Building gen_stats without COLLECT_STATS=ON!")
  endif ()

  add_executable (gen_stats misc/cpp/gen_stats.cpp)
  target_include_directories (gen_stats
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_link_libraries (gen_stats olim)
endif ()

if (${BUILD_PYTHON_BINDINGS})
  # We use a Python script to generate our interface file to reduce
  # the amount of boilerplate we need to maintain manually.
  set (GENERATE_PYOLIM_CPP_SCRIPT generate_pyolim_cpp.py)
  set (GENERATED_PYOLIM_CPP_PATH
	${CMAKE_CURRENT_BINARY_DIR}/src/pyolim.cpp)
  set (GENERATE_PYOLIM_COMMAND
	./${GENERATE_PYOLIM_CPP_SCRIPT} > ${GENERATED_PYOLIM_CPP_PATH})

  add_custom_command (
	OUTPUT src/pyolim.cpp
	COMMAND ${CMAKE_COMMAND} -E copy
	  ${CMAKE_CURRENT_SOURCE_DIR}/${GENERATE_PYOLIM_CPP_SCRIPT}
      ${GENERATE_PYOLIM_CPP_SCRIPT}
    COMMAND ${CMAKE_COMMAND} -E copy
      ${CMAKE_CURRENT_SOURCE_DIR}/pyolim_extra.cpp pyolim_extra.cpp
	COMMAND ${GENERATE_PYOLIM_COMMAND}
	DEPENDS ${GENERATE_PYOLIM_CPP_SCRIPT} ${OLIM_SRC_FILES})

  # We probably aren't doing this the best way---as it is, we end up
  # recompiling all of OLIM_SRC_FILES to do this. Ideally, we want to
  # statically link the olim target into this module.
  find_package (pybind11 QUIET)
  if (NOT ${pybind11_FOUND})
	message (STATUS "Couldn't find system pybind11 install")
	if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/external/pybind11)
	  message (FATAL_ERROR "Couldn't find pybind11 submodule")
	else ()
	  message (STATUS "Using pybind11 provided by git submodule")
      add_subdirectory (external/pybind11)
	endif ()
  endif ()
  pybind11_add_module (pyolim MODULE ${GENERATED_PYOLIM_CPP_PATH})
  target_link_libraries (pyolim PRIVATE olim)
  target_include_directories (pyolim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
endif ()

# We don't use any of these test targets when we're building our
# Release build, so save time by skipping them.
if (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Release")
  # First, take care of our C++ tests.
  find_package (Threads) # required by gtest
  set (tests
	basic_marcher
	basic_marcher_3d
    bitops
    cost_funcs
	heap
    hybrid
    marcher
    marcher_3d
    numopt
    olim4_mp0
    olim4_mp1
    olim4_rhr
	olim6_mp0
	olim6_mp1
	olim6_rhr
	olim8_mp0
	olim8_mp1
	olim8_rhr
    olim18_mp0
    olim18_mp1
	olim18_rhr
    olim26_mp0
    olim26_mp1
	olim26_rhr
    olim3d_hu_mp0
    olim3d_hu_mp1
    olim3d_hu_rhr
	updates.tetra
	updates.tri)
  foreach (test ${tests})
	add_executable (${test}.test test/${test}.test.cpp)
	target_link_libraries (${test}.test olim gtest ${CMAKE_THREAD_LIBS_INIT})
	target_include_directories (
      ${test}.test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    gtest_discover_tests (${test}.test)
  endforeach ()

  # Now, add our Python tests.
  if (BUILD_PYTHON_BINDINGS)
    set (python_tests
      pyolim)
    foreach (test ${python_tests})
      configure_file (
        ${CMAKE_CURRENT_SOURCE_DIR}/test/${test}.test
        ${CMAKE_CURRENT_BINARY_DIR}/${test}.test)
      add_test (${test}.test ${CMAKE_CURRENT_BINARY_DIR}/${test}.test)
    endforeach ()
  endif ()
endif (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Release")

enable_testing ()

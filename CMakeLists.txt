project (fastmarching)

cmake_minimum_required (VERSION 3.0)

set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DEIKONAL_DEBUG")
set (CMAKE_CXX_FLAGS_RELWITHDEBINFO
  "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DEIKONAL_DEBUG")

enable_language (Fortran)

option (OLIM8_ADJ_UPDATES
  "Have the OLIM8 fast marcher do adjacent triangle updates." OFF)
option (PRINT_UPDATES
  "Output individual updates for debugging." OFF)

# If we have a macro like PRINT_UPDATES defined, because of the way
# the preprocessor works, it's helpful to turn this warning on to make
# sure that we're using PRINT_UPDATES correctly. If we don't have it
# defined, then turning this warning on just spams us during our
# builds.
if (${PRINT_UPDATES})
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wundef")
endif ()

configure_file (src/config.hpp.in src/config.hpp)
include_directories (${PROJECT_BINARY_DIR})

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++1z")

find_package (Armadillo REQUIRED)
find_package (GSL REQUIRED)

set (RPOLY_SRC_FILES src/rpoly.f src/rpoly.cpp)

add_library (rpoly ${RPOLY_SRC_FILES})

set (FMM_SRC_FILES
  src/abstract_marcher.cpp
  src/basic_marcher.cpp
  src/basic_marcher_3d.cpp
  src/conics.cpp
  src/fmm_.cpp
  src/olim_update_rules.cpp
  src/olim_util.cpp
  src/olim4_mp0.cpp
  src/olim4_rhr.cpp
  src/olim4_rhr_lut.cpp
  src/olim8_mp0.cpp
  src/olim8_rhr.cpp
  src/qroots.cpp
  src/rootfinders.cpp
  src/solim4_mp0.cpp
  src/speed_func_cache.cpp
  src/speed_funcs.cpp
  src/test_graph_marcher.cpp)

if (${APPLE})
  add_library (fmm ${FMM_SRC_FILES})
elseif (${WIN32})
  message (FATAL_ERROR "Windows unsupported")
else ()
  add_library (fmm SHARED ${FMM_SRC_FILES})
endif ()

target_link_libraries (fmm ${ARMADILLO_LIBRARIES})
target_link_libraries (fmm rpoly)
target_link_libraries (fmm ${GSL_LIBRARIES})

# We don't use any of these test targets when we're building our
# Release build, so save time by skipping them.
if (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Release")
  add_executable (rpoly.test src/rpoly.test.cpp)
  target_link_libraries (rpoly.test rpoly)
  add_test (rpoly rpoly.test)

  set (tests
	basic_marcher
	basic_marcher_3d
	conics
	heap
	olim_update_rules
	olim4_rhr_lut
	olim6_rhr
	olim8_mp0
	olim8_mp1
	olim8_rhr
	olim18_rhr
	olim26_rhr
	qroots
	solim4_mp0
	test_graph_marcher)

  foreach (test ${tests})
	add_executable (${test}.test src/${test}.test.cpp)
	target_link_libraries (${test}.test fmm)
	add_test (
	  NAME ${test} COMMAND ${CMAKE_COMMAND}
	  -DTEST_PROG=$<TARGET_FILE:${test}.test>
	  -DSOURCEDIR=${CMAKE_CURRENT_SOURCE_DIR}
	  -P ${CMAKE_CURRENT_SOURCE_DIR}/runtest.cmake)
  endforeach ()

  set (speed_tests
	basic_marcher
	olim4_rhr_lut
	olim8_mp1)

  foreach (speed_test ${speed_tests})
	add_executable (${speed_test}.speed src/${speed_test}.speed.cpp)
	target_link_libraries (${speed_test}.speed fmm)
  endforeach ()
endif (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Release")

# See this link about maybe improving the way we're doing this:
# http://public.kitware.com/pipermail/cmake/2015-April/060350.html
set (LIBFMM_PATH "${PROJECT_BINARY_DIR}/libfmm.so")
configure_file (src/fmm.in.py ${PROJECT_BINARY_DIR}/fmm.py)

add_custom_command (
  TARGET fmm
  POST_BUILD
  COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/build_mex.sh"
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_BUILD_TYPE})

enable_testing ()

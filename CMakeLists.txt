project (fastmarching)

cmake_minimum_required (VERSION 3.9 FATAL_ERROR)

include (CMakeDependentOption)
include (GoogleTest)
include (gtest.cmake)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++1z")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DEIKONAL_DEBUG")
set (CMAKE_CXX_FLAGS_RELWITHDEBINFO
  "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DEIKONAL_DEBUG -DRELWITHDEBINFO")
set (CMAKE_CXX_VISIBILITY_PRESET hidden)

option (OLIM8_ADJ_UPDATES
  "Have the OLIM8 fast marcher do adjacent triangle updates." OFF)
option (PRINT_UPDATES "Output individual updates for debugging." OFF)
option (TRIAL_NODE_OPTIMIZATION
  "Use the 'far node optimization' during updates." ON)
option (HIERARCHICAL_UPDATES "Do hierarchical updates." ON)

option (BUILD_MEX_BINDINGS "Build MATLAB MEX bindings to libeikonal." OFF)
option (BUILD_PYTHON_BINDINGS "Build Python bindings to libeikonal." ON)
option (INCLUDE_ARMDILLO "Find and include Armadillo (for development)." OFF)
option (BUILD_SCRATCH "Build a scratch executable (for testing)." OFF)

CMAKE_DEPENDENT_OPTION (COMPILE_ALL_3D_MARCHERS
  "Compile all 256 3D marchers (this takes a while!)" OFF
  "BUILD_PYTHON_BINDINGS" ON)
CMAKE_DEPENDENT_OPTION (USE_ARMADILLO_FOR_NUMOPT
  "Use Armadillo in the implementation of numerical optimization functions." OFF
  "INCLUDE_ARMADILLO" ON)

if (INCLUDE_ARMADILLO)
  find_package (Armadillo REQUIRED)
  include_directories (${ARMADILLO_INCLUDE_DIR})
endif ()

# If we have a macro like PRINT_UPDATES defined, because of the way
# the preprocessor works, it's helpful to turn this warning on to make
# sure that we're using PRINT_UPDATES correctly. If we don't have it
# defined, then turning this warning on just spams us during our
# builds.
#
# TODO: in the future, it would be better to turn this on and off
# locally to avoid missing other errors related to undefined macros.
if (${PRINT_UPDATES})
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wundef")
endif ()

configure_file (src/config.hpp.in src/config.hpp)

set (FMM_SRC_FILES
  src/abstract_marcher.cpp
  src/basic_marcher.cpp
  src/basic_marcher_3d.cpp
  src/cost_funcs.cpp
  src/fmm.cpp
  src/numopt.cpp
  src/olim_util.cpp
  src/qroots.cpp
  src/speed_funcs.cpp)

add_library (fmm STATIC ${FMM_SRC_FILES})
target_include_directories (fmm PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
if (INCLUDE_ARMADILLO)
  target_link_libraries (fmm ${ARMADILLO_LIBRARIES})
endif ()

# This executable is just for convenience. If you want to write a
# quick test and have it build with the rest of the library, you can
# do this here.
if (BUILD_SCRATCH)
  add_executable (scratch misc/scratch.cpp)
  target_include_directories (scratch PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_link_libraries (scratch fmm)
endif ()

if (${BUILD_MEX_BINDINGS})
  set (MEX_LINKER_FLAGS "-lfmm")
  if (INCLUDE_ARMADILLO)
    set (MEX_LINKER_FLAGS "${MEX_LINKER_FLAGS} -larmadillo")
  endif ()
  add_custom_command (
	OUTPUT fmm.mexmaci64 # TODO: need to build this on multiple
						 # platforms and get the correct output target
						 # name on each.
	COMMAND
	  mex -output fmm -outdir ${CMAKE_CURRENT_BINARY_DIR}
	  -L${CMAKE_CURRENT_BINARY_DIR} ${MEX_LINKER_FLAGS} src/mex.fmm.cpp
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	MAIN_DEPENDENCY src/mex.fmm.cpp
	DEPENDS ${FMM_SRC_FILES})
  add_custom_target (mex ALL DEPENDS fmm.mexmaci64)
endif ()

if (${BUILD_PYTHON_BINDINGS})
  # We use a Python script to generate our interface file to reduce
  # the amount of boilerplate we need to maintain manually.
  set (GENERATE_PY_EIKONAL_CPP_SCRIPT generate_py_eikonal_cpp.py)
  set (GENERATED_PY_EIKONAL_CPP_PATH
	${CMAKE_CURRENT_BINARY_DIR}/src/py.eikonal.cpp)
  if (${COMPILE_ALL_3D_MARCHERS})
    set (GENERATE_PY_EIKONAL_COMMAND
	  ./${GENERATE_PY_EIKONAL_CPP_SCRIPT} --all_olim3d >
      ${GENERATED_PY_EIKONAL_CPP_PATH})
  else ()
    set (GENERATE_PY_EIKONAL_COMMAND
	  ./${GENERATE_PY_EIKONAL_CPP_SCRIPT} > ${GENERATED_PY_EIKONAL_CPP_PATH})
  endif ()

  add_custom_command (
	OUTPUT src/py.eikonal.cpp
	COMMAND ${CMAKE_COMMAND} -E copy
	  ${CMAKE_CURRENT_SOURCE_DIR}/${GENERATE_PY_EIKONAL_CPP_SCRIPT}
      ${GENERATE_PY_EIKONAL_CPP_SCRIPT}
	COMMAND ${GENERATE_PY_EIKONAL_COMMAND}
	DEPENDS ${GENERATE_PY_EIKONAL_CPP_SCRIPT} ${FMM_SRC_FILES})

  # Again, use a Python script to generate a wrapper function which
  # lets us call any 3D OLIM, specified by its tetrahedra groups
  # (i.e. marcher number).
  #
  # TODO: since it seems likely that we'll be doing some code
  # generation in the future, we should probably make these two little
  # snippets into a function.
  if (${COMPILE_ALL_3D_MARCHERS})
    set (GENERATE_PY_OLIM3D_CPP_SCRIPT generate_py_olim3d_wrapper.py)
    set (GENERATED_PY_OLIM3D_HPP_PATH
      ${CMAKE_CURRENT_BINARY_DIR}/src/py.olim3d.hpp)
    set (GENERATED_PY_OLIM3D_CPP_PATH
      ${CMAKE_CURRENT_BINARY_DIR}/src/py.olim3d.cpp)
    set (GENERATE_PY_OLIM3D_COMMAND
      ./generate_py_olim3d_wrapper.py ${GENERATED_PY_OLIM3D_HPP_PATH}
        ${GENERATED_PY_OLIM3D_CPP_PATH})
    add_custom_command (
      OUTPUT src/py.olim3d.cpp
      COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/${GENERATE_PY_OLIM3D_CPP_SCRIPT}
        ${GENERATE_PY_OLIM3D_CPP_SCRIPT}
      COMMAND ${GENERATE_PY_OLIM3D_COMMAND}
      DEPENDS ${GENERATE_PY_OLIM3D_CPP_SCRIPT} ${FMM_SRC_FILES})
  endif ()

  # We probably aren't doing this the best way---as it is, we end up
  # recompiling all of FMM_SRC_FILES to do this. Ideally, we want to
  # statically link the fmm target into this module.
  find_package (pybind11 REQUIRED)
  if (NOT ${pybind11_FOUND})
    add_subdirectory (external/pybind11)
  endif ()
  if (${COMPILE_ALL_3D_MARCHERS})
    pybind11_add_module (eikonal MODULE
      ${GENERATED_PY_EIKONAL_CPP_PATH}
      ${GENERATED_PY_OLIM3D_CPP_PATH})
  else ()
    pybind11_add_module (eikonal MODULE ${GENERATED_PY_EIKONAL_CPP_PATH})
  endif ()

  target_link_libraries (eikonal PRIVATE fmm)
  target_include_directories (eikonal PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
endif ()

# We don't use any of these test targets when we're building our
# Release build, so save time by skipping them.
if (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Release")
  find_package (Threads) # required by gtest

  set (tests
	basic_marcher
	basic_marcher_3d
    cost_funcs
	heap
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
	qroots
	update_rules.tetra_updates
	update_rules.tri_updates)

  foreach (test ${tests})
	add_executable (${test}.test test/${test}.test.cpp)
	target_link_libraries (${test}.test fmm gtest ${CMAKE_THREAD_LIBS_INIT})
	target_include_directories (
      ${test}.test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    gtest_discover_tests (${test}.test)
  endforeach ()
endif (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Release")

enable_testing ()

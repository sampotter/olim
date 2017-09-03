project (fastmarching)

cmake_minimum_required (VERSION 3.0)

enable_language (Fortran)

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

add_executable (rpoly.test src/rpoly.test.cpp)
target_link_libraries (rpoly.test rpoly)
add_test (rpoly rpoly.test)

set (tests
  basic_marcher
  basic_marcher_3d
  conics
  heap
  olim4_rhr_lut
  olim6_rhr
  olim8_mp0
  olim8_mp1
  olim8_rhr
  qroots
  solim4_mp0
  test_graph_marcher)

foreach (test ${tests})
  add_executable (${test}.test src/${test}.test.cpp)
  target_link_libraries (${test}.test fmm)
  add_test (${test} ${test}.test)
endforeach ()

set (speed_tests
  basic_marcher
  olim4_rhr_lut
  olim8_mp1)

foreach (speed_test ${speed_tests})
  add_executable (${speed_test}.speed src/${speed_test}.speed.cpp)
  target_link_libraries (${speed_test}.speed fmm)
endforeach ()

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

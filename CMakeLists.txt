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
option (BUILD_SCRATCH "Build a scratch executable (for testing)." OFF)
option (BUILD_GEN_STATS "Build the gen_stats executables." OFF)
option (LINK_PROFILE "Link Google's CPU profiler." OFF)

include (CMakeDependentOption)

if (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Release")
  include (GoogleTest)
  include (gtest.cmake)
endif ()

list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Wall -Wextra -Werror")
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
  src/olim_wrapper.cpp)

add_library (olim STATIC ${OLIM_SRC_FILES})
target_compile_features (olim PRIVATE cxx_std_17)

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
  target_compile_features (scratch PRIVATE cxx_std_17)
  target_include_directories (scratch PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_link_libraries (scratch olim)
endif ()

if (BUILD_GEN_STATS)
  if (NOT COLLECT_STATS)
    message (WARNING "Building gen_stats without COLLECT_STATS=ON!")
  endif ()
  add_executable (gen_stats misc/cpp/gen_stats.cpp)
  target_compile_features (gen_stats PRIVATE cxx_std_17)
  target_include_directories (gen_stats
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_link_libraries (gen_stats olim)
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
    range
	updates.tetra
	updates.tri
    vec)
  foreach (test ${tests})
	add_executable (${test}.test test/${test}.test.cpp)
    target_compile_features (${test}.test PRIVATE cxx_std_17)
	target_link_libraries (${test}.test olim gtest ${CMAKE_THREAD_LIBS_INIT})
	target_include_directories (
      ${test}.test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    gtest_discover_tests (${test}.test)
  endforeach ()
endif ()

enable_testing ()

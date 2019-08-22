cmake_minimum_required (VERSION 3.9)

project (olim)

enable_testing ()

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

option (CHECK_HEAP_PROP_IN_DEBUG "Validate heap continuously in debug mode." OFF)
option (BUILD_SCRATCH "Build a scratch executable (for testing)." OFF)
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
  src/eikonal.cpp
  src/eikonal/func.cpp
  src/numopt.cpp
  src/quasipot.cpp)

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
  add_executable (scratch src/scratch.cpp)
  target_compile_features (scratch PRIVATE cxx_std_17)
  target_include_directories (scratch PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_link_libraries (scratch olim)
endif ()

# We don't use any of these test targets when we're building our
# Release build, so save time by skipping them.
if (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Release")
  # First, take care of our C++ tests.
  find_package (Threads) # required by gtest
  add_subdirectory (test)
  add_subdirectory (test/eikonal)
endif ()

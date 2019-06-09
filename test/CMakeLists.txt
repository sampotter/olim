set (tests
  bitops
  heap
  hybrid
  marcher
  marcher_3d
  numopt
  ordering
  range
  vec)

foreach (test ${tests})
  add_executable (${test}.test ${test}.test.cpp)
  target_compile_features (${test}.test PRIVATE cxx_std_17)
  target_link_libraries (${test}.test olim gtest ${CMAKE_THREAD_LIBS_INIT})
  target_include_directories (${test}.test PRIVATE ${CMAKE_SOURCE_DIR}/src)
  gtest_discover_tests (${test}.test)
endforeach ()

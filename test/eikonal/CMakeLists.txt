set (tests
  fmm
  func
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
  tetra
  tri)

foreach (test ${tests})
  add_executable (${test}.test ${CMAKE_CURRENT_SOURCE_DIR}/${test}.test.cpp)
  target_compile_features (${test}.test PRIVATE cxx_std_17)
  target_link_libraries (${test}.test olim gtest ${CMAKE_THREAD_LIBS_INIT})
  target_include_directories (${test}.test PRIVATE ${CMAKE_SOURCE_DIR}/src)
  gtest_discover_tests (${test}.test)
endforeach ()

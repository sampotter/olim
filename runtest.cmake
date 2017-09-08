execute_process (
  COMMAND ${TEST_PROG}
  OUTPUT_FILE tmp.txt
  RESULT_VARIABLE HAD_ERROR)

if (HAD_ERROR)
  message (FATAL_ERROR "error running test")
endif ()

execute_process(
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/testcmp.py
  INPUT_FILE tmp.txt
  RESULT_VARIABLE TEST_FAILURE)

if (TEST_FAILURE)
  message (FATAL_ERROR "test failure")
endif ()

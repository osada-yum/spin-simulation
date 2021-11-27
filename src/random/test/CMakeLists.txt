set(TEST_LCG test_lcg.out)

add_executable(${TEST_LCG}
  test_lcg.f90
  )
target_include_directories(${TEST_LCG}
  PRIVATE
  ${CMAKE_BINARY_DIR}/src
  ${CMAKE_BINARY_DIR}/src/random
  )
target_link_libraries(${TEST_LCG}
  PRIVATE
  ${LIB_RANDOM}
  )
add_test(NAME LCG_test
  COMMAND $<TARGET_FILE:${TEST_LCG}>
  )

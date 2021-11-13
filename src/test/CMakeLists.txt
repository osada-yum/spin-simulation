set(TEST_2D-Ising    Ising2d_test.out)

add_executable(${TEST_2D-Ising}
  Ising2d_test.f90
  )
target_include_directories(${TEST_2D-Ising}
  PRIVATE
  ${CMAKE_CURRENT_BINARY_DIR}/../
  )
target_link_libraries(${TEST_2D-Ising}
  PRIVATE
  ${LIB_2D-Ising}
  )
add_test(NAME Ising2d_test
  COMMAND $<TARGET_FILE:${TEST_2D-Ising}>
  )
